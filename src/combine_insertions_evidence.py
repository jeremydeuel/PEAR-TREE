# PEAR-TREE - paired ends of aberrant retrotransposons in phylogenetic trees
#
# Copyright (C) 2025 Jeremy Deuel <jeremy.deuel@usz.ch>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Pooled per-patient junction evidence for combine_insertions (plans/tprt_hallmarks/SPEC.md).

Discovery may write, per sample, `<sample>.evidence.tsv.gz` next to `<sample>.txt.gz`
(one row per evidence read incl. mates). This module

  * loads those sidecars for the loci that survived intersect_insertions,
  * applies the SPEC "Independence rule (combine)": collapse by (sample, frag), merge
    within-sample PCR/optical duplicates from coordinates + sequence (NEVER the 0x400
    flag), keep cross-sample fragments independent unless exactly identical (flagged),
    and counts `n_independent` per junction (LEFT, RIGHT -- the poly-A end is one of them),
  * builds the indel-aware clip consensus (src/indel_consensus.py) per junction,
  * writes `<patient>.insertions.evidence.tsv.gz` and `<patient>.insertions.reads.fa.gz`.

Nothing here runs when no sidecar exists, so the legacy outputs stay byte-identical.
"""

import gzip
import os
from collections import Counter, defaultdict
from typing import Dict, List, Optional

import edlib

from indel_consensus import ClipRead, ConsensusResult, indel_aware_consensus, revcomp
from quality_seq import QualitySeq

EVIDENCE_ROLES = ("CLIP", "POLYA", "DISC", "SPAN")  # MATE rows only attach to a fragment
_ROLE_PRIORITY = {"CLIP": 0, "POLYA": 1, "DISC": 2, "SPAN": 3, "MATE": 9}
SIDES = ("LEFT", "RIGHT")

EVIDENCE_TSV_COLUMNS = [
    # SPEC columns (order fixed)
    "insertion_id", "side", "n_reads", "n_fragments", "n_independent", "n_samples", "n_mates",
    "supported", "clip_consensus", "consensus_depth", "polya_len_median", "polya_len_range",
    "beyond_polya", "beyond_polya_support",
    # additive diagnostics (appended; see SPEC)
    "polya_end", "n_duplicates", "n_cross_sample_identical", "fail_reason", "consensus_stop",
]


def sidecar_path(txt_gz: str) -> str:
    """`<sample>.txt.gz` -> its evidence sidecar.

    Rust discovery writes `<sample>.txt.gz.evidence.tsv.gz` (the cluster wrappers rename
    `<out>.<ext>` files); `<sample>.evidence.tsv.gz` is also accepted. Returns the existing
    one, else the Rust name."""
    rust = txt_gz + ".evidence.tsv.gz"
    stem = txt_gz[:-7] if txt_gz.endswith(".txt.gz") else txt_gz
    alt = stem + ".evidence.tsv.gz"
    if not os.path.exists(rust) and os.path.exists(alt):
        return alt
    return rust


def sample_name(txt_gz: str) -> str:
    b = os.path.basename(txt_gz)
    return b[:-7] if b.endswith(".txt.gz") else b


class EvidenceRow:
    __slots__ = ("sample", "locus", "side", "role", "frag", "r12", "flag", "ref", "pos", "strand",
                 "outer", "mref", "mpos", "mstrand", "tlen", "mapq", "cigar", "clip_at", "seq", "qual")

    def __init__(self, sample, d):
        self.sample = sample
        self.locus = d["locus"]
        self.side = d["side"]
        self.role = d["role"]
        self.frag = d["frag"]
        self.r12 = _int(d.get("r12"), 0)
        self.flag = _int(d.get("flag"), 0)
        self.ref = d.get("ref", "*")
        self.pos = _int(d.get("pos"), -1)
        self.strand = d.get("strand", "*")
        self.outer = _int(d.get("outer"), -1)
        self.mref = d.get("mref", "*")
        self.mpos = _int(d.get("mpos"), -1)
        self.mstrand = d.get("mstrand", "*")
        self.tlen = _int(d.get("tlen"), 0)
        self.mapq = _int(d.get("mapq"), 0)
        self.cigar = d.get("cigar", "*")
        self.clip_at = _int(d.get("clip_at"), -1)
        self.seq = d.get("seq", "")
        self.qual = d.get("qual", "")

    @property
    def mapped(self) -> bool:
        return self.ref not in ("*", "") and self.pos >= 0 and not (self.flag & 0x4)

    @property
    def mate_mapped(self) -> bool:
        return self.mref not in ("*", "") and self.mpos >= 0 and not (self.flag & 0x8)

    def quals(self) -> List[int]:
        if len(self.qual) == len(self.seq):
            return [ord(c) - 33 for c in self.qual]
        return [30] * len(self.seq)

    def outward_clip(self):
        """The junction clip of a CLIP row, oriented outward from the junction
        (first base adjacent to the junction), with qualities. RIGHT = [aligned][clip],
        LEFT = [clip][aligned] in reference-forward orientation (as in discovery)."""
        n = len(self.seq)
        at = self.clip_at
        if not (0 <= at <= n):
            at = _clip_at_from_cigar(self.cigar, self.side, n)
        q = self.quals()
        if self.side == "RIGHT":
            return self.seq[at:], q[at:]
        return revcomp(self.seq[:at]), q[:at][::-1]


def _int(x, default):
    try:
        return int(x)
    except (TypeError, ValueError):
        return default


def _clip_at_from_cigar(cigar: str, side: str, n: int) -> int:
    import re
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar or "")
    if not ops:
        return 0 if side == "LEFT" else n
    if side == "LEFT":
        return int(ops[0][0]) if ops[0][1] == "S" else 0
    return n - int(ops[-1][0]) if ops[-1][1] == "S" else n


def load_evidence(input_files: List[str], wanted_loci) -> (Dict, set):
    """Stream every existing sidecar; keep rows whose locus is in `wanted_loci`.
    Returns ({(locus, side): [EvidenceRow]}, {basename of txt.gz with a sidecar})."""
    rows = defaultdict(list)
    have = set()
    for f in input_files:
        sp = sidecar_path(f)
        if not os.path.exists(sp):
            continue
        have.add(os.path.basename(f))
        sample = sample_name(f)
        with gzip.open(sp, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            li = header.index("locus")
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) != len(header) or p[li] not in wanted_loci:
                    continue
                r = EvidenceRow(sample, dict(zip(header, p)))
                rows[(r.locus, r.side)].append(r)
    return rows, have


# ------------------------------------------------------------------ fragments

class Fragment:
    """All records of one template (sample, frag) at one junction."""
    __slots__ = ("sample", "frag", "rows", "primary", "mate")

    def __init__(self, sample, frag, rows):
        self.sample = sample
        self.frag = frag
        self.rows = rows
        ev = [r for r in rows if r.role != "MATE"]
        self.primary = min(ev, key=lambda r: (_ROLE_PRIORITY.get(r.role, 5), r.r12)) if ev else None
        self.mate = None
        if self.primary is not None:
            others = [r for r in rows if r.r12 != self.primary.r12]
            if others:
                self.mate = min(others, key=lambda r: (0 if r.role == "MATE" else 1, r.r12))

    @property
    def strand(self):
        return self.primary.strand

    @property
    def outer(self):
        return self.primary.outer

    def mate_coord(self):
        """(ref, pos, strand) of the mate, or None if the mate is unmapped/unknown.
        Uses the mate record's own `outer` when present (both fragments must then have
        one -- see _is_dup), else the primary record's mpos."""
        p = self.primary
        if not p.mate_mapped:
            return None
        return (p.mref, p.mpos, p.mstrand)

    def mate_outer(self):
        if self.mate is not None and self.mate.mapped and self.mate.outer >= 0:
            return (self.mate.ref, self.mate.outer, self.mate.strand)
        return None

    def clip_seq(self):
        p = self.primary
        if p.role == "CLIP":
            return p.outward_clip()[0].upper()
        return p.seq.upper()

    def mate_seq(self):
        return self.mate.seq.upper() if self.mate is not None else ""


def collapse_fragments(rows: List[EvidenceRow]) -> List[Fragment]:
    """SPEC rule 1: one fragment per (sample, frag). Fragments seen only as MATE rows
    (no evidence role) are not evidence and are dropped."""
    by = defaultdict(list)
    for r in rows:
        by[(r.sample, r.frag)].append(r)
    frags = [Fragment(s, f, rs) for (s, f), rs in by.items()]
    return sorted([f for f in frags if f.primary is not None], key=lambda f: (f.sample, f.outer, f.frag))


def _seq_close(a: str, b: str, max_edit: int) -> bool:
    m = min(len(a), len(b))
    if m == 0:
        return True
    a, b = a[:m], b[:m]
    if a == b:
        return True
    return edlib.align(a, b, mode="NW", task="distance", k=max_edit)["editDistance"] != -1


def _is_dup(a: Fragment, b: Fragment, tol: int, max_edit: int) -> bool:
    """SPEC rule 2 (same sample only)."""
    if a.strand != b.strand or abs(a.outer - b.outer) > tol:
        return False
    ao, bo = a.mate_outer(), b.mate_outer()
    if ao is not None and bo is not None:
        am, bm = ao, bo
    else:
        am, bm = a.mate_coord(), b.mate_coord()
    if am is not None and bm is not None:
        # outer coordinates match on both ends -> duplicate
        return am[0] == bm[0] and am[2] == bm[2] and abs(am[1] - bm[1]) <= tol
    if (am is None) != (bm is None):
        return False
    # no mate placement on either: same read outer AND junction-clip sequence within
    # max_edit; when both mate sequences are known they must agree too (a different
    # mate => a different molecule => independent).
    if not _seq_close(a.clip_seq(), b.clip_seq(), max_edit):
        return False
    ams, bms = a.mate_seq(), b.mate_seq()
    if ams and bms:
        return _seq_close(ams, bms, max_edit)
    return True


def _identical(a: Fragment, b: Fragment) -> bool:
    """SPEC rule 3: exact identity of both outer coordinates AND sequence."""
    if a.strand != b.strand or a.outer != b.outer:
        return False
    am = a.mate_outer() or a.mate_coord()
    bm = b.mate_outer() or b.mate_coord()
    if am != bm:
        return False
    if a.primary.seq.upper() != b.primary.seq.upper():
        return False
    return a.mate_seq() == b.mate_seq()


def independent_clusters(frags: List[Fragment], tol: int = 2, max_edit: int = 2):
    """Union-find over fragments. Returns (clusters, n_duplicates, n_cross_identical)."""
    parent = list(range(len(frags)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    n_dup = n_cross = 0
    by_sample = defaultdict(list)
    for i, f in enumerate(frags):
        by_sample[f.sample].append(i)
    for idx in by_sample.values():
        idx.sort(key=lambda i: frags[i].outer)
        for a_pos, i in enumerate(idx):
            for j in idx[a_pos + 1:]:
                if frags[j].outer - frags[i].outer > tol:
                    break
                if find(i) != find(j) and _is_dup(frags[i], frags[j], tol, max_edit):
                    parent[find(j)] = find(i)
                    n_dup += 1
    # cross-sample: exact identity only
    by_key = defaultdict(list)
    for i, f in enumerate(frags):
        by_key[(f.strand, f.outer)].append(i)
    for idx in by_key.values():
        for a_pos, i in enumerate(idx):
            for j in idx[a_pos + 1:]:
                if frags[i].sample != frags[j].sample and find(i) != find(j) and _identical(frags[i], frags[j]):
                    parent[find(j)] = find(i)
                    n_cross += 1
    groups = defaultdict(list)
    for i in range(len(frags)):
        groups[find(i)].append(frags[i])
    return list(groups.values()), n_dup, n_cross


# ------------------------------------------------------------------ junctions

class JunctionRecord:
    """Pooled stats + consensus for one insertion junction (one evidence.tsv row)."""

    def __init__(self, insertion_id, side):
        self.insertion_id = insertion_id
        self.side = side
        self.rows: List[EvidenceRow] = []
        self.n_reads = self.n_fragments = self.n_independent = self.n_samples = self.n_mates = 0
        self.n_duplicates = self.n_cross = 0
        self.supported = "NA"
        self.consensus: ConsensusResult = ConsensusResult()          # junction reads + mates (TSV)
        self.combined_consensus: ConsensusResult = ConsensusResult()  # junction reads only (combined.txt.gz)
        self.polya_end = 0
        self.fail_reason = ""
        self.aligned = ""   # reference part, in clip_consensus convention (uppercase)

    def tsv(self) -> str:
        c = self.consensus
        clip = c.seq.lower()
        depth = list(c.depth)
        beyond = c.beyond_polya
        if self.side == "LEFT":   # reference-forward: [clip][aligned]
            clip_cons = revcomp(clip) + self.aligned
            depth = depth[::-1]
            beyond = revcomp(beyond)
        else:                     # [aligned][clip]
            clip_cons = self.aligned + clip
        vals = [self.insertion_id, self.side, self.n_reads, self.n_fragments, self.n_independent,
                self.n_samples, self.n_mates, self.supported, clip_cons,
                ",".join(str(d) for d in depth),
                "" if c.polya_len_median is None else c.polya_len_median, c.polya_len_range,
                beyond.lower(), c.beyond_polya_support if beyond else 0,
                self.polya_end, self.n_duplicates, self.n_cross, self.fail_reason, c.stop_reason]
        return "\t".join(str(v) for v in vals) + "\n"


def evaluate_junction(insertion_id, side, rows, cfg) -> JunctionRecord:
    tol = cfg.get("dup_coord_tolerance", 2)
    max_edit = cfg.get("dup_max_edit", 2)
    min_ind = cfg.get("min_independent_fragments", 2)
    polya_min = cfg.get("polya_min_len", 8)
    rec = JunctionRecord(insertion_id, side)
    rec.rows = rows
    rec.n_mates = sum(1 for r in rows if r.role == "MATE")
    rec.n_reads = len(rows) - rec.n_mates
    frags = collapse_fragments(rows)
    rec.n_fragments = len(frags)
    rec.n_samples = len({f.sample for f in frags})
    clusters, rec.n_duplicates, rec.n_cross = independent_clusters(frags, tol, max_edit)
    rec.n_independent = len(clusters)
    # consensus input: every read, weighted 1/|cluster| so a PCR family votes once;
    # `group` = cluster so depth counts independent fragments.
    reads = []
    for gi, cl in enumerate(clusters):
        w = 1.0 / len(cl)
        for f in cl:
            for r in f.rows:
                if r.role == "CLIP":
                    s, q = r.outward_clip()
                    reads.append(ClipRead(s, q, gi, w, True))
                elif r.seq:
                    reads.append(ClipRead(r.seq, r.quals(), gi, w, False))
    rec.consensus = indel_aware_consensus(reads, min_depth=min_ind, polya_min_len=polya_min)
    if cfg.get("indel_aware_consensus", False):
        # combined.txt.gz gets the junction-read-only consensus: mate extension can run
        # through a short insertion into the flank on the far side, and that flank would
        # make the clipped-remap filter (clip maps near the breakpoint) discard a real call.
        anchored = [r for r in reads if r.anchored]
        rec.combined_consensus = (rec.consensus if len(anchored) == len(reads) else
                                  indel_aware_consensus(anchored, min_depth=min_ind, polya_min_len=polya_min))
    rec.polya_end = int(rec.consensus.polya_base is not None
                        or any(r.role == "POLYA" for r in rows)
                        or _clips_start_with_polya([r for r in rows if r.role == "CLIP"], polya_min))
    return rec


def _clips_start_with_polya(clip_rows, polya_min, within=5) -> bool:
    """Fallback poly-A-end call that works even below consensus depth (so a single-fragment
    poly-A end is still LABELLED as such when the gate drops it): >= half of the junction
    clips carry an A/T run >= polya_min starting within `within` bases of the junction."""
    if not clip_rows:
        return False
    runs = ("A" * polya_min, "T" * polya_min)
    hits = 0
    for r in clip_rows:
        head = r.outward_clip()[0].upper()[:within + polya_min]
        hits += any(0 <= head.find(x) <= within for x in runs)
    return 2 * hits >= len(clip_rows)


def _aligned_part(ins, side) -> str:
    if side == "RIGHT":
        a = getattr(ins, "right_aligned", None)
        return str(a.revcomp()).upper() if a is not None else ""
    a = getattr(ins, "left_aligned", None)
    return str(a).upper() if a is not None else ""


def apply_evidence(insertions, input_files, cfg):
    """Evaluate every junction of every insertion from the pooled sidecars.

    Returns (kept, records, failed_names, stats) or None when no sidecar exists (the
    caller then behaves exactly as before). `records` maps insertion name -> [JunctionRecord];
    gated-out insertions keep their records (supported=0) for diagnostics.
    Mutates kept insertions' clips when `indel_aware_consensus` is on and the new
    consensus is at least as long as the legacy longest-clip choice."""
    wanted = {i.name for i in insertions}
    rows, have = load_evidence(input_files, wanted)
    if not have:
        return None
    gate = bool(cfg.get("require_independent_fragments", False))
    min_ind = cfg.get("min_independent_fragments", 2)
    use_cons = bool(cfg.get("indel_aware_consensus", False))
    missing = [os.path.basename(f) for f in input_files if os.path.basename(f) not in have]
    if missing:
        print(f"WARNING: {len(missing)} discovery file(s) have no evidence sidecar; insertions "
              f"they contribute to are not gated (supported=NA): {','.join(missing[:5])}"
              f"{' ...' if len(missing) > 5 else ''}")
    kept, failed, records = [], set(), {}
    reasons = Counter()
    n_replaced = 0
    for ins in insertions:
        recs = []
        for side in SIDES:
            rec = evaluate_junction(ins.name, side, rows.get((ins.name, side), []), cfg)
            rec.aligned = _aligned_part(ins, side)
            recs.append(rec)
        records[ins.name] = recs
        eligible = all(f in have for f in ins.files)
        fails = [r for r in recs if r.n_independent < min_ind]
        for r in recs:
            if eligible:
                r.supported = 0 if r in fails else 1
                if r in fails:
                    r.fail_reason = f"n_independent<{min_ind}"
        if gate and eligible and fails:
            failed.add(ins.name)
            if all(r.n_reads == 0 for r in recs):
                reasons["no_evidence_reads"] += 1
            else:
                reasons["+".join(r.side + ("(polyA)" if r.polya_end else "") for r in fails)] += 1
            continue
        if use_cons:
            n_replaced += _replace_clips(ins, recs)
        kept.append(ins)
    print(f"evidence sidecars: {len(have)}/{len(input_files)} files; {len(insertions)} insertions evaluated, "
          f"{sum(len(v) for v in rows.values())} evidence rows")
    if gate:
        print(f"independent-fragment gate (>= {min_ind} per junction, pooled over samples): "
              f"kept {len(kept)}, dropped {len(failed)}")
        for reason, n in sorted(reasons.items()):
            print(f"  dropped: failing junction(s) {reason}: {n}")
    if use_cons:
        print(f"indel-aware consensus replaced {n_replaced} junction clip(s) in combined output")
    return kept, records, failed, {"reasons": reasons, "replaced": n_replaced}


def _replace_clips(ins, recs) -> int:
    n = 0
    for rec in recs:
        c = rec.combined_consensus
        if not c.seq:
            continue
        attr = "right_clipped" if rec.side == "RIGHT" else "left_clipped"
        old = getattr(ins, attr, None)
        if old is None:
            continue
        if len(c.seq) < len(old) or c.seq.upper() == str(old).upper():
            continue
        seq = c.seq.lower() if str(old).islower() else c.seq.upper()
        setattr(ins, attr, QualitySeq(seq, list(c.score)))
        n += 1
    return n


def write_evidence_outputs(records: Dict[str, list], names, evidence_tsv: str, reads_fa: str):
    """Write `<patient>.insertions.evidence.tsv.gz` + `<patient>.insertions.reads.fa.gz`
    for the insertions in `names` (in that order)."""
    n_rows = n_reads = 0
    with gzip.open(evidence_tsv, "wt") as t, gzip.open(reads_fa, "wt") as fa:
        t.write("\t".join(EVIDENCE_TSV_COLUMNS) + "\n")
        for name in names:
            for rec in records.get(name, ()):
                t.write(rec.tsv())
                n_rows += 1
                for r in rec.rows:
                    fa.write(f">{name}|{rec.side}|{r.role}|{r.sample}|{r.frag}|{r.r12}\n{r.seq}\n")
                    n_reads += 1
    print(f"wrote {n_rows} junction rows -> {evidence_tsv}; {n_reads} reads -> {reads_fa}")
