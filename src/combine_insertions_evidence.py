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

EVIDENCE_ROLES = ("CLIP", "POLYA", "DISC", "SPAN", "SHORT")  # MATE rows only attach to a fragment
_ROLE_PRIORITY = {"CLIP": 0, "POLYA": 1, "DISC": 2, "SPAN": 3, "SHORT": 4, "MATE": 9}
SIDES = ("LEFT", "RIGHT")

EVIDENCE_TSV_COLUMNS = [
    # SPEC columns (order fixed)
    "insertion_id", "side", "n_reads", "n_fragments", "n_independent", "n_samples", "n_mates",
    "supported", "clip_consensus", "consensus_depth", "polya_len_median", "polya_len_range",
    "beyond_polya", "beyond_polya_support",
    # additive diagnostics (appended; see SPEC)
    "polya_end", "n_duplicates", "n_cross_sample_identical", "fail_reason", "consensus_stop",
    "n_dup_coord", "n_dup_seq", "member_loci",
    "n_short_used", "n_short_rejected", "n_short_mate_inside", "n_independent_no_short",
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
    """Stream every existing sidecar; keep rows whose locus is in `wanted_loci` -- either a
    set of locus ids (any file) or a set of (txt.gz basename, locus) pairs (that file only).
    Returns ({(locus, side): [EvidenceRow]}, {basename of txt.gz with a sidecar}); with
    (file, locus) pairs the key is ((file, locus), side)."""
    rows = defaultdict(list)
    have = set()
    for f in input_files:
        sp = sidecar_path(f)
        if not os.path.exists(sp):
            continue
        fb = os.path.basename(f)
        have.add(fb)
        sample = sample_name(f)
        by_file = any(isinstance(w, tuple) for w in wanted_loci)
        want = {l for (ff, l) in wanted_loci if ff == fb} if by_file else wanted_loci
        with gzip.open(sp, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            li = header.index("locus")
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) != len(header) or p[li] not in want:
                    continue
                r = EvidenceRow(sample, dict(zip(header, p)))
                rows[((fb, r.locus) if by_file else r.locus, r.side)].append(r)
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


class DedupParams:
    """Lenient within-sample duplicate rule (SPEC "Independence rule"; 0x400 reads are already
    dropped by discovery -- this is the SECOND dedup, for PCR/optical copies markdup missed).

    tol        `dup_coord_tolerance` (5): read outer, mate outer/mpos and the junction (clip)
               position may each differ by up to tol bp
    max_edit   `dup_max_edit` (3) / `dup_max_edit_frac` (0.02): edit budget for the lenient
               sequence match = max(max_edit, ceil(frac * compared length))
    polya_min  `polya_min_len` (8): sequences are compared homopolymer-compressed and cut after
               the first A/T run >= polya_min (sequence 3' of a long poly-A is SBS junk)
    mate_min_mapq `dup_mate_min_mapq` (20): a mate record below this MAPQ (a mate inside the
               element, multi-mapped onto a random paralog) is not a placement -- the pair is then
               judged by mate SEQUENCE (E2E: 206 missed duplicates were exactly this)
    A read whose 5' end lies in the junction soft clip (LEFT clip on '+', RIGHT clip on '-') has
    an `outer` that moves with the clip length (poly-A jitter); for those the junction position
    replaces the outer-coordinate test."""
    __slots__ = ("tol", "max_edit", "max_edit_frac", "polya_min", "mate_min_mapq")

    def __init__(self, tol=5, max_edit=3, max_edit_frac=0.02, polya_min=8, mate_min_mapq=20):
        self.tol, self.max_edit, self.max_edit_frac, self.polya_min = tol, max_edit, max_edit_frac, polya_min
        self.mate_min_mapq = mate_min_mapq

    @classmethod
    def from_cfg(cls, cfg):
        return cls(cfg.get("dup_coord_tolerance", 5), cfg.get("dup_max_edit", 3),
                   cfg.get("dup_max_edit_frac", 0.02), cfg.get("polya_min_len", 8),
                   cfg.get("dup_mate_min_mapq", 20))

    def budget(self, n):
        return max(self.max_edit, int(-(-self.max_edit_frac * n // 1)))


def _hp_compress_cut(seq: str, polya_min: int) -> str:
    """Homopolymer-compressed sequence, cut right after the first A/T run >= polya_min (the
    run itself is kept as one symbol): poly-A length jitter and the low-quality sequence that
    follows a long poly-A on Illumina never count as differences."""
    out = []
    i, n = 0, len(seq)
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        out.append(seq[i])
        if seq[i] in "AT" and j - i >= polya_min:
            break
        i = j
    return "".join(out)


def _prefix_close(a: str, b: str, p: DedupParams) -> bool:
    """Homopolymer-compressed, poly-A-cut, common-prefix comparison within the edit budget."""
    a = _hp_compress_cut(a, p.polya_min)
    b = _hp_compress_cut(b, p.polya_min)
    m = min(len(a), len(b))
    if m == 0:
        return True
    a, b = a[:m], b[:m]
    if a == b:
        return True
    return edlib.align(a, b, mode="NW", task="distance", k=p.budget(m))["editDistance"] != -1


def _seq_close(a: str, b: str, p: DedupParams, shift: int = 0) -> bool:
    """Lenient sequence identity. `shift` > 0: the two reads may start up to `shift` bases
    apart (unanchored reads, e.g. an unmapped mate of a duplicate); each offset 0..shift of
    either read is tried, the comparison itself stays strict (edit budget)."""
    if not shift:            # anchored at the junction
        return _prefix_close(a, b, p)
    for k in range(shift + 1):
        if _prefix_close(a[k:], b, p) or (k and _prefix_close(a, b[k:], p)):
            return True
    return False


def _junction_pos(f: Fragment):
    """Reference position of the junction (clip) of a CLIP fragment, else None."""
    p = f.primary
    if p.role != "CLIP" or not p.mapped:
        return None
    if p.side == "LEFT":
        return p.pos
    import re
    ref_len = sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", p.cigar or ""))
    return p.pos + ref_len


def _dup_seqs(a: Fragment, b: Fragment):
    """Sequences compared by the lenient dedup: the outward junction clips when both
    fragments are CLIP fragments, else the two full primary reads (DISC / POLYA / SHORT, or a
    SHORT read against a CLIP read of the same molecule)."""
    if a.primary.role == "CLIP" and b.primary.role == "CLIP":
        return a.clip_seq(), b.clip_seq()
    return a.primary.seq.upper(), b.primary.seq.upper()


def _clip_shift(a: Fragment, b: Fragment, p: DedupParams) -> int:
    """Junction clips are anchored (shift 0); a whole-read comparison (DISC/POLYA primary) may
    be offset by the outer-coordinate tolerance."""
    return 0 if a.primary.role == "CLIP" and b.primary.role == "CLIP" else p.tol


def _outer_in_clip(f: Fragment) -> bool:
    """The primary read's 5' end is inside the junction soft clip, so its outer coordinate
    depends on the clip length (homopolymer jitter), not on the molecule."""
    p = f.primary
    return p.role == "CLIP" and ((p.side == "LEFT" and p.strand == "+") or
                                 (p.side == "RIGHT" and p.strand == "-"))


def _mate_unreliable(f: Fragment, p: DedupParams) -> bool:
    m = f.mate
    return m is not None and m.mapped and m.mapq < p.mate_min_mapq


def _is_dup(a: Fragment, b: Fragment, p: DedupParams):
    """SPEC rule 2 (same sample only). Returns '' (not a duplicate), 'coord' (mate placed on
    both: coordinates within tol + lenient clip sequence) or 'seq' (mate unplaced on both:
    read outer within tol + lenient clip AND mate sequence)."""
    tol = p.tol
    if a.strand != b.strand:
        return ""
    ja, jb = _junction_pos(a), _junction_pos(b)
    if ja is not None and jb is not None and abs(ja - jb) > tol:
        return ""
    if not (_outer_in_clip(a) and _outer_in_clip(b)) and abs(a.outer - b.outer) > tol:
        return ""
    if _mate_unreliable(a, p) or _mate_unreliable(b, p):
        am = bm = None              # multi-mapped mate: decide by sequence
    else:
        ao, bo = a.mate_outer(), b.mate_outer()
        if ao is not None and bo is not None:
            am, bm = ao, bo
        else:
            am, bm = a.mate_coord(), b.mate_coord()
    if (am is None) != (bm is None):
        return ""
    if am is not None:
        if not (am[0] == bm[0] and am[2] == bm[2] and abs(am[1] - bm[1]) <= tol):
            return ""
        return "coord" if _seq_close(*_dup_seqs(a, b), p, _clip_shift(a, b, p)) else ""
    # no mate placement on either: same read outer AND lenient junction-clip sequence; when
    # both mate sequences are known they must agree too (shift-tolerant: an unmapped mate has
    # no coordinate, but a duplicate's mate starts within tol of the original's).
    if not _seq_close(*_dup_seqs(a, b), p, _clip_shift(a, b, p)):
        return ""
    ams, bms = a.mate_seq(), b.mate_seq()
    if ams and bms and not _seq_close(ams, bms, p, shift=tol):
        return ""
    return "seq"


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


def independent_clusters(frags: List[Fragment], params: Optional[DedupParams] = None,
                         stats: Optional[dict] = None):
    """Union-find over fragments. Returns (clusters, n_duplicates, n_cross_identical);
    `stats` (optional dict) receives the split n_dup_coord / n_dup_seq."""
    p = params or DedupParams()
    tol = p.tol
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
                if find(i) == find(j):
                    continue
                kind = _is_dup(frags[i], frags[j], p)
                if kind:
                    parent[find(j)] = find(i)
                    n_dup += 1
                    if stats is not None:
                        stats["n_dup_" + kind] = stats.get("n_dup_" + kind, 0) + 1
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
        self.n_dup_coord = self.n_dup_seq = 0
        self.member_loci = ""
        self.n_short_used = self.n_short_rejected = self.n_short_mate_inside = 0
        self.n_independent_no_short = 0
        self.short_reasons = Counter()
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
                self.polya_end, self.n_duplicates, self.n_cross, self.fail_reason, c.stop_reason,
                self.n_dup_coord, self.n_dup_seq, self.member_loci or ".",
                self.n_short_used, self.n_short_rejected, self.n_short_mate_inside,
                self.n_independent_no_short]
        return "\t".join(str(v) for v in vals) + "\n"


def evaluate_junction(insertion_id, side, rows, cfg, ref_fetch=None) -> JunctionRecord:
    """Pooled stats + consensus for one junction. `ref_fetch(contig, start, end)` (0-based,
    reference-forward) is only needed for SHORT-overhang validation (`count_short_overhang`)."""
    dedup = DedupParams.from_cfg(cfg)
    min_ind = cfg.get("min_independent_fragments", 2)
    polya_min = cfg.get("polya_min_len", 8)
    use_short = bool(cfg.get("count_short_overhang", False))
    rec = JunctionRecord(insertion_id, side)
    frags = collapse_fragments(rows)
    # SHORT-only fragments (reads crossing the junction by a few bases) never feed the
    # consensus; they are validated against it afterwards and may only ADD support.
    short = [f for f in frags if f.primary.role == "SHORT"]
    main = [f for f in frags if f.primary.role != "SHORT"]
    dstats = {}
    clusters, n_dup, n_cross = independent_clusters(main, dedup, dstats)
    # consensus input: every read, weighted 1/|cluster| so a PCR family votes once;
    # `group` = cluster so depth counts independent fragments.
    reads = []
    for gi, cl in enumerate(clusters):
        w = 1.0 / len(cl)
        for f in cl:
            for r in f.rows:
                if r.role == "SHORT":
                    continue
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
    rec.n_independent_no_short = len(clusters)
    used = []
    if use_short and short:
        has_clip = any(f.primary.role == "CLIP" for f in main)
        # validation target: the junction-read clip consensus at depth >= 1 (a single CLIP
        # fragment must suffice -- the SHORT read may supply the second fragment itself)
        anchored = [r for r in reads if r.anchored]
        vcons = indel_aware_consensus(anchored, min_depth=1, polya_min_len=polya_min).seq if anchored else ""
        for f in short:
            reason = _short_overhang_check(f.primary, side, vcons, has_clip, cfg, ref_fetch)
            if reason:
                rec.short_reasons[reason] += 1
            else:
                used.append(f)
                rec.n_short_mate_inside += int(_mate_inside(f, cfg))
        rec.n_short_used, rec.n_short_rejected = len(used), len(short) - len(used)
        if used:
            dstats = {}
            clusters, n_dup, n_cross = independent_clusters(main + used, dedup, dstats)
    kept = main + used
    # rows of unused SHORT-only fragments (and their mates) are dropped from every output;
    # with count_short_overhang off this makes SHORT rows invisible (legacy behaviour)
    used_keys = {(f.sample, f.frag) for f in used}
    drop = {(f.sample, f.frag) for f in short} - used_keys
    rec.rows = [r for r in rows if (r.sample, r.frag) not in drop
                and (r.role != "SHORT" or (r.sample, r.frag) in used_keys)]
    rec.n_mates = sum(1 for r in rec.rows if r.role == "MATE")
    rec.n_reads = len(rec.rows) - rec.n_mates
    rec.n_fragments = len(kept)
    rec.n_samples = len({f.sample for f in kept})
    rec.n_duplicates, rec.n_cross = n_dup, n_cross
    rec.n_dup_coord, rec.n_dup_seq = dstats.get("n_dup_coord", 0), dstats.get("n_dup_seq", 0)
    rec.n_independent = len(clusters)
    rec.polya_end = int(rec.consensus.polya_base is not None
                        or any(r.role == "POLYA" for r in rows)
                        or _clips_start_with_polya([r for r in rows if r.role == "CLIP"], polya_min))
    return rec


def _ref_pos_at(pos: int, cigar: str, q: int) -> int:
    """Reference coordinate (0-based) aligned to read offset q (a soft-clipped / inserted base
    maps to the next reference base)."""
    import re
    r, qi = pos, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar or ""):
        n = int(n)
        if op in "M=X":
            if q < qi + n:
                return r + (q - qi)
            r += n
            qi += n
        elif op in "IS":
            if q < qi + n:
                return r
            qi += n
        elif op in "DN":
            r += n
    return r


def _short_overhang_check(r: EvidenceRow, side: str, cons: str, has_clip: bool, cfg, ref_fetch):
    """'' if a SHORT read may count as a fragment for its junction, else the rejection reason.

    1. the junction already has >= 1 full CLIP fragment (short reads only ADD support);
    2. the overhang (read bases past clip_at, outward) has >= short_overhang_min_bases (5)
       bases matching the junction's indel-aware clip consensus;
    3. >= short_overhang_min_ref_mismatch (2) of the overhang bases differ from the reference
       at the same positions, and the overhang matches the consensus better than the reference;
    4. the overhang is not a homopolymer continuing a reference homopolymer at the junction
       (poly-A slippage next to a reference A-tract)."""
    min_b = cfg.get("short_overhang_min_bases", 5)
    min_mm = cfg.get("short_overhang_min_ref_mismatch", 2)
    if not has_clip:
        return "no_clip_fragment"
    if not cons:
        return "no_consensus"
    seq = r.seq.upper()
    at = r.clip_at
    if not (0 < at < len(seq)) or not r.mapped:
        return "bad_record"
    j = _ref_pos_at(r.pos, r.cigar, at)
    if side == "RIGHT":
        over = seq[at:]
    else:
        over = revcomp(seq[:at])
    n = min(len(over), len(cons))
    if n < min_b:
        return "overhang_too_short"
    over = over[:n]
    if ref_fetch is None:
        return "no_reference"
    if side == "RIGHT":
        ref_out = ref_fetch(r.ref, j, j + n).upper()
        ref_in = ref_fetch(r.ref, max(0, j - 6), j).upper()[::-1]
    else:
        ref_out = revcomp(ref_fetch(r.ref, max(0, j - n), j).upper())
        ref_in = revcomp(ref_fetch(r.ref, j, j + 6).upper())[::-1]
    if len(ref_out) < n:
        return "no_reference"
    c = cons[:n].upper()
    m_cons = sum(1 for a, b in zip(over, c) if a == b)
    m_ref = sum(1 for a, b in zip(over, ref_out) if a == b)
    if m_cons < min_b:
        return "consensus_mismatch"
    if n - m_ref < min_mm or m_cons <= m_ref:
        return "matches_reference"
    top = max("ACGT", key=over.count)
    if over.count(top) >= 0.8 * n:
        if ref_in[:6].count(top) >= 4 or ref_out[:6].count(top) >= 4:
            return "ref_homopolymer"
    return ""


def _mate_inside(f: Fragment, cfg) -> bool:
    """The mate of a SHORT fragment lies inside the inserted element: unmapped, on another
    contig, far from the junction, or MAPQ < short_mate_min_mapq (20)."""
    p = f.primary
    m = f.mate
    if m is None:
        return not p.mate_mapped
    if not m.mapped:
        return True
    if m.ref != p.ref or abs(m.pos - p.pos) > cfg.get("short_mate_max_dist", 1000):
        return True
    return m.mapq < cfg.get("short_mate_min_mapq", 20)


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


def apply_evidence(insertions, input_files, cfg, ref_fetch=None):
    """Evaluate every junction of every insertion from the pooled sidecars.

    Returns (kept, records, failed_names, stats) or None when no sidecar exists (the
    caller then behaves exactly as before). `records` maps insertion name -> [JunctionRecord];
    gated-out insertions keep their records (supported=0) for diagnostics.
    Mutates kept insertions' clips when `indel_aware_consensus` is on and the new
    consensus is at least as long as the legacy longest-clip choice."""
    # Pooling follows combine's own cross-sample grouping: every discovery locus that
    # intersect_insertions merged into this Insertion (Insertion.member_loci, one
    # (file, locus) per contributing record) contributes its sidecar rows, so colony A's
    # chr1:100-115 and colony B's chr1:101-115 pool when combine merged them.
    members = {i.name: _member_loci(i) for i in insertions}
    wanted = {m for ms in members.values() for m in ms}
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
    short_reasons = Counter()
    if cfg.get("count_short_overhang", False) and ref_fetch is None:
        try:
            from combine_insertions_get_sequence import get_sequence as ref_fetch
        except Exception as e:     # no genome: SHORT reads are rejected ("no_reference")
            print(f"WARNING: count_short_overhang without a reference ({e}); SHORT reads ignored")
    n_replaced = 0
    for ins in insertions:
        recs = []
        open_side = getattr(ins, "open_side", None)
        if open_side is None:            # Feature-A discordant end: TYPE_RIGHT_DISC=4 / TYPE_LEFT_DISC=5
            open_side = {4: "RIGHT", 5: "LEFT"}.get(getattr(ins, "type", None))
        for side in SIDES:
            if side == open_side:
                # one-sided locus (discovery `oneside_`, Feature-A disc end, kept poly-A
                # record): the open end has no reads by construction -- only the real side is
                # gated (it still needs >= min_independent_fragments)
                continue
            pooled = [r for m in members[ins.name] for r in (rows.get((m, side)) or rows.get((m[1], side), []))]
            rec = evaluate_junction(ins.name, side, pooled, cfg, ref_fetch)
            short_reasons.update(rec.short_reasons)
            rec.aligned = _aligned_part(ins, side)
            loci = sorted({l for _, l in members[ins.name]})
            if loci != [ins.name]:
                rec.member_loci = ",".join(loci)
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
    if cfg.get("count_short_overhang", False):
        n_used = sum(r.n_short_used for recs in records.values() for r in recs)
        n_only = sum(1 for recs in records.values() for r in recs
                     if r.n_independent >= min_ind > r.n_independent_no_short)
        print(f"SHORT overhang reads: {n_used} fragment(s) used, {sum(short_reasons.values())} rejected "
              f"({', '.join(f'{k}={v}' for k, v in sorted(short_reasons.items())) or '-'}); "
              f"{n_only} junction(s) reach >= {min_ind} independent fragments only thanks to them")
    return kept, records, failed, {"reasons": reasons, "replaced": n_replaced}


def _member_loci(ins):
    """(file basename, discovery locus id) of every record merged into `ins` -- recorded by
    Insertion.__init__/__iadd__; falls back to (file, name) for objects without it."""
    ml = list(getattr(ins, "member_loci", None) or [])
    # + (file, name) for every contributing file: one-sided loci are pooled by name in
    # intersect_insertions without __iadd__ (their files list is extended instead)
    ml += [(f, ins.name) for f in getattr(ins, "files", [])]
    return list(dict.fromkeys(ml))


def allele_forward_seq(r: EvidenceRow) -> str:
    """Sequence of an evidence record in ALLELE-forward (= reference-forward at the insertion
    site) orientation, as SPEC promises for insertions.reads.fa.gz.

    CLIP/DISC/SPAN records are aligned at the site, so their stored (BAM) sequence already is
    site-forward. A MATE or POLYA record may be unmapped (stored as sequenced) or placed on a
    paralogous element copy (stored forward relative to THAT copy); its orientation in the
    allele follows from the pair geometry: in an FR pair it is opposite to its partner, whose
    strand is the record's 0x20 bit. So: site-forward = stored if (0x10 set) == (partner
    forward), else the reverse complement."""
    seq = r.seq
    if not seq or seq == "*" or r.role not in ("MATE", "POLYA") or not (r.flag & 0x1):
        return seq
    stored_rev = bool(r.flag & 0x10)
    partner_fwd = not (r.flag & 0x20)
    return seq if stored_rev == partner_fwd else revcomp(seq)


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
                    fa.write(f">{name}|{rec.side}|{r.role}|{r.sample}|{r.frag}|{r.r12}\n"
                             f"{allele_forward_seq(r)}\n")
                    n_reads += 1
    print(f"wrote {n_rows} junction rows -> {evidence_tsv}; {n_reads} reads -> {reads_fa}")
