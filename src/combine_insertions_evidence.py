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
        """Mate sequence in allele-forward orientation (a multi-mapped mate is stored relative
        to whatever paralog it landed on; two copies of one molecule can land on opposite
        strands)."""
        return allele_forward_seq(self.mate).upper() if self.mate is not None else ""

    def swapped(self):
        """The same template seen from its other read (mate as primary), or None. Two PCR
        copies of one molecule can have different junction reads (one copy's CLIP read is the
        other copy's DISC/MATE read), so the dedup also compares across that pairing."""
        if self.mate is None or self.primary is None:
            return None
        f = Fragment.__new__(Fragment)
        f.sample, f.frag, f.rows = self.sample, self.frag, self.rows
        f.primary, f.mate = self.mate, self.primary
        return f


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


def _raw_cut(seq: str, polya_min: int) -> str:
    """Raw sequence cut right after the first A/T run >= polya_min (run kept, capped at
    polya_min bases so its length jitter does not count)."""
    i, n = 0, len(seq)
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        if seq[i] in "AT" and j - i >= polya_min:
            return seq[:i + polya_min]
        i = j
    return seq


def _semi_close(a: str, b: str, budget_fn) -> bool:
    if len(a) > len(b):
        a, b = b, a
    m = len(a)
    if m == 0 or b.startswith(a):
        return True
    # semi-global: the shorter string is aligned entirely to a PREFIX of the longer one (its
    # end is free), so a read truncated earlier does not pay for the missing tail
    return edlib.align(a, b, mode="SHW", task="distance", k=budget_fn(m))["editDistance"] != -1


def _prefix_close(a: str, b: str, p: DedupParams) -> bool:
    """Lenient prefix comparison within the edit budget, two ways: homopolymer-compressed (run
    length jitter is free; a substitution may cost 2-3 RLE edits) OR raw (a substitution costs
    1; homopolymer jitter costs). Both cut after the first long A/T run."""
    if _semi_close(_hp_compress_cut(a, p.polya_min), _hp_compress_cut(b, p.polya_min), p.budget):
        return True
    return _semi_close(_raw_cut(a, p.polya_min), _raw_cut(b, p.polya_min), p.budget)


def _seq_close_from_start(a: str, b: str, p: DedupParams, shift: int = 0) -> bool:
    if not shift:            # anchored at the junction
        return _prefix_close(a, b, p)
    for k in range(shift + 1):
        if _prefix_close(a[k:], b, p) or (k and _prefix_close(a, b[k:], p)):
            return True
    return False


def _has_long_run(x: str, n: int) -> bool:
    return any(run in x for run in ("A" * n, "T" * n))


def _seq_close(a: str, b: str, p: DedupParams, shift: int = 0) -> bool:
    """Lenient sequence identity. `shift` > 0: the two reads may start up to `shift` bases
    apart (unanchored reads, e.g. an unmapped mate of a duplicate); each offset 0..shift of
    either read is tried, the comparison itself stays strict (edit budget).

    The low-quality sequence after a long poly-A lies on whichever side the SEQUENCER reached
    last, which in reference-forward / allele-forward strings can be either end. So when a long
    A/T run is present the comparison is also made from the other end (strings reversed, cut at
    the run nearest that end); the pair is close when either clean side agrees."""
    if _seq_close_from_start(a, b, p, shift):
        return True
    if _has_long_run(a, p.polya_min) or _has_long_run(b, p.polya_min):
        return _seq_close_from_start(a[::-1], b[::-1], p, shift)
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
    return allele_forward_seq(a.primary).upper(), allele_forward_seq(b.primary).upper()


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
    read outer within tol + lenient clip AND mate sequence). Also tried with b seen from its
    other read (b.swapped()): copies of one molecule may carry the evidence on different reads."""
    k = _is_dup_oriented(a, b, p)
    if k or a.strand == b.strand:
        return k
    for x, y in ((a, b.swapped()), (a.swapped(), b)):
        if x is None or y is None or not x.primary.mapped or not y.primary.mapped:
            continue
        k = _is_dup_oriented(x, y, p)
        if k:
            return k
    return ""


def _is_dup_oriented(a: Fragment, b: Fragment, p: DedupParams):
    tol = p.tol
    if a.strand != b.strand:
        return ""
    ja, jb = _junction_pos(a), _junction_pos(b)
    if ja is not None and jb is not None and abs(ja - jb) > tol:
        return ""
    # a 5' end inside the junction clip moves with homopolymer jitter as well as with the
    # fragment end: allow twice the tolerance there instead of ignoring the coordinate
    # (ignoring it merged 182 independent poly-A junction fragments in the E2E)
    otol = 2 * tol if (_outer_in_clip(a) and _outer_in_clip(b)) else tol
    if abs(a.outer - b.outer) > otol:
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
                # no pruning by the primary's outer coordinate: copies of one molecule may carry
                # the evidence on different reads (Fragment.swapped) or have a jittering outer
                # (_outer_in_clip); per junction and sample the fragment count is capped (200)
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
    # reference the overhang would read if it were reference: indel-aware (edlib infix in a
    # window), because a slipped homopolymer next to the junction shows up as a deletion in
    # an unclipped SHORT read, which shifts a per-position comparison (E2E: 156 poly-A
    # slippage junctions were rescued by exactly that before this was indel-aware)
    import re as _re
    ref_len = sum(int(k) for k, op in _re.findall(r"(\d+)([MDN=X])", r.cigar or ""))
    pad = 8
    if side == "RIGHT":
        lo, hi = j, max(j + n, r.pos + ref_len) + pad
        ref_win = ref_fetch(r.ref, lo, hi).upper()
        ref_in = ref_fetch(r.ref, max(0, j - 6), j).upper()[::-1]
    else:
        lo, hi = max(0, min(r.pos, j - n) - pad), j
        ref_win = revcomp(ref_fetch(r.ref, lo, hi).upper())
        ref_in = revcomp(ref_fetch(r.ref, j, j + 6).upper())[::-1]
    if len(ref_win) < n:
        return "no_reference"
    ref_out = ref_win[:n]
    c = cons[:n + pad].upper()
    ed_ref = edlib.align(over, ref_win, mode="HW", task="distance")["editDistance"]
    ed_cons = edlib.align(over, c, mode="SHW", task="distance")["editDistance"]
    m_cons = n - ed_cons
    if m_cons < min_b:
        return "consensus_mismatch"
    if ed_ref < min_mm or ed_cons >= ed_ref:
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


def apply_evidence(insertions, input_files, cfg, ref_fetch=None, breakpoints=None):
    """Evaluate every junction of every insertion from the pooled sidecars.

    Returns (kept, records, failed_names, stats) or None when no sidecar exists (the
    caller then behaves exactly as before). `records` maps insertion name -> [JunctionRecord];
    gated-out insertions keep their records (supported=0) for diagnostics.
    Mutates kept insertions' clips when `indel_aware_consensus` is on and the new
    consensus is at least as long as the legacy longest-clip choice.

    `breakpoints` (optional, from `discovery_breakpoints`): every per-sample discovery
    breakpoint, for the far-pair colony-consistency test (`far_pair_strict`)."""
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
    slip_on = bool(cfg.get("slippage_reject", False))
    far_on = bool(cfg.get("far_pair_strict", False))
    if (slip_on or far_on) and ref_fetch is None:
        try:
            from combine_insertions_get_sequence import get_sequence as ref_fetch
        except Exception as e:
            print(f"WARNING: slippage_reject without a reference ({e}); slippage test skipped")
            slip_on = False
    matcher = None
    if slip_on or far_on:
        from combine_insertions_tprt_filters import LibraryMatcher
        matcher = LibraryMatcher(cfg.get("rte_library") or "resources/rte_library")
    by_id = {id(i): members[i.name] for i in insertions}
    recmap = {}                     # id(insertion) -> [JunctionRecord]
    failed_objs = []
    tprt_reasons = Counter()

    def evaluate(ins, side):
        ms = by_id[id(ins)]
        allowed = getattr(ins, "member_sides", None) or {}
        pooled = [r for m in ms if side in allowed.get(m, SIDES)
                  for r in (rows.get((m, side)) or rows.get((m[1], side), []))]
        pooled = _reanchor(pooled, side, _ins_junction(ins, side))
        rec = evaluate_junction(ins.name, side, pooled, cfg, ref_fetch)
        short_reasons.update(rec.short_reasons)
        rec.aligned = _aligned_part(ins, side)
        loci = sorted({l for _, l in ms})
        if loci != [ins.name]:
            rec.member_loci = ",".join(loci)
        return rec

    for ins in insertions:
        recs = []
        open_side = _open_side(ins)
        for side in SIDES:
            if side == open_side:
                # one-sided locus (discovery `oneside_`, Feature-A disc end, kept poly-A
                # record): the open end has no reads by construction -- only the real side is
                # gated (it still needs >= min_independent_fragments)
                continue
            recs.append(evaluate(ins, side))
        recmap[id(ins)] = recs
        if far_on and open_side is None and len(recs) == 2:
            verdict = _far_pair_check(ins, recs, cfg, matcher, breakpoints, ref_fetch)
            if verdict is not None:
                reason, pside = verdict
                tprt_reasons[f"far_pair:{reason}"] += 1
                prec = [r for r in recs if r.side == pside]
                if (pside is not None and cfg.get("far_pair_split", True) and prec
                        and prec[0].n_independent >= min_ind):
                    old_name = ins.name
                    _to_one_sided(ins, pside)
                    ins.member_sides = {m: (pside,) for m in by_id[id(ins)]}
                    prec[0].insertion_id = ins.name
                    prec[0].member_loci = ",".join(sorted({l for _, l in by_id[id(ins)]} | {old_name}))
                    prec[0].fail_reason = f"split_from_far_pair:{reason}"
                    recs = prec
                    recmap[id(ins)] = recs
                    tprt_reasons["far_pair:split_to_one_sided"] += 1
                else:
                    for r in recs:
                        r.fail_reason = f"far_pair:{reason}"
                    failed_objs.append(ins)
                    reasons[f"far_pair:{reason}"] += 1
                    continue
        if slip_on:
            why = _slippage_check(ins, recs, cfg, matcher, ref_fetch)
            if why:
                tprt_reasons[why.split(":")[0]] += 1
                for r in recs:
                    r.fail_reason = why
                failed_objs.append(ins)
                reasons[why.split("(")[0]] += 1
                continue
        eligible = all(f in have for f in ins.files)
        fails = [r for r in recs if r.n_independent < min_ind]
        for r in recs:
            if eligible:
                r.supported = 0 if r in fails else 1
                if r in fails:
                    r.fail_reason = f"n_independent<{min_ind}"
        if gate and eligible and fails:
            failed_objs.append(ins)
            if all(r.n_reads == 0 for r in recs):
                reasons["no_evidence_reads"] += 1
            else:
                reasons["+".join(r.side + ("(polyA)" if r.polya_end else "") for r in fails)] += 1
            continue
        if use_cons:
            n_replaced += _replace_clips(ins, recs)
        kept.append(ins)
    # records by name (a split far pair can share its new one-sided name with another locus;
    # EvidencePool.absorb_one_sided merges such duplicates after the remap filters)
    for i in kept:
        records.setdefault(i.name, recmap[id(i)])
    kept_names = set(records)
    for i in failed_objs:
        if i.name not in kept_names:
            records.setdefault(i.name, recmap[id(i)])
            failed.add(i.name)
    if tprt_reasons:
        print("TPRT combine filters: " + ", ".join(f"{k}={v}" for k, v in sorted(tprt_reasons.items())))
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
    pool = EvidencePool(evaluate, by_id, recmap, records, cfg, use_cons)
    return kept, records, failed, {"reasons": reasons, "replaced": n_replaced, "pool": pool,
                                   "tprt": tprt_reasons}


class EvidencePool:
    """State of `apply_evidence` kept for the final one-sided merge (TPRT fuzzy merge /
    far_pair_split): run AFTER the remap filters, so a one-sided locus is folded into a
    two-sided call only when that call survived every filter (absorbing earlier lost real
    events whenever the two-sided record was later removed)."""

    def __init__(self, evaluate, by_id, recmap, records, cfg, use_cons):
        self.evaluate, self.by_id, self.recmap, self.records = evaluate, by_id, recmap, records
        self.cfg, self.use_cons = cfg, use_cons

    def absorb_one_sided(self, insertions):
        """Fold each surviving one-sided locus into a surviving insertion with the same real
        junction within `merge_tolerance_bp` (two-sided first; else another one-sided locus of the
        same side), pooling its evidence there (re-evaluated). Returns (insertions, n_absorbed)."""
        from combine_insertions_tprt_filters import clips_agree
        tol = int(self.cfg.get("merge_tolerance_bp", 0) or 0)
        if self.cfg.get("far_pair_strict", False):
            tol = max(tol, 5)
        if tol <= 0:
            return insertions, 0
        min_ind = self.cfg.get("min_independent_fragments", 2)
        alive = {id(i) for i in insertions}
        one = sorted((i for i in insertions if _open_side(i) is not None),
                     key=lambda i: (-len(i.files), i.name))
        n = 0
        for x in one:
            if id(x) not in alive:
                continue
            side = "LEFT" if _open_side(x) == "RIGHT" else "RIGHT"
            cands = [i for i in insertions if i is not x and id(i) in alive]
            t = _absorb_target(x, side, cands, tol)
            if t is None:
                continue
            xc = x.left_clipped if side == "LEFT" else x.right_clipped
            tc = t.left_clipped if side == "LEFT" else t.right_clipped
            if xc is not None and tc is not None and not clips_agree([str(tc), str(xc)]):
                continue
            old = set(self.by_id[id(t)])
            self.by_id[id(t)] = list(dict.fromkeys(self.by_id[id(t)] + self.by_id[id(x)]))
            ms = dict(getattr(t, "member_sides", None) or {})
            for m in self.by_id[id(x)]:
                if m not in old:          # x contributes only its real junction
                    ms[m] = (side,)
            t.member_sides = ms
            t.files = t.files + [f for f in x.files if f not in t.files]
            new = self.evaluate(t, side)
            recs = [new if r.side == side else r for r in self.recmap[id(t)]]
            for r in recs:
                r.supported = 1 if r.n_independent >= min_ind else 0
            self.recmap[id(t)] = recs
            self.records[t.name] = recs
            if self.use_cons:
                _replace_clips(t, [new])
            alive.discard(id(x))
            n += 1
        out = [i for i in insertions if id(i) in alive]
        names = {i.name for i in out}
        for i in insertions:
            if id(i) not in alive and i.name not in names:
                self.records.pop(i.name, None)
        return out, n


def _locus_junction(locus: str, side: str):
    """Junction coordinate of `side` from a discovery locus id (`contig:L-R`, `oneside_` /
    `polyA_` / `disc_` tokens carry a coordinate too), or None."""
    try:
        a, b = locus.rsplit(":", 1)[1].split("-")
        tok = a if side == "LEFT" else b
        return int(tok.split("_")[-1])
    except (ValueError, IndexError):
        return None


def _reanchor(rows, side, junction):
    """Rows of a fuzzy-merged member locus whose junction differs by d bp from the insertion's:
    move their `clip_at` by d so every CLIP/SHORT read's clip starts at the insertion's junction
    (a RIGHT read aligned to R+3 has the insertion's first 3 bases in its aligned part). Rows of
    the insertion's own locus are returned unchanged (exact-name pooling stays byte-identical)."""
    if junction is None:
        return rows
    out = []
    cache = {}
    for r in rows:
        if r.role in ("CLIP", "SHORT") and r.clip_at >= 0:
            jm = cache.get(r.locus)
            if jm is None:
                jm = cache[r.locus] = _locus_junction(r.locus, side)
            if jm is not None and jm != junction:
                at = r.clip_at + (junction - jm)
                if 0 <= at <= len(r.seq):
                    c = EvidenceRow.__new__(EvidenceRow)
                    for k in EvidenceRow.__slots__:
                        setattr(c, k, getattr(r, k))
                    c.clip_at = at
                    r = c
        out.append(r)
    return out


def _open_side(ins):
    """'LEFT'/'RIGHT' for a one-sided locus (the end without evidence), else None."""
    open_side = getattr(ins, "open_side", None)
    if open_side is None:            # Feature-A discordant end: TYPE_RIGHT_DISC=4 / TYPE_LEFT_DISC=5
        open_side = {4: "RIGHT", 5: "LEFT"}.get(getattr(ins, "type", None))
    return open_side


def _ins_junction(ins, side):
    return getattr(ins, "left_pos" if side == "LEFT" else "right_pos", None)


def _outward_clip(rec, ins, with_mates=False) -> str:
    """Outward clip consensus of a junction: junction reads only (combined_consensus) or with
    overlapping mates (consensus); falls back to the discovery clip when the pooled consensus
    is empty (first column below depth / disagreeing)."""
    c = rec.consensus if with_mates else (rec.combined_consensus if rec.combined_consensus.seq else rec.consensus)
    if c.seq:
        return c.seq.upper()
    old = ins.left_clipped if rec.side == "LEFT" else ins.right_clipped
    return str(old).upper() if old is not None else ""


def discovery_breakpoints(records):
    """{(contig, side): [(pos, sample)]} over every per-sample discovery record (before
    intersect): which colonies' discovery found a junction there (far-pair colony test)."""
    bp = defaultdict(list)
    for i in records:
        s = sample_name(i.files[0]) if getattr(i, "files", None) else "?"
        t = getattr(i, "type", 3)
        if t not in (2, 5) and getattr(i, "left_pos", None) is not None:      # LEFT real
            bp[(i.reference_name, "LEFT")].append((i.left_pos, s))
        if t not in (1, 4) and getattr(i, "right_pos", None) is not None:     # RIGHT real
            bp[(i.reference_name, "RIGHT")].append((i.right_pos, s))
    return bp


def _colonies(ins, rec, breakpoints, tol):
    pos = _ins_junction(ins, rec.side)
    out = {r.sample for r in rec.rows if r.role in EVIDENCE_ROLES}
    for p, s in (breakpoints or {}).get((ins.reference_name, rec.side), ()):
        if abs(p - pos) <= tol:
            out.add(s)
    return out


def _far_pair_check(ins, recs, cfg, matcher, breakpoints, ref_fetch=None):
    """None when the insertion is not a far pair or passes far_pair_verdict, else
    (reason, poly-A side or None)."""
    from combine_insertions_tprt_filters import far_geometry, far_pair_verdict
    gap = ins.right_pos - ins.left_pos
    if not far_geometry(gap, cfg):
        return None
    tol = max(int(cfg.get("merge_tolerance_bp", 0) or 0), int(cfg.get("far_pair_colony_tol", 5)))
    by = {r.side: r for r in recs}
    # junction reads only: mates of an unrelated breakpoint are reference reads and often
    # reach a nearby reference element, which would fake the element test. Candidates: the
    # pooled junction-read consensus, then the discovery clip (the consensus stops early at
    # poly-A length disagreement)
    clips = {}
    for s in by:
        old = ins.left_clipped if s == "LEFT" else ins.right_clipped
        c = [_outward_clip(by[s], ins), str(old).upper() if old is not None else ""]
        clips[s] = [x for x in dict.fromkeys(c) if x]
    n_ind = {s: by[s].n_independent for s in by}
    cols = {s: _colonies(ins, by[s], breakpoints, tol) for s in by}
    mates = {s: _inside_mates(by[s]) for s in by}
    reason, pside = far_pair_verdict(clips, n_ind, cols, matcher, cfg, mates)
    if not reason and ref_fetch is not None:
        # (f) the "poly-A tail" must not be slippage at a reference A/T tract: such a junction
        # pairs with any element-carrying breakpoint within 50 kb (the E2E's main far-pair FP)
        from combine_insertions_tprt_filters import outward_reference, slippage_junction
        line, j = outward_reference(ref_fetch, ins.reference_name, _ins_junction(ins, pside), pside)
        if slippage_junction(_outward_clip(by[pside], ins), line, j, cfg):
            reason = "polya_side_slippage"
    return None if not reason else (reason, pside)


def _inside_mates(rec, min_mapq=20, max_dist=1000):
    """Mates of a junction's evidence fragments that lie inside the insertion (unmapped, other
    contig, > max_dist away or MAPQ < min_mapq), one per fragment, oriented like the outward
    clip (RIGHT: allele-forward; LEFT: its reverse complement)."""
    ev = {(r.sample, r.frag): r for r in rec.rows if r.role in ("CLIP", "DISC")}
    out = {}
    for r in rec.rows:
        k = (r.sample, r.frag)
        if r.role != "MATE" or k not in ev or k in out or not r.seq or r.seq == "*":
            continue
        p = ev[k]
        inside = (not r.mapped or r.ref != p.ref or abs(r.pos - p.pos) > max_dist or r.mapq < min_mapq)
        if inside:
            s = allele_forward_seq(r).upper()
            out[k] = revcomp(s) if rec.side == "LEFT" else s
    return list(out.values())


def _to_one_sided(ins, real_side):
    """Reduce a two-sided insertion to a one-sided locus on `real_side` (discovery naming:
    `contig:L-oneside_L` / `contig:oneside_R-R`)."""
    c = ins.reference_name
    if real_side == "LEFT":
        ins.right_clipped = ins.right_aligned = None
        ins.right_pos = ins.left_pos
        ins.right_mates = []
        ins.type, ins.open_side = 4, "RIGHT"          # TYPE_RIGHT_DISC
        ins.name = f"{c}:{ins.left_pos}-oneside_{ins.left_pos}"
    else:
        ins.left_clipped = ins.left_aligned = None
        ins.left_pos = ins.right_pos
        ins.left_mates = []
        ins.type, ins.open_side = 5, "LEFT"           # TYPE_LEFT_DISC
        ins.name = f"{c}:oneside_{ins.right_pos}-{ins.right_pos}"


def _absorb_target(x, side, candidates, tol):
    """Insertion among `candidates` with a real `side` junction within tol of x's (full
    insertions first, then the nearest)."""
    pos = _ins_junction(x, side)
    best = None
    for t in candidates:
        if t.reference_name != x.reference_name or _open_side(t) == side:
            continue
        d = abs(_ins_junction(t, side) - pos)
        if d <= tol:
            key = (_open_side(t) is not None, d, t.name)
            if best is None or key < best[0]:
                best = (key, t)
    return best[1] if best else None


def _carries_element(rec, ins, matcher) -> bool:
    """The junction's clip hits the element / transduction-source library, or (clip too short,
    < 20 bp) >= 2 of its fragments have an inside-insertion mate that does. Mates of a slipped
    reference tract are placed reference reads and never count as inside."""
    if matcher.hit(_outward_clip(rec, ins), with_flanks=True):
        return True
    hits = 0
    for m in _inside_mates(rec):
        if matcher.hit(m, with_flanks=True) or matcher.hit(revcomp(m), with_flanks=True):
            hits += 1
            if hits >= 2:
                return True
    return False


def _slippage_check(ins, recs, cfg, matcher, ref_fetch) -> str:
    """'' or the reject reason: some junction is reference-tract slippage
    (combine_insertions_tprt_filters.slippage_junction) and no OTHER junction carries element
    or transduction-source sequence (a one-sided locus has no other junction)."""
    from combine_insertions_tprt_filters import outward_reference, slippage_junction
    slip = {}
    for r in recs:
        line, j = outward_reference(ref_fetch, ins.reference_name, _ins_junction(ins, r.side), r.side)
        # slippage only if BOTH the junction-read consensus and the consensus extended by
        # overlapping mates say so: at a real insertion in a reference tract the junction reads
        # carry the slipped poly-A (+ SBS junk), but the mates show the element beyond it
        slip[r.side] = slippage_junction(_outward_clip(r, ins), line, j, cfg)
        if slip[r.side] and r.consensus.seq and len(r.consensus.seq) > len(_outward_clip(r, ins)):
            if not slippage_junction(r.consensus.seq.upper(), line, j, cfg):
                slip[r.side] = ""
    for s, why in slip.items():
        if not why:
            continue
        # the other junction is informative only if it is not slippage itself and carries
        # element / transduction-source sequence; junction reads only (mates of a slipped
        # reference tract reach the reference Alu whose tail the tract usually is)
        others = [r for r in recs if r.side != s and not slip[r.side]]
        if not any(_carries_element(o, ins, matcher) for o in others):
            return f"slippage:{s}({why})"
    return ""


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
