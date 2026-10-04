"""Per-insertion assembly of the COVERED part of the inserted element.

Every sequence we have for an insertion -- the two junction clip consensuses (evidence TSV, or
the combined.txt.gz junction strings) and every pooled read / mate (reads.fa) -- is decomposed
into labelled segments (a "read layout"):

    REF      aligns to the insertion site on the discovery genome (or to the junction flanks)
    LOCAL    aligns to the site neighbourhood but NOT contiguously with the junction (templated)
    ELEMENT  aligns to a class consensus (L1HS, ALU_Y..., SVA_E...), with consensus coords+strand
    POLYA    A/T homopolymer run (poly-A tail; jitter-tolerant)
    FLANK3P  aligns to a source 3' flank (transduction), with source id + offset in the flank
    FLANK5P  aligns to an SVA source 5' flank
    UNKNOWN  unexplained (candidate transduction / templated / pre-mRNA / novel source)

mappy (minimap2) does the bulk; chimeric reads yield several hits (both strands), each becomes a
segment. Short gaps mappy cannot seed (15-40 bp) are rescued with edlib. Overlaps are resolved
greedily by aligned matches, REF first on ties (a site flank inside an old repeat must stay REF).

The covered element is then piled up on the chosen consensus (indel-aware: deletions vote '-',
insertions vote as strings after the preceding consensus base) and the majority sequence of each
covered interval is the covered-part consensus. Its identity to the intact elements gives the
nearest (active) element.

Orientation conventions: reads arrive reference-forward. Once the element strand is known every
layout is flipped into ELEMENT-SENSE orientation, where the 5' junction reads REF|element and the
3' junction reads element|...|POLYA|REF. All structure logic works in that frame.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field

from .sequtil import rc, polya_runs, edlib_best, edlib_path

DEFAULTS = {
    "min_segment_len": 20,       # shortest mappy segment kept
    "min_element_identity": 0.80,
    "min_ref_identity": 0.90,
    "min_flank_identity": 0.90,
    "min_gap_rescue": 15,        # unclaimed gaps >= this go to the edlib rescue
    "rescue_max_edit_frac": 0.15,
    "polya_min_run": 8,          # poly-A segment (inside reads; the hallmark threshold is in score)
    "local_window": 600,         # bp either side of the site used as the REF/LOCAL reference
    "min_element_bp": 30,        # covered element bp needed to call an element class
}


@dataclass
class Segment:
    q_st: int
    q_en: int
    kind: str
    target: str = ""
    t_st: int = -1
    t_en: int = -1
    strand: int = 0
    identity: float = 0.0
    matches: int = 0
    cigar: list | None = None

    @property
    def qlen(self):
        return self.q_en - self.q_st

    def flipped(self, n: int) -> "Segment":
        return Segment(n - self.q_en, n - self.q_st, self.kind, self.target, self.t_st, self.t_en,
                       -self.strand if self.strand else 0, self.identity, self.matches, self.cigar)


@dataclass
class ReadLayout:
    name: str
    side: str           # LEFT / RIGHT (as delivered), '' for mates without side
    role: str           # JUNCTION (clip consensus), CLIP, POLYA, MATE, DISC, SPAN
    frag_key: tuple
    seq: str
    segments: list = field(default_factory=list)
    sense: bool = False  # True once flipped into element-sense orientation

    def to_sense(self, strand: int) -> "ReadLayout":
        if strand >= 0:
            return ReadLayout(self.name, self.side, self.role, self.frag_key, self.seq,
                              list(self.segments), True)
        n = len(self.seq)
        segs = [s.flipped(n) for s in reversed(self.segments)]
        return ReadLayout(self.name, self.side, self.role, self.frag_key, rc(self.seq), segs, True)

    def kinds(self):
        return [s.kind for s in self.segments]

    def ref_left(self):
        return bool(self.segments) and self.segments[0].kind == "REF"

    def ref_right(self):
        return bool(self.segments) and self.segments[-1].kind == "REF"


@dataclass
class SiteContext:
    """What we know about the insertion site on the discovery genome."""
    title: str
    contig: str | None = None
    left_bp: int | None = None      # L: first reference base of the LEFT junction record
    right_bp: int | None = None     # R: end (exclusive) of the RIGHT junction reference
    window_start: int = 0
    window_seq: str = ""            # discovery-genome window around the site (may be '')
    left_flank: str = ""            # upper-case reference part of the LEFT junction string
    right_flank: str = ""           # upper-case reference part of the RIGHT junction string
    _aligner: object = None

    def local_reference(self):
        """(sequence, offset) used for REF/LOCAL labelling: the genome window when available,
        else the two junction flanks joined by Ns (offset None = no genome coordinates)."""
        if self.window_seq:
            return self.window_seq, self.window_start
        return self.right_flank + "N" * 30 + self.left_flank, None

    def aligner(self):
        if self._aligner is None:
            import mappy
            seq, _ = self.local_reference()
            self._aligner = mappy.Aligner(seq=seq, k=11, w=3, min_chain_score=18,
                                          min_dp_score=25, best_n=6) if len(seq) >= 30 else False
        return self._aligner or None


def _hit_identity(h):
    return h.mlen / h.blen if h.blen else 0.0


def _cut_hit_at(h, cons_end):
    """q-interval of a consensus hit truncated at the consensus element end (drops the part
    that aligned to the consensus' own poly-A tail). Returns (q_st, q_en, r_en)."""
    if h.r_en <= cons_end:
        return h.q_st, h.q_en, h.r_en
    cut = h.r_en - cons_end
    if h.strand == 1:
        return h.q_st, max(h.q_st, h.q_en - cut), cons_end
    return min(h.q_en, h.q_st + cut), h.q_en, cons_end


class Assembler:
    def __init__(self, lib, cfg: dict | None = None):
        self.lib = lib
        self.cfg = dict(DEFAULTS)
        self.cfg.update(cfg or {})
        self._edlib_targets = None

    # ------------------------------------------------------------------ layout of one read
    def layout(self, seq: str, ctx: SiteContext, name="", side="", role="", frag_key=(),
               ref_interval=None, rescue_flanks=False) -> ReadLayout:
        """Decompose one sequence into labelled segments. `ref_interval` (q_st, q_en) marks a
        known reference part (junction clip consensus: upper-case bases)."""
        c = self.cfg
        seq = seq.upper()
        n = len(seq)
        cands = []   # (priority score, Segment)
        if ref_interval is not None and ref_interval[1] > ref_interval[0]:
            a, b = ref_interval
            cands.append((3e6, Segment(a, b, "REF", "site", identity=1.0, matches=b - a)))
        else:
            al = ctx.aligner()
            if al is not None:
                for h in al.map(seq):
                    idn = _hit_identity(h)
                    if idn < c["min_ref_identity"] or h.q_en - h.q_st < c["min_segment_len"]:
                        continue
                    cands.append((2e6 + h.mlen + 2, Segment(h.q_st, h.q_en, "REF", "site", h.r_st, h.r_en,
                                                      h.strand, idn, h.mlen)))
        for h in self.lib.aligner("consensus").map(seq):
            idn = _hit_identity(h)
            if idn < c["min_element_identity"]:
                continue
            cons_end = self.lib.cons_end.get(h.ctg, h.ctg_len)
            if h.r_st >= cons_end:
                continue
            qs, qe, re_ = _cut_hit_at(h, cons_end)
            if qe - qs < c["min_segment_len"]:
                continue
            cands.append((2e6 + h.mlen * idn, Segment(qs, qe, "ELEMENT", h.ctg, h.r_st, re_, h.strand, idn,
                                                int(h.mlen * (qe - qs) / max(1, h.q_en - h.q_st)),
                                                h.cigar)))
        for kind, key in (("FLANK3P", "flanks3"), ("FLANK5P", "flanks5")):
            al = self.lib.aligner(key)
            if al is None:
                continue
            for h in al.map(seq):
                idn = _hit_identity(h)
                if idn < c["min_flank_identity"] or h.q_en - h.q_st < c["min_segment_len"]:
                    continue
                # a flank hit made of the source's own poly-A remnant is not informative
                piece = seq[h.q_st:h.q_en]
                if max(piece.count("A"), piece.count("T")) > 0.8 * len(piece):
                    continue
                cands.append((h.mlen * idn - 1, Segment(h.q_st, h.q_en, kind, h.ctg, h.r_st, h.r_en,
                                                        h.strand, idn, h.mlen)))
        # poly-A runs compete as candidates: they beat a flank/consensus hit made mostly of the
        # same homopolymer, but lose against a longer REF/ELEMENT alignment covering them
        for base in ("A", "T"):
            for a, b in polya_runs(seq, base, c["polya_min_run"]):
                cands.append((1e6 + b - a, Segment(a, b, "POLYA", base, identity=1.0, matches=b - a,
                                                   strand=1 if base == "A" else -1)))
        accepted = []
        for s in self._resolve(cands, n):
            if s.kind == "POLYA":
                if s.qlen >= c["polya_min_run"]:
                    accepted.append(s)
                continue
            piece = seq[s.q_st:s.q_en]
            if s.kind != "REF" and (s.qlen < c["min_segment_len"] or
                                    max(piece.count("A"), piece.count("T")) > 0.8 * len(piece)):
                continue      # trimmed to a sliver, or a homopolymer masquerading as a hit
            accepted.append(s)
        accepted.sort(key=lambda s: s.q_st)
        # edlib rescue of the remaining gaps
        for ga, gb in self._gaps(accepted, n, c["min_gap_rescue"]):
            seg = self._rescue(seq[ga:gb], ga, ctx, rescue_flanks)
            accepted.append(seg)
        accepted.sort(key=lambda s: s.q_st)
        # REF segments that are not contiguous with a read end are LOCAL (templated) candidates
        lay = ReadLayout(name, side, role, frag_key, seq, accepted)
        self._mark_local(lay, ctx)
        return lay

    @staticmethod
    def _gaps(segs, n, min_len):
        out = []
        pos = 0
        for s in sorted(segs, key=lambda x: x.q_st):
            if s.q_st - pos >= min_len:
                out.append((pos, s.q_st))
            pos = max(pos, s.q_en)
        if n - pos >= min_len:
            out.append((pos, n))
        return out

    @staticmethod
    def _resolve(cands, n):
        """Greedy by score; a lower-scoring candidate overlapping accepted ones is cut down to
        its largest free part (kept if still >= 15 bp). Chimeric split hits typically overlap
        by a few bp at the switch point, which this trims away."""
        acc = []
        for _, s in sorted(cands, key=lambda x: -x[0]):
            free = [(s.q_st, s.q_en)]
            for a in acc:
                nxt = []
                for fa, fb in free:
                    if a.q_en <= fa or a.q_st >= fb:
                        nxt.append((fa, fb))
                        continue
                    if fa < a.q_st:
                        nxt.append((fa, a.q_st))
                    if a.q_en < fb:
                        nxt.append((a.q_en, fb))
                free = nxt
            if not free:
                continue
            fa, fb = max(free, key=lambda x: x[1] - x[0])
            if fb - fa < (8 if s.kind == "POLYA" else 15):
                continue
            if (fa, fb) != (s.q_st, s.q_en):
                frac = (fb - fa) / max(1, s.qlen)
                t_st, t_en = s.t_st, s.t_en
                if s.t_st >= 0:     # proportional; exact coords are recomputed in the pileup
                    if s.strand >= 0:
                        t_st, t_en = s.t_st + (fa - s.q_st), s.t_en - (s.q_en - fb)
                    else:
                        t_st, t_en = s.t_st + (s.q_en - fb), s.t_en - (fa - s.q_st)
                s = Segment(fa, fb, s.kind, s.target, t_st, t_en, s.strand, s.identity,
                            int(s.matches * frac), None)
            acc.append(s)
        return acc

    def _edlib_targets_list(self):
        if self._edlib_targets is None:
            self._edlib_targets = [(n, s[:self.lib.cons_end.get(n, len(s))])
                                   for n, s in self.lib.consensus.items()]
        return self._edlib_targets

    def _rescue(self, piece, offset, ctx, rescue_flanks):
        c = self.cfg
        L = len(piece)
        mf = c["rescue_max_edit_frac"]
        best = None   # (edit_frac, Segment)
        # local site (both strands) first: templated / REF extension
        ref, off = ctx.local_reference()
        r = edlib_best(piece, ref, max_frac=0.10)
        if r is not None:
            ed, ts, te, strand = r
            best = (ed / L - 0.02, Segment(offset, offset + L, "REF", "site", ts, te, strand,
                                          1 - ed / L, L - ed))
        for name, cons in self._edlib_targets_list():
            r = edlib_best(piece, cons, max_frac=mf)
            if r is None:
                continue
            ed, ts, te, strand = r
            if best is None or ed / L < best[0]:
                best = (ed / L, Segment(offset, offset + L, "ELEMENT", name, ts, te, strand,
                                        1 - ed / L, L - ed))
        if rescue_flanks:
            for kind, fl in (("FLANK3P", self.lib.flanks3), ("FLANK5P", self.lib.flanks5)):
                for name, s in fl.items():
                    r = edlib_best(piece, s, max_frac=0.10)
                    if r is None:
                        continue
                    ed, ts, te, strand = r
                    if best is None or ed / L < best[0]:
                        best = (ed / L, Segment(offset, offset + L, kind, name, ts, te, strand,
                                                1 - ed / L, L - ed))
        if best is not None:
            return best[1]
        return Segment(offset, offset + L, "UNKNOWN")

    @staticmethod
    def _mark_local(lay, ctx):
        """REF segments that are internal to the read (not at either end) cannot be the
        junction flank: they are a site-derived template embedded in the insert -> LOCAL."""
        segs = lay.segments
        for i, s in enumerate(segs):
            if s.kind != "REF":
                continue
            at_end = (i == 0 and s.q_st <= 3) or (i == len(segs) - 1 and s.q_en >= len(lay.seq) - 3)
            if not at_end:
                s.kind = "LOCAL"

    # ------------------------------------------------------------------ whole insertion
    def assemble(self, ctx: SiteContext, junction_seqs: dict, reads: list, strand_hint=None):
        """junction_seqs: side -> (seq, ref_interval). reads: [EvidenceRead].
        Returns an AssemblyResult."""
        layouts = []
        for side, (seq, ref_iv) in junction_seqs.items():
            if seq:
                layouts.append(self.layout(seq, ctx, f"junction_{side}", side, "JUNCTION",
                                           ("consensus", side), ref_iv, rescue_flanks=True))
        for r in reads:
            layouts.append(self.layout(r.seq, ctx, f"{r.side}|{r.role}|{r.sample}|{r.frag}|{r.r12}",
                                       r.side, r.role, (r.sample, r.frag)))
        return AssemblyResult.build(self, ctx, layouts, strand_hint)


@dataclass
class AssemblyResult:
    strand: int = 0                    # element strand on the reference (+1/-1, 0 unknown)
    strand_source: str = ""
    element_class: str = ""            # L1 / ALU / SVA / '' (none)
    consensus: str = ""                # best consensus name
    element_bp: int = 0                # aligned element bases (all layouts)
    layouts: list = field(default_factory=list)          # element-sense layouts
    raw_layouts: list = field(default_factory=list)      # reference-forward layouts
    covered: list = field(default_factory=list)          # merged [(s, e)] on the consensus
    covered_seqs: list = field(default_factory=list)     # majority sequence per interval
    segments_on_cons: list = field(default_factory=list)  # (s, e, sense:bool, layout_idx)
    consensus_identity: float = 0.0
    nearest_intact: str = "."
    nearest_intact_identity: float = 0.0
    nearest_active: str = "."
    element_identity: float = 0.0
    class_bp: dict = field(default_factory=dict)

    @property
    def covered_5p(self):
        return min((s for s, _ in self.covered), default=-1)

    @property
    def covered_3p(self):
        return max((e for _, e in self.covered), default=-1)

    @classmethod
    def build(cls, asm: Assembler, ctx, raw_layouts, strand_hint=None):
        lib = asm.lib
        res = cls(raw_layouts=raw_layouts)
        # --- class / consensus by aligned matches
        per_cons = Counter()
        per_class = Counter()
        for lay in raw_layouts:
            w = 2 if lay.role == "JUNCTION" else 1
            for s in lay.segments:
                if s.kind == "ELEMENT":
                    per_cons[s.target] += s.matches * w
                    per_class[lib.cons_class.get(s.target, "OTHER")] += s.qlen
        res.class_bp = dict(per_class)
        res.element_bp = sum(per_class.values())
        if per_cons and res.element_bp >= asm.cfg["min_element_bp"]:
            best_class = max(per_class, key=per_class.get)
            res.consensus = max((n for n in per_cons if lib.cons_class.get(n) == best_class),
                                key=per_cons.get)
            res.element_class = best_class
        # --- strand
        if strand_hint:
            res.strand, res.strand_source = strand_hint
        else:
            res.strand, res.strand_source = cls._strand_from_segments(raw_layouts, res.consensus, lib)
        st = res.strand if res.strand else 1
        res.layouts = [l.to_sense(st) for l in raw_layouts]
        if res.consensus:
            res._pileup(asm)
            res._nearest(asm)
        return res

    @staticmethod
    def _strand_from_segments(layouts, consensus, lib):
        """Fallback element strand: majority (bp) strand of element segments adjacent to a
        junction REF, else of all element segments of the chosen class."""
        cls_ = lib.cons_class.get(consensus) if consensus else None
        adj = Counter()
        allc = Counter()
        for lay in layouts:
            segs = lay.segments
            for i, s in enumerate(segs):
                if s.kind != "ELEMENT" or (cls_ and lib.cons_class.get(s.target) != cls_):
                    continue
                allc[s.strand] += s.qlen
                if (i > 0 and segs[i - 1].kind == "REF") or (i + 1 < len(segs) and segs[i + 1].kind == "REF"):
                    adj[s.strand] += s.qlen
        for cnt, src in ((adj, "junction_segments"), (allc, "segments")):
            if cnt:
                s = max(cnt, key=cnt.get)
                return s, src
        return 0, "none"

    def _pileup(self, asm):
        lib = asm.lib
        cons = lib.consensus[self.consensus]
        cend = lib.cons_end.get(self.consensus, len(cons))
        cls_ = lib.cons_class.get(self.consensus)
        votes = defaultdict(Counter)        # pos -> base / '-'
        ins = defaultdict(Counter)          # pos -> inserted string after pos ('' = none)
        ident_m = ident_b = 0
        for li, lay in enumerate(self.layouts):
            for s in lay.segments:
                if s.kind != "ELEMENT" or lib.cons_class.get(s.target) != cls_:
                    continue
                piece = lay.seq[s.q_st:s.q_en]
                strand_c = s.strand
                if strand_c < 0:
                    piece = rc(piece)
                # re-place the piece on the chosen consensus (exact ops for the pileup)
                lo = max(0, s.t_st - 40) if s.target == self.consensus and s.t_st >= 0 else 0
                hi = min(cend, s.t_en + 40) if s.target == self.consensus and s.t_en >= 0 else cend
                k = max(1, int(len(piece) * 0.25))
                r = edlib_path(piece, cons[lo:hi], k=k)
                if r is None and (lo, hi) != (0, cend):
                    lo, hi = 0, cend
                    r = edlib_path(piece, cons[lo:hi], k=k)
                if r is None:
                    continue
                ed, ts, te, ops = r
                t = lo + ts
                q = 0
                # exact coordinates on the chosen consensus from here on
                s.target, s.t_st, s.t_en = self.consensus, lo + ts, lo + te
                s.identity = 1 - ed / max(1, len(piece))
                self.segments_on_cons.append((lo + ts, lo + te, strand_c > 0, li))
                ident_m += len(piece) - ed
                ident_b += max(len(piece), te - ts)
                for ln, op in ops:
                    if op == 0:
                        for kk in range(ln):
                            votes[t + kk][piece[q + kk]] += 1
                            ins[t + kk][""] += 1
                        t += ln; q += ln
                    elif op == 1:
                        if t - 1 >= 0:
                            ins[t - 1][""] -= 1
                            ins[t - 1][piece[q:q + ln]] += 1
                        q += ln
                    elif op == 2:
                        for kk in range(ln):
                            votes[t + kk]["-"] += 1
                        t += ln
        self.consensus_identity = ident_m / ident_b if ident_b else 0.0
        # merged covered intervals
        ivs = sorted((a, b) for a, b, _, _ in self.segments_on_cons)
        merged = []
        for a, b in ivs:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        self.covered = [(a, b) for a, b in merged]
        seqs = []
        for a, b in self.covered:
            out = []
            for p in range(a, b):
                v = votes.get(p)
                if v:
                    base = v.most_common(1)[0][0]
                    if base != "-":
                        out.append(base)
                iv = ins.get(p)
                if iv:
                    best, cnt = max(iv.items(), key=lambda x: x[1])
                    if best:
                        out.append(best)
            seqs.append("".join(out))
        self.covered_seqs = seqs

    def _nearest(self, asm):
        lib = asm.lib
        al = lib.aligner("intact")
        if al is None:
            return
        acc = defaultdict(lambda: [0, 0])
        for s in self.covered_seqs:
            if len(s) < 30:
                continue
            seen = set()
            for h in al.map(s):
                if h.ctg in seen:
                    continue
                seen.add(h.ctg)
                acc[h.ctg][0] += h.mlen
                acc[h.ctg][1] += h.blen
        if not acc:
            return
        def key(k):
            m, b = acc[k]
            return (m / b if b else 0, b)
        best = max(acc, key=key)
        self.nearest_intact = best
        self.nearest_intact_identity = round(key(best)[0], 4)
        act = [k for k in acc if lib.is_active_element(k)]
        if act:
            ba = max(act, key=key)
            self.nearest_active = ba
            self.element_identity = round(key(ba)[0], 4)
