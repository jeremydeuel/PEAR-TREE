"""3' transduction source lookup + the SPEC novel-source rule.

Known sources: a transduced segment (FLANK3P segment in a read layout) aligns to
`flanks_3p.fa`, whose records are 0-15 kb of genome DOWNSTREAM of each source element in element
sense. The offset of the tag's far (3') end inside the flank is the transduction endpoint (where
read-through transcription stopped, usually a downstream polyadenylation signal).

Novel sources (SPEC "Novel source rule"): a unique inserted segment that matches no library flank
is a credible novel 3' transduction source when, on the remap genome (hs1):
  * it maps uniquely (MAPQ >= 20) at >= `novel_source_min_tag_identity` (0.95) -- a diverged
    repeat-family member matching its closest genomic copy is not a placement,
  * within `novel_source_max_dist` (15 kb) DOWNSTREAM, strand-aware, of a reference L1 that is
    >= 5.5 kb long and >= 95 % identical to the L1HS consensus (or of an L1 insertion called
    elsewhere in the same cohort),
  * and the insertion carries TPRT hallmarks (checked by the caller: poly-A after the tag, TSD/EN).
Reported as TD3P_SOURCE=novel:<contig>:<start>-<end> + NOVEL_SOURCE with the source identity.

Tiers (SPEC refinement, calibrated in docs/transduction_sources.html): identity is
tools/rte_library/common.cons_identity (edlib infix both ways, robust to ragged ends) of the source
L1 to the L1HS consensus. Tier A (credible) >= `novel_source_tier_a` (0.98); tier B (reasonably
similar, L1PA2 / young L1PA3) in [`novel_source_min_identity` (0.95), 0.98) -- reported (tier in
SourceCall.tier / rte_detail novel_tier) but worth fewer score points; < 0.95 is not a source.
"""
from __future__ import annotations

import gzip
import re
from bisect import bisect_left
from collections import Counter
from dataclasses import dataclass

from .sequtil import rc


def _edlib_identity(query, target):
    """1 - edit/len(query) of the best infix (HW) alignment of query in target."""
    import edlib
    if not query or not target:
        return 0.0
    r = edlib.align(query, target, mode="HW", task="distance")
    return 1.0 - r["editDistance"] / len(query) if r["editDistance"] >= 0 else 0.0


def _load_common_cons_identity():
    import importlib.util
    import os
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rte_library", "common.py")
    try:
        spec = importlib.util.spec_from_file_location("_rte_library_common", p)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m.cons_identity
    except Exception:
        return None


_COMMON_CONS_IDENTITY = _load_common_cons_identity()


def cons_identity(seq, cons):
    """Identity of a genomic element to a class consensus, robust to ragged ends: the better of
    (consensus infix in element) and (element infix in consensus). Uses
    tools/rte_library/common.cons_identity (the definition the tiers were calibrated with) and
    falls back to an edit-distance approximation when that file is unavailable."""
    if _COMMON_CONS_IDENTITY is not None:
        return _COMMON_CONS_IDENTITY(seq, cons)
    return max(_edlib_identity(cons.upper(), seq.upper()), _edlib_identity(seq.upper(), cons.upper()))

_REGION = re.compile(r"^(.+):(\d+)-(\d+)$")

DEFAULTS = {
    "novel_source_max_dist": 15000,
    "novel_source_min_len": 5500,
    "novel_source_min_identity": 0.95,
    "novel_source_tier_a": 0.98,
    "novel_source_min_mapq": 20,
    "novel_source_min_seg": 25,
    # the tag must BE the placed genome segment, not a relative of it: a repeat-family member
    # (Alu / L1 / hAT copy) maps -- sometimes uniquely -- to the closest copy at well below
    # 100 % identity, and would make that copy's neighbourhood a fake source
    "novel_source_min_tag_identity": 0.95,
}


@dataclass
class SourceCall:
    source_id: str
    td_end: int            # offset of the transduction endpoint inside the source flank
    td_start: int          # offset of the tag start (0 = right after the element 3' end)
    n_segments: int
    novel: bool = False
    identity: float = 0.0  # novel: identity of the source L1 to the L1HS consensus
    detail: str = ""
    tier: str = ""         # novel: 'A' (>= 0.98) or 'B' (0.95-0.98); '' for known sources


def known_source(flank_segments, lib):
    """flank_segments: [Segment kind FLANK3P, element-sense frame]. Picks the source with the most
    aligned bases among sense-oriented hits (a transduction reads in element sense)."""
    per = Counter()
    ends = {}
    starts = {}
    for s in flank_segments:
        if s.strand < 0:
            continue
        sid = lib.source_for_flank(s.target)
        per[sid] += s.matches
        ends[sid] = max(ends.get(sid, -1), s.t_en)
        starts[sid] = min(starts.get(sid, 10 ** 9), s.t_st)
    if not per:
        return None
    sid = max(per, key=per.get)
    hits = [s for s in flank_segments if lib.source_for_flank(s.target) == sid and s.strand > 0]
    # strandless source (`<id>/+` / `<id>/-` flanks): the flank that matched in sense resolves
    # the source orientation -- report it
    strands = sorted({lib.flank_strand(s.target) for s in hits} - {""})
    detail = f"source_strand={','.join(strands)}" if strands else ""
    return SourceCall(sid, ends[sid], starts[sid], len(hits), detail=detail)


class L1Rmsk:
    """Young full-length L1 intervals on the remap genome from a RepeatMasker .out (or UCSC
    rmsk.txt) file: only LINE/L1 hits >= min_len are kept, so loading the whole genome's
    RepeatMasker table stays cheap."""

    def __init__(self, path, min_len=5500):
        self.by_contig = {}
        op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 15 and f[0].isdigit() and not f[0].startswith("#") and len(f) < 17:
                    # RepeatMasker .out: score div del ins contig begin end (left) strand name class/family ...
                    contig, s, e, strand, name, fam = f[4], int(f[5]) - 1, int(f[6]), f[8], f[9], f[10]
                    div = float(f[1])
                elif len(f) >= 17:
                    # UCSC rmsk.txt: bin swScore milliDiv milliDel milliIns genoName genoStart genoEnd genoLeft strand repName repClass repFamily ...
                    contig, s, e, strand, name = f[5], int(f[6]), int(f[7]), f[9], f[10]
                    fam = f"{f[11]}/{f[12]}"
                    div = int(f[2]) / 10.0
                else:
                    continue
                if not fam.startswith("LINE/L1") or e - s < min_len:
                    continue
                strand = "-" if strand in ("C", "-") else "+"
                self.by_contig.setdefault(contig, []).append((s, e, strand, name, div))
        for c in self.by_contig:
            self.by_contig[c].sort()
        self._starts = {c: [x[0] for x in v] for c, v in self.by_contig.items()}

    def upstream_of(self, contig, s, e, strand, max_dist):
        """L1s on `strand` whose 3' end lies within max_dist upstream (in that strand's sense) of
        the segment [s, e). For '+' the L1 ends before s; for '-' it starts after e."""
        lst = self.by_contig.get(contig) or self.by_contig.get(
            contig[3:] if contig.startswith("chr") else "chr" + contig) or []
        starts = self._starts.get(contig) or [x[0] for x in lst]
        out = []
        lo = bisect_left(starts, s - max_dist - 10000)
        for i in range(max(0, lo), len(lst)):
            ls, le, lstr, name, div = lst[i]
            if ls > e + max_dist:
                break
            if lstr != strand:
                continue
            if strand == "+" and le <= s + 50 and s - le <= max_dist:
                out.append((s - le, lst[i]))
            elif strand == "-" and ls >= e - 50 and ls - e <= max_dist:
                out.append((ls - e, lst[i]))
        return [x for _, x in sorted(out)]


class MappyLocator:
    """Places a sequence on the remap genome (hs1) with minimap2: path to a .mmi or FASTA."""

    def __init__(self, path):
        self.path = path
        self._al = None

    def __call__(self, seq):
        if self._al is None:
            import mappy
            self._al = mappy.Aligner(self.path, preset="sr")
        out = []
        for h in self._al.map(seq):
            if not h.is_primary:
                continue
            ctg, off = h.ctg, 0
            m = _REGION.match(ctg)          # region FASTA record `contig:start-end` (fixtures)
            if m:
                ctg, off = m.group(1), int(m.group(2))
            out.append((ctg, off + h.r_st, off + h.r_en, "+" if h.strand == 1 else "-", h.mapq,
                        h.mlen / max(1, h.blen)))
        return out


class NovelSourceFinder:
    def __init__(self, lib, cfg=None, rmsk=None, locator=None, remap_genome=None, cohort_l1=None):
        self.lib = lib
        self.cfg = dict(DEFAULTS)
        self.cfg.update(cfg or {})
        if isinstance(rmsk, str):
            rmsk = L1Rmsk(rmsk, self.cfg["novel_source_min_len"])
        self.rmsk = rmsk
        self.locator = locator
        self.genome = remap_genome
        # cohort L1 calls on the remap genome: [(contig, pos, strand)]
        self.cohort_l1 = list(cohort_l1 or [])
        self._ident_cache = {}

    def available(self):
        return self.locator is not None and (self.rmsk is not None or bool(self.cohort_l1)
                                             or bool(getattr(self.lib, "polymorphic_l1", None)))

    def _l1_identity(self, contig, s, e, strand, div):
        key = (contig, s, e)
        if key in self._ident_cache:
            return self._ident_cache[key]
        ident = None
        cons = self.lib.consensus.get("L1HS")
        if self.genome is not None and cons:
            seq = self.genome.fetch(contig, s, e)
            if seq:
                if strand == "-":
                    seq = rc(seq)
                ident = cons_identity(seq.upper(), cons)
        if ident is None:
            ident = 1.0 - div / 100.0        # RepeatMasker divergence as a proxy
        self._ident_cache[key] = ident
        return ident

    def find(self, seq):
        """seq: the transduced segment in element-sense orientation. Returns SourceCall or None."""
        c = self.cfg
        if not self.available() or len(seq) < c["novel_source_min_seg"]:
            return None
        hits = [h for h in self.locator(seq) if h[4] >= c["novel_source_min_mapq"]]
        if len(hits) != 1:
            return None
        contig, s, e, gstrand, mapq, idn = hits[0]
        if idn < c["novel_source_min_tag_identity"]:
            return None
        # source flank is downstream of the element in element sense, so the source L1 is on
        # the same genomic strand as the segment's placement and lies upstream of it.
        best = None
        if self.rmsk is not None:
            for ls, le, lstr, name, div in self.rmsk.upstream_of(contig, s, e, gstrand, c["novel_source_max_dist"]):
                ident = self._l1_identity(contig, ls, le, lstr, div)
                if ident >= c["novel_source_min_identity"]:
                    dist = (s - le) if gstrand == "+" else (ls - e)
                    best = (dist, f"{contig}:{ls}-{le}", ident, name)
                    break
        if best is None:
            for cc, pos, cstr in self.cohort_l1:
                if cc != contig or cstr != gstrand:
                    continue
                dist = (s - pos) if gstrand == "+" else (pos - e)
                if 0 <= dist <= c["novel_source_max_dist"]:
                    best = (dist, f"{contig}:{pos}", 1.0, "cohort_L1")
                    break
        poly = False
        if best is None:
            # tier B fallback: a Tubio 2014 S7 polymorphic L1 position (no strand, length or
            # sequence) in the upstream window, either orientation
            for cc, pos, pid in getattr(self.lib, "polymorphic_l1", None) or ():
                if cc != contig:
                    continue
                dist = (s - pos) if gstrand == "+" else (pos - e)
                if 0 <= dist <= c["novel_source_max_dist"] and (best is None or dist < best[0]):
                    best = (dist, f"{contig}:{pos}", 0.0, f"polymorphic_L1:{pid}")
                    poly = True
        if best is None:
            return None
        dist, src, ident, name = best
        off_end = dist + (e - s)
        tier = "B" if poly else ("A" if ident >= c["novel_source_tier_a"] else "B")
        return SourceCall(f"novel:{src}", off_end, dist, 1, True, round(ident, 4),
                          f"source={name};tag={contig}:{s}-{e}({gstrand});dist={dist};tier={tier}",
                          tier)
