"""3' transduction source lookup + the SPEC novel-source rule.

Known sources: a transduced segment (FLANK3P segment in a read layout) aligns to
`flanks_3p.fa`, whose records are 0-15 kb of genome DOWNSTREAM of each source element in element
sense. The offset of the tag's far (3') end inside the flank is the transduction endpoint (where
read-through transcription stopped, usually a downstream polyadenylation signal).

Novel sources (SPEC "Novel source rule"): a unique inserted segment that matches no library flank
is a credible novel 3' transduction source when, on the remap genome (hs1):
  * it maps uniquely (MAPQ >= 20),
  * within `novel_source_max_dist` (15 kb) DOWNSTREAM, strand-aware, of a reference L1 that is
    >= 5.5 kb long and >= 95 % identical to the L1HS consensus (or of an L1 insertion called
    elsewhere in the same cohort),
  * and the insertion carries TPRT hallmarks (checked by the caller: poly-A after the tag, TSD/EN).
Reported as TD3P_SOURCE=novel:<contig>:<start>-<end> + NOVEL_SOURCE with the source identity.
"""
from __future__ import annotations

import gzip
import re
from bisect import bisect_left
from collections import Counter
from dataclasses import dataclass

from .sequtil import rc

_REGION = re.compile(r"^(.+):(\d+)-(\d+)$")

DEFAULTS = {
    "novel_source_max_dist": 15000,
    "novel_source_min_len": 5500,
    "novel_source_min_identity": 0.95,
    "novel_source_min_mapq": 20,
    "novel_source_min_seg": 25,
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
    n = sum(1 for s in flank_segments if lib.source_for_flank(s.target) == sid and s.strand > 0)
    return SourceCall(sid, ends[sid], starts[sid], n)


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
        return self.locator is not None and (self.rmsk is not None or self.cohort_l1)

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
                import mappy
                al = mappy.Aligner(seq=cons, preset="map-ont")
                m = b = 0
                for h in al.map(seq):
                    if h.strand == 1:
                        m += h.mlen
                        b += h.blen
                if b:
                    ident = m / b
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
        if best is None:
            return None
        dist, src, ident, name = best
        off_end = dist + (e - s)
        return SourceCall(f"novel:{src}", off_end, dist, 1, True, round(ident, 4),
                          f"source={name};tag={contig}:{s}-{e}({gstrand});dist={dist}")
