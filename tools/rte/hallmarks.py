"""TPRT hallmarks at the insertion site: TSD, L1-endonuclease motif, poly-A, beyond-poly-A,
plus the two site-level artefact indicators (poly-A slippage context, fold-back clip).

Coordinate model (combined.txt.gz junction strings, reference-forward):
    LEFT  junction string = clip (lower) + reference (UPPER) starting at L
    RIGHT junction string = reference (UPPER) ending at R (exclusive) + clip (lower)
    allele = ref[..R) + INSERT + ref[L..)
so a target-site duplication is ref[L:R] (R > L), a target-site deletion is R < L.
An element in + orientation ends with its poly-A at the LEFT junction (3' end = LEFT, nick at L);
in - orientation the poly-T starts the RIGHT junction insert (3' end = RIGHT, nick at R).

EN motif: L1 ORF2p EN nicks 5'-TTTT/A-3' on the strand that primes reverse transcription.
For a + element that strand is the bottom strand at L: 7-mer = rc(ref[L-2:L+5]); for a - element it
is the top strand at R: 7-mer = ref[R-5:R+2]. Both are reported in the TTTTT/AA frame (Flasch 2019
consensus 5'-TTTTT/AA-3'); mismatches are counted against the PCAWG 5-mer TTTT|R (Rodriguez-Martin
2020: bins 0-1 / 2 / 3 / 4-5 mismatches out of five).
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from .sequtil import rc, hamming, homopolymer_at, edlib_best

_LOCUS = re.compile(r"^(.+):(-?\d+)-(-?\d+)$")


def split_junction(left_seq: str, right_seq: str):
    """(left_insert, left_flank, right_flank, right_insert) from the case-encoded junctions."""
    ls = left_seq or ""
    rs = right_seq or ""
    li_end = next((i for i, ch in enumerate(ls) if ch.isupper()), len(ls))
    ri_start = next((i for i, ch in enumerate(rs) if ch.islower()), len(rs))
    return ls[:li_end].upper(), ls[li_end:].upper(), rs[:ri_start].upper(), rs[ri_start:].upper()


def parse_locus(title: str):
    m = _LOCUS.match(title or "")
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def edge_run(insert: str, at_end: bool):
    """Homopolymer (A or T) at the REF-adjacent edge of an insert, tolerating one substitution
    after >= 4 bases. Returns (base, length)."""
    s = insert.upper()
    if not s:
        return ("", 0)
    if at_end:
        s = s[::-1]
    b = s[0]
    if b not in "AT":
        # allow one non-A/T base right at the junction (ligation / microhomology)
        if len(s) > 1 and s[1] in "AT":
            s = s[1:]
            b = s[0]
        else:
            return ("", 0)
    n = 0
    mm = 0
    for i, ch in enumerate(s):
        if ch == b:
            n = i + 1
        elif mm == 0 and n >= 4 and i + 1 < len(s) and s[i + 1] == b:
            mm = 1
        else:
            break
    return (b, n)


@dataclass
class PolyAInfo:
    strand: int = 0                # +1 poly-A at LEFT, -1 poly-T at RIGHT, 0 unknown
    source: str = "none"
    left_run: tuple = ("", 0)      # REF-adjacent homopolymer of the LEFT insert (its 3' edge)
    right_run: tuple = ("", 0)     # REF-adjacent homopolymer of the RIGHT insert (its 5' edge)
    both_sided: bool = False
    length: float = 0.0            # poly-A length on the strand-consistent side


def polya_info(left_insert, right_insert, ev_left=None, ev_right=None, min_len=10):
    """Strand-consistent poly-A. Evidence-TSV medians (indel-aware, pooled) win over the edge
    run measured on the junction consensus."""
    lr = edge_run(left_insert, at_end=True)
    rr = edge_run(right_insert, at_end=False)
    la = lr[1] if lr[0] == "A" else 0
    rt = rr[1] if rr[0] == "T" else 0
    if ev_left is not None and ev_left.polya_len_median and lr[0] == "A":
        la = max(la, ev_left.polya_len_median)
    if ev_right is not None and ev_right.polya_len_median and rr[0] == "T":
        rt = max(rt, ev_right.polya_len_median)
    info = PolyAInfo(left_run=lr, right_run=rr)
    # both junctions claim to be the 3' end (poly-A at LEFT *and* poly-T at RIGHT): the two
    # ends are orientation-incompatible -> slippage / ligation artefact signature. (A solitary
    # poly-A insertion is all-A on both clips, which is consistent and NOT flagged.)
    info.both_sided = la >= min_len and rt >= min_len
    if la >= 5 or rt >= 5:
        if la >= rt:
            info.strand, info.source, info.length = 1, "polyA_left", float(la)
        else:
            info.strand, info.source, info.length = -1, "polyT_right", float(rt)
    return info


@dataclass
class SiteInfo:
    contig: str | None = None
    L: int | None = None
    R: int | None = None
    located: bool = False
    tsd_seq: str = ""
    tsd_len: int | None = None     # >0 duplication, <0 deletion, 0 blunt, None unknown
    tsd_verified: bool = False
    en_motif: str = ""
    en_mismatches: int | None = None
    slippage: bool = False
    slippage_detail: str = ""


def _find_near(genome, contig, probe, lo, hi, want_end: bool, hint):
    """Position of `probe` in genome[lo:hi] (exact, else edlib <= 2 edits). Returns the
    coordinate of the probe start (or end if want_end), nearest to `hint`."""
    win = genome.fetch(contig, lo, hi)
    if not win or len(probe) < 12:
        return None
    hits = [m.start() for m in re.finditer(f"(?={re.escape(probe)})", win)]
    if hits:
        best = min(hits, key=lambda p: abs(lo + p - hint))
        return lo + best + (len(probe) if want_end else 0)
    r = edlib_best(probe, win, max_frac=2 / len(probe), both_strands=False)
    if r is None:
        return None
    _, ts, te, _ = r
    return lo + (te if want_end else ts)


def locate_site(title, left_flank, right_flank, genome=None, probe_len=30, slack=1500):
    si = SiteInfo()
    loc = parse_locus(title)
    if loc is not None:
        si.contig = loc[0]
    if genome is None or loc is None:
        return si
    contig, a, b = loc
    lo, hi = min(a, b) - slack, max(a, b) + slack
    if right_flank:
        R = _find_near(genome, contig, right_flank[-probe_len:], lo, hi, True, b)
        si.R = R
    if left_flank:
        L = _find_near(genome, contig, left_flank[:probe_len], lo, hi, False, a)
        si.L = L
    si.located = si.L is not None and si.R is not None
    return si


def tsd_from_flanks(left_flank, right_flank, max_len=80, min_len=4, max_mm=1):
    """Genome-free TSD: the duplicated sequence ends the RIGHT flank and starts the LEFT flank.
    Largest t with <= max_mm mismatches."""
    best = None
    for t in range(min(max_len, len(left_flank), len(right_flank)), min_len - 1, -1):
        if hamming(right_flank[-t:], left_flank[:t]) <= max_mm:
            best = t
            break
    if best is None:
        return "", None, False
    return left_flank[:best], best, True


def target_site(si: SiteInfo, left_flank, right_flank, genome=None):
    """Fill tsd_seq / tsd_len / tsd_verified on `si`."""
    if si.located and genome is not None:
        t = si.R - si.L
        si.tsd_len = t
        if t > 0:
            si.tsd_seq = genome.fetch(si.contig, si.L, si.R)
            # verify against the read-derived flanks (sample SNVs/seq errors: <= 1 mismatch)
            if len(right_flank) >= t and len(left_flank) >= t:
                si.tsd_verified = hamming(right_flank[-t:], left_flank[:t]) <= 1
            else:
                si.tsd_verified = True
        else:
            si.tsd_seq = ""
            si.tsd_verified = True
        return si
    seq, t, ok = tsd_from_flanks(left_flank, right_flank)
    if t is not None:
        si.tsd_seq, si.tsd_len, si.tsd_verified = seq, t, ok
    return si


def en_motif(si: SiteInfo, strand: int, genome=None):
    """7-mer at the nick in the TTTTT/AA frame + mismatches to TTTT|R (0..5)."""
    if genome is None or not si.located or strand == 0:
        return si
    if strand > 0:
        m = rc(genome.fetch(si.contig, si.L - 2, si.L + 5))
    else:
        m = genome.fetch(si.contig, si.R - 5, si.R + 2)
    if len(m) != 7:
        return si
    si.en_motif = f"{m[:5]}/{m[5:]}"
    si.en_mismatches = sum(1 for c in m[1:5] if c != "T") + (0 if m[5] in "AG" else 1)
    return si


def en_bin(mm):
    if mm is None:
        return "."
    if mm <= 1:
        return "0-1"
    if mm == 2:
        return "2"
    if mm == 3:
        return "3"
    return "4-5"


def slippage_context(si: SiteInfo, strand: int, genome=None, min_run=10):
    """Reference homopolymer (same base as the tail) touching the 3' breakpoint: a soft clip of
    a reference poly-A whose length slipped, not a tail (PEAR-TREE memory: poly-A slippage at
    reference Alu tails is the dominant discovery FP)."""
    if genome is None or not si.located or strand == 0:
        return si
    if strand > 0:
        p, base = si.L, "A"
    else:
        p, base = si.R, "T"
    ctx = genome.fetch(si.contig, p - 30, p + 30)
    if len(ctx) < 60:
        return si
    c = 30
    for b in (base, "A" if base == "T" else "T"):
        right = homopolymer_at(ctx, c, +1)
        left = homopolymer_at(ctx, c - 1, -1)
        for (hb, n), where in ((right, "after"), (left, "before")):
            if hb == b and n >= min_run:
                si.slippage = True
                si.slippage_detail = f"ref {b}{n} {where} 3' breakpoint"
                return si
    return si


def foldback(left_insert, left_flank, right_insert, right_flank, min_len=15, max_mm=2):
    """Clip identical to the reverse complement of the reference adjacent to the junction (a
    fold-back / hairpin chimera). RIGHT: insert starts with rc(last k ref bases); LEFT: insert
    ends with rc(first k ref bases). Returns the longest such k (0 if < min_len)."""
    best = 0
    if right_insert and right_flank:
        r = rc(right_flank)
        k = _longest_prefix_match(right_insert, r, max_mm)
        best = max(best, k)
    if left_insert and left_flank:
        r = rc(left_flank)
        k = _longest_prefix_match(left_insert[::-1], r[::-1], max_mm)
        best = max(best, k)
    return best if best >= min_len else 0


def _longest_prefix_match(a, b, max_mm):
    mm = 0
    n = 0
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            mm += 1
            if mm > max_mm:
                break
        n = i + 1
    # do not count trailing mismatches
    while n > 0 and a[n - 1] != b[n - 1]:
        n -= 1
    return n
