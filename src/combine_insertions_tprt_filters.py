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
"""Combine-level TPRT filters on the pooled junction evidence (plans/tprt_hallmarks/SPEC.md,
"Combine-level slippage reject / far-pair strictness / fuzzy merge").

All sequences here are OUTWARD from the junction (first base = the base adjacent to the
junction), exactly like `ConsensusResult.seq` and `Insertion.left_clipped/right_clipped`.

* `slippage_junction` -- a junction whose clip is (after stripping the leading copy of the
  reference homopolymer / short tandem repeat that touches the junction) either nothing
  structured (< `slippage_min_structured` bases before the next long homopolymer) or the shifted
  continuation of the outward reference. That is Illumina homopolymer / STR slippage at a
  reference tract (and STR length differences between the reads' genome and the reference):
  every slipped molecule is independent, so the >= 2 independent-fragment gate cannot remove it.
* `LibraryMatcher` -- mappy hits of a clip against the RTE consensus library (and, for the
  "other junction is informative" test, the transduction-source flanks).
* `far_pair_verdict` -- L1-mediated deletion / duplication pairs (gap < -30 or > 40) must show
  element sequence on the complex clip (sense to the consensus), a poly-A tail on the other clip,
  no conflicting element class beyond the tail, >= min_independent_fragments on both junctions
  and the same set of colonies at both breakpoints (one event)."""

import os
from typing import Dict, Iterable, List, Optional, Set, Tuple

import edlib

from indel_consensus import revcomp

_REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


# ------------------------------------------------------------------ reference repeats

def outward_reference(fetch, contig: str, junction: int, side: str, width: int = 80) -> Tuple[str, int]:
    """Reference window around a junction, oriented OUTWARD, and the junction offset in it.
    RIGHT junction B (reads aligned on [.., B), clip continues at B): ref[B-w:B+w], j = w.
    LEFT junction B (reads aligned on [B, ..), clip continues leftwards): rc(ref[B-w:B+w]), j = w.
    Returns ("", 0) when the reference is unavailable."""
    lo = max(0, junction - width)
    try:
        s = str(fetch(contig, lo, junction + width)).upper()
    except Exception:
        return "", 0
    j = junction - lo
    if side == "LEFT":
        return revcomp(s), len(s) - j
    return s, j


def repeat_at(line: str, j: int, max_period: int = 6, slack: int = 2) -> Tuple[int, int, str]:
    """Longest period-1..max_period tandem repeat of `line` that reaches within `slack` bases
    of offset j (the junction). Returns (start, end, unit) with `unit` read in the line's
    (outward) orientation at the repeat's own phase, or (j, j, "") if none spans >= 2 periods."""
    best = (j, j, "")
    n = len(line)
    for p in range(1, max_period + 1):
        for st in range(max(0, j - slack), min(n - p, j + slack) + 1):
            a = st
            while a > 0 and line[a - 1] == line[a - 1 + p]:
                a -= 1
            b = st
            while b + p < n and line[b + p] == line[b]:
                b += 1
            b += p
            if b - a >= 2 * p and b - a > best[1] - best[0] and (a <= j + slack and b >= j - slack):
                best = (a, b, line[a:a + p])
    return best


def strip_repeat(s: str, unit: str, max_mismatch_frac: float = 0.125) -> int:
    """Length of the longest prefix of `s` that continues `unit` (any phase), tolerating
    isolated sequencing errors (<= max_mismatch_frac, never two in a row, never at the end)."""
    if not unit or not s:
        return 0
    p = len(unit)
    best = 0
    for ph in range(p):
        mm = 0
        last_ok = 0
        prev_bad = False
        for i, c in enumerate(s):
            if c == unit[(i + ph) % p]:
                last_ok = i + 1
                prev_bad = False
            else:
                mm += 1
                if prev_bad or mm > max(1, int(max_mismatch_frac * (i + 1))):
                    break
                prev_bad = True
        best = max(best, last_ok)
    return best


def structured_len(s: str, polya_min: int = 8) -> int:
    """Bases of `s` before the first homopolymer run >= polya_min (the poly-A, or the junk an
    SBS read carries after a long homopolymer)."""
    i, n = 0, len(s)
    while i < n:
        k = i
        while k < n and s[k] == s[i]:
            k += 1
        if k - i >= polya_min:
            return i
        i = k
    return n


def slippage_junction(clip: str, line: str, j: int, cfg) -> str:
    """'' or the reason this junction looks like reference-tract slippage.

    clip  outward clip consensus (uppercase)
    line  outward reference window, junction at offset j (`outward_reference`)
    Rule: a reference homopolymer >= slippage_min_ref_run (8) or a period-2..6 repeat >=
    slippage_min_str_len (12) bp touches the junction; strip its continuation from the clip
    and from the outward reference; the clip is slippage when what is left has
    < slippage_min_structured (10) structured bases (`repeat_only`) or aligns (semi-global,
    edit <= max(2, 20 %)) to the start of the outward reference past the repeat
    (`repeat_shifted_reference`)."""
    if not clip or not line:
        return ""
    min_run = cfg.get("slippage_min_ref_run_combine", 8)
    min_str = cfg.get("slippage_min_str_len", 12)
    min_struct = cfg.get("slippage_min_structured", 10)
    polya_min = cfg.get("polya_min_len", 8)
    a, b, unit = repeat_at(line, j, cfg.get("slippage_max_period", 6))
    if not unit or (b - a) < (min_run if len(unit) == 1 else max(min_str, 3 * len(unit))):
        return ""
    clip = clip.upper()
    k = strip_repeat(clip, unit)
    rest = clip[k:]
    if k == 0 and rest:
        # the clip opens with its OWN long homopolymer of another base (e.g. a poly-A tail on
        # the other strand): not a continuation of this tract
        r = 0
        while r < len(rest) and rest[r] == rest[0]:
            r += 1
        if r >= polya_min and rest[0] not in unit:
            return ""
    ref_out = line[j:]
    ref_rest = ref_out[strip_repeat(ref_out, unit):]
    if structured_len(rest, polya_min) < min_struct:
        return "repeat_only"
    if len(unit) == 1:
        # post-homopolymer SBS phasing junk: the 'rest' is still dominated by the tract base
        w = rest[:40]
        if w.count(unit) >= cfg.get("slippage_junk_frac", 0.5) * len(w):
            return "repeat_junk"
    # the reads simply continue with the reference past the tract (shifted by the slipped
    # length); only the first bases are compared (sequence further out may be SBS junk)
    m = min(len(rest), 12)
    budget = max(1, m // 6)
    if len(ref_rest) >= m and edlib.align(rest[:m], ref_rest[:m + budget + 3], mode="SHW", task="distance",
                                         k=budget)["editDistance"] != -1:
        return "repeat_shifted_reference"
    return ""


# ------------------------------------------------------------------ library

def _resolve(path: str) -> str:
    if path and not os.path.isabs(path) and not os.path.exists(path):
        alt = os.path.join(_REPO, path)
        if os.path.exists(alt):
            return alt
    return path


def element_class(name: str) -> str:
    u = name.upper()
    for k in ("L1", "ALU", "SVA"):
        if u.startswith(k):
            return k
    return "FLANK"


class LibraryMatcher:
    """mappy hits of short clip sequences (>= 20 bp) against the RTE consensus library
    (`<rte_library>/consensus.fa`) and optionally the transduction-source flanks
    (`flanks_3p.fa.gz`, `flanks_5p_sva.fa.gz`). k=11 / w=3 seeds: ~90 % sensitivity at 20 bp,
    0/200 random 20-40-mers hit (E2E check)."""

    def __init__(self, library_dir: str, flanks: bool = True, min_match: int = 18):
        import mappy  # lazy: only TPRT mode needs it
        d = _resolve(library_dir)
        self.min_match = min_match
        kw = dict(k=11, w=3, min_cnt=1, min_chain_score=15, min_dp_score=20, best_n=3)
        self.elem = mappy.Aligner(os.path.join(d, "consensus.fa"), **kw)
        if not self.elem:
            raise RuntimeError(f"cannot index {os.path.join(d, 'consensus.fa')}")
        self.flank = []
        if flanks:
            for f in ("flanks_3p.fa.gz", "flanks_5p_sva.fa.gz"):
                p = os.path.join(d, f)
                if os.path.exists(p):
                    a = mappy.Aligner(p, **kw)
                    if a:
                        self.flank.append(a)
        self._cache: Dict[Tuple[str, bool], Optional[tuple]] = {}

    def hit(self, seq: str, with_flanks: bool = False):
        """Best hit (name, class, strand '+'/'-', matched bases, q_start, q_end) or None.
        Strand '+' = the sequence as given is sense to the element consensus."""
        seq = (seq or "").upper()
        if len(seq) < 20:
            return None
        key = (seq, with_flanks)
        if key in self._cache:
            return self._cache[key]
        best = None
        for al in [self.elem] + (self.flank if with_flanks else []):
            for h in al.map(seq):
                if h.mlen >= self.min_match and (best is None or h.mlen > best[3]):
                    best = (h.ctg, element_class(h.ctg) if al is self.elem else "FLANK",
                            "+" if h.strand > 0 else "-", h.mlen, h.q_st, h.q_en)
        self._cache[key] = best
        return best


# ------------------------------------------------------------------ far pairs

def leading_polyt(seq: str, n: int = 10) -> bool:
    s = (seq or "").upper()
    return len(s) >= n and s[:n].count("T") * 5 >= n * 4


def after_polyt(seq: str) -> str:
    """Clip sequence after its leading (T-dominated) poly-A tail."""
    s = (seq or "").upper()
    i = 0
    bad = 0
    while i < len(s):
        if s[i] == "T":
            bad = 0
        else:
            bad += 1
            if bad >= 2:
                i -= 1
                break
        i += 1
    return s[i:].lstrip("T")


def far_geometry(gap: int, cfg) -> bool:
    return gap < -cfg.get("far_pair_max_tsd_deletion", 30) or gap > cfg.get("far_pair_tsd_max", 40)


def colonies_consistent(a: Set[str], b: Set[str], frac: float = 0.2) -> bool:
    """Both breakpoints of one event are seen in the same colonies: they share a colony and
    differ by at most max(1, frac x |union|) colonies (a junction can be missed in a colony at
    low coverage; a TP junction paired with an unrelated breakpoint of ONE colony is not)."""
    if not (a & b):
        return False
    return len(a ^ b) <= max(1, int(frac * len(a | b)))


def far_pair_verdict(clips: Dict[str, List[str]], n_ind: Dict[str, int], colonies: Dict[str, Set[str]],
                     matcher: LibraryMatcher, cfg, inside_mates: Dict[str, List[str]] = None) -> Tuple[str, Optional[str]]:
    """('' , polya_side) if a far L1-mediated pair is credible, else (reason, polya_side).

    clips         side -> candidate outward clips, best first (pooled junction-read consensus,
                  discovery clip): the first that is long enough decides
    n_ind         side -> pooled independent fragments
    colonies      side -> colonies with a discovery breakpoint / evidence at that junction
    inside_mates  side -> mates of junction fragments that lie inside the insertion (unmapped,
                  elsewhere, MAPQ < 20), oriented like the outward clip (element sense = '+')
    (a) the complex (non-poly-A) clip hits the element library -- or, when the clip is too short
    (< 20 bp of a 5'-truncated element), >= 2 inside mates do; (b) sense to the consensus (an
    L1-mediated event is a 5'-truncated element; `far_pair_allow_antisense` admits inverted 5'
    ends); (c) a poly-A tail on the other clip (>= `far_pair_min_polya` (10) bases >= 80 % T,
    outward) and, if the sequence beyond the tail hits the library, the same element class;
    (d) >= min_independent_fragments on both junctions; (e) consistent colonies
    (`colonies_consistent`, `far_pair_colony_frac` 0.2)."""
    n_pa = cfg.get("far_pair_min_polya", 10)
    pt = {s: any(leading_polyt(c, n_pa) for c in clips.get(s, ())) for s in ("LEFT", "RIGHT")}
    pa = [s for s in pt if pt[s]]
    if len(pa) != 1:
        return ("no_polarity", None)
    pside = pa[0]
    cside = "RIGHT" if pside == "LEFT" else "LEFT"
    min_ind = cfg.get("min_independent_fragments", 2)
    if min(n_ind.get("LEFT", 0), n_ind.get("RIGHT", 0)) < min_ind:
        return ("few_fragments", pside)
    h = next((x for x in (matcher.hit(c) for c in clips.get(cside, ())) if x is not None), None)
    anti = cfg.get("far_pair_allow_antisense", False)
    if h is None:
        mh = [x for x in (matcher.hit(m) for m in (inside_mates or {}).get(cside, ())) if x is not None]
        mh = [x for x in mh if x[1] != "FLANK" and (anti or x[2] == "+")]
        if len(mh) < 2:
            return ("no_element_on_complex_clip", pside)
        h = max(mh, key=lambda x: x[3])
    if h[2] != "+" and not anti:
        return ("element_antisense", pside)
    tail = next((x for x in (matcher.hit(after_polyt(c)) for c in clips.get(pside, ())) if x is not None), None)
    if tail is not None and tail[1] != h[1]:
        return ("element_class_conflict", pside)
    if not colonies_consistent(colonies.get("LEFT", set()), colonies.get("RIGHT", set()),
                               cfg.get("far_pair_colony_frac", 0.2)):
        return ("colony_mismatch", pside)
    return ("", pside)


# ------------------------------------------------------------------ fuzzy merge helpers

def hp_compress(s: str) -> str:
    out = []
    for c in s:
        if not out or out[-1] != c:
            out.append(c)
    return "".join(out)


def polya_trimmed(seq: str, polya_min: int = 8) -> str:
    """Homopolymer-compressed clip, cut at the first A/T run >= polya_min (the poly-A tail):
    poly-A length jitter and the SBS junk a read carries after a long homopolymer never decide
    whether two colonies' clips agree. A clip that starts with its poly-A has nothing left."""
    s = str(seq).upper()
    i, n = 0, len(s)
    while i < n:
        k = i
        while k < n and s[k] == s[i]:
            k += 1
        if s[i] in "AT" and k - i >= polya_min:
            s = s[:i]
            break
        i = k
    return hp_compress(s)


def clips_agree(seqs: Iterable, min_score: float = 0.6, polya_min: int = 8, min_informative: int = 6) -> bool:
    """Clip agreement for merging records across colonies, robust to poly-A length jitter:
    clips are compared homopolymer-compressed up to their poly-A tail (`polya_trimmed`); when
    fewer than `min_informative` compressed bases remain in the shortest, they cannot disagree."""
    from sequence_checks import sequence_matching_score
    t = [polya_trimmed(s, polya_min) for s in seqs if s is not None]
    if len(t) < 2 or min(len(x) for x in t) < min_informative:
        return True
    return sequence_matching_score(t) >= min_score
