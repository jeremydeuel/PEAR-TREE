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
"""Indel-aware clip consensus (homopolymer-compressed star alignment).

Why: Illumina SBS gives homopolymer -- above all poly-A -- LENGTH jitter. The old
column-wise majority (consensus.find_consensus) left-aligns the clips, so every read
whose poly-A is one base longer/shorter is shifted for the rest of its length; the
first column after the poly-A is then ambiguous and the consensus stops there, losing
exactly the 3'-beyond-poly-A sequence (3' transductions, flank) that TPRT scoring needs.

How: work in run-length-encoded (RLE) space, where a homopolymer of any length is ONE
symbol. Length jitter disappears from the alignment problem entirely; only real
substitutions/indels of distinct bases remain, which edlib aligns exactly.

  1. RLE every read (base string + run lengths + mean run quality).
  2. Seed = a medoid among the longest anchored reads (a read whose own error splits a
     homopolymer would otherwise add spurious columns); align every read to it with edlib
     (junction reads: prefix mode, i.e. anchored at the junction; mates: infix mode).
     The target is padded with wildcard 'N's so reads longer than the seed register
     their overhang as extension columns instead of being forced onto the seed.
  3. Vote per column (A/C/G/T/gap) weighted by base quality; insertions vote per gap.
     Run length per column = MEDIAN of the reads' complete, cleanly aligned runs
     (a read's last run is censored -- it may end inside the homopolymer).
  4. Re-expand, re-RLE, iterate until stable (<= `iterations`).
  5. Final pass: walk outward from the junction and stop at the first column whose
     agreeing depth is < `min_depth` INDEPENDENT fragments (reads carry a `group`
     id -- the caller's independent-fragment cluster) or with genuine disagreement
     (best base not a strict weighted majority).

Pure Python + edlib (a declared dependency). Reads are short (<= ~150 bp) and per
junction there are at most a few hundred, so this is cheap.
"""

from dataclasses import dataclass, field
from typing import Hashable, List, Optional, Sequence, Tuple
import re

import edlib

_CIGAR = re.compile(r"(\d+)([=XIDM])")
_WILDCARD = [("N", b) for b in "ACGT"]
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class ClipRead:
    """One read's clipped sequence, oriented OUTWARD from the junction (first base =
    the base adjacent to the junction). `group` = independent-fragment id: reads of
    the same fragment / PCR-duplicate cluster share it and count once for depth.
    `anchored` False = a floating read (e.g. a mate) that may start anywhere."""
    seq: str
    qual: Sequence[int]
    group: Hashable
    weight: float = 1.0
    anchored: bool = True


@dataclass
class ConsensusResult:
    seq: str = ""                       # outward orientation
    depth: List[int] = field(default_factory=list)   # independent fragments per base
    score: List[int] = field(default_factory=list)   # phred-like vote margin per base
    stop_reason: str = "empty"
    polya_base: Optional[str] = None    # 'A' or 'T' (in outward orientation)
    polya_start: int = -1               # offset of the poly-A run in `seq`
    polya_len_median: Optional[int] = None
    polya_len_min: Optional[int] = None
    polya_len_max: Optional[int] = None
    beyond_polya: str = ""              # outward orientation (see _polya)
    beyond_start: int = -1              # offset of beyond_polya in `seq`
    beyond_polya_support: int = 0

    @property
    def polya_len_range(self) -> str:
        if self.polya_len_min is None:
            return ""
        return f"{self.polya_len_min}-{self.polya_len_max}"


# ----------------------------------------------------------------------------- RLE

def rle(seq: str, qual: Sequence[int]) -> Tuple[str, List[int], List[float]]:
    bases, lens, quals = [], [], []
    prev = None
    for b, q in zip(seq.upper(), qual):
        if b == prev:
            lens[-1] += 1
            quals[-1] += q
        else:
            bases.append(b)
            lens.append(1)
            quals.append(q)
            prev = b
    return "".join(bases), lens, [q / n for q, n in zip(quals, lens)]


def _qw(q: float) -> float:
    """quality -> vote weight; Q40 == 1.0, floor so a Q2 base still counts a little."""
    return min(max(q, 2.0), 40.0) / 40.0


# ----------------------------------------------------------------------- alignment

def _align(read_bases: str, cons_bases: str, anchored: bool):
    """Align RLE `read_bases` to RLE `cons_bases` (padded with wildcard Ns).
    Returns (ops, edit) where ops is a list of (op, read_idx, col):
      'M' read run read_idx sits on column col (match or mismatch)
      'I' read run read_idx is inserted before column col
      'D' column col is skipped by the read (read_idx = index of the next read run)
    Columns >= len(cons_bases) are extension (pad) columns."""
    target = cons_bases + "N" * len(read_bases)
    res = edlib.align(read_bases, target, mode="SHW" if anchored else "HW",
                      task="path", additionalEqualities=_WILDCARD)
    if res["editDistance"] < 0 or res["cigar"] is None:
        return None, None
    col = res["locations"][0][0]
    ri = 0
    ops = []
    for n, op in _CIGAR.findall(res["cigar"]):
        n = int(n)
        if op in "=XM":
            for _ in range(n):
                ops.append(("M", ri, col))
                ri += 1
                col += 1
        elif op == "I":
            for _ in range(n):
                ops.append(("I", ri, col))
                ri += 1
        else:  # D
            for _ in range(n):
                ops.append(("D", ri, col))
                col += 1
    return ops, res["editDistance"]


class _Col:
    __slots__ = ("w", "groups", "runs", "censored", "ins", "noins_w", "noins_groups")

    def __init__(self):
        self.w = {}            # base or '-' -> weight
        self.groups = {}       # base or '-' -> set(group)
        self.runs = {}         # base -> [complete run lengths]
        self.censored = {}     # base -> max censored run length
        self.ins = {}          # inserted RLE string BEFORE this column -> [w, groups, [lens lists]]
        self.noins_w = 0.0     # weight of reads spanning the gap before this col without insertion
        self.noins_groups = set()


def _mate_overlap_ok(ops, rb, rl, cons_bases, n_real, min_overlap):
    """A floating read is used only if it overlaps >= min_overlap real bases of the
    consensus, >= 85% of its real-column runs match, and the overlap is not a single
    homopolymer (so a poly-A mate cannot be placed anywhere)."""
    matched = mism = 0
    bp = 0
    distinct = set()
    for op, ri, col in ops:
        if op == "M" and col < n_real:
            if rb[ri] == cons_bases[col]:
                matched += 1
                bp += rl[ri]
                distinct.add(col)
            else:
                mism += 1
    if bp < min_overlap or len(distinct) < 6:
        return False, bp
    return matched / max(1, matched + mism) >= 0.85, bp


def _vote(cons_bases: str, reads, min_overlap: int):
    """Align all reads to cons_bases and accumulate votes. Returns columns list."""
    n_real = len(cons_bases)
    cols = {}

    def col_of(c):
        x = cols.get(c)
        if x is None:
            x = cols[c] = _Col()
        return x

    for r, (rb, rl, rq) in reads:
        if not rb:
            continue
        if r.anchored:
            ops, _ = _align(rb, cons_bases, True)
        else:
            best = None
            for orient in (0, 1):
                if orient == 0:
                    b, l, q = rb, rl, rq
                else:
                    b, l, q = rle(revcomp(r.seq), list(r.qual)[::-1])
                o, _ = _align(b, cons_bases, False)
                if o is None:
                    continue
                ok, bp = _mate_overlap_ok(o, b, l, cons_bases, n_real, min_overlap)
                if ok and (best is None or bp > best[0]):
                    best = (bp, o, b, l, q)
            if best is None:
                continue
            _, ops, rb, rl, rq = best
        if ops is None:
            continue
        last_ri = len(rb) - 1
        n = len(ops)
        # an op is a "clean match" if it is M on its own base (pad columns always match);
        # a run length is trusted only if every op within +-2 is a clean match, so a run
        # split by a substitution (AAAAGAAAA) never contributes its fragments' lengths.
        good = [op == "M" and (c >= n_real or rb[ri] == cons_bases[c]) for op, ri, c in ops]
        prev_col_aligned = None  # last column touched (M or D)
        for k, (op, ri, c) in enumerate(ops):
            if op == "M":
                b = rb[ri]
                w = r.weight * _qw(rq[ri])
                x = col_of(c)
                x.w[b] = x.w.get(b, 0.0) + w
                x.groups.setdefault(b, set()).add(r.group)
                clean = all(good[max(0, k - 2):k + 3])
                complete = ri != last_ri and (r.anchored or ri != 0)
                if complete and clean:
                    x.runs.setdefault(b, []).append(rl[ri])
                else:
                    x.censored[b] = max(x.censored.get(b, 0), rl[ri])
                # a read moving from column c-1 straight to c spans the gap without insertion
                if prev_col_aligned == c - 1:
                    x.noins_w += w
                    x.noins_groups.add(r.group)
                prev_col_aligned = c
            elif op == "D":
                qn = rq[ri] if ri <= last_ri else rq[last_ri]
                qp = rq[ri - 1] if ri > 0 else qn
                w = r.weight * _qw((qn + qp) / 2.0)
                x = col_of(c)
                x.w["-"] = x.w.get("-", 0.0) + w
                x.groups.setdefault("-", set()).add(r.group)
                if prev_col_aligned == c - 1:
                    x.noins_w += w
                    x.noins_groups.add(r.group)
                prev_col_aligned = c
            else:  # I before column c: collect the whole consecutive insertion
                if k > 0 and ops[k - 1][0] == "I":
                    continue
                j = k
                ins_b, ins_l = [], []
                wsum = 0.0
                while j < n and ops[j][0] == "I":
                    rj = ops[j][1]
                    ins_b.append(rb[rj])
                    ins_l.append(rl[rj])
                    wsum += _qw(rq[rj])
                    j += 1
                # only an insertion BETWEEN two aligned columns is evidence (a dangling
                # leading/trailing insertion is just unaligned read end)
                if k == 0 or j >= n or c == 0:
                    continue
                x = col_of(c)
                key = "".join(ins_b)
                e = x.ins.setdefault(key, [0.0, set(), []])
                e[0] += r.weight * wsum / len(ins_b)
                e[1].add(r.group)
                e[2].append(ins_l)
                prev_col_aligned = None  # the next M/D must not also count as no-insertion
    return cols


def _median(xs):
    s = sorted(xs)
    m = len(s)
    if m % 2:
        return s[m // 2]
    return int((s[m // 2 - 1] + s[m // 2]) / 2.0 + 0.5)


def _decide(cols, n_real, final: bool, min_depth: int):
    """Turn column votes into consensus columns.
    Returns list of (base, runlen, depth, margin, complete_runs) and a stop reason."""
    out = []
    reason = "end"
    if not cols:
        return out, "empty"
    last = max(cols)
    for c in range(0, last + 1):
        x = cols.get(c)
        if x is None:
            reason = "end"
            break
        # insertion before this column (not in the final pass: by then it was realigned)
        if not final and x.ins:
            key, (w, g, lens) = max(x.ins.items(), key=lambda kv: kv[1][0])
            if w > x.noins_w and len(g) >= 1 and out:
                for t, b in enumerate(key):
                    rl = _median([l[t] for l in lens])
                    out.append((b, rl, len(g), 1, []))
        if not x.w:
            reason = "end"
            break
        best = max(x.w, key=x.w.get)
        wb = x.w[best]
        others = sum(v for k, v in x.w.items() if k != best)
        depth = len(x.groups.get(best, ()))
        if final:
            if best == "-":
                if wb > others:
                    continue
                reason = "disagreement"
                break
            if wb <= others:
                reason = "disagreement"
                break
            if depth < min_depth:
                reason = "depth"
                break
        elif best == "-":
            continue
        runs = x.runs.get(best, [])
        rl = _median(runs) if runs else max(1, x.censored.get(best, 1))
        margin = max(1, int(round((wb - others) * 40)))
        out.append((best, rl, depth, margin, runs if runs else [x.censored.get(best, 1)]))
    return out, reason


def _expand(dec):
    """Merge adjacent same-base columns and expand to sequence/depth/score/run info."""
    merged = []
    for b, rl, d, m, runs in dec:
        if merged and merged[-1][0] == b:
            pb, prl, pd, pm, pruns = merged[-1]
            merged[-1] = (b, prl + rl, min(pd, d), min(pm, m), [])  # merged run: no clean stats
        else:
            merged.append((b, rl, d, m, runs))
    return merged


def indel_aware_consensus(reads: Sequence[ClipRead], min_depth: int = 2, iterations: int = 4,
                          polya_min_len: int = 8, mate_min_overlap: int = 20) -> ConsensusResult:
    anchored = [r for r in reads if r.anchored and len(r.seq)]
    if not anchored:
        return ConsensusResult()
    enc = [(r, rle(r.seq, list(r.qual))) for r in reads if len(r.seq)]
    cons = _pick_seed([e for r, e in enc if r.anchored])
    dec = None
    for _ in range(iterations):
        cols = _vote(cons, enc, mate_min_overlap)
        dec, _r = _decide(cols, len(cons), False, min_depth)
        new = "".join(_expand_bases(dec))
        if new == cons:
            break
        cons = new
    cols = _vote(cons, enc, mate_min_overlap)
    dec, reason = _decide(cols, len(cons), True, min_depth)
    merged = _expand(dec)
    res = ConsensusResult(stop_reason=reason)
    seq, depth, score = [], [], []
    col_off = []
    off = 0
    for b, rl, d, m, runs in merged:
        col_off.append(off)
        off += rl
        seq.append(b * rl)
        depth.extend([d] * rl)
        score.extend([min(93, m)] * rl)
    res.seq = "".join(seq)
    res.depth = depth
    res.score = score
    _polya(res, merged, col_off, polya_min_len)
    return res


def _pick_seed(encoded, n_candidates: int = 8, n_probe: int = 40):
    """Among the n_candidates longest anchored reads (RAW length -- RLE length would favour
    reads whose errors split homopolymers), take the one with the smallest summed prefix
    edit distance (RLE space) to up to n_probe anchored reads; probes longer than a
    candidate are cut to its length so coverage is not penalised."""
    enc = [e for e in encoded if e[0]]
    cands = []
    seen = set()
    for b, l, _ in sorted(enc, key=lambda e: sum(e[1]), reverse=True):
        if b not in seen:
            seen.add(b)
            cands.append(b)
        if len(cands) == n_candidates:
            break
    if len(cands) == 1:
        return cands[0]
    bases = [e[0] for e in enc]
    probe = bases if len(bases) <= n_probe else bases[::max(1, len(bases) // n_probe)][:n_probe]
    best = None
    for c in cands:
        target = c + "N" * max(len(p) for p in probe)
        tot = 0
        for p in probe:
            tot += edlib.align(p[:len(c)], target, mode="SHW", task="distance",
                               additionalEqualities=_WILDCARD)["editDistance"]
        key = (tot, -len(c))
        if best is None or key < best[0]:
            best = (key, c)
    return best[1]


def _expand_bases(dec):
    """bases (with run lengths) of the decided columns, re-expanded so re-RLE merges
    neighbours that became adjacent after a gap column was dropped."""
    s = "".join(b * rl for b, rl, _, _, _ in dec)
    return rle(s, [30] * len(s))[0]


def _polya(res: ConsensusResult, merged, col_off, polya_min_len):
    """Locate the poly-A (A run, or T run = poly-A on the other strand) in the outward
    consensus and the sequence 3' of it in ELEMENT sense:
      'A' run -> element sense == outward -> beyond = everything AFTER the run
      'T' run -> element sense is inward  -> beyond = everything BEFORE the run
    (a 'T' run at the junction therefore has an empty beyond: that side is reference)."""
    best = None
    for k, (b, rl, d, m, runs) in enumerate(merged):
        if b in "AT" and rl >= polya_min_len and (best is None or rl > merged[best][1]):
            best = k
    if best is None:
        return
    b, rl, d, m, runs = merged[best]
    res.polya_base = b
    res.polya_start = col_off[best]
    res.polya_len_median = rl
    res.polya_len_min = min(runs) if runs else rl
    res.polya_len_max = max(runs) if runs else rl
    if b == "A":
        start = col_off[best] + rl
        res.beyond_polya = res.seq[start:]
        res.beyond_start = start
        dep = res.depth[start:]
    else:
        res.beyond_polya = res.seq[:col_off[best]]
        res.beyond_start = 0
        dep = res.depth[:col_off[best]]
    res.beyond_polya_support = min(dep) if dep else 0
