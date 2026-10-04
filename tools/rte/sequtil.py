"""Small sequence helpers shared by the tools/rte modules (no heavy dependencies)."""
from __future__ import annotations

import re

_COMP = str.maketrans("ACGTNRYKMacgtnrykm", "TGCANYRMKtgcanyrmk")


def rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def polya_runs(seq: str, base: str = "A", min_len: int = 8, max_gap: int = 1):
    """Maximal homopolymer-ish runs of `base` (case-insensitive) of length >= min_len, allowing
    single non-`base` bases between sub-runs of >= 3 (Illumina substitution errors inside a
    poly-A). Returns a list of (start, end) half-open intervals on `seq`."""
    s = seq.upper()
    base = base.upper()
    raw = [(m.start(), m.end()) for m in re.finditer(f"{base}{{3,}}", s)]
    merged = []
    for st, en in raw:
        if merged and st - merged[-1][1] <= max_gap:
            merged[-1] = (merged[-1][0], en)
        else:
            merged.append((st, en))
    return [(a, b) for a, b in merged if b - a >= min_len]


def trailing_polya_start(seq: str, min_len: int = 5) -> int:
    """Index where a trailing poly-A (allowing an occasional N/non-A) starts, or len(seq)."""
    s = seq.upper()
    i = len(s)
    bad = 0
    while i > 0 and (s[i - 1] == "A" or (s[i - 1] in "N" and bad < 2)):
        if s[i - 1] != "A":
            bad += 1
        i -= 1
    return i if len(s) - i >= min_len else len(s)


def homopolymer_at(seq: str, pos: int, direction: int) -> tuple[str, int]:
    """Longest homopolymer starting at pos and extending in `direction` (+1 right, -1 left).
    Returns (base, length)."""
    s = seq.upper()
    if not (0 <= pos < len(s)):
        return ("", 0)
    b = s[pos]
    n = 0
    i = pos
    while 0 <= i < len(s) and s[i] == b:
        n += 1
        i += direction
    return (b, n)


def hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a.upper(), b.upper()) if x != y) + abs(len(a) - len(b))


def shannon(seq: str) -> float:
    from math import log2
    s = seq.upper()
    if not s:
        return 0.0
    c = {}
    for ch in s:
        c[ch] = c.get(ch, 0) + 1
    n = len(s)
    return -sum(v / n * log2(v / n) for v in c.values())


def edlib_best(query: str, target: str, max_frac: float = 0.15, both_strands: bool = True):
    """Best infix (HW) alignment of `query` inside `target` with edlib. Returns
    (edit_distance, t_start, t_end, strand) or None. Strand -1 means rc(query) matched."""
    import edlib
    if not query or not target:
        return None
    k = max(1, int(len(query) * max_frac))
    best = None
    for strand, q in ((1, query.upper()), (-1, rc(query.upper()))):
        if strand == -1 and not both_strands:
            break
        r = edlib.align(q, target.upper(), mode="HW", task="locations", k=k)
        if r["editDistance"] < 0:
            continue
        st, en = r["locations"][0]
        cand = (r["editDistance"], st, en + 1, strand)
        if best is None or cand[0] < best[0]:
            best = cand
    return best


def edlib_path(query: str, target: str, k: int = -1):
    """Global-in-query (HW) alignment returning (ed, t_start, t_end, cigar_ops) where cigar_ops is
    a list of (length, op) with op in M/I/D semantics as minimap2 integers (0=M, 1=I, 2=D)."""
    import edlib
    r = edlib.align(query.upper(), target.upper(), mode="HW", task="path", k=k)
    if r["editDistance"] < 0:
        return None
    st, en = r["locations"][0]
    ops = []
    for n, op in re.findall(r"(\d+)([=XID])", r["cigar"] or ""):
        code = {"=": 0, "X": 0, "I": 1, "D": 2}[op]
        n = int(n)
        if ops and ops[-1][1] == code:
            ops[-1] = (ops[-1][0] + n, code)
        else:
            ops.append((n, code))
    return r["editDistance"], st, en + 1, ops
