"""RepeatMasker annotation loaders (UCSC rmsk.txt table for hg38, RepeatMasker .out for hs1).

Both are normalised to `Rep` records: 0-based half-open genome interval, strand +/-, consensus
begin/end (1-based, inclusive, in consensus coordinates) and the remaining consensus bases past
the alignment end (`cons_left`), so that `cons_len = cons_end + cons_left`.
"""
import collections
import os
import pickle
import re

import numpy as np

from common import open_text

Rep = collections.namedtuple(
    "Rep", "chrom start end strand name rclass family cons_begin cons_end cons_left milli_div")

YOUNG_RE = re.compile(r"^(L1HS|L1PA\d+|L1PB\d*|AluY.*|SVA_[A-F])$")


def _ucsc_row(f):
    # bin swScore milliDiv milliDel milliIns genoName genoStart genoEnd genoLeft strand
    # repName repClass repFamily repStart repEnd repLeft id
    strand = f[9]
    rs, re_, rl = int(f[13]), int(f[14]), int(f[15])
    if strand == "+":
        cb, ce, left = rs, re_, -rl
    else:
        cb, ce, left = rl, re_, -rs
    return Rep(f[5], int(f[6]), int(f[7]), strand, f[10], f[11], f[12], cb, ce, left, int(f[2]))


def _out_row(f):
    # score div del ins chrom qbeg qend (left) strand repeat class/family beg end (left) ID [*]
    strand = "+" if f[8] == "+" else "-"
    cf = f[10].split("/")
    if strand == "+":
        cb, ce, left = int(f[11]), int(f[12]), int(f[13].strip("()"))
    else:
        left, ce, cb = int(f[11].strip("()")), int(f[12]), int(f[13])
    return Rep(f[4], int(f[5]) - 1, int(f[6]), strand, f[9], cf[0], cf[1] if len(cf) > 1 else cf[0],
               cb, ce, left, int(round(float(f[1]) * 10)))


def load(path, work, tag):
    """Returns (young_records, mask) where mask = {chrom: (starts, ends)} merged intervals of
    every repeat (used to soft-mask flanks). Cached as a pickle in `work`."""
    cache = os.path.join(work, "rmsk_%s.pkl" % tag)
    if os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(path):
        with open(cache, "rb") as fh:
            return pickle.load(fh)
    young = []
    ivs = collections.defaultdict(list)
    is_out = path.endswith(".out") or path.endswith(".out.gz")
    with open_text(path) as fh:
        for line in fh:
            if is_out:
                f = line.split()
                if len(f) < 15 or not f[0].isdigit():
                    continue
                r = _out_row(f)
            else:
                r = _ucsc_row(line.rstrip("\n").split("\t"))
            ivs[r.chrom].append((r.start, r.end))
            if YOUNG_RE.match(r.name):
                young.append(r)
    mask = {}
    for c, lst in ivs.items():
        lst.sort()
        s_out, e_out = [], []
        for s, e in lst:
            if s_out and s <= e_out[-1]:
                if e > e_out[-1]:
                    e_out[-1] = e
            else:
                s_out.append(s)
                e_out.append(e)
        mask[c] = (np.array(s_out, dtype=np.int64), np.array(e_out, dtype=np.int64))
    res = (young, mask)
    with open(cache, "wb") as fh:
        pickle.dump(res, fh, protocol=4)
    return res


def softmask(seq, chrom, start, mask):
    """Lowercase the parts of seq (genome forward strand, starting at `start`) covered by repeats."""
    if chrom not in mask:
        return seq
    S, E = mask[chrom]
    end = start + len(seq)
    i = max(0, np.searchsorted(E, start, side="right"))
    s = list(seq)
    while i < len(S) and S[i] < end:
        a = max(S[i], start) - start
        b = min(E[i], end) - start
        for k in range(a, b):
            s[k] = s[k].lower()
        i += 1
    return "".join(s)


def chain_fragments(recs, names, max_gap=150, max_overlap=200):
    """Join consecutive same-strand fragments of one element (RepeatMasker splits elements around
    small indels / insertions). `names` = set of repNames treated as one family. `max_overlap` =
    how far the next fragment may restart *inside* the previous one's consensus span (SVA VNTR
    expansions realign to the same consensus stretch: use ~1000 for SVA). Returns list of
    dicts: chrom start end strand name (longest fragment's) cons_begin cons_end cons_len milli_div
    (length-weighted) n_frag."""
    by_chrom = collections.defaultdict(list)
    for r in recs:
        if r.name in names:
            by_chrom[r.chrom].append(r)
    out = []
    for c, lst in by_chrom.items():
        lst.sort(key=lambda r: r.start)
        cur = []

        def flush():
            if not cur:
                return
            best = max(cur, key=lambda r: r.end - r.start)
            L = sum(r.end - r.start for r in cur)
            out.append(dict(
                chrom=c, start=min(r.start for r in cur), end=max(r.end for r in cur),
                strand=cur[0].strand, name=best.name,
                cons_begin=min(r.cons_begin for r in cur), cons_end=max(r.cons_end for r in cur),
                cons_len=best.cons_end + best.cons_left,
                milli_div=sum(r.milli_div * (r.end - r.start) for r in cur) / L,
                n_frag=len(cur)))

        for r in lst:
            if cur:
                p = cur[-1]
                # same strand, small genomic gap, and the consensus *continues* (next fragment
                # starts where the previous ended, +-200 bp) - so two adjacent copies are not
                # merged into one over-long element
                ok = (r.strand == p.strand and r.start - p.end <= max_gap and
                      ((r.strand == "+" and r.cons_begin > p.cons_begin
                        and p.cons_end - max_overlap <= r.cons_begin <= p.cons_end + 200) or
                       (r.strand == "-" and r.cons_end < p.cons_end
                        and p.cons_begin - 200 <= r.cons_end <= p.cons_begin + max_overlap)))
                if not ok:
                    flush()
                    cur = []
            cur.append(r)
        flush()
    return out
