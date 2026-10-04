"""Shared helpers for the RTE reference library (tools/rte_library/).

Kept dependency-light: py2bit, pyliftover, edlib (all in requirements.txt).
Coordinates are 0-based half-open internally; TSV outputs say which convention they use.
"""
import gzip
import os

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_RC)[::-1]


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path):
    """Ordered dict-like list of (name, description, seq)."""
    out = []
    name = desc = None
    buf = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, desc, "".join(buf)))
                head = line[1:].split(None, 1)
                name = head[0]
                desc = head[1] if len(head) > 1 else ""
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, desc, "".join(buf)))
    return out


def write_fasta(fh, name, seq, desc="", width=80):
    fh.write(">%s%s\n" % (name, (" " + desc) if desc else ""))
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def identity(a, b, mode="NW"):
    """Alignment identity matches/(alignment columns) via edlib (indel-aware).

    mode NW = global; HW = `a` infix inside `b` (identity over the aligned span of b).
    Returns (identity, edit_distance, (b_start, b_end_inclusive)).
    """
    import edlib
    a = a.upper()
    b = b.upper()
    r = edlib.align(a, b, mode=mode, task="path")
    cig = r["cigar"]
    m = x = ins = dele = 0
    num = ""
    for ch in cig:
        if ch.isdigit():
            num += ch
            continue
        n = int(num)
        num = ""
        if ch == "=":
            m += n
        elif ch == "X":
            x += n
        elif ch == "I":
            ins += n
        elif ch == "D":
            dele += n
    cols = m + x + ins + dele
    loc = r["locations"][0] if r["locations"] else (None, None)
    return (m / cols if cols else 0.0), r["editDistance"], loc


def cons_identity(seq, cons):
    """Identity of a genomic element to a class consensus, robust to ragged ends: the better of
    (consensus infix in element) and (element infix in consensus). Handles 5'-truncated copies
    and annotation spans that include a little flank. Used for the source / novel-source rules."""
    return max(identity(cons, seq, mode="HW")[0], identity(seq, cons, mode="HW")[0])


def ungapped_identity(a, b, mode="NW"):
    """matches / (matches + mismatches) — identity over aligned columns, ignoring indels
    (useful when one sequence has a collapsed tandem repeat, e.g. the Dfam SVA VNTR)."""
    import edlib
    r = edlib.align(a.upper(), b.upper(), mode=mode, task="path")
    m = x = 0
    num = ""
    for ch in r["cigar"]:
        if ch.isdigit():
            num += ch
            continue
        if ch == "=":
            m += int(num)
        elif ch == "X":
            x += int(num)
        num = ""
    return m / (m + x) if m + x else 0.0


def column_map(query, target):
    """Global (NW) alignment of `query` to `target`; returns list q2t where q2t[i] is the
    target index aligned to query base i (or None for an insertion in query)."""
    import edlib
    r = edlib.align(query.upper(), target.upper(), mode="NW", task="path")
    q2t = []
    qi = ti = 0
    num = ""
    for ch in r["cigar"]:
        if ch.isdigit():
            num += ch
            continue
        n = int(num)
        num = ""
        for _ in range(n):
            if ch in "=X":
                q2t.append(ti)
                qi += 1
                ti += 1
            elif ch == "I":      # extra base in query
                q2t.append(None)
                qi += 1
            elif ch == "D":      # extra base in target
                ti += 1
    return q2t


class TwoBit:
    """Thin py2bit wrapper (uppercase, clipped to contig bounds)."""

    def __init__(self, path):
        import py2bit
        self.tb = py2bit.open(path)
        self.chroms = self.tb.chroms()

    def seq(self, chrom, start, end):
        if chrom not in self.chroms:
            return ""
        n = self.chroms[chrom]
        start = max(0, start)
        end = min(n, end)
        if end <= start:
            return ""
        return self.tb.sequence(chrom, start, end).upper()


class Lifter:
    """Interval liftOver with pyliftover (chain files). Lifts both ends and checks that they
    land on one contig, same orientation and a plausible length."""

    def __init__(self, chain):
        from pyliftover import LiftOver
        self.lo = LiftOver(chain)

    def point(self, chrom, pos, same_chrom=True):
        """Lift one base. Only hits on the same chromosome name are accepted by default: the UCSC
        over.chain files also carry small non-syntenic chains that place a repeat copy on a
        paralog elsewhere in the genome (seen for SVAs: chr6 -> chr1), which is never what we
        want for an orthologous position."""
        r = self.lo.convert_coordinate(chrom, int(pos))
        if not r:
            return None
        r = sorted(r, key=lambda x: -x[3])
        for c, p, s, _ in r:
            if not same_chrom or c == chrom:
                return c, int(p), s
        return None

    def interval(self, chrom, start, end, max_len_ratio=1.5, anchor=1000):
        """Returns (chrom, start, end, strand_flip) 0-based half-open or None.
        Both ends must lift to one chromosome with one orientation and a plausible length; if
        the outside anchors (`anchor` bp beyond each end) lift, they must bracket the interval."""
        a = self.point(chrom, start)
        b = self.point(chrom, end - 1)
        # chain gaps at an edge: walk inward to the nearest liftable base (<= 10 % of the
        # interval, >= 50 bp) and extrapolate back to the edge by the same offset
        maxd = max(50, (end - start) // 10)
        if a is None:
            for d in range(5, maxd + 1, 5):
                p = self.point(chrom, start + d)
                if p:
                    a = (p[0], p[1] - d if p[2] == "+" else p[1] + d, p[2])
                    break
        if b is None:
            for d in range(5, maxd + 1, 5):
                p = self.point(chrom, end - 1 - d)
                if p:
                    b = (p[0], p[1] + d if p[2] == "+" else p[1] - d, p[2])
                    break
        if a is None or b is None:
            return None
        if a[0] != b[0] or a[2] != b[2]:
            return None
        lo, hi = sorted((a[1], b[1]))
        L = hi - lo + 1
        if L > max_len_ratio * (end - start) + 50:
            return None
        if anchor:
            al = self.point(chrom, max(0, start - anchor))
            ar = self.point(chrom, end + anchor)
            for an in (al, ar):
                if an is not None and (an[0] != a[0] or an[2] != a[2]):
                    return None
            if al is not None and ar is not None:
                lo_a, hi_a = sorted((al[1], ar[1]))
                if not (lo_a <= lo and hi <= hi_a):
                    return None
        return a[0], lo, hi + 1, a[2] == "-"

    def junction(self, chrom, pos, outward):
        """Lift the base at `pos` (or the nearest liftable one stepping `outward` = +1/-1 up to
        200 bp) — used to place the 3' junction of an element that is absent from the target."""
        for d in range(0, 201, 5):
            p = self.point(chrom, pos + outward * d)
            if p:
                return p
        return None


def load_cytobands(path):
    bands = {}
    with open_text(path) as fh:
        for line in fh:
            c, s, e, name, _ = line.rstrip("\n").split("\t")
            bands.setdefault(c, []).append((int(s), int(e), name))
    return bands


def band_of(bands, chrom, pos):
    for s, e, name in bands.get(chrom, ()):
        if s <= pos < e:
            return chrom.replace("chr", "") + name
    return "."
