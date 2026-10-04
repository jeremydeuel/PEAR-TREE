"""Sequence utilities, genome access and literature-calibrated parameter distributions.

Every distribution below names its source. Literature notes used for calibration:
  Szak 2002 (Genome Biol 3:research0052)  TSD length (median ~15 bp, 4-25 typical), poly-A
                                          length; EN motif degeneracy
  Zumalave 2026 (bioRxiv, "Synchronous L1 retrotransposition ...") solo-L1 size (median 878
                                          bp), 49% 5'-del only / 51% 5'-del + 5'-inv; twin
                                          priming: poly(dT)-primed part ~2.6x the internally
                                          primed part (median 1.51 vs 0.54 kb); 66% junction
                                          deletion (median 14.5 bp), 17% duplication (median 22);
                                          transduction termination at fixed per-source offsets
                                          (e.g. 234 bp, 64 bp downstream = poly-A signal use)
  Nam 2023 (Nature 617:540)               colon soL1R: 29.5% 5'-inversion, 0.3% fold-back
                                          inverted duplication 5' of the target site; target-site
                                          deletions (negative TSD) observed; templated / pre-mRNA
                                          co-insertions
  Ostertag & Kazazian 2001 (Genome Res)   twin priming model
  Tubio 2014 (Science 345:1251343)        3' transductions: mostly < 1 kb, up to ~12 kb;
                                          orphan transductions
  Rodriguez-Martin 2020 (Nat Genet 52:306) L1-mediated deletions/duplications, orphan TDs
  Flasch 2019 (Cell 177:837)              EN cleavage-site preference (TTTT/AA degenerate)
  Damert 2009 (Genome Res 19:1992)        SVA 5' transduction (upstream flank)
  Ewing 2013 (Genome Biol 14:R22)         processed pseudogene insertions (spliced exons)
"""
import math
import random

BASES = "ACGT"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def rnd_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def read_fasta(path):
    """Ordered dict name -> sequence (upper case kept as-is; caller may .upper())."""
    import gzip
    recs, name, buf = {}, None, []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs[name] = "".join(buf)
    return recs


def write_fasta(path, recs, width=0):
    with open(path, "w") as f:
        for name, seq in recs.items():
            if width:
                f.write(f">{name}\n")
                for i in range(0, len(seq), width):
                    f.write(seq[i:i + width] + "\n")
            else:
                f.write(f">{name}\n{seq}\n")


class Genome:
    """Read-only genome accessor over a .2bit (py2bit) or an indexed FASTA (pysam).
    Coordinates are 0-based half-open; sequence returned upper-case."""

    def __init__(self, path):
        self.path = str(path)
        if self.path.endswith(".2bit"):
            import py2bit
            self._tb = py2bit.open(self.path)
            self._fa = None
            self.lengths = dict(self._tb.chroms())
        else:
            import pysam
            self._fa = pysam.FastaFile(self.path)
            self._tb = None
            self.lengths = dict(zip(self._fa.references, self._fa.lengths))

    def fetch(self, contig, start, end, strand="+"):
        start = max(0, start)
        end = min(end, self.lengths[contig])
        if self._tb is not None:
            s = self._tb.sequence(contig, start, end).upper()
        else:
            s = self._fa.fetch(contig, start, end).upper()
        return revcomp(s) if strand == "-" else s


def homopolymer_runs(seq, min_len=6):
    """[(start, end, base)] of homopolymer runs >= min_len."""
    runs = []
    i, n = 0, len(seq)
    while i < n:
        j = i + 1
        while j < n and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            runs.append((i, j, seq[i]))
        i = j
    return runs


# ----------------------------------------------------------------------------------------
# parameter distributions
# ----------------------------------------------------------------------------------------
def lognormal_int(rng, median, sigma, lo, hi):
    v = int(round(math.exp(math.log(median) + sigma * rng.gauss(0, 1))))
    return min(max(v, lo), hi)


def sample_tsd_len(rng):
    """TSD length: peak ~15 bp, range 4-25 (Szak 2002; Nam 2023 Ext. Data Fig. 1e)."""
    if rng.random() < 0.85:
        v = int(round(rng.gauss(15, 3.0)))
    else:
        v = rng.randint(4, 25)
    return min(max(v, 4), 25)


# poly-A tail length per class. L1: median ~70, range 15-635 (Szak 2002 / long-read
# somatic L1 catalogues); Alu / SVA tails are shorter. Reads only see what fits.
POLYA_PARAMS = {
    "L1": (70, 0.75, 15, 635),
    "ALU": (30, 0.45, 10, 120),
    "SVA": (40, 0.55, 10, 200),
    "PSEUDOGENE": (60, 0.7, 15, 400),
    "POLYA_ONLY": (45, 0.6, 15, 300),
    "ORPHAN_TD": (60, 0.7, 15, 400),
}


def sample_polya_len(rng, cls="L1", scale=1.0):
    med, sig, lo, hi = POLYA_PARAMS.get(cls, POLYA_PARAMS["L1"])
    return lognormal_int(rng, med * scale, sig, lo, hi)


def make_polya(rng, n, purity=0.985):
    """Poly-A tail of length n; real tails carry occasional non-A interruptions
    (e.g. AAAAGAAAA), mostly in the distal half."""
    out = []
    for i in range(n):
        if i > 12 and rng.random() > purity:
            out.append(rng.choice("GCT"))
        else:
            out.append("A")
    return "".join(out)


def sample_l1_truncated_len(rng, full_len):
    """5'-truncated L1 insert length: 3'-clustered, median ~0.7 kb (Zumalave 2026 solo
    median 878 bp incl. full-length; Nam 2023 colon soL1R median ~0.5 kb)."""
    return lognormal_int(rng, 700, 0.85, 60, full_len - 60)


def sample_twin_priming(rng, full_len):
    """Twin-priming geometry on an L1 of length full_len (sense coords 0..full_len).

    Returns dict(c, a, b, junction, jlen):
      forward (poly-dT primed) part = L1[c:full_len]  (c = inversion breakpoint, >= 590)
      inverted (internally primed) part = L1[a:b], inserted reverse-complemented 5' of it
      junction: 'deletion' (b < c, gap jlen, 66%, median 14.5 bp) | 'duplication' (b > c,
      overlap jlen, 17%, median 22 bp) | 'clean' (b == c, 17%)        (Zumalave 2026)
      forward : inverted length ratio median ~2.3 (Zumalave: 1.51 vs 0.54 kb medians)
    """
    for _ in range(100):
        fwd = lognormal_int(rng, 1500, 0.6, 120, full_len - 590)
        c = full_len - fwd
        if c < 590:                      # breakpoint never in the first 590 bp of L1
            continue
        ratio = math.exp(math.log(2.3) + 0.6 * rng.gauss(0, 1))
        inv = max(40, int(fwd / ratio))
        u = rng.random()
        if u < 0.66:
            junction, jlen = "deletion", lognormal_int(rng, 14.5, 0.7, 1, 120)
            b = c - jlen
        elif u < 0.83:
            junction, jlen = "duplication", lognormal_int(rng, 22, 0.6, 2, 120)
            b = c + jlen
        else:
            junction, jlen = "clean", 0
            b = c
        a = b - inv
        if a < 0 or b > full_len:
            continue
        return {"c": c, "a": a, "b": b, "junction": junction, "jlen": jlen}
    # degenerate fallback (very short element)
    c = max(590, full_len // 2)
    return {"c": c, "a": max(0, c - 300), "b": c, "junction": "clean", "jlen": 0}


# L1 EN motif: bottom-strand 5'-TTTT/AA-3' == top-strand 5'-TT/AAAA-3' with the nick
# between TT and AAAA (Flasch 2019; Jurka 1997). Most real sites are degenerate.
EN_CONSENSUS = "TTAAAA"     # top strand around the nick, nick after index 2
EN_NICK_OFFSET = 2
EN_MISMATCH_DIST = [(0, 0.25), (1, 0.35), (2, 0.25), (3, 0.15)]
_TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def sample_en_mismatches(rng):
    u, acc = rng.random(), 0.0
    for k, p in EN_MISMATCH_DIST:
        acc += p
        if u <= acc:
            return k
    return EN_MISMATCH_DIST[-1][0]


def degenerate_en_motif(rng, k):
    """EN_CONSENSUS with k mismatches (mostly transitions, e.g. TTAAAA -> TCAAAA/TTAGAA)."""
    m = list(EN_CONSENSUS)
    for i in rng.sample(range(len(m)), k):
        m[i] = _TRANSITION[m[i]] if rng.random() < 0.75 else rng.choice([b for b in BASES if b != m[i]])
    return "".join(m)


def en_mismatches(hexamer):
    return sum(1 for a, b in zip(hexamer.upper(), EN_CONSENSUS) if a != b)


def sample_deletion_len(rng, lo=100, hi=20000):
    """L1-mediated deletion / duplication size: log-uniform (Gilbert 2002: bp to 24 kb;
    Rodriguez-Martin 2020 up to Mb — capped here for window-based simulation)."""
    return int(math.exp(rng.uniform(math.log(lo), math.log(hi))))
