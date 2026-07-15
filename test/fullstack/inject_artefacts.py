#!/usr/bin/env python3
"""Inject library-prep artefacts into a wgsim paired-FASTQ.

The assembly-discordance and repeat/centromere/telomere false positives come for free
from mapping hs1-derived reads to hg38; this adds the *library* artefacts that arise
before mapping (they are read-intrinsic, so a simulator that only draws from the donor
would never produce them):

  * adapter read-through (insert < read length)   -> 3' tail is the sequencing adapter
  * poly-G / two-colour dark cycle                 -> 3' tail is a G homopolymer
  * structure-specific chimera (fold-back)         -> read = seq + revcomp(seq): a palindrome
  * PCR duplicate                                  -> the pair is emitted twice (markdup flags it)

Rates are low and tunable; each modified read gets an `XA:art:<kind>` suffix in its name so
the truth of which reads are artefactual is recoverable.
"""
import argparse
import random

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAGATCGGAAGAGCGTCGTGTAGGGAAAGA"


def read_fq(fh):
    while True:
        h = fh.readline()
        if not h:
            return
        s = fh.readline().rstrip("\n")
        fh.readline()
        q = fh.readline().rstrip("\n")
        yield h.rstrip("\n"), s, q


def emit(o1, o2, n, s1, q1, s2, q2, suffix=""):
    o1.write(f"@{n}{suffix}/1\n{s1}\n+\n{q1}\n")
    o2.write(f"@{n}{suffix}/2\n{s2}\n+\n{q2}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in1", required=True); p.add_argument("--in2", required=True)
    p.add_argument("--out1", required=True); p.add_argument("--out2", required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--adapter-rate", type=float, default=0.010)
    p.add_argument("--polyg-rate", type=float, default=0.006)
    p.add_argument("--chimera-rate", type=float, default=0.004)
    p.add_argument("--dup-rate", type=float, default=0.030)
    args = p.parse_args()

    rng = random.Random(args.seed)
    f1 = open(args.in1); f2 = open(args.in2)
    o1 = open(args.out1, "w"); o2 = open(args.out2, "w")
    n_ad = n_pg = n_ch = n_du = n = 0

    for (h1, s1, q1), (h2, s2, q2) in zip(read_fq(f1), read_fq(f2)):
        n += 1
        name = h1.split()[0].lstrip("@").rstrip("/1")
        L = len(s1)
        suffix = ""
        r = rng.random()
        if r < args.adapter_rate:                       # read-through adapter on read2
            k = rng.randint(40, L - 15)
            s2 = (s2[:k] + ADAPTER)[:L]; q2 = "I" * len(s2)
            suffix = ":art_adapter"; n_ad += 1
        elif r < args.adapter_rate + args.polyg_rate:   # poly-G tail on read2
            k = rng.randint(60, L - 15)
            s2 = (s2[:k] + "G" * L)[:L]; q2 = "I" * len(s2)
            suffix = ":art_polyg"; n_pg += 1
        elif r < args.adapter_rate + args.polyg_rate + args.chimera_rate:  # fold-back chimera
            half = L // 2
            s1 = (s1[:half] + rc(s1[:half]))[:L]; q1 = "I" * len(s1)
            suffix = ":art_chimera"; n_ch += 1
        emit(o1, o2, name, s1, q1, s2, q2, suffix)
        if rng.random() < args.dup_rate:                # PCR duplicate (emit again)
            emit(o1, o2, name + "_dup", s1, q1, s2, q2, ":art_dup"); n_du += 1

    for fh in (f1, f2, o1, o2):
        fh.close()
    print(f"pairs={n} adapter={n_ad} polyg={n_pg} chimera={n_ch} dup={n_du}")


if __name__ == "__main__":
    main()
