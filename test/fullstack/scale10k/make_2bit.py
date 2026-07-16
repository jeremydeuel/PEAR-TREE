#!/usr/bin/env python3
"""Minimal .2bit writer (no faToTwoBit available). Writes the primary hg38
chromosomes so combine_insertions' py2bit-based get_sequence works locally.

2bit spec: signature 0x1A412743, per-seq records of dnaSize, N-blocks, mask-blocks
(we write none), 2-bits/base packed MSB-first with T=0,C=1,A=2,G=3.
"""
import sys
import struct
import numpy as np
import pysam

SIG = 0x1A412743
CODE = np.zeros(256, dtype=np.uint8)
CODE[ord("T")] = 0; CODE[ord("C")] = 1; CODE[ord("A")] = 2; CODE[ord("G")] = 3
IS_ACGT = np.zeros(256, dtype=bool)
for b in b"ACGT":
    IS_ACGT[b] = True


def n_blocks(is_n):
    """Return (starts, sizes) of runs where is_n is True."""
    if not is_n.any():
        return [], []
    d = np.diff(is_n.astype(np.int8))
    starts = list(np.where(d == 1)[0] + 1)
    ends = list(np.where(d == -1)[0] + 1)
    if is_n[0]:
        starts = [0] + starts
    if is_n[-1]:
        ends = ends + [len(is_n)]
    sizes = [e - s for s, e in zip(starts, ends)]
    return starts, sizes


def pack(codes):
    pad = (-len(codes)) % 4
    if pad:
        codes = np.concatenate([codes, np.zeros(pad, dtype=np.uint8)])
    q = codes.reshape(-1, 4).astype(np.uint32)
    b = (q[:, 0] << 6) | (q[:, 1] << 4) | (q[:, 2] << 2) | q[:, 3]
    return b.astype(np.uint8).tobytes()


def main():
    fa_path, out_path = sys.argv[1], sys.argv[2]
    fa = pysam.FastaFile(fa_path)
    names = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    names = [n for n in names if n in fa.references]

    records = []  # (name, dnaSize, nstarts, nsizes, packed)
    for name in names:
        seq = fa.fetch(name).upper()
        arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        is_n = ~IS_ACGT[arr]
        nstarts, nsizes = n_blocks(is_n)
        codes = CODE[arr]
        records.append((name, len(arr), nstarts, nsizes, pack(codes)))
        print(f"  {name}: {len(arr):,} bp, {len(nstarts)} N-blocks", file=sys.stderr)

    # compute offsets
    header = 16
    index = sum(1 + len(n) + 4 for n, *_ in records)
    off = header + index
    offsets = []
    for (name, dna, ns, nz, packed) in records:
        offsets.append(off)
        off += 4 + 4 + 8 * len(ns) + 4 + 0 + 4 + len(packed)

    with open(out_path, "wb") as f:
        f.write(struct.pack("<IIII", SIG, 0, len(records), 0))
        for (name, *_), o in zip(records, offsets):
            f.write(struct.pack("B", len(name)))
            f.write(name.encode("ascii"))
            f.write(struct.pack("<I", o))
        for (name, dna, ns, nz, packed) in records:
            f.write(struct.pack("<I", dna))
            f.write(struct.pack("<I", len(ns)))
            for s in ns:
                f.write(struct.pack("<I", int(s)))
            for z in nz:
                f.write(struct.pack("<I", int(z)))
            f.write(struct.pack("<I", 0))   # maskBlockCount
            f.write(struct.pack("<I", 0))   # reserved
            f.write(packed)
    print(f"wrote {out_path} ({off:,} bytes, {len(records)} seqs)", file=sys.stderr)


if __name__ == "__main__":
    main()
