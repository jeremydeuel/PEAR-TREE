#!/usr/bin/env python3
"""Emit euchromatin FP-decoy windows: contiguous blocks spread across the autosomes,
disjoint from the TP windows and the satellite/centromere compartment, so decoy coverage
is uniform (coverage_mask behaves as on real WGS). Writes fp_eu_windows.bed.

Usage: gen_fp_windows.py <hs1.fai> <tp_windows.bed> <sat_fp_windows.bed> <out.bed> [target_mb]
"""
import sys

FAI, TP, SAT, OUT = sys.argv[1:5]
TARGET_MB = float(sys.argv[5]) if len(sys.argv) > 5 else 40.0
BLOCK = 2_000_000
END_MARGIN = 3_000_000
MAIN = [f"chr{i}" for i in range(1, 23)]

lens = {}
with open(FAI) as f:
    for line in f:
        c, n = line.split()[:2]
        lens[c] = int(n)

excl = {}
for path in (TP, SAT):
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 3:
                excl.setdefault(p[0], []).append((int(p[1]), int(p[2])))


def overlaps(c, s, e):
    return any(not (e <= xs or s >= xe) for xs, xe in excl.get(c, ()))


out, total = [], 0
# two blocks per chromosome (p-arm-ish and q-arm-ish), skipping any that collide
for c in MAIN:
    if c not in lens or total >= TARGET_MB * 1e6:
        continue
    L = lens[c]
    for frac in (0.30, 0.70, 0.50, 0.85):
        if total >= TARGET_MB * 1e6:
            break
        s = int(L * frac)
        s = max(END_MARGIN, min(s, L - END_MARGIN - BLOCK))
        e = s + BLOCK
        if e > L - END_MARGIN or overlaps(c, s, e) or any(cc == c and not (e <= xs or s >= xe)
                                                          for cc, xs, xe in out):
            continue
        out.append((c, s, e))
        total += e - s

out.sort()
with open(OUT, "w") as f:
    for c, s, e in out:
        f.write(f"{c}\t{s}\t{e}\tFPEU\n")
print(f"FP euchromatin windows: {len(out)} blocks, {total/1e6:.1f} Mb")
for c, s, e in out:
    print(f"  FPEU {c}:{s}-{e}")
