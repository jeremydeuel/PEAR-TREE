#!/usr/bin/env python3
"""Emit FP compartment (centromere/satellite bins + telomere ends) and choose TP
euchromatin windows away from satellite. Writes fp_windows.bed and tp_windows.bed."""
import sys

SATBINS, FAI, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
SAT_THR = 400_000        # min satellite bp in a 1Mb bin to count as FP compartment
TELO = 500_000           # terminal window length per chromosome end
MAIN = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

lens = {}
with open(FAI) as f:
    for line in f:
        c, n = line.split()[:2]
        lens[c] = int(n)

# satellite bins -> merge consecutive
sat = {}
with open(SATBINS) as f:
    for line in f:
        c, b, v = line.split()
        if int(v) >= SAT_THR:
            sat.setdefault(c, set()).add(int(b))

fp = []  # (contig, start, end, label)
for c, bins in sat.items():
    for b in sorted(bins):
        s, e = b * 1_000_000, min((b + 1) * 1_000_000, lens[c])
        if fp and fp[-1][0] == c and fp[-1][2] == s:
            fp[-1] = (c, fp[-1][1], e, "centromere")
        else:
            fp.append((c, s, e, "centromere"))

# telomere ends (skip if already inside a satellite window)
def in_fp(c, s, e):
    return any(fc == c and not (e <= fs or s >= fe) for fc, fs, fe, _ in fp)

for c in MAIN:
    if c not in lens:
        continue
    for s, e, lab in [(0, TELO, "telo5"), (lens[c] - TELO, lens[c], "telo3")]:
        if not in_fp(c, s, e):
            fp.append((c, s, e, lab))

fp.sort()
fp_bp = sum(e - s for _, s, e, _ in fp)
with open(f"{OUT}/fp_windows.bed", "w") as f:
    for c, s, e, lab in fp:
        f.write(f"{c}\t{s}\t{e}\t{lab}\n")

# TP euchromatin windows: mid-arm 6 Mb blocks on several autosomes, each checked
# to be satellite-free (no overlap with an FP window).
CAND = [
    ("chr21", 20_000_000, 27_000_000), ("chr22", 24_000_000, 31_000_000),
    ("chr2", 100_000_000, 107_000_000), ("chr4", 80_000_000, 87_000_000),
    ("chr6", 100_000_000, 107_000_000), ("chr7", 60_000_000, 67_000_000),
]
tp = [(c, s, e) for (c, s, e) in CAND if c in lens and e <= lens[c] and not in_fp(c, s, e)]
tp_bp = sum(e - s for _, s, e in tp)
with open(f"{OUT}/tp_windows.bed", "w") as f:
    for c, s, e in tp:
        f.write(f"{c}\t{s}\t{e}\tTP\n")

print(f"FP compartment: {len(fp)} windows, {fp_bp/1e6:.1f} Mb")
print(f"TP windows: {len(tp)} windows, {tp_bp/1e6:.1f} Mb")
for c, s, e in tp:
    print(f"  TP {c}:{s}-{e}")
