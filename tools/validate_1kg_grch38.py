#!/usr/bin/env python3
"""Cross-validate the GRCh38 MEI truth set against the phase-3 hs37d5 one.

    venv/bin/python tools/validate_1kg_grch38.py

We replaced the phase-3 (2015, hs37d5, 7x) recall canary with the 1000G 30x GRCh38-native
MEI callset rather than lifting phase 3 over. That is the better source -- but "better" is a
claim, and an untested claim about a truth set is the worst kind: every downstream recall
number inherits it silently.

So: lift phase 3 to GRCh38 and check the two sets agree where they should.

MEASURED RESULT (2026-07-17), and the premise it destroyed:

  bucket           phase3   lifted    found   found%
  ALU/common           70       69       54    78.3%
  ALU/rare          12678    12676     9127    72.0%
  LINE1/common         16       16       10    62.5%
  LINE1/rare         3032     3032     1164    38.4%
  SVA/rare            833      833      381    45.7%
  TOTAL             16631    16628    10737    64.6%

This script originally FAILED at <95% of common phase-3 MEIs recovered, on the reasoning
that a common MEI is real and 30x/3202 samples cannot miss it. That reasoning is WRONG for a
cross-build comparison, and the numbers say so from two directions:

  1. GRCh38 INCORPORATED many hg19-era insertions. An Alu carried by 79% of people is a prime
     candidate to have been added to the reference — and then it is not an insertion relative
     to GRCh38 at all, so a GRCh38 callset is RIGHT to omit it. Confirmed concretely: phase-3
     ALU at hg19 1:3995268 (AF=0.79) lifts to chr1:3935209, which sits INSIDE an AluY
     annotated in hg38 at 3935087-3935379.
  2. But this explains only a MINORITY. Sampling 60 missing vs 60 found phase-3 ALUs and
     asking UCSC whether hg38 annotates an Alu at the lifted site: MISSING 13% (8/60) vs
     FOUND 5% (3/60). Real enrichment (~2.6x), small absolute share.

So ~22% of phase-3 ALUs are absent from the 30x GRCh38 callset and MOST OF THAT IS
UNEXPLAINED. Candidates not distinguished here: phase-3 false positives (7x, 2015);
liftOver imprecision (MEI insertion points sit in repeats, where liftOver is least reliable
— widening the match window 100->5000bp only moved ALU recovery 72%->80%, so this is not the
main term either); or genuinely different MELT/ensemble filtering at 30x.

WHY THE NEW SET IS STILL THE RIGHT CANARY DESPITE THIS. We use the truth set to LABEL our own
contract loci as known-real germline MEIs. What matters is that each entry IS real
(precision); completeness only changes how many anchors we get. A phase-3 record absent from
the 30x set does not create a false label — it just is not an anchor. So this disagreement
argues for reporting anchor counts honestly, NOT against using the GRCh38-native set.

WHAT YOU MUST NOT DO WITH THESE NUMBERS: do not compare any recall figure computed here
against the PD44579 table in config.discovery.grch37 ("0.6%, 7/1144"). Different donor,
different reference, and now a different truth set with a different size and composition.

This script now REPORTS and does not gate. Run it when the truth set is rebuilt or the source
changes; read the table, do not look for a PASS.

hs37d5 contigs 1..22,X,Y have IDENTICAL coordinates to hg19 chr1..chrY (hs37d5 only adds
decoys), so an hg19->hg38 liftOver is the correct transform for these records.
"""
import gzip
import sys
from collections import defaultdict

OLD = "testdata/mei/1kg.sv.vcf.gz"
NEW = "testdata/mei/1kg.mei.grch38.sites.vcf.gz"
WINDOW = 100  # bp; MEI breakpoint calls differ between callsets by tens of bp routinely

try:
    from pyliftover import LiftOver
except ImportError:
    sys.exit("pyliftover missing — run with venv/bin/python")


def info(s, key):
    for f in s.split(";"):
        if f.startswith(key + "="):
            return f[len(key) + 1:]
    return None


def cls_of(alt_or_svtype):
    u = alt_or_svtype.upper()
    if "ALU" in u:
        return "ALU"
    if "LINE1" in u or "L1" in u:
        return "LINE1"
    if "SVA" in u:
        return "SVA"
    return None


# ---- new set: class -> chrom -> sorted positions
new = defaultdict(lambda: defaultdict(list))
n_new = 0
with gzip.open(NEW, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        c = cls_of(f[4])
        if c:
            new[c][f[0]].append(int(f[1]))
            n_new += 1
for c in new:
    for ch in new[c]:
        new[c][ch].sort()
print(f"new set: {n_new} MEIs  ({', '.join(f'{k}={sum(len(v) for v in new[k].values())}' for k in sorted(new))})")

import bisect


def near(c, ch, pos):
    """Any new-set MEI of class c within +/-WINDOW of pos on ch?"""
    arr = new.get(c, {}).get(ch)
    if not arr:
        return False
    i = bisect.bisect_left(arr, pos - WINDOW)
    return i < len(arr) and arr[i] <= pos + WINDOW


print("loading hg19->hg38 chain (pyliftover downloads it on first use)...")
lo = LiftOver("hg19", "hg38")

stats = defaultdict(lambda: [0, 0, 0])  # bucket -> [total, lifted, found]
unlifted_examples, missing_examples = [], []

with gzip.open(OLD, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        svt = info(f[7], "SVTYPE") or ""
        c = cls_of(svt)
        if not c or svt.startswith("DEL"):   # DEL_ALU / DEL_LINE1 are not insertions
            continue
        af = float(info(f[7], "AF") or 0)
        bucket = f"{c}/{'common' if af >= 0.5 else 'rare'}"
        stats[bucket][0] += 1

        conv = lo.convert_coordinate("chr" + f[0], int(f[1]))
        if not conv:
            if len(unlifted_examples) < 3:
                unlifted_examples.append(f"{f[0]}:{f[1]} {svt} AF={af}")
            continue
        stats[bucket][1] += 1
        ch38, pos38 = conv[0][0], conv[0][1]
        if near(c, ch38, pos38):
            stats[bucket][2] += 1
        elif af >= 0.5 and len(missing_examples) < 5:
            missing_examples.append(f"{f[0]}:{f[1]} -> {ch38}:{pos38} {svt} AF={af:.2f}")

print()
print(f"phase-3 MEIs lifted hg19->hg38 and looked up in the new set (+/-{WINDOW}bp, same class)")
print(f"{'bucket':<16} {'phase3':>8} {'lifted':>8} {'found':>8} {'found%':>8}")
tot = [0, 0, 0]
for b in sorted(stats):
    t, l, fnd = stats[b]
    tot = [tot[0] + t, tot[1] + l, tot[2] + fnd]
    print(f"{b:<16} {t:>8} {l:>8} {fnd:>8} {100*fnd/l if l else 0:>7.1f}%")
print(f"{'TOTAL':<16} {tot[0]:>8} {tot[1]:>8} {tot[2]:>8} {100*tot[2]/tot[1] if tot[1] else 0:>7.1f}%")

if unlifted_examples:
    print(f"\nunliftable phase-3 records (a cost the GRCh38-native set simply does not have):")
    for e in unlifted_examples:
        print("  " + e)
if missing_examples:
    print("\nCOMMON phase-3 MEIs absent from the new set:")
    for e in missing_examples:
        print("  " + e)
    print("  Some are correct: GRCh38 incorporated common hg19-era insertions, so they are no")
    print("  longer insertions relative to it (verified: 1:3995268 AF=0.79 -> chr1:3935209,")
    print("  inside an hg38-annotated AluY at 3935087-3935379). Measured on 60 missing vs 60")
    print("  found ALUs, that mechanism covers 13% vs a 5% control — real, but a minority.")

com_l = sum(stats[b][1] for b in stats if b.endswith("/common"))
com_f = sum(stats[b][2] for b in stats if b.endswith("/common"))
rate = 100 * com_f / com_l if com_l else 0
print()
print("-" * 78)
print(f"{rate:.1f}% of lifted COMMON phase-3 MEIs, and "
      f"{100*tot[2]/tot[1] if tot[1] else 0:.1f}% overall, are in the GRCh38 set.")
print()
print("This is REPORTED, not gated — see the module docstring. Cross-build MEI callsets are")
print("NOT expected to agree: GRCh38 absorbed many hg19-era insertions, phase 3 was called at")
print("7x in 2015, and MEI breakpoints sit in repeats where liftOver is least reliable. A")
print("phase-3 record missing here costs us an ANCHOR, not correctness — each entry in the new")
print("set is still a real germline MEI, which is all the canary requires.")
print()
print("Treat the gap as a known unknown: ~22% of phase-3 ALUs are absent and most of that is")
print("unexplained. If a future recall number looks surprising, suspect the truth set first.")
