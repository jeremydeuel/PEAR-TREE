#!/usr/bin/env python3
"""Score the genotyping panel.

Given the lifted hg38 truth, the genotyping contract (loci actually genotyped), each
sample's genotype calls, and the final combine_genotypes matrix, report:

  1. Contract stage  -- how many true insertions reached genotyping (discovery+combine recall)
     and how many contract loci are assembly-discordance false positives.
  2. Het-recall vs depth, per element class (pseudogene broken out) -- the sensitivity curve
     the 5-100x ladder exists to measure. A trailing pseudogene row = parent-gene cross-mapping
     stealing alt reads.
  3. The 3-negative discriminator -- true loci must read wild-type in the reference-only samples.
  4. Panel outcome -- combine_genotypes must PASS every true locus and DROP every FP.

Usage: score_genotyping.py --truth truth_hg38.tsv --contract step2.genotyping.txt.gz
       --matrix matrix.csv.gz --genotypes a.genotypes.txt.gz ...
       --het-depths "100 80 60 40 20 10 5" --disc-bam het_100x.bam --window 50
"""
import argparse
import gzip
import os
import re
from collections import defaultdict

PSEUDOGENE_SYMS = {"DUX4", "MALAT1", "HNRNPA1", "CASP12", "RPL21"}
NAME_RE = re.compile(r"^([^:]+):(?:polyA_|disc_)?(\d+)-(?:polyA_|disc_)?(\d+)")
HET = {"heterozygous", "homozygous"}
PRESENT = {"heterozygous", "homozygous", "insertion"}   # colony carries the insertion (clade signal)
UNCERTAIN = {"insertion?", "wild-type?"}
NA_LIKE = {"high-coverage", "no-coverage", "error"}


def klass(family):
    return "pseudogene" if family in PSEUDOGENE_SYMS else family


def sample_stem(path):
    return re.sub(r"(\.genotypes)?(\.(txt|csv))?(\.gz)?$", "", os.path.basename(path))


PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_truth(path):
    rows = []
    skipped_alt = 0
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if r.get("status") != "scoreable":
                continue
            # A truth chain-lifted onto a non-primary scaffold (e.g. *_alt) is unscoreable:
            # reads map to the primary homolog, so no primary-contig call can ever match it.
            # Excluding it keeps recall honest (it is a lift artifact, not a discovery miss).
            if r.get("hg38_contig") not in PRIMARY:
                skipped_alt += 1
                continue
            r["L"], r["R"] = int(r["hg38_left"]), int(r["hg38_right"])
            rows.append(r)
    if skipped_alt:
        print(f"(excluded {skipped_alt} truth(s) lifted to non-primary scaffolds -- unscoreable)")
    return rows


def load_contract(path):
    names = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                names.append(line[1:].rstrip("\n"))
    return names


def load_genotypes(path):
    calls = {}
    with gzip.open(path, "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                calls[p[0]] = p[1]
    return calls


def parse_locus(name):
    m = NAME_RE.match(name)
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def classify(names, truth, window):
    """Map each contract locus -> ('TP', family) if it matches a scoreable truth, else ('FP', None)."""
    tby = defaultdict(list)
    for t in truth:
        tby[t["hg38_contig"]].append(t)
    locus_class, matched_truth = {}, set()
    for name in names:
        pl = parse_locus(name)
        if pl is None:
            locus_class[name] = ("FP", None)
            continue
        c, l, r = pl
        best, bestd = None, window + 1
        for t in tby.get(c, ()):
            d = min(abs(l - t["L"]), abs(l - t["R"]), abs(r - t["L"]), abs(r - t["R"]))
            if d <= window and d < bestd:
                best, bestd = t, d
        if best is not None:
            locus_class[name] = ("TP", best["family"])
            matched_truth.add(best["id"])
        else:
            locus_class[name] = ("FP", None)
    return locus_class, matched_truth


def pct(n, d):
    return f"{100*n/d:5.1f}%" if d else "   n/a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--contract", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--genotypes", nargs="+", required=True)
    ap.add_argument("--het-depths", default="100 80 60 40 20 10 5")
    ap.add_argument("--disc-bam", default="")
    ap.add_argument("--window", type=int, default=50)
    args = ap.parse_args()

    truth = load_truth(args.truth)
    names = load_contract(args.contract)
    locus_class, matched_truth = classify(names, truth, args.window)
    gts = {sample_stem(p): load_genotypes(p) for p in args.genotypes}
    het_depths = [int(x) for x in args.het_depths.split()]

    het_samples = sorted((s for s in gts if s.startswith("het_")),
                         key=lambda s: -int(re.search(r"_(\d+)x", s).group(1)))
    wt_samples = sorted(s for s in gts if s.startswith("wt_"))

    tp_names = [n for n in names if locus_class[n][0] == "TP"]
    fp_names = [n for n in names if locus_class[n][0] == "FP"]

    print(f"=== genotyping panel score (window +-{args.window} bp) ===\n")

    # 1. contract stage
    print("1. CONTRACT STAGE (discovery on the 100x sample + combine_insertions)")
    print(f"   scoreable truths          : {len(truth)}")
    print(f"   true insertions in contract: {len(matched_truth)}  "
          f"(recall {pct(len(matched_truth), len(truth))}) -- rest are discovery/lift misses")
    print(f"   contract loci             : {len(names)}  ({len(tp_names)} TP, {len(fp_names)} FP)")
    print(f"   FP = hs1->hg38 assembly discordance that survived combine_insertions Filter A\n")

    # 2. het-recall vs depth, per class
    classes = sorted({klass(locus_class[n][1]) for n in tp_names})
    tp_by_class = defaultdict(list)
    for n in tp_names:
        tp_by_class[klass(locus_class[n][1])].append(n)

    print("2. PRESENCE-RECALL vs DEPTH  (fraction of TP loci called het/hom/insertion per sample)")
    print("   (presence, not zygosity, is the clade signal -- 'insertion' = present, zygosity unclear)")
    header = "   " + f"{'class':12s}" + "  n   " + "".join(f"{d:>7}x" for d in het_depths)
    print(header)
    for cl in classes + ["OVERALL"]:
        loci = tp_names if cl == "OVERALL" else tp_by_class[cl]
        row = f"   {cl:12s}  {len(loci):<4d}"
        for d in het_depths:
            s = f"het_{d}x"
            calls = gts.get(s, {})
            hit = sum(1 for n in loci if calls.get(n) in PRESENT)
            row += f"{pct(hit, len(loci)):>8}"
        print(row)
    # aggregate call-type breakdown at the lowest and a mid depth (where degradation shows)
    print("\n   call-type mix at low depth (TP loci):")
    for d in het_depths:
        s = f"het_{d}x"; calls = gts.get(s, {})
        c = defaultdict(int)
        for n in tp_names:
            g = calls.get(n, "missing")
            key = ("het/hom" if g in HET else "insertion(zyg?)" if g == "insertion"
                   else "uncertain" if g in UNCERTAIN
                   else "wild-type" if g == "wild-type" else "NA" if g in NA_LIKE else g)
            c[key] += 1
        mix = "  ".join(f"{k}:{v}" for k, v in sorted(c.items(), key=lambda x: -x[1]))
        print(f"     het_{d}x : {mix}")
    print()

    # 3. the 3-negative discriminator
    print("3. NEGATIVE DISCRIMINATOR (TP loci must read wild-type in the reference-only samples)")
    print(f"   wild-type samples: {', '.join(wt_samples)}")
    allwt = 0
    for n in tp_names:
        if all(gts.get(s, {}).get(n) == "wild-type" for s in wt_samples):
            allwt += 1
    print(f"   TP loci wild-type in ALL {len(wt_samples)} negatives: {allwt}/{len(tp_names)} "
          f"({pct(allwt, len(tp_names))})")
    # contrast: FP loci should NOT read wild-type across the negatives (they carry the same
    # discordance signal in every sample) -> this is why combine_genotypes can reject them.
    fp_allwt = sum(1 for n in fp_names
                   if all(gts.get(s, {}).get(n) == "wild-type" for s in wt_samples))
    print(f"   FP loci wild-type in ALL negatives : {fp_allwt}/{len(fp_names)} "
          f"(low is good -- FPs look the same everywhere)\n")

    # 4. panel outcome
    passing = set()
    try:
        import pandas as pd
        m = pd.read_csv(args.matrix, sep=";", index_col=0)
        passing = set(m.index.astype(str))
    except Exception as e:  # noqa: BLE001
        print(f"   (could not read matrix {args.matrix}: {e})")
    tp_pass = sum(1 for n in tp_names if n in passing)
    fp_leak = sum(1 for n in fp_names if n in passing)
    print("4. PANEL OUTCOME (combine_genotypes filters)")
    print(f"   TP loci passing filters : {tp_pass}/{len(tp_names)}  ({pct(tp_pass, len(tp_names))})  [want ~all]")
    print(f"   FP loci leaked to output: {fp_leak}/{len(fp_names)}  ({pct(fp_leak, len(fp_names))})  [want 0]")
    print(f"   final matrix rows       : {len(passing)}")
    print()
    print("   verdict: PASS if TP-pass is high AND FP-leak is 0. A trailing pseudogene row in (2)")
    print("            means parent-gene cross-mapping is eroding VAF -- tune vaf bands / min_supporting_reads.")


if __name__ == "__main__":
    main()
