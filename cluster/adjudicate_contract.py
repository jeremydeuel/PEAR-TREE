#!/usr/bin/env python3
"""Score ONE contract against the frozen 9x10 truth set. Prints a one-line verdict
(TP kept/dropped, FP kept) so a threshold sweep can read results on the farm without
a download round-trip.

Usage: adjudicate_contract.py <contract.genotyping.txt.gz> <frozen_truth.tsv> [arm_label]

The truth labels (col 'label_v2'): TP, FP_tree_break, FP_nonMEI_SV, exclude. A contract
locus is a '>chr:left-right' header line. Kept = truth locus present in the contract.
"""
import csv, gzip, sys

def contract_loci(path):
    op = gzip.open if path.endswith(".gz") else open
    s = set()
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                s.add(line[1:].strip())
    return s

def truth_sets(path):
    tp, fpt, fps = set(), set(), set()
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            lab = row["label_v2"]; loc = row["locus"]
            if lab == "TP": tp.add(loc)
            elif lab == "FP_tree_break": fpt.add(loc)
            elif lab == "FP_nonMEI_SV": fps.add(loc)
    return tp, fpt, fps

def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    contract, truth = sys.argv[1], sys.argv[2]
    arm = sys.argv[3] if len(sys.argv) > 3 else contract.split("/")[-2] if "/" in contract else "contract"
    loci = contract_loci(contract)
    tp, fpt, fps = truth_sets(truth)
    tpk = len(tp & loci); tpd = len(tp) - tpk
    fptk = len(fpt & loci); fpsk = len(fps & loci)
    recall = 100 * tpk / len(tp) if tp else 0.0
    print(f"{arm}\tnloci={len(loci)}\tTPkept={tpk}\tTPdrop={tpd}\trecall={recall:.1f}%"
          f"\tFP_tree={fptk}\tFP_sv={fpsk}\tFPtot={fptk+fpsk}")

if __name__ == "__main__":
    main()
