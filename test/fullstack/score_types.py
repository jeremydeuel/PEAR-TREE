#!/usr/bin/env python3
"""Score discovery calls against a donor_types.py truth table, per insertion type.

1. Lift: each event's 300 bp hs1 flanks (flanks.fa, bwa-mapped to the discovery reference
   -> --flanks-bam) place its breakpoints in reference coordinates: lo = end of the left
   flank, hi = start of the right flank (same contig, MAPQ >= --min-mapq). One mapped flank
   is enough when the event's breakpoint span is known (hi = lo + |tsd|).
2. Match: a call `contig:L-R` matches an event when both breakpoints agree within --window
   (order-free, so deletions with L > R match); `one-sided` = only one breakpoint agrees.
3. Report recall per variant (pooled over all colony samples = patient level, and per
   sample), calls at labelled artefacts, and unexplained calls (FP). Junction-support
   counts from simulate_reads.py (<prefix>.support.tsv) are merged into --out-tsv.
"""
import argparse
import gzip
import re
from collections import defaultdict

import pysam

CALL_RE = re.compile(r"^@([^:]+):(?:polyA_|disc_)?(\d+)-(?:polyA_|disc_)?(\d+):")


def parse_calls(path):
    calls = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            m = CALL_RE.match(line)
            if m:
                calls.add((m.group(1), int(m.group(2)), int(m.group(3))))
    return calls


def lift(flanks_bam, min_mapq):
    aln = {}
    for r in pysam.AlignmentFile(flanks_bam):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < min_mapq:
            continue
        rid, side = r.query_name.rsplit("_", 1)
        aln[(rid, side)] = (r.reference_name, r.reference_start, r.reference_end, r.is_reverse)
    return aln


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth", required=True)
    p.add_argument("--flanks-bam", required=True)
    p.add_argument("--calls", required=True, help="comma list of per-sample discovery .txt.gz (S1,S2,..)")
    p.add_argument("--support", default="", help="comma list of per-sample <prefix>.support.tsv")
    p.add_argument("--window", type=int, default=30)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--out-tsv", default=None, help="per-event results (+ lifted coords, support)")
    a = p.parse_args()

    rows = []
    with open(a.truth) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    aln = lift(a.flanks_bam, a.min_mapq)
    for r in rows:
        l, rr = aln.get((r["id"], "L")), aln.get((r["id"], "R"))
        span = abs(int(r["tsd"]))
        r["ref_contig"] = r["ref_lo"] = r["ref_hi"] = None
        if l and rr and l[0] == rr[0] and not l[3] and not rr[3]:
            r["ref_contig"], r["ref_lo"], r["ref_hi"] = l[0], l[2], rr[1]
        elif l and rr and l[0] == rr[0] and l[3] and rr[3]:       # window maps reverse
            r["ref_contig"], r["ref_lo"], r["ref_hi"] = l[0], rr[2], l[1]
        elif l and not l[3]:
            r["ref_contig"], r["ref_lo"], r["ref_hi"] = l[0], l[2], l[2] + span
        elif rr and not rr[3]:
            r["ref_contig"], r["ref_lo"], r["ref_hi"] = rr[0], rr[1] - span, rr[1]
    call_sets = [parse_calls(c) for c in a.calls.split(",") if c]
    pooled = set().union(*call_sets) if call_sets else set()
    support = [defaultdict(lambda: (0, 0)) for _ in call_sets]
    for si, sp in enumerate([s for s in a.support.split(",") if s]):
        with open(sp) as f:
            next(f)
            for line in f:
                eid, side, fr, rd = line.rstrip("\n").split("\t")
                support[si][(eid, side)] = (int(fr), int(rd))

    def match(c, r):
        if r["ref_contig"] is None or c[0] != r["ref_contig"]:
            return 0
        cl, ch = sorted(c[1:])
        lo, hi = sorted((r["ref_lo"], r["ref_hi"]))
        both = abs(cl - lo) <= a.window and abs(ch - hi) <= a.window
        if both:
            return 2
        return 1 if any(abs(x - y) <= a.window for x in (cl, ch) for y in (lo, hi)) else 0

    explained = set()
    by_var = defaultdict(lambda: {"n": 0, "lifted": 0, "pooled": 0, "one_sided": 0,
                                  "per_sample": [0] * len(call_sets), "role": "TP"})
    out = []
    for r in rows:
        v = by_var[r["variant"]]
        v["n"] += 1
        v["role"] = r["role"]
        if r["ref_contig"] is None:
            r["result"] = "unlifted"
            out.append(r)
            continue
        v["lifted"] += 1
        best = 0
        for c in pooled:
            m = match(c, r)
            if m:
                explained.add(c)
            best = max(best, m)
        for si, cs in enumerate(call_sets):
            if any(match(c, r) == 2 for c in cs):
                v["per_sample"][si] += 1
        if best == 2:
            v["pooled"] += 1
        elif best == 1:
            v["one_sided"] += 1
        r["result"] = {2: "found", 1: "one_sided", 0: "missed"}[best]
        out.append(r)

    ns = len(call_sets)
    print(f"{'variant':24s} {'role':8s} {'n':>3} {'lift':>4} {'found':>5} {'1side':>5}  per-sample")
    tot = [0, 0, 0]
    for k, v in by_var.items():
        ps = " ".join(f"S{i + 1}:{x}" for i, x in enumerate(v["per_sample"]))
        print(f"{k:24s} {v['role']:8s} {v['n']:>3} {v['lifted']:>4} {v['pooled']:>5} {v['one_sided']:>5}  {ps}")
        if v["role"] == "TP":
            tot[0] += v["lifted"]; tot[1] += v["pooled"]; tot[2] += v["one_sided"]
    fp = sorted(pooled - explained)
    print(f"\nTP recall (pooled over {ns} samples, both breakpoints within {a.window} bp): "
          f"{tot[1]}/{tot[0]} lifted ({tot[1] / max(1, tot[0]):.1%}); +{tot[2]} one-sided")
    print(f"unexplained calls (FP incl. assembly-discordance / organic artefacts): {len(fp)} of {len(pooled)}")
    if a.out_tsv:
        keys = ["id", "variant", "role", "element", "structure", "tags", "strand", "hs1_contig",
                "hs1_left", "hs1_right", "tsd", "samples", "vaf_by_sample", "ref_contig", "ref_lo",
                "ref_hi", "result"]
        with open(a.out_tsv, "w") as f:
            f.write("\t".join(keys + [f"S{i + 1}_support_R_L" for i in range(ns)]) + "\n")
            for r in out:
                sup = [f"{support[i][(r['id'], 'R')][0]}/{support[i][(r['id'], 'L')][0]}" for i in range(ns)]
                f.write("\t".join(str(r.get(k)) for k in keys) + "\t" + "\t".join(sup) + "\n")


if __name__ == "__main__":
    main()
