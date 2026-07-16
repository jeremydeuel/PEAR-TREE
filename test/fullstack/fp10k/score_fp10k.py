#!/usr/bin/env python3
"""Score the FP-stress harness: TRUE-insertion recall (role=TP) and the EMERGENT false
positive count (calls landing away from every true insertion), with a breakdown of how
many FPs are explained by a planted non-MEI decoy (role=FP) vs. unlabelled.

Usage: score_fp10k.py <truth_hg38.tsv> <label>=<calls.txt.gz> [...] [--window N]
"""
import sys, gzip, re
from collections import defaultdict

CALL_RE = re.compile(r"^@([^:]+):(?:polyA_|disc_)?(\d+)-(?:polyA_|disc_)?(\d+):")


def load_truth(path):
    tp, fp = [], []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if r.get("status") != "scoreable":
                continue
            r["L"] = int(r["hg38_left"]); r["R"] = int(r["hg38_right"])
            (tp if r.get("role") == "TP" else fp).append(r)
    return tp, fp


def parse_calls(path):
    calls = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            m = CALL_RE.match(line)
            if m:
                calls.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return sorted(set(calls))


def index(rows):
    by = defaultdict(list)
    for i, t in enumerate(rows):
        by[t["hg38_contig"]].append(i)
    return by


def nearest(rows, by, c, cl, cr):
    return min((min(abs(cl - rows[ti]["L"]), abs(cl - rows[ti]["R"]),
                    abs(cr - rows[ti]["L"]), abs(cr - rows[ti]["R"]))
                for ti in by.get(c, ())), default=10**9)


def score(tp, fp, calls, window):
    tby, fby = index(tp), index(fp)
    matched = [False] * len(tp)
    # recall: 1-to-1 greedy match of calls to true insertions
    for (c, cl, cr) in calls:
        best, bd = None, window + 1
        for ti in tby.get(c, ()):
            if matched[ti]:
                continue
            d = min(abs(cl - tp[ti]["L"]), abs(cl - tp[ti]["R"]),
                    abs(cr - tp[ti]["L"]), abs(cr - tp[ti]["R"]))
            if d <= window and d < bd:
                best, bd = ti, d
        if best is not None:
            matched[best] = True
    by_cls = defaultdict(lambda: [0, 0])
    for ti, t in enumerate(tp):
        by_cls[t.get("class", "?")][1] += 1
        by_cls[t.get("class", "?")][0] += matched[ti]
    # emergent FP: call whose nearest true insertion is beyond the window
    emergent = decoy_hit = 0
    for (c, cl, cr) in calls:
        if nearest(tp, tby, c, cl, cr) > window:
            emergent += 1
            if nearest(fp, fby, c, cl, cr) <= window:
                decoy_hit += 1
    return {
        "n_tp": len(tp), "n_fp_decoy": len(fp), "n_calls": len(calls),
        "recovered": sum(matched), "recall": sum(matched) / len(tp) if tp else 0.0,
        "emergent_fp": emergent, "decoy_explained": decoy_hit,
        "by_cls": dict(by_cls),
    }


def main():
    args = list(sys.argv[1:])
    window = 50
    if "--window" in args:
        i = args.index("--window"); window = int(args[i + 1]); del args[i:i + 2]
    tp, fp = load_truth(args[0])
    print(f"true insertions (scoreable): {len(tp)}   planted decoys (scoreable): {len(fp)}"
          f"   window +/-{window} bp\n")
    for spec in args[1:]:
        label, path = spec.split("=", 1)
        calls = parse_calls(path)
        m = score(tp, fp, calls, window)
        print(f"== {label} ==")
        print(f"  TRUE recall  {m['recall']*100:6.2f}%   "
              f"({m['recovered']}/{m['n_tp']} true insertions recovered)")
        for k in sorted(m["by_cls"]):
            r, t = m["by_cls"][k]
            print(f"      {k:14s} {r:5d}/{t:<5d} {100*r/t if t else 0:5.1f}%")
        print(f"  EMERGENT FALSE POSITIVES: {m['emergent_fp']}   "
              f"(of {m['n_calls']} total calls; "
              f"{m['decoy_explained']} land on a planted decoy, "
              f"{m['emergent_fp']-m['decoy_explained']} unlabelled)")
        print()


if __name__ == "__main__":
    main()
