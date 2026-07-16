#!/usr/bin/env python3
"""Score discovery / combined calls against the chain-lifted hg38 truth.

Recall = matched scoreable truths / total scoreable truths.
FP     = calls not matching any truth (assembly-discordance false positives).
Usage: score_10k.py <truth_hg38.tsv> <label>=<calls.txt.gz> [<label>=<calls> ...] [--window N]
"""
import sys, gzip, re
from collections import defaultdict

CALL_RE = re.compile(r"^@([^:]+):(?:polyA_|disc_)?(\d+)-(?:polyA_|disc_)?(\d+):")


def load_truth(path):
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if r.get("status") != "scoreable":
                continue
            r["L"] = int(r["hg38_left"]); r["R"] = int(r["hg38_right"])
            rows.append(r)
    return rows


def parse_calls(path):
    calls = []
    with gzip.open(path, "rt") as f:
        for line in f:
            m = CALL_RE.match(line)
            if m:
                calls.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return sorted(set(calls))


def score(truth, calls, window):
    # index truth by contig for fast lookup
    tby = defaultdict(list)
    for i, t in enumerate(truth):
        tby[t["hg38_contig"]].append(i)
    matched = [False] * len(truth)
    call_hit = [False] * len(calls)
    for ci, (c, cl, cr) in enumerate(calls):
        for ti in tby.get(c, ()):
            if matched[ti]:
                continue
            t = truth[ti]
            if min(abs(cl - t["L"]), abs(cl - t["R"]),
                   abs(cr - t["L"]), abs(cr - t["R"])) <= window:
                matched[ti] = True
                call_hit[ci] = True
                break
    by_fam = defaultdict(lambda: [0, 0])
    for ti, t in enumerate(truth):
        key = f"{t['family']}_{t['variant']}"
        by_fam[key][1] += 1
        by_fam[key][0] += matched[ti]
    tp = sum(matched)
    fp = sum(1 for h in call_hit if not h)
    # genuine FP: a call whose NEAREST truth is beyond the window (reuse allowed). The
    # 1-to-1 `fp` above also counts a 2nd call landing on an already-claimed truth (a
    # duplicate call at a real insertion), which is NOT a false-positive locus.
    genuine_fp = 0
    fp_by_contig = defaultdict(int)
    for (c, cl, cr) in calls:
        d = min((min(abs(cl - truth[ti]["L"]), abs(cl - truth[ti]["R"]),
                     abs(cr - truth[ti]["L"]), abs(cr - truth[ti]["R"]))
                 for ti in tby.get(c, ())), default=10**9)
        if d > window:
            genuine_fp += 1
            fp_by_contig[c] += 1
    return {
        "n_truth": len(truth), "n_calls": len(calls), "tp": tp, "fp": fp,
        "genuine_fp": genuine_fp,
        "recall": tp / len(truth) if truth else 0.0,
        "precision": tp / len(calls) if calls else 1.0,
        "by_fam": dict(by_fam), "fp_by_contig": dict(fp_by_contig),
    }


def main():
    args = [a for a in sys.argv[1:]]
    window = 50
    if "--window" in args:
        i = args.index("--window"); window = int(args[i + 1]); del args[i:i + 2]
    truth = load_truth(args[0])
    print(f"scoreable truths: {len(truth)}  (window ±{window} bp)\n")
    for spec in args[1:]:
        label, path = spec.split("=", 1)
        calls = parse_calls(path)
        m = score(truth, calls, window)
        print(f"== {label} ==")
        print(f"  recall    {m['recall']*100:6.2f}%   ({m['tp']}/{m['n_truth']} true insertions recovered)")
        print(f"  precision {m['precision']*100:6.2f}%   ({m['tp']}/{m['n_calls']} calls)")
        print(f"  genuine false positives: {m['genuine_fp']}   "
              f"(1-to-1 scorer counts {m['fp']}; the extra are duplicate calls at true sites)")
        for k in sorted(m["by_fam"]):
            f_, t_ = m["by_fam"][k]
            print(f"      {k:26s} {f_:5d}/{t_:<5d} {100*f_/t_:5.1f}%")
        top = sorted(m["fp_by_contig"].items(), key=lambda x: -x[1])[:12]
        print(f"  FP by hg38 contig (top): " + ", ".join(f"{c}:{n}" for c, n in top))
        print()


if __name__ == "__main__":
    main()
