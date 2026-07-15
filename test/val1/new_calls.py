#!/usr/bin/env python3
"""VAL-1 orthogonal-check hook (Phase 3, requirement (b)).

Diff two discovery `.txt.gz` outputs (baseline vs a candidate toggle) and write
the calls that are NEW in the candidate as a BED + TSV. That list is what goes to
orthogonal confirmation — IGV inspection, long-read spanning, or PCR — before a
toggle that adds calls is trusted, on real data where there is no truth set.

On the synthetic set you can instead pass --truth to label each new call as a
true (matches a spike-in) or unexplained (candidate false positive) addition.
"""
import argparse
import gzip
import re

CALL_RE = re.compile(r"^@([^:]+):(\d+)-(\d+):")


def parse_calls(path):
    calls = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            m = CALL_RE.match(line)
            if m:
                calls.add((m.group(1), int(m.group(2)), int(m.group(3))))
    return calls


def load_truth(path):
    out = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            out.append((r["contig"], int(r["left"]), int(r["right"])))
    return out


def is_true(call, truth, window):
    return any(call[0] == t[0] and abs(call[1] - t[1]) <= window and abs(call[2] - t[2]) <= window
               for t in truth)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", required=True, help="baseline discovery .txt.gz")
    p.add_argument("--candidate", required=True, help="candidate (toggle) discovery .txt.gz")
    p.add_argument("--out-bed", required=True)
    p.add_argument("--truth", default=None, help="optional synthetic truth TSV to label new calls")
    p.add_argument("--window", type=int, default=3)
    args = p.parse_args()

    base = parse_calls(args.baseline)
    cand = parse_calls(args.candidate)
    new = sorted(cand - base)
    lost = sorted(base - cand)
    truth = load_truth(args.truth) if args.truth else None

    with open(args.out_bed, "w") as f:
        f.write("# new calls in candidate vs baseline — send to IGV/long-read/PCR review\n")
        for c in new:
            label = ""
            if truth is not None:
                label = "\ttrue" if is_true(c, truth, args.window) else "\tUNEXPLAINED"
            f.write(f"{c[0]}\t{c[1]}\t{c[2]}{label}\n")

    print(f"new calls (candidate not in baseline): {len(new)}  -> {args.out_bed}")
    print(f"lost calls (baseline not in candidate): {len(lost)}")
    if truth is not None and new:
        n_true = sum(is_true(c, truth, args.window) for c in new)
        print(f"of {len(new)} new: {n_true} match a spike-in, {len(new) - n_true} UNEXPLAINED (candidate FPs)")


if __name__ == "__main__":
    main()
