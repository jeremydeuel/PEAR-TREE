#!/usr/bin/env python3
"""VAL-1 toggle-matrix driver (Phase 3).

Runs the scorer over a matrix of configs and tabulates recall / precision so a
behaviour toggle can be judged against truth before it is ever defaulted on.

The matrix is `baseline` (no config) plus every `*.txt` config in --configs-dir.
Put one file per toggle for the one-at-a-time runs, and a `combined.txt` enabling
several at once for the combined-ON run the plan requires. As later phases add
config keys (SPEC-3/4 coverage gates, SENS-1/2/7 relaxations, ...), drop a config
file per new toggle into that directory — no code change needed here.

The `combined` column is the honest one for a relaxation: a toggle that lifts
recall in isolation must not wreck precision once its gate is also on.
"""
import argparse
import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run_one(bam, truth, binary, config, label, window):
    cmd = [sys.executable, os.path.join(HERE, "score.py"),
           "--bam", bam, "--truth", truth, "--binary", binary,
           "--window", str(window), "--json", "--label", label]
    if config:
        cmd += ["--config", config]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    return json.loads(out)


def low_vaf_recall(m):
    v = m["by_vaf"].get("<0.10")
    return f"{v['found']}/{v['total']}" if v else "-"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--binary", default="rust/peartree-discovery/target/release/peartree-discovery")
    p.add_argument("--configs-dir", default=os.path.join(HERE, "configs"))
    p.add_argument("--window", type=int, default=3)
    args = p.parse_args()

    runs = [("baseline", None)]
    if os.path.isdir(args.configs_dir):
        for cfg in sorted(glob.glob(os.path.join(args.configs_dir, "*.txt"))):
            runs.append((os.path.splitext(os.path.basename(cfg))[0], cfg))

    results = [run_one(args.bam, args.truth, args.binary, cfg, label, args.window)
               for label, cfg in runs]

    w = max(len(r["label"]) for r in results)
    print(f"{'config'.ljust(w)}  recall  prec.   matched  FP   lowVAF")
    print("-" * (w + 40))
    for r in results:
        print(f"{r['label'].ljust(w)}  {r['recall']:.3f}   {r['precision']:.3f}  "
              f"{r['matched']:>3}/{r['n_truth']:<3}  {r['false_pos']:>3}  {low_vaf_recall(r)}")

    base = results[0]
    print(f"\nbaseline: {base['matched']}/{base['n_truth']} truth insertions, "
          f"{base['false_pos']} false positives.")
    print("NOTE: synthetic truth only — necessary, not sufficient. Real EN/ERV/twin-priming "
          "artefacts and low-VAF specificity still need orthogonal (IGV/long-read/PCR) confirmation.")


if __name__ == "__main__":
    main()
