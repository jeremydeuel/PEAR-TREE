#!/usr/bin/env python3
"""VAL-1 scorer (Phase 3).

Run discovery on a simulated BAM under one config, match its calls against the
truth TSV, and report recall / precision overall, per element class, and per VAF
bin, plus the OBS-1 reject-counter deltas. Emits a human table and (optionally)
a JSON blob for the matrix driver.
"""
import argparse
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile

CALL_RE = re.compile(r"^@([^:]+):(\d+)-(\d+):")


def load_truth(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            vals = line.rstrip("\n").split("\t")
            r = dict(zip(header, vals))
            r["left"] = int(r["left"]); r["right"] = int(r["right"])
            r["vaf"] = float(r["vaf"])
            rows.append(r)
    return rows


def parse_calls(out_gz):
    calls = set()
    with gzip.open(out_gz, "rt") as f:
        for line in f:
            m = CALL_RE.match(line)
            if m:
                calls.add((m.group(1), int(m.group(2)), int(m.group(3))))
    return calls


def vaf_bin(v):
    if v < 0.1:
        return "<0.10"
    if v < 0.25:
        return "0.10-0.25"
    return ">=0.25"


def run_discovery(binary, bam, config, out_gz):
    cmd = [binary, "--step", "discover", "--bam", bam, "--out", out_gz]
    if config:
        cmd += ["--config", config]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def score(truth, calls, window):
    calls = list(calls)
    matched_truth = 0
    used = [False] * len(calls)
    by_class = {}
    by_vaf = {}
    for t in truth:
        cls, vb = t["class"], vaf_bin(t["vaf"])
        by_class.setdefault(cls, [0, 0]); by_vaf.setdefault(vb, [0, 0])
        by_class[cls][1] += 1; by_vaf[vb][1] += 1
        hit = None
        for i, c in enumerate(calls):
            if used[i]:
                continue
            if c[0] == t["contig"] and abs(c[1] - t["left"]) <= window and abs(c[2] - t["right"]) <= window:
                hit = i
                break
        if hit is not None:
            used[hit] = True
            matched_truth += 1
            by_class[cls][0] += 1; by_vaf[vb][0] += 1
    true_pos = sum(used)
    total_calls = len(calls)
    recall = matched_truth / len(truth) if truth else 0.0
    precision = true_pos / total_calls if total_calls else 1.0
    return {
        "n_truth": len(truth), "n_calls": total_calls,
        "matched": matched_truth, "false_pos": total_calls - true_pos,
        "recall": round(recall, 4), "precision": round(precision, 4),
        "by_class": {k: {"found": v[0], "total": v[1]} for k, v in by_class.items()},
        "by_vaf": {k: {"found": v[0], "total": v[1]} for k, v in by_vaf.items()},
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--binary", default="rust/peartree-discovery/target/release/peartree-discovery")
    p.add_argument("--config", default=None)
    p.add_argument("--window", type=int, default=3, help="breakpoint match tolerance (bp)")
    p.add_argument("--json", action="store_true", help="print JSON only")
    p.add_argument("--label", default="baseline")
    args = p.parse_args()

    truth = load_truth(args.truth)
    with tempfile.TemporaryDirectory() as tmp:
        out_gz = os.path.join(tmp, "out.txt.gz")
        run_discovery(args.binary, args.bam, args.config, out_gz)
        calls = parse_calls(out_gz)
        stats = {}
        sp = out_gz + ".stats.json"
        if os.path.exists(sp):
            with open(sp) as f:
                stats = json.load(f)
    m = score(truth, calls, args.window)
    m["label"] = args.label
    m["stats"] = stats

    if args.json:
        print(json.dumps(m))
        return
    print(f"[{args.label}]  recall {m['recall']:.3f}  precision {m['precision']:.3f}  "
          f"(matched {m['matched']}/{m['n_truth']}, calls {m['n_calls']}, FP {m['false_pos']})")
    for k, v in m["by_class"].items():
        print(f"    class {k:4s}: {v['found']}/{v['total']}")
    for k in ("<0.10", "0.10-0.25", ">=0.25"):
        if k in m["by_vaf"]:
            v = m["by_vaf"][k]
            print(f"    vaf {k:9s}: {v['found']}/{v['total']}")


if __name__ == "__main__":
    main()
