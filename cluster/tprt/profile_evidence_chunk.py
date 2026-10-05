#!/usr/bin/env python3
"""Profile combine's per-junction evidence evaluation on one evidence shard (rows.<c>.tsv.gz under
<stem>.evidence_shards/, written by apply_evidence) -- which part of evaluate_junction costs the time.

    cd <checkout> && python cluster/tprt/profile_evidence_chunk.py \
        --shard <stem>.evidence_shards/rows.0.tsv.gz --sidecar <one discovery>.txt.gz.evidence.tsv.gz [--n 300]

Evaluates the first --n (member locus, side) groups of the shard exactly as combine does for a
single-member insertion (no fuzzy-merge pooling), under cProfile; prints the top functions by
cumulative and by own time, plus the per-junction time distribution and the reads per junction.
Uses src/config.py's combine_insertions settings (run from the checkout whose config you test).
"""
import argparse
import cProfile
import gzip
import os
import pstats
import statistics
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", required=True)
    ap.add_argument("--sidecar", required=True, help="any discovery evidence sidecar (for the column header)")
    ap.add_argument("--n", type=int, default=300)
    args = ap.parse_args()
    from config import CONFIG
    import combine_insertions_evidence as ev
    cfg = CONFIG["combine_insertions"]
    with gzip.open(args.sidecar, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    groups = defaultdict(list)
    order = []
    with gzip.open(args.shard, "rt") as fh:
        for line in fh:
            fb, rest = line.split("\t", 1)
            p = rest.rstrip("\n").split("\t")
            if len(p) != len(header):
                continue
            r = ev.EvidenceRow(ev.sample_name(fb), dict(zip(header, p)))
            key = (fb, r.locus, r.side)
            if key not in groups:
                if len(order) >= args.n:
                    continue
                order.append(key)
            groups[key].append(r)
    times, sizes = [], []
    prof = cProfile.Profile()
    for key in order:
        rows = groups[key]
        t = time.perf_counter()
        prof.enable()
        ev.evaluate_junction(key[1], key[2], rows, cfg, None)
        prof.disable()
        times.append(time.perf_counter() - t)
        sizes.append(len(rows))
    print(f"{len(order)} junctions, {sum(sizes)} rows: {sum(times):.1f} s total, "
          f"median {statistics.median(times) * 1000:.0f} ms, max {max(times):.2f} s; "
          f"rows/junction median {statistics.median(sizes)}, max {max(sizes)}")
    slow = sorted(zip(times, sizes, order), reverse=True)[:5]
    for t, n, key in slow:
        print(f"  slowest: {t:.2f} s  {n} rows  {key[1]} {key[2]}")
    st = pstats.Stats(prof)
    st.sort_stats("cumulative").print_stats(25)
    st.sort_stats("tottime").print_stats(15)


if __name__ == "__main__":
    main()
