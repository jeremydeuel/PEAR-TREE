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

# a call header is @<contig>:<left>-<right>:...; either coordinate may carry a
# non-numeric partner prefix (polyA_<pos>, disc_<pos> from Feature A, or oneside_<pos>
# for a TPRT one-sided locus, which repeats the real coordinate), which we
# strip so the call is scored by its numeric breakpoint position.
CALL_RE = re.compile(r"^@([^:]+):(?:polyA_|disc_|oneside_)?(\d+)-(?:polyA_|disc_|oneside_)?(\d+):")


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


def split_roles(rows):
    """Extended (--types) truth carries role=TP|ARTEFACT; legacy truth is all TP."""
    tps = [r for r in rows if r.get("role", "TP") == "TP"]
    arts = [r for r in rows if r.get("role", "TP") == "ARTEFACT"]
    return tps, arts


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


def _match(c, t, window):
    """A call matches a truth row when both breakpoints agree within `window` (order-free,
    so target-site deletions / L1-mediated deletions with left > right still match)."""
    if c[0] != t["contig"]:
        return False
    cl, cr = sorted((c[1], c[2]))
    tl, tr = sorted((t["left"], t["right"]))
    return abs(cl - tl) <= window and abs(cr - tr) <= window


def score(truth, calls, window, artefacts=()):
    calls = list(calls)
    matched_truth = 0
    used = [False] * len(calls)
    by_class = {}
    by_vaf = {}
    partial = {}
    for t in truth:
        cls, vb = t.get("variant", t["class"]), vaf_bin(t["vaf"])
        if cls in (".", "") or cls.startswith("legacy:"):
            cls = t["class"]
        by_class.setdefault(cls, [0, 0]); by_vaf.setdefault(vb, [0, 0])
        by_class[cls][1] += 1; by_vaf[vb][1] += 1
        hit = None
        for i, c in enumerate(calls):
            if used[i]:
                continue
            if _match(c, t, window):
                hit = i
                break
        if hit is not None:
            used[hit] = True
            matched_truth += 1
            by_class[cls][0] += 1; by_vaf[vb][0] += 1
        else:
            # one-sided: some call has ONE breakpoint at either truth junction (e.g. a
            # poly-A rescue paired with the wrong partner, or an L1-mediated deletion)
            ends = (t["left"], t["right"])
            if any(c[0] == t["contig"] and any(abs(x - e) <= window for x in c[1:] for e in ends)
                   for c in calls):
                partial.setdefault(cls, 0)
                partial[cls] += 1
    true_pos = sum(used)
    # unmatched calls with one breakpoint at a TP junction (one-sided / mis-partnered calls
    # of a real event) vs genuinely stray calls
    tp_ends = {}
    for t in truth:
        tp_ends.setdefault(t["contig"], []).extend((t["left"], t["right"]))
    near_tp = sum(1 for i, c in enumerate(calls) if not used[i]
                  and any(abs(x - e) <= window for x in c[1:] for e in tp_ends.get(c[0], ())))
    # calls at labelled artefact loci (role=ARTEFACT): reported per artefact kind; they are
    # false positives (counted in false_pos as before)
    by_artefact = {}
    for t in artefacts:
        k = t.get("variant", t["class"])
        by_artefact.setdefault(k, [0, 0])
        by_artefact[k][1] += 1
        if any(not used[i] and _match(c, t, window) for i, c in enumerate(calls)):
            by_artefact[k][0] += 1
    total_calls = len(calls)
    recall = matched_truth / len(truth) if truth else 0.0
    precision = true_pos / total_calls if total_calls else 1.0
    return {
        "n_truth": len(truth), "n_calls": total_calls,
        "matched": matched_truth, "false_pos": total_calls - true_pos,
        "unmatched_at_tp_junction": near_tp,
        "recall": round(recall, 4), "precision": round(precision, 4),
        "by_class": {k: {"found": v[0], "total": v[1], "one_sided": partial.get(k, 0)}
                     for k, v in by_class.items()},
        "by_vaf": {k: {"found": v[0], "total": v[1]} for k, v in by_vaf.items()},
        "by_artefact": {k: {"called": v[0], "total": v[1]} for k, v in by_artefact.items()},
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True, help="BAM, or comma list of one patient's colony BAMs (calls pooled)")
    p.add_argument("--truth", required=True)
    p.add_argument("--binary", default="rust/peartree-discovery/target/release/peartree-discovery")
    p.add_argument("--config", default=None)
    p.add_argument("--window", type=int, default=3, help="breakpoint match tolerance (bp)")
    p.add_argument("--json", action="store_true", help="print JSON only")
    p.add_argument("--label", default="baseline")
    args = p.parse_args()

    truth, artefacts = split_roles(load_truth(args.truth))
    calls, stats = set(), {}
    with tempfile.TemporaryDirectory() as tmp:
        # --bam may list several colony BAMs of one patient (comma-separated): calls are
        # pooled (union) across samples, as the patient-level combine step would see them
        for k, bam in enumerate(args.bam.split(",")):
            out_gz = os.path.join(tmp, f"out{k}.txt.gz")
            run_discovery(args.binary, bam, args.config, out_gz)
            calls |= parse_calls(out_gz)
            sp = out_gz + ".stats.json"
            if os.path.exists(sp) and not stats:
                with open(sp) as f:
                    stats = json.load(f)
    m = score(truth, calls, args.window, artefacts)
    m["label"] = args.label
    m["stats"] = stats

    if args.json:
        print(json.dumps(m))
        return
    print(f"[{args.label}]  recall {m['recall']:.3f}  precision {m['precision']:.3f}  "
          f"(matched {m['matched']}/{m['n_truth']}, calls {m['n_calls']}, FP {m['false_pos']}, of which at a TP junction {m['unmatched_at_tp_junction']})")
    for k, v in m["by_class"].items():
        extra = f"  (+{v['one_sided']} one-sided)" if v.get("one_sided") else ""
        print(f"    class {k:4s}: {v['found']}/{v['total']}{extra}")
    for k, v in m["by_artefact"].items():
        print(f"    artefact {k}: called {v['called']}/{v['total']}")
    for k in ("<0.10", "0.10-0.25", ">=0.25"):
        if k in m["by_vaf"]:
            v = m["by_vaf"][k]
            print(f"    vaf {k:9s}: {v['found']}/{v['total']}")


if __name__ == "__main__":
    main()
