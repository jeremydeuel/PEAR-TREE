#!/usr/bin/env python3
"""Calibrate the TPRT point system against simulator truth.

    python -m tools.rte.calibrate --truth truth.tsv --annot P1.annotated.csv.gz [--json out.json]

truth.tsv (simulator): tab-separated with a header; columns
    insertion_id (or locus / id)   must match annotate's `locus` column
    element, structure, tags       expected vocabulary (plans/tprt_hallmarks/SPEC.md)
    role                           TP | ARTEFACT
annot: the annotate_v2 table with the tools/rte columns (element ... tprt_call).

Reports
  * matching / missing loci
  * element and structure confusion (truth vs called) and per-tag recall / precision
    (tag names compared before '=', e.g. TD3P_SOURCE)
  * per-feature separation: for every `feature` that appears in tprt_points, the fraction of TP
    and ARTEFACT loci carrying it, the mean contribution, and the empirical log-odds
    ln(P(f|TP)/P(f|ART)) -- a data-driven weight suggestion (Laplace-smoothed)
  * ROC of tprt_score (AUC by Mann-Whitney) + TPR/FPR at every score cut and at the current
    call thresholds
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import Counter, defaultdict


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_table(path, keys=("insertion_id", "locus", "id")):
    with _open(path) as fh:
        rows = list(csv.DictReader((l for l in fh if l.strip() and not l.startswith("##")), delimiter="\t"))
    out = {}
    for r in rows:
        k = next((r[c] for c in keys if c in r and r[c]), None)
        if k:
            out[k] = r
    return out


def tagset(s):
    if not s or s == ".":
        return set()
    return {t.split("=", 1)[0] for t in s.split(",") if t}


def parse_points(s):
    out = {}
    if not s or s == ".":
        return out
    for p in s.split(";"):
        if ":" in p:
            k, v = p.rsplit(":", 1)
            try:
                out[k] = float(v)
            except ValueError:
                pass
    return out


def auc(pos, neg):
    """Mann-Whitney AUC (ties count 1/2)."""
    if not pos or not neg:
        return float("nan")
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    rank = 0.0
    i = 0
    r_pos = 0.0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        avg = (i + 1 + j) / 2.0
        r_pos += avg * sum(1 for k in range(i, j) if allv[k][1] == 1)
        i = j
    n1, n0 = len(pos), len(neg)
    return (r_pos - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def calibrate(truth, annot, thresholds=None):
    from .score import THRESHOLDS
    th = dict(THRESHOLDS)
    th.update(thresholds or {})
    rep = {"n_truth": len(truth), "n_annot": len(annot)}
    keys = [k for k in truth if k in annot]
    rep["n_matched"] = len(keys)
    rep["missing"] = [k for k in truth if k not in annot][:50]
    conf_el = Counter()
    conf_st = Counter()
    tag_tp = Counter()
    tag_fn = Counter()
    tag_fp = Counter()
    feat = defaultdict(lambda: {"TP": 0, "ARTEFACT": 0, "sum": 0.0})
    scores = {"TP": [], "ARTEFACT": []}
    calls = Counter()
    for k in keys:
        t, a = truth[k], annot[k]
        role = (t.get("role") or "").upper()
        if t.get("element"):
            conf_el[(t["element"], a.get("element", "."))] += 1
        if t.get("structure"):
            conf_st[(t["structure"], a.get("structure", "."))] += 1
        tt, at = tagset(t.get("tags")), tagset(a.get("tags"))
        for x in tt & at:
            tag_tp[x] += 1
        for x in tt - at:
            tag_fn[x] += 1
        for x in at - tt:
            tag_fp[x] += 1
        if role not in scores:
            continue
        try:
            sc = float(a.get("tprt_score", "nan"))
        except ValueError:
            sc = float("nan")
        if not math.isnan(sc):
            scores[role].append(sc)
        calls[(role, a.get("tprt_call", "."))] += 1
        for f, v in parse_points(a.get("tprt_points")).items():
            feat[f][role] += 1
            feat[f]["sum"] += v
    n_tp, n_art = len(scores["TP"]), len(scores["ARTEFACT"])
    rep["n_tp"], rep["n_artefact"] = n_tp, n_art
    rep["element_accuracy"] = (sum(v for (t, a), v in conf_el.items() if t == a) / max(1, sum(conf_el.values())))
    rep["structure_accuracy"] = (sum(v for (t, a), v in conf_st.items() if t == a) / max(1, sum(conf_st.values())))
    rep["element_confusion"] = {f"{t}->{a}": v for (t, a), v in sorted(conf_el.items())}
    rep["structure_confusion"] = {f"{t}->{a}": v for (t, a), v in sorted(conf_st.items())}
    rep["tags"] = {x: {"recall": tag_tp[x] / max(1, tag_tp[x] + tag_fn[x]),
                       "precision": tag_tp[x] / max(1, tag_tp[x] + tag_fp[x]),
                       "tp": tag_tp[x], "fn": tag_fn[x], "fp": tag_fp[x]}
                   for x in sorted(set(tag_tp) | set(tag_fn) | set(tag_fp))}
    feats = {}
    for f, d in feat.items():
        p_tp = (d["TP"] + 1) / (n_tp + 2)
        p_art = (d["ARTEFACT"] + 1) / (n_art + 2)
        feats[f] = {"frac_TP": d["TP"] / max(1, n_tp), "frac_ARTEFACT": d["ARTEFACT"] / max(1, n_art),
                    "mean_points": d["sum"] / max(1, d["TP"] + d["ARTEFACT"]),
                    "log_odds": math.log(p_tp / p_art)}
    rep["features"] = dict(sorted(feats.items(), key=lambda x: -abs(x[1]["log_odds"])))
    rep["auc"] = auc(scores["TP"], scores["ARTEFACT"])
    cuts = sorted(set(scores["TP"] + scores["ARTEFACT"]))
    roc = []
    for c in cuts:
        tpr = sum(1 for v in scores["TP"] if v >= c) / max(1, n_tp)
        fpr = sum(1 for v in scores["ARTEFACT"] if v >= c) / max(1, n_art)
        roc.append({"cut": c, "tpr": tpr, "fpr": fpr})
    rep["roc"] = roc
    rep["thresholds"] = {name: {"cut": c,
                                "tpr": sum(1 for v in scores["TP"] if v >= c) / max(1, n_tp),
                                "fpr": sum(1 for v in scores["ARTEFACT"] if v >= c) / max(1, n_art)}
                         for name, c in th.items()}
    rep["calls"] = {f"{r}:{c}": v for (r, c), v in sorted(calls.items())}
    return rep


def format_report(rep):
    L = []
    L.append(f"matched {rep['n_matched']}/{rep['n_truth']} truth loci "
             f"({rep['n_tp']} TP, {rep['n_artefact']} ARTEFACT); annot rows {rep['n_annot']}")
    if rep["missing"]:
        L.append(f"missing (first {len(rep['missing'])}): {', '.join(rep['missing'][:10])}")
    L.append(f"element accuracy {rep['element_accuracy']:.3f}   structure accuracy {rep['structure_accuracy']:.3f}")
    for name in ("element_confusion", "structure_confusion"):
        off = {k: v for k, v in rep[name].items() if k.split("->")[0] != k.split("->")[1]}
        if off:
            L.append(f"  {name} (off-diagonal): " + ", ".join(f"{k}={v}" for k, v in off.items()))
    if rep["tags"]:
        L.append("tag            recall precision   tp   fn   fp")
        for t, d in rep["tags"].items():
            L.append(f"  {t:<14}{d['recall']:6.2f} {d['precision']:9.2f} {d['tp']:4d} {d['fn']:4d} {d['fp']:4d}")
    L.append("feature                 frac_TP frac_ART mean_pts log_odds")
    for f, d in rep["features"].items():
        L.append(f"  {f:<22}{d['frac_TP']:7.2f} {d['frac_ARTEFACT']:8.2f} {d['mean_points']:8.2f} {d['log_odds']:8.2f}")
    L.append(f"ROC AUC of tprt_score: {rep['auc']:.4f}")
    for name, d in rep["thresholds"].items():
        L.append(f"  call {name:<12} score >= {d['cut']:<5g} TPR {d['tpr']:.3f}  FPR {d['fpr']:.3f}")
    L.append("calls by role: " + ", ".join(f"{k}={v}" for k, v in rep["calls"].items()))
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--json", help="also write the full report (incl. ROC points) as JSON")
    a = ap.parse_args(argv)
    rep = calibrate(read_table(a.truth), read_table(a.annot))
    print(format_report(rep))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
