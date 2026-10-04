#!/usr/bin/env python3
"""Calibration of a tree_fit run against simulation truth.

  python tools/phylo/calibration.py --fit <tree_fit out> --truth truth.tsv --out <dir> [--label NAME]

truth.tsv: `locus  truth  carriers` with truth in {clade, private, germline, nonclade, noise,
none} (none = a called locus that is no simulated event, e.g. an hs1-vs-GRCh38 difference) and
carriers a comma list of the colonies that carry it. Loci of the fit absent from truth.tsv are
`none`.

Reports (report.md + PNGs):
  * per-colony violation calibration: among TRUE carriers, fraction with p_le_present <= alpha
    (a missing-leaf flag) vs alpha; among TRUE non-carriers of clade loci, fraction with
    p_ge_absent <= alpha (an extra-carrier flag) vs alpha. Valid p-values give <= alpha.
  * dropout calibration: true carriers binned by predicted P(n_alt = 0), observed zero-alt rate.
  * locus level: class x label per truth type; best clade == true clade; p_locus uniformity.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

ALPHAS = (0.001, 0.005, 0.01, 0.05, 0.1)


def load(fit_dir, truth_path):
    fit = pd.read_csv(os.path.join(fit_dir, "phylo_fit.tsv"), sep="\t", dtype={"locus": str})
    cells = pd.read_csv(os.path.join(fit_dir, "cells.tsv.gz"), sep="\t", dtype={"locus": str})
    truth = pd.read_csv(truth_path, sep="\t", dtype=str).drop_duplicates("locus").set_index("locus")
    fit["truth"] = fit["locus"].map(truth["truth"]).fillna("none")
    fit["true_carriers"] = fit["locus"].map(truth["carriers"]).fillna("")
    car = {l: set(str(c).split(",")) - {""} for l, c in zip(fit["locus"], fit["true_carriers"])}
    cells["truth"] = cells["locus"].map(dict(zip(fit["locus"], fit["truth"])))
    cells["true_carrier"] = [c in car.get(l, ()) for l, c in zip(cells["locus"], cells["colony"])]
    return fit, cells


def violation_calibration(cells):
    som = cells[cells["truth"].isin(["clade", "private", "nonclade"])]
    car = som[som["true_carrier"]]
    non = som[~som["true_carrier"]]
    rows = []
    for a in ALPHAS:
        rows.append((a, len(car), float((car["p_le_present"] <= a).mean()) if len(car) else np.nan,
                     len(non), float((non["p_ge_absent"] <= a).mean()) if len(non) else np.nan))
    return pd.DataFrame(rows, columns=["alpha", "n_carrier_cells", "frac_missing_flag",
                                       "n_noncarrier_cells", "frac_extra_flag"])


def dropout_calibration(cells):
    car = cells[cells["truth"].isin(["clade", "private", "nonclade"]) & cells["true_carrier"]].copy()
    bins = [0, 0.001, 0.01, 0.03, 0.1, 0.2, 0.35, 0.5, 1.0001]
    car["bin"] = pd.cut(car["p_dropout"], bins, right=False)
    g = car.groupby("bin", observed=True)
    out = pd.DataFrame({"n": g.size(), "mean_predicted": g["p_dropout"].mean(),
                        "observed_zero_alt": g["n_alt"].apply(lambda x: float((x == 0).mean()))})
    out["expected_zero"] = (g["p_dropout"].sum()).round(1)
    out["observed_zero"] = g["n_alt"].apply(lambda x: int((x == 0).sum()))
    return out.reset_index()


def locus_tables(fit):
    tab = pd.crosstab([fit["truth"], fit["class"]], fit["label"])
    sh = fit[fit["truth"] == "clade"].copy()
    sh["true_clade_found"] = [set(b.split(",")) == set(t.split(",")) for b, t in zip(sh["best_clade"], sh["true_carriers"])]
    return tab, sh


def plots(vc, dc, fit, out, label):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    files = []
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    floor = 2e-5      # log axes: a fraction of exactly 0 is drawn at the floor
    ax[0].plot(vc["alpha"], np.maximum(vc["frac_missing_flag"], floor), "o-", color="#3b6ea8",
               label="true carriers flagged missing")
    ax[0].plot(vc["alpha"], np.maximum(vc["frac_extra_flag"], floor), "s-", color="#c2553a",
               label="true non-carriers flagged extra")
    ax[0].plot([1e-3, 0.1], [1e-3, 0.1], "--", color="#888", label="nominal")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_ylim(floor / 1.5, 0.2)
    ax[0].set_xlabel("alpha"); ax[0].set_ylabel(f"fraction flagged (0 drawn at {floor:g})"); ax[0].legend(fontsize=8)
    ax[0].set_title("per-colony violation p-values")
    d = dc[dc["n"] > 0]
    ax[1].plot(np.maximum(d["mean_predicted"], floor), np.maximum(d["observed_zero_alt"], floor), "o-", color="#3b6ea8")
    for _, r in d.iterrows():
        ax[1].annotate(f"{int(r['observed_zero'])}/{int(r['n'])}", (max(r["mean_predicted"], floor),
                       max(r["observed_zero_alt"], floor)), fontsize=7, textcoords="offset points", xytext=(4, -9))
    ax[1].plot([1e-5, 0.6], [1e-5, 0.6], "--", color="#888")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("predicted P(n_alt = 0 | carrier)"); ax[1].set_ylabel("observed zero-alt rate")
    ax[1].set_title("dropout calibration (true carriers)")
    cl = fit[fit["truth"] == "clade"]["p_locus"].dropna()
    if len(cl):
        x = np.sort(cl.values)
        ax[2].plot(x, np.arange(1, len(x) + 1) / len(x), color="#3b6ea8", label=f"true clade loci (n={len(x)})")
    nc = fit[fit["truth"] == "nonclade"]["p_locus"].dropna()
    if len(nc):
        x = np.sort(nc.values)
        ax[2].plot(x, np.arange(1, len(x) + 1) / len(x), color="#c2553a", label=f"non-clade loci (n={len(x)})")
    ax[2].plot([0, 1], [0, 1], "--", color="#888")
    ax[2].set_xlabel("bootstrap p_locus"); ax[2].set_ylabel("ECDF"); ax[2].legend(fontsize=8)
    ax[2].set_title("locus-level p-value")
    fig.suptitle(label)
    fig.tight_layout()
    f = os.path.join(out, "calibration.png")
    fig.savefig(f, dpi=110)
    plt.close(fig)
    files.append(f)
    return files


def md_table(df, floatfmt="{:.4g}"):
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(floatfmt.format(v) if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    fit, cells = load(a.fit, a.truth)
    vc = violation_calibration(cells)
    dc = dropout_calibration(cells)
    tab, sh = locus_tables(fit)
    files = plots(vc, dc, fit, a.out, a.label)
    L = [f"# Phylo calibration vs truth {a.label}", "", "## Per-colony violation p-values (valid: fraction <= alpha)", "",
         md_table(vc), "", "## Dropout calibration (true carriers)", "", md_table(dc.astype({"bin": str})), "",
         "## Class x label per truth type", "", md_table(tab.reset_index()), ""]
    inf = sh[sh["class"] == "informative_shared"]
    if len(sh):
        L += [f"True clade loci (>= 2 carriers): {len(sh)}; classed informative_shared {len(inf)}; "
              f"of those best clade == truth {int(inf['true_clade_found'].sum())} "
              f"({inf['true_clade_found'].mean():.1%}), label phylo_consistent {int((inf['label'] == 'phylo_consistent').sum())}, "
              f"violating {int((inf['label'] == 'phylo_violating').sum())}; "
              f"p_locus < 0.01: {int((inf['p_locus'] < 0.01).sum())}, < 0.05: {int((inf['p_locus'] < 0.05).sum())}.", ""]
    nc = fit[(fit["truth"] == "nonclade") & (fit["class"] == "informative_shared")]
    if len(nc):
        L += [f"Non-clade loci classed informative_shared: {len(nc)}; phylo_violating {int((nc['label'] == 'phylo_violating').sum())} "
              f"({(nc['label'] == 'phylo_violating').mean():.1%}), ambiguous {int((nc['label'] == 'ambiguous').sum())}, "
              f"consistent {int((nc['label'] == 'phylo_consistent').sum())}.", ""]
    L += [f"![calibration]({os.path.basename(f)})" for f in files]
    with open(os.path.join(a.out, "report.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
