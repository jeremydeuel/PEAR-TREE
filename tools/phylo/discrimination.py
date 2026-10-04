#!/usr/bin/env python3
"""Discriminatory power of the caller's scores, with the phylogeny as the truth proxy.

  python tools/phylo/discrimination.py --fit <tree_fit out>/phylo_fit.tsv --out <dir> [--truth truth.tsv]

Among informative shared loci (>= 2 confident carriers, >= 1 confident wild-type), a locus
labelled `phylo_consistent` (log10 BF tree vs non-tree >= threshold) is the TP proxy and
`phylo_violating` the FP proxy (ambiguous loci are left out). `--include-noise` also counts
loci classed `noise` (constant sub-clonal allele fraction in every colony) as FP proxies.

Scored: `tprt_score`, `tprt_score_no_multicolony` (multi_colony is +1 on every shared locus by
construction, so it is removed for a non-circular score), the `tprt_call` ordinal, every point
feature parsed from `tprt_points` (`name:+3;name:-1`; absent = 0), numeric annotate columns,
and the legacy annotate `class` / `passes_combine_genotypes` as categorical filters. AUC
(Mann-Whitney) and average precision with bootstrap 95 % CIs; ROC / PR / score-distribution
PNGs; private loci (single colony, unverifiable by the tree) per score band.
With `--truth` (simulation), the proxy labels are also scored against the truth.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

CALL_ORDER = {"ARTEFACT_LIKE": 0, "UNCERTAIN": 1, "LIKELY_TPRT": 2, "TPRT": 3}
NUMERIC = ("polya_len", "tsd_len", "element_identity", "en_mismatches", "covered_5p", "covered_3p")
CIRCULAR = {"multi_colony"}


# ------------------------------------------------------------------ metrics


def auc(y, s):
    """Mann-Whitney AUC (ties count 1/2); NaN if a class is empty."""
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    m = np.isfinite(s)
    y, s = y[m], s[m]
    n1, n0 = y.sum(), (~y).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(s).rank().values
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def average_precision(y, s):
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    m = np.isfinite(s)
    y, s = y[m], s[m]
    if y.sum() == 0:
        return np.nan
    o = np.argsort(-s, kind="mergesort")
    y, s = y[o], s[o]
    # step-wise AP over distinct thresholds
    tp = np.cumsum(y)
    k = np.arange(1, len(y) + 1)
    last = np.r_[s[1:] != s[:-1], True]
    prec = tp[last] / k[last]
    rec = tp[last] / y.sum()
    return float(np.sum(np.diff(np.r_[0, rec]) * prec))


def roc_points(y, s):
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    m = np.isfinite(s)
    y, s = y[m], s[m]
    th = np.unique(s)[::-1]
    tpr = [(s[y] >= t).mean() for t in th]
    fpr = [(s[~y] >= t).mean() for t in th]
    return np.r_[0, fpr, 1], np.r_[0, tpr, 1]


def pr_points(y, s):
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    m = np.isfinite(s)
    y, s = y[m], s[m]
    th = np.unique(s)[::-1]
    prec = [y[s >= t].mean() for t in th]
    rec = [(s[y] >= t).mean() for t in th]
    return np.array(rec), np.array(prec)


def bootstrap_ci(y, s, fn, B=1000, seed=1):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    vals = []
    for _ in range(B):
        i = rng.integers(0, len(y), len(y))
        v = fn(y[i], s[i])
        if np.isfinite(v):
            vals.append(v)
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


# ------------------------------------------------------------------ features


def parse_points(s):
    out = {}
    if not isinstance(s, str) or s in (".", ""):
        return out
    for kv in s.split(";"):
        if ":" in kv:
            k, v = kv.rsplit(":", 1)
            try:
                out[k.strip()] = out.get(k.strip(), 0.0) + float(v)
            except ValueError:
                pass
    return out


def build_features(df):
    feats = pd.DataFrame(index=df.index)
    if "tprt_score" in df.columns:
        feats["tprt_score"] = pd.to_numeric(df["tprt_score"], errors="coerce")
    pts = df["tprt_points"].map(parse_points) if "tprt_points" in df.columns else pd.Series([{}] * len(df), index=df.index)
    names = sorted({k for d in pts for k in d})
    for k in names:
        feats[f"pt:{k}"] = [d.get(k, 0.0) for d in pts]
    if "tprt_score" in feats and "pt:multi_colony" in feats:
        feats["tprt_score_no_multicolony"] = feats["tprt_score"] - feats["pt:multi_colony"]
    if "tprt_call" in df.columns:
        feats["tprt_call_ordinal"] = df["tprt_call"].map(CALL_ORDER).astype(float)
    for c in NUMERIC:
        if c in df.columns:
            feats[c] = pd.to_numeric(df[c].replace(".", np.nan), errors="coerce")
    if "passes_combine_genotypes" in df.columns:
        feats["passes_combine_genotypes"] = df["passes_combine_genotypes"].astype(str).str.lower().eq("true").astype(float)
    return feats


# ------------------------------------------------------------------ report


def md_table(rows, cols):
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        out.append("| " + " | ".join(f"{v:.3f}" if isinstance(v, float) and np.isfinite(v) else
                                     ("NA" if isinstance(v, float) else str(v)) for v in r) + " |")
    return "\n".join(out)


def plots(out, y, feats, df, top):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    files = []
    pal = ["#3b6ea8", "#c2553a", "#4f9a5c", "#8a5fb0", "#c99a2e", "#5aa3b5", "#9a9a9a"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for k, f in enumerate(top):
        s = feats.loc[y.index, f].values
        fpr, tpr = roc_points(y.values, s)
        ax[0].plot(fpr, tpr, color=pal[k % len(pal)], label=f"{f} (AUC {auc(y.values, s):.2f})")
        rec, prec = pr_points(y.values, s)
        ax[1].plot(rec, prec, color=pal[k % len(pal)], label=f)
    ax[0].plot([0, 1], [0, 1], "--", color="#bbb")
    ax[0].set_xlabel("FPR (phylo_violating called)"); ax[0].set_ylabel("TPR (phylo_consistent called)")
    ax[0].legend(fontsize=7); ax[0].set_title("ROC, informative shared loci")
    ax[1].axhline(y.mean(), ls="--", color="#bbb")
    ax[1].set_xlabel("recall"); ax[1].set_ylabel("precision"); ax[1].set_title("precision-recall")
    fig.tight_layout()
    f = os.path.join(out, "roc_pr.png")
    fig.savefig(f, dpi=110); plt.close(fig); files.append(f)
    if "tprt_score" in feats:
        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        groups = [("phylo_consistent", (df["class"] == "informative_shared") & (df["label"] == "phylo_consistent")),
                  ("phylo_violating", (df["class"] == "informative_shared") & (df["label"] == "phylo_violating")),
                  ("private", df["class"] == "private"), ("germline", df["class"] == "germline"),
                  ("noise", df["class"] == "noise")]
        vals = feats["tprt_score"]
        bins = np.arange(np.floor(np.nanmin(vals)) if vals.notna().any() else -10,
                         (np.ceil(np.nanmax(vals)) if vals.notna().any() else 20) + 1.5, 1.0)
        for k, (g, m) in enumerate(groups):
            v = vals[m].dropna()
            if len(v):
                ax.hist(v, bins=bins, histtype="step", density=True, lw=1.6, color=pal[k], label=f"{g} (n={len(v)})")
        ax.set_xlabel("tprt_score"); ax.set_ylabel("density"); ax.legend(fontsize=8)
        ax.set_title("score distribution by phylogenetic class")
        fig.tight_layout()
        f = os.path.join(out, "score_by_class.png")
        fig.savefig(f, dpi=110); plt.close(fig); files.append(f)
    return files


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", required=True, help="phylo_fit.tsv from tree_fit.py")
    ap.add_argument("--out", required=True)
    ap.add_argument("--truth", help="simulation truth (tools/phylo/calibration.py format)")
    ap.add_argument("--include-noise", action="store_true", help="loci classed noise count as FP proxies")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--min-class-n", type=int, default=5, help="per element/structure: min loci per label")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    df = pd.read_csv(a.fit, sep="\t", dtype={"locus": str}, low_memory=False)
    feats = build_features(df)
    sel = (df["class"] == "informative_shared") & df["label"].isin(["phylo_consistent", "phylo_violating"])
    if a.include_noise:
        sel |= df["class"] == "noise"
    y = (df.loc[sel, "label"] == "phylo_consistent") & (df.loc[sel, "class"] != "noise")
    L = ["# Discrimination report (phylogeny as truth proxy)", "",
         f"Input `{a.fit}`: {len(df)} loci. Validated set: {int(sel.sum())} loci "
         f"({int(y.sum())} phylo_consistent = TP proxy, {int((~y).sum())} phylo_violating"
         f"{' + noise' if a.include_noise else ''} = FP proxy); ambiguous informative loci left out: "
         f"{int(((df['class'] == 'informative_shared') & (df['label'] == 'ambiguous')).sum())}.", "",
         "Classes: " + ", ".join(f"{k} {v}" for k, v in df["class"].value_counts().items()), ""]
    rows = []
    if y.sum() and (~y).sum():
        for f in feats.columns:
            s = feats.loc[sel, f].values.astype(float)
            if np.isfinite(s).sum() < 4 or np.nanstd(s) == 0:
                continue
            A = auc(y.values, s)
            lo, hi = bootstrap_ci(y.values, s, auc, a.bootstrap)
            ap_ = average_precision(y.values, s)
            plo, phi = bootstrap_ci(y.values, s, average_precision, a.bootstrap)
            nz = int((np.nan_to_num(s) != 0).sum()) if f.startswith("pt:") else int(np.isfinite(s).sum())
            rows.append((f, nz, A, lo, hi, ap_, plo, phi))
        rows.sort(key=lambda r: -abs(r[2] - 0.5) if np.isfinite(r[2]) else 0)
        L += ["## Per-feature discrimination (AUC > 0.5: higher in phylo_consistent loci)", "",
              f"AP baseline (prevalence of consistent): {y.mean():.3f}. CIs: {a.bootstrap} bootstrap resamples of the "
              "validated loci. `n` = loci where the feature is non-zero (points) / defined (numeric).", "",
              md_table([(r[0] + (" (circular)" if r[0].replace("pt:", "") in CIRCULAR else ""),) + tuple(r[1:]) for r in rows],
                       ["feature", "n", "AUC", "AUC lo", "AUC hi", "AP", "AP lo", "AP hi"]), ""]
        # per element / structure
        for col in ("element", "structure", "locus_kind"):
            if col not in df.columns or "tprt_score" not in feats:
                continue
            sub = []
            for v, g in df.loc[sel].groupby(col):
                yy = y.loc[g.index]
                if yy.sum() >= a.min_class_n and (~yy).sum() >= a.min_class_n:
                    s = feats.loc[g.index, "tprt_score"].values
                    lo, hi = bootstrap_ci(yy.values, s, auc, a.bootstrap)
                    sub.append((v, int(yy.sum()), int((~yy).sum()), auc(yy.values, s), lo, hi))
                else:
                    sub.append((v, int(yy.sum()), int((~yy).sum()), np.nan, np.nan, np.nan))
            L += [f"## tprt_score AUC per `{col}`", "", md_table(sub, [col, "consistent", "violating", "AUC", "lo", "hi"]), ""]
        # legacy categorical filters
        for col in ("class_annot", "class", "conclusion", "tprt_call"):
            src = "ann_class" if col == "class_annot" else col
            if src not in df.columns or src == "class":
                continue
            sub = []
            for v, g in df.loc[sel].groupby(src):
                yy = y.loc[g.index]
                sub.append((str(v)[:60], len(g), int(yy.sum()), int((~yy).sum()), float(yy.mean())))
            sub.sort(key=lambda r: -r[1])
            L += [f"## Legacy / categorical `{src}` (precision = consistent / all)", "",
                  md_table(sub[:25], [src, "n", "consistent", "violating", "precision"]), ""]
        top = [f for f in ("tprt_score", "tprt_score_no_multicolony", "tprt_call_ordinal") if f in feats]
        top += [r[0] for r in rows if r[0].startswith("pt:") and r[0] not in top][:4]
        files = plots(a.out, y, feats, df, top)
        L += [f"![{os.path.basename(f)}]({os.path.basename(f)})" for f in files] + [""]
    else:
        L += ["Too few validated loci in one of the two proxy classes for ROC / AUC.", ""]
    # private loci per score band
    if "tprt_call" in df.columns:
        grp = {"private": df["class"] == "private",
               "consistent": (df["class"] == "informative_shared") & (df["label"] == "phylo_consistent"),
               "violating": (df["class"] == "informative_shared") & (df["label"] == "phylo_violating"),
               "germline": df["class"] == "germline", "noise": df["class"] == "noise"}
        sub = []
        for band in ["TPRT", "LIKELY_TPRT", "UNCERTAIN", "ARTEFACT_LIKE"]:
            sub.append((band,) + tuple(int(((df["tprt_call"] == band) & m).sum()) for m in grp.values()))
        L += ["## Loci per score band (private loci cannot be validated by the tree; compare their band mix "
              "to the validated sets)", "", md_table(sub, ["tprt_call"] + list(grp)), ""]
        if "tprt_score" in feats:
            qs = []
            for g, m in grp.items():
                v = feats.loc[m, "tprt_score"].dropna()
                qs.append((g, len(v), float(v.median()) if len(v) else np.nan,
                           float(v.quantile(0.25)) if len(v) else np.nan, float(v.quantile(0.75)) if len(v) else np.nan))
            L += [md_table(qs, ["set", "n", "median score", "q25", "q75"]), ""]
    # proxy vs truth (simulation)
    if a.truth:
        tr = pd.read_csv(a.truth, sep="\t", dtype=str).drop_duplicates("locus").set_index("locus")
        df["truth"] = df["locus"].map(tr["truth"]).fillna("none")
        tab = pd.crosstab(df.loc[df["class"] == "informative_shared", "truth"],
                          df.loc[df["class"] == "informative_shared", "label"])
        L += ["## Proxy labels vs simulation truth (informative shared loci)", "",
              md_table([(i,) + tuple(int(v) for v in r) for i, r in tab.iterrows()], ["truth"] + list(tab.columns)), ""]
        tp = df["truth"].isin(["clade", "private", "germline"])
        fp = df["truth"].isin(["nonclade", "artefact"])
        m = sel & (tp | fp)
        if "tprt_score" in feats and m.sum():
            yt = tp[m].values
            L += [f"tprt_score AUC against TRUTH on the validated loci with an event (n={int(m.sum())}): "
                  f"{auc(yt, feats.loc[m, 'tprt_score'].values):.3f}; against the PROXY: "
                  f"{auc(y.loc[m[m].index].values, feats.loc[m, 'tprt_score'].values):.3f}.", ""]
    with open(os.path.join(a.out, "report.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    pd.DataFrame(rows, columns=["feature", "n", "auc", "auc_lo", "auc_hi", "ap", "ap_lo", "ap_hi"]).to_csv(
        os.path.join(a.out, "features.tsv"), sep="\t", index=False)
    print("\n".join(L[:40]))


if __name__ == "__main__":
    main()
