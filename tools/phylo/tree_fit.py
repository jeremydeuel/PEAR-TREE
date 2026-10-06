#!/usr/bin/env python3
"""Fit every genotyped locus of one patient onto its SNV phylogeny.

A true somatic insertion is inherited by exactly the colonies below one branch of the SNV tree;
an artefact's carriers scatter. Per locus and per branch b

    log L_b = sum_{c below b} log P(d_c | present) + sum_{c not below b} log P(d_c | absent)

with the read-vote likelihoods of tools/phylo/genotype_likelihood.py (so a clade member with
no alt reads at low depth costs little, at high depth a lot). The ROOT hypothesis is "carried
by every colony" (germline het / hom, or somatic before the MRCA). Alternatives that ignore
the tree: INDEP (each colony carries it independently with probability pi ~ U(0, 1), exact by
dynamic programming over the carrier count) and NOISE (no colony carries it; one locus-specific
alt-read rate eps ~ Beta(a, b) in every colony = the constant-allele-fraction artefact).

    BF_tree = sum_b prior_b L_b  /  (L_indep + L_noise) / 2         (reported as log10)

Input formats are auto-detected (tools/genotype2_io.py): per-colony files from the legacy
genotyper or from rust/peartree-genotype2 (numeric; votes = n_alt / n_ref, rows with status != ok
are missing data), and --genotypes as combine_genotypes' call matrix or the genotype2 joint
step's numeric P(carrier) matrix (then <P>.joint.tsv beside it is joined as joint_* columns and
cross-tabulated in summary.md: this script is the independent Python check of that Rust port).

Usage (the farm kit calls exactly this):
  python tools/phylo/tree_fit.py --genotypes P.genotypes.csv.gz --genotype-dir genotypes/ \
      --tree P.tree --annotation P.annotated.tsv --out phylo/ [--sex M|F|auto] \
      [--samples colonies.tsv] [--min-depth N]

Outputs in --out: phylo_fit.tsv (one row per locus + annotation columns), violations.tsv,
colony_params.tsv, cells.tsv.gz (per locus x colony quantities), summary.md.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.special import betaln, logsumexp

if __package__ in (None, ""):
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from tools import genotype2_io as GIO  # noqa: E402
from tools.phylo import genotype_likelihood as G  # noqa: E402
from tools.phylo import tree as T  # noqa: E402

LN10 = math.log(10.0)
CLASSES = ("informative_shared", "private", "germline", "noise", "uninformative_depth")
LABELS = ("phylo_consistent", "ambiguous", "phylo_violating")


@dataclass
class Options:
    min_depth: int = 1
    branch_prior: str = "length"      # length | uniform
    root_prior: float = 0.1
    floor_frac: float = 0.01
    bf_threshold: float = 1.0         # log10 BF: >= t consistent, <= -t violating
    alpha: float = 0.01               # per-colony violation flag
    carrier_lr: float = 2.0           # log10 L1/L0 >= this: confident carrier
    wt_lr: float = 1.0                # log10 L0/L1 >= this: confident wild-type
    noise_margin: float = 1.0         # log10: NOISE beats every tree hypothesis by this -> noise
    bootstrap: int = 200
    rounds: int = 3
    seed: int = 1


# ------------------------------------------------------------------ hypothesis scores


def branch_log_prior(br: T.Branches, mode: str, root_prior: float, floor_frac: float) -> np.ndarray:
    nb = len(br.ids)
    lens = br.lengths[1:].astype(float)
    if mode == "length" and (lens > 0).any():
        floor = floor_frac * lens[lens > 0].mean()
        w = np.maximum(lens, floor)
    else:
        w = np.ones(nb - 1)
    pr = np.empty(nb)
    pr[0] = root_prior
    pr[1:] = (1.0 - root_prior) * w / w.sum()
    return np.log(pr)


HOM_PRIOR = (8.0, 1.0)     # germline hom / hemizygous: one shared fraction ~ Beta(8, 1) (mean 0.89)


def root_scores(cl: G.CellLik, haploid: np.ndarray, alt=None, n=None):
    """log L_root per locus and its components (germline het, germline hom, somatic in all
    colonies), equal weights; het is impossible on haploid loci. The hom component is a shared
    alt fraction ~ Beta(8, 1) across colonies (closed form): a hom locus reads 0.85-0.98, not
    1 - eps, because some reference-configuration reads survive (mismaps, the other junction)."""
    het = cl.lhet.sum(1)
    hom = log_noise(alt, n, cl.mask, *HOM_PRIOR) if alt is not None else cl.lhom.sum(1)
    soma = cl.l1.sum(1)
    comps = np.vstack([het, hom, soma]).T
    w = np.where(haploid[:, None], np.array([[-np.inf, math.log(0.5), math.log(0.5)]]),
                 np.log(np.full((1, 3), 1.0 / 3)))
    with np.errstate(invalid="ignore"):
        tot = logsumexp(comps + w, axis=1)
    return tot, comps


def branch_scores(cl: G.CellLik, br: T.Branches, haploid: np.ndarray, alt=None, n=None):
    D = cl.l1 - cl.l0
    S = np.empty((D.shape[0], len(br.ids)))
    S[:, 1:] = D @ br.mask[1:].T.astype(float) + cl.l0.sum(1)[:, None]
    S[:, 0], comps = root_scores(cl, haploid, alt, n)
    return S, comps


def log_indep(cl: G.CellLik) -> np.ndarray:
    """log int_0^1 prod_c [pi P(d_c|present) + (1-pi) P(d_c|absent)] dpi, exactly: expand the
    product in the Bernstein basis pi^k (1-pi)^(C-k) by DP, integrate each term to B(k+1, C-k+1)."""
    L, C = cl.l1.shape
    E = np.full((L, C + 1), -np.inf)
    E[:, 0] = 0.0
    for c in range(C):
        new = np.full_like(E, -np.inf)
        new[:, 0] = E[:, 0] + cl.l0[:, c]
        new[:, 1:] = np.logaddexp(E[:, 1:] + cl.l0[:, c:c + 1], E[:, :-1] + cl.l1[:, c:c + 1])
        E = new
    k = np.arange(C + 1)
    return logsumexp(E + betaln(k + 1.0, C - k + 1.0)[None, :], axis=1)


def log_noise(alt, n, mask, a0, b0) -> np.ndarray:
    a = np.where(mask, alt, 0).astype(float)
    nn = np.where(mask, n, 0).astype(float)
    A, N = a.sum(1), nn.sum(1)
    return G.log_choose(nn, a).sum(1) + betaln(a0 + A, b0 + N - A) - betaln(a0, b0)


# ------------------------------------------------------------------ fitting


@dataclass
class Fit:
    S: np.ndarray
    comps: np.ndarray
    logprior: np.ndarray
    log_tree: np.ndarray
    log_indep: np.ndarray
    log_noise: np.ndarray
    post: np.ndarray
    best: np.ndarray
    bf: np.ndarray
    lr: np.ndarray          # log10 L1/L0 per cell
    cl: G.CellLik


def fit_all(alt, ref, kinds, haploid, P: G.Params, br: T.Branches, opt: Options) -> Fit:
    n = alt + ref
    cl = G.cell_likelihoods(alt, ref, kinds, haploid, P, opt.min_depth)
    S, comps = branch_scores(cl, br, haploid, alt, n)
    lp = branch_log_prior(br, opt.branch_prior, opt.root_prior, opt.floor_frac)
    J = S + lp[None, :]
    log_tree = logsumexp(J, axis=1)
    post = np.exp(J - log_tree[:, None])
    best = J.argmax(1)
    li = log_indep(cl)
    ln = log_noise(alt, n, cl.mask, P.noise_a, P.noise_b)
    bf = (log_tree - (np.logaddexp(li, ln) - math.log(2.0))) / LN10
    lr = (cl.l1 - cl.l0) / LN10
    return Fit(S, comps, lp, log_tree, li, ln, post, best, bf, lr, cl)


def classify(fit: Fit, br: T.Branches, opt: Options):
    conf_car = fit.lr >= opt.carrier_lr
    conf_wt = (-fit.lr) >= opt.wt_lr
    n_car = conf_car.sum(1)
    n_wt = conf_wt.sum(1)
    best_tree_ll = fit.S.max(1)
    noise_best = (fit.log_noise - best_tree_ll) / LN10 >= opt.noise_margin
    cls = np.empty(len(n_car), dtype=object)
    for i in range(len(n_car)):
        b = fit.best[i]
        if n_car[i] == 0:
            cls[i] = "uninformative_depth"
        elif noise_best[i]:
            cls[i] = "noise"
        elif b == 0 and n_wt[i] == 0:
            cls[i] = "germline"
        elif n_car[i] >= 2 and n_wt[i] >= 1:
            cls[i] = "informative_shared"
        elif n_car[i] == 1 and br.is_tip[b]:
            cls[i] = "private"
        else:
            cls[i] = "uninformative_depth"
    lab = np.where(fit.bf >= opt.bf_threshold, "phylo_consistent",
                   np.where(fit.bf <= -opt.bf_threshold, "phylo_violating", "ambiguous"))
    return cls, lab, conf_car, conf_wt


def update_params(alt, ref, kinds, haploid, P: G.Params, fit: Fit, br: T.Branches, opt: Options,
                  cls) -> None:
    """Purity per colony (ML over clade-member cells of confidently placed non-root loci),
    background eps per colony (pooled over the out-of-clade cells of the same loci, shrunk to
    the global rate with 1000 pseudo-reads), rho0 by ML on those cells."""
    n = alt + ref
    best = fit.best
    pbest = fit.post[np.arange(len(best)), best]
    # purity: clade members of informative SHARED loci only. A private locus is classed private
    # because its one cell is a confident carrier, i.e. it is selected on a high alt count,
    # which would bias the purity upwards (most in low-depth colonies).
    conf = ((best != 0) & (pbest >= 0.9) & (fit.bf >= opt.bf_threshold) & (cls == "informative_shared"))
    conf_eps = ((best != 0) & (pbest >= 0.9) & (fit.bf >= opt.bf_threshold)
                & np.isin(cls, ["informative_shared", "private"]))
    if conf_eps.sum() == 0:
        return
    M = br.mask[best[conf]]                      # loci x C : in clade
    a, nn = alt[conf], n[conf]
    K = P.K_matrix(np.asarray(kinds)[conf])
    hap = np.broadcast_to(np.asarray(haploid)[conf][:, None], a.shape)
    C = alt.shape[1]
    pur = np.full(C, np.nan)
    for c in range(C):
        sel = M[:, c] & (nn[:, c] > 0)
        P.n_purity_cells[c] = int(sel.sum())
        if sel.sum() >= 3:
            pur[c] = G.ml_purity(a[sel, c], nn[sel, c], K[sel, c], hap[sel, c], P)
    if np.isfinite(pur).any():
        P.purity = np.where(np.isfinite(pur), pur, np.nanmedian(pur))
    Me = br.mask[best[conf_eps]]
    ae, ne = alt[conf_eps], n[conf_eps]
    out = ~Me & (ne > 0)
    A0, N0 = (ae * out).sum(0), (ne * out).sum(0)
    g = max(float(A0.sum()) / max(float(N0.sum()), 1.0), 1e-4)
    P.eps = np.maximum((A0 + g * 1000.0) / (N0 + 1000.0), 1e-5)
    P.n_eps_cells = out.sum(0)
    r = G.ml_rho(ae[out], ne[out], np.broadcast_to(P.eps[None, :], ae.shape)[out])
    if r is not None:
        P.rho0 = max(r, 1e-4)


def bootstrap_pvalues(idx, fit: Fit, alt, ref, kinds, haploid, P: G.Params, br: T.Branches,
                      opt: Options, rng: np.random.Generator) -> np.ndarray:
    """Parametric-bootstrap p-value of 'the data are consistent with ONE branch' per locus.

    Statistic: G = sum_c max_h log P(d_c | h) - max_b log L_b (saturated minus best branch).
    Null: the locus's best branch with the fitted per-colony fractions, same informative
    depths; p = (1 + #{G_sim >= G_obs}) / (1 + B)."""
    out = np.full(len(fit.best), np.nan)
    if len(idx) == 0 or opt.bootstrap <= 0:
        return out
    B = opt.bootstrap
    n_all = alt + ref
    s1 = (1 - P.rho1) / P.rho1 if P.rho1 > 1e-9 else None
    s0 = (1 - P.rho0) / P.rho0 if P.rho0 > 1e-9 else None
    kinds = np.asarray(kinds)
    haploid = np.asarray(haploid)

    def gstat(cl, S):
        sat = np.maximum.reduce([cl.l1, cl.l0, np.where(np.isfinite(cl.lhet), cl.lhet, -np.inf), cl.lhom]).sum(1)
        return sat - S.max(1)

    def draw(f, s, n, shape):
        f = G._clip(f)
        if s is None:
            u = np.broadcast_to(f, shape)
        else:
            u = rng.beta(f * s, (1 - f) * s, size=shape)
        return rng.binomial(n, u)

    chunk = max(1, 20000 // max(B, 1))
    for st in range(0, len(idx), chunk):
        ii = np.asarray(idx[st:st + chunk])
        Lc, C = len(ii), alt.shape[1]
        n = np.broadcast_to(n_all[ii][:, None, :], (Lc, B, C))
        shape = (Lc, B, C)
        f1 = np.broadcast_to(fit.cl.f1[ii][:, None, :], shape)
        fg = np.broadcast_to(fit.cl.fg[ii][:, None, :], shape)
        eps = np.broadcast_to(P.eps[None, None, :], shape)
        b = fit.best[ii]
        inc = np.broadcast_to(br.mask[b][:, None, :], shape)
        comp = fit.comps[ii].argmax(1)                       # root: 0 het, 1 hom, 2 soma
        root = (b == 0)[:, None, None]
        a_pres = draw(f1, s1, n, shape)
        if P.lam > 0:
            sub = rng.random(shape) < P.lam
            a_sub = rng.binomial(n, rng.random(shape) * f1)
            a_pres = np.where(sub, a_sub, a_pres)
        a_abs = draw(eps, s0, n, shape)
        a_het = draw(fg, s1, n, shape)
        A_i, N_i = alt[ii].sum(1), n_all[ii].sum(1)
        f_hom = ((A_i + HOM_PRIOR[0]) / (N_i + sum(HOM_PRIOR)))[:, None, None]
        a_hom = draw(np.broadcast_to(f_hom, shape), s1, n, shape)
        a_root = np.where((comp == 0)[:, None, None], a_het,
                          np.where((comp == 1)[:, None, None], a_hom, a_pres))
        a_sim = np.where(root, a_root, np.where(inc, a_pres, a_abs))
        a2 = a_sim.reshape(Lc * B, C)
        n2 = n.reshape(Lc * B, C)
        k2 = np.repeat(kinds[ii], B)
        h2 = np.repeat(haploid[ii], B)
        cls = G.cell_likelihoods(a2, n2 - a2, k2, h2, P, opt.min_depth)
        Ss, _ = branch_scores(cls, br, h2, a2, n2)
        gs = gstat(cls, Ss).reshape(Lc, B)
        clo = G.CellLik(*(x[ii] for x in (fit.cl.l1, fit.cl.l0, fit.cl.lhet, fit.cl.lhom,
                                         fit.cl.f1, fit.cl.fg, fit.cl.mask)))
        go = gstat(clo, fit.S[ii])
        out[ii] = (1 + (gs >= go[:, None] - 1e-9).sum(1)) / (1.0 + B)
    return out


def violations(fit: Fit, br: T.Branches, alt, ref, P: G.Params, loci, colonies, cls, opt: Options):
    """Missing leaves (clade members whose reads favour absence) and extra carriers (outside
    colonies whose reads favour presence) under each locus's best branch."""
    n = alt + ref
    rows = []
    want = np.isin(cls, ["informative_shared", "private"])
    for i in np.where(want)[0]:
        b = fit.best[i]
        inc = br.mask[b]
        for c in range(len(colonies)):
            if not fit.cl.mask[i, c]:
                continue
            lr = fit.lr[i, c]
            if inc[c] and lr < 0:
                p = float(G.present_cdf(alt[i, c], n[i, c], fit.cl.f1[i, c], P.rho1, P.lam))
                pdrop = float(G.dropout_prob(n[i, c], fit.cl.f1[i, c], P))
                rows.append((loci[i], colonies[c], "missing_leaf", int(alt[i, c]), int(ref[i, c]),
                             round(float(fit.cl.f1[i, c]), 4), pdrop, p, round(lr, 3), p < opt.alpha))
            elif not inc[c] and lr > 0:
                p = float(G.absent_sf(alt[i, c], n[i, c], P.eps[c], P.rho0))
                rows.append((loci[i], colonies[c], "extra_carrier", int(alt[i, c]), int(ref[i, c]),
                             round(float(P.eps[c]), 5), np.nan, p, round(lr, 3), p < opt.alpha))
    return pd.DataFrame(rows, columns=["locus", "colony", "kind", "n_alt", "n_ref", "expected_fraction",
                                       "p_dropout", "p_chance", "log10_lr_present", "flagged"])


# ------------------------------------------------------------------ driver


def load_samples(path: Optional[str]) -> Optional[List[str]]:
    if not path:
        return None
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    col = "sample" if "sample" in df.columns else df.columns[-1]
    return [s for s in df[col].dropna().tolist() if s]


def load_annotation(path: Optional[str]) -> Optional[pd.DataFrame]:
    if not path:
        return None
    # annotate_v2 writes TAB-separated text even when named `.annotated.csv.gz` (pipeline.sh);
    # sniff the header so a real comma-separated table also works
    import gzip
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        head = fh.readline()
    sep = "\t" if "\t" in head or "," not in head else ","
    df = pd.read_csv(path, sep=sep, dtype=str, compression="infer")
    key = "locus" if "locus" in df.columns else df.columns[0]
    return df.drop_duplicates(key).set_index(key)


def run(args) -> Dict[str, object]:
    opt = Options(min_depth=args.min_depth, branch_prior=args.branch_prior, root_prior=args.root_prior,
                  bf_threshold=args.bf_threshold, alpha=args.alpha, bootstrap=args.bootstrap,
                  seed=args.seed)
    os.makedirs(args.out, exist_ok=True)
    log = []

    def say(msg):
        print(msg)
        log.append(msg)

    tree = T.read_newick(args.tree)
    tips = [t.name for t in tree.tips()]
    wanted = load_samples(args.samples)
    csv_loci = None
    csv = None
    joint = None
    if args.genotypes and os.path.exists(args.genotypes):
        # combine_genotypes call strings (gated loci) or the genotype2 numeric P(carrier) matrix
        csv = G.read_genotypes_csv(args.genotypes)
        if csv.attrs.get("format") == GIO.FMT_CALLS:
            csv_loci = set(csv.index.astype(str))
        else:
            joint = GIO.read_joint_df(GIO.joint_tsv_for(args.genotypes))
        say(f"--genotypes: {csv.attrs.get('format')} matrix, {len(csv)} loci"
            + (f"; joint table {GIO.joint_tsv_for(args.genotypes)}" if joint is not None else ""))
    if args.genotype_dir:
        counts = G.read_genotype_dir(args.genotype_dir)
    else:
        raise SystemExit("tree_fit: --genotype-dir is required: <patient>.genotypes.csv.gz (calls or "
                         "P(carrier)) carries no n_alt/n_ref, and the model needs read counts")
    say(f"--genotype-dir: {len(counts.colonies)} {counts.fmt} per-colony files"
        + (" (votes = n_alt / n_ref; n_uninf, n_art, depth unused; status != ok = missing)"
           if counts.fmt == GIO.FMT_V2 else ""))
    cols = [c for c in counts.colonies if c in set(tips) and (wanted is None or c in set(wanted))]
    dropped_geno = [c for c in counts.colonies if c not in cols]
    dropped_tips = [t for t in tips if t not in cols]
    if len(cols) < 2:
        raise SystemExit(f"tree_fit: only {len(cols)} colonies are both genotyped and in the tree "
                         f"(genotype files: {counts.colonies[:5]}..., tips: {tips[:5]}...)")
    counts = counts.subset_colonies(cols)
    tree = T.prune(tree, cols)
    cols = [t.name for t in tree.tips()]
    counts = counts.subset_colonies(cols)
    br = T.branches(tree, cols)
    say(f"colonies: {len(cols)} fitted; {len(dropped_geno)} genotyped but not in tree/samples "
        f"({', '.join(dropped_geno[:8])}); {len(dropped_tips)} tips without genotypes "
        f"({', '.join(dropped_tips[:8])})")
    say(f"branches: {len(br.ids)} (root + {int(br.is_tip.sum())} tips + {len(br.ids) - 1 - int(br.is_tip.sum())} internal)")

    loci = counts.loci
    alt, ref = counts.alt, counts.ref
    kinds = np.array([G.locus_kind(x) for x in loci])
    contigs = np.array([G.locus_contig(x) for x in loci])
    sex_info = {}
    if args.sex == "auto":
        sex, sex_info = G.infer_sex(loci, alt, ref)
        say(f"sex (auto): {sex or 'undecided -> F (diploid X)'} {sex_info}")
        sex = sex or "F"
    else:
        sex = args.sex
    haploid = np.isin(contigs, ["chrY", "Y"]) | (np.isin(contigs, ["chrX", "X"]) & (sex == "M"))
    autosomal = ~np.isin(contigs, list(G.SEX_CONTIGS))

    P = G.Params(colonies=cols, lam=args.subclonal_weight, sex=sex)
    cand = G.estimate_balance(alt, ref, kinds, autosomal, P)
    say(f"germline-het candidates: {int(cand.sum())}; K_global={P.K_global:.3f} "
        f"K_kind={ {k: round(v, 3) for k, v in P.K_kind.items()} }; rho1={P.rho1}")
    fit = None
    for r in range(opt.rounds):
        fit = fit_all(alt, ref, kinds, haploid, P, br, opt)
        cls, lab, _, cwt = classify(fit, br, opt)
        update_params(alt, ref, kinds, haploid, P, fit, br, opt, cls)
        # refine the vote odds on the loci now fitted as germline het (no confident wild-type)
        ghet = (cls == "germline") & (fit.comps.argmax(1) == 0) & autosomal
        G.estimate_balance(alt, ref, kinds, autosomal, P, cand=ghet)
        say(f"round {r + 1}: purity median {np.median(P.purity):.3f} "
            f"[{P.purity.min():.2f}-{P.purity.max():.2f}], eps median {np.median(P.eps):.5f}, rho0={P.rho0}")
    fit = fit_all(alt, ref, kinds, haploid, P, br, opt)
    cls, lab, conf_car, conf_wt = classify(fit, br, opt)

    rng = np.random.default_rng(opt.seed)
    idx = np.where(np.isin(cls, ["informative_shared", "private"]))[0]
    pboot = bootstrap_pvalues(idx, fit, alt, ref, kinds, haploid, P, br, opt, rng)
    viol = violations(fit, br, alt, ref, P, loci, cols, cls, opt)

    # ---- per-locus table
    nl = len(loci)
    order = np.argsort(-fit.post, axis=1)
    second = order[:, 1] if fit.post.shape[1] > 1 else order[:, 0]
    ar = np.arange(nl)
    n = alt + ref
    vflag = viol[viol["flagged"]] if len(viol) else viol
    nmiss = vflag[vflag["kind"] == "missing_leaf"].groupby("locus").size() if len(viol) else pd.Series(dtype=int)
    nextra = vflag[vflag["kind"] == "extra_carrier"].groupby("locus").size() if len(viol) else pd.Series(dtype=int)
    pmiss = viol[viol["kind"] == "missing_leaf"].groupby("locus")["p_chance"].min() if len(viol) else pd.Series(dtype=float)
    pextra = viol[viol["kind"] == "extra_carrier"].groupby("locus")["p_chance"].min() if len(viol) else pd.Series(dtype=float)
    root_comp = np.array(["het", "hom", "somatic_all"])[fit.comps.argmax(1)]
    df = pd.DataFrame({
        "locus": loci,
        "contig": contigs,
        "locus_kind": kinds,
        "ploidy": np.where(haploid, 1, 2),
        "n_colonies_data": fit.cl.mask.sum(1),
        "total_alt": alt.sum(1),
        "total_ref": ref.sum(1),
        "n_carriers_observed": conf_car.sum(1),
        "n_confident_wt": conf_wt.sum(1),
        "carriers": [",".join(np.array(cols)[conf_car[i]]) for i in range(nl)],
        "best_branch": np.array(br.ids)[fit.best],
        "best_clade": [",".join(br.members[b]) for b in fit.best],
        "n_clade": br.mask[fit.best].sum(1),
        "branch_length": br.lengths[fit.best],
        "post_best_branch": np.round(fit.post[ar, fit.best], 4),
        "second_branch": np.array(br.ids)[second],
        "post_second_branch": np.round(fit.post[ar, second], 4),
        "root_component": np.where(fit.best == 0, root_comp, "."),
        "log10_L_best": np.round(fit.S[ar, fit.best] / LN10, 3),
        "log10_L_root": np.round(fit.S[:, 0] / LN10, 3),
        "log10_L_indep": np.round(fit.log_indep / LN10, 3),
        "log10_L_noise": np.round(fit.log_noise / LN10, 3),
        "log10_bf_tree": np.round(fit.bf, 3),
        "n_missing_leaves": [int(nmiss.get(x, 0)) for x in loci],
        "n_extra_carriers": [int(nextra.get(x, 0)) for x in loci],
        "min_p_missing": [pmiss.get(x, np.nan) for x in loci],
        "min_p_extra": [pextra.get(x, np.nan) for x in loci],
        "p_locus": pboot,
        "class": cls,
        "label": lab,
        # one verdict per locus for downstream tools (cluster/tprt/compare_arms.py picks this
        # column first): the BF label where the tree can judge (informative_shared), the
        # constant-fraction artefact as `noise_violating`, else the class itself
        "phylo_label": np.where(cls == "informative_shared", lab,
                                np.where(cls == "noise", "noise_violating", cls)),
    })
    if csv_loci is not None:
        df["passes_combine_genotypes"] = [x in csv_loci for x in loci]
    elif csv is not None:
        # numeric matrix (genotype2 joint step): its per-colony P(carrier) over the fitted colonies,
        # and the joint step's own verdict when <P>.joint.tsv sits next to it -- the independent
        # cross-check of the Rust port against this read-vote model
        m = csv.reindex(index=loci, columns=[c for c in cols if c in csv.columns])
        df["matrix_n_carriers"] = (m >= GIO.P_CARRIER).sum(axis=1).values
        df["matrix_n_wt"] = (m <= GIO.P_ABSENT_MATRIX).sum(axis=1).values
    if joint is not None:
        j = joint.reindex(loci)
        best_j = [b if isinstance(b, str) else "" for b in j["best"]]
        df["joint_best"] = [b or "NA" for b in best_j]
        df["joint_class"] = [GIO.joint_class({"best": b, "n_carriers": n})
                             for b, n in zip(best_j, j["n_carriers"].fillna(0))]
        df["joint_n_carriers"] = j["n_carriers"].values
        df["joint_log10_bf_tree"] = j["log10_bf_tree"].values
    ann = load_annotation(args.annotation)
    if ann is not None:
        ann = ann.rename(columns={c: (f"ann_{c}" if c in df.columns else c) for c in ann.columns})
        df = df.join(ann, on="locus")
    df.to_csv(os.path.join(args.out, "phylo_fit.tsv"), sep="\t", index=False, na_rep="NA")
    viol.to_csv(os.path.join(args.out, "violations.tsv"), sep="\t", index=False, na_rep="NA")

    # ---- per-colony parameters
    med_n = np.array([np.median(n[n[:, c] > 0, c]) if (n[:, c] > 0).any() else 0 for c in range(len(cols))])
    Kt = P.K_kind.get("TSD", P.K_global) * P.b
    f_tsd = G.f_present(Kt, P.purity, False)
    pdrop_med = G.dropout_prob(med_n, f_tsd, P)
    cp = pd.DataFrame({
        "colony": cols, "purity": np.round(P.purity, 3), "purity_n_cells": P.n_purity_cells,
        "allelic_balance": np.round(P.b, 3), "epsilon": np.round(P.eps, 6), "epsilon_n_cells": P.n_eps_cells,
        "median_depth": med_n, "mean_depth": np.round(np.where(n > 0, n, np.nan).mean(0), 2) if n.size else 0,
        "f_present_tsd": np.round(f_tsd, 3), "p_dropout_median_depth": np.round(pdrop_med, 5),
        "terminal_branch_length": [br.lengths[br.index(c)] for c in cols],
    })
    cp.to_csv(os.path.join(args.out, "colony_params.tsv"), sep="\t", index=False)

    # ---- per-cell table (calibration / plotting)
    if not args.no_cells:
        L, C = alt.shape
        ii, cc = np.where(fit.cl.mask)
        inc = br.mask[fit.best[ii], cc]
        cells = pd.DataFrame({
            "locus": np.array(loci)[ii], "colony": np.array(cols)[cc],
            "n_alt": alt[ii, cc], "n_ref": ref[ii, cc],
            "f_present": np.round(fit.cl.f1[ii, cc], 4), "epsilon": P.eps[cc],
            "log10_lr_present": np.round(fit.lr[ii, cc], 3),
            "p_dropout": G.dropout_prob(n[ii, cc], fit.cl.f1[ii, cc], P),
            "p_le_present": G.present_cdf(alt[ii, cc], n[ii, cc], fit.cl.f1[ii, cc], P.rho1, P.lam),
            "p_ge_absent": G.absent_sf(alt[ii, cc], n[ii, cc], P.eps[cc], P.rho0),
            "in_best_clade": inc,
        })
        cells.to_csv(os.path.join(args.out, "cells.tsv.gz"), sep="\t", index=False, float_format="%.4g")

    write_summary(args, df, cp, viol, P, opt, br, log, sex_info)
    return {"fit": df, "colony_params": cp, "violations": viol, "params": P}


def write_summary(args, df, cp, viol, P, opt, br, log, sex_info):
    tab = pd.crosstab(df["class"], df["label"]).reindex(index=list(CLASSES), columns=list(LABELS)).fillna(0).astype(int)
    lines = ["# Phylogenetic fit summary", "",
             f"tree `{args.tree}`, genotype dir `{args.genotype_dir}`, {len(br.colonies)} colonies, "
             f"{len(df)} loci.", "", "## Run log", ""]
    lines += [f"    {x}" for x in log]
    lines += ["", "## Global parameters", "",
              f"* sex {P.sex} ({sex_info}); haploid loci: {int((df['ploidy'] == 1).sum())}",
              f"* vote odds K (per haplotype copy): global {P.K_global:.3f}; per kind "
              f"{ {k: round(v, 3) for k, v in P.K_kind.items()} } from {P.n_germline_het} germline-het loci",
              f"* overdispersion rho1 (present) {P.rho1}, rho0 (absent) {P.rho0}; subclonal weight {P.lam}",
              f"* branch prior `{opt.branch_prior}`, root prior {opt.root_prior}; BF threshold +-{opt.bf_threshold} "
              f"(log10); violation alpha {opt.alpha}; bootstrap B={opt.bootstrap}",
              "", "## Loci by class x label", "", "| class | " + " | ".join(LABELS) + " | total |",
              "|---|" + "---|" * (len(LABELS) + 1)]
    for c in tab.index:
        lines.append(f"| {c} | " + " | ".join(str(v) for v in tab.loc[c]) + f" | {tab.loc[c].sum()} |")
    sh = df[df["class"] == "informative_shared"]
    if len(sh):
        pl = sh["p_locus"].dropna()
        lines += ["", f"Informative shared loci: {len(sh)}; median clade size {sh['n_clade'].median():.0f}; "
                      f"bootstrap p < 0.01: {int((pl < 0.01).sum())} of {len(pl)}."]
    if len(viol):
        v = viol[viol["flagged"]]
        lines += ["", f"Flagged violations (p < {opt.alpha}): {int((v['kind'] == 'missing_leaf').sum())} missing "
                      f"leaves, {int((v['kind'] == 'extra_carrier').sum())} extra carriers "
                      f"(unflagged candidates: {int((~viol['flagged']).sum())})."]
    lines += ["", "## Colonies", "", "| colony | purity | balance | eps | median depth | P(dropout) at median depth |",
              "|---|---|---|---|---|---|"]
    for _, r in cp.iterrows():
        lines.append(f"| {r.colony} | {r.purity} | {r.allelic_balance} | {r.epsilon:.2e} | {r.median_depth:g} | "
                     f"{r.p_dropout_median_depth:.3g} |")
    worst = sh.sort_values("log10_bf_tree").head(15) if len(sh) else sh
    if len(worst):
        lines += ["", "## Most violating informative loci", "",
                  "| locus | carriers | best clade | log10 BF | missing | extra | p_locus |", "|---|---|---|---|---|---|---|"]
        for _, r in worst.iterrows():
            lines.append(f"| {r.locus} | {r.carriers} | {r.best_branch} ({r.n_clade}) | {r.log10_bf_tree} | "
                         f"{r.n_missing_leaves} | {r.n_extra_carriers} | {r.p_locus} |")
    if "joint_class" in df.columns:
        jt = pd.crosstab(df["class"], df["joint_class"])
        jc = [c for c in ("ROOT", "clade", "private", "INDEP", "NOISE", "NA") if c in jt.columns]
        jt = jt.reindex(index=[c for c in CLASSES if c in jt.index], columns=jc).fillna(0).astype(int)
        lines += ["", "## Cross-check: tree_fit class x genotype2 joint step", "",
                  "Rows: this read-vote model; columns: the Rust joint step (`<P>.joint.tsv`, PL-based; "
                  "clade = a branch with >= 2 carriers, private = a tip).", "",
                  "| class | " + " | ".join(jc) + " |", "|---|" + "---|" * len(jc)]
        for c in jt.index:
            lines.append(f"| {c} | " + " | ".join(str(v) for v in jt.loc[c]) + " |")
    lines += ["", "How to read: plans/tprt_hallmarks/PHYLO_EVAL.md.", ""]
    with open(os.path.join(args.out, "summary.md"), "w") as fh:
        fh.write("\n".join(lines))


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", help="<patient>.genotypes.csv.gz, auto-detected: combine_genotypes call "
                                        "strings (flags the loci passing its gates) or the genotype2 "
                                        "numeric P(carrier) matrix (adds matrix_n_carriers / matrix_n_wt, and "
                                        "joint_* columns from <patient>.joint.tsv next to it)")
    ap.add_argument("--genotype-dir", help="directory of per-colony genotype files <sample>.txt.gz, legacy "
                                           "or genotype2 numeric (auto-detected; n_alt / n_ref columns)")
    ap.add_argument("--tree", required=True, help="Newick SNV tree (tip labels = colony sample names)")
    ap.add_argument("--annotation", help="annotate_v2 table (TSV, `locus` column) joined onto phylo_fit.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sex", choices=["M", "F", "auto"], default="auto")
    ap.add_argument("--samples", help="colonies.tsv (`sample` column): restrict to these colonies")
    ap.add_argument("--min-depth", type=int, default=1,
                    help="cells with fewer informative reads are treated as missing data")
    ap.add_argument("--branch-prior", choices=["length", "uniform"], default="length")
    ap.add_argument("--root-prior", type=float, default=0.1)
    ap.add_argument("--bf-threshold", type=float, default=1.0)
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--subclonal-weight", type=float, default=0.01)
    ap.add_argument("--bootstrap", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-cells", action="store_true", help="do not write cells.tsv.gz")
    return ap


def main(argv=None):
    run(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
