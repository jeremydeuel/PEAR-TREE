#!/usr/bin/env python3
"""Count-level simulator for the phylogenetic evaluation: per-colony genotype files
(`<colony>.txt.gz`, the Rust genotyper's columns) generated from a tree, without reads.

Loci (truth in <out>/truth.tsv):
  clade     somatic insertion on a branch (drawn proportional to branch length, or uniform),
            carried by exactly the clade; per-colony expected alt-vote fraction
            f = (p/2) K / ((p/2) K + 1 - p/2), p = colony purity (U(purity range))
  germline  every colony, het (f = K / (K + 1)) or hom
  nonclade  carried by a random NON-clade subset (>= 2 carriers, >= 1 non-carrier) at the
            clonal fraction (an artefact / mis-merge that scatters over the tree)
  noise     no carrier; a locus-specific alt rate ~ U(noise range) in every colony (the
            constant-allele-fraction artefact)
Absent colonies get alt votes at the background rate eps (per colony U(eps range)); depths are
negative binomial around the colony's mean depth (a fraction of colonies is low-depth).
Counts are beta-binomial with intra-class correlation --rho.

  python tools/phylo/simulate_counts.py --tree random:12 --out sim/ --seed 3
  python tools/phylo/simulate_counts.py --tree patients/colon/PD44890/PD44890_snp_tree_with_branch_length.tree --out sim/
"""
from __future__ import annotations

import argparse
import gzip
import os
import random
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from tools.phylo import tree as T  # noqa: E402


def load_tree(spec: str, rng: random.Random) -> T.Node:
    if spec.startswith("random:"):
        return T.random_coalescent(int(spec.split(":", 1)[1]), rng)
    return T.read_newick(spec)


def nonclade_subset(rng: random.Random, colonies, clade_sets, tries: int = 1000):
    n = len(colonies)
    if n < 3:
        return None
    for _ in range(tries):
        k = rng.randint(2, n - 1)
        s = frozenset(rng.sample(colonies, k))
        if s not in clade_sets:
            return s
    return None


def bb_draw(nrng, n, f, rho):
    f = np.clip(f, 1e-6, 1 - 1e-6)
    if rho > 0:
        s = (1 - rho) / rho
        f = nrng.beta(f * s, (1 - f) * s)
    return nrng.binomial(n, f)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True, help="Newick file or random:N")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-clade", type=int, default=600)
    ap.add_argument("--n-germline", type=int, default=300)
    ap.add_argument("--n-nonclade", type=int, default=150)
    ap.add_argument("--n-noise", type=int, default=100)
    ap.add_argument("--n-chrx", type=int, default=20, help="extra germline loci on chrX")
    ap.add_argument("--sex", choices=["M", "F"], default="F")
    ap.add_argument("--branch-weight", choices=["length", "uniform"], default="uniform")
    ap.add_argument("--depth", type=float, default=15.0)
    ap.add_argument("--low-depth-frac", type=float, default=0.25)
    ap.add_argument("--low-depth-range", default="0.15,0.4", help="depth factor of low-depth colonies")
    ap.add_argument("--purity-range", default="0.7,1.0")
    ap.add_argument("--eps-range", default="0.001,0.006")
    ap.add_argument("--noise-range", default="0.05,0.25")
    ap.add_argument("--K", type=float, default=2.0, help="alt/ref vote odds per haplotype copy (TSD ~ 2)")
    ap.add_argument("--rho", type=float, default=0.01)
    ap.add_argument("--score-shift", type=float, default=2.0,
                    help="synthetic annotated.tsv: mean tprt_score gain of real insertions")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv)

    rng = random.Random(a.seed)
    nrng = np.random.default_rng(a.seed)
    tree = load_tree(a.tree, rng)
    br = T.branches(tree)
    cols = br.colonies
    C = len(cols)
    lo, hi = map(float, a.purity_range.split(","))
    purity = np.array([rng.uniform(lo, hi) for _ in cols])
    dlo, dhi = map(float, a.low_depth_range.split(","))
    dfac = np.array([rng.uniform(dlo, dhi) if rng.random() < a.low_depth_frac else rng.uniform(0.8, 1.2)
                     for _ in cols])
    elo, ehi = map(float, a.eps_range.split(","))
    eps = np.array([rng.uniform(elo, ehi) for _ in cols])
    clade_sets = {frozenset(m) for m in br.members}
    w = br.lengths[1:] if a.branch_weight == "length" else np.ones(len(br.ids) - 1)
    w = np.maximum(w, 1e-9) / np.maximum(w, 1e-9).sum()

    def depths():
        mean = a.depth * dfac
        return nrng.negative_binomial(5, 5 / (5 + mean))    # overdispersed around the mean

    K = a.K
    rows, truth = [], []
    pos = 1_000_000

    def add(kind, carriers, f_vec, contig="chr1"):
        nonlocal pos
        pos += 10_000
        name = f"{contig}:{pos}-{pos + 15}"
        n = depths()
        alt = bb_draw(nrng, n, f_vec, a.rho)
        rows.append((name, alt, n - alt))
        truth.append((name, kind, ",".join(sorted(carriers)), len(carriers)))

    h = purity / 2
    f_som = h * K / (h * K + 1 - h)
    for _ in range(a.n_clade):
        b = 1 + int(nrng.choice(len(w), p=w))
        inc = br.mask[b]
        add("clade" if inc.sum() > 1 else "private", br.members[b], np.where(inc, f_som, eps))
    for _ in range(a.n_germline):
        f = K / (K + 1) if rng.random() < 0.7 else 1 - eps
        add("germline", cols, np.broadcast_to(f, (C,)))
    for _ in range(a.n_nonclade):
        s = nonclade_subset(rng, cols, clade_sets)
        if s is None:
            break
        inc = np.array([c in s for c in cols])
        add("nonclade", s, np.where(inc, f_som, eps))
    for _ in range(a.n_noise):
        e = rng.uniform(*map(float, a.noise_range.split(",")))
        add("noise", [], np.full(C, e))
    for _ in range(a.n_chrx):
        if a.sex == "M":
            f = 1 - eps
        else:
            f = K / (K + 1) if rng.random() < 0.7 else 1 - eps
        add("germline", cols, np.broadcast_to(f, (C,)), contig="chrX")

    os.makedirs(a.out, exist_ok=True)
    gdir = os.path.join(a.out, "genotypes")
    os.makedirs(gdir, exist_ok=True)
    for j, c in enumerate(cols):
        with gzip.open(os.path.join(gdir, f"{c}.txt.gz"), "wt") as fh:
            fh.write("insertion\tgenotype\tscore_genotype\tscore_alternative\tcoverage\tn_alt\tn_ref\tn_art\n")
            for name, alt, ref in rows:
                x, y = int(alt[j]), int(ref[j])
                v = x / (x + y) if x + y else 0
                gt = ("no-coverage" if x + y == 0 else "homozygous" if v >= 0.8 else
                      "heterozygous" if v >= 0.3 and x >= 2 else "wild-type?" if v > 0.05 else "wild-type")
                fh.write(f"{name}\t{gt}\t0\t0\t{x + y}\t{x}\t{y}\t0\n")
    with open(os.path.join(a.out, "truth.tsv"), "w") as fh:
        fh.write("locus\ttruth\tcarriers\tn_carriers\n")
        for t in truth:
            fh.write("\t".join(map(str, t)) + "\n")
    # synthetic annotate_v2-like table: real insertions (clade/private/germline) score higher
    # than artefacts (nonclade/noise) by --score-shift, with two informative point features and
    # one uninformative one (for exercising discrimination.py; not a model of annotate)
    with open(os.path.join(a.out, "annotated.tsv"), "w") as fh:
        fh.write("locus\tclass\telement\tstructure\ttprt_score\ttprt_points\ttprt_call\n")
        for name, kind, car, ncar in truth:
            real = kind in ("clade", "private", "germline")
            tsd = 3 if (rng.random() < (0.85 if real else 0.4)) else 0
            pa = 1 if (rng.random() < (0.8 if real else 0.5)) else 0
            mc = 1 if ncar >= 2 else 0
            junk = rng.choice([0, 1])
            score = tsd + pa + mc + junk + rng.gauss(a.score_shift if real else 0.0, 1.5)
            pts = ";".join(f"{k}:{v:+g}" for k, v in (("tsd_4_25", tsd), ("polya_ge10", pa),
                                                      ("multi_colony", mc), ("uninformative", junk)) if v)
            call = "TPRT" if score >= 7 else "LIKELY_TPRT" if score >= 4 else "UNCERTAIN" if score >= 0 else "ARTEFACT_LIKE"
            el = rng.choice(["L1", "ALU", "SVA"])
            fh.write(f"{name}\t{'RTE' if real else 'artefact'}\t{el}\tTRUNCATED_5P\t{score:.2f}\t{pts or '.'}\t{call}\n")
    with open(os.path.join(a.out, "colonies_truth.tsv"), "w") as fh:
        fh.write("colony\tpurity\tdepth_factor\tepsilon\n")
        for j, c in enumerate(cols):
            fh.write(f"{c}\t{purity[j]:.4f}\t{dfac[j]:.3f}\t{eps[j]:.5f}\n")
    with open(os.path.join(a.out, "tree.nwk"), "w") as fh:
        fh.write(T.to_newick(tree) + "\n")
    print(f"simulated {len(rows)} loci x {C} colonies -> {a.out}")


if __name__ == "__main__":
    main()
