"""Phylogenetic evaluation (tools/phylo, plans/tprt_hallmarks/PHYLO_EVAL.md): tree parsing and
branch enumeration, read-vote likelihood math, the exact independent-presence marginal, and the
dropout / extra-carrier decisions that the method exists for.

Run:  pytest test/test_phylo_model.py
"""
import itertools
import math
import os
import random
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, REPO)
from tools.phylo import genotype_likelihood as G  # noqa: E402
from tools.phylo import tree as T  # noqa: E402
from tools.phylo import tree_fit as TF  # noqa: E402

SANGER = ("((PD1a_lo1:1402,(PD1a_lo2:1686,(PD1a_lo3:1698,PD1a_lo4:1597)100:13)18:0)20:2,PD1a_lo5:1507);")


# ------------------------------------------------------------------ tree


def test_parse_sanger_newick_with_support_and_zero_lengths():
    r = T.parse_newick(SANGER)
    assert [t.name for t in r.tips()] == ["PD1a_lo1", "PD1a_lo2", "PD1a_lo3", "PD1a_lo4", "PD1a_lo5"]
    br = T.branches(r)
    # binary tree with 5 tips: root + 2*5-2 branches
    assert len(br.ids) == 9 and br.ids[0] == "ROOT"
    assert br.mask[0].all()
    assert int(br.is_tip.sum()) == 5
    i = br.clade_index(["PD1a_lo3", "PD1a_lo4"])
    assert i is not None and br.lengths[i] == 13
    j = br.clade_index(["PD1a_lo2", "PD1a_lo3", "PD1a_lo4"])
    assert br.lengths[j] == 0.0                     # the `18:0` branch is kept, length 0


def test_prune_collapses_unary_nodes_and_sums_lengths():
    r = T.parse_newick("((A:1,B:2):3,(C:4,D:5):6);")
    p = T.prune(r, ["A", "C", "D"])
    br = T.branches(p)
    assert sorted(t.name for t in p.tips()) == ["A", "C", "D"]
    assert br.lengths[br.index("A")] == 4          # 1 + 3 (unary parent collapsed)
    assert br.clade_index(["C", "D"]) is not None
    assert len(br.ids) == 5


def test_random_coalescent_is_reproducible_and_binary():
    a = T.to_newick(T.random_coalescent(7, random.Random(3)))
    b = T.to_newick(T.random_coalescent(7, random.Random(3)))
    assert a == b
    br = T.branches(T.parse_newick(a))
    assert len(br.ids) == 2 * 7 - 1


# ------------------------------------------------------------------ likelihood math


@pytest.mark.parametrize("rho", [0.0, 0.02, 0.2])
def test_beta_binomial_pmf_normalises(rho):
    n = 17
    p = np.exp(G.bb_logpmf(np.arange(n + 1), n, 0.3, rho))
    assert abs(p.sum() - 1) < 1e-9


def test_beta_binomial_reduces_to_binomial():
    a, n, f = 3, 11, 0.4
    binom = math.comb(n, a) * f ** a * (1 - f) ** (n - a)
    assert abs(math.exp(G.bb_logpmf(a, n, f, 0.0)) - binom) < 1e-12
    assert abs(math.exp(G.bb_logpmf(a, n, f, 1e-7)) - binom) < 1e-5


def test_present_pmf_with_subclonal_component_normalises_and_cdf_matches():
    n, f1 = 12, 0.6
    pm = np.exp(G.present_logpmf(np.arange(n + 1), n, f1, 0.02, 0.05))
    assert abs(pm.sum() - 1) < 1e-9
    assert abs(G.present_cdf(3, n, f1, 0.02, 0.05) - pm[:4].sum()) < 1e-9


def test_absent_sf():
    n, eps = 30, 0.01
    pm = np.exp(G.bb_logpmf(np.arange(n + 1), n, eps, 0.0))
    assert abs(G.absent_sf(2, n, eps, 0.0) - pm[2:].sum()) < 1e-9
    assert G.absent_sf(0, n, eps, 0.0) == 1.0


def test_expected_fractions():
    # diploid clonal het, pure colony, unbiased votes: 0.5; purity 0.8 -> 0.4; haploid pure -> ~1
    assert abs(G.f_present(1.0, 1.0, False) - 0.5) < 1e-9
    assert abs(G.f_present(1.0, 0.8, False) - 0.4) < 1e-9
    assert G.f_present(1.0, 1.0, True) > 0.999
    # TSD-like vote odds K = 2: germline het reads 2/3
    assert abs(G.f_germline_het(2.0) - 2 / 3) < 1e-9


def test_log_indep_dp_equals_brute_force():
    rng = np.random.default_rng(4)
    C = 5
    l1 = rng.normal(-2, 2, (3, C))
    l0 = rng.normal(-2, 2, (3, C))
    cl = G.CellLik(l1, l0, l1, l1, l1, l1, np.ones_like(l1, bool))
    got = TF.log_indep(cl)
    for i in range(3):
        # integrate over pi numerically: prod_c (pi e^l1 + (1-pi) e^l0)
        pis = np.linspace(0, 1, 20001)
        prod = np.ones_like(pis)
        for c in range(C):
            prod *= pis * math.exp(l1[i, c]) + (1 - pis) * math.exp(l0[i, c])
        ref = math.log(np.trapezoid(prod, pis))
        assert abs(got[i] - ref) < 1e-6
        # and the subset enumeration with Beta(k+1, C-k+1) weights
        tot = 0.0
        for s in itertools.product([0, 1], repeat=C):
            k = sum(s)
            w = math.exp(math.lgamma(k + 1) + math.lgamma(C - k + 1) - math.lgamma(C + 2))
            tot += w * math.exp(sum(l1[i, c] if s[c] else l0[i, c] for c in range(C)))
        assert abs(got[i] - math.log(tot)) < 1e-9


def test_log_noise_closed_form():
    alt = np.array([[1, 2, 0]])
    n = np.array([[10, 12, 9]])
    got = TF.log_noise(alt, n, n > 0, 0.5, 2.0)[0]
    eps = np.linspace(1e-9, 1 - 1e-9, 200001)
    prior = eps ** -0.5 * (1 - eps) ** 1.0 / math.exp(math.lgamma(0.5) + math.lgamma(2) - math.lgamma(2.5))
    like = np.ones_like(eps)
    for a, m in zip(alt[0], n[0]):
        like *= math.comb(int(m), int(a)) * eps ** a * (1 - eps) ** (m - a)
    assert abs(got - math.log(np.trapezoid(like * prior, eps))) < 1e-3


# ------------------------------------------------------------------ dropout / extra carrier decisions


def _fit_one(rows, eps=0.002, purity=1.0, tree="((A:10,B:10):10,(C:10,D:10):10);"):
    """Fit loci given as {colony: (n_alt, n_ref)} on a 4-colony tree with fixed parameters
    (K = 1: a clonal het reads 0.5 at purity 1)."""
    root = T.parse_newick(tree)
    cols = [t.name for t in root.tips()]
    br = T.branches(root, cols)
    alt = np.array([[r[c][0] for c in cols] for r in rows])
    ref = np.array([[r[c][1] for c in cols] for r in rows])
    kinds = ["TSD"] * len(rows)
    hap = np.zeros(len(rows), bool)
    P = G.Params(colonies=cols, K_global=1.0, purity=np.full(len(cols), purity), eps=np.full(len(cols), eps),
                 rho0=0.001, rho1=0.01, lam=0.01)
    opt = TF.Options(bootstrap=0)
    fit = TF.fit_all(alt, ref, kinds, hap, P, br, opt)
    cls, lab, _, _ = TF.classify(fit, br, opt)
    viol = TF.violations(fit, br, alt, ref, P, [f"L{i}" for i in range(len(rows))], cols, cls, opt)
    return fit, br, cls, lab, viol, P


# A and E carry the insertion strongly; the smallest clade holding both is (A,(B,E)), so B is a
# clade member. Whether B's zero alt reads contradict the tree depends only on B's depth.
TREE5 = "((A:10,(B:10,E:10):10):10,(C:10,D:10):10);"


def _b_member(b_ref):
    return [{"A": (10, 10), "E": (11, 9), "B": (0, b_ref), "C": (0, 30), "D": (0, 25)}]


def test_low_depth_clade_member_with_no_alt_is_not_a_violation():
    fit, br, cls, lab, viol, P = _fit_one(_b_member(2), tree=TREE5)
    assert br.members[fit.best[0]] == ["A", "B", "E"]   # B is placed in the clade
    v = viol[viol.colony == "B"]
    assert len(v) == 1 and v.kind.iloc[0] == "missing_leaf"
    assert not v.flagged.iloc[0]                         # 0 of 2 reads: chance dropout
    assert v.p_dropout.iloc[0] > 0.2                     # P(0 alt | carrier, 2 reads) ~ 0.25
    assert cls[0] == "informative_shared" and lab[0] != "phylo_violating"   # 5 colonies: modest BF


def test_deep_clade_member_with_no_alt_is_a_violation():
    fit, br, cls, lab, viol, P = _fit_one(_b_member(40), tree=TREE5)
    assert cls[0] == "informative_shared"
    assert lab[0] == "phylo_violating"
    flagged = viol[viol.flagged]
    assert len(flagged) >= 1
    if br.members[fit.best[0]] == ["A", "B", "E"]:      # B is the flagged missing leaf
        v = viol[viol.colony == "B"].iloc[0]
        assert v.kind == "missing_leaf" and v.flagged and v.p_dropout < 1e-9
    # the per-cell numbers themselves
    assert G.present_cdf(0, 40, 0.5, 0.01, 0.01) < 0.01
    assert G.dropout_prob(40, 0.5, G.Params(colonies=["x"], rho1=0.01)) < 1e-9


def test_extra_carrier_with_one_read_at_high_background_is_not_a_violation():
    rows = [{"A": (10, 10), "B": (9, 12), "C": (1, 29), "D": (0, 25)}]
    fit, br, cls, lab, viol, P = _fit_one(rows, eps=0.03)
    assert br.members[fit.best[0]] == ["A", "B"]
    v = viol[viol.colony == "C"]
    assert len(v) == 0 or not v.flagged.iloc[0]       # 1/30 alt at eps 3 % is expected
    assert lab[0] != "phylo_violating" and fit.bf[0] > 0
    # the same read at a 0.1 % background, with 4 more: a real extra carrier
    rows = [{"A": (10, 10), "B": (9, 12), "C": (5, 25), "D": (0, 25), "E": (0, 30)}]
    fit, br, cls, lab, viol, P = _fit_one(rows, eps=0.001, tree="((A:10,B:10):10,(C:10,(D:10,E:10):10):10);")
    assert br.members[fit.best[0]] == ["A", "B"]
    assert lab[0] == "phylo_violating"
    v = viol[viol.colony == "C"]
    assert len(v) == 1 and v.kind.iloc[0] == "extra_carrier" and v.flagged.iloc[0]


def test_constant_fraction_everywhere_is_noise_not_germline():
    rows = [{"A": (2, 14), "B": (3, 15), "C": (2, 16), "D": (2, 13)}]
    fit, br, cls, lab, viol, P = _fit_one(rows)
    assert cls[0] == "noise" and lab[0] == "phylo_violating"


def test_germline_het_is_germline():
    rows = [{"A": (10, 9), "B": (8, 12), "C": (11, 10), "D": (9, 9)}]
    fit, br, cls, lab, viol, P = _fit_one(rows)
    assert fit.best[0] == 0 and cls[0] == "germline"


# ------------------------------------------------------------------ end to end on simulated counts


def test_tree_fit_and_discrimination_cli_on_simulated_counts(tmp_path):
    from tools.phylo import calibration, discrimination, simulate_counts
    sim = tmp_path / "sim"
    simulate_counts.main(["--tree", "random:8", "--out", str(sim), "--seed", "5", "--n-clade", "200",
                          "--n-germline", "80", "--n-nonclade", "50", "--n-noise", "30", "--n-chrx", "8"])
    out = tmp_path / "fit"
    TF.main(["--genotypes", str(tmp_path / "missing.csv.gz"), "--genotype-dir", str(sim / "genotypes"),
             "--tree", str(sim / "tree.nwk"), "--annotation", str(sim / "annotated.tsv"),
             "--out", str(out), "--sex", "auto", "--min-depth", "1", "--bootstrap", "20"])
    for f in ("phylo_fit.tsv", "violations.tsv", "colony_params.tsv", "summary.md", "cells.tsv.gz"):
        assert (out / f).exists()
    import pandas as pd
    fit = pd.read_csv(out / "phylo_fit.tsv", sep="\t")
    truth = pd.read_csv(sim / "truth.tsv", sep="\t").set_index("locus")
    fit["truth"] = fit.locus.map(truth.truth)
    sh = fit[fit["class"] == "informative_shared"]
    clade = sh[sh.truth == "clade"]
    nonc = sh[sh.truth == "nonclade"]
    assert (clade.label == "phylo_consistent").mean() > 0.75
    assert (clade.label == "phylo_violating").mean() < 0.03
    assert (nonc.label == "phylo_violating").mean() > 0.85
    assert "tprt_score" in fit.columns                       # annotation joined
    calibration.main(["--fit", str(out), "--truth", str(sim / "truth.tsv"), "--out", str(tmp_path / "cal")])
    discrimination.main(["--fit", str(out / "phylo_fit.tsv"), "--out", str(tmp_path / "disc"),
                         "--truth", str(sim / "truth.tsv"), "--bootstrap", "50"])
    rep = (tmp_path / "disc" / "report.md").read_text()
    assert "tprt_score" in rep and "private" in rep
    feats = pd.read_csv(tmp_path / "disc" / "features.tsv", sep="\t").set_index("feature")
    assert feats.loc["tprt_score", "auc"] > 0.6               # real insertions score higher by design


def _to_v2_dir(src, dst, drop_every=7):
    """Rewrite legacy per-colony files as genotype2 numeric files with the same n_alt / n_ref
    (posteriors are placeholders); every `drop_every`-th row of the first colony gets status
    no_reads (empty numerics) and is also zeroed in a legacy copy, so both fits see the same data."""
    import gzip
    from tools import genotype2_io as GIO
    os.makedirs(dst / "v2", exist_ok=True)
    os.makedirs(dst / "legacy", exist_ok=True)
    for k, f in enumerate(sorted(os.listdir(src))):
        with gzip.open(src / f, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            rows = [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh]
        with gzip.open(dst / "v2" / f, "wt") as v2, gzip.open(dst / "legacy" / f, "wt") as lg:
            v2.write("\t".join(GIO.V2_HEADER) + "\n")
            lg.write("\t".join(head) + "\n")
            for i, r in enumerate(rows):
                if k == 0 and i % drop_every == 0:
                    v2.write("\t".join([r["insertion"], "TSD", "no_reads", "0"] + [""] * 17) + "\n")
                    r = dict(r, n_alt="0", n_ref="0", genotype="no-coverage")
                else:
                    a, n = int(r["n_alt"]), int(r["n_ref"])
                    v2.write("\t".join([r["insertion"], "TSD", "ok", str(a + n + 3), str(a), str(n), "3", "0",
                                        "0", str(a), str(a), "0.5", "0.3", "0.4", "0.3", "5", "4", "5",
                                        "4", "10", "10"]) + "\n")
                lg.write("\t".join(r[h] for h in head) + "\n")


def test_tree_fit_reads_genotype2_numeric_dir_and_matrix(tmp_path):
    """v2 per-colony files feed the same read-vote model (n_alt / n_ref; status != ok = missing):
    the fit equals a legacy fit of the same counts. A numeric --genotypes matrix adds
    matrix_n_carriers / matrix_n_wt and, with <P>.joint.tsv beside it, joint_* columns + the
    cross-tab in summary.md, and never writes passes_combine_genotypes (it is not a gate)."""
    import gzip
    import pandas as pd
    from tools.phylo import simulate_counts
    sim = tmp_path / "sim"
    simulate_counts.main(["--tree", "random:6", "--out", str(sim), "--seed", "3", "--n-clade", "60",
                          "--n-germline", "30", "--n-nonclade", "15", "--n-noise", "10", "--n-chrx", "0"])
    _to_v2_dir(sim / "genotypes", tmp_path)
    loci = pd.read_csv(sim / "truth.tsv", sep="\t")["locus"].tolist()
    cols = sorted(f.split(".")[0] for f in os.listdir(sim / "genotypes"))
    with gzip.open(tmp_path / "P.genotypes.csv.gz", "wt") as fh:
        fh.write(";" + ";".join(cols) + "\n")
        for i, l in enumerate(loci):
            fh.write(l + ";" + ";".join("" if (j == 0 and i % 5 == 0) else ("0.9500" if (i + j) % 3 else "0.0100")
                                        for j in range(len(cols))) + "\n")
    with open(tmp_path / "P.joint.tsv", "w") as fh:
        fh.write("locus\tlocus_kind\tn_colonies_data\tbest\tcarriers\tn_carriers\tpost_best\tbest_tree\t"
                 "best_tree_clade\tpost_tree\tlog10_bf_tree\tlog10_L_best_tree\tlog10_L_root\t"
                 "log10_L_indep\tlog10_L_noise\n")
        for i, l in enumerate(loci):
            best, car = (("ROOT", ",".join(cols)) if i % 3 == 0 else (("NOISE", "") if i % 3 == 1 else (cols[0], cols[0])))
            fh.write(f"{l}\tTSD\t{len(cols)}\t{best}\t{car}\t{len(car.split(',')) if car else 0}\t0.9\t"
                     f"ROOT\t{','.join(cols)}\t0.9\t1.5\t-1\t-1\t-2\t-3\n")
    common = ["--tree", str(sim / "tree.nwk"), "--sex", "F", "--bootstrap", "0", "--no-cells"]
    TF.main(["--genotype-dir", str(tmp_path / "legacy"), "--out", str(tmp_path / "fit_leg")] + common)
    TF.main(["--genotypes", str(tmp_path / "P.genotypes.csv.gz"), "--genotype-dir", str(tmp_path / "v2"),
             "--out", str(tmp_path / "fit_v2")] + common)
    a = pd.read_csv(tmp_path / "fit_leg" / "phylo_fit.tsv", sep="\t").set_index("locus")
    b = pd.read_csv(tmp_path / "fit_v2" / "phylo_fit.tsv", sep="\t").set_index("locus")
    for c in ("class", "label", "best_branch", "total_alt", "total_ref", "n_colonies_data", "log10_bf_tree"):
        assert (a[c] == b.loc[a.index, c]).all(), c
    assert "passes_combine_genotypes" not in b.columns
    for c in ("matrix_n_carriers", "matrix_n_wt", "joint_best", "joint_class", "joint_n_carriers",
              "joint_log10_bf_tree"):
        assert c in b.columns, c
    assert set(b["joint_class"]) == {"ROOT", "NOISE", "private"}
    assert (b["matrix_n_carriers"] + b["matrix_n_wt"]).max() == len(cols)
    assert "Cross-check: tree_fit class x genotype2 joint step" in (tmp_path / "fit_v2" / "summary.md").read_text()
    # a directory mixing both formats is refused
    import shutil
    shutil.copy(tmp_path / "legacy" / sorted(os.listdir(tmp_path / "legacy"))[0], tmp_path / "v2" / "Zmixed.txt.gz")
    with pytest.raises(ValueError, match="mixed genotype formats"):
        G.read_genotype_dir(str(tmp_path / "v2"))


def test_simlib_tree_design_nonclade_is_never_a_clade():
    sys.path.insert(0, os.path.join(REPO, "test"))
    from simlib import phylo as PH
    rng = random.Random(9)
    d = PH.build_design("random:9", rng)
    assert d.n == 9 and len(d.clades) == 2 * 9 - 2
    cs = d.clade_set()
    for _ in range(200):
        s = PH.draw_nonclade(rng, d)
        assert s not in cs and 2 <= len(s) <= 8
    for _ in range(200):
        bid, s = PH.draw_branch(rng, d, root_frac=0.2)
        assert bid == "ROOT" and len(s) == 9 or s in cs
    assert all(0.7 <= p <= 1.0 for p in d.purity)
