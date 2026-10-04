"""Per locus x colony genotype likelihoods from the genotyper's read votes (n_alt, n_ref).

Read-vote model (see plans/tprt_hallmarks/PHYLO_EVAL.md for the derivation)
---------------------------------------------------------------------------
The Rust genotyper casts one vote per informative read: alt (crosses an insertion junction) or
ref (spans the reference junction). How many votes one haplotype copy yields depends on the
locus kind (a TSD locus has two alt junctions per reference span; far pairs and one-sided loci
are corrected by `halve_single_junction_ref` / `dup_ref_discount_min_span`), so a germline het
does not read exactly 0.5. We therefore model the per-haplotype vote odds

    K = (alt votes per alt-haplotype copy) / (ref votes per ref-haplotype copy)

and estimate it from loci that are clearly germline heterozygous (present in every colony at a
het-like allele fraction): K_kind per locus kind, times a per-colony allelic-balance factor b_c
(library / mapping bias of colony c), K_ic = K_kind(i) * b_c.

Expected alt-vote fraction of colony c (purity p_c = fraction of the colony's cells that
descend from the founding cell; the rest carry none of the colony's somatic variants):

    somatic clonal het, diploid:   f1 = (p/2) K / ((p/2) K + 1 - p/2)
    somatic, haploid (male X/Y):   f1 = p K / (p K + 1 - p)
    germline het (all cells):      fg = K / (K + 1)
    germline hom / haploid:        fh = 1 - eps_c
    absent:                        f0 = eps_c   (background alt-vote rate: mapping, slippage)

Counts are beta-binomial with intra-class correlation rho (rho1 present, rho0 absent; rho = 0
is the binomial). The present likelihood has a small "subclonal" component (weight lam) with
the fraction uniform on (0, f1), so a colony that carries the insertion at a lower clonal
fraction is not an outright contradiction:

    P(a | n, present) = (1 - lam) BB(a | n, f1, rho1) + lam * I_{f1}(a + 1, n - a + 1) / ((n + 1) f1)

where I is the regularised incomplete beta (the closed form of the integral over u ~ U(0, f1)).

Dropout: P(n_alt = 0 | present, n) is reported for a clonal carrier (beta-binomial at f1). A clade
member with 0 alt reads of 2 informative reads is expected (P ~ 0.2-0.4); with 40 reads it is not
(P ~ 1e-10). The violation p-values (`p_chance`) use the full present likelihood incl. the
subclonal component, which floors them at ~ lam / ((n + 1) f1): conservative by design.
"""
from __future__ import annotations

import gzip
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import betainc, betaln, gammaln

F_MIN, F_MAX = 1e-6, 1.0 - 1e-6
ONESIDE_TOKEN = "oneside_"
LOCUS_KINDS = ("TSD", "TSD_DELETION", "BLUNT", "L1_MED_DELETION", "L1_MED_DUPLICATION",
               "ONE_SIDED", "OTHER")
SEX_CONTIGS = {"chrX", "X", "chrY", "Y"}


# ------------------------------------------------------------------ locus names / reading


def locus_kind(name: str) -> str:
    """Geometry class of a locus name; mirrors src/combine_genotypes.locus_kind (copied: that
    module imports the pipeline CONFIG at import time)."""
    contig, sep, pos = str(name).rpartition(":")
    left, sep2, right = pos.rpartition("-")
    if not sep or not sep2:
        return "OTHER"
    if left.startswith(ONESIDE_TOKEN) != right.startswith(ONESIDE_TOKEN):
        return "ONE_SIDED"
    try:
        gap = int(right) - int(left)
    except ValueError:
        return "OTHER"
    if gap < -30:
        return "L1_MED_DELETION"
    if gap < 0:
        return "TSD_DELETION"
    if gap <= 1:
        return "BLUNT"
    if gap <= 40:
        return "TSD"
    return "L1_MED_DUPLICATION"


def locus_contig(name: str) -> str:
    return str(name).rpartition(":")[0]


def sample_stem(path: str) -> str:
    """Colony name of a per-colony genotype file (`genotypes/<SAMPLE>.txt.gz`)."""
    base = os.path.basename(path)
    return re.sub(r"(\.genotypes?)?(\.(txt|csv|tsv))?(\.gz)?$", "", base)


@dataclass
class Counts:
    """n_alt / n_ref per locus x colony (0 where a colony's file lacks the locus)."""
    loci: List[str]
    colonies: List[str]
    alt: np.ndarray            # L x C int
    ref: np.ndarray            # L x C int
    calls: Optional[pd.DataFrame] = None   # genotype strings (L x C), informational

    @property
    def depth(self) -> np.ndarray:
        return self.alt + self.ref

    def subset_colonies(self, cols: Sequence[str]) -> "Counts":
        ix = [self.colonies.index(c) for c in cols]
        calls = self.calls[list(cols)] if self.calls is not None else None
        return Counts(self.loci, list(cols), self.alt[:, ix], self.ref[:, ix], calls)


def read_genotype_file(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"insertion": str}, compression="infer")
    if "n_alt" not in df.columns or "n_ref" not in df.columns:
        raise ValueError(f"{path}: genotype file has no n_alt / n_ref columns (genotyper too old)")
    return df.drop_duplicates("insertion").set_index("insertion")


def list_genotype_files(directory: str) -> Dict[str, str]:
    out = {}
    for f in sorted(os.listdir(directory)):
        if ".tmp" in f or f.endswith((".log", ".md", ".json")):
            continue
        if not re.search(r"\.(txt|tsv)(\.gz)?$", f):
            continue
        out[sample_stem(f)] = os.path.join(directory, f)
    return out


def read_genotype_dir(directory: str, colonies: Optional[Sequence[str]] = None,
                      loci: Optional[Sequence[str]] = None) -> Counts:
    files = list_genotype_files(directory)
    if colonies is not None:
        missing = [c for c in colonies if c not in files]
        files = {c: files[c] for c in colonies if c in files}
        if missing:
            print(f"  note: {len(missing)} colonies have no genotype file: "
                  f"{', '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}")
    if not files:
        raise ValueError(f"no per-colony genotype files in {directory}")
    tables = {c: read_genotype_file(p) for c, p in files.items()}
    if loci is None:
        seen, order = set(), []
        for t in tables.values():
            for n in t.index:
                if n not in seen:
                    seen.add(n)
                    order.append(n)
        loci = order
    loci = list(loci)
    cols = list(tables)
    alt = np.zeros((len(loci), len(cols)), dtype=np.int64)
    ref = np.zeros_like(alt)
    calls = pd.DataFrame(index=loci, columns=cols, dtype=object)
    for j, c in enumerate(cols):
        t = tables[c].reindex(loci)
        alt[:, j] = t["n_alt"].fillna(0).astype(np.int64).values
        ref[:, j] = t["n_ref"].fillna(0).astype(np.int64).values
        if "genotype" in t.columns:
            calls[c] = t["genotype"].values
    return Counts(loci, cols, np.maximum(alt, 0), np.maximum(ref, 0), calls)


def read_genotypes_csv(path: str) -> pd.DataFrame:
    """`<patient>.genotypes.csv.gz` (combine_genotypes output: `;`-separated calls, loci that
    passed its gates). Carries calls only, no read counts."""
    return pd.read_csv(path, sep=";", index_col=0, compression="infer")


# ------------------------------------------------------------------ distributions


def _clip(f):
    return np.clip(f, F_MIN, F_MAX)


def log_choose(n, a):
    return gammaln(n + 1.0) - gammaln(a + 1.0) - gammaln(n - a + 1.0)


def bb_logpmf(a, n, f, rho):
    """Beta-binomial log pmf with mean f and intra-class correlation rho (binomial at rho=0).
    Broadcasts over a, n, f."""
    a = np.asarray(a, float)
    n = np.asarray(n, float)
    f = _clip(np.asarray(f, float))
    if rho <= 1e-9:
        return log_choose(n, a) + a * np.log(f) + (n - a) * np.log1p(-f)
    s = (1.0 - rho) / rho
    al, be = f * s, (1.0 - f) * s
    return log_choose(n, a) + betaln(a + al, n - a + be) - betaln(al, be)


def subclonal_pmf(a, n, f1):
    """P(a | n, fraction ~ U(0, f1)) = I_{f1}(a+1, n-a+1) / ((n+1) f1)."""
    a = np.asarray(a, float)
    n = np.asarray(n, float)
    f1 = _clip(np.asarray(f1, float))
    return betainc(a + 1.0, n - a + 1.0, f1) / ((n + 1.0) * f1)


def present_logpmf(a, n, f1, rho1, lam):
    lp = bb_logpmf(a, n, f1, rho1)
    if lam <= 0:
        return lp
    with np.errstate(divide="ignore"):
        return np.logaddexp(np.log1p(-lam) + lp, np.log(lam) + np.log(np.maximum(subclonal_pmf(a, n, f1), 1e-300)))


def present_cdf(a, n, f1, rho1, lam, cap: int = 60):
    """P(A <= a | n, present). Exact summation for a <= cap, else 1.0 (only small alt counts
    of clade members are ever tested for dropout)."""
    a = np.asarray(a, np.int64)
    n = np.asarray(n, np.int64)
    a, n, f1 = np.broadcast_arrays(a, n, np.asarray(f1, float))
    out = np.ones(a.shape)
    m = a <= cap
    if not m.any():
        return out
    am, nm, fm = a[m], n[m], f1[m]
    acc = np.zeros(am.shape)
    for j in range(int(am.max()) + 1):
        sel = am >= j
        acc[sel] += np.exp(present_logpmf(np.full(sel.sum(), j), nm[sel], fm[sel], rho1, lam))
    out[m] = np.minimum(acc, 1.0)
    return out


def absent_sf(a, n, eps, rho0):
    """P(A >= a | n, absent)."""
    a = np.asarray(a, float)
    n = np.asarray(n, float)
    eps = _clip(np.asarray(eps, float))
    a, n, eps = np.broadcast_arrays(a, n, eps)
    if rho0 <= 1e-9:
        p = stats.binom.sf(a - 1, n, eps)
    else:
        s = (1.0 - rho0) / rho0
        p = stats.betabinom.sf(a - 1, n, eps * s, (1.0 - eps) * s)
    return np.where(a <= 0, 1.0, p)


# ------------------------------------------------------------------ model parameters


@dataclass
class Params:
    colonies: List[str]
    K_kind: Dict[str, float] = field(default_factory=dict)
    K_global: float = 1.0
    b: np.ndarray = None           # per-colony allelic-balance factor
    purity: np.ndarray = None
    eps: np.ndarray = None
    rho0: float = 0.01
    rho1: float = 0.02
    lam: float = 0.01
    noise_a: float = 0.5           # locus-noise hypothesis: eps_locus ~ Beta(noise_a, noise_b)
    noise_b: float = 2.0
    sex: str = "F"
    # bookkeeping for the report
    n_germline_het: int = 0
    n_purity_cells: np.ndarray = None
    n_eps_cells: np.ndarray = None

    def __post_init__(self):
        C = len(self.colonies)
        if self.b is None:
            self.b = np.ones(C)
        if self.purity is None:
            self.purity = np.full(C, 0.85)
        if self.eps is None:
            self.eps = np.full(C, 0.003)
        if self.n_purity_cells is None:
            self.n_purity_cells = np.zeros(C, int)
        if self.n_eps_cells is None:
            self.n_eps_cells = np.zeros(C, int)

    def K_matrix(self, kinds: Sequence[str]) -> np.ndarray:
        k = np.array([self.K_kind.get(x, self.K_global) for x in kinds], float)
        return k[:, None] * self.b[None, :]


def f_present(K, purity, haploid):
    """Expected alt-vote fraction of a somatic clonal insertion (see module docstring)."""
    p = np.asarray(purity, float)
    h = np.where(np.asarray(haploid, bool), p, p / 2.0)     # alt-haplotype fraction (broadcast)
    return _clip(h * K / (h * K + 1.0 - h))


def f_germline_het(K):
    return _clip(K / (K + 1.0))


@dataclass
class CellLik:
    """Per cell (locus x colony) log-likelihoods; masked cells (no informative reads, or below
    --min-depth) are 0 in every hypothesis."""
    l1: np.ndarray       # somatic present
    l0: np.ndarray       # absent
    lhet: np.ndarray     # germline het (-inf on haploid loci)
    lhom: np.ndarray     # germline hom / hemizygous
    f1: np.ndarray
    fg: np.ndarray
    mask: np.ndarray     # True = informative


def cell_likelihoods(alt, ref, kinds, haploid, P: Params, min_depth: int = 1) -> CellLik:
    alt = np.asarray(alt)
    n = alt + np.asarray(ref)
    mask = n >= max(1, min_depth)
    K = P.K_matrix(kinds)
    hap = np.asarray(haploid, bool)[:, None]
    f1 = f_present(K, P.purity[None, :], hap)
    fg = f_germline_het(K)
    eps = P.eps[None, :]
    l1 = present_logpmf(alt, n, f1, P.rho1, P.lam)
    l0 = bb_logpmf(alt, n, eps, P.rho0)
    lhet = np.where(hap, -np.inf, bb_logpmf(alt, n, fg, P.rho1))
    lhom = bb_logpmf(alt, n, 1.0 - eps, P.rho1)
    z = ~mask
    for x in (l1, l0, lhom):
        x[z] = 0.0
    lhet = np.where(z, 0.0, lhet)
    return CellLik(l1, l0, lhet, lhom, f1, fg, mask)


def dropout_prob(n, f1, P: Params, subclonal: bool = False):
    """P(n_alt = 0 | present, n informative reads) of a CLONAL carrier (beta-binomial at f1).
    subclonal=True adds the lam component (the full present likelihood the fit uses; it floors
    the probability at ~ lam / ((n + 1) f1))."""
    z = np.zeros_like(np.asarray(n, float))
    if subclonal:
        return np.exp(present_logpmf(z, n, f1, P.rho1, P.lam))
    return np.exp(bb_logpmf(z, n, f1, P.rho1))


# ------------------------------------------------------------------ estimators


def germline_het_candidates(alt, ref, autosomal, min_cols: int = 2):
    """Loci that are clearly germline heterozygous: carried by every colony with data.

    pooled alt fraction in [0.2, 0.8]; >= max(min_cols, half the colonies) colonies with >= 3
    informative reads, ALL of them with >= 1 alt read (a somatic insertion near the root, absent
    from one colony, would otherwise drag K towards its purity-diluted fraction); >= 10 reads in
    total; autosomal. Returns a boolean locus mask."""
    n = alt + ref
    tot_a, tot_n = alt.sum(1), n.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = tot_a / tot_n
    cov3 = n >= 3
    ncov = cov3.sum(1)
    frac_alt = np.where(ncov > 0, (cov3 & (alt >= 1)).sum(1) / np.maximum(ncov, 1), 0)
    deep_zero = ((n >= 8) & (alt == 0)).any(1)
    return (autosomal & (tot_n >= 10) & (v >= 0.2) & (v <= 0.8)
            & (ncov >= max(min_cols, alt.shape[1] // 2)) & (frac_alt >= 1.0) & ~deep_zero)


def estimate_balance(alt, ref, kinds, autosomal, P: Params, min_loci_kind: int = 5,
                     rho_grid=(0.0, 0.003, 0.01, 0.02, 0.05, 0.1, 0.2), cand=None):
    """K_kind, b_c and rho1 from germline-het candidates.

    K_kind = (sum alt + 0.5) / (sum ref + 0.5) over the kind's candidates (>= min_loci_kind,
    else the pooled K over all kinds; K = 1 if there are fewer than 5 candidates at all).
    b_c = (sum_i alt_ic + 5) / (sum_i K_i ref_ic + 5): colony c's alt odds relative to what the
    kind odds predict (5 pseudo-votes shrink sparse colonies to 1). rho1 = ML on a grid over
    the candidates' cells given f = K b / (K b + 1). `cand` overrides the candidate set (tree_fit
    passes the loci it fitted as germline het after the first round)."""
    if cand is None or np.sum(cand) < 5:
        cand = germline_het_candidates(alt, ref, autosomal)
    P.n_germline_het = int(cand.sum())
    kinds = np.asarray(kinds)
    if cand.sum() < 5:
        P.K_kind, P.K_global = {}, 1.0
        P.b = np.ones(alt.shape[1])
        return cand
    a, r = alt[cand], ref[cand]
    P.K_global = float((a.sum() + 0.5) / (r.sum() + 0.5))
    P.K_kind = {}
    for k in LOCUS_KINDS:
        sel = kinds[cand] == k
        if sel.sum() >= min_loci_kind:
            P.K_kind[k] = float((a[sel].sum() + 0.5) / (r[sel].sum() + 0.5))
    Kc = np.array([P.K_kind.get(k, P.K_global) for k in kinds[cand]])
    P.b = (a.sum(0) + 5.0) / ((Kc[:, None] * r).sum(0) + 5.0)
    K = Kc[:, None] * P.b[None, :]
    f = f_germline_het(K)
    n = a + r
    m = n > 0
    best = max(rho_grid, key=lambda rho: bb_logpmf(a[m], n[m], f[m], rho).sum())
    P.rho1 = float(max(best, 1e-4))
    return cand


def ml_purity(a, n, K, haploid, P: Params, grid=None):
    """ML clonal purity of one colony from its cells in confidently placed clades."""
    if grid is None:
        grid = np.linspace(0.05, 1.0, 96)
    a = np.asarray(a, float)
    n = np.asarray(n, float)
    K = np.asarray(K, float)
    hap = np.asarray(haploid, bool)
    best, bl = None, -np.inf
    for p in grid:
        h = np.where(hap, p, p / 2.0)
        f1 = _clip(h * K / (h * K + 1.0 - h))
        ll = present_logpmf(a, n, f1, P.rho1, P.lam).sum()
        if ll > bl:
            best, bl = p, ll
    return float(best)


def ml_rho(a, n, f, grid=(0.0, 0.001, 0.003, 0.01, 0.03, 0.1)):
    a = np.asarray(a, float)
    n = np.asarray(n, float)
    f = np.asarray(f, float)
    if a.size == 0:
        return None
    return float(max(grid, key=lambda rho: bb_logpmf(a, n, f, rho).sum()))


def infer_sex(loci, alt, ref):
    """'M' / 'F' / None from germline-like chrX/chrY loci (carried by >= 90 % of colonies with
    data, >= 20 reads pooled): >= 2 het-like X loci (pooled fraction 0.2-0.8) -> F; else >= 3
    hom-like X loci (>= 0.85) or any carried chrY locus -> M; else undecided."""
    contigs = np.array([locus_contig(x) for x in loci])
    n = alt + ref
    tot_a, tot_n = alt.sum(1), n.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.where(tot_n > 0, tot_a / np.maximum(tot_n, 1), 0)
    cov = n >= 3
    carried = np.where(cov.sum(1) > 0, (cov & (alt >= 1)).sum(1) / np.maximum(cov.sum(1), 1), 0) >= 0.9
    gl = carried & (tot_n >= 20)
    isx = np.isin(contigs, ["chrX", "X"])
    isy = np.isin(contigs, ["chrY", "Y"])
    het = int((gl & isx & (v >= 0.2) & (v <= 0.8)).sum())
    hom = int((gl & isx & (v >= 0.85)).sum())
    y = int((isy & (alt.sum(1) >= 4) & ((alt >= 2).sum(1) >= 2)).sum())
    if het >= 2:
        return "F", dict(x_het=het, x_hom=hom, y_carried=y)
    if hom >= 3 or y >= 1:
        return "M", dict(x_het=het, x_hom=hom, y_carried=y)
    return None, dict(x_het=het, x_hom=hom, y_carried=y)
