//! Joint phylogenetic genotyping across the colonies of one patient (SPEC "Joint"). Owner: E.
//! Owns newick.rs too.
//!
//! Port of tools/phylo/tree_fit.py + genotype_likelihood.py. Per locus i and colony c the data
//! give two likelihoods, `l1 = log P(d_c | present)` and `l0 = log P(d_c | absent)`; missing
//! cells (no row, `no-coverage` / `high-coverage` / `error`, no informative reads) are 0 in both.
//!
//! Two likelihood sources, chosen per run:
//! * **PL mode** (every file has `pl_ref pl_het pl_hom`): `L_g = 10^(-PL_g/10)` normalised to the
//!   max; `present = ½(L_het + L_hom)`, `absent = L_ref`. ROOT = every colony present; NOISE is
//!   approximated by "absent everywhere" (Σ l0): the PLs already integrate the background alt
//!   rate, so a shared per-locus alt fraction cannot be refitted from them.
//! * **legacy mode** (any file without PL columns): the beta-binomial read-vote model of
//!   genotype_likelihood.py on `n_alt`/`n_ref`, with tree_fit.py's full parameter estimation
//!   (K per locus kind + per-colony allelic balance + ρ1 from germline-het loci; three rounds of
//!   purity / ε / ρ0 updates from confidently placed loci; sex inference for chrX/Y ploidy).
//!   ROOT is tree_fit's mixture (germline het / germline hom ~ Beta(8,1) / somatic-in-all) and
//!   NOISE the Beta(0.5, 4) marginal.
//!
//! Hypotheses: every node b of the tree (ROOT and each branch, clade = tips below), with prior
//! `root_prior` for ROOT and `1 - root_prior` spread over branches by length (floored at 1 % of
//! the mean positive length) or uniformly; INDEP (independent presence, π ~ U(0,1), exact DP over
//! the carrier count) and NOISE. `log10_bf_tree = log10(Σ_b prior_b L_b / mean(L_indep, L_noise))`.
//! The posterior over ALL hypotheses uses model weights ½ (tree) / ¼ (INDEP) / ¼ (NOISE), so the
//! posterior odds of "a tree event" vs "not" equal the BF. Per colony
//! `P(carrier) = Σ_{b ∋ c} post_b + post_INDEP · P(c present | INDEP, data)` (exact, O(C²)).

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::newick::Tree;
use crate::types::{
    GT_ERROR, GT_HETEROZYGOUS, GT_HIGH_COVERAGE, GT_HOMOZYGOUS, GT_INSERTION, GT_INSERTION_UNCERTAIN,
    GT_NO_COVERAGE, GT_WILDTYPE, GT_WILDTYPE_UNCERTAIN,
};

pub struct JointArgs {
    pub tree: String,
    /// per-colony genotype files (stem = colony id = tree tip label)
    pub genotype_files: Vec<String>,
    pub out_tsv: String,
    pub out_matrix: String,
    pub root_prior: f64,
    /// "length" | "uniform"
    pub branch_prior: String,
}

const LN10: f64 = std::f64::consts::LN_10;
const LN2: f64 = std::f64::consts::LN_2;
const NEG_INF: f64 = f64::NEG_INFINITY;
/// model weights of the three hypothesis families in the overall posterior
const LN_W_TREE: f64 = -std::f64::consts::LN_2; // ln 1/2
const LN_W_ALT: f64 = -2.0 * std::f64::consts::LN_2; // ln 1/4 each (INDEP, NOISE)
/// tree_fit.py Options.floor_frac
const FLOOR_FRAC: f64 = 0.01;

fn bad(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

// ------------------------------------------------------------------ special functions

/// ln Γ(x), Lanczos (g = 7, n = 9), reflection below 0.5. Relative error ~1e-15.
fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        let s = (std::f64::consts::PI * x).sin().abs();
        return (std::f64::consts::PI / s).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = G[0];
    for (i, g) in G.iter().enumerate().skip(1) {
        a += g / (x + i as f64);
    }
    let t = x + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

fn betaln(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

fn log_choose(n: f64, a: f64) -> f64 {
    ln_gamma(n + 1.0) - ln_gamma(a + 1.0) - ln_gamma(n - a + 1.0)
}

/// Continued fraction of the incomplete beta (modified Lentz).
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    const FPMIN: f64 = 1e-300;
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..2000 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-16 {
            break;
        }
    }
    h
}

/// Regularised incomplete beta I_x(a, b) (scipy.special.betainc).
fn betainc(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let lbt = ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (-x).ln_1p();
    if x < (a + 1.0) / (a + b + 2.0) {
        lbt.exp() * betacf(a, b, x) / a
    } else {
        1.0 - lbt.exp() * betacf(b, a, 1.0 - x) / b
    }
}

fn lae(a: f64, b: f64) -> f64 {
    if a == NEG_INF {
        return b;
    }
    if b == NEG_INF {
        return a;
    }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

fn lse(xs: impl IntoIterator<Item = f64>) -> f64 {
    let v: Vec<f64> = xs.into_iter().collect();
    let m = v.iter().cloned().fold(NEG_INF, f64::max);
    if m == NEG_INF || !m.is_finite() {
        return m;
    }
    m + v.iter().map(|x| (x - m).exp()).sum::<f64>().ln()
}

/// index of the first maximum (numpy argmax semantics; NaN never wins)
fn argmax(xs: &[f64]) -> usize {
    let mut best = 0;
    for (i, &x) in xs.iter().enumerate() {
        if x > xs[best] || (xs[best].is_nan() && !x.is_nan()) {
            best = i;
        }
    }
    best
}

// ------------------------------------------------------------------ read-vote distributions
// (tools/phylo/genotype_likelihood.py)

const F_MIN: f64 = 1e-6;
const F_MAX: f64 = 1.0 - 1e-6;

fn clipf(f: f64) -> f64 {
    f.clamp(F_MIN, F_MAX)
}

/// Beta-binomial log pmf with mean f and intra-class correlation rho (binomial at rho ≈ 0).
fn bb_logpmf(a: f64, n: f64, f: f64, rho: f64) -> f64 {
    let f = clipf(f);
    if rho <= 1e-9 {
        return log_choose(n, a) + a * f.ln() + (n - a) * (-f).ln_1p();
    }
    let s = (1.0 - rho) / rho;
    let (al, be) = (f * s, (1.0 - f) * s);
    log_choose(n, a) + betaln(a + al, n - a + be) - betaln(al, be)
}

/// P(a | n, fraction ~ U(0, f1)) = I_{f1}(a+1, n-a+1) / ((n+1) f1)
fn subclonal_pmf(a: f64, n: f64, f1: f64) -> f64 {
    let f1 = clipf(f1);
    betainc(a + 1.0, n - a + 1.0, f1) / ((n + 1.0) * f1)
}

fn present_logpmf(a: f64, n: f64, f1: f64, rho1: f64, lam: f64) -> f64 {
    let lp = bb_logpmf(a, n, f1, rho1);
    if lam <= 0.0 {
        return lp;
    }
    lae((-lam).ln_1p() + lp, lam.ln() + subclonal_pmf(a, n, f1).max(1e-300).ln())
}

fn f_present(k: f64, purity: f64, haploid: bool) -> f64 {
    let h = if haploid { purity } else { purity / 2.0 };
    clipf(h * k / (h * k + 1.0 - h))
}

fn f_germline_het(k: f64) -> f64 {
    clipf(k / (k + 1.0))
}

// ------------------------------------------------------------------ locus names

const ONESIDE_TOKEN: &str = "oneside_";
const LOCUS_KINDS: [&str; 7] = ["TSD", "TSD_DELETION", "BLUNT", "L1_MED_DELETION", "L1_MED_DUPLICATION", "ONE_SIDED", "OTHER"];

/// genotype_likelihood.locus_kind
fn locus_kind(name: &str) -> &'static str {
    let Some((_, pos)) = name.rsplit_once(':') else { return "OTHER" };
    let Some((left, right)) = pos.rsplit_once('-') else { return "OTHER" };
    if left.starts_with(ONESIDE_TOKEN) != right.starts_with(ONESIDE_TOKEN) {
        return "ONE_SIDED";
    }
    let (Ok(l), Ok(r)) = (left.parse::<i64>(), right.parse::<i64>()) else { return "OTHER" };
    match r - l {
        g if g < -30 => "L1_MED_DELETION",
        g if g < 0 => "TSD_DELETION",
        g if g <= 1 => "BLUNT",
        g if g <= 40 => "TSD",
        _ => "L1_MED_DUPLICATION",
    }
}

fn locus_contig(name: &str) -> &str {
    name.rsplit_once(':').map(|(c, _)| c).unwrap_or("")
}

// ------------------------------------------------------------------ input

#[derive(Clone, Debug)]
struct Row {
    genotype: String,
    n_alt: i64,
    n_ref: i64,
    /// [pl_ref, pl_het, pl_hom]; None when the file has no PL columns or the values are not numbers
    pl: Option<[f64; 3]>,
}

struct ColonyTable {
    has_counts: bool,
    has_pl: bool,
    rows: HashMap<String, Row>,
    order: Vec<String>,
}

fn is_na_state(g: &str) -> bool {
    g == GT_NO_COVERAGE || g == GT_HIGH_COVERAGE || g == GT_ERROR
}

/// Colony id of a per-colony genotype file: basename minus `.gz`, `.txt|.tsv|.csv`, `.genotype(s)`.
fn colony_stem(path: &str) -> String {
    let base = std::path::Path::new(path).file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
    let mut s = base.as_str();
    s = s.strip_suffix(".gz").unwrap_or(s);
    for ext in [".txt", ".tsv", ".csv"] {
        if let Some(t) = s.strip_suffix(ext) {
            s = t;
            break;
        }
    }
    for ext in [".genotypes", ".genotype"] {
        if let Some(t) = s.strip_suffix(ext) {
            s = t;
            break;
        }
    }
    s.to_string()
}

fn open_text(path: &str) -> io::Result<Box<dyn BufRead>> {
    let mut f = File::open(path).map_err(|e| io::Error::new(e.kind(), format!("cannot open {path}: {e}")))?;
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic)?;
    drop(f);
    let f = File::open(path)?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(BufReader::new(f)))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

fn parse_count(s: Option<&str>) -> i64 {
    s.and_then(|x| x.trim().parse::<f64>().ok()).filter(|x| x.is_finite()).map(|x| x.max(0.0) as i64).unwrap_or(0)
}

fn read_colony_table(path: &str) -> io::Result<ColonyTable> {
    let mut lines = open_text(path)?.lines();
    let header = lines.next().transpose()?.ok_or_else(|| bad(format!("{path}: empty genotype file")))?;
    let cols: Vec<&str> = header.trim_end_matches(['\r', '\n']).split('\t').collect();
    let find = |name: &str| cols.iter().position(|c| *c == name);
    let i_name = find("insertion").unwrap_or(0);
    let i_gt = find("genotype").unwrap_or(1);
    let (i_alt, i_ref) = (find("n_alt"), find("n_ref"));
    let i_pl = match (find("pl_ref"), find("pl_het"), find("pl_hom")) {
        (Some(a), Some(b), Some(c)) => Some([a, b, c]),
        _ => None,
    };
    let mut t = ColonyTable {
        has_counts: i_alt.is_some() && i_ref.is_some(),
        has_pl: i_pl.is_some(),
        rows: HashMap::new(),
        order: Vec::new(),
    };
    for line in lines {
        let line = line?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let Some(name) = f.get(i_name) else { continue };
        if t.rows.contains_key(*name) {
            continue; // drop_duplicates("insertion"): first wins
        }
        let pl = i_pl.and_then(|ix| {
            let v: Vec<f64> = ix.iter().filter_map(|&k| f.get(k).and_then(|x| x.trim().parse::<f64>().ok())).collect();
            (v.len() == 3 && v.iter().all(|x| x.is_finite())).then(|| [v[0], v[1], v[2]])
        });
        let row = Row {
            genotype: f.get(i_gt).map(|s| s.to_string()).unwrap_or_default(),
            n_alt: parse_count(i_alt.and_then(|k| f.get(k).copied())),
            n_ref: parse_count(i_ref.and_then(|k| f.get(k).copied())),
            pl,
        };
        t.order.push(name.to_string());
        t.rows.insert(name.to_string(), row);
    }
    Ok(t)
}

// ------------------------------------------------------------------ tree hypotheses

/// Branch prior (tree_fit.branch_log_prior), indexed by node (0 = ROOT).
fn branch_log_prior(tree: &Tree, mode: &str, root_prior: f64) -> Vec<f64> {
    let nb = tree.nodes.len();
    let lens: Vec<f64> = tree.nodes[1..].iter().map(|n| n.length).collect();
    let pos: Vec<f64> = lens.iter().cloned().filter(|&l| l > 0.0).collect();
    let w: Vec<f64> = if mode == "length" && !pos.is_empty() {
        let floor = FLOOR_FRAC * pos.iter().sum::<f64>() / pos.len() as f64;
        lens.iter().map(|&l| l.max(floor)).collect()
    } else {
        vec![1.0; nb - 1]
    };
    let tot: f64 = w.iter().sum();
    let mut lp = Vec::with_capacity(nb);
    lp.push(root_prior.ln());
    lp.extend(w.iter().map(|x| ((1.0 - root_prior) * x / tot).ln()));
    lp
}

/// Fixed per-tree quantities shared by every locus.
struct Ctx<'a> {
    tree: &'a Tree,
    ncol: usize,
    logprior: Vec<f64>,
    /// betaln(k+1, C-k+1), k = 0..=C (INDEP integral of π^k (1-π)^(C-k))
    bl_full: Vec<f64>,
    /// betaln(i+2, C-i), i = 0..C (INDEP per-colony presence functional, see `indep`)
    bl_v: Vec<f64>,
    /// node indices on the path tip -> root, per colony (the hypotheses that contain it)
    paths: Vec<Vec<usize>>,
    /// in_clade[node * C + c]
    in_clade: Vec<bool>,
}

impl<'a> Ctx<'a> {
    fn new(tree: &'a Tree, branch_prior: &str, root_prior: f64) -> Ctx<'a> {
        let c = tree.tips.len();
        let cf = c as f64;
        let bl_full = (0..=c).map(|k| betaln(k as f64 + 1.0, cf - k as f64 + 1.0)).collect();
        let bl_v = (0..c).map(|i| betaln(i as f64 + 2.0, cf - i as f64)).collect();
        let paths = (0..c)
            .map(|t| {
                let mut p = vec![tree.tip_node[t]];
                while let Some(up) = tree.nodes[*p.last().unwrap()].parent {
                    p.push(up);
                }
                p
            })
            .collect();
        let mut in_clade = vec![false; tree.nodes.len() * c];
        for (b, n) in tree.nodes.iter().enumerate() {
            for &t in &n.clade {
                in_clade[b * c + t] = true;
            }
        }
        Ctx { tree, ncol: c, logprior: branch_log_prior(tree, branch_prior, root_prior), bl_full, bl_v, paths, in_clade }
    }
}

/// INDEP: log ∫_0^1 Π_c [π P1_c + (1-π) P0_c] dπ, exactly, by DP over the carrier count
/// (tree_fit.log_indep); with `want_colony`, also log P(z_c = 1, data | INDEP) per colony.
///
/// Per colony: P(z_c=1, d) = P1_c Σ_i pre_c[i] V_c(i), where pre_c = Bernstein-like coefficients
/// of the colonies before c and V_c(i) = Σ_j suf_c[j] B(i+j+2, C-i-j) folds in the colonies after
/// c. V obeys V_{c-1}(i) = P1_c V_c(i+1) + P0_c V_c(i) with V_{C-1}(i) = B(i+2, C-i), so all
/// colonies cost O(C²) together.
fn indep(l1: &[f64], l0: &[f64], bl_full: &[f64], bl_v: &[f64], want_colony: bool) -> (f64, Vec<f64>) {
    let c = l1.len();
    let mut pre: Vec<Vec<f64>> = Vec::with_capacity(c + 1);
    pre.push(vec![0.0]);
    for m in 0..c {
        let prev = &pre[m];
        let mut row = vec![NEG_INF; m + 2];
        for (k, r) in row.iter_mut().enumerate() {
            let a = if k <= m { prev[k] + l0[m] } else { NEG_INF };
            let b = if k >= 1 { prev[k - 1] + l1[m] } else { NEG_INF };
            *r = lae(a, b);
        }
        pre.push(row);
    }
    let li = lse((0..=c).map(|k| pre[c][k] + bl_full[k]));
    if !want_colony || c == 0 {
        return (li, Vec::new());
    }
    let mut v: Vec<f64> = bl_v.to_vec();
    let mut lpres = vec![NEG_INF; c];
    for col in (0..c).rev() {
        lpres[col] = l1[col] + lse((0..=col).map(|i| pre[col][i] + v[i]));
        if col > 0 {
            v = (0..col).map(|i| lae(l1[col] + v[i + 1], l0[col] + v[i])).collect();
        }
    }
    (li, lpres)
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Best {
    Node(usize),
    Indep,
    Noise,
}

#[derive(Clone, Debug)]
struct LocusScore {
    /// log L of the ROOT hypothesis (S[0])
    s_root: f64,
    /// max over all tree hypotheses of log L (no prior)
    s_max: f64,
    /// MAP tree hypothesis (with prior), its log L and its posterior within the tree family
    best_tree: usize,
    s_best_tree: f64,
    post_within: f64,
    li: f64,
    ln: f64,
    /// log10 BF tree vs mean(INDEP, NOISE)
    bf: f64,
    best: Best,
    post_best: f64,
    post_tree: f64,
    p_carrier: Vec<f64>,
}

/// Score one locus. `l1`/`l0` per colony (missing cells 0 in both), `root_ll` = log L(ROOT),
/// `noise_ll` = log L(NOISE).
fn score_locus(ctx: &Ctx, l1: &[f64], l0: &[f64], root_ll: f64, noise_ll: f64, want_colony: bool) -> LocusScore {
    let nodes = &ctx.tree.nodes;
    let nb = nodes.len();
    let mut sum_d = vec![0.0; nb];
    for i in (0..nb).rev() {
        sum_d[i] = if nodes[i].is_tip() {
            let t = nodes[i].clade[0];
            l1[t] - l0[t]
        } else {
            nodes[i].children.iter().map(|&k| sum_d[k]).sum()
        };
    }
    let sum0: f64 = l0.iter().sum();
    let mut s = vec![0.0; nb];
    s[0] = root_ll;
    for b in 1..nb {
        s[b] = sum_d[b] + sum0;
    }
    let j: Vec<f64> = s.iter().zip(&ctx.logprior).map(|(a, b)| a + b).collect();
    let log_tree = lse(j.iter().cloned());
    let best_tree = argmax(&j);
    let (li, lpres) = indep(l1, l0, &ctx.bl_full, &ctx.bl_v, want_colony);
    let ln = noise_ll;
    let bf = (log_tree - (lae(li, ln) - LN2)) / LN10;
    let z = lse([LN_W_TREE + log_tree, LN_W_ALT + li, LN_W_ALT + ln]);
    let cand = [(Best::Node(best_tree), LN_W_TREE + j[best_tree]), (Best::Indep, LN_W_ALT + li), (Best::Noise, LN_W_ALT + ln)];
    let mut pick = cand[0];
    for x in &cand[1..] {
        if x.1 > pick.1 {
            pick = *x;
        }
    }
    let p_carrier = if want_colony {
        (0..ctx.ncol)
            .map(|c| {
                let tree_part: f64 = ctx.paths[c].iter().map(|&b| (LN_W_TREE + j[b] - z).exp()).sum();
                (tree_part + (LN_W_ALT + lpres[c] - z).exp()).min(1.0)
            })
            .collect()
    } else {
        Vec::new()
    };
    LocusScore {
        s_root: s[0],
        s_max: s.iter().cloned().fold(NEG_INF, f64::max),
        best_tree,
        s_best_tree: s[best_tree],
        post_within: (j[best_tree] - log_tree).exp(),
        li,
        ln,
        bf,
        best: pick.0,
        post_best: (pick.1 - z).exp(),
        post_tree: (LN_W_TREE + log_tree - z).exp(),
        p_carrier,
    }
}

// ------------------------------------------------------------------ PL mode

/// (log P(d|present), log P(d|absent)) from Phred-scaled genotype likelihoods.
fn pl_cell(pl: [f64; 3]) -> (f64, f64) {
    let s = -LN10 / 10.0;
    let (r, h, o) = (s * pl[0], s * pl[1], s * pl[2]);
    let m = r.max(h).max(o);
    let (r, h, o) = (r - m, h - m, o - m);
    (lae(h, o) - LN2, r)
}

// ------------------------------------------------------------------ legacy read-vote model

struct Params {
    k_kind: HashMap<&'static str, f64>,
    k_global: f64,
    b: Vec<f64>,
    purity: Vec<f64>,
    eps: Vec<f64>,
    rho0: f64,
    rho1: f64,
    lam: f64,
    noise_a: f64,
    noise_b: f64,
    n_germline_het: usize,
}

impl Params {
    fn new(c: usize) -> Params {
        Params {
            k_kind: HashMap::new(),
            k_global: 1.0,
            b: vec![1.0; c],
            purity: vec![0.85; c],
            eps: vec![0.003; c],
            rho0: 0.01,
            rho1: 0.02,
            lam: 0.01,
            noise_a: 0.5,
            noise_b: 4.0,
            n_germline_het: 0,
        }
    }
    fn k(&self, kind: &str, c: usize) -> f64 {
        self.k_kind.get(kind).copied().unwrap_or(self.k_global) * self.b[c]
    }
}

/// HOM_PRIOR of tree_fit.root_scores
const HOM_PRIOR: (f64, f64) = (8.0, 1.0);

struct Legacy {
    l: usize,
    c: usize,
    kinds: Vec<&'static str>,
    haploid: Vec<bool>,
    autosomal: Vec<bool>,
    /// flat [i * c + col]; masked cells have alt = n = 0
    alt: Vec<f64>,
    n: Vec<f64>,
}

impl Legacy {
    fn mask(&self, k: usize) -> bool {
        self.n[k] >= 1.0
    }
}

/// Per-cell likelihoods of every locus (genotype_likelihood.cell_likelihoods); masked cells 0.
struct Cells {
    l1: Vec<f64>,
    l0: Vec<f64>,
    lhet: Vec<f64>,
}

fn legacy_cells(d: &Legacy, p: &Params) -> Cells {
    let len = d.l * d.c;
    let (mut l1, mut l0, mut lhet) = (vec![0.0; len], vec![0.0; len], vec![0.0; len]);
    for i in 0..d.l {
        for c in 0..d.c {
            let k = i * d.c + c;
            if !d.mask(k) {
                continue;
            }
            let (a, n) = (d.alt[k], d.n[k]);
            let kk = p.k(d.kinds[i], c);
            let f1 = f_present(kk, p.purity[c], d.haploid[i]);
            l1[k] = present_logpmf(a, n, f1, p.rho1, p.lam);
            l0[k] = bb_logpmf(a, n, p.eps[c], p.rho0);
            lhet[k] = if d.haploid[i] { NEG_INF } else { bb_logpmf(a, n, f_germline_het(kk), p.rho1) };
        }
    }
    Cells { l1, l0, lhet }
}

/// tree_fit.log_noise for one locus
fn log_noise_locus(d: &Legacy, i: usize, a0: f64, b0: f64) -> f64 {
    let (mut sa, mut sn, mut lc) = (0.0, 0.0, 0.0);
    for c in 0..d.c {
        let k = i * d.c + c;
        if d.mask(k) {
            sa += d.alt[k];
            sn += d.n[k];
            lc += log_choose(d.n[k], d.alt[k]);
        }
    }
    lc + betaln(a0 + sa, b0 + sn - sa) - betaln(a0, b0)
}

/// tree_fit.root_scores for one locus: (log L_root, [het, hom, somatic] components)
fn root_locus(d: &Legacy, cells: &Cells, i: usize) -> (f64, [f64; 3]) {
    let r = i * d.c..(i + 1) * d.c;
    let het: f64 = cells.lhet[r.clone()].iter().sum();
    let hom = log_noise_locus(d, i, HOM_PRIOR.0, HOM_PRIOR.1);
    let soma: f64 = cells.l1[r].iter().sum();
    let comps = [het, hom, soma];
    let w = if d.haploid[i] { [NEG_INF, 0.5f64.ln(), 0.5f64.ln()] } else { [(1.0f64 / 3.0).ln(); 3] };
    (lse((0..3).map(|k| comps[k] + w[k])), comps)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Class {
    Uninformative,
    Noise,
    Germline,
    Shared,
    Private,
}

struct LegacyFit {
    cells: Cells,
    scores: Vec<LocusScore>,
    roots: Vec<(f64, [f64; 3])>,
}

fn legacy_fit(d: &Legacy, p: &Params, ctx: &Ctx, want_colony: bool) -> LegacyFit {
    let cells = legacy_cells(d, p);
    let mut scores = Vec::with_capacity(d.l);
    let mut roots = Vec::with_capacity(d.l);
    for i in 0..d.l {
        let r = i * d.c..(i + 1) * d.c;
        let root = root_locus(d, &cells, i);
        let ln = log_noise_locus(d, i, p.noise_a, p.noise_b);
        scores.push(score_locus(ctx, &cells.l1[r.clone()], &cells.l0[r], root.0, ln, want_colony));
        roots.push(root);
    }
    LegacyFit { cells, scores, roots }
}

/// tree_fit.classify (carrier_lr 2, wt_lr 1, noise_margin 1)
fn classify(d: &Legacy, fit: &LegacyFit, tree: &Tree) -> Vec<Class> {
    (0..d.l)
        .map(|i| {
            let (mut n_car, mut n_wt) = (0, 0);
            for c in 0..d.c {
                let k = i * d.c + c;
                let lr = (fit.cells.l1[k] - fit.cells.l0[k]) / LN10;
                n_car += (lr >= 2.0) as usize;
                n_wt += (-lr >= 1.0) as usize;
            }
            let sc = &fit.scores[i];
            let b = sc.best_tree;
            if n_car == 0 {
                Class::Uninformative
            } else if (sc.ln - sc.s_max) / LN10 >= 1.0 {
                Class::Noise
            } else if b == 0 && n_wt == 0 {
                Class::Germline
            } else if n_car >= 2 && n_wt >= 1 {
                Class::Shared
            } else if n_car == 1 && b != 0 && tree.nodes[b].is_tip() {
                Class::Private
            } else {
                Class::Uninformative
            }
        })
        .collect()
}

fn germline_het_candidates(d: &Legacy) -> Vec<bool> {
    let min_cov = 2usize.max(d.c / 2);
    (0..d.l)
        .map(|i| {
            let (mut ta, mut tn, mut ncov, mut ncov_alt, mut deep_zero) = (0.0, 0.0, 0usize, 0usize, false);
            for c in 0..d.c {
                let k = i * d.c + c;
                let (a, n) = (d.alt[k], d.n[k]);
                ta += a;
                tn += n;
                if n >= 3.0 {
                    ncov += 1;
                    ncov_alt += (a >= 1.0) as usize;
                }
                deep_zero |= n >= 8.0 && a == 0.0;
            }
            let v = if tn > 0.0 { ta / tn } else { f64::NAN };
            let frac_alt = if ncov > 0 { ncov_alt as f64 / ncov as f64 } else { 0.0 };
            d.autosomal[i] && tn >= 10.0 && v >= 0.2 && v <= 0.8 && ncov >= min_cov && frac_alt >= 1.0 && !deep_zero
        })
        .collect()
}

/// genotype_likelihood.estimate_balance: K_kind, b_c, rho1.
fn estimate_balance(d: &Legacy, p: &mut Params, cand: Option<Vec<bool>>) {
    let cand = match cand {
        Some(v) if v.iter().filter(|&&x| x).count() >= 5 => v,
        _ => germline_het_candidates(d),
    };
    let idx: Vec<usize> = (0..d.l).filter(|&i| cand[i]).collect();
    p.n_germline_het = idx.len();
    if idx.len() < 5 {
        p.k_kind.clear();
        p.k_global = 1.0;
        p.b = vec![1.0; d.c];
        return;
    }
    let cell = |i: usize, c: usize| {
        let k = i * d.c + c;
        (d.alt[k], d.n[k] - d.alt[k])
    };
    let (mut sa, mut sr) = (0.0, 0.0);
    for &i in &idx {
        for c in 0..d.c {
            let (a, r) = cell(i, c);
            sa += a;
            sr += r;
        }
    }
    p.k_global = (sa + 0.5) / (sr + 0.5);
    p.k_kind.clear();
    for kind in LOCUS_KINDS {
        let sel: Vec<usize> = idx.iter().cloned().filter(|&i| d.kinds[i] == kind).collect();
        if sel.len() >= 5 {
            let (mut a, mut r) = (0.0, 0.0);
            for &i in &sel {
                for c in 0..d.c {
                    let (x, y) = cell(i, c);
                    a += x;
                    r += y;
                }
            }
            p.k_kind.insert(kind, (a + 0.5) / (r + 0.5));
        }
    }
    let kc: Vec<f64> = idx.iter().map(|&i| p.k_kind.get(d.kinds[i]).copied().unwrap_or(p.k_global)).collect();
    for c in 0..d.c {
        let (mut a, mut kr) = (0.0, 0.0);
        for (q, &i) in idx.iter().enumerate() {
            let (x, y) = cell(i, c);
            a += x;
            kr += kc[q] * y;
        }
        p.b[c] = (a + 5.0) / (kr + 5.0);
    }
    let grid = [0.0, 0.003, 0.01, 0.02, 0.05, 0.1, 0.2];
    let mut best = (grid[0], NEG_INF);
    for &rho in &grid {
        let mut ll = 0.0;
        for (q, &i) in idx.iter().enumerate() {
            for c in 0..d.c {
                let k = i * d.c + c;
                if d.n[k] > 0.0 {
                    ll += bb_logpmf(d.alt[k], d.n[k], f_germline_het(kc[q] * p.b[c]), rho);
                }
            }
        }
        if ll > best.1 {
            best = (rho, ll);
        }
    }
    p.rho1 = best.0.max(1e-4);
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = v.len();
    if m % 2 == 1 { v[m / 2] } else { 0.5 * (v[m / 2 - 1] + v[m / 2]) }
}

/// tree_fit.update_params: purity per colony, background eps per colony, rho0.
fn update_params(d: &Legacy, p: &mut Params, fit: &LegacyFit, cls: &[Class], ctx: &Ctx) {
    let conf = |i: usize, eps_set: bool| {
        let s = &fit.scores[i];
        let ok_cls = if eps_set { matches!(cls[i], Class::Shared | Class::Private) } else { cls[i] == Class::Shared };
        s.best_tree != 0 && s.post_within >= 0.9 && s.bf >= 1.0 && ok_cls
    };
    let conf_eps: Vec<usize> = (0..d.l).filter(|&i| conf(i, true)).collect();
    if conf_eps.is_empty() {
        return;
    }
    let conf_pur: Vec<usize> = (0..d.l).filter(|&i| conf(i, false)).collect();
    // purity: ML on a grid over clade-member cells
    let grid: Vec<f64> = (0..96).map(|k| if k == 95 { 1.0 } else { 0.05 + k as f64 * (0.95 / 95.0) }).collect();
    let mut pur = vec![f64::NAN; d.c];
    for (c, pc) in pur.iter_mut().enumerate() {
        let cells: Vec<(f64, f64, f64, bool)> = conf_pur
            .iter()
            .filter(|&&i| ctx.in_clade[fit.scores[i].best_tree * d.c + c] && d.n[i * d.c + c] > 0.0)
            .map(|&i| (d.alt[i * d.c + c], d.n[i * d.c + c], p.k(d.kinds[i], c), d.haploid[i]))
            .collect();
        if cells.len() >= 3 {
            let mut best = (grid[0], NEG_INF);
            for &g in &grid {
                let ll: f64 = cells.iter().map(|&(a, n, k, h)| present_logpmf(a, n, f_present(k, g, h), p.rho1, p.lam)).sum();
                if ll > best.1 {
                    best = (g, ll);
                }
            }
            *pc = best.0;
        }
    }
    let mut fin: Vec<f64> = pur.iter().cloned().filter(|x| x.is_finite()).collect();
    if !fin.is_empty() {
        let med = median(&mut fin);
        p.purity = pur.iter().map(|&x| if x.is_finite() { x } else { med }).collect();
    }
    // eps: out-of-clade cells of the confidently placed loci, shrunk with 1000 pseudo-reads
    let (mut a0, mut n0) = (vec![0.0; d.c], vec![0.0; d.c]);
    let mut out_cells: Vec<(f64, f64, usize)> = Vec::new();
    for &i in &conf_eps {
        let b = fit.scores[i].best_tree;
        for c in 0..d.c {
            let k = i * d.c + c;
            if !ctx.in_clade[b * d.c + c] && d.n[k] > 0.0 {
                a0[c] += d.alt[k];
                n0[c] += d.n[k];
                out_cells.push((d.alt[k], d.n[k], c));
            }
        }
    }
    let g = (a0.iter().sum::<f64>() / n0.iter().sum::<f64>().max(1.0)).max(1e-4);
    p.eps = (0..d.c).map(|c| ((a0[c] + g * 1000.0) / (n0[c] + 1000.0)).max(1e-5)).collect();
    if !out_cells.is_empty() {
        let grid = [0.0, 0.001, 0.003, 0.01, 0.03, 0.1];
        let mut best = (grid[0], NEG_INF);
        for &rho in &grid {
            let ll: f64 = out_cells.iter().map(|&(a, n, c)| bb_logpmf(a, n, p.eps[c], rho)).sum();
            if ll > best.1 {
                best = (rho, ll);
            }
        }
        p.rho0 = best.0.max(1e-4);
    }
}

/// genotype_likelihood.infer_sex: Some('F' | 'M') or None (undecided).
fn infer_sex(loci: &[String], d: &Legacy) -> Option<char> {
    let (mut het, mut hom, mut y) = (0, 0, 0);
    for (i, name) in loci.iter().enumerate() {
        let contig = locus_contig(name);
        let (isx, isy) = (contig == "chrX" || contig == "X", contig == "chrY" || contig == "Y");
        if !isx && !isy {
            continue;
        }
        let (mut ta, mut tn, mut cov, mut cov_alt, mut alt2) = (0.0, 0.0, 0, 0, 0);
        for c in 0..d.c {
            let k = i * d.c + c;
            ta += d.alt[k];
            tn += d.n[k];
            if d.n[k] >= 3.0 {
                cov += 1;
                cov_alt += (d.alt[k] >= 1.0) as usize;
            }
            alt2 += (d.alt[k] >= 2.0) as usize;
        }
        let v = if tn > 0.0 { ta / tn } else { 0.0 };
        let carried = (if cov > 0 { cov_alt as f64 / cov as f64 } else { 0.0 }) >= 0.9;
        let gl = carried && tn >= 20.0;
        het += (gl && isx && (0.2..=0.8).contains(&v)) as usize;
        hom += (gl && isx && v >= 0.85) as usize;
        y += (isy && ta >= 4.0 && alt2 >= 2) as usize;
    }
    if het >= 2 {
        Some('F')
    } else if hom >= 3 || y >= 1 {
        Some('M')
    } else {
        None
    }
}

/// Full tree_fit.py estimation (3 rounds) then the final fit with per-colony posteriors.
fn run_legacy(loci: &[String], d: &mut Legacy, ctx: &Ctx) -> Vec<LocusScore> {
    let sex = infer_sex(loci, d);
    eprintln!("joint: sex (auto): {}", match sex { Some(s) => s.to_string(), None => "undecided -> F (diploid X)".into() });
    let male = sex == Some('M');
    for (i, name) in loci.iter().enumerate() {
        let contig = locus_contig(name);
        d.haploid[i] = contig == "chrY" || contig == "Y" || (male && (contig == "chrX" || contig == "X"));
    }
    let mut p = Params::new(d.c);
    estimate_balance(d, &mut p, None);
    eprintln!(
        "joint: germline-het candidates: {}; K_global={:.3} K_kind={:?}; rho1={}",
        p.n_germline_het, p.k_global, sorted_kinds(&p), p.rho1
    );
    for r in 0..3 {
        let fit = legacy_fit(d, &p, ctx, false);
        let cls = classify(d, &fit, ctx.tree);
        update_params(d, &mut p, &fit, &cls, ctx);
        let ghet: Vec<bool> = (0..d.l)
            .map(|i| cls[i] == Class::Germline && argmax(&fit.roots[i].1) == 0 && d.autosomal[i])
            .collect();
        estimate_balance(d, &mut p, Some(ghet));
        let mut pur = p.purity.clone();
        let mut eps = p.eps.clone();
        let (pmin, pmax) = (pur.iter().cloned().fold(f64::INFINITY, f64::min), pur.iter().cloned().fold(NEG_INF, f64::max));
        eprintln!(
            "joint: round {}: purity median {:.3} [{:.2}-{:.2}], eps median {:.5}, rho0={}",
            r + 1,
            median(&mut pur),
            pmin,
            pmax,
            median(&mut eps),
            p.rho0
        );
    }
    eprintln!(
        "joint: final K_global={:.3} K_kind={:?} from {} germline-het loci; rho1={} rho0={}",
        p.k_global, sorted_kinds(&p), p.n_germline_het, p.rho1, p.rho0
    );
    legacy_fit(d, &p, ctx, true).scores
}

fn sorted_kinds(p: &Params) -> Vec<(&'static str, f64)> {
    let mut v: Vec<(&'static str, f64)> = p.k_kind.iter().map(|(k, x)| (*k, (x * 1000.0).round() / 1000.0)).collect();
    v.sort_by(|a, b| a.0.cmp(b.0));
    v
}

// ------------------------------------------------------------------ output

/// Matrix vocabulary (SPEC "Joint"): `own` = the colony's own genotype string (None: no row).
fn matrix_call(own: Option<&str>, p_carrier: f64) -> String {
    let Some(own) = own else { return String::new() };
    if is_na_state(own) {
        return own.to_string();
    }
    if p_carrier >= 0.9 {
        if own == GT_HETEROZYGOUS || own == GT_HOMOZYGOUS { own.to_string() } else { GT_INSERTION.to_string() }
    } else if p_carrier >= 0.5 {
        GT_INSERTION_UNCERTAIN.to_string()
    } else if 1.0 - p_carrier >= 0.9 {
        GT_WILDTYPE.to_string()
    } else {
        GT_WILDTYPE_UNCERTAIN.to_string()
    }
}

/// pandas to_csv quoting (QUOTE_MINIMAL) for a `;`-separated file.
fn csv_field(s: &str) -> String {
    if s.contains([';', '"', '\n', '\r']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

fn fmt_f(x: f64, digits: usize) -> String {
    if x.is_nan() {
        "NA".into()
    } else if x == f64::INFINITY {
        "inf".into()
    } else if x == NEG_INF {
        "-inf".into()
    } else {
        let s = format!("{x:.digits$}");
        if s.starts_with('-') && s[1..].chars().all(|ch| ch == '0' || ch == '.') { s[1..].to_string() } else { s }
    }
}

fn write_out(path: &str, content: &[u8]) -> io::Result<()> {
    let f = File::create(path).map_err(|e| io::Error::new(e.kind(), format!("cannot create {path}: {e}")))?;
    if path.ends_with(".gz") {
        let mut enc = GzEncoder::new(io::BufWriter::new(f), Compression::default());
        enc.write_all(content)?;
        enc.finish()?.flush()?;
    } else {
        let mut w = io::BufWriter::new(f);
        w.write_all(content)?;
        w.flush()?;
    }
    Ok(())
}

fn best_id<'a>(ctx: &'a Ctx<'a>, b: Best) -> &'a str {
    match b {
        Best::Node(k) => &ctx.tree.nodes[k].id,
        Best::Indep => "INDEP",
        Best::Noise => "NOISE",
    }
}

fn clade_str(ctx: &Ctx, node: usize) -> String {
    ctx.tree.clade_labels(node).join(",")
}

/// Per-locus TSV and the legacy-vocabulary matrix as strings.
fn render(
    ctx: &Ctx,
    loci: &[String],
    scores: &[LocusScore],
    n_data: &[usize],
    own: &dyn Fn(usize, usize) -> Option<String>,
) -> (String, String) {
    let cols = &ctx.tree.tips;
    let mut tsv = String::from(
        "locus\tlocus_kind\tn_colonies_data\tbest\tcarriers\tn_carriers\tpost_best\tbest_tree\tbest_tree_clade\tpost_tree\tlog10_bf_tree\tlog10_L_best_tree\tlog10_L_root\tlog10_L_indep\tlog10_L_noise",
    );
    for c in cols {
        tsv.push_str("\tp_");
        tsv.push_str(c);
    }
    tsv.push('\n');
    let mut mat = String::new();
    for c in cols {
        mat.push(';');
        mat.push_str(&csv_field(c));
    }
    mat.push('\n');
    for (i, (name, s)) in loci.iter().zip(scores).enumerate() {
        let (carriers, n_car) = match s.best {
            Best::Node(b) => (clade_str(ctx, b), ctx.tree.nodes[b].clade.len()),
            _ => (String::new(), 0),
        };
        let fields = [
            name.clone(),
            locus_kind(name).to_string(),
            n_data[i].to_string(),
            best_id(ctx, s.best).to_string(),
            carriers,
            n_car.to_string(),
            fmt_f(s.post_best, 4),
            ctx.tree.nodes[s.best_tree].id.clone(),
            clade_str(ctx, s.best_tree),
            fmt_f(s.post_tree, 4),
            fmt_f(s.bf, 3),
            fmt_f(s.s_best_tree / LN10, 3),
            fmt_f(s.s_root / LN10, 3),
            fmt_f(s.li / LN10, 3),
            fmt_f(s.ln / LN10, 3),
        ];
        tsv.push_str(&fields.join("\t"));
        for p in &s.p_carrier {
            tsv.push('\t');
            tsv.push_str(&fmt_f(*p, 4));
        }
        tsv.push('\n');
        mat.push_str(&csv_field(name));
        for (c, p) in s.p_carrier.iter().enumerate() {
            mat.push(';');
            mat.push_str(&csv_field(&matrix_call(own(i, c).as_deref(), *p)));
        }
        mat.push('\n');
    }
    (tsv, mat)
}

// ------------------------------------------------------------------ driver

pub fn run(args: &JointArgs) -> io::Result<()> {
    if args.branch_prior != "length" && args.branch_prior != "uniform" {
        return Err(bad(format!("--branch-prior must be 'length' or 'uniform', got {:?}", args.branch_prior)));
    }
    if !(0.0..1.0).contains(&args.root_prior) {
        return Err(bad(format!("--root-prior must be in [0, 1), got {}", args.root_prior)));
    }
    let tree = Tree::read(&args.tree).map_err(bad)?;
    let tips = &tree.tips;

    // ---- match genotype files to tips
    let mut by_stem: HashMap<String, String> = HashMap::new();
    for f in &args.genotype_files {
        let stem = colony_stem(f);
        if let Some(prev) = by_stem.insert(stem.clone(), f.clone()) {
            return Err(bad(format!("joint: two genotype files for colony {stem}: {prev} and {f}")));
        }
    }
    let tipset: HashSet<&str> = tips.iter().map(|s| s.as_str()).collect();
    let missing_tips: Vec<&str> = tips.iter().map(|s| s.as_str()).filter(|t| !by_stem.contains_key(*t)).collect();
    let mut extra_files: Vec<&str> = args
        .genotype_files
        .iter()
        .filter(|f| !tipset.contains(colony_stem(f).as_str()))
        .map(|s| s.as_str())
        .collect();
    extra_files.sort();
    if !missing_tips.is_empty() || !extra_files.is_empty() {
        return Err(bad(format!(
            "joint: genotype files and tree tips do not match.\n  tips without a genotype file ({}): {}\n  \
             genotype files whose stem is not a tip ({}): {}",
            missing_tips.len(),
            missing_tips.join(", "),
            extra_files.len(),
            extra_files.join(", ")
        )));
    }

    // ---- read tables (in the given file order; colonies in tree tip order)
    let mut tables: Vec<Option<ColonyTable>> = (0..tips.len()).map(|_| None).collect();
    let mut loci: Vec<String> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let tip_ix: HashMap<&str, usize> = tips.iter().enumerate().map(|(i, t)| (t.as_str(), i)).collect();
    for f in &args.genotype_files {
        let t = read_colony_table(f)?;
        for n in &t.order {
            if seen.insert(n.clone()) {
                loci.push(n.clone());
            }
        }
        tables[tip_ix[colony_stem(f).as_str()]] = Some(t);
    }
    let tables: Vec<ColonyTable> = tables.into_iter().map(|t| t.unwrap()).collect();
    let pl_mode = tables.iter().all(|t| t.has_pl);
    if !pl_mode {
        if let Some(k) = tables.iter().position(|t| !t.has_counts) {
            return Err(bad(format!("joint: colony {} has neither PL nor n_alt/n_ref columns", tips[k])));
        }
        if tables.iter().any(|t| t.has_pl) {
            eprintln!("joint: WARNING: some files have PL columns, some do not -> legacy read-vote model for ALL colonies");
        }
    }
    let (l, c) = (loci.len(), tips.len());
    eprintln!(
        "joint: {c} colonies, {l} loci, {} branches; likelihoods from {}; branch prior {}, root prior {}",
        tree.nodes.len(),
        if pl_mode { "PL columns" } else { "legacy n_alt/n_ref read votes (beta-binomial)" },
        args.branch_prior,
        args.root_prior
    );
    let ctx = Ctx::new(&tree, &args.branch_prior, args.root_prior);

    let row = |i: usize, col: usize| tables[col].rows.get(&loci[i]);
    let mut n_data = vec![0usize; l];
    let scores: Vec<LocusScore> = if pl_mode {
        let mut out = Vec::with_capacity(l);
        let (mut l1, mut l0) = (vec![0.0; c], vec![0.0; c]);
        for (i, nd) in n_data.iter_mut().enumerate() {
            for col in 0..c {
                let cell = row(i, col).filter(|r| !is_na_state(&r.genotype)).and_then(|r| r.pl);
                (l1[col], l0[col]) = match cell {
                    Some(pl) => {
                        *nd += 1;
                        pl_cell(pl)
                    }
                    None => (0.0, 0.0),
                };
            }
            let root: f64 = l1.iter().sum();
            let noise: f64 = l0.iter().sum();
            out.push(score_locus(&ctx, &l1, &l0, root, noise, true));
        }
        out
    } else {
        let mut d = Legacy {
            l,
            c,
            kinds: loci.iter().map(|x| locus_kind(x)).collect(),
            haploid: vec![false; l],
            autosomal: loci.iter().map(|x| !matches!(locus_contig(x), "chrX" | "X" | "chrY" | "Y")).collect(),
            alt: vec![0.0; l * c],
            n: vec![0.0; l * c],
        };
        for (i, nd) in n_data.iter_mut().enumerate() {
            for col in 0..c {
                if let Some(r) = row(i, col).filter(|r| !is_na_state(&r.genotype)) {
                    let k = i * c + col;
                    d.alt[k] = r.n_alt as f64;
                    d.n[k] = (r.n_alt + r.n_ref) as f64;
                    *nd += (d.n[k] >= 1.0) as usize;
                }
            }
        }
        run_legacy(&loci, &mut d, &ctx)
    };

    let own = |i: usize, col: usize| row(i, col).map(|r| r.genotype.clone());
    let (tsv, mat) = render(&ctx, &loci, &scores, &n_data, &own);
    write_out(&args.out_tsv, tsv.as_bytes())?;
    write_out(&args.out_matrix, mat.as_bytes())?;
    let n_tree = scores.iter().filter(|s| matches!(s.best, Best::Node(_))).count();
    let n_noise = scores.iter().filter(|s| s.best == Best::Noise).count();
    eprintln!(
        "joint: best hypothesis: {n_tree} tree ({} ROOT), {} INDEP, {n_noise} NOISE; wrote {} and {}",
        scores.iter().filter(|s| s.best == Best::Node(0)).count(),
        l - n_tree - n_noise,
        args.out_tsv,
        args.out_matrix
    );
    Ok(())
}

// ------------------------------------------------------------------ tests

#[cfg(test)]
mod tests {
    use super::*;

    const TOY: &str = "((A:1,B:1):2,(C:1,D:1):2);";

    fn ctx_for(t: &Tree) -> Ctx<'_> {
        Ctx::new(t, "length", 0.1)
    }

    fn score_pl(ctx: &Ctx, pls: &[Option<[f64; 3]>]) -> LocusScore {
        let (mut l1, mut l0) = (vec![0.0; pls.len()], vec![0.0; pls.len()]);
        for (k, p) in pls.iter().enumerate() {
            if let Some(p) = p {
                (l1[k], l0[k]) = pl_cell(*p);
            }
        }
        score_locus(ctx, &l1, &l0, l1.iter().sum(), l0.iter().sum(), true)
    }

    const HET: Option<[f64; 3]> = Some([200.0, 0.0, 150.0]);
    const WT: Option<[f64; 3]> = Some([0.0, 60.0, 300.0]);

    #[test]
    fn special_functions() {
        assert!((ln_gamma(5.0) - 24f64.ln()).abs() < 1e-12);
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-12);
        assert!((ln_gamma(1e-6) - 13.815509980749431).abs() < 1e-9);
        assert!((betaln(3.0, 4.0) - (1.0f64 / 60.0).ln()).abs() < 1e-12);
        // I_0.3(2, 5) = scipy.special.betainc(2, 5, 0.3)
        assert!((betainc(2.0, 5.0, 0.3) - 0.579825).abs() < 1e-6);
        assert!((betainc(5.0, 2.0, 0.9) - 0.885735).abs() < 1e-6);
        // pmfs sum to one
        for &(f, rho) in &[(0.3, 0.0), (0.3, 0.05), (1e-6, 0.1), (0.999, 0.02)] {
            let s: f64 = (0..=20).map(|a| bb_logpmf(a as f64, 20.0, f, rho).exp()).sum();
            assert!((s - 1.0).abs() < 1e-9, "bb f={f} rho={rho} sum={s}");
        }
        let s: f64 = (0..=15).map(|a| present_logpmf(a as f64, 15.0, 0.4, 0.02, 0.01).exp()).sum();
        assert!((s - 1.0).abs() < 1e-9);
        assert!((lae(NEG_INF, -1.0) + 1.0).abs() < 1e-15);
        assert_eq!(locus_kind("chr1:100-112"), "TSD");
        assert_eq!(locus_kind("chr1:100-100"), "BLUNT");
        assert_eq!(locus_kind("chr1:100-90"), "TSD_DELETION");
        assert_eq!(locus_kind("chr1:100-oneside_100"), "ONE_SIDED");
        assert_eq!(locus_kind("chr1:1000-500"), "L1_MED_DELETION");
        assert_eq!(locus_kind("chr1:100-200"), "L1_MED_DUPLICATION");
    }

    #[test]
    fn pl_to_likelihood() {
        let (l1, l0) = pl_cell([0.0, 30.0, 60.0]);
        assert!(l0.abs() < 1e-15);
        assert!((l1 - (0.5 * (1e-3 + 1e-6f64)).ln()).abs() < 1e-12);
        // normalised to the max even when the min PL is not 0
        let (l1, l0) = pl_cell([20.0, 10.0, 30.0]);
        assert!((l0 - 0.1f64.ln()).abs() < 1e-12);
        assert!((l1 - (0.5 * (1.0 + 0.01f64)).ln()).abs() < 1e-12);
    }

    #[test]
    fn obvious_clade() {
        let t = Tree::parse(TOY).unwrap();
        let ctx = ctx_for(&t);
        let s = score_pl(&ctx, &[HET, HET, WT, WT]);
        assert_eq!(s.best, Best::Node(1));
        assert_eq!(t.nodes[1].id, "N1");
        assert_eq!(clade_str(&ctx, 1), "A,B");
        // 4 tips: INDEP (∫π²(1-π)² = 1/30, weight 1/4) keeps ~7 % of the mass on any pattern
        assert!(s.post_best > 0.9, "{}", s.post_best);
        assert!(s.bf > 1.0, "{}", s.bf);
        assert!(s.p_carrier[0] > 0.999 && s.p_carrier[1] > 0.999, "{:?}", s.p_carrier);
        assert!(s.p_carrier[2] < 1e-3 && s.p_carrier[3] < 1e-3, "{:?}", s.p_carrier);
        // all carriers -> ROOT; nobody -> NOISE
        assert_eq!(score_pl(&ctx, &[HET, HET, HET, HET]).best, Best::Node(0));
        assert_eq!(score_pl(&ctx, &[WT, WT, WT, WT]).best, Best::Noise);
    }

    #[test]
    fn dropout_member_still_assigned() {
        let t = Tree::parse("((A:1,B:1,E:1):2,((C:1,D:1):1,(F:1,G:1,H:1):1):2);").unwrap();
        let ctx = ctx_for(&t);
        assert_eq!(clade_str(&ctx, 1), "A,B,E");
        // B: no coverage
        let s = score_pl(&ctx, &[HET, None, HET, WT, WT, WT, WT, WT]);
        assert_eq!(s.best, Best::Node(1));
        assert!(s.p_carrier[1] > 0.95, "{:?}", s.p_carrier);
        assert_eq!(matrix_call(Some(GT_NO_COVERAGE), s.p_carrier[1]), GT_NO_COVERAGE);
        // B: low depth, weakly favouring the reference (one ref read): imputed carrier
        let s = score_pl(&ctx, &[HET, Some([0.0, 3.0, 6.0]), HET, WT, WT, WT, WT, WT]);
        assert_eq!(s.best, Best::Node(1));
        assert!(s.p_carrier[1] > 0.95, "{:?}", s.p_carrier);
        assert!(s.p_carrier[3] < 1e-3, "{:?}", s.p_carrier);
        assert_eq!(matrix_call(Some(GT_WILDTYPE_UNCERTAIN), s.p_carrier[1]), GT_INSERTION);
    }

    #[test]
    fn scattered_pattern_is_not_a_tree_event() {
        let t = Tree::parse(TOY).unwrap();
        let ctx = ctx_for(&t);
        let s = score_pl(&ctx, &[HET, WT, HET, WT]);
        assert!(matches!(s.best, Best::Indep | Best::Noise), "{:?}", s.best);
        assert!(s.bf < -5.0, "bf {}", s.bf);
        assert_eq!(s.best, Best::Indep);
        assert!(s.p_carrier[0] > 0.99 && s.p_carrier[2] > 0.99);
        assert!(s.p_carrier[1] < 0.01 && s.p_carrier[3] < 0.01);
    }

    /// brute force over all 2^C presence vectors: P(z | INDEP) = B(k+1, C-k+1)
    fn brute(l1: &[f64], l0: &[f64]) -> (f64, Vec<f64>) {
        let c = l1.len();
        let mut tot = NEG_INF;
        let mut pres = vec![NEG_INF; c];
        for z in 0u32..(1 << c) {
            let k = z.count_ones() as f64;
            let mut ll = betaln(k + 1.0, c as f64 - k + 1.0);
            for j in 0..c {
                ll += if z >> j & 1 == 1 { l1[j] } else { l0[j] };
            }
            tot = lae(tot, ll);
            for j in 0..c {
                if z >> j & 1 == 1 {
                    pres[j] = lae(pres[j], ll);
                }
            }
        }
        (tot, pres)
    }

    #[test]
    fn indep_dp_equals_brute_force() {
        let cases: [([f64; 4], [f64; 4]); 3] = [
            ([-0.3, -5.0, -1.2, 0.0], [-2.0, -0.1, -0.7, 0.0]),
            ([-40.0, -0.01, -3.0, -9.0], [-0.02, -30.0, -0.5, -0.4]),
            ([0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]),
        ];
        let t = Tree::parse(TOY).unwrap();
        let ctx = ctx_for(&t);
        for (l1, l0) in cases {
            let (li, lp) = indep(&l1, &l0, &ctx.bl_full, &ctx.bl_v, true);
            let (bt, bp) = brute(&l1, &l0);
            assert!((li - bt).abs() < 1e-10, "{li} vs {bt}");
            for j in 0..4 {
                assert!((lp[j] - bp[j]).abs() < 1e-10, "colony {j}: {} vs {}", lp[j], bp[j]);
            }
        }
        // all missing: P(present) = E[pi] = 1/2
        let (li, lp) = indep(&[0.0; 4], &[0.0; 4], &ctx.bl_full, &ctx.bl_v, true);
        assert!(li.abs() < 1e-12);
        assert!((lp[2] - 0.5f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn branch_prior_length_and_uniform() {
        let t = Tree::parse("((A:0,B:3):1,C:4);").unwrap();
        let lp = branch_log_prior(&t, "length", 0.1);
        let p: Vec<f64> = lp.iter().map(|x| x.exp()).collect();
        // nodes: ROOT, N1(1), A(0 -> floor 0.01*8/3), B(3), C(4)
        let floor = 0.01 * 8.0 / 3.0;
        let tot = 1.0 + floor + 3.0 + 4.0;
        assert!((p[0] - 0.1).abs() < 1e-12);
        assert!((p[2] - 0.9 * floor / tot).abs() < 1e-12);
        assert!((p[4] - 0.9 * 4.0 / tot).abs() < 1e-12);
        let pu: Vec<f64> = branch_log_prior(&t, "uniform", 0.1).iter().map(|x| x.exp()).collect();
        assert!((pu[3] - 0.9 / 4.0).abs() < 1e-12);
    }

    #[test]
    fn matrix_vocabulary() {
        assert_eq!(matrix_call(Some("heterozygous"), 0.95), "heterozygous");
        assert_eq!(matrix_call(Some("homozygous"), 0.9), "homozygous");
        assert_eq!(matrix_call(Some("wild-type?"), 0.97), "insertion");
        assert_eq!(matrix_call(Some("artefact"), 0.99), "insertion");
        assert_eq!(matrix_call(Some("heterozygous"), 0.7), "insertion?");
        assert_eq!(matrix_call(Some("heterozygous"), 0.5), "insertion?");
        assert_eq!(matrix_call(Some("heterozygous"), 0.05), "wild-type");
        assert_eq!(matrix_call(Some("insertion"), 0.3), "wild-type?");
        assert_eq!(matrix_call(Some("high-coverage"), 0.99), "high-coverage");
        assert_eq!(matrix_call(Some("error"), 0.01), "error");
        assert_eq!(matrix_call(Some("no-coverage"), 0.5), "no-coverage");
        assert_eq!(matrix_call(None, 0.99), "");
    }

    #[test]
    fn csv_and_number_formatting() {
        assert_eq!(csv_field("chr1:100-112"), "chr1:100-112");
        assert_eq!(csv_field("a;b"), "\"a;b\"");
        assert_eq!(csv_field("say \"x\""), "\"say \"\"x\"\"\"");
        assert_eq!(fmt_f(-0.00001, 3), "0.000");
        assert_eq!(fmt_f(1.23456, 3), "1.235");
        assert_eq!(fmt_f(NEG_INF, 3), "-inf");
        assert_eq!(colony_stem("/x/y/S1.txt.gz"), "S1");
        assert_eq!(colony_stem("PD1234b.genotypes.tsv"), "PD1234b");
    }

    fn write_gz(path: &std::path::Path, text: &str) {
        let f = File::create(path).unwrap();
        let mut e = GzEncoder::new(f, Compression::default());
        e.write_all(text.as_bytes()).unwrap();
        e.finish().unwrap();
    }

    fn read_gz(path: &std::path::Path) -> String {
        let mut s = String::new();
        MultiGzDecoder::new(File::open(path).unwrap()).read_to_string(&mut s).unwrap();
        s
    }

    #[test]
    fn end_to_end_pl_files() {
        let dir = std::env::temp_dir().join(format!("peartree_joint_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("t.nwk"), TOY).unwrap();
        let head = "insertion\tgenotype\tscore_genotype\tscore_alternative\tcoverage\tn_alt\tn_ref\tn_art\tvaf\tgq\tpl_ref\tpl_het\tpl_hom\tn_uninf\tn_disc\n";
        let het = |l: &str| format!("{l}\theterozygous\t1\t1\t20\t10\t10\t0\t0.5\t99\t200\t0\t150\t0\t0\n");
        let wt = |l: &str| format!("{l}\twild-type\t1\t1\t20\t0\t20\t0\t0\t60\t0\t60\t300\t0\t0\n");
        let noc = |l: &str| format!("{l}\tno-coverage\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
        let (l1, l2) = ("chr1:100-112", "chr1:500-515");
        let files = [
            ("A", format!("{head}{}{}", het(l1), het(l2))),
            ("B", format!("{head}{}{}", noc(l1), wt(l2))),
            ("C", format!("{head}{}{}", wt(l1), het(l2))),
            ("D", format!("{head}{}", wt(l1))), // l2 row absent
        ];
        let mut paths = Vec::new();
        for (n, txt) in &files {
            let p = dir.join(format!("{n}.txt.gz"));
            write_gz(&p, txt);
            paths.push(p.to_string_lossy().into_owned());
        }
        let out_tsv = dir.join("P.joint.tsv");
        let out_mat = dir.join("P.joint_matrix.csv.gz");
        let args = JointArgs {
            tree: dir.join("t.nwk").to_string_lossy().into_owned(),
            genotype_files: paths.clone(),
            out_tsv: out_tsv.to_string_lossy().into_owned(),
            out_matrix: out_mat.to_string_lossy().into_owned(),
            root_prior: 0.1,
            branch_prior: "length".into(),
        };
        run(&args).unwrap();
        let tsv = std::fs::read_to_string(&out_tsv).unwrap();
        let lines: Vec<&str> = tsv.lines().collect();
        assert!(lines[0].starts_with("locus\tlocus_kind\tn_colonies_data\tbest\tcarriers"));
        assert!(lines[0].ends_with("\tp_A\tp_B\tp_C\tp_D"));
        let f1: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(&f1[..6], &[l1, "TSD", "3", "N1", "A,B", "2"]);
        let f2: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(f2[3], "INDEP"); // A and C carry it, B does not: scattered
        let mat = read_gz(&out_mat);
        let m: Vec<&str> = mat.lines().collect();
        assert_eq!(m[0], ";A;B;C;D");
        assert_eq!(m[1], format!("{l1};heterozygous;no-coverage;wild-type;wild-type"));
        assert!(m[2].starts_with(&format!("{l2};heterozygous;wild-type;heterozygous;")));
        assert!(m[2].ends_with(';')); // row absent in D -> empty cell
        // unmatched tips / files are an error
        let mut bad_args = JointArgs { genotype_files: paths[..3].to_vec(), ..args };
        assert!(run(&bad_args).unwrap_err().to_string().contains("tips without a genotype file (1): D"));
        let extra = dir.join("S11.txt.gz");
        write_gz(&extra, &files[0].1);
        bad_args.genotype_files = paths.clone();
        bad_args.genotype_files.push(extra.to_string_lossy().into_owned());
        assert!(run(&bad_args).unwrap_err().to_string().contains("S11.txt.gz"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn legacy_counts_obvious_clade() {
        // read-vote model on n_alt/n_ref (no estimation: default Params)
        let t = Tree::parse(TOY).unwrap();
        let ctx = ctx_for(&t);
        let d = Legacy {
            l: 1,
            c: 4,
            kinds: vec!["TSD"],
            haploid: vec![false],
            autosomal: vec![true],
            alt: vec![9.0, 0.0, 0.0, 0.0],
            n: vec![20.0, 0.0, 18.0, 25.0],
        };
        let fit = legacy_fit(&d, &Params::new(4), &ctx, true);
        let s = &fit.scores[0];
        // B has no reads: {A,B} (length 2) beats tip A (length 1) on the prior
        assert_eq!(s.best, Best::Node(1));
        assert!(s.p_carrier[0] > 0.999);
        assert!(s.p_carrier[2] < 1e-3);
    }
}
