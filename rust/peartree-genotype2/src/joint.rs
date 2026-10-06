//! Joint phylogenetic genotyping across the colonies of one patient (SPEC "Joint phylogenetic
//! step"). Owner: E. Owns newick.rs too.
//!
//! Port of the tools/phylo/tree_fit.py hypotheses onto the NUMERIC per-colony files of
//! `--step genotype` (header `types::OUTPUT_HEADER`; columns are found by name). No legacy reader.
//!
//! Per locus i and colony c: `l0 = log P(d_c | absent) = ln 10^(-pl_absent/10)` and
//! `l1 = log P(d_c | present) = ln ½(10^(-pl_het/10) + 10^(-pl_hom/10))` (PLs normalised to the
//! best genotype). A row with `status != ok`, a row whose PLs are not numbers, or a missing row
//! contributes log 1 to every hypothesis.
//!
//! Hypotheses: every node b of the tree (ROOT = every colony present, and each branch: the tips
//! below it present, all others absent), with prior `root_prior` for ROOT and `1 - root_prior`
//! spread over branches by length (floored at 1 % of the mean positive length) or uniformly;
//! INDEP (independent presence, π ~ U(0,1), exact DP over the carrier count) and NOISE.
//! NOISE is approximated as "absent everywhere" (Σ l0): the PLs already integrate the
//! background alt rate, so the shared per-locus alt fraction of tree_fit's NOISE cannot be
//! refitted from them.
//! Genotype-error terms (PD37590, 2026-10): a germline locus present in 41 of 44 colonies has 1-3
//! colonies with 0-1 alt reads at depth 10-40 whose alt reads realign as uninformative /
//! unexplained -- a hard "absent" (PL 20-50) that costs ROOT 2-5 log10 units each and hands the
//! locus to INDEP. So per colony `P(d_c | present) = (1-ε₁) P1 + ε₁ P0` (`dropout`, default 0.02)
//! and `P(d_c | absent) = (1-ε₀) P0 + ε₀ P1` (`false_present`, default 0: a single strongly
//! present colony at an otherwise absent locus is a private event, not an error, unless asked).
//! NOISE with a shared background: when the per-colony files carry the profile columns
//! `pl_f<‰>` (one alt fraction φ shared by every colony, `cfg.noise_frac_grid`), NOISE is the
//! mean over {absent everywhere, each φ of the grid} of Π_c P(d_c | φ): mismapped paralogous
//! reads or slippage give every colony the same low alt fraction, which three-genotype PLs
//! cannot express (PD37590: 65 of 88 "clade" calls were such loci for tree_fit's noise class).
//! Files without the columns fall back to NOISE = absent everywhere.
//! `log10_bf_tree = log10(Σ_b prior_b L_b / mean(L_indep, L_noise))`. The posterior over ALL
//! hypotheses uses model weights ½ (tree) / ¼ (INDEP) / ¼ (NOISE), so the posterior odds of
//! "a tree event" vs "not" equal the BF. Per colony
//! `P(carrier) = Σ_{b ∋ c} post_b + post_INDEP · P(c present | INDEP, data)` (exact, O(C²)).
//!
//! Tree tips are matched to file stems. A file whose stem is not a tip is an error; a tip
//! without a file (farm: a colony without a BAM) is a warning: that colony is all-missing data
//! (log 1 under every hypothesis) and keeps its (empty) `p_<colony>` / matrix columns.
//!
//! Outputs: a per-locus TSV (`--out`) and a numeric `;`-separated matrix (`--matrix`, rows =
//! loci, columns = colonies in tree-tip order, cell = P(carrier), empty when the colony has no
//! usable row). Either is gzipped when its path ends in `.gz`.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::newick::Tree;
use crate::types::Status;

pub struct JointArgs {
    pub tree: String,
    /// per-colony genotype files (stem = colony id = tree tip label)
    pub genotype_files: Vec<String>,
    pub out_tsv: String,
    pub out_matrix: String,
    pub root_prior: f64,
    /// "length" | "uniform"
    pub branch_prior: String,
    /// ε₁: P(a carrier colony looks absent) -- allelic dropout / alt reads the realignment cannot place
    pub dropout: f64,
    /// ε₀: P(a non-carrier colony looks present) -- contamination / mismapping; 0 = off
    pub false_present: f64,
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

// ------------------------------------------------------------------ small numerics

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

/// ln k! for k = 0..=n
fn ln_factorials(n: usize) -> Vec<f64> {
    let mut v = Vec::with_capacity(n + 1);
    let mut acc = 0.0;
    v.push(0.0);
    for k in 1..=n {
        acc += (k as f64).ln();
        v.push(acc);
    }
    v
}

/// ln B(a, b) for positive integers: (a-1)! (b-1)! / (a+b-1)!
fn betaln_int(a: usize, b: usize, lf: &[f64]) -> f64 {
    lf[a - 1] + lf[b - 1] - lf[a + b - 1]
}

// ------------------------------------------------------------------ input

#[derive(Clone, Debug, PartialEq)]
struct Row {
    kind: String,
    /// [pl_absent, pl_het, pl_hom] of a `status == ok` row with numeric PLs; None = no data
    pl: Option<[f64; 3]>,
    /// `pl_f<‰>` profile (same Phred scale as `pl`), one per `ColonyTable::fracs`; empty when
    /// the file has no profile columns or the row has no data
    pl_frac: Vec<f64>,
}

#[derive(Debug)]
struct ColonyTable {
    rows: HashMap<String, Row>,
    order: Vec<String>,
    /// `status == ok` rows whose PL columns were not numbers (treated as missing)
    n_bad_pl: usize,
    /// shared alt fractions of the `pl_f<‰>` columns, in column order (empty: no profile)
    fracs: Vec<f64>,
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

/// Read one numeric per-colony file (header-driven: `locus`, `status`, `pl_absent`, `pl_het`,
/// `pl_hom` required, `kind` optional; other columns ignored). Duplicate loci: first row wins.
fn read_colony_table_from(path: &str, input: Box<dyn BufRead>) -> io::Result<ColonyTable> {
    let mut lines = input.lines();
    let header = lines.next().transpose()?.ok_or_else(|| bad(format!("{path}: empty genotype file")))?;
    let cols: Vec<&str> = header.trim_end_matches(['\r', '\n']).split('\t').collect();
    let find = |name: &str| cols.iter().position(|c| *c == name);
    let required = ["locus", "status", "pl_absent", "pl_het", "pl_hom"];
    let missing: Vec<&str> = required.iter().cloned().filter(|c| find(c).is_none()).collect();
    if !missing.is_empty() {
        return Err(bad(format!(
            "{path}: not a peartree-genotype2 per-colony file (missing column(s) {}); legacy genotype files are not supported",
            missing.join(", ")
        )));
    }
    let (i_locus, i_status) = (find("locus").unwrap(), find("status").unwrap());
    let i_pl = [find("pl_absent").unwrap(), find("pl_het").unwrap(), find("pl_hom").unwrap()];
    let i_kind = find("kind");
    // profile columns pl_f<per-mille>, in header order
    let mut i_frac: Vec<usize> = Vec::new();
    let mut fracs: Vec<f64> = Vec::new();
    for (k, c) in cols.iter().enumerate() {
        if let Some(d) = c.strip_prefix("pl_f") {
            if let Ok(pm) = d.parse::<u32>() {
                i_frac.push(k);
                fracs.push(pm as f64 / 1000.0);
            }
        }
    }
    let ok = Status::Ok.as_str();
    let mut t = ColonyTable { rows: HashMap::new(), order: Vec::new(), n_bad_pl: 0, fracs };
    for line in lines {
        let line = line?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let Some(name) = f.get(i_locus).filter(|s| !s.is_empty()) else { continue };
        if t.rows.contains_key(*name) {
            continue;
        }
        let pl = if f.get(i_status).map(|s| s.trim()) == Some(ok) {
            let v: Vec<f64> = i_pl.iter().filter_map(|&k| f.get(k).and_then(|x| x.trim().parse::<f64>().ok())).collect();
            let parsed = (v.len() == 3 && v.iter().all(|x| x.is_finite())).then(|| [v[0], v[1], v[2]]);
            t.n_bad_pl += parsed.is_none() as usize;
            parsed
        } else {
            None
        };
        let kind = i_kind.and_then(|k| f.get(k)).map(|s| s.trim().to_string()).unwrap_or_default();
        let pl_frac: Vec<f64> = if pl.is_some() {
            let v: Vec<f64> = i_frac.iter().filter_map(|&k| f.get(k).and_then(|x| x.trim().parse::<f64>().ok())).collect();
            if v.len() == i_frac.len() && v.iter().all(|x| x.is_finite()) { v } else { Vec::new() }
        } else {
            Vec::new()
        };
        t.order.push(name.to_string());
        t.rows.insert(name.to_string(), Row { kind, pl, pl_frac });
    }
    Ok(t)
}

fn read_colony_table(path: &str) -> io::Result<ColonyTable> {
    read_colony_table_from(path, open_text(path)?)
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
    /// ln B(k+1, C-k+1), k = 0..=C (INDEP integral of π^k (1-π)^(C-k))
    bl_full: Vec<f64>,
    /// ln B(i+2, C-i), i = 0..C (INDEP per-colony presence functional, see `indep`)
    bl_v: Vec<f64>,
    /// node indices on the path tip -> root, per colony (the hypotheses that contain it)
    paths: Vec<Vec<usize>>,
}

impl<'a> Ctx<'a> {
    fn new(tree: &'a Tree, branch_prior: &str, root_prior: f64) -> Ctx<'a> {
        let c = tree.tips.len();
        let lf = ln_factorials(c + 1);
        let bl_full = (0..=c).map(|k| betaln_int(k + 1, c - k + 1, &lf)).collect();
        let bl_v = (0..c).map(|i| betaln_int(i + 2, c - i, &lf)).collect();
        let paths = (0..c)
            .map(|t| {
                let mut p = vec![tree.tip_node[t]];
                while let Some(up) = tree.nodes[*p.last().unwrap()].parent {
                    p.push(up);
                }
                p
            })
            .collect();
        Ctx { tree, ncol: c, logprior: branch_log_prior(tree, branch_prior, root_prior), bl_full, bl_v, paths }
    }
}

/// INDEP: log ∫_0^1 Π_c [π P1_c + (1-π) P0_c] dπ, exactly, by DP over the carrier count
/// (tree_fit.log_indep), and log P(z_c = 1, data | INDEP) per colony.
///
/// Per colony: P(z_c=1, d) = P1_c Σ_i pre_c[i] V_c(i), where pre_c = Bernstein-like coefficients
/// of the colonies before c and V_c(i) = Σ_j suf_c[j] B(i+j+2, C-i-j) folds in the colonies after
/// c. V obeys V_{c-1}(i) = P1_c V_c(i+1) + P0_c V_c(i) with V_{C-1}(i) = B(i+2, C-i), so all
/// colonies cost O(C²) together.
fn indep(l1: &[f64], l0: &[f64], bl_full: &[f64], bl_v: &[f64]) -> (f64, Vec<f64>) {
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
    if c == 0 {
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
    /// log L of the ROOT hypothesis
    s_root: f64,
    /// MAP tree hypothesis (with prior) and its log L
    best_tree: usize,
    s_best_tree: f64,
    li: f64,
    ln: f64,
    /// log10 BF tree vs mean(INDEP, NOISE)
    bf: f64,
    best: Best,
    post_best: f64,
    post_tree: f64,
    p_carrier: Vec<f64>,
}

/// Score one locus from per-colony `l1`/`l0` (missing cells 0 in both) and the shared-fraction
/// profile `lf[k][c]` = log P(d_c | φ_k) (missing cells 0; `lf` empty = no profile).
/// ROOT = Σ l1 (every colony present); NOISE = mean over {Σ l0, Σ_c lf[k][c] for each k}.
fn score_locus(ctx: &Ctx, l1: &[f64], l0: &[f64], lf: &[Vec<f64>]) -> LocusScore {
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
    s[0] = l1.iter().sum();
    for b in 1..nb {
        s[b] = sum_d[b] + sum0;
    }
    let j: Vec<f64> = s.iter().zip(&ctx.logprior).map(|(a, b)| a + b).collect();
    let log_tree = lse(j.iter().cloned());
    let best_tree = argmax(&j);
    let (li, lpres) = indep(l1, l0, &ctx.bl_full, &ctx.bl_v);
    let ln = if lf.is_empty() {
        sum0
    } else {
        let mut terms: Vec<f64> = vec![sum0];
        terms.extend(lf.iter().map(|row| row.iter().sum::<f64>()));
        lse(terms.iter().cloned()) - (terms.len() as f64).ln()
    };
    let bf = (log_tree - (lae(li, ln) - LN2)) / LN10;
    let z = lse([LN_W_TREE + log_tree, LN_W_ALT + li, LN_W_ALT + ln]);
    let cand = [(Best::Node(best_tree), LN_W_TREE + j[best_tree]), (Best::Indep, LN_W_ALT + li), (Best::Noise, LN_W_ALT + ln)];
    let mut pick = cand[0];
    for x in &cand[1..] {
        if x.1 > pick.1 {
            pick = *x;
        }
    }
    let p_carrier = (0..ctx.ncol)
        .map(|c| {
            let tree_part: f64 = ctx.paths[c].iter().map(|&b| (LN_W_TREE + j[b] - z).exp()).sum();
            (tree_part + (LN_W_ALT + lpres[c] - z).exp()).min(1.0)
        })
        .collect();
    LocusScore {
        s_root: s[0],
        best_tree,
        s_best_tree: s[best_tree],
        li,
        ln,
        bf,
        best: pick.0,
        post_best: (pick.1 - z).exp(),
        post_tree: (LN_W_TREE + log_tree - z).exp(),
        p_carrier,
    }
}

/// (log P(d|present), log P(d|absent)) from Phred-scaled genotype likelihoods
/// [pl_absent, pl_het, pl_hom], normalised to the best genotype.
fn pl_cell(pl: [f64; 3]) -> (f64, f64) {
    let s = -LN10 / 10.0;
    let (r, h, o) = (s * pl[0], s * pl[1], s * pl[2]);
    let m = r.max(h).max(o);
    let (r, h, o) = (r - m, h - m, o - m);
    (lae(h, o) - LN2, r)
}

/// Mix the genotype-error terms into one colony's (log P(d|present), log P(d|absent)):
/// present = (1-ε₁) P1 + ε₁ P0, absent = (1-ε₀) P0 + ε₀ P1. ε = 0 leaves the term unchanged.
fn with_errors(l1: f64, l0: f64, dropout: f64, false_present: f64) -> (f64, f64) {
    let mix = |keep: f64, other: f64, eps: f64| {
        if eps <= 0.0 {
            keep
        } else {
            lae((1.0 - eps).ln() + keep, eps.ln() + other)
        }
    };
    (mix(l1, l0, dropout), mix(l0, l1, false_present))
}

// ------------------------------------------------------------------ output

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

/// Per-locus TSV and the numeric P(carrier) matrix as strings. `usable[i * C + c]` = colony c has
/// a `status == ok` row with numeric PLs at locus i (else its matrix cell is empty);
/// `has_file[c]` = colony c has a genotype file (else its TSV `p_<colony>` column is empty too).
fn render(
    ctx: &Ctx,
    loci: &[String],
    kinds: &[String],
    scores: &[LocusScore],
    usable: &[bool],
    has_file: &[bool],
) -> (String, String) {
    let cols = &ctx.tree.tips;
    let nc = cols.len();
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
        let row_usable = &usable[i * nc..(i + 1) * nc];
        let fields = [
            name.clone(),
            if kinds[i].is_empty() { ".".to_string() } else { kinds[i].clone() },
            row_usable.iter().filter(|&&u| u).count().to_string(),
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
        for (p, &f) in s.p_carrier.iter().zip(has_file) {
            tsv.push('\t');
            if f {
                tsv.push_str(&fmt_f(*p, 4));
            }
        }
        tsv.push('\n');
        mat.push_str(&csv_field(name));
        for (p, &u) in s.p_carrier.iter().zip(row_usable) {
            mat.push(';');
            if u {
                mat.push_str(&fmt_f(*p, 4));
            }
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
    for (name, v) in [("--dropout", args.dropout), ("--false-present", args.false_present)] {
        if !(0.0..0.5).contains(&v) {
            return Err(bad(format!("{name} must be in [0, 0.5), got {v}")));
        }
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
    if !extra_files.is_empty() {
        return Err(bad(format!(
            "joint: genotype files whose stem is not a tree tip ({}): {}",
            extra_files.len(),
            extra_files.join(", ")
        )));
    }
    if !missing_tips.is_empty() {
        eprintln!(
            "joint: WARNING: {} tree tip(s) without a genotype file (all-missing data, empty columns): {}",
            missing_tips.len(),
            missing_tips.join(", ")
        );
    }

    // ---- read tables (loci in order of first appearance over the given files; colonies in tip order)
    let mut tables: Vec<Option<ColonyTable>> = (0..tips.len()).map(|_| None).collect();
    let mut loci: Vec<String> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let tip_ix: HashMap<&str, usize> = tips.iter().enumerate().map(|(i, t)| (t.as_str(), i)).collect();
    for f in &args.genotype_files {
        let t = read_colony_table(f)?;
        if t.n_bad_pl > 0 {
            eprintln!("joint: WARNING: {f}: {} 'ok' rows without numeric PLs (treated as missing)", t.n_bad_pl);
        }
        for n in &t.order {
            if seen.insert(n.clone()) {
                loci.push(n.clone());
            }
        }
        tables[tip_ix[colony_stem(f).as_str()]] = Some(t);
    }
    let has_file: Vec<bool> = tables.iter().map(|t| t.is_some()).collect();
    // the shared-fraction grid must be the same in every file (empty = old format, no profile)
    let fracs: Vec<f64> = tables.iter().flatten().next().map(|t| t.fracs.clone()).unwrap_or_default();
    for (f, t) in args.genotype_files.iter().zip(tables.iter().filter_map(|t| t.as_ref())) {
        if t.fracs != fracs {
            return Err(bad(format!("joint: {f}: pl_f profile columns {:?} differ from the first file's {:?}", t.fracs, fracs)));
        }
    }
    let (l, c) = (loci.len(), tips.len());
    eprintln!(
        "joint: {c} colonies, {l} loci, {} hypotheses on the tree; branch prior {}, root prior {}, dropout {}, false-present {}; NOISE = {}",
        tree.nodes.len(),
        args.branch_prior,
        args.root_prior,
        args.dropout,
        args.false_present,
        if fracs.is_empty() { "absent everywhere (no pl_f profile columns)".to_string() } else { format!("absent everywhere or a shared alt fraction in {fracs:?}") }
    );
    let ctx = Ctx::new(&tree, &args.branch_prior, args.root_prior);

    let mut usable = vec![false; l * c];
    let mut kinds = vec![String::new(); l];
    let mut scores = Vec::with_capacity(l);
    let (mut l1, mut l0) = (vec![0.0; c], vec![0.0; c]);
    let mut lf: Vec<Vec<f64>> = vec![vec![0.0; c]; fracs.len()];
    let s_phred = -LN10 / 10.0;
    for (i, name) in loci.iter().enumerate() {
        for col in 0..c {
            let row = tables[col].as_ref().and_then(|t| t.rows.get(name));
            if kinds[i].is_empty() {
                if let Some(r) = row {
                    kinds[i] = r.kind.clone();
                }
            }
            (l1[col], l0[col]) = match row.and_then(|r| r.pl) {
                Some(pl) => {
                    usable[i * c + col] = true;
                    let (p1, p0) = pl_cell(pl);
                    with_errors(p1, p0, args.dropout, args.false_present)
                }
                None => (0.0, 0.0),
            };
            // profile: the pl_f values are on the pl scale (relative to the best dosage, whose pl is
            // 0, so pl_cell's normalisation is the identity); a row without a profile is 0 under every φ
            for k in 0..fracs.len() {
                lf[k][col] = match row {
                    Some(r) if r.pl.is_some() && r.pl_frac.len() == fracs.len() => s_phred * r.pl_frac[k],
                    _ => 0.0,
                };
            }
        }
        scores.push(score_locus(&ctx, &l1, &l0, &lf));
    }

    let (tsv, mat) = render(&ctx, &loci, &kinds, &scores, &usable, &has_file);
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
    let n_with = has_file.iter().filter(|&&f| f).count();
    eprintln!("joint: tips with data {n_with}, tips without {}, loci {l}", c - n_with);
    Ok(())
}

// ------------------------------------------------------------------ tests

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::OUTPUT_HEADER;

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
        score_locus(ctx, &l1, &l0, &[])
    }

    const HET: Option<[f64; 3]> = Some([200.0, 0.0, 150.0]);
    const WT: Option<[f64; 3]> = Some([0.0, 60.0, 300.0]);

    #[test]
    fn small_numerics() {
        let lf = ln_factorials(10);
        assert!((lf[5] - 120f64.ln()).abs() < 1e-12);
        assert!((betaln_int(3, 4, &lf) - (1.0f64 / 60.0).ln()).abs() < 1e-12);
        assert!((betaln_int(1, 1, &lf)).abs() < 1e-15);
        assert!((lae(NEG_INF, -1.0) + 1.0).abs() < 1e-15);
        assert!((lse([0.0, 0.0]) - LN2).abs() < 1e-15);
        assert_eq!(argmax(&[1.0, 3.0, 3.0]), 1);
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
        // B: no usable row (status != ok)
        let s = score_pl(&ctx, &[HET, None, HET, WT, WT, WT, WT, WT]);
        assert_eq!(s.best, Best::Node(1));
        assert!(s.p_carrier[1] > 0.95, "{:?}", s.p_carrier);
        // B: low depth, weakly favouring the reference (one ref read): imputed carrier
        let s = score_pl(&ctx, &[HET, Some([0.0, 3.0, 6.0]), HET, WT, WT, WT, WT, WT]);
        assert_eq!(s.best, Best::Node(1));
        assert!(s.p_carrier[1] > 0.95, "{:?}", s.p_carrier);
        assert!(s.p_carrier[3] < 1e-3, "{:?}", s.p_carrier);
    }

    #[test]
    fn scattered_pattern_is_not_a_tree_event() {
        let t = Tree::parse(TOY).unwrap();
        let ctx = ctx_for(&t);
        let s = score_pl(&ctx, &[HET, WT, HET, WT]);
        assert_eq!(s.best, Best::Indep);
        assert!(s.bf < -5.0, "bf {}", s.bf);
        assert!(s.p_carrier[0] > 0.99 && s.p_carrier[2] > 0.99);
        assert!(s.p_carrier[1] < 0.01 && s.p_carrier[3] < 0.01);
    }

    /// brute force over all 2^C presence vectors: P(z | INDEP) = B(k+1, C-k+1)
    fn brute(l1: &[f64], l0: &[f64]) -> (f64, Vec<f64>) {
        let c = l1.len();
        let lf = ln_factorials(c + 1);
        let mut tot = NEG_INF;
        let mut pres = vec![NEG_INF; c];
        for z in 0u32..(1 << c) {
            let k = z.count_ones() as usize;
            let mut ll = betaln_int(k + 1, c - k + 1, &lf);
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
    fn error_terms_rescue_root_from_one_dropout() {
        // 32 tips (balanced), germline locus: 31 colonies strongly present, one hard "absent"
        // (PL 40). INDEP's ∫π^31(1-π) = 1/(32·33) is lenient, so without an error term the one
        // dropout costs ROOT 4 log10 units and INDEP wins; at ε₁ = 0.02 it costs 1.7 and ROOT holds.
        // (With 8 tips INDEP still wins at 0.02: the prior over colony counts is flatter there.)
        fn balanced(lo: usize, hi: usize) -> String {
            if hi - lo == 1 {
                return format!("T{lo}");
            }
            let mid = (lo + hi) / 2;
            format!("({}:1,{}:1)", balanced(lo, mid), balanced(mid, hi))
        }
        let n = 32;
        let tree = Tree::parse(&format!("{};", balanced(0, n))).unwrap();
        assert_eq!(tree.tips.len(), n);
        let ctx = Ctx::new(&tree, "length", 0.1);
        let present = pl_cell([60.0, 0.0, 30.0]);
        let absent = pl_cell([0.0, 40.0, 80.0]);
        let mut l1 = vec![present.0; n];
        let mut l0 = vec![present.1; n];
        l1[3] = absent.0;
        l0[3] = absent.1;
        let plain = score_locus(&ctx, &l1, &l0, &[]);
        assert_eq!(plain.best, Best::Indep, "without an error term one dropout hands a germline locus to INDEP");
        let mixed: Vec<(f64, f64)> = (0..n).map(|c| with_errors(l1[c], l0[c], 0.02, 0.0)).collect();
        let l1m: Vec<f64> = mixed.iter().map(|x| x.0).collect();
        let l0m: Vec<f64> = mixed.iter().map(|x| x.1).collect();
        let fixed = score_locus(&ctx, &l1m, &l0m, &[]);
        assert_eq!(fixed.best, Best::Node(0), "dropout 0.02 keeps it at ROOT");
        // ε = 0 is the identity; the mixture never exceeds the larger term
        assert_eq!(with_errors(-1.0, -5.0, 0.0, 0.0), (-1.0, -5.0));
        let (a, b) = with_errors(-9.0, -0.1, 0.02, 0.0);
        assert!(a > -9.0 && a < -0.1 && (b - -0.1).abs() < 1e-12);
    }

    #[test]
    fn profile_columns_parsed_and_grid_checked() {
        let text = "locus\tkind\tstatus\tpl_absent\tpl_het\tpl_hom\tpl_f020\tpl_f050\tpl_f200\n\
                    chr1:1-13\tTSD\tok\t12\t0\t40\t-3\t5\t20\n\
                    chr1:50-50\tTSD\tno_reads\t0\t0\t0\t0\t0\t0\n";
        let t = read_colony_table_from("mem", Box::new(io::Cursor::new(text.as_bytes().to_vec()))).unwrap();
        assert_eq!(t.fracs, vec![0.02, 0.05, 0.2]);
        assert_eq!(t.rows["chr1:1-13"].pl_frac, vec![-3.0, 5.0, 20.0]);
        assert!(t.rows["chr1:50-50"].pl_frac.is_empty());
    }

    #[test]
    fn shared_fraction_noise_absorbs_a_diffuse_locus() {
        // 16 tips; every colony shows a weak alt signal (PL: absent 4, het 0 -- a 2-of-30 read
        // locus from mismapped reads), two sisters a little stronger. Without the profile the two
        // sisters make a "clade"; with a shared φ = 0.05 explaining every colony, NOISE wins.
        fn balanced(lo: usize, hi: usize) -> String {
            if hi - lo == 1 {
                return format!("T{lo}");
            }
            let mid = (lo + hi) / 2;
            format!("({}:1,{}:1)", balanced(lo, mid), balanced(mid, hi))
        }
        let n = 16;
        let tree = Tree::parse(&format!("{};", balanced(0, n))).unwrap();
        let ctx = Ctx::new(&tree, "length", 0.1);
        let s = -LN10 / 10.0;
        // weak cells: PL [absent 4, het 0, hom 40], profile at φ=0.05 is the best fit (pl_f -6)
        let weak = pl_cell([4.0, 0.0, 40.0]);
        // the two sisters T0, T1: PL [absent 15, het 0, hom 30], φ=0.05 fits less well (pl_f 6)
        let strong = pl_cell([15.0, 0.0, 30.0]);
        let mut l1 = vec![weak.0; n];
        let mut l0 = vec![weak.1; n];
        let mut lf = vec![vec![s * -6.0; n]];
        for c in 0..2 {
            l1[c] = strong.0;
            l0[c] = strong.1;
            lf[0][c] = s * 6.0;
        }
        let plain = score_locus(&ctx, &l1, &l0, &[]);
        assert!(matches!(plain.best, Best::Node(_)) || plain.best == Best::Indep, "{:?}", plain.best);
        let with = score_locus(&ctx, &l1, &l0, &lf);
        assert_eq!(with.best, Best::Noise, "a shared alt fraction explains the diffuse signal");
        assert!(with.bf < plain.bf);
        // a clean clade (two sisters het at GQ 99, everyone else absent at PL 60) is untouched by the
        // profile: the shared φ fits the 14 clean absents far worse than "absent"
        let het = pl_cell([99.0, 0.0, 60.0]);
        let abs = pl_cell([0.0, 60.0, 120.0]);
        let l1c: Vec<f64> = (0..n).map(|c| if c < 2 { het.0 } else { abs.0 }).collect();
        let l0c: Vec<f64> = (0..n).map(|c| if c < 2 { het.1 } else { abs.1 }).collect();
        let lfc = vec![(0..n).map(|c| if c < 2 { s * 60.0 } else { s * 20.0 }).collect::<Vec<f64>>()];
        let a = score_locus(&ctx, &l1c, &l0c, &[]);
        let b = score_locus(&ctx, &l1c, &l0c, &lfc);
        assert_eq!(a.best, b.best);
        assert!((a.bf - b.bf).abs() < 0.3, "{} vs {}", a.bf, b.bf);
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
            let (li, lp) = indep(&l1, &l0, &ctx.bl_full, &ctx.bl_v);
            let (bt, bp) = brute(&l1, &l0);
            assert!((li - bt).abs() < 1e-10, "{li} vs {bt}");
            for j in 0..4 {
                assert!((lp[j] - bp[j]).abs() < 1e-10, "colony {j}: {} vs {}", lp[j], bp[j]);
            }
        }
        // all missing: P(present) = E[pi] = 1/2
        let (li, lp) = indep(&[0.0; 4], &[0.0; 4], &ctx.bl_full, &ctx.bl_v);
        assert!(li.abs() < 1e-12);
        assert!((lp[2] - 0.5f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn branch_prior_length_and_uniform() {
        let t = Tree::parse("((A:0,B:3):1,C:4);").unwrap();
        let p: Vec<f64> = branch_log_prior(&t, "length", 0.1).iter().map(|x| x.exp()).collect();
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
    fn csv_and_number_formatting() {
        assert_eq!(csv_field("chr1:100-112"), "chr1:100-112");
        assert_eq!(csv_field("a;b"), "\"a;b\"");
        assert_eq!(csv_field("say \"x\""), "\"say \"\"x\"\"\"");
        assert_eq!(fmt_f(-0.00001, 3), "0.000");
        assert_eq!(fmt_f(1.23456, 3), "1.235");
        assert_eq!(fmt_f(0.99996, 4), "1.0000");
        assert_eq!(fmt_f(NEG_INF, 3), "-inf");
        assert_eq!(colony_stem("/x/y/S1.txt.gz"), "S1");
        assert_eq!(colony_stem("PD1234b.genotypes.tsv"), "PD1234b");
    }

    #[test]
    fn header_driven_parsing() {
        // columns reordered + an extra one; non-ok rows have empty p_* and zero PLs
        let text = "status\tpl_hom\textra\tlocus\tpl_het\tkind\tp_absent\tpl_absent\n\
                    ok\t150\tx\tchr1:1-13\t0\tTSD\t0.0000\t200\n\
                    no_reads\t0\tx\tchr1:50-50\t0\tBLUNT\t\t0\n\
                    high_coverage\t0\tx\tchr1:90-60\t0\tL1_MED_DELETION\t\t0\n\
                    ok\t.\tx\tchr1:200-212\t3\tTSD\t\t0\n\
                    ok\t1\tx\tchr1:1-13\t2\tTSD\t\t3\n";
        let t = read_colony_table_from("mem", Box::new(io::Cursor::new(text.as_bytes().to_vec()))).unwrap();
        assert_eq!(t.order, vec!["chr1:1-13", "chr1:50-50", "chr1:90-60", "chr1:200-212"]);
        assert_eq!(t.rows["chr1:1-13"], Row { kind: "TSD".into(), pl: Some([200.0, 0.0, 150.0]), pl_frac: vec![] }); // first wins
        assert!(t.fracs.is_empty());
        assert_eq!(t.rows["chr1:50-50"].pl, None);
        assert_eq!(t.rows["chr1:90-60"].kind, "L1_MED_DELETION");
        assert_eq!(t.rows["chr1:90-60"].pl, None);
        assert_eq!(t.rows["chr1:200-212"].pl, None);
        assert_eq!(t.n_bad_pl, 1);
        // the real header parses; a legacy header is refused
        let real = format!("{OUTPUT_HEADER}chr1:1-13\tTSD\tok\t20\t10\t10\t0\t0\t0\t5\t5\t0.5\t0\t1\t0\t200\t0\t150\t99\t300\t300\n");
        let t = read_colony_table_from("mem", Box::new(io::Cursor::new(real.into_bytes()))).unwrap();
        assert_eq!(t.rows["chr1:1-13"].pl, Some([200.0, 0.0, 150.0]));
        let legacy = "insertion\tgenotype\tscore_genotype\tscore_alternative\tcoverage\tn_alt\tn_ref\tn_art\n";
        let e = read_colony_table_from("old.txt.gz", Box::new(io::Cursor::new(legacy.as_bytes().to_vec()))).unwrap_err();
        assert!(e.to_string().contains("missing column(s) locus, status, pl_absent, pl_het, pl_hom"), "{e}");
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
    fn end_to_end_numeric_files() {
        let dir = std::env::temp_dir().join(format!("peartree_joint_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("t.nwk"), TOY).unwrap();
        let ok = |l: &str, pl: [u32; 3]| {
            format!("{l}\tTSD\tok\t20\t10\t10\t0\t0\t0\t5\t5\t0.50\t0\t1\t0\t{}\t{}\t{}\t99\t300\t300\n", pl[0], pl[1], pl[2])
        };
        let noreads = |l: &str| format!("{l}\tTSD\tno_reads\t0\t0\t0\t0\t0\t0\t0\t0\t0.00\t\t\t\t0\t0\t0\t0\t0\t0\n");
        let (het, wt) = ([200, 0, 150], [0, 60, 300]);
        let (l1, l2) = ("chr1:100-112", "chr1:500-515");
        let files = [
            ("A", format!("{OUTPUT_HEADER}{}{}", ok(l1, het), ok(l2, het))),
            ("B", format!("{OUTPUT_HEADER}{}{}", noreads(l1), ok(l2, wt))),
            ("C", format!("{OUTPUT_HEADER}{}{}", ok(l1, wt), ok(l2, het))),
            ("D", format!("{OUTPUT_HEADER}{}", ok(l1, wt))), // l2 row absent
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
            dropout: 0.0,
            false_present: 0.0,
        };
        run(&args).unwrap();
        let tsv = std::fs::read_to_string(&out_tsv).unwrap();
        let lines: Vec<&str> = tsv.lines().collect();
        assert!(lines[0].starts_with("locus\tlocus_kind\tn_colonies_data\tbest\tcarriers"));
        assert!(lines[0].ends_with("\tp_A\tp_B\tp_C\tp_D"));
        let f1: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(&f1[..6], &[l1, "TSD", "3", "N1", "A,B", "2"]);
        let f2: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(&f2[..4], &[l2, "TSD", "3", "INDEP"]); // A and C carry it, B does not: scattered
        let mat = read_gz(&out_mat);
        let m: Vec<&str> = mat.lines().collect();
        assert_eq!(m[0], ";A;B;C;D");
        let c1: Vec<&str> = m[1].split(';').collect();
        assert_eq!(c1[0], l1);
        assert_eq!(c1[2], ""); // B: no_reads -> empty cell
        assert_eq!(c1[1], f1[15]); // same number as the TSV's p_A
        assert!(c1[1].parse::<f64>().unwrap() > 0.99 && c1[3].parse::<f64>().unwrap() < 0.01);
        assert_eq!(c1[1].split('.').nth(1).unwrap().len(), 4);
        let c2: Vec<&str> = m[2].split(';').collect();
        assert_eq!(c2.len(), 5);
        assert_eq!(c2[4], ""); // D: row absent -> empty cell
        assert!(!mat.contains("wild-type") && !mat.contains("insertion"));
        // a tip without a genotype file: warning, all-missing colony, empty columns
        let mut args3 = JointArgs { genotype_files: paths[..3].to_vec(), ..args };
        run(&args3).unwrap();
        let tsv = std::fs::read_to_string(&out_tsv).unwrap();
        let r1: Vec<&str> = tsv.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(r1.len(), 19);
        assert_eq!(&r1[..6], &[l1, "TSD", "2", "N1", "A,B", "2"]);
        assert_eq!(r1[18], ""); // p_D empty
        let mat = read_gz(&out_mat);
        assert_eq!(mat.lines().next().unwrap(), ";A;B;C;D");
        assert!(mat.lines().skip(1).all(|r| r.ends_with(';')), "{mat}");
        // a genotype file whose stem is not a tip is an error
        let extra = dir.join("S11.txt.gz");
        write_gz(&extra, &files[0].1);
        args3.genotype_files = paths.clone();
        args3.genotype_files.push(extra.to_string_lossy().into_owned());
        let e = run(&args3).unwrap_err().to_string();
        assert!(e.contains("not a tree tip (1)") && e.contains("S11.txt.gz"), "{e}");
        std::fs::remove_dir_all(&dir).ok();
    }
}
