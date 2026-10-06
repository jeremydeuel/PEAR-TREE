//! Genotype likelihoods, posterior, vocabulary mapping (SPEC "Genotype model"). Owner: C.
//!
//! Model. Dosage `g ∈ {0,1,2}`; colony purity `p` on `cfg.purity_grid` (uniform weights); the
//! alt-haplotype fraction of the reads is
//!
//! ```text
//!   φ(0,p) = bg                 φ(1,p) = max(p/2, bg)          φ(2,p) = max(p, 1 - bg)
//! ```
//!
//! (`bg = cfg.bg_alt_rate`). A read with likelihoods `ll_ref`, `ll_alt` contributes
//! `ll_r(φ) = logaddexp(ln(1-φ) + ll_ref, ln φ + ll_alt)`. Everything is computed on the
//! *difference* `d = ll_alt - ll_ref`: `ll_r(φ) = ll_ref + ln((1-φ) + φ e^d)`, and `Σ ll_ref`
//! is the same for every φ, so it cancels in PL, the posterior and the VAF argmax. That keeps
//! the arithmetic well-scaled however negative the absolute alignment log-likelihoods are.
//!
//! `GL_g = logmeanexp_p Σ_r ll_r(φ(g,p))`; discordant anchors enter as `n_disc` pseudo-reads
//! with `d = cfg.disc_weight_nats` (a pseudo-read with `d = 0` contributes `ln 1 = 0`, so the
//! default weight changes nothing). PL / posterior / GQ / VAF as in the SPEC; the VAF is the
//! MLE over the real Alt/Ref reads only (pseudo-reads are not reads).

use crate::config::Config;
use crate::types::{
    Call, ReadClass, ReadObs, GT_ARTEFACT, GT_HETEROZYGOUS, GT_HOMOZYGOUS, GT_INSERTION,
    GT_INSERTION_UNCERTAIN, GT_NO_COVERAGE, GT_WILDTYPE, GT_WILDTYPE_UNCERTAIN,
};

/// Phred units per nat: `10·log10(e)`.
const PHRED_PER_NAT: f64 = 10.0 * std::f64::consts::LOG10_E;
/// A read's log-likelihood ratio is clamped to ±this many nats, so a degenerate (infinite)
/// alignment likelihood cannot turn the sums into inf − inf = NaN.
const MAX_ABS_LLR: f64 = 1000.0;
/// VAF grid resolution: φ ∈ {0, 1/VAF_STEPS, …, 1}.
const VAF_STEPS: usize = 100;

fn clamp_llr(d: f64) -> f64 {
    if d.is_nan() {
        0.0
    } else {
        d.clamp(-MAX_ABS_LLR, MAX_ABS_LLR)
    }
}

fn logaddexp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let m = a.max(b);
    m + (-(a - b).abs()).exp().ln_1p()
}

/// `ln mean exp(xs)` (−inf for an empty or all −inf input).
fn logmeanexp(xs: &[f64]) -> f64 {
    let m = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY || xs.is_empty() {
        return f64::NEG_INFINITY;
    }
    let s: f64 = xs.iter().map(|&x| (x - m).exp()).sum();
    m + (s / xs.len() as f64).ln()
}

/// `Σ_r ln((1-φ) + φ e^{d_r})` = `Σ_r ll_r(φ) − Σ_r ll_ref` for log-likelihood ratios `d`,
/// plus `n_pseudo` pseudo-reads of ratio `d_pseudo`.
fn sum_ll(d: &[f64], n_pseudo: i64, d_pseudo: f64, phi: f64) -> f64 {
    let l_ref = (1.0 - phi).ln(); // −inf at φ = 1
    let l_alt = phi.ln(); // −inf at φ = 0
    let mut s: f64 = d.iter().map(|&x| logaddexp(l_ref, l_alt + x)).sum();
    if n_pseudo > 0 && d_pseudo != 0.0 {
        s += n_pseudo as f64 * logaddexp(l_ref, l_alt + d_pseudo);
    }
    s
}

/// Alt-haplotype read fraction for dosage `g` at purity `p`.
fn phi(g: usize, p: f64, bg: f64) -> f64 {
    match g {
        0 => bg,
        1 => (p / 2.0).max(bg),
        _ => p.max(1.0 - bg),
    }
}

/// The likelihood part of a call: GL per dosage, posterior, PL, GQ, best dosage.
struct Fit {
    post: [f64; 3],
    pl: [i32; 3],
    gq: i32,
    best: usize,
}

fn fit(d: &[f64], n_disc: i64, cfg: &Config) -> Fit {
    let bg = cfg.bg_alt_rate;
    let d_pseudo = clamp_llr(cfg.disc_weight_nats);
    let mut gl = [0.0f64; 3];
    for (g, slot) in gl.iter_mut().enumerate() {
        let per_p: Vec<f64> = cfg
            .purity_grid
            .iter()
            .map(|&p| sum_ll(d, n_disc, d_pseudo, phi(g, p, bg)))
            .collect();
        *slot = logmeanexp(&per_p);
    }
    let max_gl = gl.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut pl = [0i32; 3];
    for g in 0..3 {
        // `as i32` saturates (and maps NaN to 0)
        pl[g] = (-PHRED_PER_NAT * (gl[g] - max_gl)).round().max(0.0) as i32;
    }
    // posterior ∝ e^{GL_g} · prior_g
    let mut lp = [f64::NEG_INFINITY; 3];
    for g in 0..3 {
        let pr = cfg.prior[g];
        if pr > 0.0 && pr.is_finite() {
            lp[g] = gl[g] + pr.ln();
        }
    }
    let m = lp.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let post = if m.is_finite() {
        let e = [(lp[0] - m).exp(), (lp[1] - m).exp(), (lp[2] - m).exp()];
        let z = e[0] + e[1] + e[2];
        [e[0] / z, e[1] / z, e[2] / z]
    } else {
        [1.0 / 3.0; 3]
    };
    let mut best = 0;
    for g in 1..3 {
        if post[g] > post[best] {
            best = g;
        }
    }
    // 1 - post_best as the sum of the others: no cancellation when post_best -> 1
    let rest: f64 = (0..3).filter(|&g| g != best).map(|g| post[g]).sum();
    let gq = if rest <= 0.0 { 99 } else { (-10.0 * rest.log10()).round().clamp(0.0, 99.0) as i32 };
    Fit { post, pl, gq, best }
}

/// MLE of φ on the 0..1 grid (step 1/VAF_STEPS) over the real reads; 0.0 without reads.
fn vaf_mle(d: &[f64]) -> f64 {
    if d.is_empty() {
        return 0.0;
    }
    let mut best_phi = 0.0;
    let mut best_ll = f64::NEG_INFINITY;
    for k in 0..=VAF_STEPS {
        let phi = k as f64 / VAF_STEPS as f64;
        let ll = sum_ll(d, 0, 0.0, phi);
        if ll > best_ll {
            best_ll = ll;
            best_phi = phi;
        }
    }
    best_phi
}

fn normalised_prior(cfg: &Config) -> [f64; 3] {
    let p = cfg.prior.map(|x| if x > 0.0 && x.is_finite() { x } else { 0.0 });
    let z: f64 = p.iter().sum();
    if z > 0.0 { p.map(|x| x / z) } else { [1.0 / 3.0; 3] }
}

/// Call one locus from its read observations. `n_disc` discordant anchors add
/// `cfg.disc_weight_nats` each towards alt (0 by default). Never returns `high-coverage` /
/// `error` (the driver decides those before calling).
// live once driver.rs (owner D) calls it; remove at integration
#[allow(dead_code)]
pub fn call_locus(obs: &[ReadObs], n_disc: i64, cfg: &Config) -> Call {
    let (mut n_alt, mut n_ref, mut n_art, mut n_uninf) = (0i64, 0i64, 0i64, 0i64);
    let (mut score_alt, mut score_ref) = (0i64, 0i64);
    let mut d: Vec<f64> = Vec::with_capacity(obs.len());
    for o in obs {
        match o.class {
            ReadClass::Uninformative => n_uninf += 1,
            ReadClass::Unexplained => n_art += 1,
            ReadClass::Alt => {
                let x = clamp_llr(o.llr());
                n_alt += 1;
                score_alt += (PHRED_PER_NAT * x).round() as i64;
                d.push(x);
            }
            ReadClass::Ref => {
                let x = clamp_llr(o.llr());
                n_ref += 1;
                score_ref += (PHRED_PER_NAT * -x).round() as i64;
                d.push(x);
            }
        }
    }

    let no_coverage = || Call {
        genotype: GT_NO_COVERAGE,
        score_genotype: 0,
        score_alternative: 0,
        n_alt,
        n_ref,
        n_art,
        n_uninf,
        vaf: 0.0,
        gq: 0,
        pl: [0, 0, 0],
        post: normalised_prior(cfg),
    };

    // no Alt/Ref/Unexplained reads at all
    if n_alt + n_ref + n_art == 0 {
        return no_coverage();
    }

    let f = fit(&d, n_disc, cfg);
    let vaf = vaf_mle(&d);
    let mk = |genotype: &'static str, score_genotype: i64, score_alternative: i64| Call {
        genotype,
        score_genotype,
        score_alternative,
        n_alt,
        n_ref,
        n_art,
        n_uninf,
        vaf,
        gq: f.gq,
        pl: f.pl,
        post: f.post,
    };

    // artefact-dominated locus (PL/GQ/VAF still describe the Alt/Ref reads, for the joint step)
    let total = n_alt + n_ref + n_art;
    if n_art >= cfg.min_artefact_reads && n_art as f64 >= cfg.artefact_read_fraction * total as f64 {
        return mk(GT_ARTEFACT, 0, score_alt.max(score_ref));
    }
    // Unexplained reads below the artefact rule and nothing informative: the model has no data
    // (its posterior would be the prior), so this is no-coverage as in the legacy summariser.
    if n_alt + n_ref == 0 {
        return no_coverage();
    }

    let present = f.post[1] + f.post[2];
    let wt = |g| mk(g, score_ref, score_alt);
    let ins = |g| mk(g, score_alt, score_ref);
    match f.best {
        0 if f.post[0] >= cfg.p_confident => wt(GT_WILDTYPE),
        0 => wt(GT_WILDTYPE_UNCERTAIN),
        1 if f.post[1] >= cfg.p_confident => ins(GT_HETEROZYGOUS),
        2 if f.post[2] >= cfg.p_confident => ins(GT_HOMOZYGOUS),
        _ if present >= cfg.p_present_certain => ins(GT_INSERTION),
        _ if present >= cfg.p_present_uncertain => ins(GT_INSERTION_UNCERTAIN),
        _ => wt(GT_WILDTYPE_UNCERTAIN),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(class: ReadClass, llr: f64) -> ReadObs {
        ReadObs { ll_ref: -40.0, ll_alt: -40.0 + llr, class, explained_frac: 1.0, crosses_junction: true }
    }
    fn alt(llr: f64) -> ReadObs {
        read(ReadClass::Alt, llr)
    }
    fn refr(llr: f64) -> ReadObs {
        read(ReadClass::Ref, llr)
    }
    fn reads(n_alt: usize, alt_llr: f64, n_ref: usize, ref_llr: f64) -> Vec<ReadObs> {
        let mut v = vec![alt(alt_llr); n_alt];
        v.extend(vec![refr(ref_llr); n_ref]);
        v
    }
    fn call(obs: &[ReadObs]) -> Call {
        call_locus(obs, 0, &Config::default())
    }

    #[test]
    fn clean_wild_type() {
        let c = call(&reads(0, 0.0, 8, -15.0));
        assert_eq!(c.genotype, GT_WILDTYPE);
        assert_eq!((c.n_ref, c.n_alt), (8, 0));
        assert_eq!(c.pl[0], 0);
        assert!(c.pl[1] > 10 && c.pl[2] > c.pl[1], "{:?}", c.pl);
        assert_eq!(c.vaf, 0.0);
        assert!(c.gq >= 10);
        assert_eq!(c.score_genotype, 8 * 65); // round(15 · 4.3429) = 65 per read
        assert_eq!(c.score_alternative, 0);
    }

    #[test]
    fn single_ref_read_is_uncertain() {
        // GL ∝ 0.995 (g0), mean_p(1 - p/2) = 0.575 (g1), ≈0.00375 (g2: φ ≥ 0.995)
        // → post0 = 0.995 / 1.574 ≈ 0.632 under the flat prior (not 0.67: the purity grid makes
        // a het colony look like 57.5% ref, not 50%).
        let c = call(&[refr(-15.0)]);
        assert_eq!(c.genotype, GT_WILDTYPE_UNCERTAIN);
        assert!((c.post[0] - 0.632).abs() < 0.002, "{:?}", c.post);
        assert_eq!(c.gq, 4);
    }

    #[test]
    fn balanced_het() {
        let c = call(&reads(6, 15.0, 6, -15.0));
        assert_eq!(c.genotype, GT_HETEROZYGOUS);
        assert!((c.vaf - 0.5).abs() < 1e-9);
        assert!(c.gq > 10, "gq {}", c.gq);
        assert_eq!(c.pl[1], 0);
    }

    #[test]
    fn purity_diluted_het() {
        // VAF 0.33: the purity grid down to 0.7 (φ1 = 0.35) explains it as a confident het.
        let c = call(&reads(4, 15.0, 8, -15.0));
        assert_eq!(c.genotype, GT_HETEROZYGOUS);
        assert!((c.vaf - 0.33).abs() < 1e-9, "vaf {}", c.vaf);
        assert!(c.post[1] > 0.99);
    }

    #[test]
    fn homozygous() {
        let c = call(&reads(8, 15.0, 0, 0.0));
        assert_eq!(c.genotype, GT_HOMOZYGOUS);
        assert_eq!(c.vaf, 1.0);
        assert!(c.post[2] > 0.99);
    }

    #[test]
    fn two_alt_reads_insertion() {
        // post2 ≈ 1/(1 + mean_p (p/2)^2) ≈ 0.845 < p_confident, P(present) ≈ 1
        let c = call(&reads(2, 15.0, 0, 0.0));
        assert_eq!(c.genotype, GT_INSERTION);
        assert!(c.post[2] > 0.8 && c.post[2] < 0.9, "{:?}", c.post);
    }

    #[test]
    fn single_alt_read_strong() {
        // GL ∝ 0.005·e^20 (g0) vs 0.425·e^20 (g1) vs ≈e^20 (g2): P(present) ≈ 0.9965
        let c = call(&[alt(20.0)]);
        assert_eq!(c.genotype, GT_INSERTION);
        assert!((c.post[0] - 0.0035).abs() < 0.0005, "{:?}", c.post);
    }

    #[test]
    fn single_alt_read_weak() {
        // llr 5: g0 ∝ 0.995 + 0.005·e^5 = 1.737, g1 ∝ 63.6, g2 ∝ 147.9 → P(present) ≈ 0.992.
        // bg_alt_rate = 0.005 caps a single alt read's weight against absence at ~e^5/1.74.
        let c = call(&[alt(5.0)]);
        assert_eq!(c.genotype, GT_INSERTION);
        assert!((c.post[0] - 0.0082).abs() < 0.0005, "{:?}", c.post);
    }

    #[test]
    fn low_vaf_is_not_confident_wild_type() {
        // 3 alt + 30 ref (VAF 0.09): post0 ≈ 0.78 (bg_alt_rate = 0.005 makes 3 alt reads cost
        // ~16 nats under absence, about as much as 30 ref reads cost under a 70%-pure het)
        let c = call(&reads(3, 15.0, 30, -15.0));
        assert_eq!(c.genotype, GT_WILDTYPE_UNCERTAIN);
        assert!(c.post[0] > 0.7 && c.post[0] < 0.85, "{:?}", c.post);
        assert!(c.gq < 10);
        assert!((c.vaf - 0.09).abs() < 1e-9, "vaf {}", c.vaf);
        // wild-type call: the ref side supports it
        assert_eq!(c.score_genotype, 30 * 65);
        assert_eq!(c.score_alternative, 3 * 65);
    }

    #[test]
    fn artefact_rule() {
        let mut obs = vec![read(ReadClass::Unexplained, 0.0); 3];
        obs.push(alt(15.0));
        let c = call(&obs);
        assert_eq!(c.genotype, GT_ARTEFACT);
        assert_eq!((c.n_art, c.n_alt), (3, 1));
        assert_eq!(c.score_genotype, 0);
        assert_eq!(c.score_alternative, 65);
    }

    #[test]
    fn uninformative_only_is_no_coverage() {
        let c = call(&vec![read(ReadClass::Uninformative, 1.0); 5]);
        assert_eq!(c.genotype, GT_NO_COVERAGE);
        assert_eq!(c.n_uninf, 5);
        assert_eq!((c.score_genotype, c.score_alternative, c.gq, c.pl), (0, 0, 0, [0, 0, 0]));
        assert_eq!(call(&[]).genotype, GT_NO_COVERAGE);
        // one Unexplained read (below min_artefact_reads) and nothing informative
        let c = call(&[read(ReadClass::Unexplained, 0.0)]);
        assert_eq!(c.genotype, GT_NO_COVERAGE);
        assert_eq!(c.n_art, 1);
    }

    #[test]
    fn disc_weight_moves_the_call() {
        let obs = [refr(-15.0)];
        let mut cfg = Config::default();
        // weight 0 (default): discordant anchors are counted elsewhere, change nothing here
        assert_eq!(call_locus(&obs, 5, &cfg), call_locus(&obs, 0, &cfg));
        cfg.disc_weight_nats = 10.0;
        // 1 ref read + 3 anchors at 10 nats: wild-type? -> heterozygous (the ref read rules out
        // hom, the pseudo-reads rule out absence)
        let c = call_locus(&obs, 3, &cfg);
        assert_eq!(c.genotype, GT_HETEROZYGOUS);
        assert_eq!((c.n_ref, c.n_alt), (1, 0));
        assert_eq!(c.vaf, 0.0); // pseudo-reads are not reads
        // and they cannot create coverage
        assert_eq!(call_locus(&[], 3, &cfg).genotype, GT_NO_COVERAGE);
    }

    #[test]
    fn scores_and_hand_computed_case() {
        // alt llr 10, 15; ref llr -20.  score_alt = round(43.43) + round(65.14) = 108,
        // score_ref = round(86.86) = 87.  Independent numpy reimplementation of the SPEC
        // formulas (relative to Σ ll_ref): GL = (14.407, 22.727, 19.404) → PL = (36, 0, 14),
        // post = (0.00024, 0.9650, 0.0348) → heterozygous, GQ = round(-10·log10 0.0350) = 15,
        // VAF MLE = 0.67.
        let obs = [alt(10.0), alt(15.0), refr(-20.0)];
        let c = call(&obs);
        assert_eq!((c.n_alt, c.n_ref), (2, 1));
        assert_eq!(c.genotype, GT_HETEROZYGOUS);
        assert_eq!((c.score_genotype, c.score_alternative), (108, 87));
        assert_eq!(c.pl, [HAND_PL[0], HAND_PL[1], HAND_PL[2]]);
        assert_eq!(c.gq, HAND_GQ);
        assert!((c.vaf - 0.67).abs() < 1e-9);
    }
    const HAND_PL: [i32; 3] = [36, 0, 14];
    const HAND_GQ: i32 = 15;

    /// Threshold review: `cargo test sweep -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn sweep() {
        for (label, a, llr) in [("n_ref (ref-only, llr -15)", false, -15.0), ("n_alt (alt-only, llr +15)", true, 15.0)] {
            println!("\n| {label} | genotype | GQ | post0 | post1 | post2 | PL |");
            println!("|---|---|---|---|---|---|---|");
            for n in 0..=12usize {
                let obs = if a { reads(n, llr, 0, 0.0) } else { reads(0, 0.0, n, llr) };
                let c = call(&obs);
                println!(
                    "| {n} | {} | {} | {:.4} | {:.4} | {:.4} | {} {} {} |",
                    c.genotype, c.gq, c.post[0], c.post[1], c.post[2], c.pl[0], c.pl[1], c.pl[2]
                );
            }
        }
    }
}
