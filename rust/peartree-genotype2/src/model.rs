//! Genotype likelihoods, posterior, GQ, VAF (SPEC "Genotype model"). No call strings: the
//! numbers are the output.
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
use crate::types::{AltSide, Call, ReadClass, ReadObs};

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
#[allow(dead_code)] // `best` is kept for debugging
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

/// Summarise one locus from its read observations. `n_disc` discordant anchors add
/// `cfg.disc_weight_nats` each towards alt (0 by default). The driver decides `status`
/// (`no_reads` when `obs` is empty, `high_coverage`, `error`) and formats the row.
///
/// Uninformative reads stay IN the likelihood: one carries little (|llr| below the reporting
/// threshold) but many add up -- at a far duplication every reference-junction read has
/// llr = ln(0.5) by construction and together they are the only evidence of absence. They are
/// not counted as votes. Unexplained reads (chimeras / mismaps) are counted (`n_art`) and
/// excluded from the likelihood.
pub fn call_locus(obs: &[ReadObs], n_disc: i64, cfg: &Config) -> Call {
    let (mut n_alt, mut n_ref, mut n_art, mut n_uninf) = (0i64, 0i64, 0i64, 0i64);
    let (mut n_alt_l, mut n_alt_r) = (0i64, 0i64);
    let (mut score_alt, mut score_ref) = (0i64, 0i64);
    let mut d: Vec<f64> = Vec::with_capacity(obs.len());
    for o in obs {
        match o.class {
            ReadClass::Uninformative => {
                n_uninf += 1;
                d.push(clamp_llr(o.llr()));
            }
            ReadClass::Unexplained => n_art += 1,
            ReadClass::Alt => {
                let x = clamp_llr(o.llr());
                n_alt += 1;
                match o.alt_side {
                    AltSide::Left => n_alt_l += 1,
                    AltSide::Right => n_alt_r += 1,
                    AltSide::Both => {
                        n_alt_l += 1;
                        n_alt_r += 1;
                    }
                    AltSide::None => {}
                }
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
    if d.is_empty() && n_disc == 0 {
        // nothing enters the likelihood: the posterior is the prior, PL/GQ carry no information
        return Call {
            n_alt, n_ref, n_art, n_uninf, n_alt_l, n_alt_r,
            vaf: 0.0,
            post: normalised_prior(cfg),
            pl: [0, 0, 0],
            gq: 0,
            score_alt, score_ref,
        };
    }
    let f = fit(&d, n_disc, cfg);
    Call {
        n_alt, n_ref, n_art, n_uninf, n_alt_l, n_alt_r,
        vaf: vaf_mle(&d),
        post: f.post,
        pl: f.pl,
        gq: f.gq,
        score_alt, score_ref,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(class: ReadClass, llr: f64) -> ReadObs {
        ReadObs { ll_ref: -40.0, ll_alt: -40.0 + llr, class, explained_frac: 1.0, crosses_junction: true, alt_side: AltSide::None }
    }
    fn alt(llr: f64) -> ReadObs {
        read(ReadClass::Alt, llr)
    }
    fn refr(llr: f64) -> ReadObs {
        read(ReadClass::Ref, llr)
    }
    fn reads(na: usize, la: f64, nr: usize, lr: f64) -> Vec<ReadObs> {
        let mut v = vec![alt(la); na];
        v.extend(vec![refr(lr); nr]);
        v
    }
    fn call(obs: &[ReadObs]) -> Call {
        call_locus(obs, 0, &Config::default())
    }
    fn present(c: &Call) -> f64 {
        c.post[1] + c.post[2]
    }

    #[test]
    fn clean_wild_type() {
        let c = call(&reads(0, 0.0, 8, -15.0));
        assert!(c.post[0] > 0.98, "{:?}", c.post);
        assert_eq!(c.pl[0], 0);
        assert!(c.pl[1] > 10 && c.pl[2] > 100);
        assert_eq!(c.vaf, 0.0);
        assert_eq!((c.n_alt, c.n_ref), (0, 8));
        assert_eq!(c.score_ref, 8 * 65);
        assert_eq!(c.score_alt, 0);
    }

    #[test]
    fn one_ref_read_is_weak_absence_evidence() {
        // a het alternative explains a ref read ~half the time: P(absent) ≈ 0.63 under the flat prior
        let c = call(&reads(0, 0.0, 1, -15.0));
        assert!(c.post[0] > 0.6 && c.post[0] < 0.7, "{:?}", c.post);
        assert!(c.gq < 10);
    }

    #[test]
    fn balanced_het() {
        let c = call(&reads(6, 15.0, 6, -15.0));
        assert!(c.post[1] > 0.95, "{:?}", c.post);
        assert!((c.vaf - 0.5).abs() < 1e-9);
        assert!(c.gq >= 10);
        assert_eq!(c.pl[1], 0);
    }

    #[test]
    fn purity_diluted_het_is_still_het() {
        let c = call(&reads(4, 15.0, 8, -15.0));
        assert!(c.post[1] > 0.99, "{:?}", c.post);
        assert!((c.vaf - 1.0 / 3.0).abs() < 0.01);
    }

    #[test]
    fn hom_and_presence_certain_with_unclear_zygosity() {
        let c = call(&reads(8, 15.0, 0, 0.0));
        assert!(c.post[2] > 0.99, "{:?}", c.post);
        let c2 = call(&reads(2, 15.0, 0, 0.0));
        assert!(present(&c2) > 0.999 && c2.post[2] < 0.9, "{:?}", c2.post);
    }

    #[test]
    fn one_alt_read_is_strong_presence_evidence() {
        // bg_alt_rate 0.005: a single read just above the informative threshold is ~85x likelier
        // under presence -- downstream must weigh this with depth / the joint step
        let c = call(&[alt(20.0)]);
        assert!(c.post[0] < 0.01, "{:?}", c.post);
        let c5 = call(&[alt(5.0)]);
        assert!(c5.post[0] < 0.02, "{:?}", c5.post);
    }

    #[test]
    fn low_vaf_is_ambiguous() {
        // 3 alt + 30 ref (VAF 0.09): P(absent) ≈ 0.78 -- neither a confident absence nor presence
        let c = call(&reads(3, 15.0, 30, -15.0));
        assert!(c.post[0] > 0.7 && c.post[0] < 0.85, "{:?}", c.post);
        assert!((c.vaf - 0.09).abs() < 1e-9);
        assert_eq!(c.score_ref, 30 * 65);
        assert_eq!(c.score_alt, 3 * 65);
    }

    #[test]
    fn unexplained_reads_are_counted_not_modelled() {
        let mut v = vec![read(ReadClass::Unexplained, 0.0); 3];
        v.push(alt(15.0));
        let c = call(&v);
        assert_eq!((c.n_art, c.n_alt), (3, 1));
        assert!(present(&c) > 0.99); // the one alt read is all the model sees
    }

    #[test]
    fn uninformative_zero_llr_reads_leave_the_prior() {
        let c = call(&vec![read(ReadClass::Uninformative, 0.0); 5]);
        assert_eq!(c.n_uninf, 5);
        assert!((c.post[0] - 1.0 / 3.0).abs() < 1e-6, "{:?}", c.post);
        assert!(c.gq <= 2); // the prior's own GQ: -10 log10(2/3)
    }

    #[test]
    fn far_dup_reference_reads_add_up() {
        // far duplication: every ref-junction read has llr = ln 0.5 and is Uninformative alone
        let c = call(&vec![read(ReadClass::Uninformative, 0.5f64.ln()); 12]);
        assert!(c.post[0] > 0.9, "{:?}", c.post);
        assert_eq!((c.n_alt, c.n_ref, c.n_uninf), (0, 0, 12));
    }

    #[test]
    fn empty_is_prior() {
        let c = call(&[]);
        assert_eq!(c.pl, [0, 0, 0]);
        assert_eq!(c.gq, 0);
        assert_eq!(c.vaf, 0.0);
    }

    #[test]
    fn disc_pseudo_reads_move_the_posterior_only_when_weighted() {
        let cfg = Config::default();
        let one_ref = reads(0, 0.0, 1, -15.0);
        let a = call_locus(&one_ref, 3, &cfg);
        let b = call_locus(&one_ref, 0, &cfg);
        assert_eq!(a.post, b.post);
        let w = Config { disc_weight_nats: 10.0, ..Config::default() };
        let c = call_locus(&one_ref, 3, &w);
        assert!(c.post[1] > 0.9, "{:?}", c.post);
    }

    #[test]
    fn per_junction_alt_counts() {
        let mut l = alt(15.0);
        l.alt_side = AltSide::Left;
        let mut r = alt(15.0);
        r.alt_side = AltSide::Right;
        let mut b = alt(15.0);
        b.alt_side = AltSide::Both;
        let c = call(&[l, l, r, b]);
        assert_eq!((c.n_alt, c.n_alt_l, c.n_alt_r), (4, 3, 2));
    }

    /// Prints P(absent)/GQ for n ref-only and n alt-only reads (threshold review).
    #[test]
    #[ignore]
    fn sweep() {
        for n in 0..=12 {
            let c = call(&reads(0, 0.0, n, -15.0));
            println!("ref x{n}: post {:?} gq {} pl {:?}", c.post, c.gq, c.pl);
        }
        for n in 0..=12 {
            let c = call(&reads(n, 15.0, 0, 0.0));
            println!("alt x{n}: post {:?} gq {} pl {:?}", c.post, c.gq, c.pl);
        }
    }
}
