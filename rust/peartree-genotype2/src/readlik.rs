//! Per-read likelihoods under the locus haplotypes (SPEC "Read likelihood"). Owner: B.
//!
//! Each read (reference-forward orientation, as stored in the BAM) is realigned against every
//! segment of the locus model. The band diagonal comes from the read's BAM placement through the
//! segment's genome anchors (leading soft clip subtracted); without a covering anchor the
//! alignment is unbanded. `ll_ref` / `ll_alt` = best segment of each hypothesis, with the
//! far-duplication mixture `ll_alt = logaddexp(ln(1-r) + ll_alt, ln r + ll_ref)` for r > 0.
//!
//! Missing hypotheses: a side without segments has `ll = -inf`; the read is then Uninformative
//! (or Unexplained if the side that exists explains too little of it). A model without any
//! segment yields Uninformative with `explained_frac = 0`.

use crate::align::{align_multi, Alignment};
use crate::config::Config;
use crate::types::{Hyp, LocusModel, ReadClass, ReadInput, ReadObs, Segment};

/// Realign one read against every segment of `model` and classify it.
// allow(dead_code): marks the scorer as a live root until the driver (owner D) calls it.
#[allow(dead_code)]
pub fn score_read(model: &LocusModel, read: &ReadInput, cfg: &Config) -> ReadObs {
    let n = read.seq.len();
    if n == 0 || model.segments.is_empty() {
        return ReadObs {
            ll_ref: f64::NEG_INFINITY,
            ll_alt: f64::NEG_INFINITY,
            class: ReadClass::Uninformative,
            explained_frac: 0.0,
            crosses_junction: false,
        };
    }
    let (lead, trail) = soft_clips(read.cigar);
    let mut best_ref: Option<(Alignment, &Segment)> = None;
    let mut best_alt: Option<(Alignment, &Segment)> = None;
    let mut diags: Vec<i64> = Vec::with_capacity(4);
    for seg in &model.segments {
        diagonals(seg, read, lead, trail, &mut diags);
        let a = align_multi(read.seq, read.qual, &seg.seq, &seg.qual, &diags, cfg);
        let slot = match seg.hyp {
            Hyp::Ref => &mut best_ref,
            Hyp::Alt => &mut best_alt,
        };
        if slot.as_ref().is_none_or(|(b, _)| a.ll > b.ll) {
            *slot = Some((a, seg));
        }
    }

    let ll_ref = best_ref.as_ref().map_or(f64::NEG_INFINITY, |(a, _)| a.ll);
    let ll_alt_raw = best_alt.as_ref().map_or(f64::NEG_INFINITY, |(a, _)| a.ll);
    let r = model.alt_ref_junction_fraction;
    let both = best_ref.is_some() && best_alt.is_some();
    let ll_alt = if both && r > 0.0 {
        logaddexp((1.0 - r).ln() + ll_alt_raw, r.ln() + ll_ref)
    } else {
        ll_alt_raw
    };

    let (best, seg) = match (&best_ref, &best_alt) {
        (Some(rf), Some(al)) => {
            if al.0.ll > rf.0.ll {
                *al
            } else {
                *rf
            }
        }
        (Some(x), None) | (None, Some(x)) => *x,
        (None, None) => unreachable!("non-empty segment list"),
    };
    let explained_frac = best.aligned_bases as f64 / n as f64;
    let crosses_junction = seg.junction_cols.iter().any(|&c| best.seg_start < c && c < best.seg_end);

    let llr = ll_alt - ll_ref;
    let class = if explained_frac < cfg.min_explained_frac {
        ReadClass::Unexplained
    } else if !both {
        ReadClass::Uninformative
    } else if llr > cfg.llr_informative {
        ReadClass::Alt
    } else if llr < -cfg.llr_informative {
        ReadClass::Ref
    } else {
        ReadClass::Uninformative
    };
    ReadObs { ll_ref, ll_alt, class, explained_frac, crosses_junction }
}

/// (leading, trailing) soft-clip lengths; hard clips are skipped.
fn soft_clips(cigar: &[(u8, usize)]) -> (usize, usize) {
    let lead = cigar.iter().find(|&&(op, _)| op != 5).map_or(0, |&(op, l)| if op == 4 { l } else { 0 });
    let trail = cigar.iter().rev().find(|&&(op, _)| op != 5).map_or(0, |&(op, l)| if op == 4 { l } else { 0 });
    (lead, trail)
}

/// Candidate segment indices of read base 0 from the BAM placement, via the segment anchors:
/// every anchor covering `ref_start` (query base 0 at `idx + (ref_start - ref_pos) - lead`) and
/// every anchor covering `ref_end - 1` (the last aligned query base `n-1-trail` at
/// `idx + (ref_end - 1 - ref_pos)`). Overlapping anchors (TSD: the duplicated stretch appears in
/// both genome parts of a segment) give several candidates; none -> unbanded.
fn diagonals(seg: &Segment, read: &ReadInput, lead: usize, trail: usize, out: &mut Vec<i64>) {
    out.clear();
    let covers = |a: &crate::types::Anchor, p: i64| a.ref_pos <= p && p < a.ref_pos + a.len as i64;
    for a in seg.anchors.iter().filter(|a| covers(a, read.ref_start)) {
        out.push(a.idx as i64 + (read.ref_start - a.ref_pos) - lead as i64);
    }
    let last = read.ref_end - 1;
    let q_last = read.seq.len() as i64 - 1 - trail as i64;
    for a in seg.anchors.iter().filter(|a| covers(a, last)) {
        let d = a.idx as i64 + (last - a.ref_pos) - q_last;
        if !out.contains(&d) {
            out.push(d);
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::align;
    use crate::types::{Anchor, Locus, GENOME_Q};
    use std::time::Instant;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
        fn seq(&mut self, n: usize) -> Vec<u8> {
            (0..n).map(|_| b"ACGT"[self.below(4) as usize]).collect()
        }
    }

    const F: usize = 300;
    const INS_Q: u8 = 35;

    fn other(b: u8) -> u8 {
        if b == b'A' { b'C' } else { b'A' }
    }

    fn ln_match(q: u8) -> f64 {
        ln_match2(q, GENOME_Q)
    }

    fn ln_match2(q: u8, sq: u8) -> f64 {
        let e = 10f64.powf(-(q as f64) / 10.0) + 10f64.powf(-(sq as f64) / 10.0);
        (1.0 - e).ln()
    }

    fn genome_seg(g: &[u8], lo: usize, hi: usize, label: &'static str, juncs: Vec<usize>) -> Segment {
        Segment {
            hyp: Hyp::Ref,
            label,
            seq: g[lo..hi].to_vec(),
            qual: vec![GENOME_Q; hi - lo],
            anchors: vec![Anchor { ref_pos: lo as i64, idx: 0, len: hi - lo }],
            junction_cols: juncs,
        }
    }

    fn alt_r(g: &[u8], r: usize, ins_r: &[u8]) -> Segment {
        let mut seq = g[r - F..r].to_vec();
        seq.extend_from_slice(ins_r);
        let mut qual = vec![GENOME_Q; F];
        qual.extend(std::iter::repeat_n(INS_Q, ins_r.len()));
        Segment {
            hyp: Hyp::Alt,
            label: "ALT_R",
            seq,
            qual,
            anchors: vec![Anchor { ref_pos: (r - F) as i64, idx: 0, len: F }],
            junction_cols: vec![F],
        }
    }

    fn alt_l(g: &[u8], l: usize, ins_l: &[u8]) -> Segment {
        let mut seq = ins_l.to_vec();
        seq.extend_from_slice(&g[l..l + F]);
        let mut qual = vec![INS_Q; ins_l.len()];
        qual.extend(std::iter::repeat_n(GENOME_Q, F));
        Segment {
            hyp: Hyp::Alt,
            label: "ALT_L",
            seq,
            qual,
            anchors: vec![Anchor { ref_pos: l as i64, idx: ins_l.len(), len: F }],
            junction_cols: vec![ins_l.len()],
        }
    }

    fn model(l: usize, r: usize, segments: Vec<Segment>, frac: f64, alt_full: bool) -> LocusModel {
        let locus = Locus {
            name: format!("chr1:{l}-{r}"),
            chr: "chr1".into(),
            left_pos: l as i64,
            right_pos: r as i64,
            left_open: false,
            right_open: false,
        };
        LocusModel {
            kind: locus.kind(),
            locus,
            windows: vec![(l.min(r) as i64, l.max(r) as i64)],
            breakpoints: vec![l as i64, r as i64],
            segments,
            alt_ref_junction_fraction: frac,
            alt_full,
            error: None,
        }
    }

    struct TRead {
        seq: Vec<u8>,
        qual: Vec<u8>,
        cigar: Vec<(u8, usize)>,
        ref_start: i64,
        ref_end: i64,
    }

    impl TRead {
        fn new(seq: Vec<u8>, cigar: Vec<(u8, usize)>, ref_start: usize) -> TRead {
            let rlen: usize = cigar.iter().filter(|(op, _)| matches!(op, 0 | 2 | 3 | 7 | 8)).map(|&(_, l)| l).sum();
            TRead { qual: vec![30; seq.len()], seq, cigar, ref_start: ref_start as i64, ref_end: (ref_start + rlen) as i64 }
        }
        fn input(&self) -> ReadInput<'_> {
            ReadInput {
                seq: &self.seq,
                qual: &self.qual,
                cigar: &self.cigar,
                ref_start: self.ref_start,
                ref_end: self.ref_end,
                reverse: false,
            }
        }
    }

    /// TSD locus: L = 1500, R = 1510; INS = 200 random + A×20; insR = INS[..120], insL = INS[80..].
    struct Tsd {
        g: Vec<u8>,
        ins: Vec<u8>,
        l: usize,
        r: usize,
        m: LocusModel,
    }

    fn tsd() -> Tsd {
        let mut rng = Lcg(42);
        let mut g = rng.seq(3000);
        let (l, r) = (1500usize, 1510usize);
        g[l] = b'C'; // genome after the poly-A is not an A
        g[l - 1] = b'C'; // genome before L differs from the insertion's last base (A)
        let mut ins = rng.seq(200);
        ins.extend(std::iter::repeat_n(b'A', 20));
        ins[0] = other(g[r]); // first inserted base differs from the reference after R
        ins[199] = b'G';
        let ins_r = ins[..120].to_vec();
        let ins_l = ins[80..].to_vec();
        let segs = vec![
            genome_seg(&g, l - F, r + F, "REF", vec![F, F + (r - l)]),
            alt_r(&g, r, &ins_r),
            alt_l(&g, l, &ins_l),
        ];
        let m = model(l, r, segs, 0.0, false);
        Tsd { g, ins, l, r, m }
    }

    impl Tsd {
        /// alt haplotype H = G[..R) ++ INS ++ G[L..)
        fn hap(&self) -> Vec<u8> {
            let mut h = self.g[..self.r].to_vec();
            h.extend_from_slice(&self.ins);
            h.extend_from_slice(&self.g[self.l..]);
            h
        }
        /// right-junction read with `k` inserted bases (aligned flank then soft clip)
        fn junction_r(&self, k: usize) -> TRead {
            let h = self.hap();
            let a = 151 - k;
            TRead::new(h[self.r - a..self.r + k].to_vec(), vec![(0, a), (4, k)], self.r - a)
        }
    }

    #[test]
    fn a_perfect_ref_read_is_ref() {
        let t = tsd();
        let cfg = Config::default();
        let read = TRead::new(t.g[t.l - 40..t.l + 111].to_vec(), vec![(0, 151)], t.l - 40);
        let o = score_read(&t.m, &read.input(), &cfg);
        // best alt = ALT_L: 111 aligned, the 40 bases before L clipped
        let expect = -40.0 * (ln_match(30) - 0.25f64.ln());
        assert_eq!(o.class, ReadClass::Ref, "{o:?}");
        assert!((o.llr() - expect).abs() < 0.05, "llr {} expected {}", o.llr(), expect);
        assert!((o.ll_ref - 151.0 * ln_match(30)).abs() < 0.01);
        assert_eq!(o.explained_frac, 1.0);
        assert!(o.crosses_junction);
    }

    #[test]
    fn b_junction_read_30bp_insert_is_alt() {
        let t = tsd();
        let cfg = Config::default();
        let o = score_read(&t.m, &t.junction_r(30).input(), &cfg);
        assert_eq!(o.class, ReadClass::Alt, "{o:?}");
        assert!(o.llr() > 35.0, "{o:?}");
        assert!(o.crosses_junction);
        assert_eq!(o.explained_frac, 1.0);
    }

    #[test]
    fn c_junction_read_3bp_insert_is_uninformative() {
        let t = tsd();
        let cfg = Config::default();
        let o = score_read(&t.m, &t.junction_r(3).input(), &cfg);
        assert_eq!(o.class, ReadClass::Uninformative, "{o:?}");
        assert!(o.llr() > 0.0 && o.llr() < cfg.llr_informative, "{o:?}");
    }

    #[test]
    fn d_polya_length_jitter_is_cheap() {
        let t = tsd();
        let cfg = Config::default();
        let ext = -cfg.gap_ext_phred / 10.0 * std::f64::consts::LN_10;
        let hp_open = -cfg.homopolymer_gap_open_phred / 10.0 * std::f64::consts::LN_10;
        let open = -cfg.gap_open_phred / 10.0 * std::f64::consts::LN_10;
        let read_with = |k: usize| {
            let mut s = t.ins[140..200].to_vec();
            s.extend(std::iter::repeat_n(b'A', k));
            let g = 151 - s.len();
            s.extend_from_slice(&t.g[t.l..t.l + g]);
            TRead::new(s, vec![(4, 60 + k), (0, g)], t.l)
        };
        // 80 insertion-consensus bases (q INS_Q) + 71 genome bases
        let perfect = 80.0 * ln_match2(30, INS_Q) + 71.0 * ln_match(30);
        let exact = score_read(&t.m, &read_with(20).input(), &cfg);
        assert_eq!(exact.class, ReadClass::Alt);
        assert!((exact.ll_alt - perfect).abs() < 0.01, "{exact:?}");
        for (k, gap) in [(18usize, 2usize), (23, 3)] {
            let o = score_read(&t.m, &read_with(k).input(), &cfg);
            assert_eq!(o.class, ReadClass::Alt, "k={k} {o:?}");
            assert!(o.llr() > 50.0, "k={k} {o:?}");
            let penalty = perfect - o.ll_alt;
            let hp_cost = -(hp_open + (gap - 1) as f64 * ext);
            let plain_cost = -(open + (gap - 1) as f64 * ext);
            assert!((penalty - hp_cost).abs() < 0.05, "k={k} penalty {penalty} expected {hp_cost}");
            assert!(penalty < plain_cost - 3.0);
        }
    }

    #[test]
    fn e_one_mismatch_near_junction_still_alt() {
        let t = tsd();
        let cfg = Config::default();
        for pos in [118usize, 123] {
            let mut read = t.junction_r(30);
            read.seq[pos] = other(read.seq[pos]);
            let o = score_read(&t.m, &read.input(), &cfg);
            assert_eq!(o.class, ReadClass::Alt, "pos {pos} {o:?}");
            assert!(o.llr() > 20.0, "pos {pos} {o:?}");
        }
    }

    #[test]
    fn f_unrelated_read_is_unexplained() {
        let t = tsd();
        let cfg = Config::default();
        let mut rng = Lcg(999);
        let read = TRead::new(rng.seq(151), vec![(0, 151)], t.l - 40);
        let o = score_read(&t.m, &read.input(), &cfg);
        assert_eq!(o.class, ReadClass::Unexplained, "{o:?}");
        assert!(o.explained_frac < 0.2);
    }

    #[test]
    fn g_short_insertion_alt_full() {
        let mut rng = Lcg(5);
        let mut g = rng.seq(3000);
        let (l, r) = (1500usize, 1510usize);
        let mut ins = rng.seq(40);
        ins.extend(std::iter::repeat_n(b'A', 20));
        ins[0] = other(g[r]);
        g[l - 1] = b'C';
        g[l] = b'C';
        let mut full = g[r - F..r].to_vec();
        full.extend_from_slice(&ins);
        full.extend_from_slice(&g[l..l + F]);
        let mut qual = vec![GENOME_Q; F];
        qual.extend(std::iter::repeat_n(INS_Q, ins.len()));
        qual.extend(std::iter::repeat_n(GENOME_Q, F));
        let alt = Segment {
            hyp: Hyp::Alt,
            label: "ALT_FULL",
            seq: full,
            qual,
            anchors: vec![
                Anchor { ref_pos: (r - F) as i64, idx: 0, len: F },
                Anchor { ref_pos: l as i64, idx: F + ins.len(), len: F },
            ],
            junction_cols: vec![F, F + ins.len()],
        };
        let m = model(l, r, vec![genome_seg(&g, l - F, r + F, "REF", vec![F, F + 10]), alt], 0.0, true);
        let cfg = Config::default();
        let mut s = g[r - 45..r].to_vec();
        s.extend_from_slice(&ins);
        s.extend_from_slice(&g[l..l + 46]);
        assert_eq!(s.len(), 151);
        let read = TRead::new(s, vec![(0, 45), (4, 106)], r - 45);
        let o = score_read(&m, &read.input(), &cfg);
        assert_eq!(o.class, ReadClass::Alt, "{o:?}");
        assert_eq!(o.explained_frac, 1.0);
        assert!(o.crosses_junction);
        assert!(o.llr() > 100.0, "{o:?}");
    }

    #[test]
    fn h_far_duplication_retained_ref_junction() {
        let mut rng = Lcg(77);
        let mut g = rng.seq(3000);
        let (l, r) = (1000usize, 2000usize);
        let mut ins = rng.seq(200);
        ins.extend(std::iter::repeat_n(b'A', 20));
        ins[0] = other(g[r]);
        g[l - 1] = b'C';
        let segs = vec![
            genome_seg(&g, l - F, l + F, "REF_L", vec![F]),
            genome_seg(&g, r - F, r + F, "REF_R", vec![F]),
            alt_r(&g, r, &ins[..120]),
            alt_l(&g, l, &ins[100..]),
        ];
        let m = model(l, r, segs, 0.5, false);
        assert_eq!(m.kind, crate::types::LocusKind::FarDuplication);
        let cfg = Config::default();
        let ref_read = TRead::new(g[r - 75..r + 76].to_vec(), vec![(0, 151)], r - 75);
        let o = score_read(&m, &ref_read.input(), &cfg);
        assert!((o.llr() - 0.5f64.ln()).abs() < 1e-6, "{o:?}");
        assert_eq!(o.class, ReadClass::Uninformative);
        let mut s = g[r - 121..r].to_vec();
        s.extend_from_slice(&ins[..30]);
        let alt_read = TRead::new(s, vec![(0, 121), (4, 30)], r - 121);
        let o = score_read(&m, &alt_read.input(), &cfg);
        assert_eq!(o.class, ReadClass::Alt, "{o:?}");
        assert!(o.llr() > 30.0);
    }

    #[test]
    fn missing_side_does_not_panic() {
        let t = tsd();
        let cfg = Config::default();
        let mut m = t.m.clone();
        m.segments.retain(|s| s.hyp == Hyp::Ref);
        let o = score_read(&m, &t.junction_r(30).input(), &cfg);
        assert_eq!(o.class, ReadClass::Uninformative);
        assert_eq!(o.ll_alt, f64::NEG_INFINITY);
        let mut m = t.m.clone();
        m.segments.retain(|s| s.hyp == Hyp::Alt);
        m.alt_ref_junction_fraction = 0.5;
        let o = score_read(&m, &t.junction_r(30).input(), &cfg);
        assert_eq!(o.class, ReadClass::Uninformative);
        assert_eq!(o.ll_ref, f64::NEG_INFINITY);
        m.segments.clear();
        let o = score_read(&m, &t.junction_r(30).input(), &cfg);
        assert_eq!(o.class, ReadClass::Uninformative);
    }

    /// Read drawn from `seg` starting at `p` with sparse small indels and mismatches.
    fn mutated_read(rng: &mut Lcg, seg: &[u8], p: usize, n: usize) -> (Vec<u8>, Vec<u8>) {
        let mut s = Vec::with_capacity(n);
        let mut j = p;
        while s.len() < n && j < seg.len() {
            let u = rng.below(1000);
            if u < 10 {
                for _ in 0..1 + rng.below(3) {
                    s.push(b"ACGT"[rng.below(4) as usize]);
                }
            } else if u < 20 {
                j += 1 + rng.below(3) as usize;
            } else if u < 40 {
                s.push(other(seg[j]));
                j += 1;
            } else {
                s.push(seg[j]);
                j += 1;
            }
        }
        s.truncate(n);
        let q = (0..s.len()).map(|_| 10 + rng.below(31) as u8).collect();
        (s, q)
    }

    #[test]
    fn i_banded_equals_unbanded() {
        let mut cfg = Config::default();
        cfg.realign_fallback_full = false;
        let mut rng = Lcg(2024);
        for case in 0..200 {
            let seg = rng.seq(400);
            let mut sq = vec![GENOME_Q; 400];
            for q in sq[150..250].iter_mut() {
                *q = 20 + rng.below(20) as u8;
            }
            let p = 10 + rng.below(220) as usize;
            let (s, q) = mutated_read(&mut rng, &seg, p, 151);
            let diag = p as i64 + rng.below(11) as i64 - 5;
            let b = align(&s, &q, &seg, &sq, Some(diag), &cfg);
            let u = align(&s, &q, &seg, &sq, None, &cfg);
            assert!((b.ll - u.ll).abs() < 1e-3, "case {case}: banded {b:?} unbanded {u:?}");
            assert_eq!((b.seg_start, b.seg_end), (u.seg_start, u.seg_end), "case {case}");
            assert_eq!(b.aligned_bases, u.aligned_bases, "case {case}");
        }
    }

    #[test]
    fn j_wrong_diagonal_falls_back() {
        let mut rng = Lcg(31337);
        let seg = rng.seq(400);
        let sq = vec![GENOME_Q; 400];
        let read = seg[60..211].to_vec();
        let q = vec![30u8; 151];
        let mut cfg = Config::default();
        let full = align(&read, &q, &seg, &sq, None, &cfg);
        let fb = align(&read, &q, &seg, &sq, Some(160), &cfg);
        assert!((fb.ll - full.ll).abs() < 1e-3);
        assert_eq!((fb.seg_start, fb.seg_end), (60, 211));
        cfg.realign_fallback_full = false;
        let nofb = align(&read, &q, &seg, &sq, Some(160), &cfg);
        assert!(nofb.ll < full.ll - 50.0, "{nofb:?}");
    }

    #[test]
    fn diagonal_from_anchors() {
        let t = tsd();
        // ALT_L junction read: soft clip 50 then aligned from L
        let alt_l_seg = &t.m.segments[2];
        let read = TRead::new(vec![b'A'; 151], vec![(4, 50), (0, 101)], t.l);
        let (lead, trail) = soft_clips(&read.cigar);
        let mut d = Vec::new();
        diagonals(alt_l_seg, &read.input(), lead, trail, &mut d);
        assert_eq!(d, vec![140 - 50]);
        // ALT_R read whose start is left of the anchor: placed via its aligned end
        let alt_r_seg = &t.m.segments[1];
        let read = TRead::new(vec![b'A'; 151], vec![(5, 7), (0, 141), (4, 10)], t.r - F - 20);
        let (lead, trail) = soft_clips(&read.cigar);
        assert_eq!((lead, trail), (0, 10));
        // last aligned base at ref r-F+120 -> idx 120, query index 140 -> diag -20
        diagonals(alt_r_seg, &read.input(), lead, trail, &mut d);
        assert_eq!(d, vec![-20]);
        let far = TRead::new(vec![b'A'; 151], vec![(0, 151)], 10);
        diagonals(alt_r_seg, &far.input(), 0, 0, &mut d);
        assert!(d.is_empty());
        // overlapping anchors (TSD duplicated in both genome parts): both diagonals offered
        let mut seg = alt_r_seg.clone();
        seg.anchors.push(Anchor { ref_pos: (t.r - F) as i64, idx: 200, len: 50 });
        let read = TRead::new(vec![b'A'; 151], vec![(0, 151)], t.r - F + 10);
        diagonals(&seg, &read.input(), 0, 0, &mut d);
        assert_eq!(d, vec![10, 210]);
    }

    /// `cargo test --release -- --ignored --nocapture bench_cells`
    #[test]
    #[ignore]
    fn bench_cells_per_second() {
        let mut cfg = Config::default();
        cfg.realign_fallback_full = false;
        let mut rng = Lcg(1);
        let seg = rng.seq(400);
        let sq = vec![GENOME_Q; 400];
        let reads: Vec<(Vec<u8>, Vec<u8>, usize)> = (0..64)
            .map(|_| {
                let p = rng.below(240) as usize;
                let (s, q) = mutated_read(&mut rng, &seg, p, 151);
                (s, q, p)
            })
            .collect();
        let w = 2 * cfg.band_halfwidth + 1;
        for (name, banded, iters) in [("banded", true, 40_000usize), ("unbanded", false, 8_000)] {
            let t0 = Instant::now();
            let mut acc = 0.0;
            let mut cells = 0usize;
            for it in 0..iters {
                let (s, q, p) = &reads[it % reads.len()];
                let d = if banded { Some(*p as i64) } else { None };
                acc += align(s, q, &seg, &sq, d, &cfg).ll;
                cells += if banded { s.len() * w } else { s.len() * seg.len() };
            }
            let dt = t0.elapsed().as_secs_f64();
            println!(
                "{name:>8} 151x400: {:.1} M cells/s, {:.2} us/alignment (checksum {acc:.1})",
                cells as f64 / dt / 1e6,
                dt / iters as f64 * 1e6
            );
        }
    }
}
