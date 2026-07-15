//! Port of src/breakpoint.py (Breakpoint + Breakpoint.join).

use crate::config::*;
use crate::filters::{clean_clipped_seq, find_consensus};
use crate::qseq::QualitySeq;
use crate::stats::Stats;

#[derive(Clone)]
pub struct Breakpoint {
    pub exclude: bool,
    pub is_forward: Option<bool>,
    pub is_read1: Option<bool>,
    pub query_name: Option<String>,
    pub breakpoint: i64,
    pub clipped: QualitySeq,
    pub unclipped: QualitySeq,
    pub reference_name: String,
    pub side: i32,
    pub has_mate: bool,
    pub bp_precise: bool,
    /// mapping quality of the source read (SENS-2 low-MAPQ-clip-ratio guard); 0 for
    /// a synthesised consensus breakpoint.
    pub mapq: u8,
    /// (mate_is_read1, qname)
    pub mates: Vec<(bool, String)>,
    pub mate_seqs: Vec<QualitySeq>,
    pub n_reads: usize,
}

impl Breakpoint {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        side: i32,
        reference_name: String,
        breakpoint: i64,
        query_name: Option<String>,
        clipped: QualitySeq,
        unclipped: QualitySeq,
        is_read1: Option<bool>,
        is_forward: Option<bool>,
        exclude: bool,
        mapq: u8,
    ) -> Self {
        Breakpoint {
            exclude,
            is_forward,
            is_read1,
            query_name,
            breakpoint,
            clipped,
            unclipped,
            reference_name,
            side,
            has_mate: false,
            bp_precise: true,
            mapq,
            mates: Vec::new(),
            mate_seqs: Vec::new(),
            n_reads: 1,
        }
    }
}

/// Counter.most_common()[0]: most frequent value, ties broken by first
/// insertion order. Returns (value, count).
fn most_common_first(values: &[i64]) -> (i64, usize) {
    // preserve first-seen order
    let mut order: Vec<i64> = Vec::new();
    let mut counts: Vec<usize> = Vec::new();
    for &v in values {
        if let Some(idx) = order.iter().position(|&x| x == v) {
            counts[idx] += 1;
        } else {
            order.push(v);
            counts.push(1);
        }
    }
    let mut best = 0usize;
    for i in 1..order.len() {
        if counts[i] > counts[best] {
            best = i;
        }
    }
    (order[best], counts[best])
}

/// Outcome of a polyA-tail rescue attempt. Distinguishes the three Python paths so
/// OBS-1 can count them exactly: a rescued breakpoint (`rescued_pA`), a non-polyA
/// solo (`too_few`), and a polyA whose cleaned clip fell below min_clip_len — which
/// Python returns as None **without** incrementing any counter.
enum Rescue {
    Rescued(Breakpoint),
    TooShort,
    NotPolyA,
}

/// polyA-tail rescue for a single (solo or low-support) breakpoint. Mutates and
/// returns the breakpoint if it qualifies. Mirrors the two rescue blocks in
/// Breakpoint.join. Counter bookkeeping is left to the caller (see `join`).
fn rescue_polya(mut bp: Breakpoint) -> Rescue {
    if bp.side == CLIP_LEFT && bp.clipped.pyslice(Some(-8), None).eq_bytes(b"AAAAAAAA") {
        bp.clipped = clean_clipped_seq(&bp.clipped.revcomp());
        if bp.clipped.len() < MIN_CLIP_LEN {
            return Rescue::TooShort;
        }
        if bp.is_forward == Some(false) {
            bp.has_mate = true;
            bp.mates = vec![(!(bp.is_read1.unwrap_or(false)), bp.query_name.clone().unwrap_or_default())];
        }
        Rescue::Rescued(bp)
    } else if bp.side == CLIP_RIGHT && bp.clipped.pyslice(None, Some(8)).eq_bytes(b"TTTTTTTT") {
        bp.clipped = clean_clipped_seq(&bp.clipped);
        if bp.clipped.len() < MIN_CLIP_LEN {
            return Rescue::TooShort;
        }
        if bp.is_forward == Some(true) {
            bp.has_mate = true;
            bp.mates = vec![(!(bp.is_read1.unwrap_or(false)), bp.query_name.clone().unwrap_or_default())];
        }
        Rescue::Rescued(bp)
    } else {
        Rescue::NotPolyA
    }
}

/// Port of Breakpoint.join. Consumes the group of breakpoints in a <6bp window
/// and returns a single consensus breakpoint, or None if filtered out.
pub fn join(mut breakpoints: Vec<Breakpoint>, cfg: &DiscoveryConfig, evidence_floor: usize, stats: &mut Stats) -> Option<Breakpoint> {
    if breakpoints.len() < 2 {
        let bp = breakpoints.pop().unwrap();
        let side = bp.side;
        return match rescue_polya(bp) {
            Rescue::Rescued(b) => {
                stats.side_mut(side).rescued_pa += 1;
                Some(b)
            }
            // polyA candidate whose clip cleaned below min_clip_len: Python counts nothing
            Rescue::TooShort => None,
            Rescue::NotPolyA => {
                stats.side_mut(side).too_few += 1;
                None
            }
        };
    }

    // if any breakpoint has the exclude flag, poison the whole group
    if breakpoints.iter().any(|b| b.exclude) {
        stats.side_mut(breakpoints[0].side).excluded += 1;
        return None;
    }

    let side = breakpoints[0].side;
    let reference_name = breakpoints[0].reference_name.clone();

    // SENS-2 guard: with a lowered min_mapq, drop a locus dominated by low-MAPQ
    // clipped reads (a low-quality pileup, not a real junction). Ships with the
    // relaxation so recall gains near repeats do not come with low-MAPQ artefacts.
    if let Some(max_ratio) = cfg.max_lowq_clip_ratio {
        let total = breakpoints.len();
        let lowq = breakpoints.iter().filter(|b| b.mapq < cfg.lowq_mapq_threshold).count();
        if total > 0 && (lowq as f64 / total as f64) > max_ratio {
            stats.side_mut(side).excluded += 1;
            return None;
        }
    }

    // clean/orient clipped sequences and collect qc-passing precise positions
    let mut bps: Vec<i64> = Vec::new();
    for bp in breakpoints.iter_mut() {
        if side == CLIP_LEFT {
            // left clipped is reverse complemented from now on
            bp.clipped = clean_clipped_seq(&bp.clipped.revcomp());
        } else {
            bp.clipped = clean_clipped_seq(&bp.clipped);
        }
        if side == CLIP_LEFT && bp.is_forward == Some(false) {
            bp.has_mate = true;
        }
        if side == CLIP_RIGHT && bp.is_forward == Some(true) {
            bp.has_mate = true;
        }
        if bp.bp_precise && bp.clipped.len() >= cfg.min_good_bases {
            bps.push(bp.breakpoint);
        }
    }
    if bps.is_empty() {
        stats.side_mut(side).too_few_after_filter += 1;
        return None;
    }

    let (best_bp, exact_n) = most_common_first(&bps);
    // SENS-1/OBS-3: count support within +/- evidence_window of the mode, not only
    // at the exact modal position. 0 = exact (legacy, byte-identical).
    let n = if cfg.evidence_window > 0 {
        bps.iter().filter(|&&b| (b - best_bp).abs() <= cfg.evidence_window).count()
    } else {
        exact_n
    };

    if n < evidence_floor {
        // try polyA rescue on the individual breakpoints; first hit wins
        for bp in breakpoints.into_iter() {
            let side_bp = bp.side;
            let is_a = side_bp == CLIP_LEFT && bp.clipped.pyslice(Some(-8), None).eq_bytes(b"AAAAAAAA");
            let is_t = side_bp == CLIP_RIGHT && bp.clipped.pyslice(None, Some(8)).eq_bytes(b"TTTTTTTT");
            if is_a || is_t {
                return match rescue_polya(bp) {
                    Rescue::Rescued(b) => {
                        stats.side_mut(side).rescued_pa += 1;
                        Some(b)
                    }
                    // guard above guarantees polyA, so NotPolyA is unreachable here
                    Rescue::TooShort | Rescue::NotPolyA => None,
                };
            }
        }
        stats.side_mut(side).too_few_after_filter += 1;
        return None;
    }

    // build clipped/unclipped consensus inputs
    let mut mates: Vec<(bool, String)> = Vec::new();
    let mut clipped: Vec<QualitySeq> = Vec::new();
    let mut unclipped: Vec<QualitySeq> = Vec::new();
    for bp in breakpoints.iter() {
        if bp.has_mate {
            mates.push((!(bp.is_read1.unwrap_or(false)), bp.query_name.clone().unwrap_or_default()));
        }
        if !bp.bp_precise {
            continue;
        }
        let delta_bp = bp.breakpoint - best_bp;
        let clen = bp.clipped.len() as i64;
        let ulen = bp.unclipped.len() as i64;
        if side == CLIP_RIGHT {
            if delta_bp == 0 {
                clipped.push(bp.clipped.clone());
                unclipped.push(bp.unclipped.revcomp());
            } else if delta_bp > 0 && delta_bp < clen {
                clipped.push(bp.unclipped.pyslice(Some(-delta_bp as isize), None).concat(&bp.clipped));
                unclipped.push(bp.unclipped.pyslice(None, Some(-delta_bp as isize)).revcomp());
            } else if delta_bp < 0 && delta_bp > -ulen {
                clipped.push(bp.clipped.pyslice(Some(-delta_bp as isize), None));
                unclipped.push(
                    bp.clipped
                        .pyslice(None, Some(-delta_bp as isize))
                        .revcomp()
                        .concat(&bp.unclipped.revcomp()),
                );
            }
        } else if side == CLIP_LEFT {
            clipped.push(bp.clipped.clone());
            unclipped.push(bp.unclipped.clone());
            if delta_bp == 0 {
                clipped.push(bp.clipped.clone());
                unclipped.push(bp.unclipped.clone());
            } else if delta_bp > 0 && delta_bp < ulen {
                clipped.push(bp.clipped.pyslice(Some(delta_bp as isize), None));
                unclipped.push(
                    bp.clipped
                        .pyslice(None, Some(delta_bp as isize))
                        .revcomp()
                        .concat(&bp.unclipped),
                );
            } else if delta_bp < 0 && delta_bp > -clen {
                clipped.push(
                    bp.unclipped
                        .pyslice(None, Some(-delta_bp as isize))
                        .revcomp()
                        .concat(&bp.clipped),
                );
                unclipped.push(bp.unclipped.pyslice(Some(-delta_bp as isize), None));
            }
        }
    }

    let clipped_cons = find_consensus(&clipped, cfg.consensus_tolerant);
    let unclipped_cons = find_consensus(&unclipped, cfg.consensus_tolerant);
    if clipped_cons.len() <= MIN_CLIP_LEN {
        stats.side_mut(side).clipped_failed += 1;
        return None;
    }
    if unclipped_cons.len() <= 40 {
        stats.side_mut(side).unclipped_failed += 1;
        return None;
    }
    // reject n-polymer at the breakpoint (n = 1..=4)
    for nn in 1..5usize {
        let base = unclipped_cons.pyslice(None, Some(nn as isize));
        let reps = 24 / nn;
        let query: Vec<u8> = base.seq.iter().cycle().take(base.seq.len() * reps).copied().collect();
        if unclipped_cons.pyslice(None, Some(24)).eq_bytes(&query) {
            stats.side_mut(side).polymer += 1;
            return None;
        }
    }

    let mut b = Breakpoint::new(
        side,
        reference_name,
        best_bp,
        None,
        clipped_cons,
        unclipped_cons,
        None,
        None,
        false,
        0, // synthesised consensus breakpoint: mapq unused
    );
    b.mates = mates;
    b.n_reads = clipped.len();
    stats.side_mut(side).passed += 1;
    Some(b)
}
