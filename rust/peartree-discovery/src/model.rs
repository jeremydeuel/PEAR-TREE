//! Port of src/breakpoint.py (Breakpoint + Breakpoint.join).

use crate::config::*;
use crate::evidence::{frag_hash, select_lowest, ClipLite, EvExtra};
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
    /// Feature B: mate landing sites (ref_id, 0-based pos), filled in find_mates when
    /// `splice_hallmark`/`discordant_anchor` is on. Used by the splice/pseudogene
    /// annotation; does not affect the main output.
    pub mate_dests: Vec<(Option<usize>, i64)>,
    pub n_reads: usize,
    /// distinct fragments (qname hashes) among the reads supporting the consensus
    /// position; 1 for a per-read or rescued breakpoint. Used by the one-sided floor.
    pub n_frags: usize,
    /// TPRT sidecar payload (`evidence_sidecar`); None otherwise (8 bytes).
    pub ev: Option<Box<EvExtra>>,
    /// raw SAM flag of the source read (0 for a synthesised consensus breakpoint); the
    /// sidecar uses it to re-find the record in the mate pass.
    pub flag: u16,
    /// raw RNEXT/PNEXT (0-based) of the source read: where its mate record sits
    /// (sidecar indexed fetch); -1 when unset or for a synthesised breakpoint.
    pub mref: i32,
    pub mpos: i64,
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
            mate_dests: Vec::new(),
            n_reads: 1,
            n_frags: 1,
            ev: None,
            flag: 0,
            mref: -1,
            mpos: -1,
        }
    }
}

/// SENS-8: a clipped consensus is a "pure" poly-A/T terminus when >=90% of its
/// bases are A (or >=90% T). Used to relax the clip-length floor for genuine
/// poly-A tails only.
fn is_pure_polya(qs: &QualitySeq) -> bool {
    let n = qs.len();
    if n == 0 {
        return false;
    }
    let a = qs.seq.iter().filter(|&&b| b == b'A' || b == b'a').count();
    let t = qs.seq.iter().filter(|&&b| b == b'T' || b == b't').count();
    (a.max(t) as f64 / n as f64) >= 0.9
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

/// Fragment id of a per-read breakpoint (qname hash). Mates and primary +
/// supplementary records of one template share it.
#[inline]
fn bp_frag(bp: &Breakpoint) -> u64 {
    frag_hash(bp.query_name.as_deref().unwrap_or("").as_bytes())
}

/// Number of distinct fragments in `frags`.
pub fn count_fragments(frags: &[u64]) -> usize {
    let mut v = frags.to_vec();
    v.sort_unstable();
    v.dedup();
    v.len()
}

/// The sidecar identities (qname hash, flag, breakpoint coordinate) of every read of a
/// cluster, capped deterministically at `max_evidence_reads_per_breakpoint` (lowest
/// frag hash, then flag). None when the sidecar is off.
fn merge_ev(breakpoints: &[Breakpoint], cfg: &DiscoveryConfig) -> Option<Box<EvExtra>> {
    if !cfg.evidence_sidecar {
        return None;
    }
    let lite: Vec<ClipLite> = breakpoints
        .iter()
        .map(|bp| ClipLite { frag: bp_frag(bp), flag: bp.flag, pos: bp.breakpoint, mref: bp.mref, mpos: bp.mpos })
        .collect();
    let clip_lite = select_lowest(lite, cfg.max_evidence_reads_per_breakpoint, |r| (r.frag, r.flag));
    Some(Box::new(EvExtra { clip_lite, ..Default::default() }))
}

/// Port of Breakpoint.join. Consumes the group of breakpoints in a <6bp window
/// and returns a single consensus breakpoint, or None if filtered out.
pub fn join(mut breakpoints: Vec<Breakpoint>, cfg: &DiscoveryConfig, evidence_floor: usize, stats: &mut Stats) -> Option<Breakpoint> {
    // TPRT fragment mode: the floor counts distinct fragments, and a single-fragment
    // cluster takes the normal consensus path when the floor admits it.
    let frag_mode = cfg.min_evidence_fragments_per_sample.is_some();
    if breakpoints.len() < 2 && !(frag_mode && evidence_floor <= 1) {
        let solo_ev = merge_ev(&breakpoints, cfg);
        let bp = breakpoints.pop().unwrap();
        let side = bp.side;
        return match rescue_polya(bp) {
            Rescue::Rescued(mut b) => {
                stats.side_mut(side).rescued_pa += 1;
                b.ev = solo_ev;
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
    let mut frags: Vec<u64> = Vec::new();
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
            frags.push(bp_frag(bp));
        }
    }
    if bps.is_empty() {
        stats.side_mut(side).too_few_after_filter += 1;
        return None;
    }

    let (best_bp, exact_n) = most_common_first(&bps);
    // SENS-1/OBS-3: count support within +/- evidence_window of the mode, not only
    // at the exact modal position. 0 = exact (legacy, byte-identical).
    let n = if frag_mode {
        // distinct fragments among the supporting reads (window or exact mode)
        let support: Vec<u64> = bps
            .iter()
            .zip(&frags)
            .filter(|(&b, _)| if cfg.evidence_window > 0 { (b - best_bp).abs() <= cfg.evidence_window } else { b == best_bp })
            .map(|(_, &f)| f)
            .collect();
        count_fragments(&support)
    } else if cfg.evidence_window > 0 {
        bps.iter().filter(|&&b| (b - best_bp).abs() <= cfg.evidence_window).count()
    } else {
        exact_n
    };

    if n < evidence_floor {
        // the sidecar keeps every read of the cluster on the rescued breakpoint
        let group_ev = merge_ev(&breakpoints, cfg);
        // try polyA rescue on the individual breakpoints; first hit wins
        for bp in breakpoints.into_iter() {
            let side_bp = bp.side;
            let is_a = side_bp == CLIP_LEFT && bp.clipped.pyslice(Some(-8), None).eq_bytes(b"AAAAAAAA");
            let is_t = side_bp == CLIP_RIGHT && bp.clipped.pyslice(None, Some(8)).eq_bytes(b"TTTTTTTT");
            if is_a || is_t {
                return match rescue_polya(bp) {
                    Rescue::Rescued(mut b) => {
                        stats.side_mut(side).rescued_pa += 1;
                        b.ev = group_ev;
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
    // reads actually contributing to the consensus (LEFT delta 0 is pushed twice into
    // `clipped`, so `clipped.len()` over-counts; used only in fragment mode)
    let mut n_used: usize = 0;
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
        let before = clipped.len();
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
        if clipped.len() > before {
            n_used += 1;
        }
    }

    let clipped_cons = find_consensus(&clipped, cfg.consensus_tolerant);
    let unclipped_cons = find_consensus(&unclipped, cfg.consensus_tolerant);
    // SENS-8: keep the 12 bp floor, but allow shorter clips when the clipped
    // consensus is a pure poly-A/T terminus (a genuine tail). Off = legacy floor.
    let clip_floor = if cfg.short_polya_clip && is_pure_polya(&clipped_cons) {
        cfg.short_polya_min_clip
    } else {
        MIN_CLIP_LEN
    };
    if clipped_cons.len() <= clip_floor {
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
    b.n_frags = {
        let support: Vec<u64> = bps
            .iter()
            .zip(&frags)
            .filter(|(&x, _)| if cfg.evidence_window > 0 { (x - best_bp).abs() <= cfg.evidence_window } else { x == best_bp })
            .map(|(_, &f)| f)
            .collect();
        count_fragments(&support)
    };
    // legacy n_reads double-counts LEFT reads at delta 0; it only feeds the Feature A
    // rescue (off by default), so the legacy value is kept unless in fragment mode.
    b.n_reads = if frag_mode { n_used } else { clipped.len() };
    b.ev = merge_ev(&breakpoints, cfg);
    stats.side_mut(side).passed += 1;
    Some(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// deterministic pseudo-random DNA (LCG), so consensus filters are not tripped
    fn dna(seed: u64, n: usize) -> Vec<u8> {
        let mut x = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    fn bp(side: i32, pos: i64, qname: &str, read1: bool, forward: bool) -> Breakpoint {
        let clip = dna(side as u64, 30);
        let flank = dna(10 + side as u64, 80);
        Breakpoint::new(
            side,
            "chr1".into(),
            pos,
            Some(qname.into()),
            QualitySeq::new(clip, vec![30; 30]),
            QualitySeq::new(flank, vec![30; 80]),
            Some(read1),
            Some(forward),
            false,
            60,
        )
    }

    /// RIGHT-clipped read at position 1000 sharing one junction (same clip/flank).
    fn right_bp(qname: &str, read1: bool) -> Breakpoint {
        bp(CLIP_RIGHT, 1000, qname, read1, true)
    }

    fn cfg(frag: Option<usize>) -> DiscoveryConfig {
        DiscoveryConfig { min_evidence_fragments_per_sample: frag, ..DiscoveryConfig::default() }
    }

    #[test]
    fn count_fragments_dedups() {
        assert_eq!(count_fragments(&[5, 3, 5, 5, 9]), 3);
        assert_eq!(count_fragments(&[]), 0);
    }

    #[test]
    fn both_mates_of_one_fragment_count_once() {
        let mut st = Stats::default();
        // legacy read floor: 2 reads (the two mates) pass
        let g = vec![right_bp("fragA", true), right_bp("fragA", false)];
        assert!(join(g.clone(), &cfg(None), 2, &mut st).is_some());
        // fragment floor 2: one template -> rejected
        assert!(join(g.clone(), &cfg(Some(2)), 2, &mut st).is_none());
        // fragment floor 1: passes, and n_reads counts the 2 distinct reads
        let b = join(g, &cfg(Some(1)), 1, &mut st).unwrap();
        assert_eq!(b.n_reads, 2);
    }

    #[test]
    fn primary_and_supplementary_count_once() {
        let mut st = Stats::default();
        // read1 primary + read1 supplementary (same qname) + an independent fragment
        let g = vec![right_bp("fragA", true), right_bp("fragA", true), right_bp("fragB", true)];
        assert!(join(g.clone(), &cfg(Some(2)), 2, &mut st).is_some());
        assert!(join(g.clone(), &cfg(Some(3)), 3, &mut st).is_none());
        // legacy counts 3 reads
        assert!(join(g, &cfg(None), 3, &mut st).is_some());
    }

    #[test]
    fn single_fragment_cluster_passes_only_in_fragment_mode() {
        let mut st = Stats::default();
        let g = vec![right_bp("solo", true)];
        assert!(join(g.clone(), &cfg(None), 2, &mut st).is_none()); // legacy: too_few
        let b = join(g.clone(), &cfg(Some(1)), 1, &mut st).expect("floor 1 admits one fragment");
        assert_eq!(b.n_reads, 1);
        assert!(join(g, &cfg(Some(2)), 2, &mut st).is_none());
    }

    #[test]
    fn left_n_reads_double_count_fixed_only_in_fragment_mode() {
        let mut st = Stats::default();
        let g = vec![bp(CLIP_LEFT, 500, "a", true, false), bp(CLIP_LEFT, 500, "b", true, false)];
        assert_eq!(join(g.clone(), &cfg(None), 2, &mut st).unwrap().n_reads, 4); // legacy x2
        assert_eq!(join(g, &cfg(Some(2)), 2, &mut st).unwrap().n_reads, 2);
    }
}
