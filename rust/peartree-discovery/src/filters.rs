//! Ports of src/adapter.py (is_adapter), src/sequence_checks.py
//! (clean_clipped_seq) and src/consensus.py (find_consensus).

use crate::config::{ADAPTERS, CLIP_LEFT, MIN_ADAPTERLEN_FOR_CLIP, MIN_CLIP_LEN};
use crate::qseq::{revcomp_bytes, QualitySeq};

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// Port of adapter.is_adapter. `seq` is expected to be MIN_CLIP_LEN bases (may
/// be fewer). Operates on ASCII bytes; caller passes already-oriented sequence.
pub fn is_adapter(seq_in: &[u8]) -> bool {
    let seq: Vec<u8> = seq_in.iter().map(|b| b.to_ascii_uppercase()).collect();
    if seq.len() < 2 {
        return false;
    }
    // remove all homopolymers except polyT (which corresponds to a polyA tail)
    for &b in &[b'G', b'C', b'A'] {
        if seq.len() == MIN_CLIP_LEN && seq.iter().all(|&c| c == b) {
            return true;
        }
    }
    if seq[0] != seq[1] {
        let half = MIN_CLIP_LEN / 2;
        let dimer: Vec<u8> = seq[0..2].iter().cycle().take(half * 2).copied().collect();
        if seq.len() >= half * 2 && seq[..half * 2] == dimer[..] {
            return true;
        }
    }
    let revseq = revcomp_bytes(&seq);
    for a in ADAPTERS.iter() {
        if contains(a, &seq) || contains(a, &revseq) {
            return true;
        }
    }
    false
}

/// True if `seq` is dominated by a short-period tandem repeat (period 1..=6) — a
/// homopolymer, dimer, or satellite like (CCATT)n / (GGAAT)n. Used to reject a one-sided
/// discordant call whose "element" clip is actually a low-complexity satellite array (a
/// pericentromeric/subtelomeric mismap), not a real retrotransposon: a genuine Alu/L1/
/// SVA/ERV clip is complex and matches no short period. `min_frac` of the comparable
/// positions must repeat at some period. Sequences shorter than 10 bp are left to the
/// other gates (too short to judge).
pub fn is_low_complexity(seq: &[u8], min_frac: f64) -> bool {
    let n = seq.len();
    if n < 10 {
        return false;
    }
    let up: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    for p in 1..=6usize {
        if n <= p {
            continue;
        }
        let matches = (p..n).filter(|&i| up[i] == up[i - p]).count();
        if (matches as f64) / ((n - p) as f64) >= min_frac {
            return true;
        }
    }
    false
}

/// Mean distinct-k-mer fraction (|distinct k-mers| / |k-mer positions|) over the reads
/// long enough to judge (>= 20 bp). A satellite/low-diversity read has few distinct
/// k-mers (fraction near 0.3); complex element/genomic sequence is near 0.6-0.9. Used to
/// reject a one-sided discordant call whose mate reads are a satellite array (the flank-
/// anchored mates of a real MEI are complex, so this spares real poly-A/VNTR clips —
/// which a clip-level complexity test would wrongly kill). `None` if no read qualifies.
pub fn mean_kmer_diversity(seqs: &[QualitySeq], k: usize) -> Option<f64> {
    let mut acc = 0.0f64;
    let mut n = 0usize;
    for s in seqs {
        let len = s.seq.len();
        if len < 20 || len < k + 1 {
            continue;
        }
        let total = len - k + 1;
        let mut set: std::collections::HashSet<&[u8]> = std::collections::HashSet::with_capacity(total);
        for i in 0..total {
            set.insert(&s.seq[i..i + k]);
        }
        acc += set.len() as f64 / total as f64;
        n += 1;
    }
    if n == 0 {
        None
    } else {
        Some(acc / n as f64)
    }
}

/// How far to read into the aligned side when looking for the tandem tract at the junction.
const SLIPPAGE_MAX_LOOK: usize = 64;

/// Bases of `seq` read *outward from the junction*, uppercased, capped at `n`.
/// `from_start` reads forwards from index 0, otherwise backwards from the end.
fn junction_outward(seq: &[u8], from_start: bool, n: usize) -> Vec<u8> {
    let m = seq.len().min(n);
    (0..m)
        .map(|i| {
            let b = if from_start { seq[i] } else { seq[seq.len() - 1 - i] };
            b.to_ascii_uppercase()
        })
        .collect()
}

/// The tandem repeat terminating at the junction, given the reference bases read outward
/// from it. Returns `(period, run_len)` for the period yielding the longest terminal run;
/// the smallest period wins ties, so a homopolymer reports period 1 rather than 2. A period
/// > 1 must show at least two full copies. `None` if no periodic tract reaches the junction.
fn terminal_tandem(ref_out: &[u8], max_period: usize) -> Option<(usize, usize)> {
    let mut best: Option<(usize, usize)> = None;
    for p in 1..=max_period.min(ref_out.len()) {
        let mut t = p;
        while t < ref_out.len() && ref_out[t] == ref_out[t - p] {
            t += 1;
        }
        if p > 1 && t < p * 2 {
            continue; // fewer than two copies — not a tandem
        }
        if best.is_none_or(|(_, bl)| t > bl) {
            best = Some((p, t));
        }
    }
    best
}

/// SPEC-8: true if this breakpoint looks like **replication slippage against a tandem tract
/// already in the reference**, rather than an insertion junction.
///
/// A tandem tract in the reference — canonically the poly-A tail of a reference Alu, but
/// equally a (CA)n / (TG)n / (TAAAA)n microsatellite — varies in copy number between sample
/// and reference, because polymerase slips in short repeats. The aligner cannot place the
/// surplus copies and soft-clips them, so discovery sees a clip cluster indistinguishable
/// from a real junction. The clip is simply *more of the tract that is already there*.
///
/// The test is deliberately **paired**: the reference (aligned) side of the junction must end
/// in a tandem tract of period <= `max_period` and length >= `min_ref_run`, AND the clip must
/// continue that same tract, in phase, for >= `min_clip_frac` of its bases. A clip-only
/// complexity test would wrongly kill genuine MEIs, whose 3' clip is a real poly-A tail — the
/// difference is that a real tail's poly-A is *not* also in the reference at the junction.
/// Requiring both sides is what spares them. A real MEI landing *in* an STR is spared too:
/// its clip is element body, which does not continue the tract's phase.
///
/// `max_period = 1` is the homopolymer-only behaviour; raising it to ~6 also catches STRs.
///
/// Orientation follows `src/discovery.py::add_breakpoint`:
///   `CLIP_LEFT : [clipped][unclipped]` -> tract at the **start** of `unclipped`, clip runs
///                                        leftward from the junction
///   `CLIP_RIGHT: [unclipped][clipped]` -> tract at the **end** of `unclipped`, clip runs
///                                        rightward from the junction
pub fn is_slippage_clip(
    side: i32,
    clipped: &[u8],
    unclipped: &[u8],
    min_ref_run: usize,
    min_clip_frac: f64,
    max_period: usize,
) -> bool {
    if clipped.is_empty() || unclipped.is_empty() || max_period == 0 {
        return false;
    }
    let ref_out = junction_outward(unclipped, side == CLIP_LEFT, SLIPPAGE_MAX_LOOK);
    let clip_out = junction_outward(clipped, side != CLIP_LEFT, clipped.len());
    let (period, run) = match terminal_tandem(&ref_out, max_period) {
        Some(x) => x,
        None => return false,
    };
    if run < min_ref_run.max(period * 2) {
        return false;
    }
    let unit = &ref_out[..period];
    if !unit.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
        return false;
    }
    // Continue the tract's phase across the junction: `ref_out[k]` sits k bases into the
    // reference, so the clip base k bases out sits at reference offset -(k+1), i.e. at
    // unit index (period - (k+1) % period) % period.
    let same = clip_out
        .iter()
        .enumerate()
        .filter(|(k, &c)| c == unit[(period - ((k + 1) % period)) % period])
        .count();
    (same as f64) / (clip_out.len() as f64) >= min_clip_frac
}

/// Port of sequence_checks.clean_clipped_seq (the QualitySeq version).
pub fn clean_clipped_seq(seq: &QualitySeq) -> QualitySeq {
    if seq.len() < 2 {
        return QualitySeq::empty();
    }
    // remove poly G at tail of read
    let mut i = seq.len();
    while i > 1 && seq.base(i - 1) == b'G' {
        i -= 1;
    }
    let mut seq = if i < seq.len().saturating_sub(4) {
        seq.pyslice(None, Some(i as isize))
    } else {
        seq.clone()
    };
    // remove adapter sequence, greedily matching ADAPTERS[0] anywhere.
    let adapter0 = ADAPTERS[0];
    let mut i: usize = 0;
    let mut j: usize = 0;
    while i < seq.len() {
        if seq.base(i) == adapter0[j] {
            j += 1;
        } else if j > 0 {
            // give j=0 another chance in case the adapter starts exactly here
            j = 0;
            continue;
        }
        i += 1;
        if j > 9 {
            break;
        }
    }
    if j >= MIN_ADAPTERLEN_FOR_CLIP {
        seq = seq.pyslice(None, Some((i - j) as isize));
    }
    seq
}

/// Port of consensus.find_consensus. Left-aligned quality-weighted consensus,
/// truncated at the first ambiguous position. Tie-break order is A,T,G,C.
pub fn find_consensus(seqs: &[QualitySeq], tolerant: bool) -> QualitySeq {
    if seqs.is_empty() {
        return QualitySeq::empty();
    }
    if seqs.len() == 1 {
        return seqs[0].clone();
    }
    let uppers: Vec<QualitySeq> = seqs.iter().map(|s| s.upper()).collect();
    let max_len = uppers.iter().map(|s| s.len()).max().unwrap_or(0);

    let mut consensus_seq: Vec<u8> = Vec::new();
    let mut consensus_score: Vec<u8> = Vec::new();
    for position in 0..max_len {
        // base_stat keyed A,T,G,C in that fixed order (matches Python dict order)
        let mut stat: [(u8, i64); 4] = [(b'A', 0), (b'T', 0), (b'G', 0), (b'C', 0)];
        for s in &uppers {
            if s.len() <= position {
                continue;
            }
            let base = s.base(position);
            let q = s.qual[position] as i64;
            for entry in stat.iter_mut() {
                if entry.0 == base {
                    entry.1 += q;
                    break;
                }
            }
            // bases not in A,T,G,C (e.g. N) are ignored, matching Python.
        }
        // stable sort by score descending, preserving A,T,G,C order on ties.
        let mut sorted = stat;
        sorted.sort_by(|a, b| b.1.cmp(&a.1)); // slice sort is stable
        // SENS-7: `tolerant` extends while the best base strictly beats the
        // second-best; legacy requires it to beat the sum of all others. Both
        // truncate at the first genuinely ambiguous (tied) position.
        let delta = if tolerant {
            sorted[0].1 - sorted[1].1
        } else {
            sorted[0].1 - sorted[1].1 - sorted[2].1 - sorted[3].1
        };
        if delta > 0 {
            // cap at 255: fastq() maps every score >= 93 to the same char, so this
            // is output-identical to storing the wider score (SPD-5b).
            consensus_score.push(delta.min(255) as u8);
            consensus_seq.push(sorted[0].0);
        } else {
            break; // only extract seq to the first ambiguous base
        }
    }
    QualitySeq::new(consensus_seq, consensus_score)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::CLIP_RIGHT;

    #[test]
    fn slippage_flags_extra_copies_of_a_reference_tract() {
        // CLIP_RIGHT: [unclipped][clipped]. Reference (aligned) side ends in the poly-A tail
        // of an Alu; the clip is simply more A. This is the PD44579 artefact, e.g. hg19
        // 13:21462314 where the raw clip was AAAAAAAAAAAAAAAA.
        assert!(is_slippage_clip(
            CLIP_RIGHT,
            b"AAAAAAAAAAAAAAAA",
            b"GGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAA",
            8,
            0.7,
            1
        ));
        // CLIP_LEFT: [clipped][unclipped] — reference run at the START of unclipped.
        assert!(is_slippage_clip(CLIP_LEFT, b"TTTTTTTTTTTT", b"TTTTTTTTTTTTTTTTTTGGGTCTCGCT", 8, 0.7, 1));
    }

    #[test]
    fn slippage_spares_a_real_mei_polya_clip() {
        // A genuine MEI 3' poly-A clip: the clip is poly-A, but the reference at the junction
        // is ordinary sequence. Pairing the two sides is what keeps this call.
        assert!(!is_slippage_clip(
            CLIP_RIGHT,
            b"AAAAAAAAAAAAAAAA",
            b"GCTAGCTTACGGATCCATTGCACTGGATCA",
            8,
            0.7,
            6
        ));
        // A complex element-body clip against a reference poly-A tract is also spared:
        // the clip is not the same base as the tract.
        assert!(!is_slippage_clip(
            CLIP_RIGHT,
            b"GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCA",
            b"GGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAA",
            8,
            0.7,
            6
        ));
    }

    #[test]
    fn slippage_needs_a_long_enough_reference_run() {
        // The L1 endonuclease target motif is TTTT/AA — only 4 bp. It must not trip the
        // gate, or real MEIs inserting at their canonical EN site would be lost.
        assert!(!is_slippage_clip(CLIP_RIGHT, b"AAAAAAAAAAAA", b"GCTAGCTTACGGATCCATTGCACTTTTAA", 8, 0.7, 1));
        assert!(!is_slippage_clip(CLIP_RIGHT, b"AAAA", b"", 8, 0.7, 1));
        assert!(!is_slippage_clip(CLIP_RIGHT, b"", b"AAAAAAAAAAAA", 8, 0.7, 1));
    }

    #[test]
    fn slippage_catches_strs_only_when_max_period_allows() {
        // (TG)n microsatellite, e.g. the band locus 7:63797580. Reference ends ...TGTGTGTG,
        // the clip continues the same tract in phase. Invisible at max_period 1.
        let clip = b"GTGTGTGTGTGTGT";
        let aligned = b"GACATTGGGATCCGTCATTGTGTGTGTGTGTGTGTGTGT";
        assert!(!is_slippage_clip(CLIP_RIGHT, clip, aligned, 8, 0.7, 1));
        assert!(is_slippage_clip(CLIP_RIGHT, clip, aligned, 8, 0.7, 6));
        // (TAAAA)n — period 5, the 7:96412400 pattern.
        let clip5 = b"TAAAATAAAATAAAA";
        let al5 = b"CTCCAGCCTGGGCAACTAAAATAAAATAAAATAAAATAAAA";
        assert!(is_slippage_clip(CLIP_RIGHT, clip5, al5, 8, 0.7, 6));
    }

    #[test]
    fn slippage_str_gate_is_phase_aware_not_just_composition() {
        // Same base composition as a (TG)n tract, but the clip does not continue the tract's
        // phase — it is real inserted sequence. A composition-only test would kill this.
        let aligned = b"GACATTGGGATCCGTCATTGTGTGTGTGTGTGTGTGTGT";
        assert!(!is_slippage_clip(CLIP_RIGHT, b"TTGGTTGGTTGGTTGG", aligned, 8, 0.7, 6));
        // A homopolymer is reported as period 1, never as period 2, so raising max_period
        // must not change the poly-A verdict.
        let pa_clip = b"AAAAAAAAAAAAAAAA";
        let pa_al = b"GGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAA";
        assert_eq!(
            is_slippage_clip(CLIP_RIGHT, pa_clip, pa_al, 8, 0.7, 1),
            is_slippage_clip(CLIP_RIGHT, pa_clip, pa_al, 8, 0.7, 6)
        );
    }

    #[test]
    fn low_complexity_flags_satellites_not_elements() {
        // period-5 satellite (HSat "CCATT/GGAAT") — the chrY FP flank
        assert!(is_low_complexity(b"TCCATTCCATTCCATTCCATTCCATTCCATT", 0.8));
        // (AT)n dinucleotide satellite
        assert!(is_low_complexity(b"ATATATATATATATATATATATAT", 0.8));
        // homopolymer / dimer / SVA hexamer are all short-period
        assert!(is_low_complexity(b"AAAAAAAAAAAAAAAA", 0.8));
        assert!(is_low_complexity(b"CCCTCTCCCTCTCCCTCTCCCTCT", 0.8)); // SVA VNTR hexamer
        // a genuine complex element/flank is NOT flagged
        assert!(!is_low_complexity(b"ATTAAAAAAACTGACTTCATATGAAACAACATGCAGAAAATGC", 0.8));
        assert!(!is_low_complexity(b"CTGCACTCAAAGTGCAACTTCACACAGTGAGGGAGAAACCC", 0.8));
        // too short to judge -> not flagged
        assert!(!is_low_complexity(b"ACGT", 0.8));
    }

    #[test]
    fn mate_kmer_diversity_separates_satellite_from_complex() {
        let q = |b: &[u8]| QualitySeq::new(b.to_vec(), vec![30; b.len()]);
        // satellite mate (CCATT)n — low distinct-4-mer fraction
        let sat = vec![q(b"CCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATT")];
        // complex genomic/element mate — high diversity
        let cplx = vec![q(b"TTCCTGAGGATTGTCGTCATTAAGCTTGGTGGTAAGCTTGGGCACTCAGAGTAT")];
        let sd = mean_kmer_diversity(&sat, 4).unwrap();
        let cd = mean_kmer_diversity(&cplx, 4).unwrap();
        assert!(sd < 0.4, "satellite diversity {sd} should be low");
        assert!(cd > 0.6, "complex diversity {cd} should be high");
        // reads shorter than 20 bp are not judged
        assert!(mean_kmer_diversity(&[q(b"ACGTACGT")], 4).is_none());
        assert!(mean_kmer_diversity(&[], 4).is_none());
    }

    #[test]
    fn sens7_tolerant_extends_past_isolated_ambiguity() {
        // position 2 has A(q30) vs T(q15) vs G(q15): best beats second-best (30>15)
        // but not the sum of others (30 - 15 - 15 = 0). Legacy truncates there;
        // tolerant keeps the best base and continues to the fully-agreed pos 3.
        let seqs = vec![
            QualitySeq::new(b"AAAA".to_vec(), vec![30, 30, 30, 30]),
            QualitySeq::new(b"AATA".to_vec(), vec![30, 30, 15, 30]),
            QualitySeq::new(b"AAGA".to_vec(), vec![30, 30, 15, 30]),
        ];
        assert_eq!(find_consensus(&seqs, false).len(), 2); // legacy: stops at the tie
        assert_eq!(find_consensus(&seqs, true).len(), 4); // tolerant: extends through it
    }
}
