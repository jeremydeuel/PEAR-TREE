//! Ports of src/adapter.py (is_adapter), src/sequence_checks.py
//! (clean_clipped_seq) and src/consensus.py (find_consensus).

use crate::config::{ADAPTERS, MIN_ADAPTERLEN_FOR_CLIP, MIN_CLIP_LEN};
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
