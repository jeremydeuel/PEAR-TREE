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
