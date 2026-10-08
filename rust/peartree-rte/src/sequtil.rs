//! tools/rte/sequtil.py -- small sequence helpers. FOUNDATION (implemented, unit-tested).
//!
//! Sequences are ASCII byte slices. Functions that upper-case in Python upper-case here too.

use crate::align::{self, Mode};
use crate::pyfmt::py_sum;

/// `rc(seq)`: complement via `str.maketrans("ACGTNRYKMacgtnrykm", "TGCANYRMKtgcanyrmk")`
/// (case kept; every other byte unchanged), reversed.
pub fn rc(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| comp(b)).collect()
}

#[inline]
pub fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'n' => b'n',
        b'r' => b'y',
        b'y' => b'r',
        b'k' => b'm',
        b'm' => b'k',
        o => o,
    }
}

pub fn upper(seq: &[u8]) -> Vec<u8> {
    seq.to_ascii_uppercase()
}

/// `polya_runs(seq, base, min_len, max_gap=1)`: maximal runs of `base` (>= 3, regex `B{3,}`,
/// case-insensitive) merged across gaps of <= `max_gap`, kept when >= `min_len`. Half-open.
pub fn polya_runs(seq: &[u8], base: u8, min_len: usize, max_gap: usize) -> Vec<(usize, usize)> {
    let b = base.to_ascii_uppercase();
    let mut raw = Vec::new();
    let mut i = 0;
    while i < seq.len() {
        if seq[i].to_ascii_uppercase() == b {
            let st = i;
            while i < seq.len() && seq[i].to_ascii_uppercase() == b {
                i += 1;
            }
            if i - st >= 3 {
                raw.push((st, i));
            }
        } else {
            i += 1;
        }
    }
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (st, en) in raw {
        if let Some(last) = merged.last_mut() {
            if st - last.1 <= max_gap {
                last.1 = en;
                continue;
            }
        }
        merged.push((st, en));
    }
    merged.into_iter().filter(|&(a, b)| b - a >= min_len).collect()
}

/// `trailing_polya_start(seq, min_len=5)`: index where a trailing poly-A (allowing up to 2 'N')
/// starts, or len(seq) when the run is shorter than `min_len`.
pub fn trailing_polya_start(seq: &[u8], min_len: usize) -> usize {
    let n = seq.len();
    let mut i = n;
    let mut bad = 0;
    while i > 0 {
        let c = seq[i - 1].to_ascii_uppercase();
        if c == b'A' || (c == b'N' && bad < 2) {
            if c != b'A' {
                bad += 1;
            }
            i -= 1;
        } else {
            break;
        }
    }
    if n - i >= min_len {
        i
    } else {
        n
    }
}

/// `homopolymer_at(seq, pos, direction)` -> (base (upper-case, 0 when out of range), length).
pub fn homopolymer_at(seq: &[u8], pos: i64, direction: i64) -> (u8, usize) {
    if pos < 0 || pos as usize >= seq.len() {
        return (0, 0);
    }
    let b = seq[pos as usize].to_ascii_uppercase();
    let mut n = 0;
    let mut i = pos;
    while i >= 0 && (i as usize) < seq.len() && seq[i as usize].to_ascii_uppercase() == b {
        n += 1;
        i += direction;
    }
    (b, n)
}

/// `hamming(a, b)`: mismatches over the common prefix length (case-insensitive) + |len diff|.
pub fn hamming(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b)
        .filter(|(x, y)| !x.eq_ignore_ascii_case(y))
        .count()
        + a.len().abs_diff(b.len())
}

/// `shannon(seq)` in bits; the per-character terms are summed in first-occurrence order with
/// Python's `sum()`.
pub fn shannon(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let mut order: Vec<u8> = Vec::new();
    let mut counts = [0usize; 256];
    for &c in seq {
        let u = c.to_ascii_uppercase();
        if counts[u as usize] == 0 {
            order.push(u);
        }
        counts[u as usize] += 1;
    }
    let n = seq.len() as f64;
    -py_sum(order.iter().map(|&c| {
        let p = counts[c as usize] as f64 / n;
        p * p.log2()
    }))
}

/// One `edlib_best` result: (edit distance, target start, target end (exclusive), strand).
pub type EdHit = (i32, i64, i64, i32);

/// `edlib_best(query, target, max_frac=0.15, both_strands=True)`: best infix (HW) alignment of
/// `query` (then of rc(query)) inside `target`, k = max(1, int(len(query) * max_frac)); strict `<`
/// so the + strand wins ties. None when neither strand is within k (or an input is empty).
pub fn edlib_best(query: &[u8], target: &[u8], max_frac: f64, both_strands: bool) -> Option<EdHit> {
    if query.is_empty() || target.is_empty() {
        return None;
    }
    let k = std::cmp::max(1, (query.len() as f64 * max_frac) as i64) as i32;
    let t = target.to_ascii_uppercase();
    let q = query.to_ascii_uppercase();
    let mut best: Option<EdHit> = None;
    for strand in [1, -1] {
        if strand == -1 && !both_strands {
            break;
        }
        let qq = if strand == 1 { q.clone() } else { rc(&q) };
        if let Some((ed, st, en)) = align::locate(&qq, &t, Mode::Hw, k) {
            let cand = (ed, st as i64, en as i64 + 1, strand);
            if best.is_none_or(|b| cand.0 < b.0) {
                best = Some(cand);
            }
        }
    }
    best
}

/// (edit distance, t_start, t_end (exclusive), ops [(len, code)])
pub type EdPath = (i32, i64, i64, Vec<(usize, u8)>);

/// `edlib_path(query, target, k=-1)`: (edit distance, t_start, t_end (exclusive), ops) with ops
/// `[(len, code)]`, code 0 = '='/'X', 1 = 'I', 2 = 'D' (minimap2 integers), consecutive equal
/// codes merged. None when the distance exceeds k.
pub fn edlib_path(query: &[u8], target: &[u8], k: i32) -> Option<EdPath> {
    let q = query.to_ascii_uppercase();
    let t = target.to_ascii_uppercase();
    let r = align::path(&q, &t, Mode::Hw, k, &[])?;
    let mut ops: Vec<(usize, u8)> = Vec::new();
    for &op in &r.ops {
        let code = match op {
            align::OP_MATCH | align::OP_MISMATCH => 0u8,
            align::OP_INS => 1,
            _ => 2,
        };
        match ops.last_mut() {
            Some(last) if last.1 == code => last.0 += 1,
            _ => ops.push((1, code)),
        }
    }
    Some((r.edit_distance, r.start as i64, r.end as i64 + 1, ops))
}

/// `base_fraction(seq, base)` (case-insensitive).
pub fn base_fraction(seq: &[u8], base: u8) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let b = base.to_ascii_uppercase();
    seq.iter().filter(|c| c.to_ascii_uppercase() == b).count() as f64 / seq.len() as f64
}

/// `low_complexity(seq, max_base_frac=0.7, min_entropy=1.5)`.
pub fn low_complexity(seq: &[u8], max_base_frac: f64, min_entropy: f64) -> bool {
    if seq.is_empty() {
        return true;
    }
    let s = seq.to_ascii_uppercase();
    let mx = b"ACGT".iter().map(|&b| s.iter().filter(|&&c| c == b).count()).max().unwrap();
    if mx as f64 / s.len() as f64 >= max_base_frac {
        return true;
    }
    shannon(&s) < min_entropy
}

/// `low_complexity(seq)` with the Python defaults.
pub fn low_complexity_default(seq: &[u8]) -> bool {
    low_complexity(seq, 0.7, 1.5)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_helpers_like_python() {
        // python: rc("ACGTNacgtR*") -> "*YacgtNACGT"
        assert_eq!(rc(b"ACGTNacgtR*"), b"*YacgtNACGT");
        assert_eq!(rc(b"AACG"), b"CGTT");
        // python: polya_runs("CCAAAAAAAAGAAAAAACC", "A", 8) -> [(2, 17)]
        assert_eq!(polya_runs(b"CCAAAAAAAAGAAAAAACC", b'A', 8, 1), vec![(2, 17)]);
        assert_eq!(polya_runs(b"AAAAAAAACAAA", b'a', 8, 1), vec![(0, 12)]);
        assert_eq!(polya_runs(b"AAAAAAAACCAAA", b'A', 8, 1), vec![(0, 8)]);
        // python: trailing_polya_start("ACGTAAAAAA") -> 4 ; ("ACGTAAAA") -> 8 ; ("CANAAAA") -> 1
        assert_eq!(trailing_polya_start(b"ACGTAAAAAA", 5), 4);
        assert_eq!(trailing_polya_start(b"ACGTAAAA", 5), 8);
        assert_eq!(trailing_polya_start(b"CANAAAA", 5), 1);
        assert_eq!(homopolymer_at(b"GGAAAAT", 2, 1), (b'A', 4));
        assert_eq!(homopolymer_at(b"GGAAAAT", 5, -1), (b'A', 4));
        assert_eq!(homopolymer_at(b"GG", 5, -1), (0, 0));
        assert_eq!(hamming(b"ACGT", b"acgA"), 1);
        assert_eq!(hamming(b"ACGT", b"AC"), 2);
        // python: shannon("AACG") = 1.5
        assert!((shannon(b"AACG") - 1.5).abs() < 1e-12);
        assert!(low_complexity_default(b"AAAAAAAAAC"));
        assert!(!low_complexity_default(b"ACGTTGCAAGCT"));
        assert!(low_complexity_default(b""));
    }

    #[test]
    fn edlib_wrappers_like_python() {
        // python: edlib_best("ACGTACGTTT", "GGGGAAACGTACGTTTGGG") -> (0, 6, 16, 1)
        assert_eq!(edlib_best(b"ACGTACGTTT", b"GGGGAAACGTACGTTTGGG", 0.15, true), Some((0, 6, 16, 1)));
        // reverse-complement match
        let q = b"ACGTACGTTTGCA";
        let t = [b"CCCC".as_slice(), &rc(q), b"CCCC"].concat();
        assert_eq!(edlib_best(q, &t, 0.15, true), Some((0, 4, 17, -1)));
        assert_eq!(edlib_best(q, &t, 0.15, false), None);
        // python: edlib_path("ACGTTACG", "TTACGTACGTT") -> (1, 2, 9, [(4, 0), (1, 1), (3, 0)])
        let (ed, st, en, ops) = edlib_path(b"ACGTTACG", b"TTACGTACGTT", -1).unwrap();
        assert_eq!((ed, st, en), (1, 2, 9));
        assert_eq!(ops, vec![(4, 0), (1, 1), (3, 0)]);
    }
}
