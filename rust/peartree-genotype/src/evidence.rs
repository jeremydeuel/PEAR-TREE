//! Per-read breakpoint scoring — port of src/genotyping_evidence_read.py and
//! src/genotype_qscore.py.
//!
//! A spanning read is scored on each junction it covers: the read bases on one
//! side of the breakpoint are compared, base-by-base weighted by phred quality,
//! against the reference-flank consensus (`ref`) and the inserted-junction
//! consensus (`alt`). A ±1 bp register search tolerates breakpoint imprecision.

pub const LEFT_TO_RIGHT: u8 = 1;
pub const RIGHT_TO_LEFT: u8 = 2;

// BAM CIGAR op codes grouped by what they consume (matches pysam).
#[inline]
fn consumes_query(op: u8) -> bool {
    matches!(op, 0 | 1 | 4 | 7 | 8) // M I S = X
}
#[inline]
fn consumes_ref(op: u8) -> bool {
    matches!(op, 0 | 2 | 3 | 7 | 8) // M D N = X
}
#[inline]
fn is_match_op(op: u8) -> bool {
    matches!(op, 0 | 7 | 8) // M = X
}

/// Query offset aligned to reference position `target_ref`, or None if that
/// reference base is not aligned to a query base (deletion/skip, or outside the
/// read). Port of `query_index_at_ref`: O(#cigar ops), no list materialisation.
pub fn query_index_at_ref(cigar: &[(u8, usize)], reference_start: i64, target_ref: i64) -> Option<usize> {
    let mut qpos: usize = 0;
    let mut rpos: i64 = reference_start;
    for &(op, len) in cigar {
        let len_i = len as i64;
        if is_match_op(op) && rpos <= target_ref && target_ref < rpos + len_i {
            return Some(qpos + (target_ref - rpos) as usize);
        }
        if consumes_query(op) {
            qpos += len;
        }
        if consumes_ref(op) {
            rpos += len_i;
        }
    }
    None
}

/// One register pass over the aligned read side. Returns (ref, alt, art, total_q).
///
/// `s_skip`/`ra_skip` shift the read vs. the consensus by one base: (0,0) is the
/// no-shift register, (1,0) slides the read by one, (0,1) slides the consensus.
/// `dir == RIGHT_TO_LEFT` pairs the LAST elements of each slice (Python
/// `reversed(...)`); `LEFT_TO_RIGHT` pairs from the front. Iteration length is the
/// min over all three slices after their skips, matching Python `zip` truncation.
fn pass(
    dir: u8,
    seq: &[u8],
    qual: &[u8],
    refc: &[u8],
    altc: &[u8],
    s_skip: usize,
    ra_skip: usize,
) -> (i64, i64, i64, i64) {
    let (ls, lr, la) = (seq.len(), refc.len(), altc.len());
    let n = ls.saturating_sub(s_skip).min(lr.saturating_sub(ra_skip)).min(la.saturating_sub(ra_skip));
    let (mut r, mut a, mut art, mut tq) = (0i64, 0i64, 0i64, 0i64);
    for idx in 0..n {
        let (si, ri, ai) = if dir == LEFT_TO_RIGHT {
            (s_skip + idx, ra_skip + idx, ra_skip + idx)
        } else {
            (ls - 1 - s_skip - idx, lr - 1 - ra_skip - idx, la - 1 - ra_skip - idx)
        };
        let s = seq[si];
        if s == b'N' {
            continue;
        }
        let q = qual[si] as i64;
        tq += q;
        if s == refc[ri] {
            r += q;
        } else {
            r -= q;
        }
        if s == altc[ai] {
            a += q;
        } else {
            a -= q;
        }
        if s != altc[ai] && s != refc[ri] {
            art += q;
        }
    }
    (r, a, art, tq)
}

/// Score one aligned read side against (ref, alt). Returns (ref_score, alt_score,
/// art_score). Port of `qscore`, including the perfect-match short-circuit: if the
/// no-shift register already matches one hypothesis at every scored base, it is
/// optimal and the two ±1 shift passes are skipped.
pub fn qscore(seq: &[u8], qual: &[u8], refc: &[u8], altc: &[u8], dir: u8) -> (i64, i64, i64) {
    let (r1, a1, art1, total_q) = pass(dir, seq, qual, refc, altc, 0, 0);
    let s1 = (r1, a1, art1);
    if r1 == total_q || a1 == total_q {
        return s1;
    }
    let s2 = pass(dir, seq, qual, refc, altc, 1, 0);
    let s3 = pass(dir, seq, qual, refc, altc, 0, 1);
    let mut best = s1;
    let mut best_key = r1.max(a1);
    for cand in [(s2.0, s2.1, s2.2), (s3.0, s3.1, s3.2)] {
        let ck = cand.0.max(cand.1);
        if ck > best_key {
            best = cand;
            best_key = ck;
        }
    }
    best
}

/// Score the read portion 5' of `breakpoint` (RIGHT_TO_LEFT). `Ok((0,0,0))` when the
/// breakpoint is outside the read or not aligned to a query base. `Err(())` mirrors
/// the Python crash when a consensus is missing but a spanning read reaches scoring
/// (the driver turns it into a locus-level `error`).
pub fn qleft(
    cigar: &[(u8, usize)],
    reference_start: i64,
    reference_end: i64,
    seq: &[u8],
    qual: &[u8],
    breakpoint: i64,
    refc: Option<&[u8]>,
    altc: Option<&[u8]>,
) -> Result<(i64, i64, i64), ()> {
    if breakpoint < reference_start || breakpoint > reference_end {
        return Ok((0, 0, 0));
    }
    match query_index_at_ref(cigar, reference_start, breakpoint) {
        None => Ok((0, 0, 0)),
        Some(bq) => {
            let (refc, altc) = (refc.ok_or(())?, altc.ok_or(())?);
            Ok(qscore(&seq[..bq], &qual[..bq], refc, altc, RIGHT_TO_LEFT))
        }
    }
}

/// Score the read portion 3' of `breakpoint` (LEFT_TO_RIGHT). Uses the query base
/// aligned just 5' of the breakpoint (ref == breakpoint-1), then +1, matching the
/// Python `qright` semantics.
pub fn qright(
    cigar: &[(u8, usize)],
    reference_start: i64,
    reference_end: i64,
    seq: &[u8],
    qual: &[u8],
    breakpoint: i64,
    refc: Option<&[u8]>,
    altc: Option<&[u8]>,
) -> Result<(i64, i64, i64), ()> {
    if breakpoint < reference_start || breakpoint > reference_end {
        return Ok((0, 0, 0));
    }
    match query_index_at_ref(cigar, reference_start, breakpoint - 1) {
        None => Ok((0, 0, 0)),
        Some(q) => {
            let bq = q + 1;
            let (refc, altc) = (refc.ok_or(())?, altc.ok_or(())?);
            Ok(qscore(&seq[bq..], &qual[bq..], refc, altc, LEFT_TO_RIGHT))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn query_index_simple_match() {
        // 10M starting at ref 100: ref 105 -> query offset 5.
        let cigar = [(0u8, 10usize)];
        assert_eq!(query_index_at_ref(&cigar, 100, 105), Some(5));
        assert_eq!(query_index_at_ref(&cigar, 100, 100), Some(0));
        assert_eq!(query_index_at_ref(&cigar, 100, 110), None); // past the end
    }

    #[test]
    fn query_index_with_softclip_and_del() {
        // 3S5M2D5M from ref 100: soft clip advances query only; deletion advances ref only.
        let cigar = [(4u8, 3usize), (0, 5), (2, 2), (0, 5)];
        // ref 100 -> first M, query 3 (after the 3S)
        assert_eq!(query_index_at_ref(&cigar, 100, 100), Some(3));
        // ref 104 -> last base of first M block, query 3+4 = 7
        assert_eq!(query_index_at_ref(&cigar, 100, 104), Some(7));
        // ref 105/106 fall in the 2D deletion -> no query base
        assert_eq!(query_index_at_ref(&cigar, 100, 105), None);
        // ref 107 -> first base of second M, query 3+5 = 8
        assert_eq!(query_index_at_ref(&cigar, 100, 107), Some(8));
    }

    #[test]
    fn qscore_perfect_alt_shortcircuits_and_reports_alt() {
        let seq = b"ACGT";
        let qual = [40u8, 40, 40, 40];
        let refc = b"TTTT";
        let altc = b"ACGT"; // read matches alt exactly
        let (r, a, art) = qscore(seq, &qual, refc, altc, LEFT_TO_RIGHT);
        assert_eq!(a, 160); // 4 * 40, matches alt at every base -> short-circuits
        // vs refc "TTTT": only the last base (T) matches -> +40, three mismatch -40 each
        assert_eq!(r, -80);
        assert_eq!(art, 0);
    }

    #[test]
    fn qscore_skips_n_bases() {
        let seq = b"ANGT";
        let qual = [40u8, 40, 40, 40];
        let refc = b"ANGT";
        let altc = b"CCCC";
        // 'N' at index1 is skipped entirely; the other 3 match ref -> total_q 120, r 120.
        let (r, _a, _art) = qscore(seq, &qual, refc, altc, LEFT_TO_RIGHT);
        assert_eq!(r, 120);
    }
}
