//! Sequences with qualities and the small pure sequence helpers. OWNER: P1.
//!
//! Mirrors src/quality_seq.py (QualitySeq), src/revcomp.py, src/sequence_checks.py
//! (`sequence_matching_score`) and the fuzzy-merge helpers of
//! src/combine_insertions_tprt_filters.py:316-349 (`hp_compress`, `polya_trimmed`,
//! `clips_agree`). SPEC.md §2.1, §2.2.

/// python `QualitySeq`: ASCII bases (case preserved -- clips are written lowercase into the
/// consensus FASTQ) + one quality per base. Qualities are phred (discovery FASTQ, 0..=93) or
/// consensus scores (`min(93, margin)`), so `u8` holds every value; FASTQ rendering caps at
/// chr(126) exactly like `QualitySeq.phred`.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct QualSeq {
    pub seq: Box<[u8]>,
    pub qual: Box<[u8]>,
}

impl QualSeq {
    /// Panics unless `seq.len() == qual.len()` (python assert).
    pub fn new(seq: Vec<u8>, qual: Vec<u8>) -> QualSeq {
        assert_eq!(seq.len(), qual.len(), "QualSeq: seq/qual length mismatch");
        QualSeq { seq: seq.into_boxed_slice(), qual: qual.into_boxed_slice() }
    }
    pub fn len(&self) -> usize {
        self.seq.len()
    }
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
    /// `QualitySeq.revcomp()`: revcomp of the bases ([`revcomp`]), qualities reversed.
    pub fn revcomp(&self) -> QualSeq {
        let mut qual = self.qual.to_vec();
        qual.reverse();
        QualSeq { seq: revcomp(&self.seq).into_boxed_slice(), qual: qual.into_boxed_slice() }
    }
    /// `QualitySeq.upper()`.
    pub fn upper(&self) -> QualSeq {
        QualSeq { seq: self.seq.to_ascii_uppercase().into_boxed_slice(), qual: self.qual.clone() }
    }
    /// `QualitySeq.lower()`.
    pub fn lower(&self) -> QualSeq {
        QualSeq { seq: self.seq.to_ascii_lowercase().into_boxed_slice(), qual: self.qual.clone() }
    }
    /// python slice `qs[a:b]` with python clamping semantics for in-range/over-long bounds
    /// (only non-negative bounds occur in combine).
    pub fn slice(&self, a: usize, b: usize) -> QualSeq {
        let n = self.len();
        let b = b.min(n);
        let a = a.min(b);
        QualSeq { seq: self.seq[a..b].into(), qual: self.qual[a..b].into() }
    }
    /// `self + other` (concatenation of bases and qualities).
    pub fn concat(&self, other: &QualSeq) -> QualSeq {
        let mut seq = Vec::with_capacity(self.len() + other.len());
        seq.extend_from_slice(&self.seq);
        seq.extend_from_slice(&other.seq);
        let mut qual = Vec::with_capacity(seq.capacity());
        qual.extend_from_slice(&self.qual);
        qual.extend_from_slice(&other.qual);
        QualSeq { seq: seq.into_boxed_slice(), qual: qual.into_boxed_slice() }
    }
    /// `QualitySeq.fastq(title)`: `@{title}\n{seq}\n+\n{phred}\n`, phred char =
    /// `chr(q + 33)` capped at `chr(126)`. Appends to `out`.
    pub fn fastq_into(&self, title: &str, out: &mut Vec<u8>) {
        out.reserve(title.len() + 2 * self.len() + 6);
        out.push(b'@');
        out.extend_from_slice(title.as_bytes());
        out.push(b'\n');
        out.extend_from_slice(&self.seq);
        out.extend_from_slice(b"\n+\n");
        out.extend(self.qual.iter().map(|&q| (q as u32 + 33).min(126) as u8));
        out.push(b'\n');
    }
    /// python `str(seq).islower()`: at least one cased char and no uppercase char.
    pub fn is_lower(&self) -> bool {
        let mut cased = false;
        for &c in self.seq.iter() {
            if c.is_ascii_uppercase() {
                return false;
            }
            if c.is_ascii_lowercase() {
                cased = true;
            }
        }
        cased
    }
}

/// src/revcomp.py `revcomp`: complement `ATGCatgc` -> `TACGtacg`, every other byte (N, n, ...)
/// unchanged, then reverse. (indel_consensus.revcomp additionally maps N<->N / n<->n, which is
/// the identity, so ONE function serves both.) IMPLEMENTED.
pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'a' => b't',
            b't' => b'a',
            b'g' => b'c',
            b'c' => b'g',
            x => x,
        })
        .collect()
}

/// An argument of `sequence_matching_score`: a QualitySeq (quality-weighted votes) or a plain
/// str (weight 1 per base).
#[derive(Clone, Copy, Debug)]
pub enum ScoreSeq<'a> {
    Qual(&'a QualSeq),
    Plain(&'a [u8]),
}

/// src/sequence_checks.py:89 `sequence_matching_score`. Over the first
/// `min_len = min(24, min(len))` columns: tally `bases[base.upper()] += qual` (A/T/G/C/N;
/// QualSeq -> quality, str -> 1); +1 if `sum > 0 and max/sum > 0.5` else -2; return
/// `score / min_len` (f64, correctly-rounded int/int division).
/// Python crashes (KeyError) on a base outside ACGTN and (ZeroDivisionError) when
/// `min_len == 0`; this port panics with a message in both cases (SPEC.md §2.2).
pub fn sequence_matching_score(seqs: &[ScoreSeq]) -> f64 {
    assert!(!seqs.is_empty(), "sequence_matching_score: no sequences");
        let lens = |s: &ScoreSeq| match s {
            ScoreSeq::Qual(q) => q.len(),
            ScoreSeq::Plain(p) => p.len(),
        };
        let min_len = 24.min(seqs.iter().map(lens).min().unwrap());
        assert!(min_len > 0, "sequence_matching_score: zero-length sequence (python ZeroDivisionError)");
        let mut score: i64 = 0;
        for i in 0..min_len {
            // A, T, G, C, N
            let mut bases = [0u64; 5];
            for s in seqs {
                let (b, q) = match s {
                    ScoreSeq::Qual(qs) => (qs.seq[i], qs.qual[i] as u64),
                    ScoreSeq::Plain(p) => (p[i], 1u64),
                };
                let idx = match b.to_ascii_uppercase() {
                    b'A' => 0,
                    b'T' => 1,
                    b'G' => 2,
                    b'C' => 3,
                    b'N' => 4,
                    x => panic!("sequence_matching_score: base {:?} outside ACGTN (python KeyError)", x as char),
                };
                bases[idx] += q;
            }
            let sum: u64 = bases.iter().sum();
            let max = *bases.iter().max().unwrap();
            if sum > 0 && (max as f64) / (sum as f64) > 0.5 {
                score += 1;
            } else {
                score -= 2;
            }
        }
        score as f64 / min_len as f64
}

/// tprt_filters.hp_compress: collapse runs of identical bytes.
pub fn hp_compress(s: &[u8]) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::with_capacity(s.len());
    for &c in s {
        if out.last() != Some(&c) {
            out.push(c);
        }
    }
    out
}

/// tprt_filters.polya_trimmed(seq, polya_min=8): uppercase, cut at the first A/T run of
/// length >= polya_min, then `hp_compress`.
pub fn polya_trimmed(seq: &[u8], polya_min: usize) -> Vec<u8> {
    let mut s = seq.to_ascii_uppercase();
    let (mut i, n) = (0usize, s.len());
    while i < n {
        let mut k = i;
        while k < n && s[k] == s[i] {
            k += 1;
        }
        if (s[i] == b'A' || s[i] == b'T') && k - i >= polya_min {
            s.truncate(i);
            break;
        }
        i = k;
    }
    hp_compress(&s)
}

/// tprt_filters.clips_agree(seqs, min_score=0.6, polya_min=8, min_informative=6): `polya_trimmed`
/// every present clip; True when fewer than 2 remain or the shortest has < min_informative
/// bases; else `sequence_matching_score(plain) >= min_score`. Callers pass the plain bases
/// (python stringifies QualitySeqs here).
pub fn clips_agree(seqs: &[&[u8]], min_score: f64, polya_min: usize, min_informative: usize) -> bool {
    let t: Vec<Vec<u8>> = seqs.iter().map(|s| polya_trimmed(s, polya_min)).collect();
    if t.len() < 2 || t.iter().map(|x| x.len()).min().unwrap() < min_informative {
        return true;
    }
    let refs: Vec<ScoreSeq> = t.iter().map(|x| ScoreSeq::Plain(x)).collect();
    sequence_matching_score(&refs) >= min_score
}

#[cfg(test)]
mod tests {
    use super::*;

    fn qs(s: &str, q: &[u8]) -> QualSeq {
        QualSeq::new(s.as_bytes().to_vec(), q.to_vec())
    }

    /// expected values come from python `sequence_matching_score` (QualitySeq args / plain str args)
    fn check_sms(seqs: &[(&str, &[u8])], exp_qual: f64, exp_plain: f64) {
        let q: Vec<QualSeq> = seqs.iter().map(|(s, q)| qs(s, q)).collect();
        let a: Vec<ScoreSeq> = q.iter().map(ScoreSeq::Qual).collect();
        assert_eq!(sequence_matching_score(&a), exp_qual, "qual {seqs:?}");
        let b: Vec<ScoreSeq> = seqs.iter().map(|(s, _)| ScoreSeq::Plain(s.as_bytes())).collect();
        assert_eq!(sequence_matching_score(&b), exp_plain, "plain {seqs:?}");
    }

    #[test]
    fn sms_matches_python() {
        check_sms(&[("CTAAAGACAATTAGATA", &[13, 74, 73, 81, 24, 47, 12, 70, 91, 8, 72, 7, 79, 26, 63, 87, 68]), ("CTCAGTAGAT", &[40, 43, 88, 44, 76, 63, 74, 58, 8, 11]), ("GAAAAGACATTTACATAACAGTTA", &[63, 7, 27, 36, 16, 31, 50, 50, 63, 10, 21, 57, 51, 70, 35, 17, 55, 70, 35, 90, 53, 45, 87, 48]), ("CTCNCGNCAATTACA", &[23, 33, 36, 0, 18, 53, 68, 47, 78, 72, 40, 16, 88, 65, 79])], 1.0, 0.1);
        check_sms(&[("ATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGA", &[88, 20, 66, 2, 26, 67, 46, 18, 88, 69, 3, 67, 38, 82, 11, 89, 33, 66, 46, 21, 45, 28, 68, 69, 64, 42, 81, 28, 78, 24, 30, 51, 29, 25, 66, 63, 45, 93, 3, 3]), ("ATTTTTATTACACTAATAAACAGGACTCGGGTAATTTTCA", &[78, 0, 61, 83, 44, 82, 10, 84, 15, 49, 91, 25, 61, 22, 55, 81, 42, 11, 92, 50, 59, 51, 10, 92, 20, 21, 16, 3, 19, 75, 59, 83, 18, 78, 76, 60, 84, 44, 19, 70])], 1.0, 0.625);
        check_sms(&[("CAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCG", &[8, 56, 41, 78, 64, 77, 65, 25, 88, 35, 57, 65, 68, 61, 64, 31, 89, 66, 33, 71, 25, 57, 17, 53, 15, 50, 56, 40, 9, 85, 30, 54, 9, 27, 85, 38, 15, 19, 91, 82]), ("CAAACTCCAGCGCGCTCAGTTCCATCACCCTTAGTAACCG", &[45, 40, 11, 92, 46, 2, 43, 70, 58, 56, 90, 2, 49, 42, 66, 79, 37, 65, 8, 14, 29, 13, 10, 33, 34, 5, 23, 34, 16, 54, 86, 33, 51, 19, 68, 65, 73, 63, 89, 41])], 1.0, 0.875);
        check_sms(&[("GACTAGAAGA", &[79, 16, 5, 67, 90, 30, 14, 20, 33, 6]), ("GNCTCGNAGA", &[22, 34, 44, 2, 32, 4, 1, 2, 93, 64]), ("GTCTAGAAGA", &[69, 50, 64, 39, 88, 27, 29, 43, 25, 90]), ("GAATATAAGA", &[20, 7, 10, 85, 48, 64, 85, 36, 76, 31])], 0.7, 0.7);
        check_sms(&[("GAACCCTAGGGGCAGCGCCGTATCACAAGACGTCATGAAC", &[0, 58, 8, 64, 68, 11, 84, 67, 8, 60, 32, 9, 33, 30, 93, 26, 29, 83, 58, 63, 48, 9, 61, 87, 36, 5, 78, 80, 82, 25, 9, 76, 18, 42, 32, 83, 88, 38, 79, 72]), ("GATCCGTAGGGGCA", &[37, 90, 66, 36, 59, 59, 59, 15, 70, 25, 39, 10, 60, 2]), ("gntccgtngaggctgcgtagtatgccaagactataggcac", &[46, 29, 63, 62, 50, 3, 20, 0, 62, 87, 57, 51, 38, 93, 18, 53, 44, 48, 40, 15, 42, 0, 41, 43, 50, 15, 25, 91, 1, 37, 32, 47, 8, 50, 49, 75, 9, 46, 54, 35]), ("GATCCGTAGGG", &[19, 31, 34, 55, 65, 40, 24, 47, 54, 3, 80])], 1.0, 1.0);
        check_sms(&[("CAATTCATACCTTGGGGT", &[17, 70, 24, 31, 11, 22, 43, 71, 11, 40, 30, 47, 33, 72, 25, 2, 52, 49]), ("GCCTGCGGACGTTTGAGGGCGTTA", &[62, 0, 9, 50, 67, 59, 57, 31, 13, 28, 19, 19, 66, 87, 13, 92, 89, 82, 58, 10, 70, 5, 0, 16])], 1.0, -0.6666666666666666);
        check_sms(&[("AGCGTGACC", &[2, 24, 63, 86, 82, 53, 10, 32, 29]), ("cgcgcagcgc", &[28, 62, 53, 85, 7, 76, 18, 50, 6, 27]), ("ngcgtnaagc", &[7, 23, 50, 57, 91, 40, 93, 14, 10, 21])], 0.6666666666666666, 0.3333333333333333);
        check_sms(&[("TTTNGTTTN", &[40, 35, 38, 0, 92, 76, 81, 8, 3]), ("CCTAGTGGTCAATGAATACTGGTA", &[16, 63, 23, 1, 38, 88, 19, 77, 30, 41, 40, 58, 46, 76, 10, 65, 25, 50, 20, 31, 52, 8, 83, 4])], 1.0, -1.0);
        check_sms(&[("GCTAACACACCTACC", &[12, 83, 59, 4, 13, 0, 60, 29, 57, 47, 5, 37, 29, 15, 6]), ("GCCAAGACATTTCCCTGCAGGGGG", &[85, 0, 13, 81, 76, 90, 79, 44, 27, 4, 47, 43, 18, 5, 26, 32, 4, 76, 93, 83, 26, 1, 41, 52]), ("gcttagacatttccctttacgggg", &[68, 11, 83, 20, 50, 89, 34, 52, 36, 85, 39, 53, 6, 39, 72, 45, 53, 53, 2, 46, 82, 25, 50, 93])], 1.0, 1.0);
        check_sms(&[("CATCTAATGTCCAACTAGCCGGCC", &[25, 38, 16, 5, 61, 40, 6, 77, 81, 49, 11, 91, 79, 88, 20, 81, 28, 79, 51, 78, 25, 60, 23, 72]), ("CNTCTNATGTCCCACTAGCCGGCC", &[24, 5, 71, 86, 4, 85, 41, 15, 49, 76, 58, 70, 80, 39, 83, 53, 39, 74, 31, 54, 49, 84, 47, 57])], 1.0, 0.625);
    }
    #[test]
    fn polya_and_hp_match_python() {
        assert_eq!(polya_trimmed(b"ACGTTTTTTTTTTGGA", 8), b"ACG");
        assert_eq!(hp_compress(b"ACGTTTTTTTTTTGGA"), b"ACGTGA");
        assert_eq!(polya_trimmed(b"acgtacgtAAAAAAAAAAAAcc", 8), b"ACGTACGT");
        assert_eq!(hp_compress(b"acgtacgtAAAAAAAAAAAAcc"), b"acgtacgtAc");
        assert_eq!(polya_trimmed(b"AAAAAAAAAAAA", 8), b"");
        assert_eq!(hp_compress(b"AAAAAAAAAAAA"), b"A");
        assert_eq!(polya_trimmed(b"GGCCTTAAGGCAAAAAAAGGT", 8), b"GCTAGCAGT");
        assert_eq!(hp_compress(b"GGCCTTAAGGCAAAAAAAGGT"), b"GCTAGCAGT");
        assert_eq!(polya_trimmed(b"ACGGGGTTTTTTTTTTT", 8), b"ACG");
        assert_eq!(hp_compress(b"ACGGGGTTTTTTTTTTT"), b"ACGT");
        assert_eq!(polya_trimmed(b"ACACACACGTGT", 8), b"ACACACACGTGT");
        assert_eq!(hp_compress(b"ACACACACGTGT"), b"ACACACACGTGT");
        assert_eq!(polya_trimmed(b"", 8), b"");
        assert_eq!(hp_compress(b""), b"");
        assert_eq!(polya_trimmed(b"NNNNACGTAAAAAAAA", 8), b"NACGT");
        assert_eq!(hp_compress(b"NNNNACGTAAAAAAAA"), b"NACGTA");
        assert_eq!(polya_trimmed(b"CCCCCCCCCCCCAAAAAAAAA", 8), b"C");
        assert_eq!(hp_compress(b"CCCCCCCCCCCCAAAAAAAAA"), b"CA");
        assert_eq!(polya_trimmed(b"tttttttttttttACGT", 8), b"");
        assert_eq!(hp_compress(b"tttttttttttttACGT"), b"tACGT");
        assert_eq!(polya_trimmed(b"ACGTAAAAAAATTTTTTTT", 8), b"ACGTA");
        assert_eq!(hp_compress(b"ACGTAAAAAAATTTTTTTT"), b"ACGTAT");
        assert_eq!(polya_trimmed(b"ACGTAAAAAAAT", 5), b"ACGT");
    }
    #[test]
    fn clips_agree_matches_python() {
        assert_eq!(clips_agree(&[b"TCAATTCTTCTTAACGTGATAACAGAATCA", b"TCAATTCTTCTTAACGTGATAACAGAATCA", b"CCTGCCAGGCGGTCGTCGCGGACCTCGGTC"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"AAGTAGTGGTGCGGATCCAGGGGAACCGTTAAAAAAAAAAAAAAA", b"AAGTAGTGATGCGGTTCCAGGGGAACCGTTAAAAAAAAAAAAAAA"], 0.6, 8, 6), false);
        assert_eq!(clips_agree(&[b"ATCCCAAACCTCTCGAGATATTTATCCAGCAAAAAAAAAAAAAAA"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"ACCAAAACGCAAACAAAAGCATACCCAAAAAAAAAAAAAAAAAAA"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"CTGGCGCCTCAATAGGATTATAGCGGTCTCAAAAAAAAA"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"CCGGTGCAAGCTTAATTCGTACGTACTTCCAAAAAAAAAAAAAAA", b"ACG"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"AGGTTCCTAGAGGTTAAATTGGACGTCTTC", b"AGGTTCCTAGAGGTTAAATTGGACGACTTC", b"CTCCGTTGCTGCGTGTCTAGGCGGTTTAGC"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"TAAGCGAACAGGACCCTGCCTCAGCTCATAAAAAAAAAAAAAAAA", b"TAAGCGAACAGGACCCTGCCTCAGCTCATAAAAAAAAAAAAAAAA", b"GTCCTTATTCTCTCACGTTGTGTTACGAAA"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"ATTCACTCGAGGTCGTGTGAGGGTTGGGCT", b"ATTCGCTCGACGACGTGTGAGGGTTAGGCT", b"ATGAAACTATCACATCACATAAGCGGGCTA"], 0.6, 8, 6), false);
        assert_eq!(clips_agree(&[b"ATATAATTTAATCTTAATCCATAAAACACT", b"ATATAATTTAAACTCAATCCATACAACACC", b"TTGAAAAAATGGCTAGGTTCCAGCTTTTGG"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"GAGACGTCTTTCTGAGGGTCAGCCGTGATT", b"GAGACGTCTTTCTGAGGGTGAGCCGTGATT", b"ATTCGATTAGACTGGTCCCCACGGGTCCAT"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"AGTACGAGGAAACTCGGTATCGAGCCTAAA", b"AGTACGAGAAAACTCTGTATCGAGCCTAAT", b"GCATCTCGCCCAGGAAAGTAACGACGTATG"], 0.6, 8, 6), false);
        assert_eq!(clips_agree(&[b"GTAGTTCTCCATCACCAGCTATAATGGCTAAAAAAAAAAAAAAAA"], 0.6, 8, 6), true);
        assert_eq!(clips_agree(&[b"GTCAGCATGCTAGCGTATCGCCCCCCAATG", b"GTCAGCATGCTAGCGTATCCCCCCCCAATG", b"CGCAATAGGGTAATTCGCCGACGAGTAAGC"], 0.6, 8, 6), true);
    }
    #[test]
    fn qualseq_ops_match_python() {
        {
            let q = QualSeq::new(b"cnANCggtAttGATCTn".to_vec(), [13, 71, 93, 89, 33, 76, 59, 89, 44, 93, 15, 74, 35, 15, 17, 29, 93].to_vec());
            assert_eq!(q.revcomp(), qs("nAGATCaaTaccGNTng", &[93, 29, 17, 15, 35, 74, 15, 93, 44, 89, 59, 76, 33, 89, 93, 71, 13]));
            assert_eq!(q.upper(), qs("CNANCGGTATTGATCTN", &[13, 71, 93, 89, 33, 76, 59, 89, 44, 93, 15, 74, 35, 15, 17, 29, 93]));
            assert_eq!(q.lower(), qs("cnancggtattgatctn", &[13, 71, 93, 89, 33, 76, 59, 89, 44, 93, 15, 74, 35, 15, 17, 29, 93]));
            assert_eq!(q.slice(5, 21), qs("ggtAttGATCTn", &[76, 59, 89, 44, 93, 15, 74, 35, 15, 17, 29, 93]));
            assert_eq!(q.concat(&qs("ACGT", &[1, 2, 3, 40])), qs("cnANCggtAttGATCTnACGT", &[13, 71, 93, 89, 33, 76, 59, 89, 44, 93, 15, 74, 35, 15, 17, 29, 93, 1, 2, 3, 40]));
            assert_eq!(q.is_lower(), false);
            let mut out = Vec::new();
            q.fastq_into("t:1-2:LEFT:CLIPPED", &mut out);
            assert_eq!(String::from_utf8(out).unwrap(), "@t:1-2:LEFT:CLIPPED\ncnANCggtAttGATCTn\n+\n.h~zBm\\zM~0kD02>~\n");
        }
        {
            let q = QualSeq::new(b"cGAccnntAcAa".to_vec(), [30, 0, 0, 93, 18, 0, 54, 1, 67, 41, 64, 28].to_vec());
            assert_eq!(q.revcomp(), qs("tTgTannggTCg", &[28, 64, 41, 67, 1, 54, 0, 18, 93, 0, 0, 30]));
            assert_eq!(q.upper(), qs("CGACCNNTACAA", &[30, 0, 0, 93, 18, 0, 54, 1, 67, 41, 64, 28]));
            assert_eq!(q.lower(), qs("cgaccnntacaa", &[30, 0, 0, 93, 18, 0, 54, 1, 67, 41, 64, 28]));
            assert_eq!(q.slice(1, 16), qs("GAccnntAcAa", &[0, 0, 93, 18, 0, 54, 1, 67, 41, 64, 28]));
            assert_eq!(q.concat(&qs("ACGT", &[1, 2, 3, 40])), qs("cGAccnntAcAaACGT", &[30, 0, 0, 93, 18, 0, 54, 1, 67, 41, 64, 28, 1, 2, 3, 40]));
            assert_eq!(q.is_lower(), false);
            let mut out = Vec::new();
            q.fastq_into("t:1-2:LEFT:CLIPPED", &mut out);
            assert_eq!(String::from_utf8(out).unwrap(), "@t:1-2:LEFT:CLIPPED\ncGAccnntAcAa\n+\n?!!~3!W\"dJa=\n");
        }
        {
            let q = QualSeq::new(b"gAAAnNnNtAn".to_vec(), [15, 55, 5, 39, 21, 76, 0, 75, 18, 65, 37].to_vec());
            assert_eq!(q.revcomp(), qs("nTaNnNnTTTc", &[37, 65, 18, 75, 0, 76, 21, 39, 5, 55, 15]));
            assert_eq!(q.upper(), qs("GAAANNNNTAN", &[15, 55, 5, 39, 21, 76, 0, 75, 18, 65, 37]));
            assert_eq!(q.lower(), qs("gaaannnntan", &[15, 55, 5, 39, 21, 76, 0, 75, 18, 65, 37]));
            assert_eq!(q.slice(3, 21), qs("AnNnNtAn", &[39, 21, 76, 0, 75, 18, 65, 37]));
            assert_eq!(q.concat(&qs("ACGT", &[1, 2, 3, 40])), qs("gAAAnNnNtAnACGT", &[15, 55, 5, 39, 21, 76, 0, 75, 18, 65, 37, 1, 2, 3, 40]));
            assert_eq!(q.is_lower(), false);
            let mut out = Vec::new();
            q.fastq_into("t:1-2:LEFT:CLIPPED", &mut out);
            assert_eq!(String::from_utf8(out).unwrap(), "@t:1-2:LEFT:CLIPPED\ngAAAnNnNtAn\n+\n0X&H6m!l3bF\n");
        }
        {
            let q = QualSeq::new(b"NTCtNgnnT".to_vec(), [25, 46, 70, 61, 39, 42, 65, 74, 45].to_vec());
            assert_eq!(q.revcomp(), qs("AnncNaGAN", &[45, 74, 65, 42, 39, 61, 70, 46, 25]));
            assert_eq!(q.upper(), qs("NTCTNGNNT", &[25, 46, 70, 61, 39, 42, 65, 74, 45]));
            assert_eq!(q.lower(), qs("ntctngnnt", &[25, 46, 70, 61, 39, 42, 65, 74, 45]));
            assert_eq!(q.slice(1, 10), qs("TCtNgnnT", &[46, 70, 61, 39, 42, 65, 74, 45]));
            assert_eq!(q.concat(&qs("ACGT", &[1, 2, 3, 40])), qs("NTCtNgnnTACGT", &[25, 46, 70, 61, 39, 42, 65, 74, 45, 1, 2, 3, 40]));
            assert_eq!(q.is_lower(), false);
            let mut out = Vec::new();
            q.fastq_into("t:1-2:LEFT:CLIPPED", &mut out);
            assert_eq!(String::from_utf8(out).unwrap(), "@t:1-2:LEFT:CLIPPED\nNTCtNgnnT\n+\n:Og^HKbkN\n");
        }
        let mut out = Vec::new();
        qs("ACG", &[93, 100, 255]).fastq_into("x", &mut out);
        assert_eq!(String::from_utf8(out).unwrap(), "@x\nACG\n+\n~~~\n");
    }

    #[test]
    #[should_panic]
    fn sms_zero_length_panics() {
        sequence_matching_score(&[ScoreSeq::Plain(b""), ScoreSeq::Plain(b"ACGT")]);
    }
}
