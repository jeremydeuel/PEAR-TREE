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
        todo!("P1")
    }
    /// `QualitySeq.upper()`.
    pub fn upper(&self) -> QualSeq {
        todo!("P1")
    }
    /// `QualitySeq.lower()`.
    pub fn lower(&self) -> QualSeq {
        todo!("P1")
    }
    /// python slice `qs[a:b]` with python clamping semantics for in-range/over-long bounds
    /// (only non-negative bounds occur in combine).
    pub fn slice(&self, a: usize, b: usize) -> QualSeq {
        todo!("P1")
    }
    /// `self + other` (concatenation of bases and qualities).
    pub fn concat(&self, other: &QualSeq) -> QualSeq {
        todo!("P1")
    }
    /// `QualitySeq.fastq(title)`: `@{title}\n{seq}\n+\n{phred}\n`, phred char =
    /// `chr(q + 33)` capped at `chr(126)`. Appends to `out`.
    pub fn fastq_into(&self, title: &str, out: &mut Vec<u8>) {
        todo!("P1")
    }
    /// python `str(seq).islower()`: at least one cased char and no uppercase char.
    pub fn is_lower(&self) -> bool {
        todo!("P1")
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
    todo!("P1")
}

/// tprt_filters.hp_compress: collapse runs of identical bytes.
pub fn hp_compress(s: &[u8]) -> Vec<u8> {
    todo!("P1")
}

/// tprt_filters.polya_trimmed(seq, polya_min=8): uppercase, cut at the first A/T run of
/// length >= polya_min, then `hp_compress`.
pub fn polya_trimmed(seq: &[u8], polya_min: usize) -> Vec<u8> {
    todo!("P1")
}

/// tprt_filters.clips_agree(seqs, min_score=0.6, polya_min=8, min_informative=6): `polya_trimmed`
/// every present clip; True when fewer than 2 remain or the shortest has < min_informative
/// bases; else `sequence_matching_score(plain) >= min_score`. Callers pass the plain bases
/// (python stringifies QualitySeqs here).
pub fn clips_agree(seqs: &[&[u8]], min_score: f64, polya_min: usize, min_informative: usize) -> bool {
    todo!("P1")
}
