//! Indel-aware (homopolymer-compressed star alignment) clip consensus. OWNER: P2.
//!
//! Mirrors src/indel_consensus.py line by line. SPEC.md §4.4 -- the numerics section there is
//! normative: dict-insertion-order maps for column votes, python `max()` first-wins ties,
//! Neumaier `py_sum` for `others`, `py_round` (half-even) for the margin, `py_median`,
//! edlib through `crate::align` (SHW for anchored reads / HW for floating reads, wildcard N pad).

/// python `ClipRead` (indel_consensus.py:48): one read, oriented OUTWARD from the junction.
#[derive(Clone, Debug)]
pub struct ClipRead {
    /// bases as given (rle() uppercases)
    pub seq: Vec<u8>,
    /// one quality per base (phred ints)
    pub qual: Vec<u8>,
    /// independent-fragment cluster id (python `gi`): reads of one cluster count once for depth
    pub group: u32,
    /// `1.0 / len(cluster)`
    pub weight: f64,
    /// false = floating read (a mate) that may start anywhere and either strand
    pub anchored: bool,
}

/// python `ConsensusResult` (indel_consensus.py:61). `Default` == `ConsensusResult()`
/// (empty seq, `stop_reason = "empty"`, no poly-A, beyond_* empty / -1 / 0).
#[derive(Clone, Debug, PartialEq)]
pub struct ConsensusResult {
    /// outward orientation, uppercase
    pub seq: Vec<u8>,
    /// independent fragments per base
    pub depth: Vec<u32>,
    /// `min(93, margin)` per base
    pub score: Vec<u8>,
    /// "empty" | "end" | "disagreement" | "depth"
    pub stop_reason: &'static str,
    /// b'A' / b'T'
    pub polya_base: Option<u8>,
    pub polya_start: i64,
    pub polya_len_median: Option<i64>,
    pub polya_len_min: Option<i64>,
    pub polya_len_max: Option<i64>,
    pub beyond_polya: Vec<u8>,
    pub beyond_start: i64,
    pub beyond_polya_support: i64,
}

impl Default for ConsensusResult {
    fn default() -> Self {
        ConsensusResult {
            seq: Vec::new(),
            depth: Vec::new(),
            score: Vec::new(),
            stop_reason: "empty",
            polya_base: None,
            polya_start: -1,
            polya_len_median: None,
            polya_len_min: None,
            polya_len_max: None,
            beyond_polya: Vec::new(),
            beyond_start: -1,
            beyond_polya_support: 0,
        }
    }
}

impl ConsensusResult {
    /// `polya_len_range`: "" when no poly-A, else "{min}-{max}".
    pub fn polya_len_range(&self) -> String {
        match (self.polya_len_min, self.polya_len_max) {
            (Some(a), Some(b)) => format!("{a}-{b}"),
            _ => String::new(),
        }
    }
}

/// `indel_aware_consensus(reads, min_depth, iterations=4, polya_min_len, mate_min_overlap=20)`
/// (indel_consensus.py:336). Reads are consumed in the given order -- the order is significant
/// (seed choice, vote-map insertion order, float summation order).
pub fn indel_aware_consensus(reads: &[ClipRead], min_depth: usize, polya_min_len: usize) -> ConsensusResult {
    todo!("P2: SPEC.md §4.4")
}

/// `rle(seq, qual)` (indel_consensus.py:84): uppercased run bases, run lengths, mean run
/// quality (`sum/len` as f64).
pub fn rle(seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u32>, Vec<f64>) {
    todo!("P2")
}
