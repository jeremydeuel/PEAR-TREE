//! Quality-aware affine-gap local aligner with read soft-clipping (SPEC "Read likelihood"). Owner: B.

use crate::config::Config;

/// Result of aligning one read against one segment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Alignment {
    /// natural-log likelihood of the read given the segment (best path)
    pub ll: f64,
    /// read bases aligned (not clipped)
    pub aligned_bases: usize,
    /// segment columns [start, end) covered by the aligned part
    pub seg_start: usize,
    pub seg_end: usize,
}

/// Align `seq`/`qual` (phred) against `seg_seq`/`seg_qual`. `diagonal` = expected segment index
/// of read base 0 (None -> unbanded). Implements the banded DP with full fallback.
pub fn align(seq: &[u8], qual: &[u8], seg_seq: &[u8], seg_qual: &[u8], diagonal: Option<i64>, cfg: &Config) -> Alignment {
    let _ = (seq, qual, seg_seq, seg_qual, diagonal, cfg);
    todo!("owner B: align::align")
}
