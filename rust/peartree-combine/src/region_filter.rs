//! Dense-region filter. OWNER: P1.
//!
//! Mirrors src/combine_insertions_region_filter.py `filter_dense_regions` (called with
//! bin_range=100, ins_cutoff=4 from combine_insertions.py:139-141). SPEC.md §3.4.

use crate::model::Insertion;

/// Returns (kept in input order, n_removed, n_hot_bins).
///
/// pos = right_pos if present else left_pos. bin = `int(pos / bin_range) * bin_range` (python
/// true division then truncation; positions are >= 0). Hot bins: (contig, bin) with count >=
/// cutoff. An insertion is removed iff some hot bin b on its contig satisfies
/// `pos > b - bin_range/2 and pos < b + bin_range*1.5` (float comparisons; with bin_range=100
/// these are exact integers 50 / 150). Candidate bins: b from
/// `((pos - int(bin_range*1.5)) // bin_range) * bin_range` up to (exclusive)
/// `pos + int(bin_range/2) + bin_range`, step bin_range.
pub fn filter_dense_regions(insertions: Vec<Insertion>, bin_range: i64, ins_cutoff: usize) -> (Vec<Insertion>, usize, usize) {
    todo!("P1: SPEC.md §3.4")
}
