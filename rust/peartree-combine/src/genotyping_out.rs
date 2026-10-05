//! `<stem>.genotyping.txt.gz` (the genotyping contract). OWNER: P6.
//!
//! Mirrors combine_insertions.py:298-365. SPEC.md §6.3.

use crate::context::Ctx;
use crate::model::Insertion;
use rustc_hash::FxHashSet;

/// Returns the uncompressed file text (and logs `wrote N of M, excluded K`).
/// `filter_names` = the clipped-remap filter set (python `filter_reads`, re-checked here).
pub fn genotyping_text(insertions: &[Insertion], filter_names: &FxHashSet<String>, ctx: &Ctx) -> Vec<u8> {
    todo!("P6: SPEC.md §6.3")
}
