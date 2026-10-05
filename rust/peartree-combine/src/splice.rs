//! `<stem>.combined.splice.tsv` (Feature B splice-hallmark re-keying). OWNER: P6.
//!
//! Mirrors `_splice_sidecar` / `write_combined_splice` (combine_insertions.py:32-75). SPEC.md §6.6.

use crate::context::Ctx;
use crate::model::Insertion;
use std::path::PathBuf;

/// `_splice_sidecar(combined)`: strip `.txt.gz` then append `.splice.tsv`.
pub fn splice_sidecar(combined: &str) -> PathBuf {
    PathBuf::from(format!("{}.splice.tsv", combined.strip_suffix(".txt.gz").unwrap_or(combined)))
}

/// Read `<f>.splice.tsv` of EVERY input file (all `--discovery_files`, not only accepted ones);
/// write nothing when none exists. Uncompressed output.
pub fn write_combined_splice(insertions: &[Insertion], combined: &str, ctx: &Ctx, window: i64) -> Result<(), String> {
    todo!("P6: SPEC.md §6.6")
}
