//! `<stem>.insertions.evidence.tsv.gz` + `<stem>.insertions.reads.fa.gz`. OWNER: P5.
//!
//! Mirrors `write_evidence_outputs` (src/combine_insertions_evidence.py:1307) and
//! `_evidence_paths` (combine_insertions.py:103). SPEC.md §6.4, §6.5.

use crate::context::Ctx;
use crate::evidence::EvidenceState;
use std::path::{Path, PathBuf};

/// `_evidence_paths(combined)`: strip `.combined.txt.gz` (else `.txt.gz`) ->
/// (`{stem}.insertions.evidence.tsv.gz`, `{stem}.insertions.reads.fa.gz`).
pub fn evidence_paths(combined: &str) -> (PathBuf, PathBuf) {
    let stem = combined
        .strip_suffix(".combined.txt.gz")
        .or_else(|| combined.strip_suffix(".txt.gz"))
        .unwrap_or(combined);
    (
        PathBuf::from(format!("{stem}.insertions.evidence.tsv.gz")),
        PathBuf::from(format!("{stem}.insertions.reads.fa.gz")),
    )
}

/// Write the header + one TSV row per record of each name in `names` (in order; a name without
/// records writes nothing; a repeated name is written again), and the matching reads. `names` =
/// surviving insertion names in combined.txt.gz order, then `sorted(failed)` (byte order).
pub fn write_evidence_outputs(state: &EvidenceState, names: &[String], tsv: &Path, fa: &Path, ctx: &Ctx) -> Result<(), String> {
    todo!("P5: SPEC.md §6.4/6.5")
}
