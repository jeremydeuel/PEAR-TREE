//! Discovery-file parser and the `+=` merge. OWNER: P1.
//!
//! Mirrors src/combine_insertions_insertion.py (`Insertion.__init__`, `parseFile`,
//! `__iadd__`) and the per-file import loop of combine_insertions.py:116-129. SPEC.md §2.1, §3.1.

use crate::model::{FileId, Insertion, Interner};
use std::path::Path;

/// Per-file import result (combine_insertions.py:119-129).
pub struct FileImport {
    pub file: FileId,
    /// records after the contig filter (`len(contig) < 6 and contig not in ("MT", "chrM")`),
    /// in file order. `uid` is unset (0) here.
    pub records: Vec<Insertion>,
}

/// `Insertion.parseFile(path)` + the contig filter, streaming the gzip FASTQ.
///
/// Format: 4-line FASTQ records titled `@<contig>:<start>-<end>:<SIDE>:<FIELD>` (SIDE =
/// LEFT/RIGHT; FIELD = CLIPPED / ALIGNED / CLIPPED_POLYA / MATE<n> / other). Blank lines
/// (after `strip()`) between records are skipped; seq/plus/qual lines are `strip()`ped; qual =
/// `ord(c) - 33`. Consecutive records with the same (contig, start, end) form one Insertion;
/// a repeated non-consecutive id yields a second Insertion (python behaviour). A repeated FIELD
/// within one id overwrites (python dict). MATE records are validated and DISCARDED (dead data,
/// SPEC.md "Dead data"). Type/side decoding exactly as `Insertion.__init__` (SPEC.md §2.1):
/// RIGHT real iff RIGHT:ALIGNED and RIGHT:CLIPPED present (`right_aligned = ALIGNED.revcomp()`,
/// `right_pos = int(end)`), elif end token `oneside_` (type 4, open_side RIGHT), elif `disc_`
/// (type 4), else poly-A (`right_clipped = RIGHT:CLIPPED_POLYA`, type 1); LEFT analogous
/// (`left_clipped = LEFT:CLIPPED.revcomp()`, poly-A: `LEFT:CLIPPED_POLYA.revcomp().lower()`,
/// type 2 / 5); LEFT decoding runs after RIGHT and overwrites `type` (python order).
/// `files = [file]`, `member_loci = [(file, own locus)]`.
///
/// Errors (python would raise): malformed title, missing '+', a poly-A end without its
/// CLIPPED_POLYA record, unparsable coordinate.
///
/// Memory: no per-record String; sequences are boxed slices. Called in parallel (one rayon task
/// per file) by the driver; `contigs` is the shared interner.
pub fn parse_discovery_file(path: &Path, file: FileId, contigs: &Interner) -> Result<FileImport, String> {
    todo!("P1: SPEC.md §3.1")
}

/// `Insertion.__iadd__` (combine_insertions_insertion.py:130). Only the branch reachable from
/// intersect_insertions is exercised (both FULL_INFO -- poly-A merging is dead code since the
/// poly-A loop `continue`s), but port the whole method: for each side, longer clipped replaces,
/// longer aligned replaces (strictly longer); mates are not stored; then
/// `files += other.files`, `member_loci += other.member_loci`. Name/positions unchanged in the
/// reachable branch.
pub fn merge_into(target: &mut Insertion, other: Insertion) {
    todo!("P1: SPEC.md §3.3 step 2")
}
