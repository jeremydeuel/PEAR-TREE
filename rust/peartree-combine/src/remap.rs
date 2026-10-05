//! Consensus FASTQ writing, the two bowtie2 remaps and their filters. OWNER: P6.
//!
//! Mirrors combine_insertions.py:78-100 (`_far_flank_trimmed`), 160-261. SPEC.md §5.
//! bowtie2 / samtools stay subprocesses; BAMs are read back as SAM text via `samtools view`.

use crate::context::Ctx;
use crate::liftover::LiftOver;
use crate::model::Insertion;
use crate::seq::QualSeq;
use rustc_hash::FxHashSet;
use std::path::Path;

/// The consensus FASTQ text written to `<stem>.fq.gz` (first remap input) and to
/// `<stem>.combined.txt.gz` (final output): five passes over `insertions` in order --
/// FULL_INFO: `{name}:L` left_consensus then `{name}:R` right_consensus; RIGHT_POLYA: `:L`;
/// LEFT_POLYA: `:R`; RIGHT_DISC: `:L`; LEFT_DISC: `:R`. SPEC.md §6.1.
pub fn consensus_fastq(insertions: &[Insertion], ctx: &Ctx) -> Vec<u8> {
    todo!("P6: SPEC.md §6.1")
}

/// Second remap input (overwrites `<stem>.fq.gz`): per insertion in order,
/// `(lc.fastq("{name}:L") if lc else "") + (rc.fastq("{name}:R") if rc else "")` with
/// (lc, rc) = far-flank-trimmed clips when `trim_far_flank_before_remap`, else the raw
/// left_clipped / right_clipped. SPEC.md §6.2.
pub fn clip_fastq(insertions: &[Insertion], ctx: &Ctx) -> Vec<u8> {
    todo!("P6: SPEC.md §6.2")
}

/// `_far_flank_trimmed(ins, side, probe=20, min_insert=10)`. `side` b'L' / b'R'.
pub fn far_flank_trimmed(ins: &Insertion, side: u8) -> Option<QualSeq> {
    todo!("P6")
}

/// gzip-write `data` to `path` (compression level `level`; python uses 1 for `.fq.gz`, the
/// gzip default 9 for outputs -- only the decompressed bytes are specified).
pub fn write_gz(path: &Path, data: &[u8], level: u32) -> Result<(), String> {
    todo!("P6")
}

/// Run `cmd` through `/bin/sh -c` (python `os.system`), echoing nothing; a non-zero exit is
/// logged, not fatal (python ignores os.system's return value).
pub fn sh(cmd: &str) {
    todo!("P6")
}

/// First remap command string (byte-identical to python, combine_insertions.py:178):
/// `{bowtie2} {fq} -x {index} --end-to-end --sensitive --threads {threads} --qc-filter |
/// {samtools} view -F 4 -b -o {bam}`. Skipped by the caller when `bam` exists.
pub fn bowtie2_end_to_end_cmd(ctx: &Ctx, fq: &str, bam: &str) -> String {
    todo!("P6")
}

/// Second remap command (combine_insertions.py:225): `{bowtie2} {fq} -k 1000 -x {index2}
/// --local --very-fast --threads {threads} --qc-filter | {samtools} view -F 2308 -b -o {bam}`.
pub fn bowtie2_local_cmd(ctx: &Ctx, fq: &str, bam: &str) -> String {
    todo!("P6")
}

/// One SAM record as needed by the filters.
#[derive(Clone, Debug)]
pub struct SamRec {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    /// 0-based leftmost position
    pub pos0: i64,
    /// (len, op) pairs
    pub cigar: Vec<(u32, u8)>,
    /// AS:i tag
    pub as_tag: Option<i64>,
}

/// Stream `samtools view <bam>` and parse each record. SPEC.md §5.1.
pub fn read_bam(ctx: &Ctx, bam: &Path) -> Result<Vec<SamRec>, String> {
    todo!("P6")
}

/// Clean-remap filter (combine_insertions.py:187-201): names (`qname[:-2]`) of mapped records
/// whose longest CIGAR I < clean_remap_max_insertion and AS (missing -> -999) >= clean_remap_min_as.
pub fn clean_remap_names(recs: &[SamRec], ctx: &Ctx) -> FxHashSet<String> {
    todo!("P6: SPEC.md §5.1")
}

/// Clipped-remap filter (combine_insertions.py:227-258): primary mapped records only; parse
/// `qname.split(":")` = (contig, "L-R", side), prefix "chr" to the contig unless it starts with
/// "chr" (and has >= 3 chars); pos = R for side 'R' else L (python int()); lift
/// `reference_start` (forward) / `reference_end` (reverse, exclusive end) of the hit; any result
/// on the (prefixed) contig with |pos - lifted| < 1000 -> name `qname[:-2]`. SPEC.md §5.2.
pub fn clipped_remap_names(recs: &[SamRec], lo: &LiftOver) -> Result<FxHashSet<String>, String> {
    todo!("P6: SPEC.md §5.2")
}
