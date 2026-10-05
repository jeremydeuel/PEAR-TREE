//! Per-junction evaluation and the evidence.tsv row. OWNER: P3.
//!
//! Mirrors `JunctionRecord`, `evaluate_junction`, `_short_overhang_check`, `_mate_inside`,
//! `_clips_start_with_polya` of src/combine_insertions_evidence.py:523-755. SPEC.md §4.3, §6.4.

use crate::config::Config;
use crate::consensus::ConsensusResult;
use crate::evidence::row::EvidenceRow;
use crate::genome::RefFetch;
use crate::model::{InputFile, Interner, Side};

/// `supported` column: "NA" until apply_evidence decides, then 0 / 1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Supported {
    Na,
    No,
    Yes,
}

/// Where a detached record's reads live (shard store; python `fa_ref`). Opaque to everyone
/// but evidence::store.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadsRef {
    pub chunk: u32,
    pub offset: u64,
    pub len: u32,
}

/// python `JunctionRecord` -- one evidence.tsv row (+ its reads for reads.fa.gz).
#[derive(Clone, Debug)]
pub struct JunctionRecord {
    pub insertion_id: String,
    pub side: Side,
    /// kept rows (unused SHORT fragments removed). None once detached to the store.
    pub rows: Option<Vec<EvidenceRow>>,
    pub reads_ref: Option<ReadsRef>,
    pub n_reads: u32,
    pub n_fragments: u32,
    pub n_independent: u32,
    pub n_samples: u32,
    pub n_mates: u32,
    pub n_duplicates: u32,
    pub n_cross: u32,
    pub n_dup_coord: u32,
    pub n_dup_seq: u32,
    /// "" or the comma-joined sorted member locus ids ("." in the TSV when empty)
    pub member_loci: String,
    pub n_short_used: u32,
    pub n_short_rejected: u32,
    pub n_short_mate_inside: u32,
    pub n_independent_no_short: u32,
    /// SHORT rejection reasons (log only); insertion-ordered (reason, count)
    pub short_reasons: Vec<(&'static str, u32)>,
    pub supported: Supported,
    /// junction reads + mates (TSV)
    pub consensus: ConsensusResult,
    /// junction reads only (combined.txt.gz); `ConsensusResult::default()` unless
    /// `indel_aware_consensus`
    pub combined_consensus: ConsensusResult,
    pub polya_end: bool,
    pub fail_reason: String,
    /// reference part in clip_consensus convention, uppercase (`_aligned_part`)
    pub aligned: Vec<u8>,
}

impl JunctionRecord {
    /// `JunctionRecord.tsv()` (evidence.py:546): the 26 EVIDENCE_TSV_COLUMNS joined by "\t" +
    /// "\n". Byte-exact format in SPEC.md §6.4.
    pub fn tsv(&self) -> String {
        todo!("P3: SPEC.md §6.4")
    }

    /// FASTA text of this record's rows (`_render_reads(name, rec)`, evidence.py:1301):
    /// `>{name}|{side}|{role}|{sample}|{frag}|{r12}\n{allele_forward_seq}\n` per row.
    /// Panics if detached (rows None).
    pub fn render_reads(&self, name: &str, files: &[InputFile]) -> String {
        todo!("P3: SPEC.md §6.5")
    }
}

/// `evaluate_junction(insertion_id, side, rows, cfg, ref_fetch)` (evidence.py:569).
/// `rows` = the pooled (and re-anchored) rows, in pooling order; they are moved into the record
/// (minus dropped SHORT fragments). `ref_fetch` is only used with `count_short_overhang`.
/// `rec.aligned` / `rec.member_loci` are left empty (the caller sets them).
pub fn evaluate_junction(
    insertion_id: &str,
    side: Side,
    rows: Vec<EvidenceRow>,
    cfg: &Config,
    ref_fetch: Option<&dyn RefFetch>,
    contigs: &Interner,
    files: &[InputFile],
) -> JunctionRecord {
    todo!("P3: SPEC.md §4.3")
}
