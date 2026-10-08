//! tools/rte/annotator.py -- orchestration (RteAnnotator) + the cross-locus passes.
//!
//! STATUS: InsertionInput + the small helpers = FOUNDATION; `Resources::load`, `RteAnnotator::
//! annotate` / `annotate_key` / `annotate_all`, `_premrna_fn`, `_pseudogene_structure_fn`,
//! `_supported_tails`, `_beyond_polya`, `cohort_source_pass`, `recurrence_pass`, `_score` =
//! INTEGRATION stage (PORT_PLAN.md "WP-INT"), after the module packages land.
//!
//! Per-insertion data flow (python `_annotate_key` -> `annotate`):
//!   LocusData (stream.rs) -> pool combine reads + GT reads -> cap_reads(rte_max_reads)
//!   -> split_junction / polya_info / locate_site / target_site (hallmarks)
//!   -> SiteContext windows (genome) -> Assembler.assemble -> supported tails / polya_reads
//!   -> en_motif / slippage_context -> ExonJunctionIndex.find (if candidate genes)
//!   -> structure.classify (NovelSourceFinder, pre-mRNA fn, pseudogene structure fn)
//!   -> site-level tags -> beyond-poly-A -> ScoreInput -> score -> RteRecord
//!   (+ a combine-reads-only rerun for `gt_changed` when the locus has GT reads and
//!   rte_gt_compare).
//! Cross-locus (python `annotate_all`, needs every record first; records are small, kept in
//! memory in input order):
//!   1. cohort_source_pass: only when genome_2bit == remap_2bit and a locator exists; collects
//!      L1 TPRT/LIKELY_TPRT sites, then RE-ANNOTATES (re-loads the locus reads from the store)
//!      every TD3P record without a TD3P_SOURCE= tag with cohort_l1 set; an error keeps the old
//!      record.
//!   2. recurrence_pass: signature (element, structure, int(j5)//5) over L1/ALU/SVA records not
//!      FULL_LENGTH / 5P_UNRESOLVED with detail j5 (or fwd_start); count > rte_recurrence_max ->
//!      score_input.recurrent = True, detail["recurrence"] = count, re-score.
//!
//! Errors: any failure annotating one locus -> RteRecord(key) with detail["error"] =
//! f"{type(e).__name__}: {e}"[:200] (python never lets one locus kill the table; in Rust:
//! catch_unwind around the locus or Result plumbing).

use crate::config::RteConfig;
use crate::record::RteRecord;
use crate::stream::{LocusData, LocusLoader};

/// annotator.InsertionInput (what annotate_v2 hands over per insertion; SPEC.md "Input").
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InsertionInput {
    /// locus name (= insertion id)
    pub title: String,
    /// combined.txt.gz LEFT junction string (lower = clip, UPPER = reference)
    pub left_seq: Vec<u8>,
    pub right_seq: Vec<u8>,
    /// pseudogene candidate genes (InsertionInput.from_legacy order)
    pub pseudogene_genes: Vec<String>,
    /// annotate_v2 element_class(ins.conclusion())
    pub legacy_class: Option<String>,
    /// annotate_v2 `_sv_subtype(allow_rte=True)` rank (only `sv[0]` is read: rank 2 =
    /// intrachromosomal partner)
    pub sv_rank: Option<i64>,
}

/// `_richer_junction(evidence, combined)`: the string with more lower-case (clipped) bases; the
/// evidence consensus on a tie / when combined is empty.
pub fn richer_junction(evidence: &[u8], combined: &[u8]) -> Vec<u8> {
    if evidence.is_empty() {
        return combined.to_vec();
    }
    let lc = |s: &[u8]| s.iter().filter(|c| c.is_ascii_lowercase()).count();
    if !combined.is_empty() && lc(combined) > lc(evidence) {
        combined.to_vec()
    } else {
        evidence.to_vec()
    }
}

/// `default_sidecars(insertions_file)` -> (evidence.tsv.gz, reads.fa.gz)
pub fn default_sidecars(insertions_file: &str) -> (String, String) {
    let mut base = insertions_file;
    for suf in [".combined.txt.gz", ".txt.gz"] {
        if let Some(b) = base.strip_suffix(suf) {
            base = b;
            break;
        }
    }
    (format!("{base}.insertions.evidence.tsv.gz"), format!("{base}.insertions.reads.fa.gz"))
}

/// `default_gt_reads(insertions_file)`
pub fn default_gt_reads(insertions_file: &str) -> String {
    default_sidecars(insertions_file).1.replace(".insertions.reads.fa.gz", ".insertions.genotype_reads.fa.gz")
}

/// Everything loaded once per run (library, genomes, rmsk, locator, exon index, gene model).
/// WP-INT: owns what RteAnnotator borrows.
pub struct Resources {
    pub lib: crate::library::RteLibrary,
}

impl Resources {
    /// Load every configured resource (python RteAnnotator.__init__). WP-INT.
    pub fn load(_cfg: &RteConfig) -> Result<Resources, String> {
        todo!("WP-INT: port RteAnnotator.__init__ (library, genomes, NovelSourceFinder, exon index, gene model)")
    }
}

/// annotator.RteAnnotator
pub struct RteAnnotator<'a> {
    pub cfg: &'a RteConfig,
    pub res: &'a Resources,
}

impl<'a> RteAnnotator<'a> {
    pub fn new(cfg: &'a RteConfig, res: &'a Resources) -> RteAnnotator<'a> {
        RteAnnotator { cfg, res }
    }

    /// `annotate(inp, ev)`: one insertion from its (already pooled) evidence. WP-INT.
    pub fn annotate(&self, _inp: &InsertionInput, _ev: &crate::inputs::InsertionEvidence) -> RteRecord {
        todo!("WP-INT: port RteAnnotator.annotate")
    }

    /// `_annotate_key(key, inp)` with the locus's data. WP-INT.
    pub fn annotate_key(&self, _inp: &InsertionInput, _data: LocusData) -> RteRecord {
        todo!("WP-INT: port RteAnnotator._annotate_key / annotate")
    }

    /// `annotate_all(inputs)`: per-locus annotation (parallel, bounded, input order) + the
    /// cohort-source and recurrence passes. Returns one record per input, in input order. WP-INT.
    pub fn annotate_all(&self, _inputs: &[InsertionInput], _loader: &LocusLoader) -> Vec<RteRecord> {
        todo!("WP-INT: port RteAnnotator.annotate_all / cohort_source_pass / recurrence_pass")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers_like_python() {
        assert_eq!(richer_junction(b"", b"aaGG"), b"aaGG");
        assert_eq!(richer_junction(b"aGG", b"aaGG"), b"aaGG");
        assert_eq!(richer_junction(b"aaGG", b"ccGG"), b"aaGG");
        assert_eq!(richer_junction(b"aaGG", b""), b"aaGG");
        assert_eq!(
            default_sidecars("/x/P1.combined.txt.gz"),
            ("/x/P1.insertions.evidence.tsv.gz".to_string(), "/x/P1.insertions.reads.fa.gz".to_string())
        );
        assert_eq!(default_gt_reads("P1.txt.gz"), "P1.insertions.genotype_reads.fa.gz");
    }
}
