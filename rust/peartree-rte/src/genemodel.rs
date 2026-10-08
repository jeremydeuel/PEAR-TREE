//! The part of annotate_v2.GeneModel (tools/annotate_v2.py:83-300) that the pre-mRNA lookup
//! (annotator._premrna_fn) uses: `_load`, `_resolve`, `_candidates`, `_genic_feature`,
//! `_splice_class` (windows from CONFIG['annotate'] splice_* / promoter_*; `_candidates` reaches
//! `promoter_up`, only the returned FEATURE KEYWORD is used by _premrna_fn).
//!
//! STATUS: WORK PACKAGE "WP-TD" (small; bundled with transduction / pseudogene).
//!
//! Golden: events `premrna` (in: seq, site; out: label|null) recorded from the annotator's
//! pre-mRNA callback, plus `gene_candidates` (in: contig, p; out: [[gs, ge, name, strand], ...]).

use crate::config::GeneModelCfg;
use std::path::Path;

/// One gene: (start, end, name, strand, merged exons).
#[derive(Clone, Debug, PartialEq)]
pub struct Gene {
    pub start: i64,
    pub end: i64,
    pub name: String,
    pub strand: String,
    pub exons: Vec<(i64, i64)>,
}

/// annotate_v2.GeneModel (subset). Windows default like python.
pub struct GeneModel {
    pub donor_window: i64,
    pub acceptor_window: i64,
    pub ppt_window: i64,
    pub branch_window: i64,
    pub prom_up: i64,
    /// contig -> genes sorted by (start, end, name, ...) like python's tuple sort
    pub genes: rustc_hash::FxHashMap<String, Vec<Gene>>,
    pub maxspan: rustc_hash::FxHashMap<String, i64>,
}

impl GeneModel {
    /// `GeneModel(path, cfg)`. WP-TD.
    pub fn open(_path: &Path, _cfg: &GeneModelCfg) -> Result<GeneModel, String> {
        todo!("WP-TD: port annotate_v2.GeneModel._load")
    }

    /// `_candidates(contig, p)`: None when the contig is absent. WP-TD.
    pub fn candidates(&self, _contig: &str, _p: i64) -> Option<Vec<&Gene>> {
        todo!("WP-TD: port annotate_v2.GeneModel._candidates")
    }

    /// `_genic_feature(exons, strand, p)` -> (rank, keyword, label). WP-TD.
    pub fn genic_feature(&self, _exons: &[(i64, i64)], _strand: &str, _p: i64) -> (i64, &'static str, &'static str) {
        todo!("WP-TD: port annotate_v2.GeneModel._genic_feature")
    }
}
