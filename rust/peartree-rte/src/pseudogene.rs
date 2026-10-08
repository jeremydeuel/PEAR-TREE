//! tools/rte/pseudogene.py -- processed-pseudogene proof (exon-exon junction read).
//!
//! STATUS: WORK PACKAGE "WP-TD" (PORT_PLAN.md; small, bundled with transduction).
//!
//! Python functions to port: load_exons_by_gene (27; merged per gene, `sorted()` of
//! (contig, start, end) tuples), load_gene_strands (43), _merge (57), ExonJunctionIndex.mrna (76),
//! .structure (86), .cores (110), .find (129). Caches (`_cores`, `_mrna`) must be thread-safe
//! (the annotator runs loci in parallel): a Mutex'd map or OnceLock per gene.
//!
//! Golden: events `exon_find` (in: genes, seqs; out: [[label, name], ...]), `exon_structure`
//! (in: genes, seqs, tol; out: str|null), with the exon track / remap genome ids.

use crate::config::PseudogeneCfg;
use crate::genome::Genome;
use rustc_hash::FxHashMap;
use std::path::Path;

/// gene -> merged exons [(contig, start, end)] (0-based half-open).
pub type ExonsByGene = FxHashMap<String, Vec<(String, i64, i64)>>;

/// `load_exons_by_gene(path)`. WP-TD.
pub fn load_exons_by_gene(_path: &Path) -> Result<ExonsByGene, String> {
    todo!("WP-TD: port pseudogene.load_exons_by_gene")
}

/// `load_gene_strands(path)` -> gene -> '+'/'-' (first seen). WP-TD.
pub fn load_gene_strands(_path: &Path) -> Result<FxHashMap<String, char>, String> {
    todo!("WP-TD: port pseudogene.load_gene_strands")
}

/// pseudogene.ExonJunctionIndex
pub struct ExonJunctionIndex<'a> {
    pub cfg: PseudogeneCfg,
    pub exons: ExonsByGene,
    pub genome: Option<&'a dyn Genome>,
    pub strands: FxHashMap<String, char>,
}

impl ExonJunctionIndex<'_> {
    /// `structure(genes, five_prime_seqs, tol=15, probe=25)`. WP-TD.
    pub fn structure(&self, _genes: &[String], _five_prime_seqs: &[Vec<u8>], _tol: i64) -> Option<String> {
        todo!("WP-TD: port ExonJunctionIndex.structure")
    }

    /// `find(genes, seqs)` -> [(junction label, name)]. WP-TD.
    pub fn find(&self, _genes: &[String], _seqs: &[(String, Vec<u8>)]) -> Vec<(String, String)> {
        todo!("WP-TD: port ExonJunctionIndex.find")
    }
}
