//! RTE consensus-library matcher (python `LibraryMatcher`, mappy). OWNER: P4.
//!
//! Mirrors src/combine_insertions_tprt_filters.py:171-229. SPEC.md §4.6 "LibraryMatcher".
//! Index options exactly as python: `k=11, w=3, min_cnt=1, min_chain_score=15,
//! min_dp_score=20 (map_opt.min_dp_max), best_n=3`, no preset, alignment on (mappy always sets
//! MM_F_CIGAR). Element index = `<lib>/consensus.fa`; flank indices = `flanks_3p.fa.gz`,
//! `flanks_5p_sva.fa.gz` when present (skipped when missing or empty).

use std::path::Path;

/// `element_class(name)`: "L1" / "ALU" / "SVA" by case-insensitive prefix, else "FLANK".
pub fn element_class(name: &str) -> &'static str {
    let u = name.to_ascii_uppercase();
    for k in ["L1", "ALU", "SVA"] {
        if u.starts_with(k) {
            return k;
        }
    }
    "FLANK"
}

/// python hit tuple `(name, class, strand, mlen, q_st, q_en)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LibHit {
    pub ctg: String,
    /// element_class(ctg) for the element index, "FLANK" for flank indices
    pub class: &'static str,
    /// b'+' = the query as given is sense to the consensus
    pub strand: u8,
    pub mlen: i64,
    pub q_st: i64,
    pub q_en: i64,
}

/// Library matcher interface (a trait so unit tests can stub it). `Sync`: shared by rayon
/// workers.
pub trait Matcher: Sync {
    /// `LibraryMatcher.hit(seq, with_flanks)`: uppercase; len < 20 -> None; over the element
    /// index then (with_flanks) each flank index, in that order, over each index's hits in
    /// minimap2 order: keep the first hit with the strictly largest `mlen >= min_match (18)`.
    /// (The python result cache does not change results; cache or not.)
    fn hit(&self, seq: &[u8], with_flanks: bool) -> Option<LibHit>;
}

/// The minimap2-backed matcher (feature `mappy`).
pub struct LibraryMatcher {
    // P4: element aligner, flank aligners, min_match
}

impl LibraryMatcher {
    /// `LibraryMatcher(library_dir, flanks=True, min_match=18)`. Error when consensus.fa cannot be
    /// indexed (python RuntimeError).
    pub fn open(library_dir: &Path) -> Result<LibraryMatcher, String> {
        todo!("P4")
    }
}

impl Matcher for LibraryMatcher {
    fn hit(&self, seq: &[u8], with_flanks: bool) -> Option<LibHit> {
        todo!("P4")
    }
}
