//! Fragments and the lenient within-sample dedup (independent_clusters). OWNER: P3.
//!
//! Mirrors `Fragment`, `collapse_fragments`, `DedupParams`, `_hp_compress_cut`, `_raw_cut`,
//! `_semi_close`, `_prefix_close`, `_seq_close_from_start`, `_has_long_run`, `_seq_close`,
//! `_junction_pos`, `_dup_seqs`, `_clip_shift`, `_outer_in_clip`, `_mate_unreliable`, `_is_dup`,
//! `_is_dup_oriented`, `_identical`, `independent_clusters` of
//! src/combine_insertions_evidence.py:178-519. SPEC.md §4.2.
//!
//! NOT a gate: the clusters decide the consensus vote weights (1/|cluster|), the consensus depth
//! unit (cluster id = `group`) and the reported `n_independent` / `supported` columns, all of which
//! reach the outputs (SPEC.md §0 "n_independent decision").

use crate::config::Config;
use crate::evidence::row::EvidenceRow;
use crate::model::{InputFile, Interner};

/// All rows of one template (sample, frag) at one junction. Holds indices into the junction's
/// row slice (`rows[i]`), so a Fragment borrows nothing.
#[derive(Clone, Debug)]
pub struct Fragment {
    /// FileId of the sample
    pub file: u32,
    pub frag: Box<str>,
    /// indices of this fragment's rows, in row order
    pub rows: Vec<usize>,
    /// `primary`: among non-MATE rows, min by (role priority, r12) -- first wins on ties
    pub primary: usize,
    /// `mate`: among rows with r12 != primary.r12, min by (0 if MATE else 1, r12) -- first wins
    pub mate: Option<usize>,
}

impl Fragment {
    /// `Fragment.swapped()`: same rows with primary and mate exchanged (None without a mate).
    pub fn swapped(&self) -> Option<Fragment> {
        todo!("P3")
    }
}

/// `collapse_fragments(rows)`: group by (sample, frag) in first-appearance order, build
/// Fragments, drop those without a primary, then STABLE sort by (sample NAME string, primary
/// outer, frag string). `files` gives the sample names.
pub fn collapse_fragments(rows: &[EvidenceRow], files: &[InputFile]) -> Vec<Fragment> {
    todo!("P3")
}

/// python `DedupParams` (`from_cfg`).
#[derive(Clone, Copy, Debug)]
pub struct DedupParams {
    pub tol: i64,
    pub max_edit: i64,
    pub max_edit_frac: f64,
    pub polya_min: usize,
    pub mate_min_mapq: i64,
}

impl DedupParams {
    pub fn from_cfg(cfg: &Config) -> DedupParams {
        DedupParams {
            tol: cfg.dup_coord_tolerance,
            max_edit: cfg.dup_max_edit,
            max_edit_frac: cfg.dup_max_edit_frac,
            polya_min: cfg.polya_min_len,
            mate_min_mapq: cfg.dup_mate_min_mapq,
        }
    }
    /// `budget(n)` -> `crate::pyfmt::dedup_budget`.
    pub fn budget(&self, n: usize) -> i64 {
        crate::pyfmt::dedup_budget(self.max_edit, self.max_edit_frac, n)
    }
}

/// Dedup statistics (`stats` dict: n_dup_coord / n_dup_seq).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DedupStats {
    pub n_dup_coord: u32,
    pub n_dup_seq: u32,
}

/// Result of `independent_clusters`.
#[derive(Clone, Debug, Default)]
pub struct Clusters {
    /// clusters as indices into the `frags` argument; clusters ordered by their smallest member
    /// index, members ascending (python `groups` dict order -- independent of union-find roots)
    pub clusters: Vec<Vec<usize>>,
    pub n_dup: u32,
    pub n_cross: u32,
    pub stats: DedupStats,
}

/// `independent_clusters(frags, params, stats)` (evidence.py:472): within-sample pairs
/// (samples in first-appearance order of `frags`; each sample's indices STABLE-sorted by primary
/// outer; all pairs i<j in that order, skipping already-connected pairs, `_is_dup` -> union,
/// n_dup += 1, stats by kind), then cross-sample exact identity (`by_key[(strand, outer)]` groups
/// in first-appearance order, pairs in index order, different sample, not connected,
/// `_identical` -> union, n_cross += 1).
pub fn independent_clusters(frags: &[Fragment], rows: &[EvidenceRow], p: &DedupParams, contigs: &Interner) -> Clusters {
    todo!("P3: SPEC.md §4.2")
}

/// `_is_dup(a, b, p)` -> "" / "coord" / "seq" (as an enum-free &'static str for parity).
pub fn is_dup(a: &Fragment, b: &Fragment, rows: &[EvidenceRow], p: &DedupParams, contigs: &Interner) -> &'static str {
    todo!("P3")
}

/// `_seq_close(a, b, p, shift)` (exposed for unit tests against the python).
pub fn seq_close(a: &[u8], b: &[u8], p: &DedupParams, shift: i64) -> bool {
    todo!("P3")
}
