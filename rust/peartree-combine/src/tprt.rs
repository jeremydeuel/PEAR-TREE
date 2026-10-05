//! Combine-level TPRT filter primitives. OWNER: P4.
//!
//! Mirrors src/combine_insertions_tprt_filters.py:45-311 (reference repeats, slippage rule,
//! far-pair verdict). The fuzzy-merge helpers (hp_compress / polya_trimmed / clips_agree) live
//! in `crate::seq` (P1). All sequences OUTWARD from the junction. SPEC.md §4.6.

use crate::config::Config;
use crate::genome::RefFetch;
use crate::library::Matcher;
use crate::model::Side;

/// `outward_reference(fetch, contig, junction, side, width=80)` -> (line, j).
/// lo = max(0, junction - width); s = fetch(contig, lo, junction + width) uppercased;
/// j = junction - lo; LEFT -> (revcomp(s), len(s) - j) else (s, j). NB python returns
/// ("", 0) only on an exception, which `get_sequence` never raises (it returns ""): an empty
/// fetch therefore yields ("", j) / ("", -j) -- callers only test `not line`.
pub fn outward_reference(fetch: &dyn RefFetch, contig: &str, junction: i64, side: Side, width: i64) -> (Vec<u8>, i64) {
    todo!("P4")
}

/// `repeat_at(line, j, max_period=6, slack=2)` -> (start, end, unit).
pub fn repeat_at(line: &[u8], j: i64, max_period: usize, slack: i64) -> (i64, i64, Vec<u8>) {
    todo!("P4")
}

/// `strip_repeat(s, unit, max_mismatch_frac=0.125)`.
pub fn strip_repeat(s: &[u8], unit: &[u8], max_mismatch_frac: f64) -> usize {
    todo!("P4")
}

/// `structured_len(s, polya_min=8)`.
pub fn structured_len(s: &[u8], polya_min: usize) -> usize {
    todo!("P4")
}

/// `slippage_junction(clip, line, j, cfg)` -> "" | "repeat_only" | "repeat_junk" |
/// "repeat_shifted_reference".
pub fn slippage_junction(clip: &[u8], line: &[u8], j: i64, cfg: &Config) -> &'static str {
    todo!("P4: SPEC.md §4.6")
}

/// `leading_polyt(seq, n=10)`.
pub fn leading_polyt(seq: &[u8], n: usize) -> bool {
    todo!("P4")
}

/// `after_polyt(seq)`.
pub fn after_polyt(seq: &[u8]) -> Vec<u8> {
    todo!("P4")
}

/// `far_geometry(gap, cfg)`: gap < -far_pair_max_tsd_deletion or gap > far_pair_tsd_max.
pub fn far_geometry(gap: i64, cfg: &Config) -> bool {
    todo!("P4")
}

/// `colonies_consistent(a, b, frac)`: share a colony and `|a ^ b| <= max(1, int(frac*|a|b|))`.
/// Colonies are sample-name ids (FileId).
pub fn colonies_consistent(a: &[u32], b: &[u32], frac: f64) -> bool {
    todo!("P4")
}

/// Per-side inputs of `far_pair_verdict`. Index 0 = LEFT, 1 = RIGHT (python dict keys; a side
/// missing from the python dict is an empty Vec here, which python's `.get(s, ())` equals).
pub struct FarPairInput {
    /// candidate outward clips, best first, deduplicated, non-empty
    pub clips: [Vec<Vec<u8>>; 2],
    /// colonies (sorted, unique) per side
    pub colonies: [Vec<u32>; 2],
    /// `_inside_mates` per side
    pub inside_mates: [Vec<Vec<u8>>; 2],
}

/// `far_pair_verdict(clips, colonies, matcher, cfg, inside_mates, n_ind=None)` -> (reason, polya
/// side). The `n_ind` / "few_fragments" branch belongs to the dropped pooled gate and is NOT
/// ported (python passes n_ind=None whenever require_independent_fragments is off).
pub fn far_pair_verdict(inp: &FarPairInput, matcher: &dyn Matcher, cfg: &Config) -> (&'static str, Option<Side>) {
    todo!("P4: SPEC.md §4.6")
}
