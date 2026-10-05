//! Pooled per-patient junction evidence (TPRT sidecars). OWNER: P5 (this file, store.rs,
//! output.rs); row.rs / dedup.rs / junction.rs: P3; filters.rs: P4.
//!
//! Mirrors `apply_evidence`, `_make_evaluate`, `_judge_chunk`, `EvidencePool.absorb_one_sided`,
//! `_reanchor`, `_member_loci`, `_absorb_target`, `_replace_clips` of
//! src/combine_insertions_evidence.py:766-1243, 1255-1298. SPEC.md §4 (overview §4.0).
//!
//! The pooled independent-fragment gate (`require_independent_fragments`, fail reason
//! `n_independent<k`) is NOT ported (SPEC.md §0); `supported` is still reported.

pub mod dedup;
pub mod filters;
pub mod junction;
pub mod output;
pub mod row;
pub mod store;

use crate::context::Ctx;
use crate::evidence::filters::Breakpoints;
use crate::evidence::junction::JunctionRecord;
use crate::evidence::row::EvidenceRow;
use crate::evidence::store::Store;
use crate::model::{FileId, Insertion, Member, Side};
use rustc_hash::FxHashMap;
use std::path::Path;

/// State kept from apply_evidence for absorb_one_sided and write_evidence_outputs (python
/// `(kept, records, failed, {"pool", "store", ...})`).
pub struct EvidenceState {
    /// python `records`: insertion name -> uid whose record list is current for that name
    /// (setdefault for kept then failed; absorb overwrites; absorbed names popped)
    pub records: FxHashMap<String, u32>,
    /// python `recmap` (`id(ins)` -> [JunctionRecord]), indexed by uid
    pub recmap: Vec<Option<Vec<JunctionRecord>>>,
    /// python `failed` names (written sorted after the survivors)
    pub failed: Vec<String>,
    /// python `by_id` (`_member_loci(ins)` at apply time; absorb extends it), indexed by uid
    pub members: Vec<Vec<Member>>,
    /// sidecar rows + detached reads
    pub store: Store,
    /// basenames with a sidecar (python `have`)
    pub have: Vec<FileId>,
}

/// `apply_evidence(insertions, accepted_files, cfg, breakpoints, threads, shard_dir)`.
///
/// Returns `(insertions, None)` unchanged when no accepted file has a sidecar (python returns
/// None -> legacy path). Otherwise `(kept, Some(state))`: kept insertions in input order (with
/// far-pair splits / clip replacements applied), failed (far_pair / slippage) recorded in
/// `state.failed`. Chunks are evaluated in parallel (rayon, `ctx.threads`), results merged in
/// insertion order -- output identical for any thread count / chunking. SPEC.md §4.0.
pub fn apply_evidence(
    insertions: Vec<Insertion>,
    accepted: &[FileId],
    ctx: &Ctx,
    breakpoints: Option<&Breakpoints>,
    shard_dir: &Path,
) -> Result<(Vec<Insertion>, Option<EvidenceState>), String> {
    todo!("P5: SPEC.md §4.0")
}

impl EvidenceState {
    /// `EvidencePool.absorb_one_sided(insertions)` (evidence.py:1193), run by the driver after
    /// the clipped-remap filter. Returns (surviving insertions, n_absorbed). SPEC.md §4.7.
    /// Must not be O(one-sided x all): index candidates per (contig) by junction position; the
    /// selection key `(target is one-sided, |d|, target NAME string)` makes the result
    /// independent of scan order.
    pub fn absorb_one_sided(&mut self, insertions: Vec<Insertion>, ctx: &Ctx) -> (Vec<Insertion>, usize) {
        todo!("P5: SPEC.md §4.7")
    }
}

/// `_member_loci(ins)`: member_loci + (f, ins.name) for f in files, deduplicated in order.
pub fn member_loci(ins: &Insertion) -> Vec<Member> {
    todo!("P5")
}

/// `_reanchor(rows, side, junction)` (evidence.py:1257): CLIP/SHORT rows with clip_at >= 0 whose
/// own locus junction (`LocusKey::junction(side)`) differs from `junction` get
/// `clip_at += junction - jm` when the result stays within `0..=len(seq)`; others unchanged.
/// `junction` None -> unchanged.
pub fn reanchor(rows: Vec<EvidenceRow>, side: Side, junction: Option<i64>) -> Vec<EvidenceRow> {
    todo!("P5")
}

/// `_replace_clips(ins, recs)` -> number replaced: for each rec with non-empty
/// combined_consensus.seq and an existing clip of rec.side: skip if shorter than the old clip or
/// equal ignoring case; else the clip becomes QualSeq(seq lowercased if the old clip
/// `is_lower()` else uppercased, combined_consensus.score).
pub fn replace_clips(ins: &mut Insertion, recs: &[JunctionRecord]) -> usize {
    todo!("P5")
}
