//! Evidence-level TPRT filters (far pair, slippage) and their helpers. OWNER: P4.
//!
//! Mirrors `discovery_breakpoints`, `_colonies`, `_far_pair_check`, `_inside_mates`,
//! `_to_one_sided`, `_carries_element`, `_slippage_check`, `_outward_clip`, `_aligned_part` of
//! src/combine_insertions_evidence.py:1294-1452 and 358-363. SPEC.md §4.5, §4.6.

use crate::config::Config;
use crate::evidence::junction::JunctionRecord;
use crate::genome::RefFetch;
use crate::library::Matcher;
use crate::model::{ContigId, FileId, Insertion, Interner, Side};
use rustc_hash::FxHashMap;

/// `discovery_breakpoints(records)`: {(contig, side): [(pos, sample)]} over EVERY accepted
/// per-sample discovery record before intersect, in record order. LEFT when type not in (2, 5)
/// and left_pos is set; RIGHT when type not in (1, 4) and right_pos is set. Sample = files[0].
pub type Breakpoints = FxHashMap<(ContigId, Side), Vec<(i64, FileId)>>;

pub fn discovery_breakpoints(records: &[Insertion]) -> Breakpoints {
    todo!("P4")
}

/// `_aligned_part(ins, side)`: RIGHT -> `str(right_aligned.revcomp()).upper()`; LEFT ->
/// `str(left_aligned).upper()`; "" when absent.
pub fn aligned_part(ins: &Insertion, side: Side) -> Vec<u8> {
    todo!("P4")
}

/// `_outward_clip(rec, ins, with_mates=False)`: with_mates -> rec.consensus, else
/// combined_consensus if its seq is non-empty else consensus; non-empty -> its seq uppercased;
/// else the insertion's clipped seq of rec.side uppercased; else "".
pub fn outward_clip(rec: &JunctionRecord, ins: &Insertion, with_mates: bool) -> Vec<u8> {
    todo!("P4")
}

/// `_inside_mates(rec, min_mapq=20, max_dist=1000)`: needs `rec.rows` (not detached).
pub fn inside_mates(rec: &JunctionRecord, contigs: &Interner) -> Vec<Vec<u8>> {
    todo!("P4")
}

/// `_colonies(ins, rec, breakpoints, tol)`: samples of rec's evidence-role rows plus every
/// discovery breakpoint sample within tol of the insertion's junction. Sorted unique FileIds.
pub fn colonies(ins: &Insertion, rec: &JunctionRecord, bps: Option<&Breakpoints>, tol: i64) -> Vec<u32> {
    todo!("P4")
}

/// `_far_pair_check(ins, recs, cfg, matcher, breakpoints, ref_fetch)` -> None (not a far pair /
/// passes) or Some((reason, polya side)). `recs` = [LEFT, RIGHT] records (both present: the
/// caller checks `open_side is None and len(recs) == 2`). Reasons are the python strings
/// ("no_polarity", "no_element_on_complex_clip", "element_antisense", "element_class_conflict",
/// "colony_mismatch", "polya_side_slippage"); "few_fragments" is never produced (gate dropped).
pub fn far_pair_check(
    ins: &Insertion,
    recs: &[JunctionRecord],
    cfg: &Config,
    matcher: &dyn Matcher,
    bps: Option<&Breakpoints>,
    ref_fetch: Option<&dyn RefFetch>,
    contigs: &Interner,
) -> Option<(String, Option<Side>)> {
    todo!("P4: SPEC.md §4.5")
}

/// `_to_one_sided(ins, real_side)`: drop the other side (clipped/aligned None, mates cleared),
/// LEFT real: right_pos = left_pos, type 4, open_side RIGHT, name `c:{L}-oneside_{L}`;
/// RIGHT real: left_pos = right_pos, type 5, open_side LEFT, name `c:oneside_{R}-{R}`.
pub fn to_one_sided(ins: &mut Insertion, real_side: Side) {
    todo!("P4")
}

/// `_carries_element(rec, ins, matcher)`.
pub fn carries_element(rec: &JunctionRecord, ins: &Insertion, matcher: &dyn Matcher, contigs: &Interner) -> bool {
    todo!("P4")
}

/// `_slippage_check(ins, recs, cfg, matcher, ref_fetch)` -> "" or `slippage:{SIDE}({why})`.
pub fn slippage_check(
    ins: &Insertion,
    recs: &[JunctionRecord],
    cfg: &Config,
    matcher: &dyn Matcher,
    ref_fetch: &dyn RefFetch,
    contigs: &Interner,
) -> String {
    todo!("P4: SPEC.md §4.5")
}
