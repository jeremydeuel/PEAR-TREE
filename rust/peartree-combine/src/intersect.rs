//! Pool per-sample discovery records by locus. OWNER: P1.
//!
//! Mirrors src/combine_insertions_intersect_insertions.py (`intersect_insertions` and helpers).
//! SPEC.md §3.3 -- read it for the exact ordering rules: the OUTPUT ORDER of this function is
//! the order of every later output file, and it is python dict insertion order.

use crate::config::Config;
use crate::model::{FileId, Insertion, InputFile, Interner, Side};

/// `intersect_insertions(insertions)` with the three settings taken from `cfg`
/// (`keep_polya_one_sided`, `merge_tolerance_bp`, `polya_aware_clip_agreement`).
///
/// Input: every accepted record, files in command-line order, records in file order.
/// Output: surviving insertions in python `full_insertions.values()` order (SPEC.md §3.3), with
/// `uid` = output index (0..n).
pub fn intersect_insertions(records: Vec<Insertion>, cfg: &Config, contigs: &Interner, files: &[InputFile]) -> Vec<Insertion> {
    todo!("P1: SPEC.md §3.3")
}

/// `_shift_tolerant_agree(a, b, tol, polya_aware)` (intersect:57): None on either side -> true;
/// polya_trimmed (polya_min 8) both when polya_aware else uppercase; shorter = x; len(x) < 6 ->
/// true; `edlib HW distance(x[:20], y[:20+tol+4], k=max(1, len(q)//4)) != -1`.
pub fn shift_tolerant_agree(a: Option<&[u8]>, b: Option<&[u8]>, tol: i64, polya_aware: bool) -> bool {
    todo!("P1")
}

/// `_fuzzy_clusters(keys, weight, tol)` (intersect:92) over full-insertion keys
/// (contig, left_pos, right_pos): order by (-weight, contig NAME string, L, R); bucket width
/// `max(tol, 1)` on L; a key joins the lowest-index cluster whose REPRESENTATIVE has
/// |dL| <= tol and |dR| <= tol among buckets b-1, b, b+1 of the key's L bucket (python floor
/// division -- positions are non-negative), else starts a cluster (registered in its own
/// bucket only). Returns clusters as index lists into `keys`, representative first.
pub fn fuzzy_clusters(keys: &[(u32, i64, i64)], weight: &[usize], tol: i64, contigs: &Interner) -> Vec<Vec<usize>> {
    todo!("P1")
}

/// `_absorb(target, other, sides)` (intersect:116): member loci of `other` (its member_loci +
/// (f, other.name) for f in other.files, deduplicated in order) not yet in target.member_loci are
/// appended and get `member_sides[m] = sides` (target.member_sides starts from a copy of the
/// existing dict or empty); files of `other` not in target.files are appended. (Mates are not
/// stored.)
pub fn absorb(target: &mut Insertion, other: &Insertion, sides: &[Side]) {
    todo!("P1")
}

/// helper used by absorb / evidence: `(f, name)` member for each file id of `ins`.
pub fn own_members(ins: &Insertion) -> Vec<(FileId, crate::model::LocusKey)> {
    ins.files.iter().map(|&f| (f, ins.locus())).collect()
}
