//! Thin FFI over the vendored edlib 1.2.7 (vendor/edlib, compiled by build.rs).
//!
//! IMPLEMENTED (architect, phase 1) -- every other module calls edlib only through here.
//!
//! Mirrors python-edlib 1.3.9 `edlib.align(query, target, mode, task, k, additionalEqualities)`
//! for ASCII inputs (python-edlib passes ASCII `str` straight through as bytes, so the C call is
//! identical):
//!   * `distance()`  == `edlib.align(q, t, mode=M, task="distance", k=K)["editDistance"]`
//!     (-1 when the distance exceeds `k`; `k = -1` = unlimited).
//!   * `path()`      == `edlib.align(q, t, mode=M, task="path", ...)` -> edit distance,
//!     `locations[0][0]` (start in target of the first optimal location) and the alignment
//!     path as raw edit ops (one per column; run-length encoding them gives exactly the
//!     EXTENDED cigar python-edlib returns: 0 '=' match, 1 'I' insertion to target (query
//!     base consumed), 2 'D' deletion from target (target base consumed), 3 'X' mismatch).
//!
//! SPEC.md "edlib" lists every call site with its mode/task/k/equalities.

use std::os::raw::{c_char, c_int, c_uchar};

#[repr(C)]
#[derive(Clone, Copy)]
struct EdlibEqualityPair {
    first: c_char,
    second: c_char,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct EdlibAlignConfig {
    k: c_int,
    mode: c_int, // EdlibAlignMode: NW=0, SHW=1, HW=2
    task: c_int, // EdlibAlignTask: DISTANCE=0, LOC=1, PATH=2
    additional_equalities: *const EdlibEqualityPair,
    additional_equalities_length: c_int,
}

#[repr(C)]
struct EdlibAlignResult {
    status: c_int,
    edit_distance: c_int,
    end_locations: *mut c_int,
    start_locations: *mut c_int,
    num_locations: c_int,
    alignment: *mut c_uchar,
    alignment_length: c_int,
    alphabet_length: c_int,
}

extern "C" {
    fn edlibAlign(
        query: *const c_char,
        query_length: c_int,
        target: *const c_char,
        target_length: c_int,
        config: EdlibAlignConfig,
    ) -> EdlibAlignResult;
    fn edlibFreeAlignResult(result: EdlibAlignResult);
}

/// edlib alignment mode (python strings "NW" / "SHW" / "HW").
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// global
    Nw,
    /// prefix: the end of the TARGET is free (python "SHW")
    Shw,
    /// infix: both ends of the target are free (python "HW")
    Hw,
}

impl Mode {
    fn c(self) -> c_int {
        match self {
            Mode::Nw => 0,
            Mode::Shw => 1,
            Mode::Hw => 2,
        }
    }
}

/// edit op codes of `PathResult::ops` (edlib EDLIB_EDOP_*)
pub const OP_MATCH: u8 = 0;
/// query base consumed, no target base ('I' in the extended cigar)
pub const OP_INS: u8 = 1;
/// target base consumed, no query base ('D')
pub const OP_DEL: u8 = 2;
pub const OP_MISMATCH: u8 = 3;

/// `[("N","A"),("N","C"),("N","G"),("N","T")]` -- indel_consensus._WILDCARD.
pub const WILDCARD_N: &[(u8, u8)] = &[(b'N', b'A'), (b'N', b'C'), (b'N', b'G'), (b'N', b'T')];

/// Result of a `task="path"` alignment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathResult {
    pub edit_distance: i32,
    /// python `res["locations"][0][0]`: 0-based start of the first optimal alignment in target
    pub start: i32,
    /// one op per alignment column (OP_MATCH / OP_INS / OP_DEL / OP_MISMATCH)
    pub ops: Vec<u8>,
}

fn run(query: &[u8], target: &[u8], mode: Mode, task: c_int, k: i32, eq: &[(u8, u8)]) -> EdlibAlignResult {
    let pairs: Vec<EdlibEqualityPair> = eq
        .iter()
        .map(|&(a, b)| EdlibEqualityPair { first: a as c_char, second: b as c_char })
        .collect();
    let cfg = EdlibAlignConfig {
        k,
        mode: mode.c(),
        task,
        additional_equalities: if pairs.is_empty() { std::ptr::null() } else { pairs.as_ptr() },
        additional_equalities_length: pairs.len() as c_int,
    };
    // SAFETY: query/target are valid for their lengths for the duration of the call; edlib
    // copies nothing past the call; the equality array outlives the call.
    let res = unsafe {
        edlibAlign(
            query.as_ptr() as *const c_char,
            query.len() as c_int,
            target.as_ptr() as *const c_char,
            target.len() as c_int,
            cfg,
        )
    };
    assert_eq!(res.status, 0, "edlib returned an error status");
    res
}

/// `edlib.align(query, target, mode=mode, task="distance", k=k, additionalEqualities=eq)["editDistance"]`.
/// Returns -1 when the best distance is > `k` (k < 0 = unlimited).
pub fn distance(query: &[u8], target: &[u8], mode: Mode, k: i32, eq: &[(u8, u8)]) -> i32 {
    let res = run(query, target, mode, 0, k, eq);
    let d = res.edit_distance;
    // SAFETY: result came from edlibAlign and is freed exactly once.
    unsafe { edlibFreeAlignResult(res) };
    d
}

/// `edlib.align(query, target, mode=mode, task="path", additionalEqualities=eq)` (k = -1).
/// `None` when python would see `editDistance < 0 or cigar is None` (indel_consensus._align).
pub fn path(query: &[u8], target: &[u8], mode: Mode, eq: &[(u8, u8)]) -> Option<PathResult> {
    let res = run(query, target, mode, 2, -1, eq);
    let out = if res.edit_distance < 0 || res.alignment.is_null() || res.num_locations < 1 {
        None
    } else {
        // SAFETY: edlib allocated `alignment_length` ops and `num_locations` start locations.
        let ops = unsafe { std::slice::from_raw_parts(res.alignment, res.alignment_length as usize) }.to_vec();
        let start = if res.start_locations.is_null() { 0 } else { unsafe { *res.start_locations } };
        Some(PathResult { edit_distance: res.edit_distance, start, ops })
    };
    // SAFETY: freed exactly once.
    unsafe { edlibFreeAlignResult(res) };
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn infix_distance() {
        assert_eq!(distance(b"ACGT", b"GGACGTGG", Mode::Hw, -1, &[]), 0);
        assert_eq!(distance(b"ACGT", b"GGACGAGG", Mode::Hw, -1, &[]), 1);
        assert_eq!(distance(b"ACGTACGT", b"TTTTTTTT", Mode::Hw, 1, &[]), -1);
    }

    #[test]
    fn prefix_path_with_wildcards() {
        // reads longer than the seed register as extension over the N pad
        let p = path(b"ACGTA", b"ACGNNNNN", Mode::Shw, WILDCARD_N).unwrap();
        assert_eq!(p.edit_distance, 0);
        assert_eq!(p.start, 0);
        assert_eq!(p.ops, vec![0, 0, 0, 0, 0]);
    }
}
