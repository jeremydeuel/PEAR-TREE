//! 2bit genome access. OWNER: P1.
//!
//! Mirrors src/combine_insertions_get_sequence.py (`get_sequence`, py2bit semantics).
//! SPEC.md §2.3.

use std::path::Path;

/// Reference access used by the genotyping output, the slippage / far-pair filters and the
/// SHORT-overhang check (python `ref_fetch`). Must be `Sync`: called from rayon workers.
pub trait RefFetch: Sync {
    /// python `get_sequence(seqname, start, end)` -> UPPERCASE bases, or empty on any failure.
    fn fetch(&self, seqname: &str, start: i64, end: i64) -> Vec<u8>;
}

/// A .2bit file (UCSC format v0, both endiannesses, 32-bit offsets; 64-bit offset variant
/// optional). Reads with positional reads (`FileExt::read_at`) so one handle is shared by all
/// threads; only the header/index (+ per-sequence N-block and length tables, loaded lazily or at
/// open) is held in memory. Soft-mask blocks are ignored (py2bit default `storeMasked=False`
/// returns uppercase); N blocks yield 'N'.
pub struct Genome {
    // P1: file handle, seqnames in FILE INDEX ORDER (python `GENOME.chroms().keys()` order --
    // the substring fallback takes the first match in that order), lengths, offsets, N blocks
}

impl Genome {
    pub fn open(path: &Path) -> Result<Genome, String> {
        todo!("P1")
    }

    /// Sequence names in 2bit index order.
    pub fn seqnames(&self) -> Vec<&str> {
        todo!("P1")
    }
}

impl RefFetch for Genome {
    /// Name mapping (get_sequence.py:37-55): if `seqname` not present: 'MT' -> 'chrM'; else if
    /// `chr{seqname}` present use it; else if len(seqname) > 2: strip a trailing '.1', then the
    /// FIRST name (index order) that contains it as a substring. Still absent -> "".
    /// `start >= end` -> "". py2bit then: `end` is clamped to the sequence length; `start < 0`
    /// or `start >= clamped end` -> error -> "" (verified against py2bit, SPEC.md §2.3).
    fn fetch(&self, seqname: &str, start: i64, end: i64) -> Vec<u8> {
        todo!("P1")
    }
}
