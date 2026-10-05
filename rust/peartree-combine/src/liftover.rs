//! UCSC chain-file liftover with pyliftover 0.4.1 semantics. OWNER: P6.
//!
//! SPEC.md §5.3. pyliftover: `LiftOver(path)` opens gzip iff the path ends with `.gz`
//! (case-insensitive); each `chain` header (score, tName, tSize, tStrand(+ required), tStart,
//! tEnd, qName, qSize, qStrand, qStart, qEnd[, id]) is followed by `size dt dq` lines and a final
//! `size` line; blocks `(sfrom, sfrom+size, tfrom)`. Lines starting with '#', '\n', '\r' between
//! chains are skipped. `convert_coordinate(chrom, pos)` (strand '+'):
//!   * None if `chrom` has no chain;
//!   * else for every block with `start <= pos < end`: `tpos = tfrom + (pos - start)`, if the
//!     chain's target strand is '-': `tpos = target_size - 1 - tpos`; result
//!     `(target_name, tpos, target_strand, score)`, sorted by score descending (stable).
//! combine only asks "does ANY result on the breakpoint's contig lie within 1000 bp", so result
//! order is irrelevant there.

use std::path::Path;

pub struct LiftOver {
    // P6: per source contig: blocks sorted by start (+ an interval index -- blocks of different
    // chains may overlap), chain table (target name, size, strand, score)
}

/// One `convert_coordinate` result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Lifted<'a> {
    pub chrom: &'a str,
    pub pos: i64,
    pub strand: u8,
    pub score: i64,
}

impl LiftOver {
    pub fn open(path: &Path) -> Result<LiftOver, String> {
        todo!("P6")
    }

    pub fn convert_coordinate(&self, chrom: &str, pos: i64) -> Option<Vec<Lifted<'_>>> {
        todo!("P6")
    }
}
