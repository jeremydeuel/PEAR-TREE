//! Reference genome access: FASTA (+ `.fai`) or UCSC 2bit, behind one trait. Owner: A.
//! `fetch` returns UPPERCASE bases, `N` for coordinates outside the contig.

use std::io;

pub trait RefSeq: Send {
    /// `genome[start, end)` of `chr`, 0-based half-open, uppercase; positions outside the contig
    /// are filled with `N` so callers never need to clamp. Unknown contig -> Err.
    fn fetch(&mut self, chr: &str, start: i64, end: i64) -> io::Result<Vec<u8>>;
    fn contig_len(&self, chr: &str) -> Option<i64>;
}

/// Open `path` as FASTA (needs `<path>.fai`) or 2bit (by extension `.2bit`).
pub fn open_reference(path: &str) -> io::Result<Box<dyn RefSeq>> {
    let _ = path;
    todo!("owner A: refseq::open_reference")
}
