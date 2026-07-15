//! Port of src/polyABreakpoint.py (PolyABreakpoint).

use crate::config::*;
use crate::filters::clean_clipped_seq;
use crate::qseq::QualitySeq;
use crate::read::BamRead;

pub struct PolyABreakpoint {
    pub polya: bool,
    pub qname: String,
    pub is_forward: bool,
    pub is_read1: bool,
    pub reference_name: Option<String>,
    pub breakpoint: Option<i64>,
    pub clipped: Option<QualitySeq>,
    pub clip: i32,
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// (firstInside, firstOutside, polyA_len) — longest homopolymer run of `base`.
fn find_parts(seq: &[u8], polya: bool) -> (Option<i64>, Option<i64>, i64) {
    let base = if polya { b'A' } else { b'T' };
    let mut current_start: i64 = 0;
    let mut current_end: i64 = 0;
    let mut in_run = false;
    let mut polya_len: i64 = 0;
    let mut first_inside: Option<i64> = None;
    let mut first_outside: Option<i64> = None;

    let consider = |cs: i64, ce: i64, pl: &mut i64, fi: &mut Option<i64>, fo: &mut Option<i64>| {
        if ce - cs + 1 > *pl {
            if polya {
                *fi = Some(cs);
                *fo = Some(ce);
            } else {
                *fi = Some(ce);
                *fo = Some(cs);
            }
            *pl = ce - cs + 1;
        }
    };

    for (i, &b) in seq.iter().enumerate() {
        let i = i as i64;
        if b == base {
            if in_run {
                current_end = i;
            } else {
                in_run = true;
                current_start = i;
                current_end = i;
            }
        } else {
            consider(current_start, current_end, &mut polya_len, &mut first_inside, &mut first_outside);
            current_end = 0;
            current_start = 0;
            in_run = false;
        }
    }
    consider(current_start, current_end, &mut polya_len, &mut first_inside, &mut first_outside);
    (first_inside, first_outside, polya_len)
}

impl PolyABreakpoint {
    fn new(read: &BamRead, polya: bool) -> PolyABreakpoint {
        let seq = QualitySeq::new(read.seq.clone(), read.qual.clone());
        let (first_inside, first_outside, polya_len) = find_parts(&read.seq, polya);
        debug_assert!(polya_len >= POLYA_CUTOFF as i64);

        let clip = if read.is_forward() ^ polya { CLIP_RIGHT } else { CLIP_LEFT };
        let mut clipped: Option<QualitySeq> = None;

        if polya {
            if let Some(fi) = first_inside {
                if fi > 6 {
                    clipped = Some(if read.is_read1 ^ read.is_forward() {
                        clean_clipped_seq(&seq.pyslice(None, Some(fi as isize)).revcomp())
                    } else {
                        clean_clipped_seq(&seq.pyslice(None, Some(fi as isize))).revcomp()
                    });
                }
            }
        } else {
            let fo = first_outside.unwrap_or(0);
            clipped = Some(if read.is_read1 ^ read.is_forward() {
                clean_clipped_seq(&seq.pyslice(Some((fo + 1) as isize), None).revcomp()).revcomp()
            } else {
                clean_clipped_seq(&seq.pyslice(Some((fo + 1) as isize), None))
            });
        }

        PolyABreakpoint {
            polya,
            qname: read.query_name.clone(),
            is_forward: read.is_forward(),
            is_read1: read.is_read1,
            reference_name: None,
            breakpoint: None,
            clipped,
            clip,
        }
    }

    pub fn has_clipped(&self) -> bool {
        matches!(&self.clipped, Some(c) if !c.is_empty())
    }

    /// Port of PolyABreakpoint.setMate.
    pub fn set_mate(&mut self, mate: &BamRead) {
        self.reference_name = mate.reference_name.clone();
        if mate.mapq < MIN_MAPQ {
            return;
        }
        if !mate.mapped {
            self.breakpoint = None;
            self.reference_name = None;
            return;
        }
        self.breakpoint = Some(if self.is_forward ^ self.polya ^ mate.is_forward() {
            mate.reference_start
        } else {
            mate.reference_end
        });
    }

    /// Port of PolyABreakpoint.findPolyA. Returns a PolyABreakpoint with a valid
    /// clipped sequence, or None.
    pub fn find_polya(read: &BamRead) -> Option<PolyABreakpoint> {
        if read.is_proper_pair {
            return None;
        }
        if !read.mate_is_mapped {
            return None;
        }
        let polya_seq = [b'A'; POLYA_CUTOFF];
        let polyt_seq = [b'T'; POLYA_CUTOFF];
        let b = if contains(&read.seq, &polya_seq) {
            if read.is_reverse {
                if read.is_read2 { Some(PolyABreakpoint::new(read, true)) } else { None }
            } else if read.is_read1 {
                Some(PolyABreakpoint::new(read, true))
            } else {
                None
            }
        } else if contains(&read.seq, &polyt_seq) {
            if read.is_reverse {
                if read.is_read1 { Some(PolyABreakpoint::new(read, false)) } else { None }
            } else if read.is_read2 {
                Some(PolyABreakpoint::new(read, false))
            } else {
                None
            }
        } else {
            None
        };
        let b = b?;
        if !b.has_clipped() {
            return None;
        }
        Some(b)
    }
}
