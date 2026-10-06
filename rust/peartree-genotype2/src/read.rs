//! Minimal alignment-record decoding for genotyping. Works over any
//! `sam::alignment::Record` (BAM and CRAM). Cheap fields (flags, mapq, positions, mate
//! fields, full CIGAR) decode in `decode_light_into`; the heavy per-read fields (sequence,
//! qualities) decode on demand only for the reads that survive the gates.

use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::Record as AlignmentRecord;
use noodles_sam::Header;
use std::io;

/// noodles CIGAR `Kind` -> BAM op code (matches pysam's integer ops).
#[inline]
fn kind_code(k: Kind) -> u8 {
    match k {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}

/// The cheap, always-decoded part of a record.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LightRec {
    pub reference_sequence_id: Option<usize>,
    pub reference_start: i64, // 0-based, -1 if absent
    pub reference_end: i64,   // 0-based, exclusive (pysam reference_end)
    pub mapq: u8,
    pub is_paired: bool,
    pub is_unmapped: bool,
    pub is_reverse: bool,
    pub is_first: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub is_qcfail: bool,
    pub is_duplicate: bool,
    pub mate_unmapped: bool,
    pub mate_reverse: bool,
    pub mate_reference_sequence_id: Option<usize>,
    pub mate_start: i64, // 0-based, -1 if absent
    pub tlen: i64,
    /// BAM CIGAR as (op code, len)
    pub cigar: Vec<(u8, usize)>,
}

impl LightRec {
    /// Leading soft-clip length (after any hard clip).
    pub fn leading_softclip(&self) -> usize {
        self.cigar.iter().find(|&&(op, _)| op != 5).filter(|&&(op, _)| op == 4).map(|&(_, l)| l).unwrap_or(0)
    }
    /// Trailing soft-clip length (before any hard clip).
    pub fn trailing_softclip(&self) -> usize {
        self.cigar.iter().rev().find(|&&(op, _)| op != 5).filter(|&&(op, _)| op == 4).map(|&(_, l)| l).unwrap_or(0)
    }
}

/// Decode the cheap fields plus the full CIGAR into `out` (its CIGAR buffer is reused, so a
/// stream of records allocates nothing per record).
pub fn decode_light_into(record: &dyn AlignmentRecord, header: &Header, out: &mut LightRec) -> io::Result<()> {
    let flags: Flags = record.flags()?;
    out.reference_sequence_id = record.reference_sequence_id(header).transpose()?;
    out.reference_start = record.alignment_start().transpose()?.map(|p| usize::from(p) as i64 - 1).unwrap_or(-1);

    let mut ref_len: i64 = 0;
    out.cigar.clear();
    for op in record.cigar().iter() {
        let op = op?;
        let kind = op.kind();
        let len = op.len();
        out.cigar.push((kind_code(kind), len));
        if matches!(kind, Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch) {
            ref_len += len as i64;
        }
    }
    out.reference_end = out.reference_start + ref_len;

    out.is_unmapped = flags.is_unmapped();
    out.mapq = if !out.is_unmapped { record.mapping_quality().transpose()?.map(|m| m.get()).unwrap_or(255) } else { 0 };
    out.is_paired = flags.is_segmented();
    out.is_reverse = flags.is_reverse_complemented();
    out.is_first = flags.is_first_segment();
    out.is_secondary = flags.is_secondary();
    out.is_supplementary = flags.is_supplementary();
    out.is_qcfail = flags.is_qc_fail();
    out.is_duplicate = flags.is_duplicate();
    out.mate_unmapped = flags.is_mate_unmapped();
    out.mate_reverse = flags.is_mate_reverse_complemented();
    out.mate_reference_sequence_id = record.mate_reference_sequence_id(header).transpose()?;
    out.mate_start = record.mate_alignment_start().transpose()?.map(|p| usize::from(p) as i64 - 1).unwrap_or(-1);
    out.tlen = record.template_length()? as i64;
    Ok(())
}

/// Read name (qname) bytes, empty if absent. Borrowed: no allocation.
#[inline]
pub fn name(record: &dyn AlignmentRecord) -> &[u8] {
    record.name().map(|n| -> &[u8] { n.as_ref() }).unwrap_or(&[])
}

/// Read sequence bytes (uppercase ACGTN), forward/reference orientation, into `out`.
/// Empty for SEQ '*'.
pub fn seq_into(record: &dyn AlignmentRecord, out: &mut Vec<u8>) {
    out.clear();
    out.extend(record.sequence().iter());
}

/// Per-base phred qualities (0-93) into `out`, as stored (BAM and CRAM `RecordBuf` both hold raw
/// phred, not +33). Left empty when absent — including BAM's "missing" encoding (first byte 0xFF).
pub fn qual_into(record: &dyn AlignmentRecord, out: &mut Vec<u8>) -> io::Result<()> {
    out.clear();
    for q in record.quality_scores().iter() {
        out.push(q?);
    }
    if out.first() == Some(&0xFF) {
        out.clear();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn softclip_lengths_skip_hardclips() {
        let r = LightRec { cigar: vec![(5, 3), (4, 25), (0, 100), (4, 7), (5, 2)], ..Default::default() };
        assert_eq!(r.leading_softclip(), 25);
        assert_eq!(r.trailing_softclip(), 7);
        let r = LightRec { cigar: vec![(0, 150)], ..Default::default() };
        assert_eq!((r.leading_softclip(), r.trailing_softclip()), (0, 0));
    }
}
