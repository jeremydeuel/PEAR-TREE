//! Minimal alignment-record decoding for genotyping. Works over any
//! `sam::alignment::Record` (BAM here; CRAM would slot in the same way). Cheap
//! fields (flags, mapq, positions, full CIGAR) decode in `decode_light`; the heavy
//! per-read fields (name, sequence, qualities) decode on demand for the reads that
//! survive the gates, mirroring the Python genotyper which only scores survivors.

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

pub struct LightRec {
    pub reference_sequence_id: Option<usize>,
    pub reference_start: i64, // 0-based
    pub reference_end: i64,   // 0-based, exclusive (pysam reference_end)
    pub mapq: u8,
    pub is_unmapped: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub is_qcfail: bool,
    pub is_duplicate: bool,
    pub cigar: Vec<(u8, usize)>,
}

/// Decode the cheap fields plus the full CIGAR (needed for the breakpoint walk).
pub fn decode_light(record: &dyn AlignmentRecord, header: &Header) -> io::Result<LightRec> {
    let flags: Flags = record.flags()?;
    let mapped = !flags.is_unmapped();
    let reference_sequence_id = record.reference_sequence_id(header).transpose()?;
    let reference_start = record
        .alignment_start()
        .transpose()?
        .map(|p| usize::from(p) as i64 - 1)
        .unwrap_or(-1);

    let mut ref_len: i64 = 0;
    let mut cigar: Vec<(u8, usize)> = Vec::new();
    for op in record.cigar().iter() {
        let op = op?;
        let kind = op.kind();
        let len = op.len();
        cigar.push((kind_code(kind), len));
        if matches!(
            kind,
            Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
        ) {
            ref_len += len as i64;
        }
    }
    let reference_end = reference_start + ref_len;

    let mapq = if mapped {
        record.mapping_quality().transpose()?.map(|m| m.get()).unwrap_or(255)
    } else {
        0
    };

    Ok(LightRec {
        reference_sequence_id,
        reference_start,
        reference_end,
        mapq,
        is_unmapped: flags.is_unmapped(),
        is_secondary: flags.is_secondary(),
        is_supplementary: flags.is_supplementary(),
        is_qcfail: flags.is_qc_fail(),
        is_duplicate: flags.is_duplicate(),
        cigar,
    })
}

/// Read name (qname). Empty string if absent.
pub fn name(record: &dyn AlignmentRecord) -> Vec<u8> {
    record
        .name()
        .map(|n| {
            let b: &[u8] = n.as_ref();
            b.to_vec()
        })
        .unwrap_or_default()
}

/// Read sequence bytes (uppercase ACGTN), forward/reference orientation. Empty for SEQ '*'.
pub fn seq(record: &dyn AlignmentRecord) -> Vec<u8> {
    record.sequence().iter().collect()
}

/// Per-base phred qualities. Empty if absent.
pub fn qual(record: &dyn AlignmentRecord) -> Vec<u8> {
    record.quality_scores().iter().collect::<io::Result<Vec<u8>>>().unwrap_or_default()
}

/// Resolve a contig name from the header by reference id.
pub fn reference_name(header: &Header, id: usize) -> Option<Vec<u8>> {
    header
        .reference_sequences()
        .get_index(id)
        .map(|(name, _)| {
            let b: &[u8] = name.as_ref();
            b.to_vec()
        })
}
