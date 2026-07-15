//! A decoded BAM record with exactly the fields the discovery logic uses,
//! decoupled from noodles. Flag semantics mirror pysam.

use noodles_bam::Record;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::Header;
use std::io;

const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_READ1: u16 = 0x40;
const FLAG_READ2: u16 = 0x80;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QCFAIL: u16 = 0x200;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

pub struct BamRead {
    pub reference_name: Option<String>,
    pub mapped: bool,
    pub reference_start: i64, // 0-based
    pub reference_end: i64,   // 0-based, exclusive
    pub mapq: u8,
    pub is_read1: bool,
    pub is_read2: bool,
    pub is_reverse: bool,
    pub is_secondary: bool,
    pub is_qcfail: bool,
    pub is_duplicate: bool,
    pub is_supplementary: bool,
    pub is_proper_pair: bool,
    pub mate_is_mapped: bool,
    pub query_name: String,
    pub seq: Vec<u8>,
    pub qual: Vec<i32>,
    pub has_cigar: bool,
    pub left_is_soft: bool,
    pub left_len: usize,
    pub right_is_soft: bool,
    pub right_len: usize,
    pub sa: Option<String>,
}

impl BamRead {
    #[inline]
    pub fn is_forward(&self) -> bool {
        !self.is_reverse
    }

    pub fn from_record(record: &Record, header: &Header) -> io::Result<BamRead> {
        let flags = u16::from(record.flags());
        let mapped = flags & FLAG_UNMAPPED == 0;

        let refs = header.reference_sequences();
        let rid = record.reference_sequence_id().transpose()?;
        let reference_name = rid
            .and_then(|id| refs.get_index(id))
            .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).into_owned());

        let start1 = record.alignment_start().transpose()?.map(usize::from);
        let reference_start = start1.map(|s| (s - 1) as i64).unwrap_or(-1);

        // cigar: first/last op + reference span
        let mut ref_len: i64 = 0;
        let mut first: Option<(Kind, usize)> = None;
        let mut last: Option<(Kind, usize)> = None;
        for op in record.cigar().iter() {
            let op = op?;
            let kind = op.kind();
            let len = op.len();
            if first.is_none() {
                first = Some((kind, len));
            }
            last = Some((kind, len));
            if matches!(
                kind,
                Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
            ) {
                ref_len += len as i64;
            }
        }
        let has_cigar = first.is_some();
        let reference_end = reference_start + ref_len;
        let (left_is_soft, left_len) = first.map(|(k, l)| (k == Kind::SoftClip, l)).unwrap_or((false, 0));
        let (right_is_soft, right_len) = last.map(|(k, l)| (k == Kind::SoftClip, l)).unwrap_or((false, 0));

        let seq: Vec<u8> = record.sequence().iter().collect();
        let qual: Vec<i32> = record.quality_scores().as_ref().iter().map(|&q| q as i32).collect();

        let query_name = record
            .name()
            .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
            .unwrap_or_default();

        let sa = match record.data().get(&Tag::OTHER_ALIGNMENTS).transpose()? {
            Some(value) => Some(value_to_string(&value)),
            None => None,
        };

        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);

        Ok(BamRead {
            reference_name,
            mapped,
            reference_start,
            reference_end,
            mapq,
            is_read1: flags & FLAG_READ1 != 0,
            is_read2: flags & FLAG_READ2 != 0,
            is_reverse: flags & FLAG_REVERSE != 0,
            is_secondary: flags & FLAG_SECONDARY != 0,
            is_qcfail: flags & FLAG_QCFAIL != 0,
            is_duplicate: flags & FLAG_DUPLICATE != 0,
            is_supplementary: flags & FLAG_SUPPLEMENTARY != 0,
            is_proper_pair: flags & FLAG_PROPER_PAIR != 0,
            mate_is_mapped: flags & FLAG_MATE_UNMAPPED == 0,
            query_name,
            seq,
            qual,
            has_cigar,
            left_is_soft,
            left_len,
            right_is_soft,
            right_len,
            sa,
        })
    }
}

fn value_to_string(value: &noodles_sam::alignment::record::data::field::Value) -> String {
    use noodles_sam::alignment::record::data::field::Value;
    match value {
        Value::String(s) => String::from_utf8_lossy(s.as_ref()).into_owned(),
        other => format!("{:?}", other),
    }
}
