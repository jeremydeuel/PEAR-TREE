//! A decoded alignment record exposing exactly the fields the discovery logic uses.
//!
//! Backed by `&dyn sam::alignment::Record`, so the same view works over a BAM record
//! and a CRAM record (both implement the trait). Cheap fields (flags, mapq, reference
//! id, cigar-derived positions/clips) are decoded eagerly in `from_record`; the heavy
//! fields — sequence, qualities, name, SA/XA tags, contig name — are decoded on demand
//! via accessor methods (SPD-1 lazy decode). Flag semantics mirror pysam.

use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::Record as AlignmentRecord;
use noodles_sam::Header;
use std::io;

pub struct BamRead<'a> {
    record: &'a dyn AlignmentRecord,
    header: &'a Header,
    /// index into the header's reference sequences (SPD-5a: compare by id, resolve
    /// the contig name string only when it changes)
    pub reference_sequence_id: Option<usize>,
    pub mapped: bool,
    pub reference_start: i64, // 0-based
    pub reference_end: i64,   // 0-based, exclusive
    /// mate coordinate (RNEXT/PNEXT) for the SPD-4 mate fetch
    pub mate_ref_id: Option<usize>,
    pub mate_pos: i64, // 0-based, -1 if unset
    /// raw SAM FLAG (sidecar: kept verbatim, dup bit included for diagnostics)
    pub flag: u16,
    /// mate strand (0x20)
    pub mate_is_reverse: bool,
    /// SAM TLEN
    pub tlen: i32,
    /// soft-clip lengths adjacent to the alignment (hard clips skipped), for the
    /// sidecar `outer` coordinate. Unlike `left_len`/`right_len` these look past a
    /// leading/trailing hard clip.
    pub lead_soft: usize,
    pub trail_soft: usize,
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
    pub has_cigar: bool,
    pub left_is_soft: bool,
    pub left_len: usize,
    pub right_is_soft: bool,
    pub right_len: usize,
}

impl<'a> BamRead<'a> {
    #[inline]
    pub fn is_forward(&self) -> bool {
        !self.is_reverse
    }

    /// Decode only the cheap fields. Sequence/qualities/name/tags/contig-name are
    /// left to the accessor methods below. Works for any `sam::alignment::Record`
    /// (BAM or CRAM).
    pub fn from_record(record: &'a dyn AlignmentRecord, header: &'a Header) -> io::Result<BamRead<'a>> {
        let flags: Flags = record.flags()?;
        let mapped = !flags.is_unmapped();

        let reference_sequence_id = record.reference_sequence_id(header).transpose()?;

        let reference_start = record.alignment_start().transpose()?.map(|p| usize::from(p) as i64 - 1).unwrap_or(-1);

        // cigar: first/last op + reference span
        let mut ref_len: i64 = 0;
        let mut first: Option<(Kind, usize)> = None;
        let mut last: Option<(Kind, usize)> = None;
        // soft clips adjacent to the alignment, skipping hard clips (sidecar `outer`)
        let mut lead_soft: usize = 0;
        let mut trail_soft: usize = 0;
        let mut seen_aligned = false;
        for op in record.cigar().iter() {
            let op = op?;
            let kind = op.kind();
            let len = op.len();
            if first.is_none() {
                first = Some((kind, len));
            }
            last = Some((kind, len));
            match kind {
                Kind::SoftClip => {
                    if seen_aligned {
                        trail_soft += len;
                    } else {
                        lead_soft += len;
                    }
                }
                Kind::HardClip => {}
                _ => {
                    seen_aligned = true;
                    trail_soft = 0;
                }
            }
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

        // An unmapped read's MAPQ is meaningless; pysam/samtools/BAM all report 0 for it.
        // noodles-bam returns Some(0) but noodles-cram returns None here, and defaulting
        // None to 255 would let an unmapped read spuriously pass the MAPQ gate on the CRAM
        // path (it would then be dropped for lacking a CIGAR, silently losing polyA
        // evidence). Force 0 for unmapped reads so BAM and CRAM decode identically.
        let mapq = if mapped {
            record.mapping_quality().transpose()?.map(|m| m.get()).unwrap_or(255)
        } else {
            0
        };

        let mate_ref_id = record.mate_reference_sequence_id(header).transpose()?;
        let mate_pos = record.mate_alignment_start().transpose()?.map(|p| usize::from(p) as i64 - 1).unwrap_or(-1);
        let tlen = record.template_length()?;

        Ok(BamRead {
            record,
            header,
            reference_sequence_id,
            mapped,
            reference_start,
            reference_end,
            mate_ref_id,
            mate_pos,
            flag: flags.bits(),
            mate_is_reverse: flags.is_mate_reverse_complemented(),
            tlen,
            lead_soft,
            trail_soft,
            mapq,
            is_read1: flags.is_first_segment(),
            is_read2: flags.is_last_segment(),
            is_reverse: flags.is_reverse_complemented(),
            is_secondary: flags.is_secondary(),
            is_qcfail: flags.is_qc_fail(),
            is_duplicate: flags.is_duplicate(),
            is_supplementary: flags.is_supplementary(),
            is_proper_pair: flags.is_properly_segmented(),
            mate_is_mapped: !flags.is_mate_unmapped(),
            has_cigar,
            left_is_soft,
            left_len,
            right_is_soft,
            right_len,
        })
    }

    /// Unclipped 5' end of this read on the reference (0-based, inclusive): the
    /// leftmost base incl. leading soft clip for a forward read, the rightmost base
    /// incl. trailing soft clip for a reverse read. -1 if unmapped.
    pub fn outer(&self) -> i64 {
        if !self.mapped || self.reference_start < 0 {
            return -1;
        }
        if self.is_reverse {
            self.reference_end - 1 + self.trail_soft as i64
        } else {
            self.reference_start - self.lead_soft as i64
        }
    }

    /// 1 / 2 for first / last segment, 0 if neither (unpaired).
    pub fn r12(&self) -> u8 {
        if self.is_read1 {
            1
        } else if self.is_read2 {
            2
        } else {
            0
        }
    }

    // --- lazy heavy fields (decode on demand) ---

    /// CIGAR as a SAM string (`*` if empty). Sidecar only.
    pub fn cigar_string(&self) -> String {
        let mut s = String::new();
        for op in self.record.cigar().iter() {
            let Ok(op) = op else { break };
            let c = match op.kind() {
                Kind::Match => 'M',
                Kind::Insertion => 'I',
                Kind::Deletion => 'D',
                Kind::Skip => 'N',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                Kind::Pad => 'P',
                Kind::SequenceMatch => '=',
                Kind::SequenceMismatch => 'X',
            };
            s.push_str(&op.len().to_string());
            s.push(c);
        }
        if s.is_empty() {
            s.push('*');
        }
        s
    }

    /// Offset in the stored sequence (soft clips included) of the base aligned at reference
    /// position `refpos` (0-based). Before the alignment: counted back into the leading
    /// soft clip; at/after the alignment end: counted forward into the trailing soft clip
    /// (`refpos == reference_end` gives the first trailing-clip base, the RIGHT-junction
    /// convention); inside a deletion: the next aligned base. Clamped to [0, len].
    pub fn query_offset_at(&self, refpos: i64) -> i64 {
        let len = self.record_len() as i64;
        if refpos < self.reference_start {
            return (self.lead_soft as i64 - (self.reference_start - refpos)).clamp(0, len);
        }
        if refpos >= self.reference_end {
            return (len - self.trail_soft as i64 + (refpos - self.reference_end)).clamp(0, len);
        }
        let (mut q, mut r) = (0i64, self.reference_start);
        for op in self.record.cigar().iter() {
            let Ok(op) = op else { break };
            let n = op.len() as i64;
            match op.kind() {
                Kind::SoftClip | Kind::Insertion => q += n,
                Kind::HardClip | Kind::Pad => {}
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if refpos < r + n {
                        return q + (refpos - r);
                    }
                    q += n;
                    r += n;
                }
                Kind::Deletion | Kind::Skip => {
                    if refpos < r + n {
                        return q;
                    }
                    r += n;
                }
            }
        }
        q.clamp(0, len)
    }

    /// Stored sequence length (no decode).
    pub fn record_len(&self) -> usize {
        self.record.sequence().len()
    }

    /// Raw query-name bytes (no allocation). Used for the qname hash.
    pub fn name_bytes(&self) -> &[u8] {
        self.record.name().map(|n| -> &[u8] { n.as_ref() }).unwrap_or(&[])
    }

    pub fn seq(&self) -> Vec<u8> {
        self.record.sequence().iter().collect()
    }

    pub fn qual(&self) -> Vec<u8> {
        self.record.quality_scores().iter().collect::<io::Result<Vec<u8>>>().unwrap_or_default()
    }

    pub fn query_name(&self) -> String {
        self.record
            .name()
            .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
            .unwrap_or_default()
    }

    /// Resolve the contig name from the header by reference id (SPD-5a). None for
    /// unmapped reads.
    pub fn reference_name(&self) -> Option<String> {
        self.reference_sequence_id.and_then(|id| {
            self.header
                .reference_sequences()
                .get_index(id)
                .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).into_owned())
        })
    }

    pub fn sa(&self) -> Option<String> {
        self.tag_string(Tag::OTHER_ALIGNMENTS)
    }

    /// Placement (reference id, 0-based POS) of this read's PRIMARY record, taken from the
    /// first SA entry (bwa lists the primary first on a supplementary record). None when
    /// the SA tag is absent/unparseable or names a contig not in the header.
    pub fn sa_primary_loc(&self) -> Option<(i32, i64)> {
        let sa = self.sa()?;
        let mut f = sa.split(';').next()?.split(',');
        let (name, pos) = (f.next()?, f.next()?.parse::<i64>().ok()?);
        let id = self.header.reference_sequences().get_index_of(name.as_bytes())?;
        Some((id as i32, pos - 1))
    }

    pub fn xa(&self) -> Option<String> {
        self.tag_string(Tag::from([b'X', b'A']))
    }

    /// Mate mapping quality (`MQ` integer tag, added by `samtools fixmate`). None if the
    /// tag is absent or not an integer. Used by mate-anchored rescue.
    pub fn mq(&self) -> Option<u8> {
        match self.record.data().get(&Tag::from([b'M', b'Q'])) {
            Some(Ok(value)) => value_to_u8(&value),
            _ => None,
        }
    }

    fn tag_string(&self, tag: Tag) -> Option<String> {
        match self.record.data().get(&tag) {
            Some(Ok(value)) => Some(value_to_string(&value)),
            _ => None,
        }
    }
}

fn value_to_string(value: &Value) -> String {
    match value {
        Value::String(s) => String::from_utf8_lossy(s.as_ref()).into_owned(),
        other => format!("{:?}", other),
    }
}

/// Extract an unsigned mapping-quality from any BAM integer aux value, clamped to u8.
fn value_to_u8(value: &Value) -> Option<u8> {
    let n: i64 = match value {
        Value::UInt8(v) => *v as i64,
        Value::Int8(v) => *v as i64,
        Value::UInt16(v) => *v as i64,
        Value::Int16(v) => *v as i64,
        Value::UInt32(v) => *v as i64,
        Value::Int32(v) => *v as i64,
        _ => return None,
    };
    Some(n.clamp(0, 255) as u8)
}
