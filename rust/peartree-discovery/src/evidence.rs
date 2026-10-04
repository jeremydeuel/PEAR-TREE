//! TPRT-hallmark evidence sidecar (`<out>.evidence.tsv.gz`, see
//! plans/tprt_hallmarks/SPEC.md). One row per read record supporting an emitted
//! locus: the junction-clipped reads (CLIP), poly-A reads placed by their mate
//! (POLYA), discordant anchors near the junction (DISC) and the mates of those
//! (MATE). Everything here is inert unless `evidence_sidecar` is on, so the main
//! `.txt.gz` output never depends on it.
//!
//! Memory: a record is kept as an `EvRec` (raw seq/qual/CIGAR + a few ints) only for
//! breakpoints that can still be emitted; fragments are identified by a 64-bit FNV-1a
//! hash of the qname (`frag`) plus the r12 bit, never by the qname string.

use std::collections::{BTreeSet, VecDeque};
use std::io::{self, Write};

use crate::read::BamRead;

/// FNV-1a 64-bit hash of the query name. Stable across runs/platforms/versions (the
/// sidecar `frag` column must link a read to its mate across files written by the
/// same binary, and combine compares it within a sample).
#[inline]
pub fn frag_hash(name: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in name {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    Clip,
    PolyA,
    Mate,
    Disc,
    /// short-overhang read: crosses the junction by only 1..`short_overhang_max` bases
    /// (soft clip < MIN_CLIP_LEN, or aligned through) — collected, never judged here
    Short,
}

impl Role {
    pub fn as_str(self) -> &'static str {
        match self {
            Role::Clip => "CLIP",
            Role::PolyA => "POLYA",
            Role::Mate => "MATE",
            Role::Disc => "DISC",
            Role::Short => "SHORT",
        }
    }
}

/// One read record as written to the sidecar. `ref_id`/`mref` are header indices
/// (-1 = unmapped / unset); names are resolved at output time.
#[derive(Clone, Debug)]
pub struct EvRec {
    pub role: Role,
    pub frag: u64,
    pub flag: u16,
    pub ref_id: i32,
    pub pos: i64,
    pub outer: i64,
    pub mref: i32,
    pub mpos: i64,
    pub tlen: i32,
    pub mapq: u8,
    pub clip_at: i32,
    pub cigar: Box<str>,
    /// read bases as stored in the BAM (= reference-forward orientation)
    pub seq: Box<[u8]>,
    /// raw phred scores (not +33); empty / 0xFF-filled when absent
    pub qual: Box<[u8]>,
}

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_READ1: u16 = 0x40;
pub const FLAG_READ2: u16 = 0x80;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

impl EvRec {
    pub fn from_read(read: &BamRead<'_>, role: Role, clip_at: i32) -> EvRec {
        let mapped = read.mapped && read.reference_start >= 0;
        let mate_mapped = read.flag & FLAG_PAIRED != 0 && read.mate_is_mapped && read.mate_pos >= 0;
        EvRec {
            role,
            frag: frag_hash(read.name_bytes()),
            flag: read.flag,
            ref_id: if mapped { read.reference_sequence_id.map_or(-1, |i| i as i32) } else { -1 },
            pos: if mapped { read.reference_start } else { -1 },
            outer: read.outer(),
            mref: if mate_mapped { read.mate_ref_id.map_or(-1, |i| i as i32) } else { -1 },
            mpos: if mate_mapped { read.mate_pos } else { -1 },
            tlen: read.tlen,
            mapq: read.mapq,
            clip_at,
            cigar: read.cigar_string().into_boxed_str(),
            seq: read.seq().into_boxed_slice(),
            qual: read.qual().into_boxed_slice(),
        }
    }

    #[inline]
    pub fn r12(&self) -> u8 {
        r12_of(self.flag)
    }

}

#[inline]
pub fn r12_of(flag: u16) -> u8 {
    if flag & FLAG_READ1 != 0 {
        1
    } else if flag & FLAG_READ2 != 0 {
        2
    } else {
        0
    }
}

/// Compact discordant-anchor observation collected during the extract scan (no
/// sequence). Promoted to a DISC `EvRec` in the mate pass if it lies next to a
/// breakpoint that can still be emitted.
#[derive(Clone, Copy, Debug)]
pub struct DiscLite {
    pub frag: u64,
    pub flag: u16,
    pub ref_id: i32,
    pub start: i64,
    pub end: i64,
}

impl DiscLite {
    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.flag & FLAG_REVERSE != 0
    }
    #[inline]
    pub fn r12(&self) -> u8 {
        r12_of(self.flag)
    }
}

/// Compact identity of one junction-clipped read of a cluster (no sequence): enough to
/// re-find its record in the mate pass — qname hash, raw flag (r12, strand,
/// supplementary) and its breakpoint coordinate (alignment start for a LEFT clip, end
/// for a RIGHT clip).
#[derive(Clone, Copy, Debug)]
pub struct ClipLite {
    pub frag: u64,
    pub flag: u16,
    pub pos: i64,
}

impl ClipLite {
    #[inline]
    pub fn r12(&self) -> u8 {
        r12_of(self.flag)
    }
    #[inline]
    pub fn is_primary(&self) -> bool {
        self.flag & (FLAG_SUPPLEMENTARY | 0x100) == 0
    }
    #[inline]
    pub fn mate_r12(&self) -> u8 {
        match self.r12() {
            1 => 2,
            2 => 1,
            _ => 0,
        }
    }
}

/// Per-breakpoint sidecar payload. Boxed behind an `Option` on `Breakpoint` so the
/// default (sidecar off) path pays 8 bytes per breakpoint. The `*_lite` lists are
/// filled at extract/join time; the full records are captured in the mate pass, and
/// only for breakpoints that can still be emitted.
#[derive(Clone, Debug, Default)]
pub struct EvExtra {
    /// the cluster's clipped reads (capped at `max_evidence_reads_per_breakpoint`)
    pub clip_lite: Vec<ClipLite>,
    /// discordant anchors selected for this breakpoint
    pub disc_lite: Vec<DiscLite>,
    /// short-overhang reads selected for this breakpoint (`short_overhang_evidence`)
    pub short_lite: Vec<ShortLite>,
    pub clip: Vec<EvRec>,
    pub disc: Vec<EvRec>,
    pub short: Vec<EvRec>,
    pub mates: Vec<EvRec>,
}

/// Compact identity of a primary high-MAPQ read that may cross a junction by only a few
/// bases (short-overhang candidate), collected during the extract scan. Soft-clip lengths
/// are the ones adjacent to the alignment (hard clips skipped).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShortLite {
    pub frag: u64,
    pub flag: u16,
    pub ref_id: i32,
    pub start: i64,
    pub end: i64,
    pub lead_soft: u32,
    pub trail_soft: u32,
}

impl ShortLite {
    #[inline]
    pub fn r12(&self) -> u8 {
        r12_of(self.flag)
    }
    #[inline]
    pub fn mate_r12(&self) -> u8 {
        match self.r12() {
            1 => 2,
            2 => 1,
            _ => 0,
        }
    }
}

/// Is `r` a short-overhang read for a junction at `b`? LEFT junction (reads start at b,
/// clip to the left): (a) a leading soft clip of 1..min_clip-1 bases with the alignment
/// starting within ±`window` of b, or (b) no leading soft clip and the alignment starting
/// 1..`max` bases before b while still ending past b. RIGHT junction mirrored on the
/// alignment end. Strand-agnostic (the junction side is a reference property).
pub fn short_candidate(r: &ShortLite, b: i64, left: bool, window: i64, max: i64, min_clip: u32) -> bool {
    if left {
        let a = r.lead_soft >= 1 && r.lead_soft < min_clip && (r.start - b).abs() <= window;
        let bb = r.lead_soft == 0 && r.start < b && b - r.start <= max && r.end > b;
        a || bb
    } else {
        let a = r.trail_soft >= 1 && r.trail_soft < min_clip && (r.end - b).abs() <= window;
        let bb = r.trail_soft == 0 && r.end > b && r.end - b <= max && r.start < b;
        a || bb
    }
}

/// Extract-scan collector for short-overhang candidates (one contig at a time). Every
/// eligible read enters a start-ordered look-back buffer; a junction-clipped read (clip
/// >= MIN_CLIP_LEN) marks its breakpoint "hot". A buffered read is kept only if a hot
/// position of the matching side lies within its candidate window when it leaves the
/// buffer (all clipped reads that could define such a window have been seen by then:
/// `horizon` > the longest aligned span). Memory: O(reads in `horizon` bp + kept).
#[derive(Default)]
pub struct ShortCollector {
    buf: VecDeque<ShortLite>,
    hot_l: BTreeSet<i64>,
    hot_r: BTreeSet<i64>,
    kept: Vec<ShortLite>,
    pub window: i64,
    pub max: i64,
    pub horizon: i64,
}

impl ShortCollector {
    pub fn new(window: i64, max: i64) -> Self {
        ShortCollector { window, max, horizon: 1500, ..Default::default() }
    }

    fn is_hot(&self, r: &ShortLite) -> bool {
        let (w, m) = (self.window, self.max);
        // LEFT: (a) b in [start-w, start+w]; (b) b in [start+1, start+m]
        let l = self.hot_l.range(r.start - w..=r.start + w.max(m)).next().is_some();
        // RIGHT: (a) b in [end-w, end+w]; (b) b in [end-m, end-1]
        let rr = self.hot_r.range(r.end - w.max(m)..=r.end + w).next().is_some();
        l || rr
    }

    fn evict_before(&mut self, pos: i64) {
        while let Some(front) = self.buf.front() {
            if front.start >= pos {
                break;
            }
            let r = self.buf.pop_front().unwrap();
            if self.is_hot(&r) {
                self.kept.push(r);
            }
        }
        let cut = pos - 2 * self.horizon;
        while let Some(&x) = self.hot_l.iter().next() {
            if x >= cut {
                break;
            }
            self.hot_l.remove(&x);
        }
        while let Some(&x) = self.hot_r.iter().next() {
            if x >= cut {
                break;
            }
            self.hot_r.remove(&x);
        }
    }

    /// Add one eligible read (coordinate-sorted input).
    pub fn push(&mut self, r: ShortLite) {
        self.evict_before(r.start - self.horizon);
        self.buf.push_back(r);
    }

    pub fn hot(&mut self, b: i64, left: bool) {
        if left {
            self.hot_l.insert(b);
        } else {
            self.hot_r.insert(b);
        }
    }

    /// End of contig: flush the buffer and return the kept candidates (start-sorted).
    pub fn finish(&mut self) -> Vec<ShortLite> {
        self.evict_before(i64::MAX);
        self.hot_l.clear();
        self.hot_r.clear();
        let mut v = std::mem::take(&mut self.kept);
        v.sort_by_key(|r| (r.start, r.frag, r.flag));
        v
    }
}

/// Compact identity of a poly-A read, kept from the extract scan: its placement and
/// its mate's (the mate places the poly-A breakpoint, so `mref`/`mpos` predict where
/// the poly-A end can be emitted).
#[derive(Clone, Copy, Debug)]
pub struct PaLite {
    pub flag: u16,
    pub ref_id: i32,
    pub pos: i64,
    pub mref: i32,
    pub mpos: i64,
}

/// Poly-A read payload captured in the mate pass: the poly-A read itself and its
/// (anchoring) primary mate.
#[derive(Clone, Debug, Default)]
pub struct PaEv {
    pub read: Option<EvRec>,
    pub mate: Option<EvRec>,
}

/// Deterministic capped subset: keep the `cap` items with the smallest `key`
/// (ties keep input order), returned in their ORIGINAL relative order. Independent
/// of thread scheduling and of input order when keys are distinct.
pub fn select_lowest<T, K: Ord + Copy>(items: Vec<T>, cap: usize, key: impl Fn(&T) -> K) -> Vec<T> {
    if items.len() <= cap {
        return items;
    }
    let mut idx: Vec<usize> = (0..items.len()).collect();
    idx.sort_by_key(|&i| (key(&items[i]), i));
    let mut keep = vec![false; items.len()];
    for &i in idx.iter().take(cap) {
        keep[i] = true;
    }
    items.into_iter().zip(keep).filter_map(|(it, k)| if k { Some(it) } else { None }).collect()
}

/// A request to capture one record in the mate pass.
#[derive(Clone, Copy, Debug)]
pub struct MateReq {
    pub hash: u64,
    /// index into the target list (final LEFT / final RIGHT breakpoints, or poly-A reads)
    pub idx: u32,
    /// expected placement of the wanted record (self-requests only; -1 = unmapped)
    pub ref_id: i32,
    pub pos: i64,
    /// r12 of the wanted record
    pub r12: u8,
    pub kind: ReqKind,
    pub target: Target,
    /// raw flag bits the record must carry for ClipSelf (0x800 supplementary)
    pub supp: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReqKind {
    /// the primary mate of an evidence read (matched by qname hash + r12; a 64-bit
    /// collision needs ~1e9 records x 1e6 requests to reach 1e-4, so no placement check)
    Mate,
    /// a discordant anchor read itself (primary; placement checked)
    DiscSelf,
    /// a junction-clipped read itself (primary or supplementary as recorded; its
    /// breakpoint coordinate checked: start for LEFT, end for RIGHT)
    ClipSelf,
    /// a poly-A read itself (primary; placement checked, or unmapped)
    PolySelf,
    /// a short-overhang read itself (primary; placement checked)
    ShortSelf,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Target {
    Left,
    Right,
    PolyA,
}

impl MateReq {
    /// Does `read` (secondary/QC-fail/dup already filtered) satisfy this request? The
    /// qname hash is matched by the caller; this checks r12, primary/supplementary and,
    /// for self-requests, the recorded placement (which also rules out hash collisions).
    pub fn matches(&self, read: &BamRead<'_>) -> bool {
        if read.r12() != self.r12 {
            return false;
        }
        let same_ref = || read.reference_sequence_id.map(|i| i as i32) == Some(self.ref_id);
        match self.kind {
            ReqKind::Mate => !read.is_supplementary,
            ReqKind::DiscSelf | ReqKind::ShortSelf => {
                !read.is_supplementary && read.mapped && same_ref() && read.reference_start == self.pos
            }
            ReqKind::PolySelf => {
                !read.is_supplementary
                    && if self.ref_id < 0 { !read.mapped } else { read.mapped && same_ref() && read.reference_start == self.pos }
            }
            ReqKind::ClipSelf => {
                read.is_supplementary == self.supp
                    && read.mapped
                    && same_ref()
                    && self.pos == if self.target == Target::Left { read.reference_start } else { read.reference_end }
            }
        }
    }
}

pub const HEADER: &str = "locus\tside\trole\tfrag\tr12\tflag\tref\tpos\tstrand\touter\tmref\tmpos\tmstrand\ttlen\tmapq\tcigar\tclip_at\tseq\tqual\n";

/// Format one sidecar row (SPEC column order). `names` maps header ids to contig names.
pub fn write_row<W: Write>(w: &mut W, locus: &str, side: &str, r: &EvRec, names: &[String]) -> io::Result<()> {
    let name = |id: i32| -> &str {
        if id < 0 {
            "*"
        } else {
            names.get(id as usize).map(|s| s.as_str()).unwrap_or("*")
        }
    };
    let mapped = r.ref_id >= 0 && r.flag & FLAG_UNMAPPED == 0;
    let strand = if !mapped {
        "*"
    } else if r.flag & FLAG_REVERSE != 0 {
        "-"
    } else {
        "+"
    };
    let mstrand = if r.mref < 0 {
        "*"
    } else if r.flag & FLAG_MATE_REVERSE != 0 {
        "-"
    } else {
        "+"
    };
    let seq: &str = if r.seq.is_empty() { "*" } else { std::str::from_utf8(&r.seq).unwrap_or("*") };
    let qual: String = if r.qual.is_empty() || r.qual.iter().all(|&q| q == 0xFF) {
        "*".to_string()
    } else {
        r.qual.iter().map(|&q| (q.min(93) + 33) as char).collect()
    };
    writeln!(
        w,
        "{locus}\t{side}\t{}\t{:016x}\t{}\t{}\t{}\t{}\t{strand}\t{}\t{}\t{}\t{mstrand}\t{}\t{}\t{}\t{}\t{seq}\t{qual}",
        r.role.as_str(),
        r.frag,
        r.r12(),
        r.flag,
        if mapped { name(r.ref_id) } else { "*" },
        if mapped { r.pos } else { -1 },
        if mapped { r.outer } else { -1 },
        name(r.mref),
        if r.mref < 0 { -1 } else { r.mpos },
        r.tlen,
        r.mapq,
        r.cigar,
        r.clip_at,
    )
}

/// Open sidecar writer + header names, threaded through `Discovery::output`.
pub struct Sidecar<'a> {
    pub w: &'a mut dyn Write,
    pub names: Vec<String>,
    pub stats: SidecarStats,
}

impl Sidecar<'_> {
    pub fn row(&mut self, locus: &str, side: &str, r: &EvRec) -> io::Result<()> {
        self.stats.count(r.role);
        write_row(&mut self.w, locus, side, r, &self.names)
    }
}

/// Run-level sidecar counters, printed to stderr.
#[derive(Default, Debug)]
pub struct SidecarStats {
    pub loci: u64,
    pub rows_clip: u64,
    pub rows_polya: u64,
    pub rows_mate: u64,
    pub rows_disc: u64,
    pub rows_short: u64,
    /// SHORT rows per emitted breakpoint side (Bp ends only)
    pub short_per_bp: Vec<u32>,
    /// MATE rows per emitted breakpoint side (Bp ends only), for the distribution
    pub mates_per_bp: Vec<u32>,
}

impl SidecarStats {
    pub fn count(&mut self, role: Role) {
        match role {
            Role::Clip => self.rows_clip += 1,
            Role::PolyA => self.rows_polya += 1,
            Role::Mate => self.rows_mate += 1,
            Role::Disc => self.rows_disc += 1,
            Role::Short => self.rows_short += 1,
        }
    }

    pub fn summary(&mut self) -> String {
        let total = self.rows_clip + self.rows_polya + self.rows_mate + self.rows_disc + self.rows_short;
        self.mates_per_bp.sort_unstable();
        let n = self.mates_per_bp.len();
        let (med, max, mean) = if n == 0 {
            (0, 0, 0.0)
        } else {
            let s: u64 = self.mates_per_bp.iter().map(|&x| x as u64).sum();
            (self.mates_per_bp[n / 2], self.mates_per_bp[n - 1], s as f64 / n as f64)
        };
        let mut out = format!(
            "evidence sidecar: {} loci, {} rows (CLIP {}, POLYA {}, MATE {}, DISC {}); mates/breakpoint-side: mean {:.2} median {} max {} over {} sides",
            self.loci, total, self.rows_clip, self.rows_polya, self.rows_mate, self.rows_disc, mean, med, max, n
        );
        if self.rows_short > 0 || !self.short_per_bp.is_empty() {
            self.short_per_bp.sort_unstable();
            let k = self.short_per_bp.len();
            let (smed, smax) = if k == 0 { (0, 0) } else { (self.short_per_bp[k / 2], self.short_per_bp[k - 1]) };
            let with = self.short_per_bp.iter().filter(|&&x| x > 0).count();
            out.push_str(&format!(
                "; SHORT {} rows, per breakpoint-side median {} max {} ({} of {} sides with >= 1)",
                self.rows_short, smed, smax, with, k
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::RecordBuf;
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    use noodles_sam::Header;
    use std::num::NonZeroUsize;

    fn header() -> Header {
        Header::builder()
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(NonZeroUsize::new(100_000).unwrap()))
            .build()
    }

    /// 10S20M5S read at 1-based start 1001 (0-based 1000), 35 bases.
    fn rec(flags: u16) -> RecordBuf {
        let seq: Vec<u8> = b"ACGTACGTAC".iter().chain(b"GGGGGCCCCCAAAAATTTTT").chain(b"ACGTA").copied().collect();
        RecordBuf::builder()
            .set_name("read/1")
            .set_flags(Flags::from_bits_truncate(flags))
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(1001).unwrap())
            .set_cigar(vec![Op::new(Kind::SoftClip, 10), Op::new(Kind::Match, 20), Op::new(Kind::SoftClip, 5)].into())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(1301).unwrap())
            .set_template_length(-321)
            .set_sequence(seq.into())
            .set_quality_scores(vec![30u8; 35].into())
            .build()
    }

    #[test]
    fn outer_coordinate_forward_and_reverse_with_soft_clips() {
        let h = header();
        // forward: 5' end = start - leading soft clip = 1000 - 10
        let r = rec(0x1 | 0x40 | 0x20);
        let b = BamRead::from_record(&r, &h).unwrap();
        assert_eq!(b.outer(), 990);
        assert!(b.mate_is_reverse);
        assert_eq!(b.tlen, -321);
        // reverse: 5' end = last aligned base + trailing soft clip = 1019 + 5
        let r = rec(0x1 | 0x80 | 0x10);
        let b = BamRead::from_record(&r, &h).unwrap();
        assert_eq!(b.reference_end, 1020);
        assert_eq!(b.outer(), 1024);
        assert!(!b.mate_is_reverse);
    }

    #[test]
    fn outer_skips_hard_clip() {
        let h = header();
        let r = RecordBuf::builder()
            .set_flags(Flags::from_bits_truncate(0x10))
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(501).unwrap())
            .set_cigar(vec![Op::new(Kind::HardClip, 7), Op::new(Kind::SoftClip, 3), Op::new(Kind::Match, 10), Op::new(Kind::SoftClip, 4), Op::new(Kind::HardClip, 9)].into())
            .set_sequence(vec![b'A'; 17].into())
            .build();
        let b = BamRead::from_record(&r, &h).unwrap();
        assert_eq!(b.lead_soft, 3);
        assert_eq!(b.trail_soft, 4);
        assert_eq!(b.outer(), 500 + 10 - 1 + 4);
        assert_eq!(b.cigar_string(), "7H3S10M4S9H");
    }

    #[test]
    fn sidecar_row_formatting() {
        let h = header();
        let names = vec!["chr1".to_string()];
        let r = rec(0x1 | 0x80 | 0x10 | 0x400);
        let b = BamRead::from_record(&r, &h).unwrap();
        let ev = EvRec::from_read(&b, Role::Clip, 30);
        let mut out = Vec::new();
        write_row(&mut out, "chr1:1020-1035", "RIGHT", &ev, &names).unwrap();
        let line = String::from_utf8(out).unwrap();
        let f: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(f.len(), 19);
        assert_eq!(f[0], "chr1:1020-1035");
        assert_eq!(f[1], "RIGHT");
        assert_eq!(f[2], "CLIP");
        assert_eq!(f[3], format!("{:016x}", frag_hash(b"read/1")));
        assert_eq!(f[4], "2");
        assert_eq!(f[5], (0x1 | 0x80 | 0x10 | 0x400).to_string()); // dup bit kept
        assert_eq!(&f[6..10], &["chr1", "1000", "-", "1024"]);
        assert_eq!(&f[10..13], &["chr1", "1300", "+"]);
        assert_eq!(f[13], "-321");
        assert_eq!(f[15], "10S20M5S");
        assert_eq!(f[16], "30");
        assert_eq!(f[17], "ACGTACGTACGGGGGCCCCCAAAAATTTTTACGTA"); // stored = reference-forward
        assert_eq!(f[18], "?".repeat(35)); // phred 30 + 33
    }

    #[test]
    fn sidecar_row_unmapped_mate_and_record() {
        let names = vec!["chr1".to_string()];
        let h = header();
        // unmapped record whose mate is mapped: ref '*', pos -1, strand '*', seq kept
        let r = RecordBuf::builder()
            .set_name("q")
            .set_flags(Flags::from_bits_truncate(0x1 | 0x4 | 0x80 | 0x20))
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(1001).unwrap())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(1001).unwrap())
            .set_sequence(b"ACGTTT".to_vec().into())
            .set_quality_scores(vec![40u8; 6].into())
            .build();
        let b = BamRead::from_record(&r, &h).unwrap();
        let ev = EvRec::from_read(&b, Role::Mate, -1);
        let mut out = Vec::new();
        write_row(&mut out, "chr1:1-2", "LEFT", &ev, &names).unwrap();
        let line = String::from_utf8(out).unwrap();
        let f: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(&f[6..10], &["*", "-1", "*", "-1"]);
        assert_eq!(&f[10..13], &["chr1", "1000", "-"]);
        assert_eq!(f[15], "*");
        assert_eq!(f[17], "ACGTTT");
        // a record whose mate is unmapped: mref '*', mpos -1, mstrand '*'
        let r = RecordBuf::builder()
            .set_name("q")
            .set_flags(Flags::from_bits_truncate(0x1 | 0x8 | 0x40))
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(1001).unwrap())
            .set_cigar(vec![Op::new(Kind::Match, 6)].into())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(1001).unwrap())
            .set_sequence(b"ACGTTT".to_vec().into())
            .build();
        let b = BamRead::from_record(&r, &h).unwrap();
        let ev = EvRec::from_read(&b, Role::Clip, 3);
        let mut out = Vec::new();
        write_row(&mut out, "chr1:1-2", "LEFT", &ev, &names).unwrap();
        let line = String::from_utf8(out).unwrap();
        let f: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(&f[10..13], &["*", "-1", "*"]);
        assert_eq!(f[18], "*"); // no qualities
    }

    #[test]
    fn select_lowest_is_deterministic_and_order_independent() {
        let items: Vec<u64> = (0..500u64).map(|i| frag_hash(format!("q{i}").as_bytes())).collect();
        let a = select_lowest(items.clone(), 50, |&h| h);
        let mut shuffled = items.clone();
        shuffled.reverse();
        shuffled.rotate_left(137);
        let b = select_lowest(shuffled, 50, |&h| h);
        assert_eq!(a.len(), 50);
        let mut sa = a.clone();
        sa.sort();
        let mut sb = b.clone();
        sb.sort();
        assert_eq!(sa, sb, "same subset regardless of input order");
        // the subset is exactly the 50 smallest hashes
        let mut all = items.clone();
        all.sort();
        assert_eq!(sa, all[..50].to_vec());
        // original relative order preserved
        let pos: Vec<usize> = a.iter().map(|h| items.iter().position(|x| x == h).unwrap()).collect();
        assert!(pos.windows(2).all(|w| w[0] < w[1]));
        // under the cap: untouched
        assert_eq!(select_lowest(vec![3, 1, 2], 5, |&x| x), vec![3, 1, 2]);
    }

    fn sl(frag: u64, flag: u16, start: i64, end: i64, lead: u32, trail: u32) -> ShortLite {
        ShortLite { frag, flag, ref_id: 0, start, end, lead_soft: lead, trail_soft: trail }
    }

    #[test]
    fn short_candidate_modes_both_sides_both_strands() {
        let b = 1000;
        for flag in [0x1u16 | 0x40, 0x1 | 0x10 | 0x80] {
            // LEFT (a): 5 bp leading clip, alignment starts at b (+-window)
            assert!(short_candidate(&sl(1, flag, 1000, 1140, 5, 0), b, true, 3, 20, 12));
            assert!(short_candidate(&sl(1, flag, 1002, 1140, 5, 0), b, true, 3, 20, 12));
            assert!(!short_candidate(&sl(1, flag, 1005, 1140, 5, 0), b, true, 3, 20, 12)); // outside window
            assert!(!short_candidate(&sl(1, flag, 1000, 1140, 12, 0), b, true, 3, 20, 12)); // a CLIP read
            // LEFT (b): unclipped, aligned 1..20 bp past b into the clip side, far-side clip ok
            assert!(short_candidate(&sl(1, flag, 990, 1140, 0, 0), b, true, 3, 20, 12));
            assert!(short_candidate(&sl(1, flag, 980, 1100, 0, 30), b, true, 3, 20, 12));
            assert!(!short_candidate(&sl(1, flag, 979, 1130, 0, 0), b, true, 3, 20, 12)); // 21 bp: plain reference read
            assert!(!short_candidate(&sl(1, flag, 1000, 1150, 0, 0), b, true, 3, 20, 12)); // does not cross
            // RIGHT (a): 7 bp trailing clip, alignment ends at b
            assert!(short_candidate(&sl(1, flag, 860, 1000, 0, 7), b, false, 3, 20, 12));
            assert!(short_candidate(&sl(1, flag, 860, 997, 0, 7), b, false, 3, 20, 12));
            assert!(!short_candidate(&sl(1, flag, 860, 996, 0, 7), b, false, 3, 20, 12));
            // RIGHT (b): unclipped, alignment ends 1..20 bp past b
            assert!(short_candidate(&sl(1, flag, 870, 1015, 0, 0), b, false, 3, 20, 12));
            assert!(short_candidate(&sl(1, flag, 870, 1020, 40, 0), b, false, 3, 20, 12));
            assert!(!short_candidate(&sl(1, flag, 870, 1021, 0, 0), b, false, 3, 20, 12));
            assert!(!short_candidate(&sl(1, flag, 870, 1000, 0, 0), b, false, 3, 20, 12)); // ends exactly at b
        }
    }

    #[test]
    fn short_collector_keeps_only_reads_near_hot_junctions() {
        let mut c = ShortCollector::new(3, 20);
        c.push(sl(1, 0x41, 500, 650, 0, 0)); // nowhere near
        c.push(sl(2, 0x41, 990, 1140, 0, 0)); // LEFT (b) for b=1000 — hot marked later
        c.push(sl(3, 0x41, 1000, 1140, 30, 0)); // the clipped read itself
        c.hot(1000, true);
        c.push(sl(4, 0x41, 1001, 1150, 4, 0)); // LEFT (a)
        c.push(sl(5, 0x51, 1900, 2050, 0, 6)); // RIGHT (a) for b=2050, marked later
        c.hot(2050, false);
        c.push(sl(6, 0x41, 9000, 9150, 0, 0));
        let kept: Vec<u64> = c.finish().iter().map(|r| r.frag).collect();
        assert_eq!(kept, vec![2, 3, 4, 5]);
    }

    #[test]
    fn frag_hash_is_stable_fnv1a() {
        // FNV-1a 64 reference values
        assert_eq!(frag_hash(b""), 0xcbf29ce484222325);
        assert_eq!(frag_hash(b"a"), 0xaf63dc4c8601ec8c);
    }
}
