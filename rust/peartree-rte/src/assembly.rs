//! tools/rte/assembly.py -- per-insertion read layouts + covered-element consensus.
//!
//! STATUS: types FOUNDATION (implemented); algorithms = WORK PACKAGE "WP-ASM" (PORT_PLAN.md).
//!
//! Python functions to port (assembly.py line refs at eaa2718):
//!   Assembler.layout (178) + _gaps (258), _resolve (270), _edlib_targets_list (306),
//!   _rescue (312), _mark_local (347), _local_is_element (373), _element_is_local (392),
//!   _merge_ref_runs (423), _smooth_polya (450), _wide_local (520), assemble (548);
//!   AssemblyResult.build (590), _terminal_3p (635), _strand_from_polya (676),
//!   _strand_from_segments (705), _pileup (730), _nearest (807); _cut_hit_at (159).
//!
//! Golden: golden events `assemble` (in: ctx, junction_seqs, reads, strand_hint, cfg; out: the
//! whole AssemblyResult incl. raw_layouts and sense layouts). NOTE python aliasing: for
//! strand >= 0, `to_sense` copies the segment LIST but shares the Segment OBJECTS, so `_pileup`'s
//! in-place update of ELEMENT segments (target / t_st / t_en / identity) shows up in
//! `raw_layouts` too; for strand < 0 `flipped()` makes new objects and raw_layouts keep the
//! mappy coordinates. The golden `out.raw_layouts` captures that; reproduce it.

use crate::config::AssemblyCfg;
use crate::inputs::EvidenceRead;
use crate::library::RteLibrary;
use crate::mm::{Aligner, MapOpts};
use crate::sequtil::rc;
use std::sync::OnceLock;

/// Segment kinds (python strings).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SegKind {
    Ref,
    Local,
    Element,
    PolyA,
    Flank3p,
    Flank5p,
    Unknown,
}

impl SegKind {
    pub fn as_str(self) -> &'static str {
        match self {
            SegKind::Ref => "REF",
            SegKind::Local => "LOCAL",
            SegKind::Element => "ELEMENT",
            SegKind::PolyA => "POLYA",
            SegKind::Flank3p => "FLANK3P",
            SegKind::Flank5p => "FLANK5P",
            SegKind::Unknown => "UNKNOWN",
        }
    }
    pub fn parse(s: &str) -> Option<SegKind> {
        Some(match s {
            "REF" => SegKind::Ref,
            "LOCAL" => SegKind::Local,
            "ELEMENT" => SegKind::Element,
            "POLYA" => SegKind::PolyA,
            "FLANK3P" => SegKind::Flank3p,
            "FLANK5P" => SegKind::Flank5p,
            "UNKNOWN" => SegKind::Unknown,
            _ => return None,
        })
    }
}

/// assembly.Segment (half-open query interval; target coordinates -1 when unknown).
#[derive(Clone, Debug, PartialEq)]
pub struct Segment {
    pub q_st: i64,
    pub q_en: i64,
    pub kind: SegKind,
    /// consensus / flank / "site" / "wide" / "A" / "T" (POLYA base) / ""
    pub target: String,
    pub t_st: i64,
    pub t_en: i64,
    /// +1 / -1 / 0
    pub strand: i32,
    pub identity: f64,
    pub matches: i64,
    /// mappy cigar of an ELEMENT hit (stored, never read downstream)
    pub cigar: Option<Vec<(u32, u32)>>,
}

impl Segment {
    pub fn new(q_st: i64, q_en: i64, kind: SegKind) -> Segment {
        Segment { q_st, q_en, kind, target: String::new(), t_st: -1, t_en: -1, strand: 0, identity: 0.0, matches: 0, cigar: None }
    }
    pub fn qlen(&self) -> i64 {
        self.q_en - self.q_st
    }
    /// `flipped(n)`: the same segment on the reverse-complemented read of length n.
    pub fn flipped(&self, n: i64) -> Segment {
        Segment {
            q_st: n - self.q_en,
            q_en: n - self.q_st,
            strand: if self.strand != 0 { -self.strand } else { 0 },
            ..self.clone()
        }
    }
}

/// assembly.counts_as_fragment: False for genotype2 extra-pass reads (GT_* roles).
pub fn counts_as_fragment(lay: &ReadLayout) -> bool {
    !lay.role.starts_with("GT_")
}

/// assembly.ReadLayout
#[derive(Clone, Debug, PartialEq)]
pub struct ReadLayout {
    pub name: String,
    /// LEFT / RIGHT / '' (as delivered)
    pub side: String,
    /// JUNCTION (clip consensus), CLIP, POLYA, MATE, DISC, SPAN, GT_*
    pub role: String,
    /// (sample, frag); ("consensus", side) for a junction consensus
    pub frag_key: (String, String),
    /// upper case
    pub seq: Vec<u8>,
    pub segments: Vec<Segment>,
    /// True once flipped into element-sense orientation
    pub sense: bool,
}

impl ReadLayout {
    /// `to_sense(strand)`
    pub fn to_sense(&self, strand: i32) -> ReadLayout {
        if strand >= 0 {
            return ReadLayout { sense: true, ..self.clone() };
        }
        let n = self.seq.len() as i64;
        ReadLayout {
            segments: self.segments.iter().rev().map(|s| s.flipped(n)).collect(),
            seq: rc(&self.seq),
            sense: true,
            ..self.clone()
        }
    }
    pub fn ref_left(&self) -> bool {
        self.segments.first().is_some_and(|s| s.kind == SegKind::Ref)
    }
    pub fn ref_right(&self) -> bool {
        self.segments.last().is_some_and(|s| s.kind == SegKind::Ref)
    }
    pub fn piece(&self, s: &Segment) -> &[u8] {
        &self.seq[s.q_st.max(0) as usize..(s.q_en.max(0) as usize).min(self.seq.len())]
    }
}

/// assembly.SiteContext. The aligners are built lazily (python `aligner()` / `wide_aligner()`).
#[derive(Debug, Default)]
pub struct SiteContext {
    pub title: String,
    pub contig: Option<String>,
    /// L: first reference base of the LEFT junction record
    pub left_bp: Option<i64>,
    /// R: end (exclusive) of the RIGHT junction reference
    pub right_bp: Option<i64>,
    pub window_start: i64,
    pub window_seq: Vec<u8>,
    pub left_flank: Vec<u8>,
    pub right_flank: Vec<u8>,
    pub wide_start: i64,
    pub wide_seq: Vec<u8>,
    aligner: OnceLock<Option<Aligner>>,
    wide_aligner: OnceLock<Option<Aligner>>,
}

impl Clone for SiteContext {
    fn clone(&self) -> Self {
        SiteContext {
            title: self.title.clone(),
            contig: self.contig.clone(),
            left_bp: self.left_bp,
            right_bp: self.right_bp,
            window_start: self.window_start,
            window_seq: self.window_seq.clone(),
            left_flank: self.left_flank.clone(),
            right_flank: self.right_flank.clone(),
            wide_start: self.wide_start,
            wide_seq: self.wide_seq.clone(),
            aligner: OnceLock::new(),
            wide_aligner: OnceLock::new(),
        }
    }
}

impl SiteContext {
    pub fn new(title: &str, contig: Option<String>, left_bp: Option<i64>, right_bp: Option<i64>, left_flank: Vec<u8>, right_flank: Vec<u8>) -> SiteContext {
        SiteContext { title: title.to_string(), contig, left_bp, right_bp, left_flank, right_flank, ..Default::default() }
    }

    /// `local_reference()`: (sequence, offset) -- the genome window when available, else the
    /// two junction flanks joined by 30 N (offset None = no genome coordinates).
    pub fn local_reference(&self) -> (Vec<u8>, Option<i64>) {
        if !self.window_seq.is_empty() {
            return (self.window_seq.clone(), Some(self.window_start));
        }
        let mut s = self.right_flank.clone();
        s.extend(std::iter::repeat_n(b'N', 30));
        s.extend_from_slice(&self.left_flank);
        (s, None)
    }

    /// `aligner()`: mappy index of the local reference (MapOpts::LOCAL), None when < 30 bp.
    pub fn aligner(&self) -> Option<&Aligner> {
        self.aligner
            .get_or_init(|| {
                let (seq, _) = self.local_reference();
                if seq.len() >= 30 {
                    Aligner::from_seq(&seq, MapOpts::LOCAL)
                } else {
                    None
                }
            })
            .as_ref()
    }

    /// `wide_aligner()`: mappy `preset="sr"` index of the wide window, None when < 100 bp.
    pub fn wide_aligner(&self) -> Option<&Aligner> {
        self.wide_aligner
            .get_or_init(|| if self.wide_seq.len() >= 100 { Aligner::from_seq(&self.wide_seq, MapOpts::SR) } else { None })
            .as_ref()
    }
}

/// assembly.AssemblyResult
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AssemblyResult {
    /// element strand on the reference (+1 / -1, 0 unknown)
    pub strand: i32,
    pub strand_source: String,
    /// L1 / ALU / SVA / '' (none)  (python `element_class`)
    pub element_class: String,
    /// best consensus name ('' none)
    pub consensus: String,
    pub element_bp: i64,
    /// element-sense layouts
    pub layouts: Vec<ReadLayout>,
    /// reference-forward layouts
    pub raw_layouts: Vec<ReadLayout>,
    /// merged covered intervals on the consensus
    pub covered: Vec<(i64, i64)>,
    pub covered_seqs: Vec<Vec<u8>>,
    /// (s, e, sense, layout index)
    pub segments_on_cons: Vec<(i64, i64, bool, usize)>,
    pub consensus_identity: f64,
    pub nearest_intact: String,
    pub nearest_intact_identity: f64,
    pub nearest_active: String,
    pub element_identity: f64,
    /// class -> bp, python dict order
    pub class_bp: Vec<(String, i64)>,
}

impl AssemblyResult {
    pub fn empty() -> AssemblyResult {
        AssemblyResult { nearest_intact: ".".into(), nearest_active: ".".into(), ..Default::default() }
    }
    /// `covered_5p` (-1 when nothing is covered)
    pub fn covered_5p(&self) -> i64 {
        self.covered.iter().map(|c| c.0).min().unwrap_or(-1)
    }
    /// `covered_3p` (-1 when nothing is covered)
    pub fn covered_3p(&self) -> i64 {
        self.covered.iter().map(|c| c.1).max().unwrap_or(-1)
    }
}

/// One junction string handed to the assembler: (sequence, reference interval (q_st, q_en)).
pub type JunctionSeq = (Vec<u8>, (i64, i64));

/// assembly.Assembler
pub struct Assembler<'a> {
    pub lib: &'a RteLibrary,
    pub cfg: AssemblyCfg,
}

impl<'a> Assembler<'a> {
    pub fn new(lib: &'a RteLibrary, cfg: AssemblyCfg) -> Assembler<'a> {
        Assembler { lib, cfg }
    }

    /// `layout(seq, ctx, name, side, role, frag_key, ref_interval, rescue_flanks)`. WP-ASM.
    #[allow(clippy::too_many_arguments)]
    pub fn layout(
        &self,
        _seq: &[u8],
        _ctx: &SiteContext,
        _name: &str,
        _side: &str,
        _role: &str,
        _frag_key: (String, String),
        _ref_interval: Option<(i64, i64)>,
        _rescue_flanks: bool,
    ) -> ReadLayout {
        todo!("WP-ASM: port assembly.Assembler.layout")
    }

    /// `assemble(ctx, junction_seqs, reads, strand_hint)`: junction_seqs in python dict order
    /// (LEFT before RIGHT as the annotator inserts them). WP-ASM.
    pub fn assemble(
        &self,
        _ctx: &SiteContext,
        _junction_seqs: &[(String, JunctionSeq)],
        _reads: &[EvidenceRead],
        _strand_hint: Option<(i32, String)>,
    ) -> AssemblyResult {
        todo!("WP-ASM: port assembly.Assembler.assemble + AssemblyResult.build")
    }
}

/// `AssemblyResult._strand_from_polya(layouts, min_run=10)` -> (strand, [tail lengths]); also
/// used by the annotator (polya_reads). WP-ASM.
pub fn strand_from_polya(_layouts: &[ReadLayout], _min_run: i64) -> (i32, Vec<i64>) {
    todo!("WP-ASM: port AssemblyResult._strand_from_polya")
}
