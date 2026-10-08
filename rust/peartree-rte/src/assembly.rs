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
use crate::library::AlignerKind;
use crate::mm::Hit;
use crate::pyfmt::py_round;
use crate::sequtil::{base_fraction, edlib_best, edlib_path, polya_runs, rc};
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

// ----------------------------------------------------------------------------- helpers

/// python `seq[a:b]` (clamped; empty when a >= b).
fn sl(seq: &[u8], a: i64, b: i64) -> &[u8] {
    let n = seq.len() as i64;
    let a = a.clamp(0, n);
    let b = b.clamp(0, n);
    if a >= b {
        &[]
    } else {
        &seq[a as usize..b as usize]
    }
}

/// `max(piece.count("A"), piece.count("T")) > 0.8 * len(piece)`
fn homopolymer_like(piece: &[u8]) -> bool {
    let a = piece.iter().filter(|&&c| c == b'A').count();
    let t = piece.iter().filter(|&&c| c == b'T').count();
    (a.max(t) as f64) > 0.8 * piece.len() as f64
}

/// A small insertion-ordered counter (python `Counter` / dict with first-max semantics).
struct OrdCount<K: PartialEq>(Vec<(K, i64)>);

impl<K: PartialEq> OrdCount<K> {
    fn new() -> Self {
        OrdCount(Vec::new())
    }
    fn add(&mut self, k: K, v: i64) {
        bump(&mut self.0, k, v);
    }
    /// `max(d, key=d.get)` -> FIRST maximal key
    fn argmax(&self) -> Option<&K> {
        first_max(&self.0).map(|e| &e.0)
    }
}

fn bump<K: PartialEq>(v: &mut Vec<(K, i64)>, k: K, d: i64) {
    match v.iter_mut().find(|e| e.0 == k) {
        Some(e) => e.1 += d,
        None => v.push((k, d)),
    }
}

/// python `max(items, key=count)`: the FIRST entry with the maximal count.
fn first_max<K>(v: &[(K, i64)]) -> Option<&(K, i64)> {
    let mut best: Option<&(K, i64)> = None;
    for e in v {
        if best.is_none_or(|b| e.1 > b.1) {
            best = Some(e);
        }
    }
    best
}

/// `_cut_hit_at(h, cons_end)` -> (q_st, q_en, r_en)
fn cut_hit_at(h: &Hit, cons_end: i64) -> (i64, i64, i64) {
    if h.r_en <= cons_end {
        return (h.q_st, h.q_en, h.r_en);
    }
    let cut = h.r_en - cons_end;
    if h.strand == 1 {
        return (h.q_st, h.q_st.max(h.q_en - cut), cons_end);
    }
    (h.q_en.min(h.q_st + cut), h.q_en, cons_end)
}

#[allow(clippy::too_many_arguments)]
fn seg_full(q_st: i64, q_en: i64, kind: SegKind, target: &str, t_st: i64, t_en: i64, strand: i32, identity: f64, matches: i64) -> Segment {
    Segment { q_st, q_en, kind, target: target.to_string(), t_st, t_en, strand, identity, matches, cigar: None }
}

/// `(i == 0 and s.q_st <= 3) or (i == len(segs) - 1 and s.q_en >= len(seq) - 3)`
fn at_read_end(i: usize, nsegs: usize, s: &Segment, seqlen: usize) -> bool {
    (i == 0 && s.q_st <= 3) || (i + 1 == nsegs && s.q_en >= seqlen as i64 - 3)
}

/// stable `sorted(segs, key=lambda s: s.q_st)`
fn sort_q(segs: &mut [Segment]) {
    segs.sort_by_key(|s| s.q_st);
}

/// (edit distance, name, t_st, t_en, strand)
type NamedEd<'n> = (i32, &'n str, i64, i64, i32);

impl<'a> Assembler<'a> {
    pub fn new(lib: &'a RteLibrary, cfg: AssemblyCfg) -> Assembler<'a> {
        Assembler { lib, cfg }
    }

    /// `_edlib_targets_list()`: (name, consensus cut at cons_end), consensus file order.
    fn edlib_targets(&self) -> impl Iterator<Item = (&str, &[u8])> {
        self.lib.consensus.iter().map(|(n, s)| {
            let end = self.lib.cons_end.get(n).copied().unwrap_or(s.len()).min(s.len());
            (n, &s[..end])
        })
    }

    /// `layout(seq, ctx, name, side, role, frag_key, ref_interval, rescue_flanks)`.
    #[allow(clippy::too_many_arguments)]
    pub fn layout(
        &self,
        seq: &[u8],
        ctx: &SiteContext,
        name: &str,
        side: &str,
        role: &str,
        frag_key: (String, String),
        ref_interval: Option<(i64, i64)>,
        rescue_flanks: bool,
    ) -> ReadLayout {
        let c = &self.cfg;
        let seq = seq.to_ascii_uppercase();
        let n = seq.len() as i64;
        let mut cands: Vec<(f64, Segment)> = Vec::new();
        match ref_interval {
            Some((a, b)) if b > a => {
                cands.push((3e6, seg_full(a, b, SegKind::Ref, "site", -1, -1, 0, 1.0, b - a)));
            }
            _ => {
                if let Some(al) = ctx.aligner() {
                    for h in al.map(&seq) {
                        let idn = h.identity();
                        if idn < c.min_ref_identity || h.q_en - h.q_st < c.min_segment_len {
                            continue;
                        }
                        cands.push((2e6 + h.mlen as f64 + 2.0, seg_full(h.q_st, h.q_en, SegKind::Ref, "site", h.r_st, h.r_en, h.strand, idn, h.mlen)));
                    }
                }
            }
        }
        if let Some(al) = self.lib.aligner(AlignerKind::Consensus) {
            for h in al.map(&seq) {
                let idn = h.identity();
                if idn < c.min_element_identity {
                    continue;
                }
                let cons_end = self.lib.cons_end.get(&*h.ctg).map(|&e| e as i64).unwrap_or(h.ctg_len);
                if h.r_st >= cons_end {
                    continue;
                }
                let (qs, qe, re_) = cut_hit_at(&h, cons_end);
                if qe - qs < c.min_segment_len {
                    continue;
                }
                let matches = ((h.mlen * (qe - qs)) as f64 / (h.q_en - h.q_st).max(1) as f64) as i64;
                let mut sg = seg_full(qs, qe, SegKind::Element, &h.ctg, h.r_st, re_, h.strand, idn, matches);
                sg.cigar = Some(h.cigar.clone());
                cands.push((2e6 + h.mlen as f64 * idn, sg));
            }
        }
        for (kind, key) in [(SegKind::Flank3p, AlignerKind::Flanks3), (SegKind::Flank5p, AlignerKind::Flanks5)] {
            let Some(al) = self.lib.aligner(key) else { continue };
            for h in al.map(&seq) {
                let idn = h.identity();
                if idn < c.min_flank_identity || h.q_en - h.q_st < c.min_segment_len {
                    continue;
                }
                // a flank hit made of the source's own poly-A remnant is not informative
                if homopolymer_like(sl(&seq, h.q_st, h.q_en)) {
                    continue;
                }
                cands.push((h.mlen as f64 * idn - 1.0, seg_full(h.q_st, h.q_en, kind, &h.ctg, h.r_st, h.r_en, h.strand, idn, h.mlen)));
            }
        }
        // poly-A runs compete as candidates
        for base in [b'A', b'T'] {
            for (a, b) in polya_runs(&seq, base, c.polya_min_run.max(0) as usize, 1) {
                let (a, b) = (a as i64, b as i64);
                let st = if base == b'A' { 1 } else { -1 };
                cands.push((1e6 + b as f64 - a as f64, seg_full(a, b, SegKind::PolyA, if base == b'A' { "A" } else { "T" }, -1, -1, st, 1.0, b - a)));
            }
        }
        let mut accepted: Vec<Segment> = Vec::new();
        for s in Self::resolve(cands) {
            if s.kind == SegKind::PolyA {
                if s.qlen() >= c.polya_min_run {
                    accepted.push(s);
                }
                continue;
            }
            if s.kind != SegKind::Ref && (s.qlen() < c.min_segment_len || homopolymer_like(sl(&seq, s.q_st, s.q_en))) {
                continue; // trimmed to a sliver, or a homopolymer masquerading as a hit
            }
            accepted.push(s);
        }
        sort_q(&mut accepted);
        // edlib rescue of the remaining gaps
        for (ga, gb) in Self::gaps(&accepted, n, c.min_gap_rescue) {
            let seg = self.rescue(sl(&seq, ga, gb), ga, ctx, rescue_flanks);
            accepted.push(seg);
        }
        sort_q(&mut accepted);
        let mut lay = ReadLayout { name: name.to_string(), side: side.to_string(), role: role.to_string(), frag_key, seq, segments: accepted, sense: false };
        Self::mark_local(&mut lay);
        self.local_is_element(&mut lay);
        self.element_is_local(&mut lay, ctx);
        self.smooth_polya(&mut lay);
        self.wide_local(&mut lay, ctx);
        lay
    }

    /// `_gaps(segs, n, min_len)`
    fn gaps(segs: &[Segment], n: i64, min_len: i64) -> Vec<(i64, i64)> {
        let mut sorted: Vec<&Segment> = segs.iter().collect();
        sorted.sort_by_key(|s| s.q_st);
        let mut out = Vec::new();
        let mut pos = 0i64;
        for s in sorted {
            if s.q_st - pos >= min_len {
                out.push((pos, s.q_st));
            }
            pos = pos.max(s.q_en);
        }
        if n - pos >= min_len {
            out.push((pos, n));
        }
        out
    }

    /// `_resolve(cands, n)`: greedy by score (stable sort over insertion order); a candidate
    /// overlapping accepted ones is cut to its largest free part (FIRST largest).
    fn resolve(mut cands: Vec<(f64, Segment)>) -> Vec<Segment> {
        cands.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut acc: Vec<Segment> = Vec::new();
        for (_, s) in cands {
            let mut free = vec![(s.q_st, s.q_en)];
            for a in &acc {
                let mut nxt = Vec::new();
                for &(fa, fb) in &free {
                    if a.q_en <= fa || a.q_st >= fb {
                        nxt.push((fa, fb));
                        continue;
                    }
                    if fa < a.q_st {
                        nxt.push((fa, a.q_st));
                    }
                    if a.q_en < fb {
                        nxt.push((a.q_en, fb));
                    }
                }
                free = nxt;
            }
            if free.is_empty() {
                continue;
            }
            let mut best = free[0];
            for &f in &free[1..] {
                if f.1 - f.0 > best.1 - best.0 {
                    best = f;
                }
            }
            let (fa, fb) = best;
            if fb - fa < if s.kind == SegKind::PolyA { 8 } else { 15 } {
                continue;
            }
            let s = if (fa, fb) != (s.q_st, s.q_en) {
                let frac = (fb - fa) as f64 / s.qlen().max(1) as f64;
                let (mut t_st, mut t_en) = (s.t_st, s.t_en);
                if s.t_st >= 0 {
                    // proportional; exact coords are recomputed in the pileup
                    if s.strand >= 0 {
                        t_st = s.t_st + (fa - s.q_st);
                        t_en = s.t_en - (s.q_en - fb);
                    } else {
                        t_st = s.t_st + (s.q_en - fb);
                        t_en = s.t_en - (fa - s.q_st);
                    }
                }
                Segment { q_st: fa, q_en: fb, kind: s.kind, target: s.target, t_st, t_en, strand: s.strand, identity: s.identity, matches: (s.matches as f64 * frac) as i64, cigar: None }
            } else {
                s
            };
            acc.push(s);
        }
        acc
    }

    /// `_rescue(piece, offset, ctx, rescue_flanks)`
    fn rescue(&self, piece: &[u8], offset: i64, ctx: &SiteContext, rescue_flanks: bool) -> Segment {
        let c = &self.cfg;
        let l = piece.len() as i64;
        let lf = l as f64;
        let mf = c.rescue_max_edit_frac;
        let mut best: Option<(f64, Segment)> = None;
        // local site (both strands) first: templated / REF extension
        let (refseq, _off) = ctx.local_reference();
        if let Some((ed, ts, te, strand)) = edlib_best(piece, &refseq, 0.10, true) {
            let edf = ed as f64;
            best = Some((edf / lf - 0.02, seg_full(offset, offset + l, SegKind::Ref, "site", ts, te, strand, 1.0 - edf / lf, l - ed as i64)));
        }
        for (name, cons) in self.edlib_targets() {
            let Some((ed, ts, te, strand)) = edlib_best(piece, cons, mf, true) else { continue };
            let edf = ed as f64;
            if best.as_ref().is_none_or(|b| edf / lf < b.0) {
                best = Some((edf / lf, seg_full(offset, offset + l, SegKind::Element, name, ts, te, strand, 1.0 - edf / lf, l - ed as i64)));
            }
        }
        if rescue_flanks {
            for (kind, fl) in [(SegKind::Flank3p, &self.lib.flanks3), (SegKind::Flank5p, &self.lib.flanks5)] {
                for (name, s) in fl.iter() {
                    let Some((ed, ts, te, strand)) = edlib_best(piece, s, 0.10, true) else { continue };
                    let edf = ed as f64;
                    if best.as_ref().is_none_or(|b| edf / lf < b.0) {
                        best = Some((edf / lf, seg_full(offset, offset + l, kind, name, ts, te, strand, 1.0 - edf / lf, l - ed as i64)));
                    }
                }
            }
        }
        match best {
            Some((_, s)) => s,
            None => Segment::new(offset, offset + l, SegKind::Unknown),
        }
    }

    /// `_mark_local(lay, ctx)`
    fn mark_local(lay: &mut ReadLayout) {
        Self::merge_ref_runs(lay, 20);
        let n = lay.segments.len();
        let seqlen = lay.seq.len();
        for i in 0..n {
            let s = &mut lay.segments[i];
            if s.kind == SegKind::Ref && !at_read_end(i, n, s, seqlen) {
                s.kind = SegKind::Local;
            }
        }
        if !(lay.side == "LEFT" || lay.side == "RIGHT") || n < 2 {
            return;
        }
        let (a, b) = (&lay.segments[0], &lay.segments[n - 1]);
        if a.kind == SegKind::Ref && b.kind == SegKind::Ref && n == 2 && b.q_st - a.q_en <= 10 && a.t_st >= 0 && b.t_st >= 0 {
            let contiguous = a.strand == b.strand && if a.strand >= 0 { (b.t_st - a.t_en).abs() <= 10 } else { (a.t_st - b.t_en).abs() <= 10 };
            if !contiguous {
                let idx = if lay.side == "RIGHT" { n - 1 } else { 0 };
                lay.segments[idx].kind = SegKind::Local;
            }
        }
    }

    /// `_local_is_element(lay)`
    fn local_is_element(&self, lay: &mut ReadLayout) {
        for si in 0..lay.segments.len() {
            let s = &lay.segments[si];
            if s.kind != SegKind::Local || s.qlen() < 30 {
                continue; // a short piece matches somewhere in 6 kb of consensus by chance
            }
            let piece = sl(&lay.seq, s.q_st, s.q_en);
            let mut best: Option<NamedEd> = None;
            for (name, cons) in self.edlib_targets() {
                if let Some(r) = edlib_best(piece, cons, 0.10, true) {
                    if best.is_none_or(|b| r.0 < b.0) {
                        best = Some((r.0, name, r.1, r.2, r.3));
                    }
                }
            }
            if let Some((ed, name, ts, te, strand)) = best {
                let plen = piece.len() as i64;
                let s = &mut lay.segments[si];
                s.kind = SegKind::Element;
                s.target = name.to_string();
                s.t_st = ts;
                s.t_en = te;
                s.strand = strand;
                s.identity = 1.0 - ed as f64 / plen.max(1) as f64;
                s.matches = plen - ed as i64;
            }
        }
    }

    /// `_element_is_local(lay, ctx)`
    fn element_is_local(&self, lay: &mut ReadLayout, ctx: &SiteContext) {
        let c = &self.cfg;
        let (refseq, off) = ctx.local_reference();
        if refseq.len() < 30 {
            return;
        }
        let n = lay.segments.len();
        let seqlen = lay.seq.len();
        for i in 0..n {
            let s = &lay.segments[i];
            if s.kind != SegKind::Element || s.qlen() < c.min_segment_len {
                continue;
            }
            let piece = sl(&lay.seq, s.q_st, s.q_en);
            let Some((ed, ts, te, strand)) = edlib_best(piece, &refseq, 1.0 - c.element_local_min_identity, true) else { continue };
            let plen = piece.len() as i64;
            let idn = 1.0 - ed as f64 / plen.max(1) as f64;
            if idn < c.element_local_min_identity || idn < s.identity + c.element_local_margin {
                continue;
            }
            let end = at_read_end(i, n, s, seqlen);
            let s = &mut lay.segments[i];
            s.kind = if end { SegKind::Ref } else { SegKind::Local };
            s.target = "site".to_string();
            (s.t_st, s.t_en) = if off.is_some() { (ts, te) } else { (-1, -1) };
            s.strand = strand;
            s.identity = idn;
            s.matches = plen - ed as i64;
        }
    }

    /// `_merge_ref_runs(lay, max_gap=20)`
    fn merge_ref_runs(lay: &mut ReadLayout, max_gap: i64) {
        let mut segs = lay.segments.clone();
        sort_q(&mut segs);
        let mut out = Vec::new();
        let mut i = 0;
        while i < segs.len() {
            let a = &segs[i];
            let mut j = i + 1;
            if a.kind == SegKind::Ref && a.t_st >= 0 {
                while j < segs.len() && segs[j].kind != SegKind::Ref && segs[j].q_en - a.q_en <= max_gap {
                    j += 1;
                }
                if j < segs.len() && segs[j].kind == SegKind::Ref && segs[j].t_st >= 0 && segs[j].strand == a.strand {
                    let b = &segs[j];
                    let qgap = b.q_st - a.q_en;
                    let rgap = if a.strand >= 0 { b.t_st - a.t_en } else { a.t_st - b.t_en };
                    if qgap <= max_gap && (-5..=max_gap).contains(&rgap) && (qgap - rgap).abs() <= max_gap {
                        let m = seg_full(a.q_st, b.q_en, SegKind::Ref, &a.target, a.t_st.min(b.t_st), a.t_en.max(b.t_en), a.strand, a.identity.min(b.identity), a.matches + b.matches);
                        segs.splice(i..=j, [m]);
                        continue;
                    }
                }
            }
            out.push(segs[i].clone());
            i += 1;
        }
        lay.segments = out;
    }

    /// `_smooth_polya(lay)`
    fn smooth_polya(&self, lay: &mut ReadLayout) {
        let c = &self.cfg;
        let mut segs = lay.segments.clone();
        sort_q(&mut segs);
        if !segs.iter().any(|x| x.kind == SegKind::PolyA) {
            return; // lay.segments stays as it was (python returns before reassigning)
        }
        let n0 = segs.len();
        let flank_idx: Option<usize> = match lay.side.as_str() {
            "RIGHT" => Some(0),
            "LEFT" => Some(n0 - 1),
            _ => None,
        };
        let unsided = lay.side.is_empty();
        let junction = lay.role == "JUNCTION";
        // python closure: `len(segs)` is the CURRENT list (late binding), flank_idx the initial
        let is_flank = |i: usize, x: &Segment, len: usize| -> bool {
            if x.kind != SegKind::Ref {
                return false;
            }
            if junction || flank_idx == Some(i) {
                return true;
            }
            unsided && (i == 0 || i + 1 == len)
        };
        loop {
            let mut changed = false;
            'outer: for i in 0..segs.len() {
                if segs[i].kind != SegKind::PolyA {
                    continue;
                }
                let x = segs[i].clone();
                let base: String = if !x.target.is_empty() {
                    x.target.clone()
                } else if x.strand >= 0 {
                    "A".into()
                } else {
                    "T".into()
                };
                let len = segs.len() as i64;
                for j in [i as i64 - 1, i as i64 + 1] {
                    if !(0 <= j && j < len) {
                        continue;
                    }
                    let ju = j as usize;
                    let y = &segs[ju];
                    if y.kind == SegKind::PolyA || is_flank(ju, y, segs.len()) || (y.kind == SegKind::Element && y.qlen() > 30) {
                        continue;
                    }
                    let piece = sl(&lay.seq, y.q_st, y.q_en);
                    let k = j + (j - i as i64); // the segment beyond y
                    let sandwiched = 0 <= k && k < len && segs[k as usize].kind == SegKind::PolyA && segs[k as usize].target == base;
                    let terminal = j == 0 || j == len - 1;
                    let rich = base_fraction(piece, base.as_bytes().first().copied().unwrap_or(0)) >= c.polya_absorb_base_frac;
                    if (sandwiched && (y.qlen() <= c.polya_absorb_max_gap || rich)) || (terminal && rich) {
                        let mut lo = x.q_st.min(y.q_st);
                        let mut hi = x.q_en.max(y.q_en);
                        if sandwiched {
                            lo = lo.min(segs[k as usize].q_st);
                            hi = hi.max(segs[k as usize].q_en);
                        }
                        let merged = seg_full(lo, hi, SegKind::PolyA, &base, -1, -1, x.strand, 1.0, hi - lo);
                        let ku = if sandwiched { Some(k as usize) } else { None };
                        let mut next: Vec<Segment> = segs.iter().enumerate().filter(|(t, _)| *t != i && *t != ju && Some(*t) != ku).map(|(_, s)| s.clone()).collect();
                        next.push(merged);
                        sort_q(&mut next);
                        segs = next;
                        changed = true;
                        break 'outer;
                    }
                }
            }
            if !changed {
                break;
            }
        }
        // a gap of <= max_gap bp between two same-base runs (no segment in between)
        let mut out: Vec<Segment> = Vec::new();
        for s in segs {
            if let Some(p) = out.last_mut() {
                if s.kind == SegKind::PolyA && p.kind == SegKind::PolyA && p.target == s.target && s.q_st - p.q_en <= c.polya_absorb_max_gap {
                    *p = seg_full(p.q_st, s.q_en, SegKind::PolyA, &p.target, -1, -1, p.strand, 1.0, s.q_en - p.q_st);
                    continue;
                }
            }
            out.push(s);
        }
        lay.segments = out;
    }

    /// `_wide_local(lay, ctx)`
    fn wide_local(&self, lay: &mut ReadLayout, ctx: &SiteContext) {
        let c = &self.cfg;
        if ctx.wide_seq.is_empty() {
            return;
        }
        let cand: Vec<usize> = (0..lay.segments.len()).filter(|&i| lay.segments[i].kind == SegKind::Unknown && lay.segments[i].qlen() >= c.wide_min_bp).collect();
        if cand.is_empty() {
            return;
        }
        let Some(al) = ctx.wide_aligner() else { return };
        let hits: Vec<Hit> = al.map(&lay.seq).into_iter().filter(|h| h.identity() >= c.min_ref_identity).collect();
        for si in cand {
            let s = &lay.segments[si];
            let mut best: Option<(i64, &Hit)> = None;
            for h in &hits {
                let ov = s.q_en.min(h.q_en) - s.q_st.max(h.q_st);
                if ov as f64 >= 0.8 * s.qlen() as f64 && best.is_none_or(|b| ov > b.0) {
                    best = Some((ov, h));
                }
            }
            let Some((_, h)) = best else { continue };
            let s = &mut lay.segments[si];
            s.kind = SegKind::Local;
            s.target = "wide".to_string();
            s.t_st = ctx.wide_start + h.r_st;
            s.t_en = ctx.wide_start + h.r_en;
            s.strand = h.strand;
            s.identity = h.identity();
            s.matches = h.mlen;
        }
    }

    /// `assemble(ctx, junction_seqs, reads, strand_hint)`: junction_seqs in python dict order
    /// (LEFT before RIGHT as the annotator inserts them).
    pub fn assemble(&self, ctx: &SiteContext, junction_seqs: &[(String, JunctionSeq)], reads: &[EvidenceRead], strand_hint: Option<(i32, String)>) -> AssemblyResult {
        let mut layouts = Vec::with_capacity(junction_seqs.len() + reads.len());
        for (side, (seq, ref_iv)) in junction_seqs {
            if !seq.is_empty() {
                layouts.push(self.layout(seq, ctx, &format!("junction_{side}"), side, "JUNCTION", ("consensus".to_string(), side.clone()), Some(*ref_iv), true));
            }
        }
        for r in reads {
            let name = format!("{}|{}|{}|{}|{}", r.side, r.role, r.sample, r.frag, r.r12);
            layouts.push(self.layout(&r.seq(), ctx, &name, &r.side, &r.role, (r.sample.to_string(), r.frag.to_string()), None, false));
        }
        AssemblyResult::build(self, layouts, strand_hint)
    }
}

/// `AssemblyResult._strand_from_polya(layouts, min_run=10)` -> (strand, [tail lengths]); also
/// used by the annotator (polya_reads). GT_* reads never vote (`counts_as_fragment`).
pub fn strand_from_polya(layouts: &[ReadLayout], min_run: i64) -> (i32, Vec<i64>) {
    // index 0 = +1, 1 = -1 (only emptiness of the vote sets matters)
    let mut votes = [0usize; 2];
    let mut lens: [Vec<i64>; 2] = [Vec::new(), Vec::new()];
    for lay in layouts {
        if !counts_as_fragment(lay) {
            continue;
        }
        let segs = &lay.segments;
        for (i, s) in segs.iter().enumerate() {
            if s.kind != SegKind::PolyA || s.qlen() < min_run {
                continue;
            }
            let base = if !s.target.is_empty() {
                s.target.as_str()
            } else if s.strand >= 0 {
                "A"
            } else {
                "T"
            };
            if base == "A" && i + 1 < segs.len() && segs[i + 1].kind == SegKind::Ref {
                votes[0] += 1;
                lens[0].push(s.qlen());
            } else if base == "T" && i > 0 && segs[i - 1].kind == SegKind::Ref {
                votes[1] += 1;
                lens[1].push(s.qlen());
            }
        }
    }
    for (k, st) in [(0usize, 1i32), (1, -1)] {
        if votes[k] > 0 && votes[1 - k] == 0 {
            return (st, std::mem::take(&mut lens[k]));
        }
    }
    (0, Vec::new())
}

impl AssemblyResult {
    /// `AssemblyResult.build(asm, ctx, raw_layouts, strand_hint)`
    fn build(asm: &Assembler, raw_layouts: Vec<ReadLayout>, strand_hint: Option<(i32, String)>) -> AssemblyResult {
        let lib = asm.lib;
        let mut res = AssemblyResult::empty();
        // --- class / consensus by aligned matches
        let mut per_cons: OrdCount<String> = OrdCount::new();
        let mut per_class: OrdCount<String> = OrdCount::new();
        for lay in &raw_layouts {
            let w = if lay.role == "JUNCTION" { 2 } else { 1 };
            for s in &lay.segments {
                if s.kind == SegKind::Element {
                    per_cons.add(s.target.clone(), s.matches * w);
                    per_class.add(lib.class_of(&s.target).unwrap_or("OTHER").to_string(), s.qlen());
                }
            }
        }
        res.class_bp = per_class.0.clone();
        res.element_bp = per_class.0.iter().map(|e| e.1).sum();
        let (term_strand, term) = Self::terminal_3p(&raw_layouts, lib, &asm.cfg);
        if !per_cons.0.is_empty() && res.element_bp >= asm.cfg.min_element_bp {
            let best_class = per_class.argmax().unwrap().clone();
            let mut best: Option<&(String, i64)> = None;
            for e in &per_cons.0 {
                if lib.class_of(&e.0) == Some(best_class.as_str()) && best.is_none_or(|b| e.1 > b.1) {
                    best = Some(e);
                }
            }
            let Some(best) = best else { panic!("ValueError: max() iterable argument is empty") };
            res.consensus = best.0.clone();
            res.element_class = best_class;
        } else if term_strand != 0 {
            // a piece ending AT the consensus 3' terminus next to the tail establishes the class
            let mut best: OrdCount<String> = OrdCount::new();
            for sg in &term {
                best.add(sg.target.clone(), sg.matches);
            }
            res.consensus = best.argmax().unwrap().clone();
            res.element_class = lib.class_of(&res.consensus).unwrap_or("").to_string();
        }
        // --- strand
        if term_strand != 0 && strand_hint.as_ref().is_none_or(|h| h.0 != term_strand) {
            res.strand = term_strand;
            res.strand_source = "element_3p_end".into();
        } else if let Some((st, src)) = strand_hint {
            res.strand = st;
            res.strand_source = src;
        } else {
            let (st, src) = Self::strand_from_segments(&raw_layouts, &res.consensus, lib);
            res.strand = st;
            res.strand_source = src.into();
        }
        let st = if res.strand != 0 { res.strand } else { 1 };
        res.layouts = raw_layouts.iter().map(|l| l.to_sense(st)).collect();
        res.raw_layouts = raw_layouts;
        if !res.consensus.is_empty() {
            res.pileup(asm);
            if st >= 0 {
                // python aliasing: to_sense(+) shares the Segment OBJECTS with raw_layouts, so
                // the pileup's in-place updates show up there too (not for strand < 0)
                for (r, l) in res.raw_layouts.iter_mut().zip(&res.layouts) {
                    r.segments.clone_from(&l.segments);
                }
            }
            res.nearest(asm);
        }
        res
    }

    /// `_terminal_3p(layouts, lib, cfg)` -> (strand or 0, element segments of that strand)
    fn terminal_3p(layouts: &[ReadLayout], lib: &RteLibrary, cfg: &AssemblyCfg) -> (i32, Vec<Segment>) {
        let tol = cfg.terminal_tolerance;
        let mbp = cfg.terminal_min_bp;
        let mid = cfg.terminal_min_identity;
        let rgap = cfg.terminal_max_ref_gap;
        let mut votes = [0usize; 2]; // [+1, -1]
        let mut segs_by: [Vec<Segment>; 2] = [Vec::new(), Vec::new()];
        let term_ok = |e: &Segment, strand: i32| -> bool {
            if e.kind != SegKind::Element || e.strand != strand || e.qlen() < mbp || e.identity < mid {
                return false;
            }
            let cend = lib.cons_end.get(&e.target).map(|&x| x as i64).unwrap_or(0);
            cend > 0 && e.t_en >= cend - tol
        };
        for lay in layouts {
            if !counts_as_fragment(lay) {
                continue;
            }
            let sg = &lay.segments;
            for i in 0..sg.len().saturating_sub(2) {
                let (a, b, c) = (&sg[i], &sg[i + 1], &sg[i + 2]);
                let bt = if b.target.is_empty() { "A" } else { b.target.as_str() };
                if b.kind == SegKind::PolyA && bt == "A" && c.kind == SegKind::Ref && term_ok(a, 1) && b.q_st - a.q_en <= 5 && c.q_st - b.q_en <= rgap {
                    votes[0] += 1;
                    segs_by[0].push(a.clone());
                } else if a.kind == SegKind::Ref && b.kind == SegKind::PolyA && b.target == "T" && term_ok(c, -1) && b.q_st - a.q_en <= rgap && c.q_st - b.q_en <= 5 {
                    votes[1] += 1;
                    segs_by[1].push(c.clone());
                }
            }
        }
        for (k, st) in [(0usize, 1i32), (1, -1)] {
            if votes[k] > 0 && votes[1 - k] == 0 {
                return (st, std::mem::take(&mut segs_by[k]));
            }
        }
        (0, Vec::new())
    }

    /// `_strand_from_segments(layouts, consensus, lib)`
    fn strand_from_segments(layouts: &[ReadLayout], consensus: &str, lib: &RteLibrary) -> (i32, &'static str) {
        let (pst, _) = strand_from_polya(layouts, 10);
        if pst != 0 {
            return (pst, "polya_reads");
        }
        let cls_ = if consensus.is_empty() { None } else { lib.class_of(consensus) }.filter(|c| !c.is_empty());
        let mut adj: OrdCount<i32> = OrdCount::new();
        let mut allc: OrdCount<i32> = OrdCount::new();
        for lay in layouts {
            let segs = &lay.segments;
            for (i, s) in segs.iter().enumerate() {
                if s.kind != SegKind::Element || cls_.is_some_and(|c| lib.class_of(&s.target) != Some(c)) {
                    continue;
                }
                allc.add(s.strand, s.qlen());
                if (i > 0 && segs[i - 1].kind == SegKind::Ref) || (i + 1 < segs.len() && segs[i + 1].kind == SegKind::Ref) {
                    adj.add(s.strand, s.qlen());
                }
            }
        }
        for (cnt, src) in [(&adj, "junction_segments"), (&allc, "segments")] {
            if let Some(&s) = cnt.argmax() {
                return (s, src);
            }
        }
        (0, "none")
    }

    /// `_pileup(asm)`: re-place every ELEMENT piece of the chosen class on the chosen consensus
    /// (in place on `self.layouts`), majority sequence per covered interval. Memory is bounded
    /// by the consensus length (votes per consensus position).
    fn pileup(&mut self, asm: &Assembler) {
        let lib = asm.lib;
        let Some(cons) = lib.consensus.get(&self.consensus) else { panic!("KeyError: '{}'", self.consensus) };
        let cend = lib.cons_end.get(&self.consensus).copied().unwrap_or(cons.len()) as i64;
        let cls_ = lib.class_of(&self.consensus);
        let npos = cons.len();
        // pos -> [(base | '-', count)] / [(inserted string ('' = none), count)], insertion order
        let mut votes: Vec<Vec<(u8, i64)>> = vec![Vec::new(); npos];
        let mut ins: Vec<Vec<(Vec<u8>, i64)>> = vec![Vec::new(); npos];
        let (mut ident_m, mut ident_b) = (0i64, 0i64);
        for li in 0..self.layouts.len() {
            for si in 0..self.layouts[li].segments.len() {
                let lay = &self.layouts[li];
                let s = &lay.segments[si];
                if s.kind != SegKind::Element || lib.class_of(&s.target) != cls_ {
                    continue;
                }
                let mut piece = sl(&lay.seq, s.q_st, s.q_en).to_vec();
                let strand_c = s.strand;
                if strand_c < 0 {
                    piece = rc(&piece);
                }
                // re-place the piece on the chosen consensus (exact ops for the pileup)
                let same = s.target == self.consensus;
                let mut lo = if same && s.t_st >= 0 { (s.t_st - 40).max(0) } else { 0 };
                let hi = if same && s.t_en >= 0 { cend.min(s.t_en + 40) } else { cend };
                let k = std::cmp::max(1, (piece.len() as f64 * 0.25) as i64) as i32;
                let mut r = edlib_path(&piece, sl(cons, lo, hi), k);
                if r.is_none() && (lo, hi) != (0, cend) {
                    lo = 0;
                    r = edlib_path(&piece, sl(cons, 0, cend), k);
                }
                let Some((ed, ts, te, ops)) = r else { continue };
                let plen = piece.len() as i64;
                let mut t = lo + ts;
                let mut q = 0i64;
                {
                    // exact coordinates on the chosen consensus from here on
                    let s = &mut self.layouts[li].segments[si];
                    s.target = self.consensus.clone();
                    s.t_st = lo + ts;
                    s.t_en = lo + te;
                    s.identity = 1.0 - ed as f64 / plen.max(1) as f64;
                }
                self.segments_on_cons.push((lo + ts, lo + te, strand_c > 0, li));
                ident_m += plen - ed as i64;
                ident_b += plen.max(te - ts);
                for (ln, op) in ops {
                    let ln = ln as i64;
                    match op {
                        0 => {
                            for kk in 0..ln {
                                let p = (t + kk) as usize;
                                bump(&mut votes[p], piece[(q + kk) as usize], 1);
                                bump(&mut ins[p], Vec::new(), 1);
                            }
                            t += ln;
                            q += ln;
                        }
                        1 => {
                            if t >= 1 {
                                let p = (t - 1) as usize;
                                bump(&mut ins[p], Vec::new(), -1);
                                bump(&mut ins[p], sl(&piece, q, q + ln).to_vec(), 1);
                            }
                            q += ln;
                        }
                        _ => {
                            for kk in 0..ln {
                                bump(&mut votes[(t + kk) as usize], b'-', 1);
                            }
                            t += ln;
                        }
                    }
                }
            }
        }
        self.consensus_identity = if ident_b != 0 { ident_m as f64 / ident_b as f64 } else { 0.0 };
        // merged covered intervals
        let mut ivs: Vec<(i64, i64)> = self.segments_on_cons.iter().map(|x| (x.0, x.1)).collect();
        ivs.sort();
        let mut merged: Vec<(i64, i64)> = Vec::new();
        for (a, b) in ivs {
            match merged.last_mut() {
                Some(m) if a <= m.1 => m.1 = m.1.max(b),
                _ => merged.push((a, b)),
            }
        }
        self.covered = merged;
        let mut seqs = Vec::with_capacity(self.covered.len());
        for &(a, b) in &self.covered {
            let mut out: Vec<u8> = Vec::new();
            for p in a.max(0)..b.min(npos as i64) {
                // Counter.most_common(1) / max(items) -> FIRST maximal entry
                if let Some(&(base, _)) = first_max(&votes[p as usize]) {
                    if base != b'-' {
                        out.push(base);
                    }
                }
                if let Some((best, _)) = first_max(&ins[p as usize]) {
                    out.extend_from_slice(best);
                }
            }
            seqs.push(out);
        }
        self.covered_seqs = seqs;
    }

    /// `_nearest(asm)`: identity of the covered parts to the intact elements (per contig the
    /// first hit only; key (identity, blen), FIRST max).
    fn nearest(&mut self, asm: &Assembler) {
        let lib = asm.lib;
        let Some(al) = lib.aligner(AlignerKind::Intact) else { return };
        let mut acc: Vec<(String, i64, i64)> = Vec::new();
        for s in &self.covered_seqs {
            if s.len() < 30 {
                continue;
            }
            let hits = al.map(s);
            let mut seen: Vec<&str> = Vec::new();
            for h in &hits {
                if seen.contains(&&*h.ctg) {
                    continue;
                }
                seen.push(&h.ctg);
                match acc.iter_mut().find(|e| e.0 == *h.ctg) {
                    Some(e) => {
                        e.1 += h.mlen;
                        e.2 += h.blen;
                    }
                    None => acc.push((h.ctg.to_string(), h.mlen, h.blen)),
                }
            }
        }
        if acc.is_empty() {
            return;
        }
        fn key(e: &(String, i64, i64)) -> (f64, i64) {
            (if e.2 != 0 { e.1 as f64 / e.2 as f64 } else { 0.0 }, e.2)
        }
        fn best_of<'x>(items: impl Iterator<Item = &'x (String, i64, i64)>) -> Option<(String, f64)> {
            let mut best: Option<&(String, i64, i64)> = None;
            for e in items {
                let (ki, kb) = key(e);
                if best.is_none_or(|b| {
                    let (bi, bb) = key(b);
                    ki > bi || (ki == bi && kb > bb)
                }) {
                    best = Some(e);
                }
            }
            best.map(|b| (b.0.clone(), key(b).0))
        }
        let (bi, bid) = best_of(acc.iter()).unwrap();
        self.nearest_intact = bi;
        self.nearest_intact_identity = py_round(bid, 4);
        if let Some((ba, bad)) = best_of(acc.iter().filter(|e| lib.is_active_element(&e.0))) {
            self.nearest_active = ba;
            self.element_identity = py_round(bad, 4);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::golden;
    use serde_json::Value;
    use std::path::{Path, PathBuf};

    fn events() -> Vec<Value> {
        let mut v = golden::load_events(&golden::golden_dir().join("pytest.jsonl.gz"));
        if let Some(m) = golden::e2e_manifest() {
            v.extend(golden::load_events(Path::new(m["golden"].as_str().unwrap())));
        }
        v.into_iter().filter(|e| e["kind"] == "assemble").collect()
    }

    fn first_diff(got: &AssemblyResult, want: &AssemblyResult) -> String {
        macro_rules! chk {
            ($f:ident) => {
                if got.$f != want.$f {
                    return format!("{}: got {:?} want {:?}", stringify!($f), got.$f, want.$f);
                }
            };
        }
        chk!(strand);
        chk!(strand_source);
        chk!(element_class);
        chk!(consensus);
        chk!(element_bp);
        chk!(class_bp);
        for (nm, g, w) in [("raw", &got.raw_layouts, &want.raw_layouts), ("sense", &got.layouts, &want.layouts)] {
            if g.len() != w.len() {
                return format!("{nm} layouts len {} vs {}", g.len(), w.len());
            }
            for (i, (a, b)) in g.iter().zip(w).enumerate() {
                if a != b {
                    let f = |l: &ReadLayout| -> Vec<_> { l.segments.iter().map(|s| (s.q_st, s.q_en, s.kind.as_str(), s.target.clone(), s.t_st, s.t_en, s.strand, s.identity, s.matches)).collect() };
                    return format!("{nm} layout {i} {} (seq eq {})\n got  {:?}\n want {:?}", a.name, a.seq == b.seq, f(a), f(b));
                }
            }
        }
        chk!(segments_on_cons);
        chk!(covered);
        chk!(covered_seqs);
        chk!(consensus_identity);
        chk!(nearest_intact);
        chk!(nearest_intact_identity);
        chk!(nearest_active);
        chk!(element_identity);
        "?".into()
    }

    /// Field-level report over every recorded `assemble` call (diagnostics; the gate is
    /// tests/golden_modules.rs::golden_assemble).
    #[test]
    #[ignore = "diagnostic"]
    fn assemble_report() {
        let mut libs: Vec<(PathBuf, RteLibrary)> = Vec::new();
        let (mut ok, mut bad) = (0, 0);
        for e in events() {
            let dir = golden::lib_dir(&e).unwrap_or_else(|| golden::repo_path("test/fixtures/rte_library"));
            if !libs.iter().any(|(k, _)| *k == dir) {
                libs.push((dir.clone(), RteLibrary::open(&dir.to_string_lossy(), None).unwrap()));
            }
            let lib = &libs.iter().find(|(k, _)| *k == dir).unwrap().1;
            let cfg = AssemblyCfg::from(e["cfg"].as_object()).unwrap();
            let asm = Assembler::new(lib, cfg);
            let ctx = golden::ctx(&e["in"]["ctx"]);
            let js: Vec<(String, JunctionSeq)> = e["in"]["junction_seqs"].as_array().unwrap().iter().map(|x| (golden::s(&x[0]), (golden::b(&x[1]), (golden::i(&x[2][0]), golden::i(&x[2][1]))))).collect();
            let reads: Vec<_> = e["in"]["reads"].as_array().unwrap().iter().map(golden::read).collect();
            let hint = e["in"]["strand_hint"].as_array().map(|h| (golden::i(&h[0]) as i32, golden::s(&h[1])));
            let mut got = asm.assemble(&ctx, &js, &reads, hint);
            golden::strip_cigars(&mut got);
            let want = golden::assembly(&e["out"]);
            if got == want {
                ok += 1;
            } else {
                bad += 1;
                if bad <= 8 {
                    eprintln!("MISMATCH {} seq {}: {}", e["case"], e["seq"], first_diff(&got, &want));
                }
            }
        }
        eprintln!("assemble: {ok} ok, {bad} mismatching");
        assert_eq!(bad, 0);
    }
}
