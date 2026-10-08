//! tools/rte/transduction.py -- 3' transduction source lookup + the novel-source rule.
//!
//! STATUS: SourceCall, known_source, MappyLocator = FOUNDATION (implemented: structure.rs needs
//! known_source, and the locator is a thin mm.rs wrapper). L1Rmsk, cons_identity and
//! NovelSourceFinder = WORK PACKAGE "WP-TD" (PORT_PLAN.md).
//!
//! Python functions to port: L1Rmsk.__init__ / upstream_of (112-170), cons_identity (59) =
//! tools/rte_library/common.cons_identity -> common.identity (edlib HW path, columns = M+X+I+D,
//! identity = matches / columns; max of both directions), NovelSourceFinder.available /
//! _l1_identity / find (192-277).
//!
//! Golden: events `novel_find` (in: seq, the locator's answer, cohort_l1, cfg, rmsk/genome
//! paths; out: SourceCall|null), `cons_identity` (in: seq, cons; out: float),
//! `l1rmsk_upstream_of` (in: contig, s, e, strand, max_dist; out: rows).

use crate::assembly::Segment;
use crate::config::TransductionCfg;
use crate::genome::{parse_region, Genome};
use crate::library::RteLibrary;
use crate::mm::{Aligner, MapOpts};
use crate::structure::SourceFinder;
use std::path::Path;

/// transduction.SourceCall
#[derive(Clone, Debug, PartialEq)]
pub struct SourceCall {
    pub source_id: String,
    /// offset of the transduction endpoint inside the source flank
    pub td_end: i64,
    /// offset of the tag start (0 = right after the element 3' end)
    pub td_start: i64,
    pub n_segments: i64,
    pub novel: bool,
    /// novel: identity of the source L1 to the L1HS consensus (rounded to 4)
    pub identity: f64,
    pub detail: String,
    /// novel: "A" / "B"; "" for known sources
    pub tier: String,
}

/// `known_source(flank_segments, lib)`: the source with the most aligned bases among
/// sense-oriented FLANK3P hits (python Counter + `max(per, key=per.get)` = first maximal in
/// insertion order).
pub fn known_source(flank_segments: &[&Segment], lib: &RteLibrary) -> Option<SourceCall> {
    let mut per: Vec<(String, i64)> = Vec::new();
    let mut ends: Vec<(String, i64)> = Vec::new();
    let mut starts: Vec<(String, i64)> = Vec::new();
    for s in flank_segments {
        if s.strand < 0 {
            continue;
        }
        let sid = lib.source_for_flank(&s.target);
        match per.iter_mut().find(|(k, _)| *k == sid) {
            Some(e) => e.1 += s.matches,
            None => per.push((sid.clone(), s.matches)),
        }
        match ends.iter_mut().find(|(k, _)| *k == sid) {
            Some(e) => e.1 = e.1.max(s.t_en),
            None => ends.push((sid.clone(), (-1i64).max(s.t_en))),
        }
        match starts.iter_mut().find(|(k, _)| *k == sid) {
            Some(e) => e.1 = e.1.min(s.t_st),
            None => starts.push((sid, 1_000_000_000i64.min(s.t_st))),
        }
    }
    let mut best: Option<&(String, i64)> = None;
    for e in &per {
        if best.is_none_or(|b| e.1 > b.1) {
            best = Some(e);
        }
    }
    let sid = best?.0.clone();
    let hits: Vec<&&Segment> = flank_segments.iter().filter(|s| lib.source_for_flank(&s.target) == sid && s.strand > 0).collect();
    let mut strands: Vec<String> = hits.iter().map(|s| lib.flank_strand(&s.target)).filter(|x| !x.is_empty()).collect();
    strands.sort();
    strands.dedup();
    let detail = if strands.is_empty() { String::new() } else { format!("source_strand={}", strands.join(",")) };
    let get = |v: &Vec<(String, i64)>| v.iter().find(|(k, _)| *k == sid).map(|e| e.1).unwrap();
    Some(SourceCall {
        td_end: get(&ends),
        td_start: get(&starts),
        n_segments: hits.len() as i64,
        source_id: sid,
        novel: false,
        identity: 0.0,
        detail,
        tier: String::new(),
    })
}

/// One `MappyLocator` placement: (contig, start, end, '+'/'-', mapq, identity).
#[derive(Clone, Debug, PartialEq)]
pub struct LocatorHit {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub strand: char,
    pub mapq: i64,
    pub identity: f64,
}

/// `locator(seq)` -> primary placements on the remap genome.
pub trait Locator: Sync {
    fn locate(&self, seq: &[u8]) -> Vec<LocatorHit>;
}

/// transduction.MappyLocator: minimap2 `preset="sr"` on a .mmi / FASTA; primary hits only; a
/// region-FASTA record `contig:start-end` is shifted to genome coordinates.
pub struct MappyLocator {
    al: Aligner,
}

impl MappyLocator {
    pub fn open(path: &Path) -> Option<MappyLocator> {
        Aligner::from_path(path, MapOpts::SR).map(|al| MappyLocator { al })
    }
}

impl Locator for MappyLocator {
    fn locate(&self, seq: &[u8]) -> Vec<LocatorHit> {
        let mut out = Vec::new();
        for h in self.al.map(seq) {
            if !h.is_primary {
                continue;
            }
            let (ctg, off) = match parse_region(&h.ctg) {
                Some((c, s, _)) => (c, s),
                None => (h.ctg.to_string(), 0),
            };
            out.push(LocatorHit {
                contig: ctg,
                start: off + h.r_st,
                end: off + h.r_en,
                strand: if h.strand == 1 { '+' } else { '-' },
                mapq: h.mapq as i64,
                identity: h.mlen as f64 / h.blen.max(1) as f64,
            });
        }
        out
    }
}

/// One young L1: (start, end, strand '+'/'-', name, divergence %).
pub type L1Row = (i64, i64, char, String, f64);

/// transduction.L1Rmsk: young full-length L1s (LINE/L1, >= min_len) per contig. WP-TD.
pub struct L1Rmsk {
    /// contig -> sorted rows
    pub by_contig: rustc_hash::FxHashMap<String, Vec<L1Row>>,
}

impl L1Rmsk {
    /// Parse a RepeatMasker .out or UCSC rmsk.txt (plain or gzipped). WP-TD.
    pub fn open(_path: &Path, _min_len: i64) -> Result<L1Rmsk, String> {
        todo!("WP-TD: port transduction.L1Rmsk.__init__")
    }

    /// `upstream_of(contig, s, e, strand, max_dist)`. WP-TD.
    pub fn upstream_of(&self, _contig: &str, _s: i64, _e: i64, _strand: char, _max_dist: i64) -> Vec<L1Row> {
        todo!("WP-TD: port transduction.L1Rmsk.upstream_of")
    }
}

/// `cons_identity(seq, cons)` (tools/rte_library/common.cons_identity). WP-TD.
pub fn cons_identity(_seq: &[u8], _cons: &[u8]) -> f64 {
    todo!("WP-TD: port tools/rte_library/common.cons_identity")
}

/// transduction.NovelSourceFinder. `cohort_l1` is set by the annotator's cohort pass.
pub struct NovelSourceFinder<'a> {
    pub lib: &'a RteLibrary,
    pub cfg: TransductionCfg,
    pub rmsk: Option<L1Rmsk>,
    pub locator: Option<Box<dyn Locator + 'a>>,
    pub genome: Option<&'a dyn Genome>,
    /// cohort L1 calls on the remap genome: (contig, pos, '+'/'-')
    pub cohort_l1: Vec<(String, i64, char)>,
    /// `_l1_identity` cache ((contig, s, e) -> identity); shared across threads
    pub ident_cache: std::sync::Mutex<rustc_hash::FxHashMap<(String, i64, i64), f64>>,
}

impl NovelSourceFinder<'_> {
    /// `available()`. WP-TD.
    pub fn available(&self) -> bool {
        todo!("WP-TD: port NovelSourceFinder.available")
    }
}

impl SourceFinder for NovelSourceFinder<'_> {
    /// `find(seq)`. WP-TD.
    fn find(&self, _seq: &[u8]) -> Option<SourceCall> {
        todo!("WP-TD: port NovelSourceFinder.find")
    }
}
