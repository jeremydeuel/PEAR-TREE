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
use crate::align::{self, Mode};
use crate::mm::{Aligner, MapOpts};
use crate::pyfmt::py_round;
use crate::sequtil::rc;
use crate::structure::SourceFinder;
use std::io::BufRead;
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

fn row_cmp(a: &L1Row, b: &L1Row) -> std::cmp::Ordering {
    a.0.cmp(&b.0)
        .then(a.1.cmp(&b.1))
        .then(a.2.cmp(&b.2))
        .then_with(|| a.3.cmp(&b.3))
        .then_with(|| a.4.partial_cmp(&b.4).unwrap_or(std::cmp::Ordering::Equal))
}

/// transduction.L1Rmsk: young full-length L1s (LINE/L1, >= min_len) per contig.
pub struct L1Rmsk {
    /// contig -> sorted rows
    pub by_contig: rustc_hash::FxHashMap<String, Vec<L1Row>>,
}

fn pyint(s: &str) -> Result<i64, String> {
    s.trim().parse::<i64>().map_err(|_| format!("invalid literal for int(): '{s}'"))
}

fn pyfloat(s: &str) -> Result<f64, String> {
    s.trim().parse::<f64>().map_err(|_| format!("could not convert string to float: '{s}'"))
}

impl L1Rmsk {
    /// Parse a RepeatMasker .out or UCSC rmsk.txt (plain or gzipped).
    pub fn open(path: &Path, min_len: i64) -> Result<L1Rmsk, String> {
        let rd = crate::inputs::open_text(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        let mut by_contig: rustc_hash::FxHashMap<String, Vec<L1Row>> = Default::default();
        for line in rd.lines() {
            let line = line.map_err(|e| format!("{}: {e}", path.display()))?;
            let f: Vec<&str> = line.split_whitespace().collect();
            let (contig, s, e, strand, name, fam, div);
            if f.len() >= 15 && f[0].bytes().all(|b| b.is_ascii_digit()) && f.len() < 17 {
                // RepeatMasker .out: score div del ins contig begin end (left) strand name class/family ...
                contig = f[4];
                s = pyint(f[5])? - 1;
                e = pyint(f[6])?;
                strand = f[8];
                name = f[9];
                fam = f[10].to_string();
                div = pyfloat(f[1])?;
            } else if f.len() >= 17 {
                // UCSC rmsk.txt: bin swScore milliDiv ... genoName genoStart genoEnd genoLeft strand repName repClass repFamily
                contig = f[5];
                s = pyint(f[6])?;
                e = pyint(f[7])?;
                strand = f[9];
                name = f[10];
                fam = format!("{}/{}", f[11], f[12]);
                div = pyint(f[2])? as f64 / 10.0;
            } else {
                continue;
            }
            if !fam.starts_with("LINE/L1") || e - s < min_len {
                continue;
            }
            let strand = if strand == "C" || strand == "-" { '-' } else { '+' };
            by_contig.entry(contig.to_string()).or_default().push((s, e, strand, name.to_string(), div));
        }
        for v in by_contig.values_mut() {
            v.sort_by(row_cmp);
        }
        Ok(L1Rmsk { by_contig })
    }

    /// `upstream_of(contig, s, e, strand, max_dist)`: L1s on `strand` whose 3' end lies within
    /// max_dist upstream (in that strand's sense) of the segment [s, e).
    pub fn upstream_of(&self, contig: &str, s: i64, e: i64, strand: char, max_dist: i64) -> Vec<L1Row> {
        let alt = match contig.strip_prefix("chr") {
            Some(x) => x.to_string(),
            None => format!("chr{contig}"),
        };
        let Some(lst) = self.by_contig.get(contig).or_else(|| self.by_contig.get(&alt)) else { return Vec::new() };
        let target = s - max_dist - 10000;
        // bisect_left over the start column
        let lo = lst.partition_point(|r| r.0 < target);
        let mut out: Vec<(i64, &L1Row)> = Vec::new();
        for r in &lst[lo..] {
            let (ls, le, lstr) = (r.0, r.1, r.2);
            if ls > e + max_dist {
                break;
            }
            if lstr != strand {
                continue;
            }
            if strand == '+' && le <= s + 50 && s - le <= max_dist {
                out.push((s - le, r));
            } else if strand == '-' && ls >= e - 50 && ls - e <= max_dist {
                out.push((ls - e, r));
            }
        }
        out.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| row_cmp(a.1, b.1)));
        out.into_iter().map(|(_, r)| r.clone()).collect()
    }
}

/// common.identity(a, b, mode="HW")[0]: matches / (M + X + I + D) of the edlib path alignment of
/// `a` as an infix of `b` (both upper-cased).
fn hw_identity(a: &[u8], b: &[u8]) -> f64 {
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let Some(r) = align::path(&a, &b, Mode::Hw, -1, &[]) else { return 0.0 };
    let m = r.ops.iter().filter(|&&op| op == align::OP_MATCH).count();
    let cols = r.ops.len();
    if cols > 0 {
        m as f64 / cols as f64
    } else {
        0.0
    }
}

/// `cons_identity(seq, cons)` (tools/rte_library/common.cons_identity): the better of the
/// consensus as an infix of the element and the element as an infix of the consensus.
pub fn cons_identity(seq: &[u8], cons: &[u8]) -> f64 {
    let a = hw_identity(cons, seq);
    let b = hw_identity(seq, cons);
    if b > a {
        b
    } else {
        a
    }
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
    /// `available()`.
    pub fn available(&self) -> bool {
        self.locator.is_some() && (self.rmsk.is_some() || !self.cohort_l1.is_empty() || !self.lib.polymorphic_l1.is_empty())
    }

    fn l1_identity(&self, contig: &str, s: i64, e: i64, strand: char, div: f64) -> f64 {
        let key = (contig.to_string(), s, e);
        if let Some(v) = self.ident_cache.lock().unwrap().get(&key) {
            return *v;
        }
        let mut ident: Option<f64> = None;
        let cons = self.lib.consensus.get("L1HS").filter(|c| !c.is_empty());
        if let (Some(g), Some(cons)) = (self.genome, cons) {
            let mut seq = g.fetch(contig, s, e);
            if !seq.is_empty() {
                if strand == '-' {
                    seq = rc(&seq);
                }
                ident = Some(cons_identity(&seq.to_ascii_uppercase(), cons));
            }
        }
        let ident = ident.unwrap_or(1.0 - div / 100.0);
        self.ident_cache.lock().unwrap().insert(key, ident);
        ident
    }
}

impl SourceFinder for NovelSourceFinder<'_> {
    /// `find(seq)`: seq = the transduced segment in element-sense orientation.
    fn find(&self, seq: &[u8]) -> Option<SourceCall> {
        let c = &self.cfg;
        if !self.available() || (seq.len() as i64) < c.novel_source_min_seg {
            return None;
        }
        let hits: Vec<LocatorHit> = self.locator.as_ref()?.locate(seq).into_iter().filter(|h| h.mapq >= c.novel_source_min_mapq).collect();
        if hits.len() != 1 {
            return None;
        }
        let h = &hits[0];
        let (contig, s, e, gstrand) = (h.contig.as_str(), h.start, h.end, h.strand);
        if h.identity < c.novel_source_min_tag_identity {
            return None;
        }
        // (dist, src, ident, name)
        let mut best: Option<(i64, String, f64, String)> = None;
        if let Some(rmsk) = &self.rmsk {
            for (ls, le, lstr, name, div) in rmsk.upstream_of(contig, s, e, gstrand, c.novel_source_max_dist) {
                let ident = self.l1_identity(contig, ls, le, lstr, div);
                if ident >= c.novel_source_min_identity {
                    let dist = if gstrand == '+' { s - le } else { ls - e };
                    best = Some((dist, format!("{contig}:{ls}-{le}"), ident, name));
                    break;
                }
            }
        }
        if best.is_none() {
            for (cc, pos, cstr) in &self.cohort_l1 {
                if cc != contig || *cstr != gstrand {
                    continue;
                }
                let dist = if gstrand == '+' { s - pos } else { pos - e };
                if (0..=c.novel_source_max_dist).contains(&dist) {
                    best = Some((dist, format!("{contig}:{pos}"), 1.0, "cohort_L1".to_string()));
                    break;
                }
            }
        }
        let mut poly = false;
        if best.is_none() {
            // tier B fallback: a polymorphic L1 position in the upstream window, either orientation
            for (cc, pos, pid) in &self.lib.polymorphic_l1 {
                if cc != contig {
                    continue;
                }
                let dist = if gstrand == '+' { s - pos } else { pos - e };
                if (0..=c.novel_source_max_dist).contains(&dist) && best.as_ref().is_none_or(|b| dist < b.0) {
                    best = Some((dist, format!("{contig}:{pos}"), 0.0, format!("polymorphic_L1:{pid}")));
                    poly = true;
                }
            }
        }
        let (dist, src, ident, name) = best?;
        let off_end = dist + (e - s);
        let tier = if poly || ident < c.novel_source_tier_a { "B" } else { "A" };
        Some(SourceCall {
            source_id: format!("novel:{src}"),
            td_end: off_end,
            td_start: dist,
            n_segments: 1,
            novel: true,
            identity: py_round(ident, 4),
            detail: format!("source={name};tag={contig}:{s}-{e}({gstrand});dist={dist};tier={tier}"),
            tier: tier.to_string(),
        })
    }
}
