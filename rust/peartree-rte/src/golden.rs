//! Golden-event support for tests (golden/make_golden.py output). FOUNDATION (implemented).
//!
//! A golden file is JSON lines; each line is one event recorded at a python module boundary:
//!   {"kind": K, "case": <pytest node id or e2e locus>, "lib": <library dir>|null,
//!    "genome": <genome id>|null, "remap": <genome id>|null, "cfg": {...CONFIG['annotate']...},
//!    "in": {...}, "out": ...}
//! plus `{"kind": "genome", "id": .., "path": ..}` / `{"kind": "genome", "id": .., "regions":
//! [[contig, offset, seq], ...]}` lines defining the genomes the other events refer to.
//! SPEC.md "Golden harness" lists every kind and its in/out fields; the converters below turn
//! the JSON into the crate's types (the field names are python's).

use crate::assembly::{AssemblyResult, ReadLayout, SegKind, Segment, SiteContext};
use crate::hallmarks::{PolyAInfo, SiteInfo};
use crate::genome::{open_genome, FastaGenome, Genome};
use crate::inputs::{EvidenceRead, JunctionEvidence};
use crate::packed::PackedSeq;
use crate::record::{Detail, RteRecord};
use crate::score::ScoreInput;
use crate::structure::{SourceFinder, StructureCall};
use crate::transduction::SourceCall;
use serde_json::Value;
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// Directory of the committed golden data (override with PEARTREE_RTE_GOLDEN).
pub fn golden_dir() -> PathBuf {
    std::env::var_os("PEARTREE_RTE_GOLDEN")
        .map(PathBuf::from)
        .unwrap_or_else(|| Path::new(env!("CARGO_MANIFEST_DIR")).join("golden/data"))
}

/// A recorded path: relative ones are relative to the repository root (recorder.rel).
pub fn repo_path(p: &str) -> PathBuf {
    let pb = PathBuf::from(p);
    if pb.is_absolute() {
        pb
    } else {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../..").join(pb)
    }
}

/// The library directory of an event (`lib` field), None for a test stand-in (`{"mock": ..}`).
pub fn lib_dir(ev: &Value) -> Option<PathBuf> {
    ev["lib"].as_str().map(repo_path)
}

/// Every event of a golden file (plain or .gz).
pub fn load_events(path: &Path) -> Vec<Value> {
    let rdr = crate::inputs::open_text(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    rdr.lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(&l).unwrap())
        .collect()
}

/// Events of one kind from every `*.jsonl[.gz]` in `golden_dir()` (empty when none).
pub fn events_of(kind: &str) -> Vec<Value> {
    let mut out = Vec::new();
    let Ok(rd) = std::fs::read_dir(golden_dir()) else { return out };
    let mut files: Vec<PathBuf> = rd.filter_map(|e| e.ok().map(|e| e.path())).collect();
    files.sort();
    for f in files {
        let n = f.to_string_lossy().to_string();
        if n.ends_with(".jsonl") || n.ends_with(".jsonl.gz") {
            out.extend(load_events(&f).into_iter().filter(|e| e["kind"] == kind));
        }
    }
    out
}

pub fn s(v: &Value) -> String {
    v.as_str().unwrap_or("").to_string()
}
pub fn b(v: &Value) -> Vec<u8> {
    v.as_str().unwrap_or("").as_bytes().to_vec()
}
pub fn oi(v: &Value) -> Option<i64> {
    v.as_i64().or_else(|| v.as_f64().map(|f| f as i64))
}
pub fn i(v: &Value) -> i64 {
    oi(v).unwrap_or(0)
}
pub fn f(v: &Value) -> f64 {
    v.as_f64().unwrap_or(0.0)
}
pub fn bo(v: &Value) -> bool {
    v.as_bool().unwrap_or(false)
}
pub fn strs(v: &Value) -> Vec<String> {
    v.as_array().map(|a| a.iter().map(s).collect()).unwrap_or_default()
}

pub fn segment(v: &Value) -> Segment {
    Segment {
        q_st: i(&v["q"][0]),
        q_en: i(&v["q"][1]),
        kind: SegKind::parse(v["kind"].as_str().unwrap()).unwrap(),
        target: s(&v["target"]),
        t_st: i(&v["t"][0]),
        t_en: i(&v["t"][1]),
        strand: i(&v["strand"]) as i32,
        identity: f(&v["identity"]),
        matches: i(&v["matches"]),
        cigar: None,
    }
}

pub fn layout(v: &Value) -> ReadLayout {
    ReadLayout {
        name: s(&v["name"]),
        side: s(&v["side"]),
        role: s(&v["role"]),
        frag_key: (s(&v["frag"][0]), s(&v["frag"][1])),
        seq: b(&v["seq"]),
        segments: v["segs"].as_array().unwrap().iter().map(segment).collect(),
        sense: bo(&v["sense"]),
    }
}

pub fn ctx(v: &Value) -> SiteContext {
    let mut c = SiteContext::new(&s(&v["title"]), v["contig"].as_str().map(String::from), oi(&v["left_bp"]), oi(&v["right_bp"]), b(&v["left_flank"]), b(&v["right_flank"]));
    c.window_start = i(&v["window_start"]);
    c.window_seq = b(&v["window_seq"]);
    c.wide_start = i(&v["wide_start"]);
    c.wide_seq = b(&v["wide_seq"]);
    c
}

pub fn assembly(v: &Value) -> AssemblyResult {
    let lays = |k: &str| v[k].as_array().unwrap().iter().map(layout).collect::<Vec<_>>();
    AssemblyResult {
        strand: i(&v["strand"]) as i32,
        strand_source: s(&v["strand_source"]),
        element_class: s(&v["element_class"]),
        consensus: s(&v["consensus"]),
        element_bp: i(&v["element_bp"]),
        layouts: lays("layouts"),
        raw_layouts: lays("raw_layouts"),
        covered: v["covered"].as_array().unwrap().iter().map(|x| (i(&x[0]), i(&x[1]))).collect(),
        covered_seqs: v["covered_seqs"].as_array().unwrap().iter().map(b).collect(),
        segments_on_cons: v["segments_on_cons"].as_array().unwrap().iter().map(|x| (i(&x[0]), i(&x[1]), bo(&x[2]), i(&x[3]) as usize)).collect(),
        consensus_identity: f(&v["consensus_identity"]),
        nearest_intact: s(&v["nearest_intact"]),
        nearest_intact_identity: f(&v["nearest_intact_identity"]),
        nearest_active: s(&v["nearest_active"]),
        element_identity: f(&v["element_identity"]),
        class_bp: v["class_bp"].as_object().unwrap().iter().map(|(k, x)| (k.clone(), i(x))).collect(),
    }
}

/// `[side, role, sample, frag, r12, seq]`
pub fn read(v: &Value) -> EvidenceRead {
    EvidenceRead {
        side: s(&v[0]).into(),
        role: s(&v[1]).into(),
        sample: s(&v[2]).into(),
        frag: s(&v[3]).into(),
        r12: s(&v[4]).into(),
        seq: PackedSeq::encode(&b(&v[5])),
    }
}

pub fn junction(v: &Value) -> JunctionEvidence {
    JunctionEvidence {
        side: s(&v["side"]),
        n_reads: i(&v["n_reads"]),
        n_fragments: i(&v["n_fragments"]),
        n_independent: i(&v["n_independent"]),
        n_samples: i(&v["n_samples"]),
        n_mates: i(&v["n_mates"]),
        supported: i(&v["supported"]),
        clip_consensus: s(&v["clip_consensus"]),
        polya_len_median: f(&v["polya_len_median"]),
        polya_len_range: s(&v["polya_len_range"]),
        beyond_polya: s(&v["beyond_polya"]),
        beyond_polya_support: i(&v["beyond_polya_support"]),
        cross_sample_identical: i(&v["cross_sample_identical"]),
        n_short_used: i(&v["n_short_used"]),
        n_short_mate_inside: i(&v["n_short_mate_inside"]),
    }
}

pub fn source_call(v: &Value) -> Option<SourceCall> {
    if v.is_null() {
        return None;
    }
    Some(SourceCall {
        source_id: s(&v["source_id"]),
        td_end: i(&v["td_end"]),
        td_start: i(&v["td_start"]),
        n_segments: i(&v["n_segments"]),
        novel: bo(&v["novel"]),
        identity: f(&v["identity"]),
        detail: s(&v["detail"]),
        tier: s(&v["tier"]),
    })
}

pub fn structure_call(v: &Value) -> StructureCall {
    StructureCall {
        element: s(&v["element"]),
        structure: s(&v["structure"]),
        tags: strs(&v["tags"]),
        detail: Detail::from_json(&v["detail"]),
        source: source_call(&v["source"]),
        j5_class: s(&v["j5_class"]),
        j3_class: s(&v["j3_class"]),
        j5_pos: oi(&v["j5_pos"]),
        inv_p1: oi(&v["inv_p1"]),
        three_prime_truncated: bo(&v["three_prime_truncated"]),
        three_prime_short: bo(&v["three_prime_short"]),
        has_polya_3p: bo(&v["has_polya_3p"]),
        td3p_seq: b(&v["td3p_seq"]),
    }
}

pub fn score_input(v: &Value) -> ScoreInput {
    ScoreInput {
        element: s(&v["element"]),
        structure: s(&v["structure"]),
        tags: strs(&v["tags"]),
        tsd_len: oi(&v["tsd_len"]),
        tsd_verified: bo(&v["tsd_verified"]),
        polya_len: f(&v["polya_len"]),
        polya_both_sides: bo(&v["polya_both_sides"]),
        slippage: bo(&v["slippage"]),
        beyond_polya_len: i(&v["beyond_polya_len"]),
        beyond_polya_support: i(&v["beyond_polya_support"]),
        en_mismatches: oi(&v["en_mismatches"]),
        ends_concordant: bo(&v["ends_concordant"]),
        td_source_matches_5p: bo(&v["td_source_matches_5p"]),
        element_identity: f(&v["element_identity"]),
        inactive_only: bo(&v["inactive_only"]),
        inv_p1: oi(&v["inv_p1"]),
        junctions_supported: i(&v["junctions_supported"]),
        n_samples: i(&v["n_samples"]),
        foldback: bo(&v["foldback"]),
        recurrent: bo(&v["recurrent"]),
        cross_sample_identical: bo(&v["cross_sample_identical"]),
        novel_tier: s(&v["novel_tier"]),
    }
}

pub fn record(v: &Value) -> RteRecord {
    let mut r = RteRecord::new(&s(&v["insertion_id"]));
    r.element = s(&v["element"]);
    r.structure = s(&v["structure"]);
    r.tags = strs(&v["tags"]);
    r.covered_5p = oi(&v["covered_5p"]);
    r.covered_3p = oi(&v["covered_3p"]);
    r.covered_intervals = v["covered_intervals"].as_array().map(|a| a.iter().map(|x| (i(&x[0]), i(&x[1]))).collect()).unwrap_or_default();
    r.consensus = s(&v["consensus"]);
    r.element_identity = v["element_identity"].as_f64();
    r.nearest_active = s(&v["nearest_active"]);
    r.tsd_seq = s(&v["tsd_seq"]);
    r.tsd_len = oi(&v["tsd_len"]);
    r.en_motif = s(&v["en_motif"]);
    r.en_mismatches = oi(&v["en_mismatches"]);
    r.polya_len = v["polya_len"].as_f64();
    r.beyond_polya = s(&v["beyond_polya"]);
    r.beyond_polya_support = i(&v["beyond_polya_support"]);
    r.strand = i(&v["strand"]) as i32;
    r.tprt_score = f(&v["tprt_score"]);
    r.tprt_points = s(&v["tprt_points"]);
    r.tprt_call = s(&v["tprt_call"]);
    r.detail = Detail::from_json(&v["detail"]);
    r.score_input = if v["score_input"].is_null() { None } else { Some(score_input(&v["score_input"])) };
    r.site = (v["site"][0].as_str().map(String::from), oi(&v["site"][1]), oi(&v["site"][2]));
    r.gt_reads = i(&v["gt_reads"]);
    r.gt_changed = s(&v["gt_changed"]);
    r
}

/// Genomes defined by `genome` events, by id.
#[derive(Default)]
pub struct Genomes {
    by_id: rustc_hash::FxHashMap<String, Box<dyn Genome>>,
    missing: rustc_hash::FxHashSet<String>,
}

impl Genomes {
    /// Register every `genome` event of `events`. A path that does not exist locally is noted
    /// (tests needing it skip).
    pub fn from_events(events: &[Value]) -> Genomes {
        let mut g = Genomes::default();
        for e in events.iter().filter(|e| e["kind"] == "genome") {
            let id = s(&e["id"]);
            if let Some(p) = e["path"].as_str() {
                let p = repo_path(p);
                if p.exists() {
                    g.by_id.insert(id, open_genome(Some(&p.to_string_lossy())).unwrap().unwrap());
                } else {
                    g.missing.insert(id);
                }
            } else {
                let mut fg = FastaGenome::default();
                for r in e["regions"].as_array().unwrap() {
                    let (c, off, seq) = (s(&r[0]), i(&r[1]), b(&r[2]));
                    let name = if off == 0 { c } else { format!("{c}:{off}-{}", off + seq.len() as i64) };
                    fg.add(&name, &seq);
                }
                g.by_id.insert(id, Box::new(fg));
            }
        }
        g
    }
    /// Event field value (genome id or null) -> genome. Err(()) when the genome is missing
    /// locally (the caller skips the event).
    #[allow(clippy::result_unit_err)]
    pub fn get(&self, id: &Value) -> Result<Option<&dyn Genome>, ()> {
        match id.as_str() {
            None => Ok(None),
            Some(k) if self.missing.contains(k) => Err(()),
            Some(k) => Ok(Some(self.by_id.get(k).ok_or(())?.as_ref())),
        }
    }
}

/// A SourceFinder replaying the answers recorded with a `classify` event
/// (`in.novel_answers` = [[seq, SourceCall|null], ...]). Panics on an unrecorded query: the port
/// asked something python did not.
pub struct ReplayFinder(pub Vec<(Vec<u8>, Option<SourceCall>)>);

impl ReplayFinder {
    pub fn from_json(v: &Value) -> ReplayFinder {
        ReplayFinder(v.as_array().map(|a| a.iter().map(|x| (b(&x[0]), source_call(&x[1]))).collect()).unwrap_or_default())
    }
}

impl SourceFinder for ReplayFinder {
    fn find(&self, seq: &[u8]) -> Option<SourceCall> {
        match self.0.iter().find(|(q, _)| q.as_slice() == seq) {
            Some((_, a)) => a.clone(),
            None => panic!("novel_finder.find({:?}) was not called by python", String::from_utf8_lossy(seq)),
        }
    }
}

/// [[seq, answer|null], ...] -> lookup fn (panics on an unrecorded query).
pub fn replay_strs(v: &Value) -> impl Fn(&[u8]) -> Option<String> + Sync {
    let answers: Vec<(Vec<u8>, Option<String>)> =
        v.as_array().map(|a| a.iter().map(|x| (b(&x[0]), x[1].as_str().map(String::from))).collect()).unwrap_or_default();
    move |q: &[u8]| match answers.iter().find(|(s, _)| s.as_slice() == q) {
        Some((_, a)) => a.clone(),
        None => panic!("callback({:?}) was not called by python", String::from_utf8_lossy(q)),
    }
}

/// hallmarks.SiteInfo (`vars(si)`: contig, L, R, located, tsd_*, en_*, slippage*)
pub fn site_info(v: &Value) -> SiteInfo {
    SiteInfo {
        contig: v["contig"].as_str().map(String::from),
        l: oi(&v["L"]),
        r: oi(&v["R"]),
        located: bo(&v["located"]),
        tsd_seq: s(&v["tsd_seq"]),
        tsd_len: oi(&v["tsd_len"]),
        tsd_verified: bo(&v["tsd_verified"]),
        en_motif: s(&v["en_motif"]),
        en_mismatches: oi(&v["en_mismatches"]),
        slippage: bo(&v["slippage"]),
        slippage_detail: s(&v["slippage_detail"]),
    }
}

/// `["A", 12]` / `["", 0]` -> (b'A', 12) / (0, 0)
pub fn run(v: &Value) -> (u8, i64) {
    (v[0].as_str().and_then(|x| x.bytes().next()).unwrap_or(0), i(&v[1]))
}

/// hallmarks.PolyAInfo
pub fn polya(v: &Value) -> PolyAInfo {
    PolyAInfo {
        strand: i(&v["strand"]) as i32,
        source: s(&v["source"]),
        left_run: run(&v["left_run"]),
        right_run: run(&v["right_run"]),
        both_sided: bo(&v["both_sided"]),
        length: f(&v["length"]),
    }
}

/// The AssemblyResult a `classify` event refers to (`in.assembly` inline, or `in.assembly_ref` =
/// the `seq` of an `assemble` event among `events`).
pub fn classify_assembly(ev: &Value, events: &[Value]) -> AssemblyResult {
    if let Some(r) = ev["in"]["assembly_ref"].as_i64() {
        let a = events
            .iter()
            .find(|e| e["kind"] == "assemble" && e["seq"].as_i64() == Some(r) && e["case"] == ev["case"])
            .expect("assemble event of assembly_ref");
        assembly(&a["out"])
    } else {
        assembly(&ev["in"]["assembly"])
    }
}

/// Cigars are not recorded (never read downstream): drop them before comparing.
pub fn strip_cigars(a: &mut AssemblyResult) {
    for l in a.layouts.iter_mut().chain(a.raw_layouts.iter_mut()) {
        for sg in l.segments.iter_mut() {
            sg.cigar = None;
        }
    }
}

/// `{k: v}` -> pairs (score weights / thresholds overrides; null -> empty).
pub fn pairs(v: &Value) -> Vec<(String, f64)> {
    match v {
        Value::Object(o) => o.iter().map(|(k, x)| (k.clone(), f(x))).collect(),
        _ => Vec::new(),
    }
}

/// Positional-or-keyword argument (index `idx`, name `name`) of a hallmark event.
pub fn arg<'a>(ev: &'a Value, idx: usize, name: &str) -> Option<&'a Value> {
    ev["args"].as_array().and_then(|a| a.get(idx)).or_else(|| ev["kwargs"].get(name))
}

/// The e2e manifest written by `make_golden.py e2e` (path: PEARTREE_RTE_E2E, default the
/// foundation session's scratchpad). None when absent: e2e tests skip.
pub fn e2e_manifest() -> Option<Value> {
    let p = std::env::var("PEARTREE_RTE_E2E").unwrap_or_else(|_| {
        "/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/352d40f9-4660-4fb5-a308-8d4d14cd55e7/scratchpad/golden/e2e.manifest.json".into()
    });
    let txt = std::fs::read_to_string(&p).ok()?;
    serde_json::from_str(&txt).ok()
}
