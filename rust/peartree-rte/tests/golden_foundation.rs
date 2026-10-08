//! Golden tests of the FOUNDATION modules (always on). Golden data: golden/data/*.jsonl.gz
//! (pytest suite, committed) and the e2e_phylo run (scratchpad; skipped when absent).

use peartree_rte::golden::{self, events_of, load_events};
use peartree_rte::inputs::EvidenceRead;
use peartree_rte::library::RteLibrary;
use peartree_rte::stream::{cap_reads, LocusLoader};
use rustc_hash::FxHashSet;
use std::path::Path;

fn names(rs: &[EvidenceRead]) -> Vec<String> {
    rs.iter().map(|r| format!("{}|{}|{}|{}|{}", r.side, r.role, r.sample, r.frag, r.r12)).collect()
}

#[test]
fn golden_present_and_parses() {
    let n_asm = events_of("assemble").len();
    assert!(n_asm > 50, "golden/data/pytest.jsonl.gz missing or empty ({n_asm} assemble events)");
    for e in events_of("assemble") {
        golden::ctx(&e["in"]["ctx"]);
        for r in e["in"]["reads"].as_array().unwrap() {
            golden::read(r);
        }
        golden::assembly(&e["out"]);
    }
    let all = load_events(&golden::golden_dir().join("pytest.jsonl.gz"));
    for e in events_of("classify") {
        golden::classify_assembly(&e, &all);
        golden::structure_call(&e["out"]);
    }
    for e in events_of("score") {
        golden::score_input(&e["in"]["score_input"]);
    }
    let g = golden::Genomes::from_events(&all);
    for e in all.iter().filter(|e| e["kind"].as_str().is_some_and(|k| k.starts_with("hallmark."))) {
        for a in e["args"].as_array().unwrap() {
            if a.get("genome").is_some() {
                assert!(g.get(&a["genome"]).is_ok());
            }
        }
    }
}

/// record.rs `row()` / `detail_string()` / `_f` against every python record in the golden data.
#[test]
fn record_rows_match_python() {
    let mut n = 0;
    for kind in ["annotate", "annotate_key", "final"] {
        for e in events_of(kind) {
            let rec = golden::record(&e["out"]);
            let want: Vec<String> = e["out"]["row"].as_array().unwrap().iter().map(golden::s).collect();
            assert_eq!(rec.row(), want, "{} {}", e["case"], e["out"]["insertion_id"]);
            n += 1;
        }
    }
    if let Some(m) = golden::e2e_manifest() {
        for e in load_events(Path::new(m["golden"].as_str().unwrap())) {
            if matches!(e["kind"].as_str(), Some("annotate" | "annotate_key" | "final")) {
                let rec = golden::record(&e["out"]);
                let want: Vec<String> = e["out"]["row"].as_array().unwrap().iter().map(golden::s).collect();
                assert_eq!(rec.row(), want, "e2e {}", e["out"]["insertion_id"]);
                n += 1;
            }
        }
    }
    assert!(n > 100, "{n} records checked");
}

/// transduction::known_source (foundation) against python.
#[test]
fn known_source_matches_python() {
    let mut libs: Vec<(String, RteLibrary)> = Vec::new();
    for e in events_of("known_source") {
        let dir = golden::lib_dir(&e).unwrap();
        let key = dir.to_string_lossy().to_string();
        if !libs.iter().any(|(k, _)| *k == key) {
            libs.push((key.clone(), RteLibrary::open(&key, None).unwrap()));
        }
        let lib = &libs.iter().find(|(k, _)| *k == key).unwrap().1;
        let segs: Vec<_> = e["in"]["segments"].as_array().unwrap().iter().map(golden::segment).collect();
        let refs: Vec<&_> = segs.iter().collect();
        assert_eq!(peartree_rte::transduction::known_source(&refs, lib), golden::source_call(&e["out"]), "{}", e["case"]);
    }
}

/// Streaming store + cap against python on the e2e_phylo sidecars: for every python `annotate`
/// call, the pooled read list (combine reads, then GT reads when python pooled them) capped with
/// rte_max_reads must give exactly python's capped read names, in order; junction rows equal.
#[test]
fn e2e_stream_and_cap_match_python() {
    let Some(m) = golden::e2e_manifest() else {
        eprintln!("skipping: no e2e manifest");
        return;
    };
    let p = |k: &str| m[k].as_str().map(std::path::PathBuf::from);
    let inputs = peartree_rte::io::read_inputs(&p("inputs").unwrap()).unwrap();
    assert_eq!(inputs.len() as i64, m["n_inputs"].as_i64().unwrap());
    let wanted: FxHashSet<String> = inputs.iter().map(|i| i.title.clone()).collect();
    let tmp = std::env::temp_dir();
    let loader = LocusLoader::open(p("evidence").as_deref(), p("reads").as_deref(), p("gt_reads").as_deref(), &wanted, &tmp).unwrap();
    let cfg = peartree_rte::config::RteConfig::load(&p("config").unwrap()).unwrap();
    let mut checked = 0;
    let mut pooled = 0;
    for e in load_events(&p("golden").unwrap()).into_iter().filter(|e| e["kind"] == "annotate") {
        let key = e["in"]["input"]["locus"].as_str().unwrap();
        let ev = &e["in"]["evidence"];
        if ev.is_null() {
            continue;
        }
        let d = loader.load(key).unwrap();
        let want_j: Vec<_> = ev["junctions"].as_array().unwrap().iter().map(golden::junction).collect();
        assert_eq!(d.evidence.junctions, want_j, "{key} junctions");
        let n = ev["n_reads"].as_u64().unwrap() as usize;
        let mut reads = d.evidence.reads.clone();
        if n != reads.len() {
            reads.extend(d.gt_reads.iter().cloned());
            pooled += 1;
        }
        assert_eq!(reads.len(), n, "{key}: read count (combine {} + gt {})", d.evidence.reads.len(), d.gt_reads.len());
        if let Some(full) = ev["reads"].as_array() {
            let want: Vec<EvidenceRead> = full.iter().map(golden::read).collect();
            assert_eq!(names(&reads), names(&want), "{key}: pooled read order");
            for (a, b) in reads.iter().zip(&want) {
                assert_eq!(a.seq(), b.seq(), "{key}: read sequence");
            }
        }
        let capped = cap_reads(reads, cfg.max_reads);
        let want: Vec<String> = ev["capped"].as_array().unwrap().iter().map(golden::s).collect();
        assert_eq!(names(&capped), want, "{key}: capped reads");
        checked += 1;
    }
    eprintln!("e2e: {checked} annotate calls checked ({pooled} pooled with genotype reads)");
    assert!(checked > 200 && pooled > 50);
}
