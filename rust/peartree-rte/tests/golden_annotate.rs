//! WP-INT golden: `RteAnnotator::annotate` per insertion against python's recorded `annotate`
//! events (pytest golden + the e2e golden when its manifest is present). Each event carries the
//! InsertionInput, the junction evidence and every read; the novel-source and pre-mRNA answers
//! are replayed from the event's `classify` child (python's callbacks may have been test
//! stand-ins), the exon track from its `exon_find` / `exon_structure` children.

use peartree_rte::annotator::RteAnnotator;
use peartree_rte::config::RteConfig;
use peartree_rte::golden::{self, load_events, Genomes, ReplayFinder};
use peartree_rte::inputs::InsertionEvidence;
use peartree_rte::library::RteLibrary;
use peartree_rte::pseudogene::ExonsByGene;
use rustc_hash::FxHashMap;
use serde_json::Value;
use std::path::{Path, PathBuf};

fn all_events() -> Vec<Value> {
    let mut v = load_events(&golden::golden_dir().join("pytest.jsonl.gz"));
    if let Some(m) = golden::e2e_manifest() {
        let mut e2e = load_events(Path::new(m["golden"].as_str().unwrap()));
        golden::namespace_genomes(&mut e2e, "e2e:");
        // the e2e seq numbers restart at 1: offset them so parent links stay unique
        for e in e2e.iter_mut() {
            for k in ["seq", "parent", "begin"] {
                if let Some(x) = e[k].as_i64() {
                    e[k] = Value::from(x + 10_000_000);
                }
            }
        }
        v.extend(e2e);
    }
    v
}

#[test]
fn golden_annotate() {
    let all = all_events();
    let g = Genomes::from_events(&all);
    let mut children: FxHashMap<i64, Vec<&Value>> = FxHashMap::default();
    for e in &all {
        if let Some(p) = e["parent"].as_i64() {
            children.entry(p).or_default().push(e);
        }
    }
    let mut libs: Vec<(PathBuf, RteLibrary)> = Vec::new();
    let (mut n, mut skipped) = (0, 0);
    for e in all.iter().filter(|e| e["kind"] == "annotate") {
        let (Ok(genome), Ok(remap)) = (g.get(&e["genome"]), g.get(&e["remap"])) else {
            skipped += 1;
            continue;
        };
        let Some(dir) = golden::lib_dir(e) else {
            skipped += 1;
            continue;
        };
        if !libs.iter().any(|(k, _)| *k == dir) {
            libs.push((dir.clone(), RteLibrary::open(&dir.to_string_lossy(), None).unwrap()));
        }
        let lib = &libs.iter().find(|(k, _)| *k == dir).unwrap().1;
        let cfg = RteConfig::from_json(&e["cfg"]).unwrap();
        let kids = children.get(&e["begin"].as_i64().unwrap()).cloned().unwrap_or_default();
        let classify = kids.iter().find(|c| c["kind"] == "classify").copied();
        let finder = classify.map(|c| ReplayFinder::from_json(&c["in"]["novel_answers"])).unwrap_or(ReplayFinder(Vec::new()));
        let premrna = golden::replay_strs(&classify.map(|c| c["in"]["premrna_answers"].clone()).unwrap_or(Value::Null));
        let has_premrna = classify.is_some_and(|c| golden::bo(&c["in"]["has_premrna"]));
        let (mut exons, mut strands) = (ExonsByGene::default(), FxHashMap::default());
        if let Some(x) = kids.iter().find(|c| c["kind"] == "exon_find" || c["kind"] == "exon_structure") {
            for (k, v) in x["in"]["exons"].as_object().unwrap() {
                exons.insert(k.clone(), v.as_array().unwrap().iter().map(|x| (golden::s(&x[0]), golden::i(&x[1]), golden::i(&x[2]))).collect());
            }
            for (k, v) in x["in"]["strands"].as_object().unwrap() {
                strands.insert(k.clone(), golden::s(v).chars().next().unwrap());
            }
        }
        let mut ann = RteAnnotator::with_parts(&cfg, lib, genome, remap, None, None, exons, strands, None);
        ann.source_finder_override = Some(&finder);
        if has_premrna {
            ann.premrna_override = Some(&premrna);
        }
        let inp = peartree_rte::io::input_from_json(&e["in"]["input"]).unwrap();
        let mut ev = InsertionEvidence::new(&inp.title);
        if let Some(x) = e["in"]["evidence"].as_object() {
            for j in x["junctions"].as_array().unwrap() {
                ev.set_junction(golden::junction(j));
            }
            ev.reads = x["reads"].as_array().unwrap().iter().map(golden::read).collect();
        }
        let got = ann.annotate(&inp, &ev);
        let want = golden::record(&e["out"]);
        assert_eq!(got, want, "{} seq {}", e["case"], e["seq"]);
        assert_eq!(got.row(), want.row());
        n += 1;
    }
    eprintln!("golden_annotate: {n} annotate events matched, {skipped} skipped (genome not local / mock library)");
    assert!(n > 100);
}
