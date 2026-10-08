//! Golden tests of the WORK PACKAGES (PORT_PLAN.md). Each is `#[ignore]` until its package lands;
//! the package owner removes the `#[ignore]` line of ITS tests only (this file is shared -- edit
//! nothing else here; add extra tests in your own module's `#[cfg(test)]`).
//!
//! Run one: `cargo test --release --test golden_modules -- --ignored golden_hallmarks`

use peartree_rte::config::{AssemblyCfg, PseudogeneCfg, StructureCfg, TransductionCfg};
use peartree_rte::golden::{self, load_events, Genomes};
use peartree_rte::library::RteLibrary;
use serde_json::Value;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

fn all_events() -> Vec<Value> {
    let mut v = load_events(&golden::golden_dir().join("pytest.jsonl.gz"));
    if let Some(m) = golden::e2e_manifest() {
        if std::env::var_os("PEARTREE_RTE_SKIP_E2E").is_none() {
            v.extend(load_events(Path::new(m["golden"].as_str().unwrap())));
        }
    }
    v
}

/// One RteLibrary per directory (fixture / resources / snapshot), opened once.
#[derive(Default)]
struct Libs(Vec<(PathBuf, RteLibrary)>);

impl Libs {
    fn get(&mut self, ev: &Value) -> Option<&RteLibrary> {
        let dir = match golden::lib_dir(ev) {
            Some(d) => d,
            None if ev["lib"].get("mock").is_some() => golden::repo_path("test/fixtures/rte_library"),
            None => return None,
        };
        if !self.0.iter().any(|(k, _)| *k == dir) {
            let mut lib = RteLibrary::open(&dir.to_string_lossy(), None).unwrap();
            if let Some(m) = ev["lib"].get("mock") {
                // the test stand-in: fixture consensus + its own polymorphic_l1 list
                lib.polymorphic_l1 = m["polymorphic_l1"]
                    .as_array()
                    .map(|a| a.iter().map(|x| (golden::s(&x[0]), golden::i(&x[1]), golden::s(&x[2]))).collect())
                    .unwrap_or_default();
            }
            self.0.push((dir.clone(), lib));
        }
        self.0.iter().find(|(k, _)| *k == dir).map(|(_, l)| l)
    }
}

fn opt_obj(v: &Value) -> Option<&serde_json::Map<String, Value>> {
    v.as_object()
}

// ------------------------------------------------------------------------------- WP-HALL
#[test]
fn golden_hallmarks() {
    use peartree_rte::hallmarks as H;
    let all = all_events();
    let g = Genomes::from_events(&all);
    let b = |v: Option<&Value>| v.map(golden::b).unwrap_or_default();
    let genome = |v: Option<&Value>| v.map(|x| g.get(&x["genome"])).unwrap_or(Ok(None));
    let mut n = 0;
    for e in &all {
        let Some(kind) = e["kind"].as_str().and_then(|k| k.strip_prefix("hallmark.")) else { continue };
        let ctx = format!("{} {} seq {}", kind, e["case"], e["seq"]);
        let out = &e["out"];
        match kind {
            "split_junction" => {
                let (a, bb, c, d) = H::split_junction(&b(golden::arg(e, 0, "left_seq")), &b(golden::arg(e, 1, "right_seq")));
                let want: Vec<Vec<u8>> = out.as_array().unwrap().iter().map(golden::b).collect();
                assert_eq!(vec![a, bb, c, d], want, "{ctx}");
            }
            "parse_locus" => {
                let got = H::parse_locus(golden::arg(e, 0, "title").and_then(|v| v.as_str()).unwrap_or(""));
                let want = if out.is_null() { None } else { Some((golden::s(&out[0]), golden::i(&out[1]), golden::i(&out[2]))) };
                assert_eq!(got, want, "{ctx}");
            }
            "edge_run" => {
                let at_end = golden::arg(e, 1, "at_end").map(golden::bo).unwrap();
                assert_eq!(H::edge_run(&b(golden::arg(e, 0, "insert")), at_end), golden::run(out), "{ctx}");
            }
            "polya_info" => {
                let j = |i: usize, k: &str| golden::arg(e, i, k).filter(|v| !v.is_null()).map(|v| golden::junction(&v["junction"]));
                let (jl, jr) = (j(2, "ev_left"), j(3, "ev_right"));
                let min_len = golden::arg(e, 4, "min_len").map(golden::f).unwrap_or(10.0);
                let got = H::polya_info(&b(golden::arg(e, 0, "left_insert")), &b(golden::arg(e, 1, "right_insert")), jl.as_ref(), jr.as_ref(), min_len);
                assert_eq!(got, golden::polya(out), "{ctx}");
            }
            "locate_site" => {
                if e["kwargs"].as_object().is_some_and(|k| !k.is_empty()) || e["args"].as_array().unwrap().len() > 4 {
                    continue; // non-default probe_len / slack: not in the Rust signature
                }
                let Ok(gn) = genome(golden::arg(e, 3, "genome")) else { continue };
                let got = H::locate_site(golden::arg(e, 0, "title").and_then(|v| v.as_str()).unwrap(), &b(golden::arg(e, 1, "left_flank")), &b(golden::arg(e, 2, "right_flank")), gn);
                assert_eq!(got, golden::site_info(out), "{ctx}");
            }
            "tsd_from_flanks" => {
                if e["kwargs"].as_object().is_some_and(|k| !k.is_empty()) || e["args"].as_array().unwrap().len() > 2 {
                    continue;
                }
                let got = H::tsd_from_flanks(&b(golden::arg(e, 0, "left_flank")), &b(golden::arg(e, 1, "right_flank")));
                assert_eq!(got, (golden::s(&out[0]), golden::oi(&out[1]), golden::bo(&out[2])), "{ctx}");
            }
            "target_site" | "en_motif" | "slippage_context" => {
                let mut si = golden::site_info(&golden::arg(e, 0, "si").unwrap()["site"]);
                if kind == "target_site" {
                    let Ok(gn) = genome(golden::arg(e, 3, "genome")) else { continue };
                    H::target_site(&mut si, &b(golden::arg(e, 1, "left_flank")), &b(golden::arg(e, 2, "right_flank")), gn);
                } else {
                    if kind == "slippage_context" && e["args"].as_array().unwrap().len() > 3 {
                        continue; // non-default min_run
                    }
                    let Ok(gn) = genome(golden::arg(e, 2, "genome")) else { continue };
                    let strand = golden::arg(e, 1, "strand").map(golden::i).unwrap() as i32;
                    if kind == "en_motif" {
                        H::en_motif(&mut si, strand, gn);
                    } else {
                        H::slippage_context(&mut si, strand, gn);
                    }
                }
                assert_eq!(si, golden::site_info(out), "{ctx}");
            }
            "en_bin" => assert_eq!(H::en_bin(golden::arg(e, 0, "mm").and_then(golden::oi)), golden::s(out), "{ctx}"),
            "foldback" => {
                if e["args"].as_array().unwrap().len() > 4 || e["kwargs"].as_object().is_some_and(|k| !k.is_empty()) {
                    continue;
                }
                let got = H::foldback(&b(golden::arg(e, 0, "left_insert")), &b(golden::arg(e, 1, "left_flank")), &b(golden::arg(e, 2, "right_insert")), &b(golden::arg(e, 3, "right_flank")));
                assert_eq!(got, golden::i(out), "{ctx}");
            }
            other => panic!("unknown hallmark event {other}"),
        }
        n += 1;
    }
    assert!(n > 1000, "{n} hallmark events");
}

#[test]
fn golden_score() {
    let mut n = 0;
    for e in all_events().iter().filter(|e| e["kind"] == "score") {
        let si = golden::score_input(&e["in"]["score_input"]);
        let (sc, pts, call) = peartree_rte::score::score(&si, &golden::pairs(&e["in"]["weights"]), &golden::pairs(&e["in"]["thresholds"]));
        assert_eq!((sc, pts.as_str(), call.as_str()), (golden::f(&e["out"][0]), e["out"][1].as_str().unwrap(), e["out"][2].as_str().unwrap()), "{}", e["case"]);
        n += 1;
    }
    assert!(n > 100);
}

// -------------------------------------------------------------------------------- WP-ASM
#[test]
#[ignore = "WP-ASM"]
fn golden_assemble() {
    let mut libs = Libs::default();
    let mut n = 0;
    for e in all_events().iter().filter(|e| e["kind"] == "assemble") {
        let lib = libs.get(e).unwrap();
        let cfg = AssemblyCfg::from(opt_obj(&e["cfg"])).unwrap();
        let asm = peartree_rte::assembly::Assembler::new(lib, cfg);
        let ctx = golden::ctx(&e["in"]["ctx"]);
        let js: Vec<(String, peartree_rte::assembly::JunctionSeq)> = e["in"]["junction_seqs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| (golden::s(&x[0]), (golden::b(&x[1]), (golden::i(&x[2][0]), golden::i(&x[2][1])))))
            .collect();
        let reads: Vec<_> = e["in"]["reads"].as_array().unwrap().iter().map(golden::read).collect();
        let hint = e["in"]["strand_hint"].as_array().map(|h| (golden::i(&h[0]) as i32, golden::s(&h[1])));
        let mut got = asm.assemble(&ctx, &js, &reads, hint);
        golden::strip_cigars(&mut got);
        let want = golden::assembly(&e["out"]);
        assert_eq!(got, want, "{} seq {}", e["case"], e["seq"]);
        n += 1;
    }
    assert!(n > 100);
}

// ------------------------------------------------------------------------------ WP-STRUCT
#[test]
#[ignore = "WP-STRUCT"]
fn golden_classify() {
    use peartree_rte::structure::{classify, PseudogeneArg, SourceFinder};
    let all = all_events();
    let mut libs = Libs::default();
    let mut n = 0;
    for e in all.iter().filter(|e| e["kind"] == "classify") {
        let lib = libs.get(e).unwrap();
        let res = golden::classify_assembly(e, &all);
        let ctx = (!e["in"]["ctx"].is_null()).then(|| golden::ctx(&e["in"]["ctx"]));
        let cfg = StructureCfg::from(opt_obj(&e["cfg"])).unwrap();
        let finder = golden::ReplayFinder::from_json(&e["in"]["novel_answers"]);
        let premrna = golden::replay_strs(&e["in"]["premrna_answers"]);
        let pg = &e["in"]["pseudogene"];
        let genes = golden::strs(&pg["genes"]);
        let hits: Vec<(String, String)> = pg["hits"].as_array().map(|a| a.iter().map(|h| (golden::s(&h[0]), golden::s(&h[1]))).collect()).unwrap_or_default();
        // pseudogene structure fn: python's answers in call order
        let answers: Vec<Option<String>> = e["in"]["pg_structure_answers"].as_array().unwrap().iter().map(|a| a[1].as_str().map(String::from)).collect();
        let k = Mutex::new(0usize);
        let pgs = |_l: &[peartree_rte::assembly::ReadLayout]| {
            let mut k = k.lock().unwrap();
            *k += 1;
            answers[*k - 1].clone()
        };
        let pg_arg = PseudogeneArg { genes: &genes, hits: &hits, structure_fn: pg["has_structure_fn"].as_bool().unwrap_or(false).then_some(&pgs as _) };
        let got = classify(
            &res,
            lib,
            ctx.as_ref(),
            &cfg,
            golden::bo(&e["in"]["has_novel"]).then_some(&finder as &dyn SourceFinder),
            golden::bo(&e["in"]["has_premrna"]).then_some(&premrna as _),
            (!pg.is_null()).then_some(&pg_arg),
            e["in"]["legacy_class"].as_str(),
        );
        assert_eq!(got, golden::structure_call(&e["out"]), "{} seq {}", e["case"], e["seq"]);
        n += 1;
    }
    assert!(n > 100);
}

// --------------------------------------------------------------------------------- WP-TD
struct ReplayLocator(Vec<(Vec<u8>, Vec<peartree_rte::transduction::LocatorHit>)>);

impl peartree_rte::transduction::Locator for ReplayLocator {
    fn locate(&self, seq: &[u8]) -> Vec<peartree_rte::transduction::LocatorHit> {
        self.0.iter().find(|(q, _)| q.as_slice() == seq).map(|(_, h)| h.clone()).expect("locator query python did not make")
    }
}

#[test]
#[ignore = "WP-TD"]
fn golden_novel_find() {
    use peartree_rte::structure::SourceFinder;
    use peartree_rte::transduction::{L1Rmsk, LocatorHit, NovelSourceFinder};
    let all = all_events();
    let g = Genomes::from_events(&all);
    let mut libs = Libs::default();
    let mut n = 0;
    for e in all.iter().filter(|e| e["kind"] == "novel_find") {
        let Ok(remap) = g.get(&e["remap"]) else { continue };
        let cfg = TransductionCfg::from(opt_obj(&e["cfg"])).unwrap();
        let lib = libs.get(e).unwrap();
        let answers = e["in"]["locator_answers"].as_array().unwrap().iter().map(|a| {
            let hits = a[1].as_array().unwrap().iter().map(|h| LocatorHit {
                contig: golden::s(&h[0]),
                start: golden::i(&h[1]),
                end: golden::i(&h[2]),
                strand: golden::s(&h[3]).chars().next().unwrap(),
                mapq: golden::i(&h[4]),
                identity: golden::f(&h[5]),
            });
            (golden::b(&a[0]), hits.collect())
        });
        let has_locator = !e["in"]["locator_answers"].as_array().unwrap().is_empty() || golden::bo(&e["in"]["available"]);
        let rmsk = e["in"]["rmsk"].as_str().map(|p| L1Rmsk::open(&golden::repo_path(p), cfg.novel_source_min_len).unwrap());
        let f = NovelSourceFinder {
            lib,
            cfg: cfg.clone(),
            rmsk,
            locator: has_locator.then(|| Box::new(ReplayLocator(answers.collect())) as Box<dyn peartree_rte::transduction::Locator>),
            genome: remap,
            cohort_l1: e["in"]["cohort_l1"].as_array().unwrap().iter().map(|c| (golden::s(&c[0]), golden::i(&c[1]), golden::s(&c[2]).chars().next().unwrap())).collect(),
            ident_cache: Default::default(),
        };
        assert_eq!(f.find(&golden::b(&e["in"]["seq"])), golden::source_call(&e["out"]), "{} seq {}", e["case"], e["seq"]);
        n += 1;
    }
    assert!(n > 5);
}

#[test]
#[ignore = "WP-TD"]
fn golden_cons_identity() {
    for e in all_events().iter().filter(|e| e["kind"] == "cons_identity") {
        let got = peartree_rte::transduction::cons_identity(&golden::b(&e["in"]["seq"]), &golden::b(&e["in"]["cons"]));
        assert_eq!(got, golden::f(&e["out"]), "{}", e["case"]);
    }
}

#[test]
#[ignore = "WP-TD"]
fn golden_exon_junctions() {
    use peartree_rte::pseudogene::ExonJunctionIndex;
    let all = all_events();
    let g = Genomes::from_events(&all);
    let mut n = 0;
    for e in all.iter().filter(|e| e["kind"] == "exon_find" || e["kind"] == "exon_structure") {
        let Ok(remap) = g.get(&e["remap"]) else { continue };
        let exons = e["in"]["exons"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(k, v)| (k.clone(), v.as_array().unwrap().iter().map(|x| (golden::s(&x[0]), golden::i(&x[1]), golden::i(&x[2]))).collect()))
            .collect();
        let strands = e["in"]["strands"].as_object().unwrap().iter().map(|(k, v)| (k.clone(), golden::s(v).chars().next().unwrap())).collect();
        let idx = ExonJunctionIndex { cfg: PseudogeneCfg::from(opt_obj(&e["cfg"])).unwrap(), exons, genome: remap, strands };
        let genes = golden::strs(&e["in"]["genes"]);
        if e["kind"] == "exon_find" {
            let seqs: Vec<(String, Vec<u8>)> = e["in"]["seqs"].as_array().unwrap().iter().map(|x| (golden::s(&x[0]), golden::b(&x[1]))).collect();
            let want: Vec<(String, String)> = e["out"].as_array().unwrap().iter().map(|h| (golden::s(&h[0]), golden::s(&h[1]))).collect();
            assert_eq!(idx.find(&genes, &seqs), want, "{}", e["case"]);
        } else {
            if golden::i(&e["in"]["probe"]) != 25 {
                continue;
            }
            let seqs: Vec<Vec<u8>> = e["in"]["seqs"].as_array().unwrap().iter().map(golden::b).collect();
            assert_eq!(idx.structure(&genes, &seqs, golden::i(&e["in"]["tol"])), e["out"].as_str().map(String::from), "{}", e["case"]);
        }
        n += 1;
    }
    assert!(n > 5);
}

// -------------------------------------------------------------------------------- WP-INT
/// End to end: the binary over the e2e_phylo sidecars == python's final rows (+ gt columns).
#[test]
#[ignore = "WP-INT"]
fn golden_e2e_binary() {
    let Some(m) = golden::e2e_manifest() else {
        eprintln!("skipping: no e2e manifest");
        return;
    };
    let out = std::env::temp_dir().join(format!("peartree-rte-e2e-{}.tsv", std::process::id()));
    let mut cmd = std::process::Command::new(env!("CARGO_BIN_EXE_peartree-rte"));
    cmd.args(["annotate", "--config", m["config"].as_str().unwrap(), "--inputs", m["inputs"].as_str().unwrap(), "--out", out.to_str().unwrap()]);
    for (flag, key) in [("--evidence", "evidence"), ("--reads", "reads"), ("--gt-reads", "gt_reads")] {
        if let Some(p) = m[key].as_str() {
            cmd.args([flag, p]);
        }
    }
    let st = cmd.status().unwrap();
    assert!(st.success());
    let txt = std::fs::read_to_string(&out).unwrap();
    let mut rows = txt.lines().skip(1).map(|l| l.split('\t').map(String::from).collect::<Vec<_>>());
    let finals: Vec<Value> = load_events(Path::new(m["golden"].as_str().unwrap())).into_iter().filter(|e| e["kind"] == "final").collect();
    assert_eq!(finals.len() as i64, m["n_inputs"].as_i64().unwrap());
    for e in &finals {
        let got = rows.next().unwrap();
        let rec = &e["out"];
        let mut want: Vec<String> = vec![golden::s(&rec["insertion_id"])];
        want.extend(rec["row"].as_array().unwrap().iter().map(golden::s));
        want.push(golden::i(&rec["gt_reads"]).to_string());
        want.push(rec["gt_changed"].as_str().filter(|s| !s.is_empty()).unwrap_or(".").to_string());
        assert_eq!(&got[..want.len()], &want[..], "{}", want[0]);
    }
    std::fs::remove_file(&out).ok();
}
