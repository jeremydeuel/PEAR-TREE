//! WP-TD parity tests that have no pytest/e2e golden: the gene model, the exon-track loaders and
//! `L1Rmsk` (all three file layouts). Expectations in `golden/data/td/expect.json` are produced
//! by `golden/make_td_expect.py`, which runs the PYTHON implementations (annotate_v2.GeneModel,
//! pseudogene.load_*, transduction.L1Rmsk) on the tracks in the same directory.

use peartree_rte::config::GeneModelCfg;
use peartree_rte::genemodel::GeneModel;
use peartree_rte::pseudogene::{load_exons_by_gene, load_gene_strands};
use peartree_rte::transduction::L1Rmsk;
use serde_json::Value;
use std::path::PathBuf;

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden/data/td")
}

fn expect() -> Value {
    serde_json::from_str(&std::fs::read_to_string(dir().join("expect.json")).unwrap()).unwrap()
}

fn i(v: &Value) -> i64 {
    v.as_i64().unwrap()
}

fn cfg1() -> GeneModelCfg {
    GeneModelCfg { splice_donor_window: 6, splice_acceptor_window: 3, splice_ppt_window: 17, splice_branch_window: 45, promoter_up: 2000 }
}

fn cfg2(e: &Value) -> GeneModelCfg {
    let c = &e["cfg2"];
    GeneModelCfg {
        splice_donor_window: i(&c["splice_donor_window"]),
        splice_acceptor_window: i(&c["splice_acceptor_window"]),
        splice_ppt_window: i(&c["splice_ppt_window"]),
        splice_branch_window: i(&c["splice_branch_window"]),
        promoter_up: i(&c["promoter_up"]),
    }
}

fn feat(t: (i64, &str, &str)) -> Value {
    serde_json::json!([t.0, t.1, t.2])
}

#[test]
fn gene_model_matches_python() {
    let e = expect();
    let gm = GeneModel::open(&dir().join("track.tsv"), &cfg1()).unwrap();
    let gm2 = GeneModel::open(&dir().join("track.tsv"), &cfg2(&e)).unwrap();

    // loaded genes: order, merged exons, strands
    for (c, want) in e["genes"].as_object().unwrap() {
        let got = &gm.genes[c];
        assert_eq!(got.len(), want.as_array().unwrap().len(), "{c}");
        for (g, w) in got.iter().zip(want.as_array().unwrap()) {
            let exons: Vec<(i64, i64)> = w[4].as_array().unwrap().iter().map(|x| (i(&x[0]), i(&x[1]))).collect();
            assert_eq!((g.start, g.end, g.name.as_str(), g.strand.as_str(), g.exons.clone()), (i(&w[0]), i(&w[1]), w[2].as_str().unwrap(), w[3].as_str().unwrap(), exons), "{c}");
        }
        assert_eq!(gm.maxspan[c], i(&e["maxspan"][c]));
    }
    assert_eq!(gm.genes.len(), e["genes"].as_object().unwrap().len());

    for r in e["resolve"].as_array().unwrap() {
        assert_eq!(gm.resolve(r[0].as_str().unwrap()).as_deref(), r[1].as_str());
    }

    // _candidates: membership AND order
    let mut n = 0;
    for r in e["candidates"].as_array().unwrap() {
        let got = gm.candidates(r[0].as_str().unwrap(), i(&r[1]));
        match (&got, &r[2]) {
            (None, Value::Null) => {}
            (Some(g), Value::Array(w)) => {
                let g: Vec<(i64, i64, &str, &str)> = g.iter().map(|g| (g.start, g.end, g.name.as_str(), g.strand.as_str())).collect();
                let w: Vec<(i64, i64, &str, &str)> = w.iter().map(|x| (i(&x[0]), i(&x[1]), x[2].as_str().unwrap(), x[3].as_str().unwrap())).collect();
                assert_eq!(g, w, "{r}");
            }
            _ => panic!("candidates mismatch {r}: {got:?}"),
        }
        n += 1;
    }
    assert!(n > 2000);

    // _genic_feature with both window sets, over the candidate grid and a dense intron walk
    let find = |m: &GeneModel, name: &str, strand: &str, p: i64| {
        let g = m.genes.values().flatten().find(|g| g.name == name && g.strand == strand).unwrap();
        feat(m.genic_feature(&g.exons, strand, p))
    };
    for (key, m) in [("features", &gm), ("features2", &gm2)] {
        for r in e[key].as_array().unwrap() {
            assert_eq!(find(m, r[0].as_str().unwrap(), r[1].as_str().unwrap(), i(&r[2])), r[3], "{key} {r}");
        }
    }
    for r in e["walk"].as_array().unwrap() {
        let (name, strand, p) = (r[0].as_str().unwrap(), r[1].as_str().unwrap(), i(&r[2]));
        assert_eq!(find(&gm, name, strand, p), r[3], "walk {r}");
        assert_eq!(find(&gm2, name, strand, p), r[4], "walk2 {r}");
    }
    for r in e["splice"].as_array().unwrap() {
        assert_eq!(feat(gm.splice_class(i(&r[0]), i(&r[1]))), r[2], "{r}");
        assert_eq!(feat(gm2.splice_class(i(&r[0]), i(&r[1]))), r[3], "{r}");
    }
}

#[test]
fn exon_track_loaders_match_python() {
    let e = expect();
    let got = load_exons_by_gene(&dir().join("track.tsv")).unwrap();
    let want = e["exons_by_gene"].as_object().unwrap();
    assert_eq!(got.len(), want.len());
    for (g, w) in want {
        let w: Vec<(String, i64, i64)> = w.as_array().unwrap().iter().map(|x| (x[0].as_str().unwrap().to_string(), i(&x[1]), i(&x[2]))).collect();
        assert_eq!(got[g], w, "{g}");
    }
    let strands = load_gene_strands(&dir().join("track.tsv")).unwrap();
    let want = e["gene_strands"].as_object().unwrap();
    assert_eq!(strands.len(), want.len());
    for (g, w) in want {
        assert_eq!(strands[g], w.as_str().unwrap().chars().next().unwrap(), "{g}");
    }
}

#[test]
fn l1_rmsk_matches_python() {
    let e = expect();
    for (name, file) in [("out", "rmsk.out"), ("gz", "rmsk.out.gz"), ("ucsc", "rmsk_ucsc.txt")] {
        let r = L1Rmsk::open(&dir().join(file), 5500).unwrap();
        let w = &e["rmsk"][name];
        let rows = w["rows"].as_object().unwrap();
        assert_eq!(r.by_contig.len(), rows.len(), "{name}");
        for (c, v) in rows {
            let want: Vec<_> = v.as_array().unwrap().iter().map(|x| (i(&x[0]), i(&x[1]), x[2].as_str().unwrap().chars().next().unwrap(), x[3].as_str().unwrap().to_string(), x[4].as_f64().unwrap())).collect();
            assert_eq!(r.by_contig[c], want, "{name} {c}");
        }
        let qs = w["queries"].as_array().unwrap();
        for (key, dist) in [("up", 15000), ("up_small", 3000)] {
            for (q, want) in qs.iter().zip(w[key].as_array().unwrap()) {
                let got = r.upstream_of(q[0].as_str().unwrap(), i(&q[1]), i(&q[2]), q[3].as_str().unwrap().chars().next().unwrap(), dist);
                let want: Vec<_> = want.as_array().unwrap().iter().map(|x| (i(&x[0]), i(&x[1]), x[2].as_str().unwrap().chars().next().unwrap(), x[3].as_str().unwrap().to_string(), x[4].as_f64().unwrap())).collect();
                assert_eq!(got, want, "{name} {key} {q}");
            }
        }
        // the dataset must actually exercise hits
        let nonempty = w["up"].as_array().unwrap().iter().filter(|x| !x.as_array().unwrap().is_empty()).count();
        assert!(nonempty > 10, "{name}: only {nonempty} non-empty");
    }
}
