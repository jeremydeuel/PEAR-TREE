//! peartree-rte command line (SPEC.md "Interface").
//!
//!   peartree-rte annotate --config CFG.json --inputs IN.jsonl[.gz] --out OUT.tsv[.gz]
//!                         [--evidence P.insertions.evidence.tsv.gz] [--reads P.insertions.reads.fa.gz]
//!                         [--gt-reads P.insertions.genotype_reads.fa.gz] [--tmp-dir DIR]
//!                         [--threads N] [--chunk N]
//!   peartree-rte scan --reads FASTA[.gz] [--gt] [--tmp-dir DIR]
//!       one streaming pass + every locus read back (memory / grouping diagnostics)

#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use peartree_rte::stream::{map_loci, LocusLoader, ReadStore, StoreOptions};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

const USAGE: &str = "usage:\n  peartree-rte annotate --config CFG.json --inputs IN.jsonl --out OUT.tsv \
[--evidence TSV] [--reads FA] [--gt-reads FA] [--tmp-dir DIR] [--threads N] [--chunk N]\n  \
peartree-rte scan --reads FA [--gt] [--tmp-dir DIR]";

fn parse_flags(args: &[String], flags: &[&str], switches: &[&str]) -> Result<(HashMap<String, String>, Vec<String>), String> {
    let mut kv = HashMap::new();
    let mut sw = Vec::new();
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if switches.contains(&a.as_str()) {
            sw.push(a.clone());
        } else if flags.contains(&a.as_str()) {
            let v = it.next().ok_or_else(|| format!("{a} needs a value"))?;
            kv.insert(a.clone(), v.clone());
        } else {
            return Err(format!("unknown argument {a:?}\n{USAGE}"));
        }
    }
    Ok((kv, sw))
}

fn tmp_dir(kv: &HashMap<String, String>, default_near: &Path) -> PathBuf {
    kv.get("--tmp-dir").map(PathBuf::from).unwrap_or_else(|| {
        default_near.parent().filter(|p| !p.as_os_str().is_empty()).map(Path::to_path_buf).unwrap_or_else(|| PathBuf::from("."))
    })
}

fn scan(args: &[String]) -> Result<(), String> {
    let (kv, sw) = parse_flags(args, &["--reads", "--tmp-dir"], &["--gt"])?;
    let reads = PathBuf::from(kv.get("--reads").ok_or("scan: --reads is required")?);
    let t0 = Instant::now();
    let opts = StoreOptions {
        role_prefix: sw.contains(&"--gt".to_string()).then(|| "GT_".to_string()),
        tmp_dir: tmp_dir(&kv, &reads),
    };
    let store = ReadStore::build(&reads, None, &opts).map_err(|e| format!("{}: {e}", reads.display()))?;
    let st = &store.stats;
    eprintln!(
        "[scan] {} records ({} kept), {} loci, {} runs, {} loci in >1 run, max {} reads/locus, spill {} bytes, {:.1}s",
        st.records_seen, st.records_kept, st.loci, st.runs, st.multi_run_loci, st.max_reads_per_locus, st.spill_bytes,
        t0.elapsed().as_secs_f64()
    );
    // read every locus back (in parallel chunks) to exercise the per-locus path
    let t1 = Instant::now();
    let loci = store_loci(&store);
    let counts = map_loci(&loci, 256, |l| {
        let r = store.reads(l).expect("spill read");
        (r.len(), r.iter().map(|x| x.seq.len()).sum::<usize>())
    });
    let n: usize = counts.iter().map(|c| c.0).sum();
    let bases: usize = counts.iter().map(|c| c.1).sum();
    eprintln!("[scan] read back {n} reads / {bases} bases of {} loci in {:.1}s", loci.len(), t1.elapsed().as_secs_f64());
    Ok(())
}

fn store_loci(store: &ReadStore) -> Vec<String> {
    let mut v = store.loci();
    v.sort();
    v
}

fn annotate(args: &[String]) -> Result<(), String> {
    let (kv, _) = parse_flags(
        args,
        &["--config", "--inputs", "--out", "--evidence", "--reads", "--gt-reads", "--tmp-dir", "--threads", "--chunk"],
        &[],
    )?;
    let need = |k: &str| kv.get(k).cloned().ok_or_else(|| format!("annotate: {k} is required\n{USAGE}"));
    let cfg = peartree_rte::config::RteConfig::load(Path::new(&need("--config")?))?;
    let out = PathBuf::from(need("--out")?);
    if let Some(n) = kv.get("--threads") {
        let n: usize = n.parse().map_err(|_| "--threads: not a number")?;
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().map_err(|e| e.to_string())?;
    }
    let t0 = Instant::now();
    let inputs = peartree_rte::io::read_inputs(Path::new(&need("--inputs")?))?;
    let wanted = inputs.iter().map(|i| i.title.clone()).collect();
    let p = |k: &str| kv.get(k).map(PathBuf::from);
    let loader = LocusLoader::open(p("--evidence").as_deref(), p("--reads").as_deref(), p("--gt-reads").as_deref(), &wanted, &tmp_dir(&kv, &out))
        .map_err(|e| format!("reading the sidecars: {e}"))?;
    eprintln!(
        "[rte] {} insertions; evidence rows for {}; {} reads indexed ({} loci, {} in >1 run); genotype reads: {}; {:.1}s",
        inputs.len(),
        loader.evidence.len(),
        loader.reads.stats.records_kept,
        loader.reads.stats.loci,
        loader.reads.stats.multi_run_loci,
        if loader.has_gt_reads { loader.gt.stats.records_kept.to_string() } else { "no file".into() },
        t0.elapsed().as_secs_f64()
    );
    let res = peartree_rte::annotator::Resources::load(&cfg)?;
    let mut ann = peartree_rte::annotator::RteAnnotator::new(&cfg, &res);
    if let Some(c) = kv.get("--chunk") {
        ann.chunk = c.parse::<usize>().map_err(|_| "--chunk: not a number")?.max(1);
    }
    let t1 = Instant::now();
    let records = ann.annotate_all(&inputs, &loader)?;
    let n_err = records.iter().filter(|r| r.detail.contains("error")).count();
    eprintln!(
        "[rte] annotated {} insertions in {:.1}s ({} TPRT, {} errors{})",
        records.len(),
        t1.elapsed().as_secs_f64(),
        records.iter().filter(|r| r.tprt_call == "TPRT").count(),
        n_err,
        if loader.has_gt_reads {
            format!(
                "; genotype reads used at {}, changing a call at {}",
                records.iter().filter(|r| r.gt_reads != 0).count(),
                records.iter().filter(|r| !r.gt_changed.is_empty()).count()
            )
        } else {
            String::new()
        }
    );
    let rows = inputs.iter().zip(&records).map(|(i, r)| peartree_rte::io::output_row(&i.title, Some(r)));
    peartree_rte::io::write_output(&out, rows)?;
    eprintln!("[rte] wrote {} rows -> {} ({:.1}s total)", records.len(), out.display(), t0.elapsed().as_secs_f64());
    Ok(())
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let r = match argv.first().map(|s| s.as_str()) {
        Some("annotate") => annotate(&argv[1..]),
        Some("scan") => scan(&argv[1..]),
        _ => Err(USAGE.to_string()),
    };
    if let Err(e) = r {
        eprintln!("peartree-rte: {e}");
        std::process::exit(2);
    }
}
