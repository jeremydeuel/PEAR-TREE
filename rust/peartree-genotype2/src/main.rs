//! peartree-genotype2 — realignment-based, phylogeny-aware genotyper. See plans/genotype_v2/SPEC.md.
//!
//!   --step genotype        --bam <bam|cram> --insertions <contract> [--combined <P.combined.txt.gz>]
//!                          --reference <fa|2bit> --out <out.txt.gz> [--threads N] [--config <file>]
//!   --step genotype_batch  --manifest <input<TAB>output per line> --insertions .. --reference .. [..]
//!   --step joint           --tree <newick> (--genotype-dir <dir> | --genotypes f1 f2 ..)
//!                          --out <P.joint.tsv> --matrix <P.joint_matrix.csv.gz>
//!                          [--root-prior 0.1] [--branch-prior length|uniform] [--dropout 0.02] [--false-present 0] [--noise-max-frac 1.0]

mod align;
mod config;
mod contract;
mod driver;
mod haplotype;
mod joint;
mod model;
mod newick;
mod output;
mod read;
mod readlik;
mod refseq;
mod source;
mod types;

use std::fs::File;
use std::io::{self, BufRead, BufWriter};

use flate2::write::GzEncoder;
use flate2::Compression;

use config::Config;

fn usage() -> ! {
    eprintln!(
        "usage:\n  \
         peartree-genotype2 --step genotype --bam <bam|cram> --insertions <contract.txt.gz> \
         [--combined <P.combined.txt.gz>] --reference <ref.fa|ref.2bit> --out <out.txt.gz> \
         [--threads N] [--config <file>]\n  \
         peartree-genotype2 --step genotype_batch --manifest <samples.tsv> --insertions <contract> \
         [--combined ..] --reference <ref> [--threads N] [--config <file>]\n  \
         peartree-genotype2 --step joint --tree <newick> (--genotype-dir <dir> | --genotypes f1 f2 ..) \
         --out <P.joint.tsv> --matrix <P.joint_matrix.csv.gz> [--root-prior 0.1] [--branch-prior length|uniform] [--dropout 0.02] [--false-present 0] [--noise-max-frac 1.0]"
    );
    std::process::exit(1);
}

fn die(msg: impl std::fmt::Display) -> ! {
    eprintln!("{msg}");
    std::process::exit(1);
}

fn genotype_to_file(
    loci: &[(types::Locus, contract::ContractSides)],
    combined: Option<&std::collections::HashMap<String, contract::ContractSides>>,
    input: &str,
    output: &str,
    reference: &str,
    cfg: &Config,
    threads: usize,
) -> io::Result<()> {
    let file = File::create(output)?;
    let encoder = GzEncoder::new(BufWriter::new(file), Compression::default());
    let mut writer = BufWriter::new(encoder);
    driver::run(loci, combined, input, reference, cfg, threads, &mut writer)?;
    writer.into_inner()?.finish()?;
    Ok(())
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let mut step: Option<String> = None;
    let mut bam: Option<String> = None;
    let mut out: Option<String> = None;
    let mut insertions: Option<String> = None;
    let mut combined: Option<String> = None;
    let mut manifest: Option<String> = None;
    let mut reference: Option<String> = None;
    let mut config_path: Option<String> = None;
    let mut tree: Option<String> = None;
    let mut genotype_dir: Option<String> = None;
    let mut genotype_files: Vec<String> = Vec::new();
    let mut matrix: Option<String> = None;
    let mut root_prior: f64 = 0.1;
    let mut dropout: f64 = 0.02;
    let mut false_present: f64 = 0.0;
    let mut noise_max_frac: f64 = 1.0;
    let mut branch_prior: String = "length".to_string();
    let mut threads: usize = std::env::var("PEARTREE_THREADS").ok().and_then(|s| s.parse().ok()).unwrap_or(1);

    let mut i = 1;
    while i < args.len() {
        let next = |i: usize| args.get(i + 1).cloned().unwrap_or_else(|| usage());
        match args[i].as_str() {
            "--step" | "-s" => { step = Some(next(i)); i += 2; }
            "--bam" | "-d" => { bam = Some(next(i)); i += 2; }
            "--out" | "-o" => { out = Some(next(i)); i += 2; }
            "--insertions" | "-i" => { insertions = Some(next(i)); i += 2; }
            "--combined" => { combined = Some(next(i)); i += 2; }
            "--manifest" | "-m" => { manifest = Some(next(i)); i += 2; }
            "--reference" | "-T" => { reference = Some(next(i)); i += 2; }
            "--threads" | "-@" => { threads = next(i).parse().unwrap_or_else(|_| usage()); i += 2; }
            "--config" | "-c" => { config_path = Some(next(i)); i += 2; }
            "--tree" => { tree = Some(next(i)); i += 2; }
            "--genotype-dir" => { genotype_dir = Some(next(i)); i += 2; }
            "--genotypes" => {
                i += 1;
                while i < args.len() && !args[i].starts_with("--") {
                    genotype_files.push(args[i].clone());
                    i += 1;
                }
            }
            "--matrix" => { matrix = Some(next(i)); i += 2; }
            "--root-prior" => { root_prior = next(i).parse().unwrap_or_else(|_| usage()); i += 2; }
            "--branch-prior" => { branch_prior = next(i); i += 2; }
            "--dropout" => { dropout = next(i).parse().unwrap_or_else(|_| usage()); i += 2; }
            "--false-present" => { false_present = next(i).parse().unwrap_or_else(|_| usage()); i += 2; }
            "--noise-max-frac" => { noise_max_frac = next(i).parse().unwrap_or_else(|_| usage()); i += 2; }
            "-h" | "--help" => usage(),
            other => die(format!("unknown argument {other}")),
        }
    }
    let threads = threads.max(1);

    match step.as_deref() {
        Some("joint") => {
            let Some(tree) = tree else { usage() };
            if let Some(dir) = genotype_dir {
                let mut v: Vec<String> = std::fs::read_dir(&dir)?
                    .filter_map(|e| e.ok())
                    .map(|e| e.path().to_string_lossy().into_owned())
                    .filter(|p| p.ends_with(".txt.gz"))
                    .collect();
                v.sort();
                genotype_files.extend(v);
            }
            if genotype_files.is_empty() {
                die("joint: no genotype files (--genotype-dir or --genotypes)");
            }
            let (Some(out_tsv), Some(out_matrix)) = (out, matrix) else { usage() };
            eprintln!("PEAR-TREE joint phylogenetic genotyping (rust): {} colonies", genotype_files.len());
            joint::run(&joint::JointArgs { tree, genotype_files, out_tsv, out_matrix, root_prior, branch_prior, dropout, false_present, noise_max_frac })?;
            eprintln!("done.");
            return Ok(());
        }
        Some("genotype") | Some("genotype_batch") => {}
        _ => usage(),
    }

    let Some(insertions) = insertions else { usage() };
    let Some(reference) = reference else { die("--reference <fa|2bit> is required (haplotype flanks; CRAM decoding)") };
    for p in [&insertions, &reference].into_iter().chain(combined.iter()) {
        if !std::path::Path::new(p).exists() {
            die(format!("input file {p} does not exist"));
        }
    }
    let cfg = Config::load(config_path.as_deref()).unwrap_or_else(|e| die(format!("config error: {e}")));

    eprintln!("PEAR-TREE genotyping v2 (rust): realignment genotyper");
    eprintln!("effective min_mapq: {} (clipped reads >= {}), reads_for_high_coverage: {}, flank: {}",
        cfg.min_mapq, cfg.min_mapq_clipped, cfg.reads_for_high_coverage, cfg.flank);

    let loci = contract::load_loci(&insertions)?;
    let n_one = loci.iter().filter(|(l, _)| l.is_one_sided()).count();
    eprintln!("contract: {} loci ({n_one} one-sided) from {insertions}", loci.len());
    let combined_map = match &combined {
        Some(p) => {
            let m = contract::load_combined(p)?;
            let hit = loci.iter().filter(|(l, _)| m.contains_key(&l.name)).count();
            eprintln!("combined consensus: {} records, {hit}/{} contract loci covered ({p})", m.len(), loci.len());
            Some(m)
        }
        None => {
            eprintln!("WARNING: no --combined file: junction consensus limited to the contract's 12 bp per side");
            None
        }
    };

    match step.as_deref() {
        Some("genotype") => {
            let (Some(bam), Some(out)) = (bam, out) else { usage() };
            driver::check_input(&bam).unwrap_or_else(|e| die(e));
            eprintln!("input: {bam} -> {out}  ({threads} thread(s))");
            genotype_to_file(&loci, combined_map.as_ref(), &bam, &out, &reference, &cfg, threads)?;
        }
        Some("genotype_batch") => {
            let Some(manifest) = manifest else { usage() };
            let mut samples: Vec<(String, String)> = Vec::new();
            let mf = File::open(&manifest).unwrap_or_else(|e| die(format!("cannot read manifest {manifest}: {e}")));
            for (lineno, line) in io::BufReader::new(mf).lines().enumerate() {
                let line = line?;
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let mut it = line.split('\t');
                let (Some(inp), Some(outp), None) = (it.next(), it.next(), it.next()) else {
                    die(format!("{manifest}:{}: expected 'input<TAB>output', got {line:?}", lineno + 1));
                };
                samples.push((inp.trim().to_string(), outp.trim().to_string()));
            }
            if samples.is_empty() {
                die(format!("manifest {manifest} lists no samples"));
            }
            for (inp, _) in &samples {
                driver::check_input(inp).unwrap_or_else(|e| die(e));
            }
            eprintln!("batch: {} samples, CONSECUTIVE, {threads} thread(s)/sample", samples.len());
            for (k, (inp, outp)) in samples.iter().enumerate() {
                eprintln!("[{}/{}] {inp} -> {outp}", k + 1, samples.len());
                genotype_to_file(&loci, combined_map.as_ref(), inp, outp, &reference, &cfg, threads)?;
            }
        }
        _ => usage(),
    }
    eprintln!("done.");
    Ok(())
}
