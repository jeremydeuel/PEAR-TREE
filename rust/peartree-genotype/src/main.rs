//! PEAR-TREE genotyping step, Rust port. Drop-in replacement for
//!   python src/main.py --step genotype --bam <bam> --out <out.txt.gz> --insertions <ins.genotyping.txt.gz>
//! Emits the same 8-column gzip table. Byte-identical (decompressed) to the Python
//! oracle under the generic src/config.py genotyping defaults.
//!
//! Two steps:
//!   --step genotype        one alignment file -> one output
//!   --step genotype_batch  a manifest of many files, genotyped CONSECUTIVELY
//!                          (one file at a time, `--threads` cores spread across
//!                          that file's loci). Bounds memory + open handles on the
//!                          cluster: 4 cores over 20 files never opens 20 readers.
//!
//! BAM and CRAM input (indexed: `.bai`/`.csi`, or `.crai` + `--reference`).

mod config;
mod evidence;
mod genotype;
mod insertion;
mod read;
mod source;

use std::fs::File;
use std::io::{self, BufRead, BufWriter};

use flate2::write::GzEncoder;
use flate2::Compression;

use config::GenotypingConfig;
use insertion::Insertion;

fn usage() -> ! {
    eprintln!(
        "usage:\n  \
         peartree-genotype --step genotype --bam <bam|cram> --insertions <ins.genotyping.txt.gz> \
         --out <out.txt.gz> [--threads N] [--config <file>] [--reference <ref.fa> (CRAM)]\n  \
         peartree-genotype --step genotype_batch --manifest <samples.tsv: input<TAB>output per line> \
         --insertions <ins.genotyping.txt.gz> [--threads N] [--config <file>] [--reference <ref.fa> (CRAM)]"
    );
    std::process::exit(1);
}

/// True when a coordinate index sits next to `path` (`.bai`/`.csi` for BAM, `.crai` for CRAM).
fn has_index(path: &str) -> bool {
    if path.ends_with(".cram") {
        std::path::Path::new(&format!("{path}.crai")).exists()
    } else {
        std::path::Path::new(&format!("{path}.bai")).exists()
            || std::path::Path::new(&format!("{path}.csi")).exists()
    }
}

/// Validate one input alignment file up front (exists, indexed, reference for CRAM).
fn check_input(path: &str, reference: Option<&str>) -> Result<(), String> {
    if !std::path::Path::new(path).exists() {
        return Err(format!("input file {path} does not exist"));
    }
    if !has_index(path) {
        let want = if path.ends_with(".cram") { ".crai" } else { ".bai/.csi" };
        return Err(format!("no index found next to {path} (expected {want}); run `samtools index`"));
    }
    if path.ends_with(".cram") && reference.is_none() {
        return Err(format!("CRAM input {path} requires a reference FASTA (--reference <ref.fa>)"));
    }
    Ok(())
}

/// Open a gzip writer, run genotyping for one file, finish the stream.
fn genotype_to_file(
    insertions: &[Insertion],
    input: &str,
    output: &str,
    reference: Option<&str>,
    cfg: &GenotypingConfig,
    threads: usize,
) -> io::Result<()> {
    let file = File::create(output)?;
    let encoder = GzEncoder::new(BufWriter::new(file), Compression::default());
    let mut writer = BufWriter::new(encoder);
    genotype::run(insertions, input, reference, cfg, threads, &mut writer)?;
    writer.into_inner()?.finish()?;
    Ok(())
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let mut bam: Option<String> = None;
    let mut out: Option<String> = None;
    let mut insertions_path: Option<String> = None;
    let mut manifest: Option<String> = None;
    let mut reference: Option<String> = None;
    let mut step: Option<String> = None;
    let mut config_path: Option<String> = None;
    let mut threads: usize = std::env::var("PEARTREE_THREADS").ok().and_then(|s| s.parse().ok()).unwrap_or(0);

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--step" | "-s" => { step = args.get(i + 1).cloned(); i += 2; }
            "--bam" | "-d" => { bam = args.get(i + 1).cloned(); i += 2; }
            "--out" | "-o" => { out = args.get(i + 1).cloned(); i += 2; }
            "--insertions" | "-i" => { insertions_path = args.get(i + 1).cloned(); i += 2; }
            "--manifest" | "-m" => { manifest = args.get(i + 1).cloned(); i += 2; }
            "--reference" | "-T" => { reference = args.get(i + 1).cloned(); i += 2; }
            "--threads" | "-@" => { threads = args.get(i + 1).and_then(|s| s.parse().ok()).unwrap_or(threads); i += 2; }
            "--config" | "-c" => { config_path = args.get(i + 1).cloned(); i += 2; }
            _ => { i += 1; }
        }
    }

    let step = step.as_deref();
    if step != Some("genotype") && step != Some("genotype_batch") {
        eprintln!("this binary implements --step genotype and --step genotype_batch");
        usage();
    }

    let Some(insertions_path) = insertions_path else { usage() };
    if !std::path::Path::new(&insertions_path).exists() {
        eprintln!("input insertions file {insertions_path} does not exist!");
        std::process::exit(1);
    }

    let cfg = match GenotypingConfig::load(config_path.as_deref()) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("config error: {e}");
            std::process::exit(1);
        }
    };

    if threads == 0 {
        threads = std::thread::available_parallelism().map(|n| n.get().min(8)).unwrap_or(1);
    }
    threads = threads.max(1);
    let reference = reference.as_deref();

    eprintln!("PEAR-TREE genotyping (rust)");
    eprintln!("effective min_mapq: {}, reads_for_high_coverage: {}", cfg.min_mapq, cfg.reads_for_high_coverage);

    // Parse the insertion contract ONCE, shared by every sample.
    let insertions = genotype::load_contract(&insertions_path)?;
    let n_one_sided = insertions.iter().filter(|i| i.is_one_sided()).count();
    eprintln!("contract: {} loci ({n_one_sided} one-sided) from {insertions_path}", insertions.len());

    match step {
        Some("genotype") => {
            let (Some(bam), Some(out)) = (bam, out) else { usage() };
            if let Err(e) = check_input(&bam, reference) {
                eprintln!("{e}");
                std::process::exit(1);
            }
            eprintln!("input: {bam} -> {out}  ({threads} thread(s))");
            genotype_to_file(&insertions, &bam, &out, reference, &cfg, threads)?;
        }
        Some("genotype_batch") => {
            let Some(manifest) = manifest else { usage() };
            // Parse the manifest: one "input<TAB>output" per line (# / blank ignored).
            let mut samples: Vec<(String, String)> = Vec::new();
            let mf = File::open(&manifest).unwrap_or_else(|e| {
                eprintln!("cannot read manifest {manifest}: {e}");
                std::process::exit(1);
            });
            for (lineno, line) in io::BufReader::new(mf).lines().enumerate() {
                let line = line?;
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let mut it = line.split('\t');
                let (Some(inp), Some(outp), None) = (it.next(), it.next(), it.next()) else {
                    eprintln!("{manifest}:{}: expected 'input<TAB>output', got {line:?}", lineno + 1);
                    std::process::exit(1);
                };
                samples.push((inp.trim().to_string(), outp.trim().to_string()));
            }
            if samples.is_empty() {
                eprintln!("manifest {manifest} lists no samples");
                std::process::exit(1);
            }
            // Validate every input up front so a bad path fails before hours of work.
            for (inp, _) in &samples {
                if let Err(e) = check_input(inp, reference) {
                    eprintln!("{e}");
                    std::process::exit(1);
                }
            }
            eprintln!("batch: {} samples, CONSECUTIVE, {threads} thread(s)/sample", samples.len());
            // Files are processed ONE AT A TIME (cluster-friendly: bounded memory +
            // open handles); the threads parallelise loci WITHIN each file.
            for (k, (inp, outp)) in samples.iter().enumerate() {
                eprintln!("[{}/{}] {inp} -> {outp}", k + 1, samples.len());
                genotype_to_file(&insertions, inp, outp, reference, &cfg, threads)?;
            }
        }
        _ => usage(),
    }

    eprintln!("done.");
    Ok(())
}
