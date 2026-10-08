//! peartree-genotype2 — realignment-based, phylogeny-aware genotyper. See plans/genotype_v2/SPEC.md.
//!
//!   --step genotype        --bam <bam|cram> --insertions <contract> [--combined <P.combined.txt.gz>]
//!                          --reference <fa|2bit> --out <out.txt.gz> [--threads N] [--config <file>]
//!                          [--members <P.members.tsv.gz | P.insertions.reads.fa.gz>] [--sample <id>]
//!                          (with `gt_extra_reads = true`: also <out>.extra_reads.fa.gz, extra.rs)
//!   --step genotype_batch  --manifest <input<TAB>output per line> --insertions .. --reference .. [..]
//!   --step joint           --tree <newick> (--genotype-dir <dir> | --genotypes f1 f2 ..)
//!                          --out <P.joint.tsv> --matrix <P.joint_matrix.csv.gz>
//!                          [--root-prior 0.1] [--branch-prior length|uniform] [--dropout 0.02] [--false-present 0] [--noise-max-frac 1.0] [--ref-bias off|auto|<b>] [--zygosity colony|locus]

mod align;
mod config;
mod contract;
mod driver;
mod extra;
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
         [--threads N] [--config <file>] [--members <members.tsv.gz|reads.fa.gz>] [--sample <id>]\n  \
         peartree-genotype2 --step genotype_batch --manifest <samples.tsv> --insertions <contract> \
         [--combined ..] --reference <ref> [--threads N] [--config <file>]\n  \
         peartree-genotype2 --step joint --tree <newick> (--genotype-dir <dir> | --genotypes f1 f2 ..) \
         --out <P.joint.tsv> --matrix <P.joint_matrix.csv.gz> [--root-prior 0.1] [--branch-prior length|uniform] [--dropout 0.02] [--false-present 0] [--noise-max-frac 1.0] [--ref-bias off|auto|<b>] [--zygosity colony|locus]"
    );
    std::process::exit(1);
}

fn die(msg: impl std::fmt::Display) -> ! {
    eprintln!("{msg}");
    std::process::exit(1);
}

/// `members` = (path, sample) of the extra pass (`gt_extra_reads`; None = off): the sidecar
/// `<output>.extra_reads.fa.gz` is written next to the output, streamed like the rows.
#[allow(clippy::too_many_arguments)]
fn genotype_to_file(
    loci: &[(types::Locus, contract::ContractSides)],
    combined: Option<&std::collections::HashMap<String, contract::ContractSides>>,
    input: &str,
    output: &str,
    reference: &str,
    cfg: &Config,
    threads: usize,
    members: Option<(&str, &str)>,
) -> io::Result<()> {
    let file = File::create(output)?;
    let encoder = GzEncoder::new(BufWriter::new(file), Compression::default());
    let mut writer = BufWriter::new(encoder);
    match members {
        None => driver::run(loci, combined, input, reference, cfg, threads, &mut writer, None)?,
        Some((path, sample)) => {
            let set = extra::load_members(path, sample, cfg.gt_extra_max_member_frac)?;
            let in_contract = loci.iter().filter(|(l, _)| set.mine.contains(&l.name)).count();
            let germ = loci.iter().filter(|(l, _)| set.germline.contains(&l.name)).count();
            let side = extra::sidecar_path(output);
            eprintln!(
                "extra reads: {sample} is a discovery member of {in_contract}/{} contract loci ({path}); the others with \
                 ALT support -> {side}",
                loci.len()
            );
            match set.n_colonies {
                Some(n) => eprintln!(
                    "extra reads: {germ} contract loci discovered in > {} of the {n} colonies are skipped as germline",
                    cfg.gt_extra_max_member_frac
                ),
                None => eprintln!("WARNING: {path} has no '#colonies N' line: no germline skip in the extra pass"),
            }
            let xenc = GzEncoder::new(BufWriter::new(File::create(&side)?), Compression::default());
            let mut xw = BufWriter::new(xenc);
            let args = driver::ExtraArgs { members: &set, sample, out: &mut xw };
            driver::run(loci, combined, input, reference, cfg, threads, &mut writer, Some(args))?;
            xw.into_inner()?.finish()?;
        }
    }
    writer.into_inner()?.finish()?;
    Ok(())
}

/// The extra pass's (members file, sample) when `gt_extra_reads` is on and `--members` given.
fn extra_members<'a>(cfg: &Config, members: Option<&'a str>, sample: &'a str) -> Option<(&'a str, &'a str)> {
    if !cfg.gt_extra_reads {
        return None;
    }
    match members {
        Some(m) => Some((m, sample)),
        None => {
            eprintln!("WARNING: gt_extra_reads is on but no --members file: the extra pass is skipped");
            None
        }
    }
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let mut step: Option<String> = None;
    let mut bam: Option<String> = None;
    let mut out: Option<String> = None;
    let mut insertions: Option<String> = None;
    let mut combined: Option<String> = None;
    let mut members: Option<String> = None;
    let mut sample: Option<String> = None;
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
    let mut ref_bias = joint::RefBias::Off;
    let mut zygosity = joint::Zygosity::Colony;
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
            "--members" => { members = Some(next(i)); i += 2; }
            "--sample" => { sample = Some(next(i)); i += 2; }
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
            "--zygosity" => { zygosity = joint::Zygosity::parse(&next(i)).unwrap_or_else(|e| die(e)); i += 2; }
            "--ref-bias" => { ref_bias = joint::RefBias::parse(&next(i)).unwrap_or_else(|e| die(e)); i += 2; }
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
            joint::run(&joint::JointArgs { tree, genotype_files, out_tsv, out_matrix, root_prior, branch_prior, dropout, false_present, noise_max_frac, ref_bias, zygosity })?;
            eprintln!("done.");
            return Ok(());
        }
        Some("genotype") | Some("genotype_batch") => {}
        _ => usage(),
    }

    let Some(insertions) = insertions else { usage() };
    let Some(reference) = reference else { die("--reference <fa|2bit> is required (haplotype flanks; CRAM decoding)") };
    for p in [&insertions, &reference].into_iter().chain(combined.iter()).chain(members.iter()) {
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
            let sample = sample.unwrap_or_else(|| extra::sample_from_out(&out));
            let xm = extra_members(&cfg, members.as_deref(), &sample);
            genotype_to_file(&loci, combined_map.as_ref(), &bam, &out, &reference, &cfg, threads, xm)?;
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
                // batch: the sample is always the output's stem (--sample names one colony)
                let smp = extra::sample_from_out(outp);
                let xm = extra_members(&cfg, members.as_deref(), &smp);
                genotype_to_file(&loci, combined_map.as_ref(), inp, outp, &reference, &cfg, threads, xm)?;
            }
        }
        _ => usage(),
    }
    eprintln!("done.");
    Ok(())
}
