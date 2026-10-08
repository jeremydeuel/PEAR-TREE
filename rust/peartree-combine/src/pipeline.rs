//! Stage driver -- python `combine_insertions()` (combine_insertions.py:113-365) and the
//! `--step combine_insertions` branch of src/main.py:148-188. OWNER: P6.
//! SPEC.md §1 (CLI, outputs) and §3-§6 (stages in order).

use crate::config::Config;
use crate::context::Ctx;
use crate::evidence::filters::discovery_breakpoints;
use crate::evidence::output::{evidence_paths, write_evidence_outputs};
use crate::evidence::apply_evidence;
use crate::genome::Genome;
use crate::genotyping_out::genotyping_text;
use crate::insertion::parse_discovery_file;
use crate::intersect::intersect_insertions;
use crate::liftover::LiftOver;
use crate::model::{FileId, InputFile, Insertion, Interner};
use crate::region_filter::filter_dense_regions;
use crate::remap::{
    bowtie2_end_to_end_cmd, bowtie2_local_cmd, clean_remap_names_bam, clip_fastq, clipped_remap_names_bam, consensus_fastq, sh, write_gz,
};
use crate::splice::write_combined_splice;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::path::{Path, PathBuf};

/// Parsed command line (SPEC.md §1.1).
#[derive(Clone, Debug)]
pub struct Args {
    pub discovery_files: Vec<String>,
    /// `--out` with a trailing `.gz` stripped (python main.py)
    pub out_stem: String,
    pub threads: usize,
    pub config: PathBuf,
}

const LONG_OPTS: [&str; 6] = ["step", "out", "discovery_files", "discovery-files", "threads", "config"];

impl Args {
    /// Accepts the python spelling (`--step combine_insertions` ignored / validated,
    /// `--discovery_files F...`, `--out STEM`, `--threads N`) plus `--config PATH`
    /// (default `$PEARTREE_CONFIG`, else `src/config.py` relative to the cwd). Multi-value
    /// `--discovery_files` consumes arguments until the next `--option`. Also accepts
    /// `--discovery-files`.
    pub fn parse(argv: &[String]) -> Result<Args, String> {
        let mut step: Option<String> = None;
        let mut out: Option<String> = None;
        let mut files: Vec<String> = Vec::new();
        let mut files_given = false;
        let mut threads: usize = 1;
        let mut config: Option<String> = None;

        let is_opt = |s: &str| s.len() > 1 && s.starts_with('-') && s.parse::<f64>().is_err();
        let mut i = 0;
        while i < argv.len() {
            let tok = &argv[i];
            i += 1;
            if !is_opt(tok) {
                return Err(format!("unexpected argument {tok:?}"));
            }
            let (key, inline): (String, Option<String>) = match tok.split_once('=') {
                Some((k, v)) if k.starts_with("--") => (k.to_string(), Some(v.to_string())),
                _ => (tok.clone(), None),
            };
            let canon: &str = if let Some(long) = key.strip_prefix("--") {
                if LONG_OPTS.contains(&long) {
                    LONG_OPTS[LONG_OPTS.iter().position(|o| *o == long).unwrap()]
                } else {
                    // argparse accepts unambiguous prefixes
                    let m: Vec<&str> = LONG_OPTS.iter().copied().filter(|o| o.starts_with(long)).collect();
                    let canon_set: FxHashSet<&str> = m.iter().map(|o| if *o == "discovery-files" { "discovery_files" } else { *o }).collect();
                    match (m.len(), canon_set.len()) {
                        (0, _) => return Err(format!("unrecognized option {key}")),
                        (_, 1) => *canon_set.iter().next().unwrap(),
                        _ => return Err(format!("ambiguous option {key}")),
                    }
                }
            } else {
                match key.as_str() {
                    "-s" => "step",
                    "-o" => "out",
                    "-f" => "discovery_files",
                    "-@" => "threads",
                    "-c" => "config",
                    _ => return Err(format!("unrecognized option {key}")),
                }
            };
            let canon = if canon == "discovery-files" { "discovery_files" } else { canon };
            // collect the option's value(s)
            let mut vals: Vec<String> = Vec::new();
            if let Some(v) = inline {
                vals.push(v);
            }
            if canon == "discovery_files" {
                while i < argv.len() && !is_opt(&argv[i]) {
                    vals.push(argv[i].clone());
                    i += 1;
                }
            } else if vals.is_empty() {
                if i < argv.len() && !is_opt(&argv[i]) {
                    vals.push(argv[i].clone());
                    i += 1;
                }
                if vals.is_empty() {
                    return Err(format!("option {key} expects one argument"));
                }
            }
            match canon {
                "step" => step = Some(vals.remove(0)),
                "out" => out = Some(vals.remove(0)),
                "threads" => {
                    threads = vals[0].trim().parse::<usize>().map_err(|_| format!("invalid --threads value {:?}", vals[0]))?;
                }
                "config" => config = Some(vals.remove(0)),
                _ => {
                    files_given = true;
                    files.extend(vals);
                }
            }
        }
        if let Some(s) = &step {
            if s != "combine_insertions" {
                return Err(format!("--step {s}: this binary only implements combine_insertions"));
            }
        }
        if !files_given || files.is_empty() || out.is_none() {
            return Err("usage: peartree-combine [--step combine_insertions] [--config CONFIG] --discovery_files F [F ...] --out STEM [--threads N]".into());
        }
        let mut stem = out.unwrap();
        if stem.ends_with(".gz") {
            stem.truncate(stem.len() - 3);
        }
        let config = config
            .or_else(|| std::env::var("PEARTREE_CONFIG").ok().filter(|s| !s.is_empty()))
            .unwrap_or_else(|| "src/config.py".to_string());
        Ok(Args { discovery_files: files, out_stem: stem, threads, config: PathBuf::from(config) })
    }
}

fn drop_filtered(insertions: &mut Vec<Insertion>, set: &FxHashSet<String>, contigs: &Interner) {
    if !set.is_empty() {
        insertions.retain(|i| !set.contains(&i.name(contigs)));
    }
}

/// Run the whole step. Stage order (SPEC.md §3-§6):
///  1. load config, open genome (error if missing), validate inputs (python main.py checks:
///     files exist, samtools/bowtie2 executables exist, `{bowtie2_index}.1.bt2` exists,
///     threads <= cpu count), reject duplicate basenames;
///  2. parse discovery files in parallel, import loop / exclusion in command-line order (§3.1);
///  3. far_pair_strict -> discovery_breakpoints over all accepted records (§3.2);
///  4. intersect_insertions (§3.3); 5. filter_dense_regions(100, 4) (§3.4);
///  6. evidence (if any accepted sidecar) -> apply_evidence (§4);
///  7. write `<stem>.fq.gz` consensus FASTQ, bowtie2 end-to-end -> `<stem>.bam` unless it exists,
///     clean-remap filter (§5.1);
///  8. unless `<stem>.insertionsonly.bam` exists: write clip FASTQ to `<stem>.fq.gz`, bowtie2
///     local; clipped-remap filter with liftover (§5.2);
///  9. evidence -> absorb_one_sided (§4.7);
/// 10. write `<stem>.combined.txt.gz` (§6.1); 11. `<stem>.combined.splice.tsv` (§6.6);
/// 12. evidence -> write evidence TSV + reads FASTA (§6.4-6.5), cleanup shard dir;
/// 13. write `<stem>.genotyping.txt.gz` (§6.3).


/// Maximum |TSD / target-site deletion| (bp) of a two-sided locus (Jeremy, 2026-10-08).
pub const MAX_SITE_GAP: i64 = 120;

/// A locus passes unless BOTH name tokens are plain junction coordinates more than MAX_SITE_GAP
/// apart (`polyA_` = a mate position, `disc_` / `oneside_` = no second junction: no TSD).
pub fn site_gap_ok(i: &Insertion) -> bool {
    use crate::model::TokKind::Pos;
    !(i.name_start.kind == Pos && i.name_end.kind == Pos && (i.name_end.pos - i.name_start.pos).abs() > MAX_SITE_GAP)
}

pub fn run(args: &Args) -> Result<(), String> {
    crate::diag::start();
    // ---- 1. config, validation, genome
    let cfg = Config::load(&args.config)?;
    let threads = args.threads;
    let cpus = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    if threads > cpus {
        return Err(format!("this machine only has {cpus} CPUs, do not run this script with more threads, you have requested {threads}."));
    }
    if threads == 0 {
        return Err("--threads must be >= 1".into());
    }
    for f in &args.discovery_files {
        if !Path::new(f).exists() {
            return Err(format!("input insertions file {f} does not exist!"));
        }
    }
    if !Path::new(&cfg.samtools_executable).exists() {
        return Err(format!("samtools not found in {}, specify this in the config file", cfg.samtools_executable));
    }
    if !Path::new(&cfg.bowtie2_executable).exists() {
        return Err(format!("bowtie2 not found in {}, specify this in the config file", cfg.bowtie2_executable));
    }
    if !Path::new(&format!("{}.1.bt2", cfg.bowtie2_index)).exists() {
        return Err(format!("bowtie2 index not found in {}, generate this using `{}-build [REFERENCE_GENOME.fasta] {}`", cfg.bowtie2_index, cfg.bowtie2_executable, cfg.bowtie2_index));
    }
    let files: Vec<InputFile> = args.discovery_files.iter().map(|p| InputFile::new(p)).collect();
    {
        let mut seen: FxHashSet<&str> = FxHashSet::default();
        for f in &files {
            if !seen.insert(f.basename.as_str()) {
                return Err(format!("duplicate input basename {:?}: python keys files by basename, rename one of them", f.basename));
            }
        }
    }
    let genome = Genome::open(Path::new(&cfg.genome_2bit))?;
    let ctx = Ctx { cfg, contigs: Interner::new(), files, genome, threads };
    let stem = args.out_stem.as_str();
    let fq = format!("{stem}.fq.gz");
    let bam = format!("{stem}.bam");
    let clipped_bam = format!("{stem}.insertionsonly.bam");
    let combined = format!("{stem}.combined.txt.gz");
    let genotyping = format!("{stem}.genotyping.txt.gz");
    let cfgc = &ctx.cfg;

    // ---- 2. import
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().map_err(|e| format!("thread pool: {e}"))?;
    let imports: Vec<Result<crate::insertion::FileImport, String>> = pool.install(|| {
        ctx.files
            .par_iter()
            .enumerate()
            .map(|(k, f)| parse_discovery_file(Path::new(&f.path), k as FileId, &ctx.contigs))
            .collect()
    });
    let mut all: Vec<Insertion> = Vec::new();
    let mut accepted: Vec<FileId> = Vec::new();
    for (k, imp) in imports.into_iter().enumerate() {
        let imp = imp.map_err(|e| format!("{}: {e}", ctx.files[k].path))?;
        println!("File {}: imported {} insertions", ctx.files[k].basename, imp.records.len());
        if imp.records.len() as i64 > cfgc.exclude_files_with_many_insertions {
            println!("removed file {} since it contains too many insertions.", ctx.files[k].path);
        } else {
            all.extend(imp.records);
            accepted.push(k as FileId);
        }
    }
    // HARD rule: a TSD / target-site deletion is never longer than MAX_SITE_GAP (120 bp). Discovery
    // refuses wider windows since 2026-10-08, but older per-colony files (e.g. PD51635, run with
    // L1-mediated far pairing) still carry 14-22 kb "pairs": drop them here, before anything else.
    let n0 = all.len();
    all.retain(site_gap_ok);
    if all.len() < n0 {
        println!("site gap: removed {} loci with |TSD/target-site deletion| > {MAX_SITE_GAP} bp", n0 - all.len());
    }
    println!("intersecting insertions from {} files...", ctx.files.len());
    crate::diag::memlog("import");

    // ---- 3. per-sample discovery breakpoints (before intersect merges records)
    let breakpoints = if cfgc.far_pair_strict { Some(discovery_breakpoints(&all)) } else { None };

    // ---- 4/5. intersect + dense-region filter (inputs are consumed)
    let insertions = intersect_insertions(all, cfgc, &ctx.contigs, &ctx.files);
    let (bin_range, ins_cutoff) = (100i64, 4usize);
    let (mut insertions, removed, regions) = filter_dense_regions(insertions, bin_range, ins_cutoff);
    println!(
        "filtering regions with very high insertion rate of {ins_cutoff} or higher per {bin_range} bases ,removed {removed} insertions in {regions} regions, {} insertions are remaining",
        insertions.len()
    );

    crate::diag::memlog("intersect + dense filter");
    // ---- 6. evidence
    let shard_dir = PathBuf::from(format!("{stem}.evidence_shards"));
    let (ins2, mut evidence) = apply_evidence(insertions, &accepted, &ctx, breakpoints.as_ref(), &shard_dir)?;
    insertions = ins2;
    drop(breakpoints);
    crate::diag::memlog("apply_evidence");

    // ---- 7. consensus FASTQ + end-to-end remap + clean-remap filter
    println!("writing summarised insertions fasta file {fq}");
    write_gz(Path::new(&fq), &consensus_fastq(&insertions, &ctx), 1)?;
    if !Path::new(&bam).exists() {
        println!("running bowtie2 {} with index {}", cfgc.bowtie2_executable, cfgc.bowtie2_index);
        sh(&bowtie2_end_to_end_cmd(&ctx, &fq, &bam));
    }
    let filter_reads = clean_remap_names_bam(&ctx, Path::new(&bam))?;
    println!(
        "detected {} insertions where at least one end maps cleanly (no inserted block) to the reference genome, removing these (since they can not be chimeric)...",
        filter_reads.len()
    );
    drop_filtered(&mut insertions, &filter_reads, &ctx.contigs);

    crate::diag::memlog("clean remap");
    // ---- 8. clipped-part local remap (the chain loads while bowtie2 runs) + filter
    println!("now re-mapping in local mode all clipped parts of reads");
    let lo_path = cfgc.bowtie2_index2_lo.clone();
    let lo: Result<LiftOver, String> = std::thread::scope(|s| {
        let loader = s.spawn(|| LiftOver::open(Path::new(&lo_path)));
        if !Path::new(&clipped_bam).exists() {
            let res = write_gz(Path::new(&fq), &clip_fastq(&insertions, &ctx), 1);
            if let Err(e) = res {
                return Err(e);
            }
            println!("running bowtie2 {} with index {}", cfgc.bowtie2_executable, cfgc.bowtie2_index2);
            sh(&bowtie2_local_cmd(&ctx, &fq, &clipped_bam));
        }
        loader.join().map_err(|_| "liftover loader panicked".to_string())?
    });
    let lo = lo?;
    println!("removing all reads where one of the clipped ends maps within 1000bp of the breakpoint.");
    let filter_reads = clipped_remap_names_bam(&ctx, Path::new(&clipped_bam), &lo)?;
    drop(lo);
    println!("detected {} insertions where the clipped part maps near the breakpoint. Removing these", filter_reads.len());
    drop_filtered(&mut insertions, &filter_reads, &ctx.contigs);

    crate::diag::memlog("clipped remap");
    // ---- 9. fold surviving one-sided loci
    if let Some(state) = evidence.as_mut() {
        let (ins2, n_abs) = state.absorb_one_sided(insertions, &ctx);
        insertions = ins2;
        if n_abs > 0 {
            println!("folded {n_abs} one-sided loci into a surviving call of the same junction");
        }
    }

    crate::diag::memlog("absorb_one_sided");
    // ---- 10. combined.txt.gz
    write_gz(Path::new(&combined), &consensus_fastq(&insertions, &ctx), 9)?;
    let n_written = insertions.iter().filter(|i| !filter_reads.contains(&i.name(&ctx.contigs))).count();
    println!("wrote {n_written} insertions to insertions.txt.gz");

    // ---- 11. splice sidecar
    write_combined_splice(&insertions, &combined, &ctx, 25)?;

    // ---- 12. evidence outputs
    if let Some(mut state) = evidence {
        let (tsv, fa) = evidence_paths(&combined);
        let mut names: Vec<String> = insertions.iter().map(|i| i.name(&ctx.contigs)).collect();
        let mut failed = state.failed.clone();
        failed.sort();
        names.extend(failed);
        write_evidence_outputs(&state, &names, &tsv, &fa, &ctx)?;
        state.store.cleanup();
    }

    crate::diag::memlog("evidence outputs");
    // ---- 13. genotyping contract
    let text = genotyping_text(&insertions, &filter_reads, &ctx);
    write_gz(Path::new(&genotyping), &text, 9)?;
    crate::diag::memlog("genotyping");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn a(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parse_python_spelling() {
        let r = Args::parse(&a(&["--step", "combine_insertions", "--discovery_files", "a.txt.gz", "b.txt.gz", "--out", "x/stem.gz", "--threads", "4"])).unwrap();
        assert_eq!(r.discovery_files, vec!["a.txt.gz", "b.txt.gz"]);
        assert_eq!(r.out_stem, "x/stem");
        assert_eq!(r.threads, 4);
    }

    #[test]
    fn parse_variants() {
        let r = Args::parse(&a(&["--discovery-files", "a", "b", "-o", "s", "--config=/c.json", "-@", "2"])).unwrap();
        assert_eq!(r.discovery_files, vec!["a", "b"]);
        assert_eq!(r.config, PathBuf::from("/c.json"));
        assert_eq!(r.threads, 2);
        let r = Args::parse(&a(&["--out", "s", "--disc", "a", "--thr", "3"])).unwrap();
        assert_eq!((r.discovery_files.len(), r.threads, r.out_stem.as_str()), (1, 3, "s"));
        // files up to the next option, then more options
        let r = Args::parse(&a(&["-f", "a", "b", "c", "--out", "s"])).unwrap();
        assert_eq!(r.discovery_files.len(), 3);
        assert!(Args::parse(&a(&["--out", "s"])).is_err());
        assert!(Args::parse(&a(&["-f", "a"])).is_err());
        assert!(Args::parse(&a(&["--step", "discover", "-f", "a", "-o", "s"])).is_err());
        assert!(Args::parse(&a(&["-f", "a", "-o", "s", "--bogus", "1"])).is_err());
    }
}
