//! PEAR-TREE discovery step, Rust port. Drop-in replacement for
//!   python src/main.py --step discover --bam <bam> --out <out.txt.gz>
//! Emits the same custom FASTQ-based .txt.gz breakpoint format.

mod config;
mod coverage;
mod discovery;
mod evidence;
mod exons;
mod filters;
mod intervals;
mod mem;
mod model;
mod polya;
mod qseq;
mod read;
mod stats;

use std::fs::File;
use std::io::{self, BufWriter, Write};

use flate2::write::GzEncoder;
use flate2::Compression;

use config::DiscoveryConfig;
use discovery::Discovery;

// jemalloc global allocator (default feature; disable with `--no-default-features`).
// glibc's malloc holds freed memory from discovery's millions of tiny per-read/
// per-breakpoint allocations at the process high-water mark: measured 33 GB RSS for
// ~1 GB live data on the deepest colonies, which is the rare 16 GB OOM. jemalloc
// packs the same run to ~2.5 GB.
#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

// Bake in the tuning that produced the 2.5 GB result, so a plain `bsub ... binary`
// gets it without every array task needing MALLOC_CONF in its environment (a runtime
// MALLOC_CONF env var still overrides this if set). `background_thread:true` runs
// page purging on a background thread so the aggressive `*_decay_ms:0` (return dirty
// and muzzy pages to the OS immediately) costs little on the allocation hot path.
#[cfg(feature = "jemalloc")]
#[allow(non_upper_case_globals)]
#[export_name = "malloc_conf"]
pub static malloc_conf: &[u8] = b"background_thread:true,dirty_decay_ms:0,muzzy_decay_ms:0\0";

fn usage() -> ! {
    eprintln!("usage: peartree-discovery --step discover --bam <bam|cram> --out <out.txt.gz> [--threads N] [--config <file>] [--reference <ref.fa> (required for CRAM)]");
    std::process::exit(1);
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let mut bam: Option<String> = None;
    let mut out: Option<String> = None;
    let mut step: Option<String> = None;
    let mut config_path: Option<String> = None;
    let mut reference: Option<String> = None;
    // 0 = auto (resolved from available parallelism after arg parsing). An explicit
    // `--threads`/`-@` or PEARTREE_BAM_THREADS overrides, including `--threads 1`.
    let mut threads: usize = std::env::var("PEARTREE_BAM_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--step" | "-s" => { step = args.get(i + 1).cloned(); i += 2; }
            "--bam" | "-d" => { bam = args.get(i + 1).cloned(); i += 2; }
            "--out" | "-o" => { out = args.get(i + 1).cloned(); i += 2; }
            "--threads" | "-@" => { threads = args.get(i + 1).and_then(|s| s.parse().ok()).unwrap_or(threads); i += 2; }
            "--config" | "-c" => { config_path = args.get(i + 1).cloned(); i += 2; }
            "--reference" | "-T" => { reference = args.get(i + 1).cloned(); i += 2; }
            _ => { i += 1; }
        }
    }

    let step = step.as_deref().unwrap_or("discover").to_string();
    if step != "discover" && step != "coverage-median" {
        eprintln!("unknown step '{step}'; this binary implements --step discover | coverage-median");
        usage();
    }
    let Some(bam) = bam else { usage() };

    if !std::path::Path::new(&bam).exists() {
        eprintln!("input bam file {bam} does not exist!");
        std::process::exit(1);
    }

    let mut config = match DiscoveryConfig::load(config_path.as_deref()) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("config error: {e}");
            std::process::exit(1);
        }
    };

    // SPD-6: resolve the thread count (auto = min(8, available cores)) and route it to
    // the right parallelism knob. On indexed BAM, `--threads N` drives contig-level
    // parallelism (SPD-3, the big win); `discovery()` falls back to the single-threaded
    // scan with N-way BGZF decode when there is no index, and CRAM stays single-threaded
    // until native CRAM parallelism lands. An explicit `contig_threads` in a config file
    // wins. Output is byte-identical either way.
    if threads == 0 {
        threads = std::thread::available_parallelism().map(|n| n.get().min(8)).unwrap_or(1);
    }
    threads = threads.max(1);
    // SPD-6: on indexed BAM, `--threads N` drives contig-level parallel extract (SPD-3);
    // the mate pass and CRAM stay single-threaded (with N-way BGZF decode for BAM). An
    // explicit `contig_threads` in a config file wins. Output is byte-identical either way.
    if !bam.ends_with(".cram") && config.contig_threads <= 1 && threads > 1 {
        config.contig_threads = threads;
    }

    // `--step coverage-median`: run only the coverage pre-pass on the full BAM and print
    // `<bam>\t<median>`. Used to precompute the true genome median that region-slice runs
    // pin via PEARTREE_COVERAGE_MEDIAN. No output file, no exclude/rm/exon/rte setup.
    if step == "coverage-median" {
        let mut d = Discovery::new(bam.clone(), threads, config, None, None);
        d.set_reference_path(reference);
        let med = d.compute_coverage_median()?;
        println!("{bam}\t{med}");
        return Ok(());
    }
    let Some(out) = out else { usage() };

    // SPEC-5: build the exclude-BED interval index up front so a bad path fails fast.
    let exclude = match config.exclude_bed.as_deref() {
        Some(p) => match intervals::IntervalIndex::from_bed(p) {
            Ok(ix) => Some(ix),
            Err(e) => {
                eprintln!("cannot read exclude_bed {p}: {e}");
                std::process::exit(1);
            }
        },
        None => None,
    };

    // SPEC-7: build the young-RepeatMasker mask (divergence-gated) when enabled.
    let rm_mask = if config.rm_self_mask {
        match config.rm_track.as_deref() {
            Some(p) => match intervals::IntervalIndex::from_repeatmasker(p, config.rm_divergence_max) {
                Ok(ix) => Some(ix),
                Err(e) => {
                    eprintln!("cannot read rm_track {p}: {e}");
                    std::process::exit(1);
                }
            },
            None => {
                eprintln!("rm_self_mask is on but rm_track is not set");
                std::process::exit(1);
            }
        }
    } else {
        None
    };

    eprintln!("PEAR-TREE discovery (rust)");
    eprintln!("input bam: {bam}, output file: {out}");
    // min_mapq is the parameter most likely to differ from production (generic
    // config = 40, config_hs/mm = 60); report the effective value for the record.
    eprintln!("effective min_mapq: {}", config.min_mapq);
    match &config.contig_allowlist {
        Some(set) => eprintln!("contig allowlist: {} contigs", set.len()),
        None => eprintln!("contig allowlist: none (legacy len<=5 + not-MT filter)"),
    }
    // Report the slippage gate's EFFECTIVE state. Every setting that changes what is emitted
    // should be visible in this banner: an A/B of this gate once returned "cuts 0.0%" because
    // the binary predated the feature and ignored the key, and nothing in the output said so.
    // A filter you cannot see is a filter you cannot trust you ran.
    if config.slippage_filter {
        eprintln!(
            "slippage filter: ON (min_ref_run {}, min_clip_frac {}, max_period {})",
            config.slippage_min_ref_run, config.slippage_min_clip_frac, config.slippage_max_period
        );
    } else {
        eprintln!("slippage filter: OFF");
    }
    if config.exclude_bed.is_some() {
        eprintln!("exclude-bed: {}", config.exclude_bed.as_deref().unwrap());
    }
    if config.coverage_mask {
        eprintln!("coverage mask: ON (> {}x median)", config.coverage_mask_multiplier);
    }
    if config.adaptive_evidence {
        eprintln!("adaptive evidence floor: ON");
    }
    if config.rm_self_mask {
        eprintln!("RM self-mask: ON (divergence <= {})", config.rm_divergence_max);
    }
    if config.contig_threads > 1 {
        let indexed = std::path::Path::new(&format!("{bam}.bai")).exists()
            || std::path::Path::new(&format!("{bam}.csi")).exists();
        if indexed {
            eprintln!("contig parallelism: {} threads (SPD-3 extract)", config.contig_threads);
        } else {
            eprintln!(
                "contig parallelism requested ({} threads) but no .bai/.csi index found — \
                 falling back to single-threaded scan with {}-way BGZF decode",
                config.contig_threads, threads
            );
        }
    }
    if config.discordant_anchor {
        eprintln!(
            "discordant anchoring: ON (>= {} reads{})",
            config.discordant_min_reads,
            if config.discordant_rte_only {
                format!(", RTE-origin gate >= {}", config.discordant_rte_min)
            } else if config.discordant_rte_track.is_some() {
                ", RTE-origin labelled".to_string()
            } else {
                String::new()
            }
        );
    }

    // D3: build the discordant mate-origin RTE index (reuses the SPEC-7 RM loader).
    let discordant_rte = match (config.discordant_anchor, config.discordant_rte_track.as_deref()) {
        (true, Some(p)) => match intervals::IntervalIndex::from_repeatmasker(p, config.discordant_rte_divergence_max) {
            Ok(ix) => Some(ix),
            Err(e) => {
                eprintln!("cannot read discordant_rte_track {p}: {e}");
                std::process::exit(1);
            }
        },
        _ => None,
    };
    if config.discordant_rte_only && discordant_rte.is_none() {
        eprintln!("discordant_rte_only is on but discordant_rte_track is not set");
        std::process::exit(1);
    }

    // D5: build the exon model for the splice / processed-pseudogene annotation.
    let exon_model = match (config.splice_hallmark, config.exon_annotation.as_deref()) {
        (true, Some(p)) => match exons::GeneModel::load(p) {
            Ok(m) => Some(m),
            Err(e) => {
                eprintln!("cannot read exon_annotation {p}: {e}");
                std::process::exit(1);
            }
        },
        (true, None) => {
            eprintln!("splice_hallmark is on but exon_annotation is not set");
            std::process::exit(1);
        }
        _ => None,
    };

    // CRAM input requires a reference FASTA (to decode read sequences).
    if bam.ends_with(".cram") && reference.is_none() {
        eprintln!("CRAM input requires a reference FASTA: pass --reference <ref.fa> (names matching the CRAM @SQ)");
        std::process::exit(1);
    }

    // TPRT-hallmark keys: report the effective state (every emission-changing key must be
    // visible in the banner).
    let sidecar_on = config.evidence_sidecar;
    if config.ignore_dup_flag {
        eprintln!("dup flag: IGNORED (0x400 reads kept in clip path, mate pass and sidecar)");
    }
    if let Some(n) = config.min_evidence_fragments_per_sample {
        eprintln!("evidence floor: {n} distinct fragment(s) per breakpoint (fragment mode)");
    }
    if config.drop_dup_in_polya_path {
        eprintln!("dup flag: also dropped in the low-MAPQ poly-A path");
    }
    if config.extra_pairing() || config.clip_slippage_junction_spare > 0 {
        eprintln!(
            "pairing modes: target-site deletion <= {} bp, blunt {}, L1-mediated span {} bp, one-sided {} (>= {} fragments, poly-A only {}), SPEC-8b junction spare {} bp",
            config.max_target_site_deletion,
            config.allow_blunt_pairs,
            config.max_l1_mediated_span,
            config.one_sided_loci,
            config.one_sided_min_fragments,
            config.one_sided_require_polya,
            config.clip_slippage_junction_spare
        );
    }
    if config.short_overhang_evidence {
        if sidecar_on {
            eprintln!(
                "short-overhang evidence: ON (window +-{} bp, overhang <= {} bp, cap {} per breakpoint side)",
                config.short_overhang_window, config.short_overhang_max, config.max_short_per_breakpoint
            );
        } else {
            eprintln!("warning: short_overhang_evidence has no effect without evidence_sidecar");
        }
    }
    if sidecar_on {
        eprintln!(
            "evidence sidecar: ON (fetch_all_mates {}, caps: {} reads / {} mates per breakpoint side, disc span {} bp)",
            config.fetch_all_mates, config.max_evidence_reads_per_breakpoint, config.max_mates_per_breakpoint, config.sidecar_disc_span
        );
    } else if config.fetch_all_mates {
        eprintln!("warning: fetch_all_mates has no effect without evidence_sidecar");
    }

    let mut d = Discovery::new(bam, threads, config, exclude, rm_mask);
    d.set_discordant_rte(discordant_rte);
    d.set_exon_model(exon_model);
    d.set_reference_path(reference);
    d.discovery()?;
    mem::phase("after discovery (extract+find_mates+cluster)");

    let file = File::create(&out)?;
    let encoder = GzEncoder::new(BufWriter::new(file), Compression::default());
    let mut writer = BufWriter::new(encoder);
    let mut hallmarks: Vec<u8> = Vec::new();
    // TPRT: per-read evidence sidecar `<out>.evidence.tsv.gz` (only when enabled). The
    // `.txt.gz` writer above is untouched by it.
    let ev_path = format!("{out}.evidence.tsv.gz");
    let mut ev_writer = if sidecar_on {
        let f = File::create(&ev_path)?;
        let mut w = BufWriter::new(GzEncoder::new(BufWriter::new(f), Compression::default()));
        w.write_all(evidence::HEADER.as_bytes())?;
        Some(w)
    } else {
        None
    };
    let names = if sidecar_on { d.reference_names()? } else { Vec::new() };
    let mut sidecar = ev_writer
        .as_mut()
        .map(|w| evidence::Sidecar { w, names, stats: Default::default() });
    d.output(&mut writer, &mut hallmarks, sidecar.as_mut())?;
    mem::phase("after output");
    // Feature A: append discordant-anchored calls (no-op unless discordant_anchor).
    d.discordant_rescue(&mut writer, sidecar.as_mut())?;
    writer.into_inner()?.finish()?;
    mem::phase("after rescue+flush");
    if let Some(mut sc) = sidecar {
        eprintln!("{}", sc.stats.summary());
    }
    if let Some(w) = ev_writer {
        w.into_inner()?.finish()?.flush()?;
        eprintln!("evidence: {ev_path}");
    }

    // OBS-1: reject-counter sidecar next to the output.
    let stats_path = format!("{out}.stats.json");
    std::fs::write(&stats_path, d.stats_json())?;
    eprintln!("stats: {stats_path}");

    // SENS-5: hallmark annotation sidecar (only when enabled).
    if !hallmarks.is_empty() {
        let hm_path = format!("{out}.hallmarks.tsv");
        std::fs::write(&hm_path, &hallmarks)?;
        eprintln!("hallmarks: {hm_path}");
    }

    // D5: splice / processed-pseudogene annotation sidecar (only when enabled).
    let splice = d.splice_annotate()?;
    mem::phase("after splice_annotate");
    if !splice.is_empty() {
        let sp_path = format!("{out}.splice.tsv");
        std::fs::write(&sp_path, &splice)?;
        eprintln!("splice: {sp_path}");
    }

    eprintln!("done.");
    Ok(())
}
