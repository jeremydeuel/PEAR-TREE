//! PEAR-TREE discovery step, Rust port. Drop-in replacement for
//!   python src/main.py --step discover --bam <bam> --out <out.txt.gz>
//! Emits the same custom FASTQ-based .txt.gz breakpoint format.

mod config;
mod discovery;
mod filters;
mod model;
mod polya;
mod qseq;
mod read;

use std::fs::File;
use std::io::{self, BufWriter};

use flate2::write::GzEncoder;
use flate2::Compression;

use discovery::Discovery;

fn usage() -> ! {
    eprintln!("usage: peartree-discovery --step discover --bam <bam> --out <out.txt.gz> [--threads N]");
    std::process::exit(1);
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let mut bam: Option<String> = None;
    let mut out: Option<String> = None;
    let mut step: Option<String> = None;
    let mut threads: usize = std::env::var("PEARTREE_BAM_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--step" | "-s" => { step = args.get(i + 1).cloned(); i += 2; }
            "--bam" | "-d" => { bam = args.get(i + 1).cloned(); i += 2; }
            "--out" | "-o" => { out = args.get(i + 1).cloned(); i += 2; }
            "--threads" | "-@" => { threads = args.get(i + 1).and_then(|s| s.parse().ok()).unwrap_or(threads); i += 2; }
            _ => { i += 1; }
        }
    }

    if step.as_deref() != Some("discover") {
        eprintln!("this binary only implements --step discover");
        usage();
    }
    let (Some(bam), Some(out)) = (bam, out) else { usage() };

    if !std::path::Path::new(&bam).exists() {
        eprintln!("input bam file {bam} does not exist!");
        std::process::exit(1);
    }

    eprintln!("PEAR-TREE discovery (rust)");
    eprintln!("input bam: {bam}, output file: {out}");

    let mut d = Discovery::new(bam, threads);
    d.discovery()?;

    let file = File::create(&out)?;
    let encoder = GzEncoder::new(BufWriter::new(file), Compression::default());
    let mut writer = BufWriter::new(encoder);
    d.output(&mut writer)?;
    writer.into_inner()?.finish()?;

    eprintln!("done.");
    Ok(())
}
