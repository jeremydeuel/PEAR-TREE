//! Per-colony genotyping driver (SPEC "Driver / I/O"). Owner: D. Owns source.rs / read.rs too.

use crate::config::Config;
use crate::contract::ContractSides;
use crate::types::Locus;
use std::io::{self, Write};

/// Genotype `loci` (contract order; the driver sorts internally) against one alignment file,
/// writing rows to `writer` (header included). `reference_path` is the genome (FASTA/2bit) used
/// for haplotypes AND for CRAM decoding. `combined` = the per-locus consensus from `--combined`
/// (None -> contract fallback with one warning).
pub fn run<W: Write>(
    loci: &[(Locus, ContractSides)],
    combined: Option<&std::collections::HashMap<String, ContractSides>>,
    input: &str,
    reference_path: &str,
    cfg: &Config,
    threads: usize,
    writer: &mut W,
) -> io::Result<()> {
    let _ = (loci, combined, input, reference_path, cfg, threads, writer);
    todo!("owner D: driver::run")
}

/// Validate an input alignment file up front (exists, indexed).
pub fn check_input(path: &str) -> Result<(), String> {
    let _ = path;
    todo!("owner D: driver::check_input")
}
