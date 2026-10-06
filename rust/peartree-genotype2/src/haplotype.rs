//! Per-locus haplotype model construction (SPEC "Per-locus haplotype model"). Owner: A.

use crate::config::Config;
use crate::contract::ContractSides;
use crate::refseq::RefSeq;
use crate::types::{Locus, LocusModel};
use std::io;

/// Build the segments / windows / breakpoints for one locus. `sides` = the best available
/// consensus (combined if present, else contract 12-bp). Never panics on bad input: returns
/// `LocusModel::error(..)` for a locus that cannot be modelled (missing consensus, contig not in
/// the reference, ...). I/O errors from the reference propagate.
pub fn build_model(locus: &Locus, sides: &ContractSides, reference: &mut dyn RefSeq, cfg: &Config) -> io::Result<LocusModel> {
    let _ = (locus, sides, reference, cfg);
    todo!("owner A: haplotype::build_model")
}
