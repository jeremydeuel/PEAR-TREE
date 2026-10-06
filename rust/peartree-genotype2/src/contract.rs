//! Contract + combined-consensus parsing. Owner: A.
//!
//! * `load_loci`: `<patient>.genotyping[.tprt].txt.gz` -> ordered loci (+ the 12-bp fallback
//!   consensus per side).
//! * `load_combined`: `<patient>.combined.txt.gz` -> per locus, per side `JunctionConsensus`
//!   (case-split: uppercase = reference flank, lowercase = inserted part; qualities phred+33).

use crate::types::{JunctionConsensus, Locus};
use std::collections::HashMap;
use std::io;

/// A contract entry: the locus plus whatever consensus the contract itself carries (12 bp).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ContractSides {
    pub left: Option<JunctionConsensus>,
    pub right: Option<JunctionConsensus>,
}

pub fn load_loci(path: &str) -> io::Result<Vec<(Locus, ContractSides)>> {
    let _ = path;
    todo!("owner A: contract::load_loci")
}

/// Per locus name: (left, right) consensus from the combined file. Missing side -> None.
pub fn load_combined(path: &str) -> io::Result<HashMap<String, ContractSides>> {
    let _ = path;
    todo!("owner A: contract::load_combined")
}

/// Parse a locus name `<contig>:<L>-<R>` (either end may carry `oneside_`); split from the right
/// so contig names containing ':' or '-' work.
pub fn parse_locus_name(name: &str) -> io::Result<Locus> {
    let _ = name;
    todo!("owner A: contract::parse_locus_name")
}
