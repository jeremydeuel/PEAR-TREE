//! Joint phylogenetic genotyping across the colonies of one patient (SPEC "Joint"). Owner: E.
//! Owns newick.rs too.

use std::io;

pub struct JointArgs {
    pub tree: String,
    /// per-colony genotype files (stem = colony id = tree tip label)
    pub genotype_files: Vec<String>,
    pub out_tsv: String,
    pub out_matrix: String,
    pub root_prior: f64,
    /// "length" | "uniform"
    pub branch_prior: String,
}

pub fn run(args: &JointArgs) -> io::Result<()> {
    let _ = args;
    todo!("owner E: joint::run")
}
