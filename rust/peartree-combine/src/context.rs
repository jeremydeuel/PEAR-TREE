//! Run-wide shared state. IMPLEMENTED (architect).

use crate::config::Config;
use crate::genome::Genome;
use crate::model::{InputFile, Interner};

/// Everything the stages share (read-only after startup, `Sync`).
pub struct Ctx {
    pub cfg: Config,
    pub contigs: Interner,
    /// `--discovery_files` in command-line order (FileId = index)
    pub files: Vec<InputFile>,
    /// opened at startup like python (combine_insertions imports
    /// combine_insertions_get_sequence, which opens genome_2bit at import time)
    pub genome: Genome,
    pub threads: usize,
}
