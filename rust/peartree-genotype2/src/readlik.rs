//! Per-read likelihoods under the locus haplotypes (SPEC "Read likelihood"). Owner: B.

use crate::config::Config;
use crate::types::{LocusModel, ReadInput, ReadObs};

/// Realign one read against every segment of `model` and classify it.
pub fn score_read(model: &LocusModel, read: &ReadInput, cfg: &Config) -> ReadObs {
    let _ = (model, read, cfg);
    todo!("owner B: readlik::score_read")
}
