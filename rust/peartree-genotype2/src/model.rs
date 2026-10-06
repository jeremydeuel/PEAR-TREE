//! Genotype likelihoods, posterior, vocabulary mapping (SPEC "Genotype model"). Owner: C.

use crate::config::Config;
use crate::types::{Call, ReadObs};

/// Call one locus from its read observations. `n_disc` discordant anchors add
/// `cfg.disc_weight_nats` each towards alt (0 by default). Never returns `high-coverage` /
/// `error` (the driver decides those before calling).
pub fn call_locus(obs: &[ReadObs], n_disc: i64, cfg: &Config) -> Call {
    let _ = (obs, n_disc, cfg);
    todo!("owner C: model::call_locus")
}
