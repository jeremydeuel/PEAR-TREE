//! tools/rte/score.py -- the transparent additive TPRT point system.
//!
//! STATUS: ScoreInput FOUNDATION; WEIGHTS / THRESHOLDS / score = WORK PACKAGE "WP-HALL".
//!
//! Port notes: `total = round(sum(v for _, v in pts), 2)` uses python's `sum()` (pyfmt::py_sum,
//! Neumaier on 3.12) and `round` (pyfmt::py_round); points string `f"{n}:{v:+g}"` =
//! `pyfmt::py_g(v, 6, true)`; weights/threshold overrides from config `rte_score`
//! (`RteConfig::score_weights` / `score_thresholds`). A weight of 0 never adds a point (`if v:`).
//!
//! Golden: events `score` (in: score_input, weights, thresholds; out: [score, points, call]).

/// score.ScoreInput (field names identical to python).
#[derive(Clone, Debug, PartialEq)]
pub struct ScoreInput {
    pub element: String,
    pub structure: String,
    pub tags: Vec<String>,
    pub tsd_len: Option<i64>,
    pub tsd_verified: bool,
    pub polya_len: f64,
    pub polya_both_sides: bool,
    pub slippage: bool,
    pub beyond_polya_len: i64,
    pub beyond_polya_support: i64,
    pub en_mismatches: Option<i64>,
    pub ends_concordant: bool,
    pub td_source_matches_5p: bool,
    pub element_identity: f64,
    pub inactive_only: bool,
    pub inv_p1: Option<i64>,
    pub junctions_supported: i64,
    pub n_samples: i64,
    pub foldback: bool,
    pub recurrent: bool,
    pub cross_sample_identical: bool,
    pub novel_tier: String,
}

impl Default for ScoreInput {
    fn default() -> Self {
        ScoreInput {
            element: "UNKNOWN".into(),
            structure: "5P_UNRESOLVED".into(),
            tags: Vec::new(),
            tsd_len: None,
            tsd_verified: false,
            polya_len: 0.0,
            polya_both_sides: false,
            slippage: false,
            beyond_polya_len: 0,
            beyond_polya_support: 0,
            en_mismatches: None,
            ends_concordant: false,
            td_source_matches_5p: false,
            element_identity: 0.0,
            inactive_only: false,
            inv_p1: None,
            junctions_supported: 0,
            n_samples: 0,
            foldback: false,
            recurrent: false,
            cross_sample_identical: false,
            novel_tier: String::new(),
        }
    }
}

/// `score(si, weights, thresholds)` -> (score, points string, call). Overrides are applied on
/// top of WEIGHTS / THRESHOLDS like python's `dict.update`. WP-HALL.
pub fn score(_si: &ScoreInput, _weights: &[(String, f64)], _thresholds: &[(String, f64)]) -> (f64, String, String) {
    todo!("WP-HALL: port score.WEIGHTS / THRESHOLDS / score")
}
