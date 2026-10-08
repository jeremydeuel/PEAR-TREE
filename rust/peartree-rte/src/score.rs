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

use crate::pyfmt::{py_g, py_round, py_sum};

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

/// score.WEIGHTS (python dict order).
pub const WEIGHTS: &[(&str, f64)] = &[
    ("tsd_4_25", 3.0),
    ("tsd_other_1_50", 1.0),
    ("tsd_deletion_le20", 1.0),
    ("polya_ge10", 2.0),
    ("polya_5_9", 0.5),
    ("beyond_polya", 2.0),
    ("beyond_polya_2frag", 1.5),
    ("en_0_1", 2.0),
    ("en_2", 1.0),
    ("en_3", 0.5),
    ("ends_concordant", 1.0),
    ("td_source_matches_5p", 2.0),
    ("active_identity_ge98", 1.0),
    ("twin_priming_ge590", 1.0),
    ("twin_priming_lt590", -1.5),
    ("exon_junction", 3.0),
    ("novel_source", 1.0),
    ("novel_source_tier_b", 0.5),
    ("junction_supported", 0.5),
    ("multi_colony", 1.0),
    ("polya_both_sides", -4.0),
    ("slippage_context", -4.0),
    ("chimeric_ends", -6.0),
    ("tsd_gt50", -6.0),
    ("inactive_only", -2.0),
    ("foldback_clip", -3.0),
    ("recurrence", -2.0),
    ("cross_sample_identical", -2.0),
    ("no_element_no_tail", -1.0),
    ("en_independent", -2.0),
];

/// score.THRESHOLDS
pub const THRESHOLDS: &[(&str, f64)] = &[("TPRT", 7.0), ("LIKELY_TPRT", 4.0), ("UNCERTAIN", 0.0)];

fn lookup(defaults: &[(&str, f64)], over: &[(String, f64)], name: &str) -> f64 {
    // dict.update: the last override of a key wins
    if let Some((_, v)) = over.iter().rev().find(|(k, _)| k == name) {
        return *v;
    }
    defaults.iter().find(|(k, _)| *k == name).map(|(_, v)| *v).expect("unknown score weight")
}

/// `score(si, weights, thresholds)` -> (score, points string, call). Overrides are applied on
/// top of WEIGHTS / THRESHOLDS like python's `dict.update`.
pub fn score(si: &ScoreInput, weights: &[(String, f64)], thresholds: &[(String, f64)]) -> (f64, String, String) {
    let mut pts: Vec<(&'static str, f64)> = Vec::new();
    let mut add = |name: &'static str, n: Option<f64>| {
        let v = n.unwrap_or_else(|| lookup(WEIGHTS, weights, name));
        if v != 0.0 {
            pts.push((name, v));
        }
    };
    let has = |t: &str| si.tags.iter().any(|x| x == t);

    let t = si.tsd_len;
    let l1dup = has("L1_MED_DUPLICATION")
        && si.en_mismatches.is_some_and(|m| m <= 2)
        && !si.slippage
        && t.is_some_and(|t| t > 150 || si.n_samples >= 2);
    if let Some(t) = t {
        if t > 50 && !l1dup {
            add("tsd_gt50", None);
        } else if t > 50 {
            // L1-mediated duplication: no penalty, no points
        } else if (4..=25).contains(&t) && si.tsd_verified {
            add("tsd_4_25", None);
        } else if t > 0 && si.tsd_verified {
            add("tsd_other_1_50", None);
        } else if (-20..0).contains(&t) {
            add("tsd_deletion_le20", None);
        }
    }
    if si.polya_len >= 10.0 {
        add("polya_ge10", None);
    } else if si.polya_len >= 5.0 {
        add("polya_5_9", None);
    }
    if si.beyond_polya_len >= 10 && si.polya_len >= 5.0 {
        add("beyond_polya", None);
        if si.beyond_polya_support >= 2 {
            add("beyond_polya_2frag", None);
        }
    }
    if let Some(mm) = si.en_mismatches {
        if !si.slippage {
            if mm <= 1 {
                add("en_0_1", None);
            } else if mm == 2 {
                add("en_2", None);
            } else if mm == 3 {
                add("en_3", None);
            }
        }
    }
    if si.td_source_matches_5p {
        add("td_source_matches_5p", None);
    } else if si.ends_concordant {
        add("ends_concordant", None);
    }
    if si.element_identity >= 0.98 {
        add("active_identity_ge98", None);
    }
    if si.element == "L1" && si.structure.starts_with("INVERTED_5P") {
        if let Some(p1) = si.inv_p1 {
            add(if p1 >= 590 { "twin_priming_ge590" } else { "twin_priming_lt590" }, None);
        }
    }
    if has("EXON_JUNCTION") {
        add("exon_junction", None);
    }
    if has("EN_INDEPENDENT") {
        add("en_independent", None);
    }
    if has("NOVEL_SOURCE") {
        add(if si.novel_tier == "B" { "novel_source_tier_b" } else { "novel_source" }, None);
    }
    if si.junctions_supported != 0 {
        let w = lookup(WEIGHTS, weights, "junction_supported");
        add("junction_supported", Some(w * si.junctions_supported.min(3) as f64));
    }
    if si.n_samples >= 2 {
        add("multi_colony", None);
    }
    // artefacts
    if si.polya_both_sides {
        add("polya_both_sides", None);
    }
    if si.slippage {
        add("slippage_context", None);
    }
    if has("CHIMERIC_ENDS") {
        add("chimeric_ends", None);
    }
    if si.inactive_only {
        add("inactive_only", None);
    }
    if si.foldback {
        add("foldback_clip", None);
    }
    if si.recurrent {
        add("recurrence", None);
    }
    if si.cross_sample_identical && si.n_samples <= 1 {
        add("cross_sample_identical", None);
    }
    if (si.element == "UNKNOWN" || si.element == "NON_TPRT") && si.polya_len < 5.0 {
        add("no_element_no_tail", None);
    }
    let total = py_round(py_sum(pts.iter().map(|(_, v)| *v)), 2);
    let th = |name: &str| lookup(THRESHOLDS, thresholds, name);
    let mut call = if total >= th("TPRT") {
        "TPRT"
    } else if total >= th("LIKELY_TPRT") {
        "LIKELY_TPRT"
    } else if total >= th("UNCERTAIN") {
        "UNCERTAIN"
    } else {
        "ARTEFACT_LIKE"
    };
    // hard artefact signatures cap the call regardless of the other points
    if (call == "TPRT" || call == "LIKELY_TPRT")
        && pts.iter().any(|(n, _)| matches!(*n, "chimeric_ends" | "tsd_gt50" | "polya_both_sides"))
    {
        call = "UNCERTAIN";
    }
    let points = if pts.is_empty() {
        ".".to_string()
    } else {
        pts.iter().map(|(n, v)| format!("{}:{}", n, py_g(*v, 6, true))).collect::<Vec<_>>().join(";")
    };
    (total, points, call.to_string())
}
