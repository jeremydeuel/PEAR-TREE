//! CONFIG['annotate'] keys read by tools/rte, typed, with the Python defaults. FOUNDATION
//! (implemented).
//!
//! Input: the JSON dump of CONFIG['annotate'] that annotate_v2 writes (tools/rte/rust_bridge.py
//! `write_config_json`; callables are resolved or dropped there). Keys tools/rte does not read
//! are ignored. Sub-dicts (`rte_assembly`, `rte_structure`, `rte_transduction`,
//! `rte_pseudogene`, `rte_score`) override the module DEFAULTS key by key, exactly like
//! python's `dict(DEFAULTS); update(cfg)`.
//!
//! Integer-valued keys must hold integral numbers (a float like 600.0 is accepted); python would
//! also accept 600.5 for some of them -- not supported, rejected with an error.

use serde_json::{Map, Value};
use std::path::Path;

fn num(m: &Map<String, Value>, k: &str, d: f64) -> Result<f64, String> {
    match m.get(k) {
        None | Some(Value::Null) => Ok(d),
        Some(Value::Bool(b)) => Ok(if *b { 1.0 } else { 0.0 }),
        Some(v) => v.as_f64().ok_or_else(|| format!("config key {k:?}: expected a number, got {v}")),
    }
}

fn int(m: &Map<String, Value>, k: &str, d: i64) -> Result<i64, String> {
    let x = num(m, k, d as f64)?;
    if x.fract() != 0.0 {
        return Err(format!("config key {k:?}: expected an integer, got {x}"));
    }
    Ok(x as i64)
}

fn opt_str(m: &Map<String, Value>, k: &str) -> Option<String> {
    match m.get(k) {
        Some(Value::String(s)) if !s.is_empty() => Some(s.clone()),
        _ => None,
    }
}

fn sub<'a>(m: &'a Map<String, Value>, k: &str) -> Result<Option<&'a Map<String, Value>>, String> {
    match m.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Object(o)) => Ok(Some(o)),
        Some(v) => Err(format!("config key {k:?}: expected a dict, got {v}")),
    }
}

/// assembly.py DEFAULTS (+ the `.get` defaults of `_terminal_3p`). `rte_assembly` overrides.
#[derive(Clone, Debug, PartialEq)]
pub struct AssemblyCfg {
    pub min_segment_len: i64,
    pub min_element_identity: f64,
    pub min_ref_identity: f64,
    pub min_flank_identity: f64,
    pub min_gap_rescue: i64,
    pub rescue_max_edit_frac: f64,
    pub polya_min_run: i64,
    pub local_window: i64,
    pub min_element_bp: i64,
    pub polya_absorb_max_gap: i64,
    pub polya_absorb_base_frac: f64,
    pub wide_min_bp: i64,
    pub element_local_min_identity: f64,
    pub element_local_margin: f64,
    pub terminal_tolerance: i64,
    pub terminal_min_bp: i64,
    pub terminal_min_identity: f64,
    pub terminal_max_ref_gap: i64,
}

impl AssemblyCfg {
    pub fn from(m: Option<&Map<String, Value>>) -> Result<AssemblyCfg, String> {
        let e = Map::new();
        let m = m.unwrap_or(&e);
        Ok(AssemblyCfg {
            min_segment_len: int(m, "min_segment_len", 20)?,
            min_element_identity: num(m, "min_element_identity", 0.80)?,
            min_ref_identity: num(m, "min_ref_identity", 0.90)?,
            min_flank_identity: num(m, "min_flank_identity", 0.90)?,
            min_gap_rescue: int(m, "min_gap_rescue", 15)?,
            rescue_max_edit_frac: num(m, "rescue_max_edit_frac", 0.15)?,
            polya_min_run: int(m, "polya_min_run", 8)?,
            local_window: int(m, "local_window", 600)?,
            min_element_bp: int(m, "min_element_bp", 30)?,
            polya_absorb_max_gap: int(m, "polya_absorb_max_gap", 20)?,
            polya_absorb_base_frac: num(m, "polya_absorb_base_frac", 0.6)?,
            wide_min_bp: int(m, "wide_min_bp", 30)?,
            element_local_min_identity: num(m, "element_local_min_identity", 0.95)?,
            element_local_margin: num(m, "element_local_margin", 0.03)?,
            terminal_tolerance: int(m, "terminal_tolerance", 5)?,
            terminal_min_bp: int(m, "terminal_min_bp", 20)?,
            terminal_min_identity: num(m, "terminal_min_identity", 0.9)?,
            terminal_max_ref_gap: int(m, "terminal_max_ref_gap", 12)?,
        })
    }
}

/// structure.py DEFAULTS. `rte_structure` overrides (`full_length_tolerance` is replaced whole).
#[derive(Clone, Debug, PartialEq)]
pub struct StructureCfg {
    /// class -> tolerance; `.get(cls_, 50)`
    pub full_length_tolerance: Vec<(String, i64)>,
    pub truncated_3p_tolerance: i64,
    pub en_independent_3p_tolerance: i64,
    pub switch_cluster_bp: i64,
    pub foldback_tolerance: i64,
    pub min_unknown_bp: i64,
    pub templated_max_dist: i64,
    pub premrna_min_bp: i64,
    pub td_min_bp: i64,
    pub td_min_flank_bp: i64,
    pub td_short_flank_identity: f64,
    pub td_max_masked_frac: f64,
    pub td_short_flank_max_at: f64,
    pub td_min_fragments: i64,
    pub templated_min_bp: i64,
    pub templated_min_identity: f64,
    pub templated_min_fragments: i64,
    pub switch_min_seg: i64,
    pub pseudogene_full_length_tol: i64,
    pub terminal_inv_tolerance: i64,
    pub inverted_tail_min: i64,
    pub inverted_tail_max_gap: i64,
}

impl StructureCfg {
    pub fn from(m: Option<&Map<String, Value>>) -> Result<StructureCfg, String> {
        let e = Map::new();
        let m = m.unwrap_or(&e);
        let flt = match m.get("full_length_tolerance") {
            None | Some(Value::Null) => vec![("L1".into(), 60), ("ALU".into(), 20), ("SVA".into(), 60)],
            Some(Value::Object(o)) => {
                let mut v = Vec::new();
                for (k, x) in o {
                    let f = x.as_f64().ok_or_else(|| format!("full_length_tolerance[{k}]: not a number"))?;
                    v.push((k.clone(), f as i64));
                }
                v
            }
            Some(v) => return Err(format!("full_length_tolerance: expected a dict, got {v}")),
        };
        Ok(StructureCfg {
            full_length_tolerance: flt,
            truncated_3p_tolerance: int(m, "truncated_3p_tolerance", 60)?,
            en_independent_3p_tolerance: int(m, "en_independent_3p_tolerance", 20)?,
            switch_cluster_bp: int(m, "switch_cluster_bp", 50)?,
            foldback_tolerance: int(m, "foldback_tolerance", 30)?,
            min_unknown_bp: int(m, "min_unknown_bp", 20)?,
            templated_max_dist: int(m, "templated_max_dist", 250)?,
            premrna_min_bp: int(m, "premrna_min_bp", 30)?,
            td_min_bp: int(m, "td_min_bp", 30)?,
            td_min_flank_bp: int(m, "td_min_flank_bp", 30)?,
            td_short_flank_identity: num(m, "td_short_flank_identity", 0.95)?,
            td_max_masked_frac: num(m, "td_max_masked_frac", 0.5)?,
            td_short_flank_max_at: num(m, "td_short_flank_max_at", 0.5)?,
            td_min_fragments: int(m, "td_min_fragments", 2)?,
            templated_min_bp: int(m, "templated_min_bp", 20)?,
            templated_min_identity: num(m, "templated_min_identity", 0.90)?,
            templated_min_fragments: int(m, "templated_min_fragments", 2)?,
            switch_min_seg: int(m, "switch_min_seg", 20)?,
            pseudogene_full_length_tol: int(m, "pseudogene_full_length_tol", 15)?,
            terminal_inv_tolerance: int(m, "terminal_inv_tolerance", 5)?,
            inverted_tail_min: int(m, "inverted_tail_min", 10)?,
            inverted_tail_max_gap: int(m, "inverted_tail_max_gap", 12)?,
        })
    }

    pub fn full_length_tol(&self, cls: &str) -> i64 {
        self.full_length_tolerance.iter().find(|(k, _)| k == cls).map_or(50, |(_, v)| *v)
    }
}

/// transduction.py DEFAULTS. `rte_transduction` overrides.
#[derive(Clone, Debug, PartialEq)]
pub struct TransductionCfg {
    pub novel_source_max_dist: i64,
    pub novel_source_min_len: i64,
    pub novel_source_min_identity: f64,
    pub novel_source_tier_a: f64,
    pub novel_source_min_mapq: i64,
    pub novel_source_min_seg: i64,
    pub novel_source_min_tag_identity: f64,
}

impl TransductionCfg {
    pub fn from(m: Option<&Map<String, Value>>) -> Result<TransductionCfg, String> {
        let e = Map::new();
        let m = m.unwrap_or(&e);
        Ok(TransductionCfg {
            novel_source_max_dist: int(m, "novel_source_max_dist", 15000)?,
            novel_source_min_len: int(m, "novel_source_min_len", 5500)?,
            novel_source_min_identity: num(m, "novel_source_min_identity", 0.95)?,
            novel_source_tier_a: num(m, "novel_source_tier_a", 0.98)?,
            novel_source_min_mapq: int(m, "novel_source_min_mapq", 20)?,
            novel_source_min_seg: int(m, "novel_source_min_seg", 25)?,
            novel_source_min_tag_identity: num(m, "novel_source_min_tag_identity", 0.95)?,
        })
    }
}

/// pseudogene.py DEFAULTS. `rte_pseudogene` overrides.
#[derive(Clone, Debug, PartialEq)]
pub struct PseudogeneCfg {
    pub exon_junction_overhang: i64,
    pub exon_junction_max_edits: i64,
    pub exon_junction_min_intron: i64,
    pub exon_junction_skip: i64,
}

impl PseudogeneCfg {
    pub fn from(m: Option<&Map<String, Value>>) -> Result<PseudogeneCfg, String> {
        let e = Map::new();
        let m = m.unwrap_or(&e);
        Ok(PseudogeneCfg {
            exon_junction_overhang: int(m, "exon_junction_overhang", 20)?,
            exon_junction_max_edits: int(m, "exon_junction_max_edits", 2)?,
            exon_junction_min_intron: int(m, "exon_junction_min_intron", 30)?,
            exon_junction_skip: int(m, "exon_junction_skip", 1)?,
        })
    }
}

/// annotate_v2 GeneModel windows used by the pre-mRNA lookup (`_candidates` reaches
/// `promoter_up`; `_genic_feature` labels with the splice windows).
#[derive(Clone, Debug, PartialEq)]
pub struct GeneModelCfg {
    pub splice_donor_window: i64,
    pub splice_acceptor_window: i64,
    pub splice_ppt_window: i64,
    pub splice_branch_window: i64,
    pub promoter_up: i64,
}

/// Everything tools/rte reads from CONFIG['annotate'].
#[derive(Clone, Debug, PartialEq)]
pub struct RteConfig {
    /// `rte_library` (REQUIRED; relative = cwd, else the repository root -- library.rs)
    pub rte_library: String,
    pub genome_2bit: Option<String>,
    pub remap_2bit: Option<String>,
    pub remap_index: Option<String>,
    pub remap_rmsk: Option<String>,
    pub exon_annotation: Option<String>,
    pub rte_exon_annotation: Option<String>,
    /// annotate_v2 `gene_model` track (pre-mRNA lookup); None = no gene model
    pub gene_model: Option<String>,
    pub gene_model_cfg: GeneModelCfg,
    /// library.py `young_consensus_regex` (None = the built-in DEFAULT_YOUNG rule)
    pub young_consensus_regex: Option<String>,
    pub gt_compare: bool,
    pub max_reads: usize,
    pub premrna_window: i64,
    pub wide_window: i64,
    pub max_target_site_deletion: i64,
    pub polya_min_fragments: i64,
    pub recurrence_max: i64,
    pub assembly: AssemblyCfg,
    pub structure: StructureCfg,
    pub transduction: TransductionCfg,
    pub pseudogene: PseudogeneCfg,
    /// `rte_score['weights']` overrides (key order as given)
    pub score_weights: Vec<(String, f64)>,
    /// `rte_score['thresholds']` overrides
    pub score_thresholds: Vec<(String, f64)>,
}

fn pairs(m: Option<&Map<String, Value>>, what: &str) -> Result<Vec<(String, f64)>, String> {
    let mut v = Vec::new();
    if let Some(m) = m {
        for (k, x) in m {
            v.push((k.clone(), x.as_f64().ok_or_else(|| format!("rte_score.{what}.{k}: not a number"))?));
        }
    }
    Ok(v)
}

impl RteConfig {
    pub fn from_json(v: &Value) -> Result<RteConfig, String> {
        let m = v.as_object().ok_or("config: expected a JSON object (CONFIG['annotate'])")?;
        let rte_library = opt_str(m, "rte_library").ok_or("config: rte_library is required")?;
        let score = sub(m, "rte_score")?;
        Ok(RteConfig {
            rte_library,
            genome_2bit: opt_str(m, "genome_2bit"),
            remap_2bit: opt_str(m, "remap_2bit"),
            remap_index: opt_str(m, "remap_index"),
            remap_rmsk: opt_str(m, "remap_rmsk"),
            exon_annotation: opt_str(m, "exon_annotation"),
            rte_exon_annotation: opt_str(m, "rte_exon_annotation"),
            gene_model: opt_str(m, "gene_model"),
            gene_model_cfg: GeneModelCfg {
                splice_donor_window: int(m, "splice_donor_window", 6)?,
                splice_acceptor_window: int(m, "splice_acceptor_window", 3)?,
                splice_ppt_window: int(m, "splice_ppt_window", 17)?,
                splice_branch_window: int(m, "splice_branch_window", 45)?,
                promoter_up: int(m, "promoter_up", 2000)?,
            },
            young_consensus_regex: opt_str(m, "young_consensus_regex"),
            gt_compare: match m.get("rte_gt_compare") {
                None | Some(Value::Null) => true,
                Some(Value::Bool(b)) => *b,
                Some(x) => x.as_f64().map(|f| f != 0.0).unwrap_or(true),
            },
            max_reads: int(m, "rte_max_reads", 400)?.max(0) as usize,
            premrna_window: int(m, "rte_premrna_window", 1_000_000)?,
            wide_window: int(m, "rte_wide_window", 10_000)?,
            max_target_site_deletion: int(m, "rte_max_target_site_deletion", 30)?,
            polya_min_fragments: int(m, "rte_polya_min_fragments", 2)?,
            recurrence_max: int(m, "rte_recurrence_max", 3)?,
            assembly: AssemblyCfg::from(sub(m, "rte_assembly")?)?,
            structure: StructureCfg::from(sub(m, "rte_structure")?)?,
            transduction: TransductionCfg::from(sub(m, "rte_transduction")?)?,
            pseudogene: PseudogeneCfg::from(sub(m, "rte_pseudogene")?)?,
            score_weights: pairs(score.and_then(|s| s.get("weights")).and_then(|w| w.as_object()), "weights")?,
            score_thresholds: pairs(score.and_then(|s| s.get("thresholds")).and_then(|w| w.as_object()), "thresholds")?,
        })
    }

    pub fn load(path: &Path) -> Result<RteConfig, String> {
        let txt = std::fs::read_to_string(path).map_err(|e| format!("cannot read config {}: {e}", path.display()))?;
        let v: Value = serde_json::from_str(&txt).map_err(|e| format!("config {}: {e}", path.display()))?;
        RteConfig::from_json(&v)
    }

    /// The exon track used for the exon-exon junction cores: `rte_exon_annotation or
    /// exon_annotation`.
    pub fn exon_track(&self) -> Option<&str> {
        self.rte_exon_annotation.as_deref().or(self.exon_annotation.as_deref())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_and_overrides() {
        let c = RteConfig::from_json(&serde_json::json!({
            "rte_library": "resources/rte_library",
            "genome_2bit": "/x/hg38.2bit",
            "remap_2bit": null,
            "rte_assembly": {"local_window": 300},
            "rte_structure": {"full_length_tolerance": {"L1": 10}},
            "rte_score": {"weights": {"polya_ge10": 3}, "thresholds": {"TPRT": 8.5}},
            "unrelated": [1, 2]
        }))
        .unwrap();
        assert_eq!(c.assembly.local_window, 300);
        assert_eq!(c.assembly.min_segment_len, 20);
        assert_eq!(c.structure.full_length_tol("L1"), 10);
        assert_eq!(c.structure.full_length_tol("ALU"), 50);
        assert_eq!(c.max_reads, 400);
        assert!(c.gt_compare);
        assert_eq!(c.remap_2bit, None);
        assert_eq!(c.score_weights, vec![("polya_ge10".to_string(), 3.0)]);
        assert_eq!(c.score_thresholds, vec![("TPRT".to_string(), 8.5)]);
        assert!(RteConfig::from_json(&serde_json::json!({})).is_err());
        assert!(RteConfig::from_json(&serde_json::json!({"rte_library": "x", "rte_max_reads": 1.5})).is_err());
    }
}
