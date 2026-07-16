//! Genotyping configuration, mirroring the `genotyping` block of src/config.py.
//!
//! Every field defaults to the *generic* src/config.py value, so a run with no
//! `--config` file is byte-identical to `python src/main.py --step genotype` under
//! the generic config. Production configs (config_hs.py / config_mm.py) raise
//! `min_mapq` to 60 — pin it via `--config` on real-WGS runs.

/// Genotype call vocabulary. These strings are the on-disk contract consumed by
/// combine_genotypes.py; keep them byte-identical to src/genotyping_insertion.py.
pub const GT_ARTEFACT: &str = "artefact";
pub const GT_WILDTYPE: &str = "wild-type";
pub const GT_HETEROZYGOUS: &str = "heterozygous";
pub const GT_HOMOZYGOUS: &str = "homozygous";
pub const GT_INSERTION: &str = "insertion";
pub const GT_INSERTION_UNCERTAIN: &str = "insertion?";
pub const GT_WILDTYPE_UNCERTAIN: &str = "wild-type?";
pub const GT_NO_COVERAGE: &str = "no-coverage";
pub const GT_HIGH_COVERAGE: &str = "high-coverage";
pub const GT_ERROR: &str = "error";

#[derive(Clone, Debug)]
pub struct GenotypingConfig {
    pub min_mapq: u8,
    pub min_score_for_call: i64,
    pub min_supporting_reads: i64,
    pub min_reads_for_zygosity: i64,
    pub recover_low_coverage_presence: bool,
    pub reads_for_high_coverage: i64,
    pub art_min_score: i64,
    pub vaf_wildtype_max: f64,
    pub vaf_het_min: f64,
    pub vaf_hom_min: f64,
    pub artefact_read_fraction: f64,
    pub min_artefact_reads: i64,
    pub double_alt_is_artefact: bool,
}

impl Default for GenotypingConfig {
    fn default() -> Self {
        // Generic src/config.py['genotyping'] defaults.
        GenotypingConfig {
            min_mapq: 40,
            min_score_for_call: 6,
            min_supporting_reads: 2,
            min_reads_for_zygosity: 6,
            recover_low_coverage_presence: true,
            reads_for_high_coverage: 180,
            art_min_score: 60,
            vaf_wildtype_max: 0.10,
            vaf_het_min: 0.30,
            vaf_hom_min: 0.85,
            artefact_read_fraction: 0.5,
            min_artefact_reads: 2,
            double_alt_is_artefact: true,
        }
    }
}

fn parse_num<T: std::str::FromStr>(v: &str) -> Result<T, String> {
    v.parse().map_err(|_| format!("not a valid number: '{v}'"))
}

fn parse_bool(v: &str) -> Result<bool, String> {
    match v {
        "true" | "True" | "1" => Ok(true),
        "false" | "False" | "0" => Ok(false),
        _ => Err(format!("not a valid bool: '{v}'")),
    }
}

impl GenotypingConfig {
    /// Start from defaults, overlay a `key = value` file (if given), then env
    /// overrides (which win). `#` starts a comment. Same format as the discovery
    /// crate's loader, so a single config file can carry both sections.
    pub fn load(path: Option<&str>) -> Result<GenotypingConfig, String> {
        let mut cfg = GenotypingConfig::default();
        if let Some(p) = path {
            let text = std::fs::read_to_string(p).map_err(|e| format!("cannot read config {p}: {e}"))?;
            for (i, raw) in text.lines().enumerate() {
                let line = raw.split('#').next().unwrap_or("").trim();
                if line.is_empty() {
                    continue;
                }
                let (key, val) = line
                    .split_once('=')
                    .ok_or_else(|| format!("{p}:{}: expected 'key = value'", i + 1))?;
                cfg.set(key.trim(), val.trim()).map_err(|e| format!("{p}:{}: {e}", i + 1))?;
            }
        }
        if let Ok(v) = std::env::var("PEARTREE_MIN_MAPQ") {
            cfg.min_mapq = v.trim().parse().map_err(|_| format!("PEARTREE_MIN_MAPQ: not a valid u8: '{v}'"))?;
        }
        Ok(cfg)
    }

    fn set(&mut self, key: &str, val: &str) -> Result<(), String> {
        match key {
            "min_mapq" => self.min_mapq = parse_num(val)?,
            "min_score_for_call" => self.min_score_for_call = parse_num(val)?,
            "min_supporting_reads" => self.min_supporting_reads = parse_num(val)?,
            "min_reads_for_zygosity" => self.min_reads_for_zygosity = parse_num(val)?,
            "recover_low_coverage_presence" => self.recover_low_coverage_presence = parse_bool(val)?,
            "reads_for_high_coverage" => self.reads_for_high_coverage = parse_num(val)?,
            "art_min_score" => self.art_min_score = parse_num(val)?,
            "vaf_wildtype_max" => self.vaf_wildtype_max = parse_num(val)?,
            "vaf_het_min" => self.vaf_het_min = parse_num(val)?,
            "vaf_hom_min" => self.vaf_hom_min = parse_num(val)?,
            "artefact_read_fraction" => self.artefact_read_fraction = parse_num(val)?,
            "min_artefact_reads" => self.min_artefact_reads = parse_num(val)?,
            "double_alt_is_artefact" => self.double_alt_is_artefact = parse_bool(val)?,
            other => eprintln!("warning: ignoring unknown config key '{other}'"),
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_generic_config() {
        let c = GenotypingConfig::default();
        assert_eq!(c.min_mapq, 40);
        assert_eq!(c.min_supporting_reads, 2);
        assert_eq!(c.min_reads_for_zygosity, 6);
        assert!(c.recover_low_coverage_presence);
        assert_eq!(c.reads_for_high_coverage, 180);
        assert!((c.vaf_hom_min - 0.85).abs() < 1e-12);
    }

    #[test]
    fn set_overrides() {
        let mut c = GenotypingConfig::default();
        c.set("min_mapq", "60").unwrap();
        c.set("reads_for_high_coverage", "250").unwrap();
        c.set("recover_low_coverage_presence", "false").unwrap();
        assert_eq!(c.min_mapq, 60);
        assert_eq!(c.reads_for_high_coverage, 250);
        assert!(!c.recover_low_coverage_presence);
    }
}
