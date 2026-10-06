//! Configuration (`key = value` file, `#` comments; same loader format as the other crates).
//! Unknown keys warn and are ignored (so a legacy `config.genotype.*` file can be passed).

#[derive(Clone, Debug)]
pub struct Config {
    // ---- read gates (driver) ----
    pub min_mapq: u8,
    /// a read below `min_mapq` is still used when its MAPQ >= this AND it carries a soft clip
    /// >= `min_clip_for_lowmapq` at a junction-facing end (its MAPQ reflects only the flank)
    pub min_mapq_clipped: u8,
    pub min_clip_for_lowmapq: usize,
    pub reads_for_high_coverage: i64,
    pub lenient_dedup: bool,
    pub io_buffer_bytes: usize,
    pub heartbeat_every: usize,
    // ---- discordant anchors ----
    pub disc_max_tlen: i64,
    pub disc_span: i64,
    /// log-likelihood contribution (nats, towards alt) of one discordant anchor; 0 = count only
    pub disc_weight_nats: f64,
    // ---- haplotype construction ----
    pub flank: i64,
    pub merge_overlap_min: usize,
    pub merge_overlap_max_mismatch_frac: f64,
    pub dup_retained_min_span: i64,
    // ---- aligner ----
    pub base_q_min: u8,
    pub base_q_max: u8,
    pub gap_open_phred: f64,
    pub gap_ext_phred: f64,
    pub homopolymer_min_len: usize,
    pub homopolymer_gap_open_phred: f64,
    pub clip_prob: f64,
    pub band_halfwidth: usize,
    pub realign_fallback_full: bool,
    pub fallback_slack_nats: f64,
    pub llr_informative: f64,
    pub min_explained_frac: f64,
    // ---- genotype model ----
    pub bg_alt_rate: f64,
    pub purity_grid: Vec<f64>,
    pub prior: [f64; 3],
    pub p_confident: f64,
    pub p_present_certain: f64,
    pub p_present_uncertain: f64,
    pub min_artefact_reads: i64,
    pub artefact_read_fraction: f64,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            min_mapq: 60,
            min_mapq_clipped: 20,
            min_clip_for_lowmapq: 20,
            reads_for_high_coverage: 180,
            lenient_dedup: false,
            io_buffer_bytes: 4 << 20,
            heartbeat_every: 1000,
            disc_max_tlen: 1000,
            disc_span: 300,
            disc_weight_nats: 0.0,
            flank: 300,
            merge_overlap_min: 20,
            merge_overlap_max_mismatch_frac: 0.1,
            dup_retained_min_span: 150,
            base_q_min: 2,
            base_q_max: 40,
            gap_open_phred: 30.0,
            gap_ext_phred: 10.0,
            homopolymer_min_len: 4,
            homopolymer_gap_open_phred: 13.0,
            clip_prob: 0.25,
            band_halfwidth: 32,
            realign_fallback_full: true,
            fallback_slack_nats: 20.0,
            llr_informative: 4.6,
            min_explained_frac: 0.8,
            bg_alt_rate: 0.005,
            purity_grid: vec![1.0, 0.9, 0.8, 0.7],
            prior: [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            p_confident: 0.9,
            p_present_certain: 0.9,
            p_present_uncertain: 0.5,
            min_artefact_reads: 2,
            artefact_read_fraction: 0.5,
        }
    }
}

fn num<T: std::str::FromStr>(v: &str) -> Result<T, String> {
    v.parse().map_err(|_| format!("not a valid number: '{v}'"))
}

fn boolean(v: &str) -> Result<bool, String> {
    match v {
        "true" | "True" | "1" => Ok(true),
        "false" | "False" | "0" => Ok(false),
        _ => Err(format!("not a valid bool: '{v}'")),
    }
}

fn list(v: &str) -> Result<Vec<f64>, String> {
    v.split(',').map(|s| num::<f64>(s.trim())).collect()
}

impl Config {
    pub fn load(path: Option<&str>) -> Result<Config, String> {
        let mut cfg = Config::default();
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
        cfg.validate()?;
        Ok(cfg)
    }

    pub fn set(&mut self, key: &str, val: &str) -> Result<(), String> {
        match key {
            "min_mapq" => self.min_mapq = num(val)?,
            "min_mapq_clipped" => self.min_mapq_clipped = num(val)?,
            "min_clip_for_lowmapq" => self.min_clip_for_lowmapq = num(val)?,
            "reads_for_high_coverage" => self.reads_for_high_coverage = num(val)?,
            "lenient_dedup" => self.lenient_dedup = boolean(val)?,
            "io_buffer_bytes" => self.io_buffer_bytes = num(val)?,
            "heartbeat_every" => self.heartbeat_every = num(val)?,
            "disc_max_tlen" => self.disc_max_tlen = num(val)?,
            "disc_span" => self.disc_span = num(val)?,
            "disc_weight_nats" => self.disc_weight_nats = num(val)?,
            "flank" => self.flank = num(val)?,
            "merge_overlap_min" => self.merge_overlap_min = num(val)?,
            "merge_overlap_max_mismatch_frac" => self.merge_overlap_max_mismatch_frac = num(val)?,
            "dup_retained_min_span" => self.dup_retained_min_span = num(val)?,
            "base_q_min" => self.base_q_min = num(val)?,
            "base_q_max" => self.base_q_max = num(val)?,
            "gap_open_phred" => self.gap_open_phred = num(val)?,
            "gap_ext_phred" => self.gap_ext_phred = num(val)?,
            "homopolymer_min_len" => self.homopolymer_min_len = num(val)?,
            "homopolymer_gap_open_phred" => self.homopolymer_gap_open_phred = num(val)?,
            "clip_prob" => self.clip_prob = num(val)?,
            "band_halfwidth" => self.band_halfwidth = num(val)?,
            "realign_fallback_full" => self.realign_fallback_full = boolean(val)?,
            "fallback_slack_nats" => self.fallback_slack_nats = num(val)?,
            "llr_informative" => self.llr_informative = num(val)?,
            "min_explained_frac" => self.min_explained_frac = num(val)?,
            "bg_alt_rate" => self.bg_alt_rate = num(val)?,
            "purity_grid" => self.purity_grid = list(val)?,
            "prior" => {
                let l = list(val)?;
                if l.len() != 3 {
                    return Err("prior needs 3 comma-separated values".into());
                }
                self.prior = [l[0], l[1], l[2]];
            }
            "p_confident" => self.p_confident = num(val)?,
            "p_present_certain" => self.p_present_certain = num(val)?,
            "p_present_uncertain" => self.p_present_uncertain = num(val)?,
            "min_artefact_reads" => self.min_artefact_reads = num(val)?,
            "artefact_read_fraction" => self.artefact_read_fraction = num(val)?,
            other => eprintln!("warning: ignoring unknown config key '{other}'"),
        }
        Ok(())
    }

    fn validate(&self) -> Result<(), String> {
        if self.purity_grid.is_empty() || self.purity_grid.iter().any(|&p| !(p > 0.0 && p <= 1.0)) {
            return Err("purity_grid must be non-empty values in (0, 1]".into());
        }
        if !(self.clip_prob > 0.0 && self.clip_prob < 1.0) {
            return Err("clip_prob must be in (0, 1)".into());
        }
        if self.flank < 50 {
            return Err("flank must be >= 50".into());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_and_overrides() {
        let mut c = Config::default();
        assert_eq!(c.min_mapq, 60);
        c.set("min_mapq", "20").unwrap();
        c.set("purity_grid", "1.0, 0.8").unwrap();
        c.set("prior", "0.98,0.015,0.005").unwrap();
        assert_eq!(c.min_mapq, 20);
        assert_eq!(c.purity_grid, vec![1.0, 0.8]);
        assert!((c.prior[1] - 0.015).abs() < 1e-12);
        assert!(c.validate().is_ok());
        c.set("clip_prob", "1.5").unwrap();
        assert!(c.validate().is_err());
    }
}
