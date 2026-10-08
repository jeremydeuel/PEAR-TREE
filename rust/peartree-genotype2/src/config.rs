//! Configuration (`key = value` file, `#` comments; same loader format as the other crates).
//! Unknown keys warn and are ignored (so a legacy `config.genotype.*` file can be passed).

/// Upper bound of a reference-bias value (b > 1 = ALT reads over-captured; allowed, rarely real).
pub const MAX_REF_BIAS: f64 = 4.0;

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
    pub io_fill_bytes: usize,
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
    /// shared alt-read fractions φ at which the per-locus log-likelihood profile `pl_f<‰>` is
    /// reported (the joint step's NOISE hypothesis: one φ shared by every colony)
    pub noise_frac_grid: Vec<f64>,
    /// reference bias `b` (relative capture / assignment efficiency of ALT vs REF reads) applied
    /// to the het alt fraction: φ(1,p) = h·b / (h·b + 1 − h), h = p/2. 1 = no correction (the
    /// pre-correction model, byte-identical output)
    pub ref_bias: f64,
    /// per-locus-kind overrides of `ref_bias` (`ref_bias_kind = TSD:0.95,L1_MED_DELETION:0.6`)
    pub ref_bias_kind: Vec<(String, f64)>,
    /// bias values `b` at which the per-locus het likelihood profile `pl_het_b<‰>` is reported
    /// (empty = no columns); the joint step's `--ref-bias auto|<b>` plugs an estimated `b` in
    pub ref_bias_grid: Vec<f64>,
    pub prior: [f64; 3],
    // ---- extra evidence for undiscovered carriers (extra.rs; default off) ----
    /// collect, at loci this colony shows ALT support for but did not discover, its ALT junction
    /// reads + discordant anchors and their inside mates into `<out>.extra_reads.fa.gz`
    /// (needs `--members`); never changes the genotype rows
    pub gt_extra_reads: bool,
    /// ALT gate: at least this many reads realigned as Alt (LLR >= `llr_informative`)
    pub gt_extra_min_alt: usize,
    /// germline skip: no extra pass at a locus discovered in MORE than this fraction of the
    /// patient's colonies (`#colonies N` of the members table)
    pub gt_extra_max_member_frac: f64,
    /// discordant anchors: last base within this many bp before R / first base after L
    pub gt_extra_disc_span: i64,
    /// discordant anchors: MAPQ floor (the anchor sits in the flank; its mate is in the element)
    pub gt_extra_anchor_mapq: u8,
    /// per locus: at most this many discordant anchors (and inside-mate fetches)
    pub gt_extra_max_mates: usize,
    /// per locus: at most this many ALT junction reads written
    pub gt_extra_max_reads: usize,
    /// per colony: at most this many inside-mate fetches (random access; bounds the extra I/O)
    pub gt_extra_max_mate_fetches: usize,
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
            io_fill_bytes: 256 << 10,
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
            noise_frac_grid: vec![0.01, 0.02, 0.05, 0.1, 0.2],
            ref_bias: 1.0,
            ref_bias_kind: Vec::new(),
            ref_bias_grid: Vec::new(),
            prior: [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            gt_extra_reads: false,
            gt_extra_min_alt: 1,
            gt_extra_max_member_frac: 0.5,
            gt_extra_disc_span: 500,
            gt_extra_anchor_mapq: 20,
            gt_extra_max_mates: 10,
            gt_extra_max_reads: 40,
            gt_extra_max_mate_fetches: 3000,
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

/// `KIND:value,KIND:value` (empty string = no overrides)
fn kind_list(v: &str) -> Result<Vec<(String, f64)>, String> {
    v.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            let (k, x) = s.split_once(':').ok_or_else(|| format!("expected KIND:value, got '{s}'"))?;
            Ok((k.trim().to_string(), num::<f64>(x.trim())?))
        })
        .collect()
}

fn list(v: &str) -> Result<Vec<f64>, String> {
    v.split(',').map(|s| num::<f64>(s.trim())).collect()
}

impl Config {
    pub fn load(path: Option<&str>) -> Result<Config, String> {
        let mut cfg = Config::default();
        if let Some(p) = path {
            cfg.apply_file(p, 0)?;
        }
        if let Ok(v) = std::env::var("PEARTREE_MIN_MAPQ") {
            cfg.min_mapq = v.trim().parse().map_err(|_| format!("PEARTREE_MIN_MAPQ: not a valid u8: '{v}'"))?;
        }
        cfg.validate()?;
        Ok(cfg)
    }

    /// Apply one config file. `include = <path>` (relative to the including file's directory)
    /// applies another file at that point, so a variant config can be "the base + a few keys".
    fn apply_file(&mut self, p: &str, depth: usize) -> Result<(), String> {
        if depth > 8 {
            return Err(format!("{p}: config includes nested too deeply"));
        }
        let text = std::fs::read_to_string(p).map_err(|e| format!("cannot read config {p}: {e}"))?;
        for (i, raw) in text.lines().enumerate() {
            let line = raw.split('#').next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }
            let (key, val) = line.split_once('=').ok_or_else(|| format!("{p}:{}: expected 'key = value'", i + 1))?;
            let (key, val) = (key.trim(), val.trim());
            if key == "include" {
                let base = std::path::Path::new(p).parent().unwrap_or(std::path::Path::new("."));
                let inc = base.join(val);
                self.apply_file(&inc.to_string_lossy(), depth + 1).map_err(|e| format!("{p}:{}: {e}", i + 1))?;
                continue;
            }
            self.set(key, val).map_err(|e| format!("{p}:{}: {e}", i + 1))?;
        }
        Ok(())
    }

    pub fn set(&mut self, key: &str, val: &str) -> Result<(), String> {
        match key {
            "min_mapq" => self.min_mapq = num(val)?,
            "min_mapq_clipped" => self.min_mapq_clipped = num(val)?,
            "min_clip_for_lowmapq" => self.min_clip_for_lowmapq = num(val)?,
            "reads_for_high_coverage" => self.reads_for_high_coverage = num(val)?,
            "lenient_dedup" => self.lenient_dedup = boolean(val)?,
            "io_buffer_bytes" => self.io_buffer_bytes = num(val)?,
            "io_fill_bytes" => self.io_fill_bytes = num(val)?,
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
            "noise_frac_grid" => self.noise_frac_grid = list(val)?,
            "ref_bias" => self.ref_bias = num(val)?,
            "ref_bias_kind" => self.ref_bias_kind = kind_list(val)?,
            "ref_bias_grid" => self.ref_bias_grid = if val.trim().is_empty() { Vec::new() } else { list(val)? },
            "gt_extra_reads" => self.gt_extra_reads = boolean(val)?,
            "gt_extra_min_alt" => self.gt_extra_min_alt = num(val)?,
            "gt_extra_max_member_frac" => self.gt_extra_max_member_frac = num(val)?,
            "gt_extra_disc_span" => self.gt_extra_disc_span = num(val)?,
            "gt_extra_anchor_mapq" => self.gt_extra_anchor_mapq = num(val)?,
            "gt_extra_max_mates" => self.gt_extra_max_mates = num(val)?,
            "gt_extra_max_reads" => self.gt_extra_max_reads = num(val)?,
            "gt_extra_max_mate_fetches" => self.gt_extra_max_mate_fetches = num(val)?,
            "prior" => {
                let l = list(val)?;
                if l.len() != 3 {
                    return Err("prior needs 3 comma-separated values".into());
                }
                self.prior = [l[0], l[1], l[2]];
            }
            other => eprintln!("warning: ignoring unknown config key '{other}'"),
        }
        Ok(())
    }

    /// The reference bias `b` used for loci of `kind` (`LocusKind::as_str`).
    pub fn ref_bias_for(&self, kind: &str) -> f64 {
        self.ref_bias_kind.iter().find(|(k, _)| k == kind).map(|(_, b)| *b).unwrap_or(self.ref_bias)
    }

    fn validate(&self) -> Result<(), String> {
        if self.purity_grid.is_empty() || self.purity_grid.iter().any(|&p| !(p > 0.0 && p <= 1.0)) {
            return Err("purity_grid must be non-empty values in (0, 1]".into());
        }
        if self.noise_frac_grid.iter().any(|&p| !(p > 0.0 && p < 1.0)) {
            return Err("noise_frac_grid values must be in (0, 1)".into());
        }
        let ok_bias = |b: f64| b > 0.0 && b <= MAX_REF_BIAS;
        if !ok_bias(self.ref_bias) || self.ref_bias_kind.iter().any(|(_, b)| !ok_bias(*b)) {
            return Err(format!("ref_bias / ref_bias_kind values must be in (0, {MAX_REF_BIAS}]"));
        }
        if self.ref_bias_grid.iter().any(|&b| !ok_bias(b)) || self.ref_bias_grid.windows(2).any(|w| w[0] >= w[1]) {
            return Err(format!("ref_bias_grid values must be increasing, in (0, {MAX_REF_BIAS}]"));
        }
        const KINDS: [&str; 6] = ["TSD", "BLUNT", "TSD_DELETION", "L1_MED_DELETION", "L1_MED_DUPLICATION", "ONE_SIDED"];
        if let Some((k, _)) = self.ref_bias_kind.iter().find(|(k, _)| !KINDS.contains(&k.as_str())) {
            return Err(format!("ref_bias_kind: unknown locus kind '{k}' (one of {})", KINDS.join(", ")));
        }
        if !(self.clip_prob > 0.0 && self.clip_prob < 1.0) {
            return Err("clip_prob must be in (0, 1)".into());
        }
        if self.gt_extra_min_alt == 0 {
            return Err("gt_extra_min_alt must be >= 1 (the ALT gate)".into());
        }
        if !(self.gt_extra_max_member_frac > 0.0 && self.gt_extra_max_member_frac <= 1.0) {
            return Err("gt_extra_max_member_frac must be in (0, 1] (1 = no germline skip)".into());
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
        c.set("ref_bias", "0.8").unwrap();
        c.set("ref_bias_kind", "L1_MED_DELETION:0.6, TSD:0.95").unwrap();
        c.set("ref_bias_grid", "0.4,0.6,0.8,1.0").unwrap();
        assert!(c.validate().is_ok());
        assert_eq!(c.ref_bias_for("TSD"), 0.95);
        assert_eq!(c.ref_bias_for("L1_MED_DELETION"), 0.6);
        assert_eq!(c.ref_bias_for("BLUNT"), 0.8);
        c.set("ref_bias_kind", "TSDX:0.9").unwrap();
        assert!(c.validate().is_err());
        c.set("ref_bias_kind", "").unwrap();
        c.set("ref_bias_grid", "0.8,0.6").unwrap();
        assert!(c.validate().is_err());
        c.set("ref_bias_grid", "").unwrap();
        assert!(c.validate().is_ok());
        c.set("clip_prob", "1.5").unwrap();
        assert!(c.validate().is_err());
    }

    #[test]
    fn include_applies_the_base_first() {
        let dir = std::env::temp_dir().join(format!("peartree_cfg_inc_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("base.cfg"), "min_mapq = 30\nref_bias = 0.9\n").unwrap();
        std::fs::write(dir.join("variant.cfg"), "include = base.cfg\nref_bias_grid = 0.5,1.0\nref_bias = 0.8\n").unwrap();
        let c = Config::load(Some(&dir.join("variant.cfg").to_string_lossy())).unwrap();
        assert_eq!((c.min_mapq, c.ref_bias), (30, 0.8));
        assert_eq!(c.ref_bias_grid, vec![0.5, 1.0]);
        std::fs::write(dir.join("loop.cfg"), "include = loop.cfg\n").unwrap();
        assert!(Config::load(Some(&dir.join("loop.cfg").to_string_lossy())).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
