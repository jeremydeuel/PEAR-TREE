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
    /// TPRT one-sided loci (`contig:L-oneside_L` / `contig:oneside_R-R`): score only the
    /// real junction; the open end carries no consensus and is never scored. Off (default)
    /// -> such a locus has no consensus on its open side and is reported as `error`
    /// (legacy contracts never contain one-sided loci, so their output is unchanged).
    pub one_sided_loci: bool,
    /// Far Bp+Bp pairs (L1-mediated deletion `R < L` / duplication `R > L`, up to 50 kb):
    /// when `|R - L|` exceeds this, the two breakpoints are depth-gated and fetched as two
    /// separate 1-bp windows instead of one `[min, max]` window (which would count every read
    /// of the deleted/duplicated span as depth -> `high-coverage`, and let a non-spanning
    /// mate claim the fragment's qname before the spanning one). 0 = off (legacy).
    pub split_breakpoint_span: i64,
    /// One-sided loci only: a read whose alignment ends, on the OPEN side, in a soft clip of
    /// at least `one_sided_open_min_clip` bases within this many bp of the real breakpoint is a
    /// junction read of the missing end (its aligned part runs over the real breakpoint in
    /// reference configuration) and is skipped. Two-sided genotyping scores such a read alt
    /// on its own side; with the open side unscored it would be a false ref vote, pulling a
    /// het towards VAF 1/3. 0 = off. Spans the TSD (<= 40) and target-site deletion (<= 30).
    pub one_sided_open_window: i64,
    pub one_sided_open_min_clip: i64,
    /// Far L1-mediated duplications (`R - L` >= this, i.e. longer than a read): the
    /// alt haplotype `ref[..R) + element + ref[L..)` still carries BOTH reference junctions
    /// (the duplicated copy), so every alt-haplotype molecule also yields about one
    /// reference-configuration read; a het reads VAF ~1/3, a hom ~1/2. Discount one ref read
    /// per alt read (`n_ref -= min(n_ref, n_alt)`, ref score scaled alike) before the VAF
    /// bands. Wild-type colonies (n_alt = 0) are unchanged. 0 = off.
    pub dup_ref_discount_min_span: i64,
    /// Single-junction evidence (one-sided loci, split far pairs): count each reference vote
    /// as HALF (`n_ref -> ceil(n_ref / 2)`, ref score halved). The VAF bands are calibrated on
    /// a TSD locus, where alt reads come from TWO junctions but one reference read spans both
    /// (het VAF = 2f/(2f+1), f = junction-read yield). A one-sided locus has alt from one
    /// junction per reference span, a far pair has two disjoint reference spans: both read
    /// f/(1+f) (E2E: ~0.30-0.38 for true hets) unless the reference side is halved. Applied
    /// before the duplication discount. false = off.
    pub halve_single_junction_ref: bool,
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
            one_sided_loci: false,
            split_breakpoint_span: 0,
            one_sided_open_window: 0,
            one_sided_open_min_clip: 5,
            dup_ref_discount_min_span: 0,
            halve_single_junction_ref: false,
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
            "one_sided_loci" => self.one_sided_loci = parse_bool(val)?,
            "split_breakpoint_span" => self.split_breakpoint_span = parse_num(val)?,
            "one_sided_open_window" => self.one_sided_open_window = parse_num(val)?,
            "one_sided_open_min_clip" => self.one_sided_open_min_clip = parse_num(val)?,
            "dup_ref_discount_min_span" => self.dup_ref_discount_min_span = parse_num(val)?,
            "halve_single_junction_ref" => self.halve_single_junction_ref = parse_bool(val)?,
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
        // TPRT keys default off (legacy byte-identical output)
        assert!(!c.one_sided_loci);
        assert_eq!(c.split_breakpoint_span, 0);
        assert_eq!(c.one_sided_open_window, 0);
        assert_eq!(c.dup_ref_discount_min_span, 0);
        assert!(!c.halve_single_junction_ref);
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
        c.set("one_sided_loci", "true").unwrap();
        c.set("split_breakpoint_span", "50").unwrap();
        assert!(c.one_sided_loci);
        assert_eq!(c.split_breakpoint_span, 50);
    }
}
