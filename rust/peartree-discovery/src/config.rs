//! Discovery-relevant configuration, mirroring the `discovery` and `adapters`
//! sections of src/config.py. (Stage 2 will make these load from a file.)

use rustc_hash::FxHashSet;

pub const MIN_MAPQ: u8 = 40;
pub const MIN_CLIP_LEN: usize = 12;
pub const MIN_EVIDENCE_READS_PER_BREAKPOINT: usize = 2;
pub const MIN_ADAPTERLEN_FOR_CLIP: usize = 4;
pub const MIN_GOOD_BASES: usize = 10;
pub const EXCLUDE_SAME_CONTIG_SUPPLEMENTARY: i64 = 1000;

pub const POLYA_CUTOFF: usize = 12;

// Clustering + TSD-pairing windows used in cleanup()/output(). Extracted verbatim
// from discovery.py (OBS-2); values unchanged, so output stays byte-identical.
// OBS-3 will revisit CLUSTER_WINDOW vs TSD_MAX jointly with SENS-1.
//
// single-linkage cluster window: two breakpoints closer than this join one group
pub const CLUSTER_WINDOW: i64 = 6;
// TSD size bounds when pairing a left with a right breakpoint (right - left):
// below TSD_MIN it is not a TSD; above TSD_MAX the gap is too large to pair.
pub const TSD_MIN: i64 = 2;
pub const TSD_MAX: i64 = 40;
// polyA-rescue proximity to the anchor breakpoint: a polyA within POLYA_NEAR_DIST
// is too close (skipped); it rescues only when strictly inside POLYA_FAR_DIST.
pub const POLYA_NEAR_DIST: i64 = 12;
pub const POLYA_FAR_DIST: i64 = 120;

// drop clipped reads whose XA/SA shows the whole read maps contiguously
// elsewhere (not a real junction). On by default; set PEARTREE_KEEP_FULLMAP=1
// to disable (e.g. to reproduce pre-filter output).
pub fn reject_fully_mapping_reads() -> bool {
    !matches!(std::env::var("PEARTREE_KEEP_FULLMAP").as_deref(), Ok("1") | Ok("true") | Ok("True"))
}

// clip side constants (match the Python ints)
pub const CLIP_RIGHT: i32 = 1;
pub const CLIP_LEFT: i32 = 2;

pub const ADAPTERS: [&[u8]; 4] = [
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    b"AGATCGGAAAGCACACGTCTGAACTCCAGTCA",
    b"AGATCGGAAAGCGTCGTGTAGGGAAAGAGTGT",
];

/// Runtime discovery configuration (OBS-4). Mirrors the tunable parameters of the
/// `discovery` block in src/config.py. Every field defaults to the module constant
/// above, so a run with no `--config` file and no env override is byte-identical to
/// the pre-OBS-4 port.
///
/// ⚠ The default `min_mapq` is 40, matching the *generic* src/config.py. Production
/// configs `config_hs.py` / `config_mm.py` use **60**. Pin `min_mapq` (via `--config`
/// or `PEARTREE_MIN_MAPQ`) on any real-WGS run, or the port is only faithful to the
/// generic config.
#[derive(Clone, Debug)]
pub struct DiscoveryConfig {
    pub min_mapq: u8,
    pub min_evidence_reads_per_breakpoint: usize,
    pub min_good_bases: usize,
    pub exclude_same_contig_supplementary: i64,
    pub cluster_window: i64,
    pub tsd_min: i64,
    pub tsd_max: i64,
    pub polya_near_dist: i64,
    pub polya_far_dist: i64,
    pub reject_fully_mapping_reads: bool,
    /// SPEC-5/SENS-4: explicit primary-assembly allowlist. When `Some`, a contig is
    /// processed iff its name is in the set (replacing the `len(name) > 5` + MT/chrM
    /// heuristic). When `None`, the legacy heuristic applies, so output is unchanged.
    pub contig_allowlist: Option<FxHashSet<String>>,
    /// SPEC-5: path to a BED file of regions to drop breakpoints in. `None` = no-op.
    pub exclude_bed: Option<String>,
}

impl Default for DiscoveryConfig {
    fn default() -> Self {
        DiscoveryConfig {
            min_mapq: MIN_MAPQ,
            min_evidence_reads_per_breakpoint: MIN_EVIDENCE_READS_PER_BREAKPOINT,
            min_good_bases: MIN_GOOD_BASES,
            exclude_same_contig_supplementary: EXCLUDE_SAME_CONTIG_SUPPLEMENTARY,
            cluster_window: CLUSTER_WINDOW,
            tsd_min: TSD_MIN,
            tsd_max: TSD_MAX,
            polya_near_dist: POLYA_NEAR_DIST,
            polya_far_dist: POLYA_FAR_DIST,
            reject_fully_mapping_reads: reject_fully_mapping_reads(),
            contig_allowlist: None,
            exclude_bed: None,
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

impl DiscoveryConfig {
    /// Build config: start from defaults, overlay a `key = value` file (if given),
    /// then overlay environment overrides (which always win). `#` starts a comment.
    pub fn load(path: Option<&str>) -> Result<DiscoveryConfig, String> {
        let mut cfg = DiscoveryConfig::default();
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
        cfg.apply_env()?;
        Ok(cfg)
    }

    fn set(&mut self, key: &str, val: &str) -> Result<(), String> {
        match key {
            "min_mapq" => self.min_mapq = parse_num(val)?,
            "min_evidence_reads_per_breakpoint" => self.min_evidence_reads_per_breakpoint = parse_num(val)?,
            "min_good_bases" => self.min_good_bases = parse_num(val)?,
            "exclude_same_contig_supplementary" => self.exclude_same_contig_supplementary = parse_num(val)?,
            "cluster_window" => self.cluster_window = parse_num(val)?,
            "tsd_min" => self.tsd_min = parse_num(val)?,
            // accept the Python name `max_bp_window` as an alias for tsd_max
            "tsd_max" | "max_bp_window" => self.tsd_max = parse_num(val)?,
            "polya_near_dist" => self.polya_near_dist = parse_num(val)?,
            "polya_far_dist" => self.polya_far_dist = parse_num(val)?,
            "reject_fully_mapping_reads" => self.reject_fully_mapping_reads = parse_bool(val)?,
            // comma-separated inline allowlist, e.g. "1,2,...,X,Y,MT"
            "contig_allowlist" => {
                self.contig_allowlist = Some(val.split(',').map(|s| s.trim()).filter(|s| !s.is_empty()).map(String::from).collect())
            }
            // allowlist from a file, one contig name per line (# comments allowed)
            "contig_allowlist_file" => {
                let text = std::fs::read_to_string(val).map_err(|e| format!("cannot read contig_allowlist_file {val}: {e}"))?;
                self.contig_allowlist = Some(
                    text.lines()
                        .map(|l| l.split('#').next().unwrap_or("").trim())
                        .filter(|s| !s.is_empty())
                        .map(String::from)
                        .collect(),
                )
            }
            "exclude_bed" => self.exclude_bed = Some(val.to_string()),
            other => eprintln!("warning: ignoring unknown config key '{other}'"),
        }
        Ok(())
    }

    fn apply_env(&mut self) -> Result<(), String> {
        if let Ok(v) = std::env::var("PEARTREE_MIN_MAPQ") {
            self.min_mapq = v.trim().parse().map_err(|_| format!("PEARTREE_MIN_MAPQ: not a valid u8: '{v}'"))?;
        }
        // PEARTREE_KEEP_FULLMAP=1 keeps full-mapping reads (env wins over the file).
        if matches!(std::env::var("PEARTREE_KEEP_FULLMAP").as_deref(), Ok("1") | Ok("true") | Ok("True")) {
            self.reject_fully_mapping_reads = false;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_constants() {
        let c = DiscoveryConfig::default();
        assert_eq!(c.min_mapq, MIN_MAPQ);
        assert_eq!(c.min_evidence_reads_per_breakpoint, MIN_EVIDENCE_READS_PER_BREAKPOINT);
        assert_eq!(c.cluster_window, CLUSTER_WINDOW);
        assert_eq!(c.tsd_max, TSD_MAX);
        assert_eq!(c.polya_far_dist, POLYA_FAR_DIST);
    }

    #[test]
    fn set_overrides_and_alias() {
        let mut c = DiscoveryConfig::default();
        c.set("min_mapq", "60").unwrap();
        c.set("max_bp_window", "50").unwrap(); // Python-name alias for tsd_max
        c.set("reject_fully_mapping_reads", "false").unwrap();
        assert_eq!(c.min_mapq, 60);
        assert_eq!(c.tsd_max, 50);
        assert!(!c.reject_fully_mapping_reads);
    }

    #[test]
    fn unknown_key_is_ignored_not_errored() {
        let mut c = DiscoveryConfig::default();
        assert!(c.set("no_such_key", "3").is_ok());
    }

    #[test]
    fn bad_value_errors() {
        let mut c = DiscoveryConfig::default();
        assert!(c.set("min_mapq", "notnum").is_err());
        assert!(c.set("reject_fully_mapping_reads", "maybe").is_err());
    }

    #[test]
    fn inline_allowlist_parses() {
        let mut c = DiscoveryConfig::default();
        assert!(c.contig_allowlist.is_none());
        c.set("contig_allowlist", "1, 2 ,X,, NC_000014.9").unwrap();
        let set = c.contig_allowlist.unwrap();
        assert_eq!(set.len(), 4); // empty entry between the commas is dropped
        assert!(set.contains("1"));
        assert!(set.contains("X"));
        assert!(set.contains("NC_000014.9"));
        assert!(!set.contains("3"));
    }
}
