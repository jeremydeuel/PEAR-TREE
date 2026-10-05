//! Configuration. OWNER: P1.
//!
//! The Python step reads `CONFIG` from `src/config.py` (a copy of a `cluster/config.py.*`).
//! The Rust binary takes `--config <path>`:
//!   * `*.json`: the JSON dump of the whole `CONFIG` dict;
//!   * anything else (e.g. `src/config.py`): executed by a python subprocess
//!     (`$PEARTREE_PYTHON`, else `python3`) that prints
//!     `json.dumps(runpy.run_path(path)["CONFIG"], default=lambda o: None)` (lambdas in the
//!     `annotate` section become null). SPEC.md §1.2 "Config".
//!
//! Every key the Python step consumes, with the Python default (`cfg.get(key, default)`), is a
//! field below. Keys the Python reads with `CONFIG[...]` (no default) are REQUIRED here too.
//! Truthiness follows Python: `bool(value)` for flags, `int(v or 0)` where Python does that.

use std::path::{Path, PathBuf};

/// Values of `CONFIG['combine_insertions']` + the two keys combine reads from other sections.
#[derive(Clone, Debug)]
pub struct Config {
    // ---- required (python `CONFIG['combine_insertions'][k]`)
    pub genome_2bit: String,
    pub exclude_files_with_many_insertions: i64,
    pub samtools_executable: String,
    pub bowtie2_executable: String,
    pub bowtie2_index: String,
    pub bowtie2_index2: String,
    pub bowtie2_index2_lo: String,
    /// `CONFIG['genotyping']['max_bases']` (required)
    pub genotyping_max_bases: usize,

    // ---- remap filters
    /// `clean_remap_max_insertion`, default `CONFIG['discovery']['min_clip_len']` (then required)
    pub clean_remap_max_insertion: i64,
    /// `clean_remap_min_as` (-15)
    pub clean_remap_min_as: i64,
    /// `trim_far_flank_before_remap` (False)
    pub trim_far_flank_before_remap: bool,

    // ---- intersect (combine_insertions_intersect_insertions.py:178-183)
    /// `keep_polya_one_sided` (False)
    pub keep_polya_one_sided: bool,
    /// `int(merge_tolerance_bp or 0)` (0)
    pub merge_tolerance_bp: i64,
    /// `polya_aware_clip_agreement` (False)
    pub polya_aware_clip_agreement: bool,

    // ---- evidence (combine_insertions_evidence.py)
    /// `require_independent_fragments` (False). DROPPED in the port: true is a startup error.
    pub require_independent_fragments: bool,
    /// `min_independent_fragments` (2): `supported` threshold + consensus depth floor
    pub min_independent_fragments: usize,
    /// `indel_aware_consensus` (False)
    pub indel_aware_consensus: bool,
    /// `dup_coord_tolerance` (5)
    pub dup_coord_tolerance: i64,
    /// `dup_max_edit` (3)
    pub dup_max_edit: i64,
    /// `dup_max_edit_frac` (0.02)
    pub dup_max_edit_frac: f64,
    /// `polya_min_len` (8) -- dedup cut, consensus poly-A, slippage `polya_min`
    pub polya_min_len: usize,
    /// `dup_mate_min_mapq` (20)
    pub dup_mate_min_mapq: i64,
    /// `count_short_overhang` (False)
    pub count_short_overhang: bool,
    /// `short_overhang_min_bases` (5)
    pub short_overhang_min_bases: usize,
    /// `short_overhang_min_ref_mismatch` (2)
    pub short_overhang_min_ref_mismatch: i64,
    /// `short_mate_max_dist` (1000)
    pub short_mate_max_dist: i64,
    /// `short_mate_min_mapq` (20)
    pub short_mate_min_mapq: i64,

    // ---- TPRT combine filters (combine_insertions_tprt_filters.py)
    /// `slippage_reject` (False)
    pub slippage_reject: bool,
    /// `far_pair_strict` (False) -- also switches on discovery_breakpoints in the driver
    pub far_pair_strict: bool,
    /// `far_pair_split` (True)
    pub far_pair_split: bool,
    /// `cfg.get("rte_library") or "resources/rte_library"`, resolved (see `resolve_rte_library`)
    pub rte_library: String,
    /// `slippage_min_ref_run_combine` (8)
    pub slippage_min_ref_run_combine: usize,
    /// `slippage_min_str_len` (12)
    pub slippage_min_str_len: usize,
    /// `slippage_min_structured` (10)
    pub slippage_min_structured: usize,
    /// `slippage_max_period` (6)
    pub slippage_max_period: usize,
    /// `slippage_junk_frac` (0.5)
    pub slippage_junk_frac: f64,
    /// `far_pair_max_tsd_deletion` (30)
    pub far_pair_max_tsd_deletion: i64,
    /// `far_pair_tsd_max` (40)
    pub far_pair_tsd_max: i64,
    /// `far_pair_min_polya` (10)
    pub far_pair_min_polya: usize,
    /// `far_pair_allow_antisense` (False)
    pub far_pair_allow_antisense: bool,
    /// `far_pair_colony_frac` (0.2)
    pub far_pair_colony_frac: f64,
    /// `far_pair_colony_tol` (5) (python `int(...)`)
    pub far_pair_colony_tol: i64,

    /// directory of the config file (python repo-root fallback for relative paths)
    pub config_dir: Option<PathBuf>,
}

use serde_json::Value;

/// Python truthiness of a JSON value.
fn truthy(v: &Value) -> bool {
    match v {
        Value::Null => false,
        Value::Bool(b) => *b,
        Value::Number(n) => n.as_f64().map(|x| x != 0.0).unwrap_or(true),
        Value::String(s) => !s.is_empty(),
        Value::Array(a) => !a.is_empty(),
        Value::Object(o) => !o.is_empty(),
    }
}

/// `v` as a Python number (bools are ints), None for non-numbers.
fn as_num(v: &Value) -> Option<f64> {
    match v {
        Value::Bool(b) => Some(*b as i64 as f64),
        Value::Number(n) => n.as_f64(),
        _ => None,
    }
}

struct Sect<'a>(&'a serde_json::Map<String, Value>);

impl<'a> Sect<'a> {
    /// key present and not null (a null value is treated like an absent key)
    fn get(&self, k: &str) -> Option<&'a Value> {
        self.0.get(k).filter(|v| !v.is_null())
    }
    fn req_str(&self, k: &str) -> Result<String, String> {
        match self.get(k) {
            Some(Value::String(s)) => Ok(s.clone()),
            Some(_) => Err(format!("config: combine_insertions.{k} must be a string")),
            None => Err(format!("config: required key combine_insertions.{k} is missing")),
        }
    }
    fn int(&self, k: &str, d: i64) -> Result<i64, String> {
        match self.get(k) {
            None => Ok(d),
            Some(v) => as_num(v).map(|x| x.trunc() as i64).ok_or_else(|| format!("config: combine_insertions.{k} must be a number")),
        }
    }
    fn req_int(&self, k: &str) -> Result<i64, String> {
        match self.get(k) {
            None => Err(format!("config: required key combine_insertions.{k} is missing")),
            Some(_) => self.int(k, 0),
        }
    }
    fn uint(&self, k: &str, d: usize) -> Result<usize, String> {
        let v = self.int(k, d as i64)?;
        if v < 0 {
            return Err(format!("config: combine_insertions.{k} must be >= 0"));
        }
        Ok(v as usize)
    }
    fn float(&self, k: &str, d: f64) -> Result<f64, String> {
        match self.get(k) {
            None => Ok(d),
            Some(v) => as_num(v).ok_or_else(|| format!("config: combine_insertions.{k} must be a number")),
        }
    }
    fn flag(&self, k: &str, d: bool) -> bool {
        self.get(k).map(truthy).unwrap_or(d)
    }
}

impl Config {
    /// Load from `--config` (see module doc). Missing required keys are an error naming the key.
    pub fn load(path: &Path) -> Result<Config, String> {
        let dir = path.parent().map(|p| if p.as_os_str().is_empty() { PathBuf::from(".") } else { p.to_path_buf() });
        let text = if path.extension().map(|e| e == "json").unwrap_or(false) {
            std::fs::read_to_string(path).map_err(|e| format!("cannot read config {}: {e}", path.display()))?
        } else {
            let py = std::env::var("PEARTREE_PYTHON").ok().filter(|s| !s.is_empty()).unwrap_or_else(|| "python3".to_string());
            let script = "import json, runpy, sys; print(json.dumps(runpy.run_path(sys.argv[1])[\"CONFIG\"], default=lambda o: None))";
            let out = std::process::Command::new(&py)
                .arg("-c")
                .arg(script)
                .arg(path)
                .output()
                .map_err(|e| format!("cannot run {py} to read config {}: {e}", path.display()))?;
            if !out.status.success() {
                return Err(format!("reading config {} with {py} failed: {}", path.display(), String::from_utf8_lossy(&out.stderr)));
            }
            String::from_utf8(out.stdout).map_err(|e| format!("config output is not UTF-8: {e}"))?
        };
        let v: Value = serde_json::from_str(&text).map_err(|e| format!("config {} is not valid JSON: {e}", path.display()))?;
        Config::from_json(&v, dir)
    }

    /// Build from an already-parsed JSON value of the whole CONFIG dict (unit tests use this).
    pub fn from_json(v: &serde_json::Value, config_dir: Option<PathBuf>) -> Result<Config, String> {
        let root = v.as_object().ok_or("config: top level must be an object")?;
        let ci = match root.get("combine_insertions") {
            Some(Value::Object(o)) => Sect(o),
            _ => return Err("config: required section combine_insertions is missing".into()),
        };
        let genotyping_max_bases = match root.get("genotyping").and_then(|g| g.get("max_bases")).filter(|v| !v.is_null()) {
            Some(v) => as_num(v).ok_or("config: genotyping.max_bases must be a number")?.trunc().max(0.0) as usize,
            None => return Err("config: required key genotyping.max_bases is missing".into()),
        };
        let clean_remap_max_insertion = match ci.get("clean_remap_max_insertion") {
            Some(_) => ci.int("clean_remap_max_insertion", 0)?,
            None => match root.get("discovery").and_then(|d| d.get("min_clip_len")).filter(|v| !v.is_null()).and_then(as_num) {
                Some(x) => x.trunc() as i64,
                None => return Err("config: combine_insertions.clean_remap_max_insertion is missing and so is discovery.min_clip_len".into()),
            },
        };
        let require_independent_fragments = ci.flag("require_independent_fragments", false);
        if require_independent_fragments {
            return Err("config: require_independent_fragments is not supported by the Rust port (the pooled fragment gate was dropped; see SPEC.md section 0)".into());
        }
        let rte_library = match ci.get("rte_library") {
            Some(Value::String(s)) if !s.is_empty() => s.clone(),
            Some(x) if truthy(x) => return Err("config: combine_insertions.rte_library must be a string".into()),
            _ => "resources/rte_library".to_string(),
        };
        Ok(Config {
            genome_2bit: ci.req_str("genome_2bit")?,
            exclude_files_with_many_insertions: ci.req_int("exclude_files_with_many_insertions")?,
            samtools_executable: ci.req_str("samtools_executable")?,
            bowtie2_executable: ci.req_str("bowtie2_executable")?,
            bowtie2_index: ci.req_str("bowtie2_index")?,
            bowtie2_index2: ci.req_str("bowtie2_index2")?,
            bowtie2_index2_lo: ci.req_str("bowtie2_index2_lo")?,
            genotyping_max_bases,
            clean_remap_max_insertion,
            clean_remap_min_as: ci.int("clean_remap_min_as", -15)?,
            trim_far_flank_before_remap: ci.flag("trim_far_flank_before_remap", false),
            keep_polya_one_sided: ci.flag("keep_polya_one_sided", false),
            // python: int(v or 0)
            merge_tolerance_bp: match ci.get("merge_tolerance_bp") {
                Some(x) if truthy(x) => ci.int("merge_tolerance_bp", 0)?,
                _ => 0,
            },
            polya_aware_clip_agreement: ci.flag("polya_aware_clip_agreement", false),
            require_independent_fragments,
            min_independent_fragments: ci.uint("min_independent_fragments", 2)?,
            indel_aware_consensus: ci.flag("indel_aware_consensus", false),
            dup_coord_tolerance: ci.int("dup_coord_tolerance", 5)?,
            dup_max_edit: ci.int("dup_max_edit", 3)?,
            dup_max_edit_frac: ci.float("dup_max_edit_frac", 0.02)?,
            polya_min_len: ci.uint("polya_min_len", 8)?,
            dup_mate_min_mapq: ci.int("dup_mate_min_mapq", 20)?,
            count_short_overhang: ci.flag("count_short_overhang", false),
            short_overhang_min_bases: ci.uint("short_overhang_min_bases", 5)?,
            short_overhang_min_ref_mismatch: ci.int("short_overhang_min_ref_mismatch", 2)?,
            short_mate_max_dist: ci.int("short_mate_max_dist", 1000)?,
            short_mate_min_mapq: ci.int("short_mate_min_mapq", 20)?,
            slippage_reject: ci.flag("slippage_reject", false),
            far_pair_strict: ci.flag("far_pair_strict", false),
            far_pair_split: ci.flag("far_pair_split", true),
            rte_library,
            slippage_min_ref_run_combine: ci.uint("slippage_min_ref_run_combine", 8)?,
            slippage_min_str_len: ci.uint("slippage_min_str_len", 12)?,
            slippage_min_structured: ci.uint("slippage_min_structured", 10)?,
            slippage_max_period: ci.uint("slippage_max_period", 6)?,
            slippage_junk_frac: ci.float("slippage_junk_frac", 0.5)?,
            far_pair_max_tsd_deletion: ci.int("far_pair_max_tsd_deletion", 30)?,
            far_pair_tsd_max: ci.int("far_pair_tsd_max", 40)?,
            far_pair_min_polya: ci.uint("far_pair_min_polya", 10)?,
            far_pair_allow_antisense: ci.flag("far_pair_allow_antisense", false),
            far_pair_colony_frac: ci.float("far_pair_colony_frac", 0.2)?,
            far_pair_colony_tol: ci.int("far_pair_colony_tol", 5)?,
            config_dir,
        })
    }

    /// tprt_filters._resolve: a relative `rte_library` that does not exist relative to the cwd is
    /// tried relative to the repository root (python: `<dir of src/>/..`). Rust candidates, in
    /// order: cwd, `<config_dir>/..` (config lives in src/), `$PT_ROOT`, the binary's
    /// `../../../..` (rust/peartree-combine/target/release). First existing wins; else as given.
    pub fn resolve_rte_library(&self) -> PathBuf {
        let given = PathBuf::from(&self.rte_library);
        if self.rte_library.is_empty() || given.is_absolute() || given.exists() {
            return given;
        }
        let mut roots: Vec<PathBuf> = Vec::new();
        if let Some(d) = &self.config_dir {
            roots.push(d.join(".."));
        }
        if let Ok(r) = std::env::var("PT_ROOT") {
            if !r.is_empty() {
                roots.push(PathBuf::from(r));
            }
        }
        if let Ok(exe) = std::env::current_exe() {
            let mut p = exe;
            for _ in 0..5 {
                // exe -> release -> target -> peartree-combine -> rust -> repo
                p = match p.parent() {
                    Some(q) => q.to_path_buf(),
                    None => break,
                };
            }
            roots.push(p);
        }
        for r in roots {
            let alt = r.join(&given);
            if alt.exists() {
                return alt;
            }
        }
        given
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn base() -> serde_json::Value {
        json!({
            "combine_insertions": {
                "genome_2bit": "g.2bit", "exclude_files_with_many_insertions": 5000,
                "samtools_executable": "samtools", "bowtie2_executable": "bowtie2",
                "bowtie2_index": "i", "bowtie2_index2": "i2", "bowtie2_index2_lo": "lo"
            },
            "genotyping": {"max_bases": 50},
            "discovery": {"min_clip_len": 12}
        })
    }

    #[test]
    fn defaults_match_python() {
        let c = Config::from_json(&base(), None).unwrap();
        assert_eq!(c.clean_remap_max_insertion, 12);
        assert_eq!(c.clean_remap_min_as, -15);
        assert!(!c.keep_polya_one_sided && !c.far_pair_strict && c.far_pair_split);
        assert_eq!(c.merge_tolerance_bp, 0);
        assert_eq!((c.min_independent_fragments, c.polya_min_len, c.dup_coord_tolerance), (2, 8, 5));
        assert_eq!((c.dup_max_edit, c.dup_max_edit_frac, c.dup_mate_min_mapq), (3, 0.02, 20));
        assert_eq!((c.slippage_min_ref_run_combine, c.slippage_min_str_len, c.slippage_min_structured, c.slippage_max_period), (8, 12, 10, 6));
        assert_eq!((c.far_pair_max_tsd_deletion, c.far_pair_tsd_max, c.far_pair_min_polya, c.far_pair_colony_tol), (30, 40, 10, 5));
        assert_eq!((c.far_pair_colony_frac, c.slippage_junk_frac), (0.2, 0.5));
        assert_eq!(c.rte_library, "resources/rte_library");
        assert_eq!(c.genotyping_max_bases, 50);
    }

    #[test]
    fn overrides_and_truthiness() {
        let mut v = base();
        let ci = v["combine_insertions"].as_object_mut().unwrap();
        ci.insert("merge_tolerance_bp".into(), json!(7));
        ci.insert("keep_polya_one_sided".into(), json!(1));
        ci.insert("slippage_reject".into(), json!("yes"));
        ci.insert("far_pair_split".into(), json!(0));
        ci.insert("rte_library".into(), json!(""));
        ci.insert("clean_remap_max_insertion".into(), json!(30));
        ci.insert("far_pair_colony_tol".into(), json!(3.9));
        let c = Config::from_json(&v, None).unwrap();
        assert_eq!(c.merge_tolerance_bp, 7);
        assert!(c.keep_polya_one_sided && c.slippage_reject && !c.far_pair_split);
        assert_eq!(c.rte_library, "resources/rte_library");
        assert_eq!(c.clean_remap_max_insertion, 30);
        assert_eq!(c.far_pair_colony_tol, 3);
    }

    #[test]
    fn errors() {
        let mut v = base();
        v["combine_insertions"].as_object_mut().unwrap().remove("bowtie2_index");
        assert!(Config::from_json(&v, None).unwrap_err().contains("bowtie2_index"));
        let mut v = base();
        v["combine_insertions"]["require_independent_fragments"] = json!(true);
        assert!(Config::from_json(&v, None).is_err());
        let mut v = base();
        v.as_object_mut().unwrap().remove("discovery");
        assert!(Config::from_json(&v, None).is_err());
    }

    #[test]
    fn loads_json_file_and_resolves_library() {
        let dir = std::env::temp_dir().join(format!("pt_combine_p1_cfg_{}", std::process::id()));
        std::fs::create_dir_all(dir.join("src")).unwrap();
        std::fs::create_dir_all(dir.join("lib_unique_p1")).unwrap();
        let p = dir.join("src").join("config.json");
        std::fs::write(&p, base().to_string()).unwrap();
        let mut c = Config::load(&p).unwrap();
        c.rte_library = "lib_unique_p1".into();
        assert_eq!(c.resolve_rte_library(), dir.join("src").join("..").join("lib_unique_p1"));
        c.rte_library = "does_not_exist_p1".into();
        assert_eq!(c.resolve_rte_library(), PathBuf::from("does_not_exist_p1"));
        std::fs::remove_dir_all(&dir).ok();
    }
}
