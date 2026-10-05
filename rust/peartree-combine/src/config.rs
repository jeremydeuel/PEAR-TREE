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

impl Config {
    /// Load from `--config` (see module doc). Missing required keys are an error naming the key.
    pub fn load(path: &Path) -> Result<Config, String> {
        todo!("P1: SPEC.md §1.2")
    }

    /// Build from an already-parsed JSON value of the whole CONFIG dict (unit tests use this).
    pub fn from_json(v: &serde_json::Value, config_dir: Option<PathBuf>) -> Result<Config, String> {
        todo!("P1: SPEC.md §1.2 -- every key + default listed on the fields above")
    }

    /// tprt_filters._resolve: a relative `rte_library` that does not exist relative to the cwd is
    /// tried relative to the repository root (python: `<dir of src/>/..`). Rust candidates, in
    /// order: cwd, `<config_dir>/..` (config lives in src/), `$PT_ROOT`, the binary's
    /// `../../../..` (rust/peartree-combine/target/release). First existing wins; else as given.
    pub fn resolve_rte_library(&self) -> PathBuf {
        todo!("P1")
    }
}
