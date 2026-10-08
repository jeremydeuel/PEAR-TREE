//! tools/rte/hallmarks.py -- TSD, EN motif, poly-A, slippage, fold-back.
//!
//! STATUS: types FOUNDATION; functions = WORK PACKAGE "WP-HALL" (PORT_PLAN.md).
//!
//! Python functions to port (hallmarks.py at eaa2718): split_junction (27), parse_locus (37),
//! edge_run (44), polya_info (84), _find_near (124) (python `re.finditer` of the exact probe ->
//! all overlapping occurrences, nearest to the hint, first on ties; else edlib_best + strand only),
//! locate_site (140), tsd_from_flanks (158), target_site (172), en_motif (194), en_bin (209),
//! slippage_context (221), foldback (247), _longest_prefix_match (262).
//!
//! Golden: events `hallmark.<function name>` (in: the python arguments, genome as a genome id;
//! out: the return value / the mutated SiteInfo).

use crate::genome::Genome;
use crate::inputs::JunctionEvidence;

/// (left_insert, left_flank, right_flank, right_insert), all upper case.
pub type SplitJunction = (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>);

/// `split_junction(left_seq, right_seq)`. WP-HALL.
pub fn split_junction(_left_seq: &[u8], _right_seq: &[u8]) -> SplitJunction {
    todo!("WP-HALL: port hallmarks.split_junction")
}

/// `parse_locus(title)` -> (contig, a, b) for `contig:a-b` (a, b may be negative). WP-HALL.
pub fn parse_locus(_title: &str) -> Option<(String, i64, i64)> {
    todo!("WP-HALL: port hallmarks.parse_locus")
}

/// `edge_run(insert, at_end)` -> (base (0 = none), length). WP-HALL.
pub fn edge_run(_insert: &[u8], _at_end: bool) -> (u8, i64) {
    todo!("WP-HALL: port hallmarks.edge_run")
}

/// hallmarks.PolyAInfo
#[derive(Clone, Debug, PartialEq)]
pub struct PolyAInfo {
    /// +1 poly-A at LEFT, -1 poly-T at RIGHT, 0 unknown
    pub strand: i32,
    pub source: String,
    /// REF-adjacent homopolymer of the LEFT insert (base 0 = none)
    pub left_run: (u8, i64),
    pub right_run: (u8, i64),
    pub both_sided: bool,
    pub length: f64,
}

impl Default for PolyAInfo {
    fn default() -> Self {
        PolyAInfo { strand: 0, source: "none".into(), left_run: (0, 0), right_run: (0, 0), both_sided: false, length: 0.0 }
    }
}

/// `polya_info(left_insert, right_insert, ev_left, ev_right, min_len=10)`. WP-HALL.
pub fn polya_info(
    _left_insert: &[u8],
    _right_insert: &[u8],
    _ev_left: Option<&JunctionEvidence>,
    _ev_right: Option<&JunctionEvidence>,
    _min_len: f64,
) -> PolyAInfo {
    todo!("WP-HALL: port hallmarks.polya_info")
}

/// hallmarks.SiteInfo
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SiteInfo {
    pub contig: Option<String>,
    pub l: Option<i64>,
    pub r: Option<i64>,
    pub located: bool,
    pub tsd_seq: String,
    /// > 0 duplication, < 0 deletion, 0 blunt, None unknown
    pub tsd_len: Option<i64>,
    pub tsd_verified: bool,
    pub en_motif: String,
    pub en_mismatches: Option<i64>,
    pub slippage: bool,
    pub slippage_detail: String,
}

/// `locate_site(title, left_flank, right_flank, genome, probe_len=30, slack=1500)`. WP-HALL.
pub fn locate_site(_title: &str, _left_flank: &[u8], _right_flank: &[u8], _genome: Option<&dyn Genome>) -> SiteInfo {
    todo!("WP-HALL: port hallmarks.locate_site")
}

/// `tsd_from_flanks(left_flank, right_flank, max_len=80, min_len=4, max_mm=1)` ->
/// (seq, len or None, verified). WP-HALL.
pub fn tsd_from_flanks(_left_flank: &[u8], _right_flank: &[u8]) -> (String, Option<i64>, bool) {
    todo!("WP-HALL: port hallmarks.tsd_from_flanks")
}

/// `target_site(si, left_flank, right_flank, genome)` (fills tsd_* in place). WP-HALL.
pub fn target_site(_si: &mut SiteInfo, _left_flank: &[u8], _right_flank: &[u8], _genome: Option<&dyn Genome>) {
    todo!("WP-HALL: port hallmarks.target_site")
}

/// `en_motif(si, strand, genome)` (in place). WP-HALL.
pub fn en_motif(_si: &mut SiteInfo, _strand: i32, _genome: Option<&dyn Genome>) {
    todo!("WP-HALL: port hallmarks.en_motif")
}

/// `en_bin(mm)`. WP-HALL.
pub fn en_bin(_mm: Option<i64>) -> &'static str {
    todo!("WP-HALL: port hallmarks.en_bin")
}

/// `slippage_context(si, strand, genome, min_run=10)` (in place). WP-HALL.
pub fn slippage_context(_si: &mut SiteInfo, _strand: i32, _genome: Option<&dyn Genome>) {
    todo!("WP-HALL: port hallmarks.slippage_context")
}

/// `foldback(left_insert, left_flank, right_insert, right_flank, min_len=15, max_mm=2)`. WP-HALL.
pub fn foldback(_left_insert: &[u8], _left_flank: &[u8], _right_insert: &[u8], _right_flank: &[u8]) -> i64 {
    todo!("WP-HALL: port hallmarks.foldback")
}
