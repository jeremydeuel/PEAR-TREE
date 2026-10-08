//! tools/rte/structure.py -- element / 5' structure / tags from the element-sense layouts.
//!
//! STATUS: types + callback interfaces FOUNDATION; algorithms = WORK PACKAGE "WP-STRUCT".
//!
//! Python functions to port (structure.py at eaa2718): classify (106), _first_element_after_ref
//! (87), _last_element_before_ref (97), _chain (405), _inversion (431), _element_after (476),
//! _unexplained_tail (482), _template_abs (521), _foldback_5p (533), _local_templates (564),
//! _sense_switch (604), _inverted_tail_5p (622), _flank_masked_frac (636), _source_class_ok
//! (646), _at_rich (657), _td5_at_junction (664), _explained_by_consensus (674),
//! _ref_at_breakpoint (684).
//!
//! Python identity semantics to keep: `_local_templates(exclude={id(x) ...})` excludes the fold-
//! back segments by OBJECT identity -> use (layout index, segment index) here; `e in f5`
//! (structure.py:308) is dataclass EQUALITY (all fields) -> `Segment: PartialEq`.
//!
//! Callbacks (so this package does not depend on transduction / pseudogene / the gene model):
//! the annotator passes a [`SourceFinder`] (NovelSourceFinder), a pre-mRNA labeller and the
//! pseudogene arguments; the golden `classify` events record every callback call with its answer
//! so the port is testable alone (golden.rs `ReplayFinder` / `replay_premrna`).
//!
//! Golden: events `classify` (in: assembly, ctx, cfg, legacy_class, pseudogene genes/hits,
//! recorded callback answers; out: StructureCall).

use crate::assembly::{AssemblyResult, ReadLayout, SiteContext};
use crate::config::StructureCfg;
use crate::library::RteLibrary;
use crate::record::Detail;
use crate::transduction::SourceCall;

/// structure.StructureCall
#[derive(Clone, Debug, PartialEq)]
pub struct StructureCall {
    pub element: String,
    pub structure: String,
    pub tags: Vec<String>,
    pub detail: Detail,
    pub source: Option<SourceCall>,
    pub j5_class: String,
    pub j3_class: String,
    pub j5_pos: Option<i64>,
    pub inv_p1: Option<i64>,
    pub three_prime_truncated: bool,
    pub three_prime_short: bool,
    pub has_polya_3p: bool,
    pub td3p_seq: Vec<u8>,
}

impl Default for StructureCall {
    fn default() -> Self {
        StructureCall {
            element: "UNKNOWN".into(),
            structure: "5P_UNRESOLVED".into(),
            tags: Vec::new(),
            detail: Detail::default(),
            source: None,
            j5_class: String::new(),
            j3_class: String::new(),
            j5_pos: None,
            inv_p1: None,
            three_prime_truncated: false,
            three_prime_short: false,
            has_polya_3p: false,
            td3p_seq: Vec::new(),
        }
    }
}

impl StructureCall {
    /// `add(tag)`: append when not present
    pub fn add(&mut self, tag: &str) {
        if !self.tags.iter().any(|t| t == tag) {
            self.tags.push(tag.to_string());
        }
    }
}

/// `novel_finder.find(seq)` (transduction.NovelSourceFinder implements it).
pub trait SourceFinder: Sync {
    fn find(&self, seq: &[u8]) -> Option<SourceCall>;
}

/// `premrna(seq) -> label or None` (annotator._premrna_fn).
pub type PremrnaFn<'a> = &'a (dyn Fn(&[u8]) -> Option<String> + Sync);

/// `pseudogene_structure_fn(layouts) -> 'FULL_LENGTH' | 'TRUNCATED_5P' | None`.
pub type PgStructureFn<'a> = &'a (dyn Fn(&[ReadLayout]) -> Option<String> + Sync);

/// python's `pseudogene=(candidate_genes, exon_junction_hits, structure_fn)` argument.
pub struct PseudogeneArg<'a> {
    pub genes: &'a [String],
    /// [(junction label, read name)]
    pub hits: &'a [(String, String)],
    pub structure_fn: Option<PgStructureFn<'a>>,
}

/// `classify(res, lib, ctx, cfg, novel_finder, premrna, pseudogene, legacy_class)`. WP-STRUCT.
#[allow(clippy::too_many_arguments)]
pub fn classify(
    _res: &AssemblyResult,
    _lib: &RteLibrary,
    _ctx: Option<&SiteContext>,
    _cfg: &StructureCfg,
    _novel_finder: Option<&dyn SourceFinder>,
    _premrna: Option<PremrnaFn>,
    _pseudogene: Option<&PseudogeneArg>,
    _legacy_class: Option<&str>,
) -> StructureCall {
    todo!("WP-STRUCT: port structure.classify")
}
