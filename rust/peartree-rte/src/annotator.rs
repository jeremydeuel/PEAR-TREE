//! tools/rte/annotator.py -- orchestration (RteAnnotator) + the cross-locus passes. WP-INT.
//!
//! Per-insertion data flow (python `_annotate_key` -> `annotate`):
//!   LocusData (stream.rs) -> pool combine reads + GT reads -> cap_reads(rte_max_reads)
//!   -> split_junction / polya_info / locate_site / target_site (hallmarks)
//!   -> SiteContext windows (genome) -> Assembler.assemble -> supported tails / polya_reads
//!   -> en_motif / slippage_context -> ExonJunctionIndex.find (if candidate genes)
//!   -> structure.classify (NovelSourceFinder, pre-mRNA fn, pseudogene structure fn)
//!   -> site-level tags -> beyond-poly-A -> ScoreInput -> score -> RteRecord
//!   (+ a combine-reads-only rerun for `gt_changed` when the locus has GT reads and
//!   rte_gt_compare).
//! Cross-locus (python `annotate_all`, needs every record first; records are small, kept in
//! memory in input order):
//!   1. cohort_source_pass: only when genome_2bit == remap_2bit and a locator exists; collects
//!      L1 TPRT/LIKELY_TPRT sites, then RE-ANNOTATES (re-loads the locus reads from the store)
//!      every TD3P record without a TD3P_SOURCE= tag with cohort_l1 set; an error keeps the old
//!      record.
//!   2. recurrence_pass: signature (element, structure, int(j5)//5) over L1/ALU/SVA records not
//!      FULL_LENGTH / 5P_UNRESOLVED with detail j5 (or fwd_start); count > rte_recurrence_max ->
//!      score_input.recurrent = True, detail["recurrence"] = count, re-score.
//!
//! Errors: any failure annotating one locus -> RteRecord(key) with detail["error"] =
//! f"{type(e).__name__}: {e}"[:200]. In Rust a per-locus panic is caught (`catch_unwind`); the
//! ported modules panic with python's exception text where python raises (e.g.
//! `ValueError: max() iterable argument is empty`), anything else is reported as
//! `RustPanic: <message>`. An I/O error reading the spilled reads is NOT a locus error: it fails
//! the run.

use crate::assembly::{counts_as_fragment, strand_from_polya, AssemblyResult, Assembler, JunctionSeq, ReadLayout, SegKind, SiteContext};
use crate::config::RteConfig;
use crate::genemodel::GeneModel;
use crate::genome::{open_genome, Genome};
use crate::hallmarks::{edge_run, en_motif, foldback, locate_site, parse_locus, polya_info, slippage_context, split_junction, target_site, PolyAInfo, SiteInfo};
use crate::inputs::{EvidenceRead, InsertionEvidence, JunctionEvidence, GT_ROLE_PREFIX};
use crate::library::RteLibrary;
use crate::mm::{Aligner, MapOpts};
use crate::pseudogene::{load_exons_by_gene, load_gene_strands, ExonJunctionIndex, ExonsByGene};
use crate::record::{gt_changes, RteRecord};
use crate::score::{score, ScoreInput};
use crate::stream::{cap_reads, map_loci, LocusData, LocusLoader};
use crate::structure::{classify, PremrnaFn, PseudogeneArg, SourceFinder};
use crate::transduction::{L1Rmsk, Locator, LocatorHit, MappyLocator, NovelSourceFinder};
use rustc_hash::FxHashMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// annotator.InsertionInput (what annotate_v2 hands over per insertion; SPEC.md "Input").
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InsertionInput {
    /// locus name (= insertion id)
    pub title: String,
    /// combined.txt.gz LEFT junction string (lower = clip, UPPER = reference)
    pub left_seq: Vec<u8>,
    pub right_seq: Vec<u8>,
    /// pseudogene candidate genes (InsertionInput.from_legacy order)
    pub pseudogene_genes: Vec<String>,
    /// annotate_v2 element_class(ins.conclusion())
    pub legacy_class: Option<String>,
    /// annotate_v2 `_sv_subtype(allow_rte=True)` rank (only `sv[0]` is read: rank 2 =
    /// intrachromosomal partner)
    pub sv_rank: Option<i64>,
}

/// `_richer_junction(evidence, combined)`: the string with more lower-case (clipped) bases; the
/// evidence consensus on a tie / when combined is empty.
pub fn richer_junction(evidence: &[u8], combined: &[u8]) -> Vec<u8> {
    if evidence.is_empty() {
        return combined.to_vec();
    }
    let lc = |s: &[u8]| s.iter().filter(|c| c.is_ascii_lowercase()).count();
    if !combined.is_empty() && lc(combined) > lc(evidence) {
        combined.to_vec()
    } else {
        evidence.to_vec()
    }
}

/// `default_sidecars(insertions_file)` -> (evidence.tsv.gz, reads.fa.gz)
pub fn default_sidecars(insertions_file: &str) -> (String, String) {
    let mut base = insertions_file;
    for suf in [".combined.txt.gz", ".txt.gz"] {
        if let Some(b) = base.strip_suffix(suf) {
            base = b;
            break;
        }
    }
    (format!("{base}.insertions.evidence.tsv.gz"), format!("{base}.insertions.reads.fa.gz"))
}

/// `default_gt_reads(insertions_file)`
pub fn default_gt_reads(insertions_file: &str) -> String {
    default_sidecars(insertions_file).1.replace(".insertions.reads.fa.gz", ".insertions.genotype_reads.fa.gz")
}

/// `title_gap(title)`: R - L of a numeric locus name `contig:L-R` (None otherwise).
pub fn title_gap(title: &str) -> Option<i64> {
    parse_locus(title).map(|(_, l, r)| r - l)
}

/// python `MappyLocator(path)`: the remap index is only loaded at the first query (an hs1
/// minimap2 index is several GB; python builds it lazily too).
pub struct LazyLocator {
    path: PathBuf,
    al: OnceLock<Option<MappyLocator>>,
}

impl LazyLocator {
    pub fn new(path: &Path) -> LazyLocator {
        LazyLocator { path: path.to_path_buf(), al: OnceLock::new() }
    }
}

impl Locator for LazyLocator {
    fn locate(&self, seq: &[u8]) -> Vec<LocatorHit> {
        match self.al.get_or_init(|| MappyLocator::open(&self.path)) {
            Some(l) => l.locate(seq),
            None => panic!("OSError: cannot load the minimap2 index {}", self.path.display()),
        }
    }
}

/// Everything loaded once per run (python `RteAnnotator.__init__`): library, discovery / remap
/// genomes, novel-source locator + rmsk, exon track, gene model.
pub struct Resources {
    pub lib: RteLibrary,
    pub genome: Option<Box<dyn Genome>>,
    pub remap: Option<Box<dyn Genome>>,
    pub locator: Option<Box<dyn Locator>>,
    pub rmsk: Option<L1Rmsk>,
    pub exons: ExonsByGene,
    pub strands: FxHashMap<String, char>,
    pub gene_model: Option<GeneModel>,
}

impl Resources {
    /// Load every configured resource (python RteAnnotator.__init__ + annotate_v2's
    /// `read_gene_model`, whose GeneModel annotate_v2 hands to the annotator).
    pub fn load(cfg: &RteConfig) -> Result<Resources, String> {
        let lib = RteLibrary::open(&cfg.rte_library, cfg.young_consensus_regex.as_deref())?;
        let genome = open_genome(cfg.genome_2bit.as_deref())?;
        let remap = open_genome(cfg.remap_2bit.as_deref())?;
        let locator: Option<Box<dyn Locator>> = match cfg.remap_index.as_deref().filter(|s| !s.is_empty()) {
            Some(p) => {
                if !Path::new(p).exists() {
                    return Err(format!("remap_index {p}: no such file"));
                }
                Some(Box::new(LazyLocator::new(Path::new(p))))
            }
            None => None,
        };
        let rmsk = match cfg.remap_rmsk.as_deref().filter(|s| !s.is_empty()) {
            Some(p) => Some(L1Rmsk::open(Path::new(p), cfg.transduction.novel_source_min_len)?),
            None => None,
        };
        let (mut exons, mut strands) = (ExonsByGene::default(), FxHashMap::default());
        if let Some(t) = cfg.exon_track().filter(|t| !t.is_empty() && Path::new(t).exists()) {
            exons = load_exons_by_gene(Path::new(t))?;
            strands = load_gene_strands(Path::new(t))?;
        }
        let gene_model = match cfg.gene_model.as_deref().filter(|s| !s.is_empty()) {
            Some(p) => Some(GeneModel::open(Path::new(p), &cfg.gene_model_cfg)?),
            None => None,
        };
        Ok(Resources { lib, genome, remap, locator, rmsk, exons, strands, gene_model })
    }
}

/// annotator.RteAnnotator
pub struct RteAnnotator<'a> {
    pub cfg: &'a RteConfig,
    pub lib: &'a RteLibrary,
    pub assembler: Assembler<'a>,
    pub genome: Option<&'a dyn Genome>,
    pub novel: NovelSourceFinder<'a>,
    pub exon_index: ExonJunctionIndex<'a>,
    pub gene_model: Option<&'a GeneModel>,
    /// loci in flight per parallel chunk (`--chunk`)
    pub chunk: usize,
    /// golden replays only: answer classify's novel-source / pre-mRNA queries from a recording
    pub source_finder_override: Option<&'a dyn SourceFinder>,
    pub premrna_override: Option<PremrnaFn<'a>>,
}

/// One locus result of the parallel map: the record, or a fatal (I/O) error.
type LocusResult = Result<RteRecord, String>;

impl<'a> RteAnnotator<'a> {
    pub fn new(cfg: &'a RteConfig, res: &'a Resources) -> RteAnnotator<'a> {
        RteAnnotator::with_parts(
            cfg,
            &res.lib,
            res.genome.as_deref(),
            res.remap.as_deref(),
            res.locator.as_deref(),
            res.rmsk.as_ref(),
            res.exons.clone(),
            res.strands.clone(),
            res.gene_model.as_ref(),
        )
    }

    /// The annotator from borrowed parts (`new` = from a [`Resources`]; tests pass their own).
    #[allow(clippy::too_many_arguments)]
    pub fn with_parts(
        cfg: &'a RteConfig,
        lib: &'a RteLibrary,
        genome: Option<&'a dyn Genome>,
        remap: Option<&'a dyn Genome>,
        locator: Option<&'a dyn Locator>,
        rmsk: Option<&'a L1Rmsk>,
        exons: ExonsByGene,
        strands: FxHashMap<String, char>,
        gene_model: Option<&'a GeneModel>,
    ) -> RteAnnotator<'a> {
        RteAnnotator {
            cfg,
            lib,
            assembler: Assembler::new(lib, cfg.assembly.clone()),
            genome,
            novel: NovelSourceFinder {
                lib,
                cfg: cfg.transduction.clone(),
                rmsk,
                locator,
                genome: remap,
                cohort_l1: Vec::new(),
                ident_cache: Default::default(),
            },
            exon_index: ExonJunctionIndex { cfg: cfg.pseudogene.clone(), exons, genome: remap, strands },
            gene_model,
            chunk: 64,
            source_finder_override: None,
            premrna_override: None,
        }
    }

    // ------------------------------------------------------------------ pre-mRNA helper
    /// `_premrna_fn(site)`: None without a discovery genome / gene model / site.
    fn premrna(&self, site: &SiteInfo) -> Option<Premrna<'a>> {
        let (genome, gm, contig) = (self.genome?, self.gene_model?, site.contig.as_ref()?);
        let anchor = site.l.or(site.r)?;
        let win = self.cfg.premrna_window;
        Some(Premrna { genome, gm, contig: contig.clone(), lo: (anchor - win).max(0), anchor, win, al: OnceLock::new() })
    }

    // ------------------------------------------------------------------ one insertion
    /// `annotate(inp, ev)`: one insertion from its (already pooled) evidence.
    pub fn annotate(&self, inp: &InsertionInput, ev: &InsertionEvidence) -> RteRecord {
        let cfg = self.cfg;
        let mut rec = RteRecord::new(&inp.title);
        let (jl, jr) = (ev.junction("LEFT"), ev.junction("RIGHT"));
        let ev_left: &[u8] = jl.map_or(&[][..], |j| j.clip_consensus.as_bytes());
        let ev_right: &[u8] = jr.map_or(&[][..], |j| j.clip_consensus.as_bytes());
        let left_str = richer_junction(ev_left, &inp.left_seq);
        let right_str = richer_junction(ev_right, &inp.right_seq);
        let (li, lf, rf, ri) = split_junction(&left_str, &right_str);
        let pa0 = polya_info(&li, &ri, jl, jr, 10.0);
        let genome = self.genome;
        let mut site = locate_site(&inp.title, &lf, &rf, genome);
        target_site(&mut site, &lf, &rf, genome);
        let mut ctx = SiteContext::new(&inp.title, site.contig.clone(), site.l, site.r, lf.clone(), rf.clone());
        if let Some(g) = genome {
            let pts: Vec<i64> = [site.l, site.r].into_iter().flatten().collect();
            if !pts.is_empty() {
                let contig = site.contig.as_deref().unwrap_or("");
                let (mn, mx) = (*pts.iter().min().unwrap(), *pts.iter().max().unwrap());
                let w = cfg.assembly.local_window;
                ctx.window_start = (mn - w).max(0);
                ctx.window_seq = g.fetch(contig, ctx.window_start, mx + w);
                let ww = cfg.wide_window;
                if ww > w {
                    ctx.wide_start = (mn - ww).max(0);
                    ctx.wide_seq = g.fetch(contig, ctx.wide_start, mx + ww);
                }
            }
        }
        let mut junction_seqs: Vec<(String, JunctionSeq)> = Vec::new();
        if !left_str.is_empty() {
            junction_seqs.push(("LEFT".into(), (left_str.clone(), (li.len() as i64, left_str.len() as i64))));
        }
        if !right_str.is_empty() {
            junction_seqs.push(("RIGHT".into(), (right_str.clone(), (0, rf.len() as i64))));
        }
        let reads = cap_reads(ev.reads.clone(), cfg.max_reads);
        rec.gt_reads = reads.iter().filter(|r| r.role.starts_with(GT_ROLE_PREFIX)).count() as i64;
        let hint = (pa0.strand != 0).then(|| (pa0.strand, pa0.source.clone()));
        let asm = self.assembler.assemble(&ctx, &junction_seqs, &reads, hint);
        let strand = asm.strand;
        // scored poly-A: the tail length reached by >= rte_polya_min_fragments fragments
        let hint_left = if pa0.left_run.0 == b'A' { pa0.left_run.1 } else { 0 };
        let hint_right = if pa0.right_run.0 == b'T' { pa0.right_run.1 } else { 0 };
        let (sup_a, sup_t) = self.supported_tails(&asm, ev_left, ev_right, jl, jr);
        let mut pa = PolyAInfo { strand: 0, left_run: pa0.left_run, right_run: pa0.right_run, ..Default::default() };
        pa.both_sided = sup_a >= 10.0 && sup_t >= 10.0;
        if sup_a >= 5.0 || sup_t >= 5.0 {
            if sup_a >= sup_t {
                (pa.strand, pa.source, pa.length) = (1, "polyA_left".into(), sup_a);
            } else {
                (pa.strand, pa.source, pa.length) = (-1, "polyT_right".into(), sup_t);
            }
        }
        if pa.strand != 0 && strand != 0 && pa.strand != strand {
            pa.strand = strand;
            pa.source = asm.strand_source.clone();
            pa.length = if strand < 0 { sup_t } else { sup_a };
            pa.both_sided = false;
        }
        let hint_len = if (if strand != 0 { strand } else { pa.strand }) < 0 { hint_right } else { hint_left };
        let polya_unsupported = (hint_len >= 5 && (hint_len as f64) > pa.length).then_some(hint_len);
        let mut polya_reads: Option<i64> = None;
        if pa.strand == 0 && polya_unsupported.is_none() && asm.strand_source == "polya_reads" {
            let (st, mut lens) = strand_from_polya(&asm.raw_layouts, 10);
            if !lens.is_empty() && st == strand {
                lens.sort();
                polya_reads = Some(lens[lens.len() / 2]);
            }
        }
        en_motif(&mut site, strand, genome);
        slippage_context(&mut site, strand, genome);
        // pseudogene proof
        let mut pg_hits: Vec<(String, String)> = Vec::new();
        if !inp.pseudogene_genes.is_empty() {
            let mut seqs: Vec<(String, Vec<u8>)> =
                junction_seqs.iter().map(|(k, v)| (format!("junction_{k}"), v.0.to_ascii_uppercase())).collect();
            seqs.extend(reads.iter().map(|r| (format!("{}|{}|{}", r.role, r.sample, r.frag), r.seq())));
            pg_hits = self.exon_index.find(&inp.pseudogene_genes, &seqs);
        }
        let premrna = self.premrna(&site);
        let premrna_call = |s: &[u8]| premrna.as_ref().and_then(|p| p.call(s));
        let premrna_fn: Option<PremrnaFn> = match self.premrna_override {
            Some(f) => Some(f),
            None => premrna.is_some().then_some(&premrna_call as PremrnaFn),
        };
        let tol = cfg.structure.pseudogene_full_length_tol;
        let genes = &inp.pseudogene_genes;
        let pg_structure = |layouts: &[ReadLayout]| -> Option<String> {
            let mut order: Vec<&ReadLayout> = layouts.iter().collect();
            order.sort_by_key(|l| l.role != "JUNCTION");
            let seqs: Vec<Vec<u8>> = order
                .iter()
                .filter(|l| l.segments.len() >= 2 && l.segments[0].kind == SegKind::Ref)
                .map(|l| l.seq[(l.segments[0].q_en.max(0) as usize).min(l.seq.len())..].to_vec())
                .collect();
            self.exon_index.structure(genes, &seqs, tol)
        };
        let pg_arg = PseudogeneArg { genes, hits: &pg_hits, structure_fn: (!genes.is_empty()).then_some(&pg_structure as _) };
        let finder: &dyn SourceFinder = self.source_finder_override.unwrap_or(&self.novel);
        let mut call =
            classify(&asm, self.lib, Some(&ctx), &cfg.structure, Some(finder), premrna_fn, Some(&pg_arg), inp.legacy_class.as_deref());
        let rte_elem = matches!(call.element.as_str(), "L1" | "ALU" | "SVA");
        // ---- site-level tags (locus-name geometry gap = R - L)
        let gap = site.tsd_len.or_else(|| title_gap(&inp.title));
        if gap.is_some() && site.tsd_len.is_none() {
            site.tsd_len = gap;
        }
        let md = cfg.max_target_site_deletion;
        if let Some(g) = gap {
            if -md <= g && g < 0 {
                call.add("TSD_DELETION");
            }
        }
        let has = |tags: &[String], t: &str| tags.iter().any(|x| x == t);
        if rte_elem {
            let sv_intra = inp.sv_rank == Some(2) && !has(&call.tags, "PREMRNA_COINSERT") && !has(&call.tags, "TEMPLATED_LOCAL");
            let polarised = pa.length >= 10.0 && !pa.both_sided && !has(&call.tags, "CHIMERIC_ENDS");
            let g0 = gap.unwrap_or(0);
            if gap.is_some_and(|g| g < -md) || (sv_intra && g0 <= 0) {
                call.add("L1_MED_DELETION");
            } else if (sv_intra && g0 > 0) || (gap.is_some_and(|g| g > 40) && polarised) {
                call.add("L1_MED_DUPLICATION");
            }
            if pa.length < 10.0 && site.tsd_len.is_none_or(|t| t <= 1) && call.three_prime_short && call.structure != "FULL_LENGTH" {
                call.add("EN_INDEPENDENT");
            }
        }
        // ---- beyond poly-A on the 3' side
        let side3 = if strand >= 0 { "LEFT" } else { "RIGHT" };
        let ev3 = ev.junction(side3);
        let (mut bseq, mut bsup) = beyond_polya(&asm);
        if let Some(e3) = ev3.filter(|e| !e.beyond_polya.is_empty()) {
            if bseq.is_empty() {
                bseq = e3.beyond_polya.as_bytes().to_vec();
            }
            bsup = bsup.max(e3.beyond_polya_support);
        }
        // ---- support
        let supported = ev.junctions.iter().filter(|j| j.supported != 0).count() as i64;
        let n_samples = ev.junctions.iter().map(|j| j.n_samples).max().unwrap_or(0);
        let csi = ev.junctions.iter().any(|j| j.cross_sample_identical != 0);
        let fb = foldback(&li, &lf, &ri, &rf);
        let mut td_match = false;
        if let Some(src) = call.source.as_ref().filter(|s| !s.novel) {
            let se = self.lib.source_element(&src.source_id);
            td_match = se.as_deref() == Some(asm.nearest_intact.as_str()) && call.j5_class == "L1";
        }
        // ---- record
        rec.element = call.element.clone();
        rec.structure = call.structure.clone();
        rec.tags = call.tags.clone();
        rec.detail = call.detail.clone();
        rec.consensus = asm.consensus.clone();
        rec.strand = strand;
        if !asm.covered.is_empty() {
            rec.covered_5p = Some(asm.covered_5p());
            rec.covered_3p = Some(asm.covered_3p());
            rec.covered_intervals = asm.covered.clone();
        }
        if asm.nearest_active != "." {
            rec.element_identity = Some(asm.element_identity);
            rec.nearest_active = asm.nearest_active.clone();
        }
        rec.site = (site.contig.clone(), site.l, site.r);
        rec.tsd_seq = site.tsd_seq.clone();
        rec.tsd_len = site.tsd_len;
        rec.en_motif = site.en_motif.clone();
        rec.en_mismatches = site.en_mismatches;
        rec.polya_len = Some(pa.length);
        rec.beyond_polya = String::from_utf8_lossy(&bseq).into_owned();
        rec.beyond_polya_support = bsup;
        if site.slippage {
            rec.detail.set("slippage", site.slippage_detail.as_str());
        }
        if let Some(v) = polya_reads {
            rec.detail.set("polya_reads", v);
        }
        if let Some(v) = polya_unsupported {
            rec.detail.set("polya_unsupported", v);
        }
        if fb != 0 {
            rec.detail.set("foldback_bp", fb);
        }
        let n_short: i64 = ev.junctions.iter().map(|j| j.n_short_used).sum();
        if n_short != 0 {
            rec.detail.set("short_used", n_short);
            rec.detail.set("short_mate_inside", ev.junctions.iter().map(|j| j.n_short_mate_inside).sum::<i64>());
        }
        rec.score_input = Some(ScoreInput {
            element: call.element.clone(),
            structure: call.structure.clone(),
            tags: call.tags.clone(),
            tsd_len: site.tsd_len,
            tsd_verified: site.tsd_verified,
            polya_len: pa.length,
            polya_both_sides: pa.both_sided,
            slippage: site.slippage,
            beyond_polya_len: bseq.len() as i64,
            beyond_polya_support: bsup,
            en_mismatches: site.en_mismatches,
            ends_concordant: !call.j5_class.is_empty() && call.j5_class == call.j3_class,
            td_source_matches_5p: td_match,
            element_identity: asm.element_identity,
            inactive_only: !asm.consensus.is_empty() && !self.lib.is_young(&asm.consensus),
            inv_p1: call.inv_p1,
            junctions_supported: supported,
            n_samples,
            foldback: fb != 0,
            recurrent: false,
            cross_sample_identical: csi,
            novel_tier: call.source.as_ref().map(|s| s.tier.clone()).unwrap_or_default(),
        });
        self.rescore(&mut rec);
        rec
    }

    /// `_score(rec)`
    fn rescore(&self, rec: &mut RteRecord) {
        if let Some(si) = &rec.score_input {
            let (s, p, c) = score(si, &self.cfg.score_weights, &self.cfg.score_thresholds);
            (rec.tprt_score, rec.tprt_points, rec.tprt_call) = (s, p, c);
        }
    }

    /// `_supported_tails(asm, ev_left, ev_right, jl, jr)` -> (A-run at the LEFT insert's REF
    /// edge, T-run at the RIGHT one) reached by >= rte_polya_min_fragments fragments.
    fn supported_tails(&self, asm: &AssemblyResult, ev_left: &[u8], ev_right: &[u8], jl: Option<&JunctionEvidence>, jr: Option<&JunctionEvidence>) -> (f64, f64) {
        let k = self.cfg.polya_min_fragments.max(1) as usize;
        let (eli, _elf, _erf, eri) = split_junction(ev_left, ev_right);
        let (lr, rr) = (edge_run(&eli, true), edge_run(&eri, false));
        let mut sup_a: f64 = if lr.0 == b'A' { lr.1 as f64 } else { 0.0 };
        let mut sup_t: f64 = if rr.0 == b'T' { rr.1 as f64 } else { 0.0 };
        if let Some(j) = jl {
            if j.polya_len_median != 0.0 && lr.0 == b'A' && j.n_independent >= k as i64 {
                sup_a = sup_a.max(j.polya_len_median);
            }
        }
        if let Some(j) = jr {
            if j.polya_len_median != 0.0 && rr.0 == b'T' && j.n_independent >= k as i64 {
                sup_t = sup_t.max(j.polya_len_median);
            }
        }
        // per fragment: the longest A-run | REF (A) / REF | T-run (T)
        let mut per: [FxHashMap<&(String, String), i64>; 2] = [FxHashMap::default(), FxHashMap::default()];
        for lay in &asm.raw_layouts {
            if lay.role == "JUNCTION" || !counts_as_fragment(lay) {
                continue;
            }
            let segs = &lay.segments;
            for (i, sg) in segs.iter().enumerate() {
                if sg.kind != SegKind::PolyA {
                    continue;
                }
                let base = if !sg.target.is_empty() { sg.target.as_str() } else if sg.strand >= 0 { "A" } else { "T" };
                let d = if base == "A" && i + 1 < segs.len() && segs[i + 1].kind == SegKind::Ref {
                    &mut per[0]
                } else if base == "T" && i > 0 && segs[i - 1].kind == SegKind::Ref {
                    &mut per[1]
                } else {
                    continue;
                };
                let e = d.entry(&lay.frag_key).or_insert(0);
                *e = (*e).max(sg.qlen());
            }
        }
        for (b, d) in per.iter().enumerate() {
            let mut lens: Vec<i64> = d.values().copied().collect();
            lens.sort_unstable_by(|a, b| b.cmp(a));
            let kth = if lens.len() >= k { lens[k - 1] as f64 } else { 0.0 };
            if b == 0 {
                sup_a = sup_a.max(kth);
            } else {
                sup_t = sup_t.max(kth);
            }
        }
        (sup_a, sup_t)
    }

    /// `_annotate_key(key, inp)` with the locus's data: pooled with the genotype reads when the
    /// locus has any (then `gt_reads` / `gt_changed` on the record).
    pub fn annotate_key(&self, inp: &InsertionInput, data: LocusData) -> RteRecord {
        let LocusData { evidence, gt_reads } = data;
        if gt_reads.is_empty() {
            return self.annotate(inp, &evidence);
        }
        let mut pooled = InsertionEvidence::new(&inp.title);
        pooled.junctions = evidence.junctions.clone();
        pooled.reads = evidence.reads.iter().cloned().chain(gt_reads).collect::<Vec<EvidenceRead>>();
        let mut rec = self.annotate(inp, &pooled);
        drop(pooled);
        if self.cfg.gt_compare {
            rec.gt_changed = gt_changes(&self.annotate(inp, &evidence), &rec);
        }
        rec
    }

    /// One locus, a panic turned into python's error record; Err only for a failed read-back.
    fn annotate_locus(&self, inp: &InsertionInput, loader: &LocusLoader) -> LocusResult {
        let data = loader.load(&inp.title).map_err(|e| format!("reading the reads of {}: {e}", inp.title))?;
        Ok(match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| self.annotate_key(inp, data))) {
            Ok(r) => r,
            Err(p) => error_record(&inp.title, &panic_text(p.as_ref())),
        })
    }

    /// `annotate_all(inputs)`: per-locus annotation (parallel in bounded chunks, input order) +
    /// the cohort-source and recurrence passes. One record per input, in input order.
    pub fn annotate_all(&mut self, inputs: &[InsertionInput], loader: &LocusLoader) -> Result<Vec<RteRecord>, String> {
        let this = &*self;
        let mut recs: Vec<RteRecord> = map_loci(inputs, self.chunk, |inp| this.annotate_locus(inp, loader)).into_iter().collect::<Result<_, _>>()?;
        self.cohort_source_pass(inputs, loader, &mut recs)?;
        self.recurrence_pass(&mut recs);
        Ok(recs)
    }

    /// `cohort_source_pass`: an unsourced 3' transduction may come from an L1 called elsewhere in
    /// the cohort (only when the discovery genome IS the remap genome).
    pub fn cohort_source_pass(&mut self, inputs: &[InsertionInput], loader: &LocusLoader, recs: &mut [RteRecord]) -> Result<(), String> {
        let (g, r) = (self.cfg.genome_2bit.as_deref().unwrap_or(""), self.cfg.remap_2bit.as_deref().unwrap_or(""));
        if g.is_empty() || g != r || self.novel.locator.is_none() {
            return Ok(());
        }
        let mut l1 = Vec::new();
        for rec in recs.iter() {
            if rec.element == "L1" && (rec.tprt_call == "TPRT" || rec.tprt_call == "LIKELY_TPRT") {
                if let Some(c) = rec.site.0.as_ref().filter(|c| !c.is_empty()) {
                    let pos = if rec.strand >= 0 { rec.site.1 } else { rec.site.2 };
                    if let Some(p) = pos {
                        l1.push((c.clone(), p, if rec.strand >= 0 { '+' } else { '-' }));
                    }
                }
            }
        }
        let redo: Vec<usize> = (0..recs.len())
            .filter(|&k| recs[k].tags.iter().any(|t| t == "TD3P") && !recs[k].tags.iter().any(|t| t.starts_with("TD3P_SOURCE=")))
            .collect();
        if l1.is_empty() || redo.is_empty() {
            return Ok(());
        }
        self.novel.cohort_l1 = l1;
        let this = &*self;
        let again: Vec<Option<RteRecord>> = map_loci(&redo, self.chunk, |&k| {
            let data = loader.load(&inputs[k].title).ok()?;
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| this.annotate_key(&inputs[k], data))).ok()
        });
        for (k, r) in redo.into_iter().zip(again) {
            if let Some(r) = r {
                recs[k] = r; // python: an exception keeps the old record
            }
        }
        Ok(())
    }

    /// `recurrence_pass`: one 5' truncation signature at more than `rte_recurrence_max` loci.
    pub fn recurrence_pass(&self, recs: &mut [RteRecord]) {
        let mx = self.cfg.recurrence_max;
        let mut count: FxHashMap<(String, String, i64), i64> = FxHashMap::default();
        let mut keys: Vec<(usize, (String, String, i64))> = Vec::new();
        for (k, r) in recs.iter().enumerate() {
            let j5 = r.detail.get("j5").or_else(|| r.detail.get("fwd_start"));
            if !matches!(r.element.as_str(), "L1" | "ALU" | "SVA") || r.structure == "FULL_LENGTH" || r.structure == "5P_UNRESOLVED" {
                continue;
            }
            let Some(j5) = j5 else { continue };
            let s = (r.element.clone(), r.structure.clone(), py_int_value(j5).div_euclid(5));
            *count.entry(s.clone()).or_insert(0) += 1;
            keys.push((k, s));
        }
        for (k, s) in keys {
            let n = count[&s];
            if n > mx && recs[k].score_input.is_some() {
                recs[k].score_input.as_mut().unwrap().recurrent = true;
                recs[k].detail.set("recurrence", n);
                self.rescore(&mut recs[k]);
            }
        }
    }
}

/// python `int(v)` of a detail value (int / float truncation / numeric string).
fn py_int_value(v: &crate::record::DetailValue) -> i64 {
    use crate::record::DetailValue::*;
    match v {
        Int(i) => *i,
        Float(f) => crate::pyfmt::py_int(*f),
        Str(s) => s.trim().parse::<i64>().unwrap_or_else(|_| panic!("ValueError: invalid literal for int() with base 10: '{s}'")),
    }
}

/// `_beyond_polya(asm)`: sequence beyond the 3' poly-A (element sense: ...X | POLYA | REF) and
/// the number of independent fragments carrying >= 10 bp of it.
fn beyond_polya(asm: &AssemblyResult) -> (Vec<u8>, i64) {
    let mut best: Vec<u8> = Vec::new();
    let mut frags: Vec<&(String, String)> = Vec::new();
    for lay in &asm.layouts {
        if !counts_as_fragment(lay) {
            continue;
        }
        let segs = &lay.segments;
        let n = segs.len();
        if n < 3 || segs[n - 1].kind != SegKind::Ref || segs[n - 2].kind != SegKind::PolyA {
            continue;
        }
        let x = &segs[n - 3];
        if matches!(x.kind, SegKind::PolyA | SegKind::Ref | SegKind::Local) || x.qlen() < 10 {
            continue;
        }
        let s = lay.piece(x).to_vec();
        if lay.role == "JUNCTION" {
            if s.len() > best.len() {
                best = s;
            }
        } else {
            if !frags.contains(&&lay.frag_key) {
                frags.push(&lay.frag_key);
            }
            if best.is_empty() {
                best = s;
            }
        }
    }
    (best, frags.len() as i64)
}

/// The `_premrna_fn(site)` closure state: the gene-model window around the site, indexed once
/// (lazily, at the first query) with minimap2 `sr`.
struct Premrna<'a> {
    genome: &'a dyn Genome,
    gm: &'a GeneModel,
    contig: String,
    lo: i64,
    anchor: i64,
    win: i64,
    al: OnceLock<Option<Aligner>>,
}

impl Premrna<'_> {
    fn call(&self, seq: &[u8]) -> Option<String> {
        let al = self
            .al
            .get_or_init(|| {
                let r = self.genome.fetch(&self.contig, self.lo, self.anchor + self.win);
                if r.is_empty() {
                    None
                } else {
                    Aligner::from_seq(&r, MapOpts::SR)
                }
            })
            .as_ref()?;
        for h in al.map(seq) {
            if (h.mlen as f64) / (h.blen.max(1) as f64) < 0.9 {
                continue;
            }
            let p = self.lo + (h.r_st + h.r_en).div_euclid(2);
            if (p - self.anchor).abs() <= 250 {
                continue; // a local template, not a distal pre-mRNA
            }
            for g in self.gm.candidates(&self.contig, p).unwrap_or_default() {
                if g.start <= p && p < g.end {
                    let feat = self.gm.genic_feature(&g.exons, &g.strand, p).1;
                    return Some(format!("{}:{}@{}:{}", g.name, feat, self.contig, p));
                }
            }
        }
        None
    }
}

/// The text of a caught panic.
fn panic_text(p: &(dyn std::any::Any + Send)) -> String {
    if let Some(s) = p.downcast_ref::<&str>() {
        s.to_string()
    } else if let Some(s) = p.downcast_ref::<String>() {
        s.clone()
    } else {
        "panic".to_string()
    }
}

/// python's error record: `RteRecord(key)` + `detail["error"] = f"{type(e).__name__}: {e}"[:200]`.
/// Panics that carry a python exception text (`ValueError: ...`) keep it; others are
/// `RustPanic: <message>`.
pub fn error_record(key: &str, msg: &str) -> RteRecord {
    let head = msg.split(':').next().unwrap_or("");
    let pyish = !head.is_empty() && head.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') && head.ends_with("Error") && msg.contains(": ");
    let full = if pyish { msg.to_string() } else { format!("RustPanic: {msg}") };
    let mut rec = RteRecord::new(key);
    rec.detail.set("error", full.chars().take(200).collect::<String>().as_str());
    rec
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers_like_python() {
        assert_eq!(richer_junction(b"", b"aaGG"), b"aaGG");
        assert_eq!(richer_junction(b"aGG", b"aaGG"), b"aaGG");
        assert_eq!(richer_junction(b"aaGG", b"ccGG"), b"aaGG");
        assert_eq!(richer_junction(b"aaGG", b""), b"aaGG");
        assert_eq!(
            default_sidecars("/x/P1.combined.txt.gz"),
            ("/x/P1.insertions.evidence.tsv.gz".to_string(), "/x/P1.insertions.reads.fa.gz".to_string())
        );
        assert_eq!(default_gt_reads("P1.txt.gz"), "P1.insertions.genotype_reads.fa.gz");
        assert_eq!(title_gap("chr1:100-115"), Some(15));
        assert_eq!(title_gap("chr1:polyA_3"), None);
    }

    #[test]
    fn error_record_like_python() {
        let r = error_record("x", "ValueError: max() iterable argument is empty");
        assert_eq!(r.detail.get("error").unwrap().py_str(), "ValueError: max() iterable argument is empty");
        let r = error_record("x", "index out of bounds: the len is 3 but the index is 7");
        assert_eq!(r.detail.get("error").unwrap().py_str(), "RustPanic: index out of bounds: the len is 3 but the index is 7");
        let r = error_record("x", &"KeyError: 'a'".repeat(40));
        assert_eq!(r.detail.get("error").unwrap().py_str().chars().count(), 200);
        assert_eq!(r.row()[16], format!("error={}", &"KeyError: 'a'".repeat(40)[..200]));
    }
}
