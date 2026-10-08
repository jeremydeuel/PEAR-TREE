//! tools/rte/structure.py -- element / 5' structure / tags from the element-sense layouts.
//!
//! STATUS: types + callback interfaces FOUNDATION; algorithms ported by WP-STRUCT (golden
//! `classify`: 109 pytest + 336 e2e events identical).
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

use crate::assembly::{counts_as_fragment, AssemblyResult, ReadLayout, SegKind, Segment, SiteContext};
use crate::config::StructureCfg;
use crate::library::{Ordered, RteLibrary};
use crate::record::Detail;
use crate::sequtil::{edlib_best, low_complexity_default};
use crate::transduction::{known_source, SourceCall};

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

// ------------------------------------------------------------------------------ small helpers

/// python `s[st:en]` (negative indices count from the end, clamped, empty when st >= en).
fn py_slice(s: &[u8], st: i64, en: i64) -> &[u8] {
    let n = s.len() as i64;
    let norm = |i: i64| if i < 0 { (i + n).max(0) } else { i.min(n) };
    let (a, b) = (norm(st), norm(en));
    if a >= b {
        &[]
    } else {
        &s[a as usize..b as usize]
    }
}

/// python `needle in hay` for strings.
fn contains(hay: &[u8], needle: &[u8]) -> bool {
    needle.is_empty() || hay.windows(needle.len()).any(|w| w == needle)
}

/// `lib.cons_class.get(target) == cls_`
fn class_is(lib: &RteLibrary, target: &str, cls: &str) -> bool {
    lib.class_of(target) == Some(cls)
}

/// python `Counter` + `max(c, key=c.get)`: first maximal key in insertion order.
fn first_max<K: Clone>(items: &[(K, i64)]) -> Option<K> {
    let mut best: Option<&(K, i64)> = None;
    for e in items {
        if best.is_none_or(|b| e.1 > b.1) {
            best = Some(e);
        }
    }
    best.map(|b| b.0.clone())
}

fn bump<K: PartialEq>(items: &mut Vec<(K, i64)>, k: K, v: i64) {
    match items.iter_mut().find(|(x, _)| *x == k) {
        Some(e) => e.1 += v,
        None => items.push((k, v)),
    }
}

/// python set.add on an insertion-ordered Vec.
fn add_unique<T: PartialEq>(v: &mut Vec<T>, x: T) {
    if !v.contains(&x) {
        v.push(x);
    }
}

/// `_first_element_after_ref(lay)`: (element segment, [segments between REF and it]).
fn first_element_after_ref(lay: &ReadLayout) -> (Option<&Segment>, Vec<&Segment>) {
    let mut extras = Vec::new();
    for s in lay.segments.iter().skip(1) {
        if s.kind == SegKind::Element {
            return (Some(s), extras);
        }
        extras.push(s);
    }
    (None, extras)
}

/// `_last_element_before_ref(lay)` (only the element is used by classify).
fn last_element_before_ref(lay: &ReadLayout) -> Option<&Segment> {
    let n = lay.segments.len();
    lay.segments[..n.saturating_sub(1)].iter().rev().find(|s| s.kind == SegKind::Element)
}

// ------------------------------------------------------------------------------ classify

/// `classify(res, lib, ctx, cfg, novel_finder, premrna, pseudogene, legacy_class)`.
#[allow(clippy::too_many_arguments)]
pub fn classify(
    res: &AssemblyResult,
    lib: &RteLibrary,
    ctx: Option<&SiteContext>,
    cfg: &StructureCfg,
    novel_finder: Option<&dyn SourceFinder>,
    premrna: Option<PremrnaFn>,
    pseudogene: Option<&PseudogeneArg>,
    legacy_class: Option<&str>,
) -> StructureCall {
    let c = cfg;
    let mut call = StructureCall::default();
    let cls: &str = &res.element_class;
    let has_cls = !cls.is_empty();
    let cons: &str = &res.consensus;
    let cend: i64 = if cons.is_empty() { 0 } else { lib.cons_end.get(cons).map_or(0, |&v| v as i64) };
    let layouts = &res.layouts;

    // indices into `layouts`; junction consensus first, then reads (stable)
    let mut five: Vec<usize> =
        (0..layouts.len()).filter(|&k| layouts[k].ref_left() && layouts[k].segments.len() > 1).collect();
    let mut three: Vec<usize> =
        (0..layouts.len()).filter(|&k| layouts[k].ref_right() && layouts[k].segments.len() > 1).collect();
    five.sort_by_key(|&k| layouts[k].role != "JUNCTION");
    three.sort_by_key(|&k| layouts[k].role != "JUNCTION");

    // ---------------------------------------------------------------- 5' junction
    let mut j5: Option<&Segment> = None;
    let mut j5_extras: Vec<&Segment> = Vec::new();
    let mut votes: Vec<(bool, i64)> = Vec::new();
    let mut cands: Vec<(&Segment, Vec<&Segment>)> = Vec::new();
    for &k in &five {
        let lay = &layouts[k];
        let (seg, extras) = first_element_after_ref(lay);
        let Some(seg) = seg else { continue };
        if !ref_at_breakpoint(&lay.segments[0], ctx, 30) {
            continue; // REF | element read from elsewhere in the window (a reference copy)
        }
        if lib.class_of(&seg.target) != Some(cls) && has_cls {
            cands.push((seg, extras)); // chimera candidate: no vote
            continue;
        }
        cands.push((seg, extras));
        bump(&mut votes, seg.strand > 0, if lay.role == "JUNCTION" { 3 } else { 1 });
    }
    if cands.is_empty() {
        for &k in &five {
            let lay = &layouts[k];
            let extras: Vec<&Segment> = lay.segments[1..]
                .iter()
                .filter(|x| matches!(x.kind, SegKind::Flank5p | SegKind::Local | SegKind::Flank3p))
                .collect();
            if extras.is_empty() {
                continue;
            }
            if let Some(seg) = chain(layouts, extras[0], true, cls, lib) {
                bump(&mut votes, seg.strand > 0, 1);
                cands.push((seg, extras));
                break;
            }
        }
    }
    if !cands.is_empty() {
        let pick: Vec<&(&Segment, Vec<&Segment>)> = match first_max(&votes) {
            Some(want) => cands
                .iter()
                .filter(|x| (x.0.strand > 0) == want && (!has_cls || class_is(lib, &x.0.target, cls)))
                .collect(),
            None => cands.iter().collect(),
        };
        let chosen = pick.first().copied().unwrap_or(&cands[0]);
        j5 = Some(chosen.0);
        j5_extras = chosen.1.clone();
        call.j5_class = lib.class_of(&chosen.0.target).unwrap_or("").to_string();
    }

    // ---------------------------------------------------------------- 3' junction
    let mut j3: Option<&Segment> = three.iter().find_map(|&k| last_element_before_ref(&layouts[k]));
    if j3.is_none() && !three.is_empty() {
        let segs = &layouts[three[0]].segments;
        if let Some(fl) = segs[..segs.len() - 1].iter().find(|x| x.kind == SegKind::Flank3p) {
            j3 = chain(layouts, fl, false, cls, lib);
        }
    }
    call.j3_class = j3.map_or(String::new(), |s| lib.class_of(&s.target).unwrap_or("").to_string());
    call.has_polya_3p = three.iter().any(|&k| {
        let s = &layouts[k].segments;
        s.len() >= 2 && s[s.len() - 2].kind == SegKind::PolyA
    });
    if let Some(j3s) = j3 {
        if has_cls && class_is(lib, &j3s.target, cls) {
            call.three_prime_truncated = j3s.strand > 0 && j3s.t_en < cend - c.truncated_3p_tolerance;
            call.three_prime_short = j3s.strand > 0 && j3s.t_en < cend - c.en_independent_3p_tolerance;
            call.detail.set("j3", j3s.t_en);
        }
    }

    // ---------------------------------------------------------------- structure
    if let Some(j5s) = j5.filter(|s| has_cls && class_is(lib, &s.target, cls)) {
        let tol = c.full_length_tol(cls);
        if j5s.strand > 0 {
            call.j5_pos = Some(j5s.t_st);
            let in_hex = cls == "SVA" && lib.landmark_at(&j5s.target, j5s.t_st).to_lowercase() == "hexamer";
            call.structure = if j5s.t_st <= tol || in_hex { "FULL_LENGTH" } else { "TRUNCATED_5P" }.into();
            call.detail.set("j5", j5s.t_st);
            if let Some((a, b)) = sense_switch(res, j5s, lib, cls, c) {
                call.structure = "INVERTED_5P_SWITCH".into();
                call.detail.set("switch", format!("{a}>{b}"));
            }
        } else {
            call.j5_pos = Some(j5s.t_en);
            call.structure = "INVERTED_5P".into();
            inversion(&mut call, res, j5s, lib, c);
        }
        let lm = lib.landmark_at(cons, call.j5_pos.unwrap_or(-1));
        if lm != "." {
            call.detail.set("j5_feature", lm);
        }
    } else if has_cls {
        call.structure = "5P_UNRESOLVED".into();
        if let Some(j3s) = j3 {
            if j3s.strand > 0
                && call.has_polya_3p
                && j3s.t_en >= cend - c.terminal_inv_tolerance
                && inverted_tail_5p(&five, layouts, c)
            {
                // minimal twin priming: REF | poly-T at the 5' junction, terminus + poly-A at 3'
                call.structure = "INVERTED_5P".into();
                call.j5_pos = Some(cend);
                call.detail.set("inv", "polyA");
                call.detail.set("inv_junction", "unresolved");
            }
        }
    }
    if cls == "SVA" && call.structure == "5P_UNRESOLVED" && td5_at_junction(&five, layouts) {
        call.structure = "FULL_LENGTH".into();
        call.detail.set("j5", "td5p");
    }

    // ---------------------------------------------------------------- segments by kind
    let mut flank3_any = false;
    let mut flank5: Vec<&Segment> = Vec::new();
    let mut unknown: Vec<(usize, &Segment)> = Vec::new();
    for (k, lay) in layouts.iter().enumerate() {
        for s in &lay.segments {
            match s.kind {
                SegKind::Flank3p => flank3_any = true,
                SegKind::Flank5p => flank5.push(s),
                SegKind::Unknown if s.qlen() >= c.min_unknown_bp => unknown.push((k, s)),
                _ => {}
            }
        }
    }

    // ---------------------------------------------------------------- 3' transduction
    let side3 = if res.strand >= 0 { "LEFT" } else { "RIGHT" };
    let mut flank3_ok: Vec<&Segment> = Vec::new();
    let mut side3_sources: Vec<String> = Vec::new();
    let mut long_sources: Vec<String> = Vec::new();
    for lay in layouts {
        for (i, s) in lay.segments.iter().enumerate() {
            if s.kind != SegKind::Flank3p || element_after(lay, i, cls, lib) {
                continue;
            }
            if s.qlen() < 60 && explained_by_consensus(py_slice(&lay.seq, s.q_st, s.q_en), lib, 0.15) {
                continue; // an element end that also sits (as a repeat) in some flank
            }
            if flank_masked_frac(&lib.flanks3, s) >= c.td_max_masked_frac {
                continue; // the hit lies in a soft-masked repeat of the flank
            }
            let sid = lib.source_for_flank(&s.target);
            if !source_class_ok(lib, &sid, cls) {
                continue;
            }
            flank3_ok.push(s);
            let nxt = lay.segments.get(i + 1);
            let tail_pos = nxt.is_some_and(|n| n.kind == SegKind::PolyA && n.q_st - s.q_en <= 3);
            let short_ok = tail_pos
                && s.identity >= c.td_short_flank_identity
                && !at_rich(
                    py_slice(lib.flanks3.get(&s.target).map_or(&[][..], |v| v.as_slice()), s.t_st, s.t_en),
                    c.td_short_flank_max_at,
                );
            if s.qlen() >= c.td_min_flank_bp || short_ok {
                add_unique(&mut long_sources, sid.clone());
            }
            if lay.side == side3 || lay.side.is_empty() || j5.is_none() {
                add_unique(&mut side3_sources, sid);
            }
        }
    }
    flank3_ok.retain(|s| {
        let sid = lib.source_for_flank(&s.target);
        side3_sources.contains(&sid) && long_sources.contains(&sid)
    });
    let mut src = if flank3_ok.is_empty() { None } else { known_source(&flank3_ok, lib) };
    let mut td_seq = unexplained_tail(layouts, cls, lib, c, &mut call);
    if src.is_none() && !td_seq.is_empty() {
        if let Some(f) = novel_finder {
            src = f.find(&td_seq);
        }
    }
    let mut premrna_lab = String::new();
    if src.is_none() && !td_seq.is_empty() && td_seq.len() as i64 >= c.premrna_min_bp {
        if let Some(p) = premrna {
            // host-gene sequence between element and poly-A: pre-mRNA co-insertion
            premrna_lab = p(&td_seq).unwrap_or_default();
            if !premrna_lab.is_empty() {
                call.add("PREMRNA_COINSERT");
                call.detail.set("premrna", premrna_lab.as_str());
                td_seq = Vec::new();
            }
        }
    }
    let has_td = src.is_some() || (!td_seq.is_empty() && has_cls);
    if has_td && (has_cls || src.is_some()) {
        call.add("TD3P");
    }
    if let Some(sc) = &src {
        call.source = Some(sc.clone());
        call.add(&format!("TD3P_SOURCE={}", sc.source_id));
        call.detail.set("td_end", sc.td_end);
        if sc.novel {
            call.add("NOVEL_SOURCE");
            call.detail.set("source_identity", sc.identity);
            if !sc.tier.is_empty() {
                call.detail.set("novel_tier", sc.tier.as_str());
            }
        } else if !sc.detail.is_empty() {
            call.detail.set("td_flank", sc.detail.as_str());
        }
    }
    call.td3p_seq = td_seq.clone();

    // ---------------------------------------------------------------- SVA 5' transduction
    if cls == "SVA" {
        let f5: Vec<&Segment> = flank5
            .iter()
            .copied()
            .filter(|s| s.strand > 0 && flank_masked_frac(&lib.flanks5, s) < c.td_max_masked_frac)
            .collect();
        // `e in f5` is dataclass EQUALITY (Segment: PartialEq)
        if !f5.is_empty()
            || j5_extras.iter().any(|e| e.kind == SegKind::Unknown && e.qlen() >= 30)
            || j5_extras.iter().any(|e| e.kind == SegKind::Flank5p && e.qlen() >= 30 && f5.contains(e))
        {
            call.add("TD5P");
            if !f5.is_empty() {
                let mut best: Vec<(String, i64)> = Vec::new();
                for s5 in &f5 {
                    bump(&mut best, lib.source_for_flank(&s5.target), s5.matches);
                }
                if let Some(sid) = first_max(&best) {
                    call.detail.set("td5_source", sid.as_str());
                    call.add(&format!("TD5P_SOURCE={sid}"));
                }
            }
        }
    }

    // ---------------------------------------------------------------- templated / pre-mRNA
    let fb_segs = foldback_5p(&five, layouts, c, 12, 15);
    if !fb_segs.is_empty() {
        call.add("FOLDBACK_INVDUP_5P");
        let mut keys: Vec<&(String, String)> = Vec::new();
        for (fk, _) in &fb_segs {
            add_unique(&mut keys, *fk);
        }
        call.detail.set("foldback_frags", keys.len() as i64);
    }
    // python excludes by object identity -> (layout index, segment index)
    let exclude: Vec<(usize, usize)> = fb_segs.iter().map(|(_, id)| *id).collect();
    let tmpl = local_templates(layouts, ctx, c, &exclude);
    if let Some((d, n, _iv)) = tmpl.templated {
        call.add("TEMPLATED_LOCAL");
        if let Some(d) = d {
            call.detail.set("templated_dist", d);
        }
        call.detail.set("templated_frags", n);
    }
    if let Some((d, _n, iv)) = tmpl.distal {
        if premrna_lab.is_empty() {
            // a distal local template = co-inserted local pre-mRNA (d and iv are always set here)
            let contig = ctx.and_then(|x| x.contig.clone()).unwrap_or_else(|| "None".into());
            let dd = d.map_or("None".to_string(), |x| x.to_string());
            let (a, b) = iv.unwrap_or((0, 0));
            premrna_lab = format!("local:{contig}:{a}-{b}(d={dd})");
            call.add("PREMRNA_COINSERT");
            call.detail.set("premrna", premrna_lab.as_str());
        }
    }
    if let Some(p) = premrna {
        if premrna_lab.is_empty() {
            for (k, s) in &unknown {
                if s.qlen() < c.premrna_min_bp {
                    continue;
                }
                let seq = py_slice(&layouts[*k].seq, s.q_st, s.q_en);
                if !td_seq.is_empty() && contains(&td_seq, seq) {
                    continue;
                }
                if let Some(lab) = p(seq).filter(|l| !l.is_empty()) {
                    call.add("PREMRNA_COINSERT");
                    call.detail.set("premrna", lab);
                    break;
                }
            }
        }
    }

    // ---------------------------------------------------------------- chimeric ends
    // ALU at the 5' end of an SVA = its Alu-like domain: compatible, not a chimera
    let compatible = call.j5_class == "ALU" && call.j3_class == "SVA";
    if !call.j5_class.is_empty() && !call.j3_class.is_empty() && call.j5_class != call.j3_class && !compatible {
        call.add("CHIMERIC_ENDS");
        let ends = format!("{}/{}", call.j5_class, call.j3_class);
        call.detail.set("ends", ends);
    }

    // ---------------------------------------------------------------- element
    let (pg_genes, pg_hits, pg_structure): (&[String], &[(String, String)], Option<PgStructureFn>) = match pseudogene {
        Some(p) => (p.genes, p.hits, p.structure_fn),
        None => (&[], &[], None),
    };
    let nonpolya_unknown: i64 = unknown.iter().map(|(_, s)| s.qlen()).sum();
    if !pg_genes.is_empty() && !pg_hits.is_empty() {
        // an exon-exon junction read of a candidate gene proves a spliced mRNA
        call.element = "PSEUDOGENE".into();
        call.add("EXON_JUNCTION");
        call.detail.set("exon_junction", pg_hits[0].0.as_str());
        if has_cls {
            call.detail.set("rte_in_mrna", cls);
        }
        const DROP: [&str; 9] = [
            "TD3P",
            "TD3P_SOURCE",
            "NOVEL_SOURCE",
            "TD5P",
            "TD5P_SOURCE",
            "CHIMERIC_ENDS",
            "FOLDBACK_INVDUP_5P",
            "PREMRNA_COINSERT",
            "TEMPLATED_LOCAL",
        ];
        call.tags.retain(|t| !DROP.contains(&t.split('=').next().unwrap_or("")));
        call.detail.pop("premrna");
        call.source = None;
        let st = pg_structure.and_then(|f| f(layouts)).filter(|s| !s.is_empty());
        call.structure = st.unwrap_or_else(|| "5P_UNRESOLVED".into());
        return call;
    }
    if has_cls {
        call.element = match cls {
            "L1" | "ALU" | "SVA" => cls,
            _ => "UNKNOWN",
        }
        .into();
    } else if !pg_genes.is_empty() {
        call.add("PSEUDOGENE_CANDIDATE");
    }
    if call.element == "UNKNOWN" && !has_cls {
        if src.is_some() && call.has_polya_3p {
            call.element = "ORPHAN_TD".into();
        } else if call.has_polya_3p && nonpolya_unknown < c.min_unknown_bp && !flank3_any && pg_genes.is_empty() {
            call.element = "POLYA_ONLY".into();
        } else if matches!(
            legacy_class,
            Some(
                "non_RTE_SV"
                    | "SV_DELETION"
                    | "SV_DUPLICATION"
                    | "SV_INVERSION"
                    | "microsatellite"
                    | "templated_insertion"
            )
        ) {
            call.element = "NON_TPRT".into();
        }
    }
    if matches!(call.element.as_str(), "POLYA_ONLY" | "ORPHAN_TD" | "PSEUDOGENE" | "NON_TPRT" | "UNKNOWN") && !has_cls {
        call.structure = "5P_UNRESOLVED".into();
    }
    call
}

/// `_chain(layouts, piece, after, cls_, lib)`: a segment of the same kind/target/strand as
/// `piece` joined (<= 10 bp, skipping short POLYA) to an element segment of class `cls`.
fn chain<'a>(layouts: &'a [ReadLayout], piece: &Segment, after: bool, cls: &str, lib: &RteLibrary) -> Option<&'a Segment> {
    for lay in layouts {
        let segs = &lay.segments;
        let n = segs.len() as i64;
        for (i, s) in segs.iter().enumerate() {
            if s.kind != piece.kind || s.target != piece.target || s.strand != piece.strand {
                continue;
            }
            if s.kind == SegKind::Local && (s.t_en < piece.t_st - 20 || s.t_st > piece.t_en + 20) {
                continue; // a different local template
            }
            let d: i64 = if after { 1 } else { -1 };
            let mut j = i as i64 + d;
            while 0 <= j && j < n && segs[j as usize].kind == SegKind::PolyA && segs[j as usize].qlen() <= 30 {
                j += d;
            }
            if !(0 <= j && j < n) {
                continue;
            }
            let e = &segs[j as usize];
            let near = &segs[(j - d) as usize];
            let gap = if !after { near.q_st - e.q_en } else { e.q_st - near.q_en };
            if e.kind == SegKind::Element && gap <= 10 && (cls.is_empty() || class_is(lib, &e.target, cls)) {
                return Some(e);
            }
        }
    }
    None
}

/// `_inversion(call, res, j5, lib, c)`: twin-priming geometry from switch points.
fn inversion(call: &mut StructureCall, res: &AssemblyResult, j5: &Segment, lib: &RteLibrary, c: &StructureCfg) {
    let cls: &str = &res.element_class;
    let p2 = j5.t_en;
    let mut switches: Vec<i64> = Vec::new();
    let mut inner: Option<(i64, i64)> = None;
    for lay in &res.layouts {
        let el: Vec<&Segment> =
            lay.segments.iter().filter(|s| s.kind == SegKind::Element && class_is(lib, &s.target, cls)).collect();
        for w in el.windows(2) {
            let (a, b) = (w[0], w[1]);
            if (a.strand > 0) == (b.strand > 0) {
                continue;
            }
            if b.q_st - a.q_en > 10 {
                continue; // not adjacent in the read
            }
            if a.strand < 0 && b.strand > 0 {
                switches.push(b.t_st);
                if inner.is_none() {
                    inner = Some((a.t_st, b.t_st));
                }
            } else {
                switches.push(a.t_en);
            }
        }
    }
    let sense_min = res.segments_on_cons.iter().filter(|x| x.2).map(|x| x.0).min();
    let p1 = inner.map(|iv| iv.1).or(sense_min);
    let mut pts = switches;
    pts.sort();
    pts.dedup();
    let mut clusters: Vec<i64> = Vec::new();
    for p in pts {
        if clusters.last().is_none_or(|&l| p - l > c.switch_cluster_bp) {
            clusters.push(p);
        }
    }
    if clusters.len() >= 2 {
        call.structure = "INVERTED_5P_SWITCH".into();
    }
    let anti_lo = res.segments_on_cons.iter().filter(|x| !x.2).map(|x| x.0).min();
    let inv = format!("{}-{p2}", anti_lo.map_or("?".to_string(), |v| v.to_string()));
    call.detail.set("inv", inv);
    if let Some(p1) = p1 {
        call.inv_p1 = Some(p1);
        call.detail.set("fwd_start", p1);
        let d = p1 - p2;
        let j = match d.cmp(&0) {
            std::cmp::Ordering::Greater => format!("del{d}"),
            std::cmp::Ordering::Less => format!("dup{}", -d),
            std::cmp::Ordering::Equal => "blunt".to_string(),
        };
        call.detail.set("inv_junction", j);
        call.detail.set("inv_exact", inner.is_some() as i64);
    } else {
        call.detail.set("inv_junction", "unresolved");
    }
    if let Some((a, b)) = inner {
        if (a - b).abs() <= c.foldback_tolerance {
            call.add("FOLDBACK_INVDUP_5P");
        }
    }
}

/// `_element_after(lay, i, cls_, lib)`
fn element_after(lay: &ReadLayout, i: usize, cls: &str, lib: &RteLibrary) -> bool {
    lay.segments[i + 1..]
        .iter()
        .any(|t| t.kind == SegKind::Element && (cls.is_empty() || class_is(lib, &t.target, cls)))
}

/// `_unexplained_tail(layouts, cls_, lib, c, call)`: GT_* reads never count (counts_as_fragment).
fn unexplained_tail(
    layouts: &[ReadLayout],
    cls: &str,
    lib: &RteLibrary,
    c: &StructureCfg,
    call: &mut StructureCall,
) -> Vec<u8> {
    let mut frags: Vec<&(String, String)> = Vec::new();
    let mut best: &[u8] = &[];
    let mut best_j: &[u8] = &[];
    for lay in layouts {
        if !counts_as_fragment(lay) {
            continue;
        }
        let segs = &lay.segments;
        for i in 0..segs.len().saturating_sub(1) {
            let (x, nxt) = (&segs[i], &segs[i + 1]);
            if x.kind != SegKind::Unknown || nxt.kind != SegKind::PolyA || nxt.q_st - x.q_en > 5 {
                continue;
            }
            if let Some(after) = segs.get(i + 2) {
                if after.kind != SegKind::Ref && !(after.kind == SegKind::Unknown && i + 3 == segs.len()) {
                    continue; // an A-run inside genomic sequence, not the poly-A tail
                }
            }
            if x.qlen() < c.td_min_bp {
                continue;
            }
            let piece = py_slice(&lay.seq, x.q_st, x.q_en);
            if low_complexity_default(piece) {
                continue;
            }
            if let Some(prev) = segs[..i].iter().rev().find(|t| t.kind == SegKind::Element) {
                if !cls.is_empty() && !class_is(lib, &prev.target, cls) {
                    continue;
                }
            }
            add_unique(&mut frags, &lay.frag_key);
            if lay.role == "JUNCTION" {
                if piece.len() > best_j.len() {
                    best_j = piece;
                }
            } else if piece.len() > best.len() {
                best = piece;
            }
        }
    }
    if !frags.is_empty() {
        call.detail.set("td_frags", frags.len() as i64);
    }
    if (frags.len() as i64) < c.td_min_fragments {
        return Vec::new();
    }
    if !best_j.is_empty() { best_j } else { best }.to_vec()
}

/// `_template_abs(seg, ctx)`
fn template_abs(seg: &Segment, ctx: Option<&SiteContext>) -> Option<(i64, i64)> {
    if seg.t_st < 0 {
        return None;
    }
    if seg.target == "wide" {
        return Some((seg.t_st, seg.t_en));
    }
    let ctx = ctx?;
    if ctx.window_seq.is_empty() {
        return None;
    }
    Some((ctx.window_start + seg.t_st, ctx.window_start + seg.t_en))
}

type FragKey = (String, String);

/// `_foldback_5p(five_layouts, c, max_gap, min_len)` -> [(frag_key, (layout index, segment index))]
fn foldback_5p<'a>(
    five: &[usize],
    layouts: &'a [ReadLayout],
    c: &StructureCfg,
    max_gap: i64,
    min_len: i64,
) -> Vec<(&'a FragKey, (usize, usize))> {
    let mut hits = Vec::new();
    for &k in five {
        let lay = &layouts[k];
        if !counts_as_fragment(lay) {
            continue;
        }
        let segs = &lay.segments;
        if segs.len() < 2 || segs[0].kind != SegKind::Ref || segs[0].strand == 0 || segs[0].t_st < 0 {
            continue;
        }
        let r = &segs[0];
        let x = &segs[1];
        if x.kind != SegKind::Local || x.target != "site" || x.qlen() < min_len || x.q_st - r.q_en > max_gap {
            continue;
        }
        if x.strand == r.strand || x.t_st < 0 {
            continue;
        }
        let ok = if r.strand > 0 {
            let j = r.t_en;
            j - max_gap <= x.t_en && x.t_en <= j + 2
        } else {
            let j = r.t_st;
            j - 2 <= x.t_st && x.t_st <= j + max_gap
        };
        if ok {
            hits.push((&lay.frag_key, (k, 1usize)));
        }
    }
    let mut keys: Vec<&FragKey> = Vec::new();
    for (fk, _) in &hits {
        add_unique(&mut keys, *fk);
    }
    if (keys.len() as i64) < c.templated_min_fragments {
        return Vec::new();
    }
    hits
}

/// (distance or None, n_fragments, interval or None)
type Template = (Option<i64>, i64, Option<(i64, i64)>);

/// (template interval, distance, fragment)
type TemplateCand<'a> = (Option<(i64, i64)>, Option<i64>, &'a FragKey);

#[derive(Default)]
struct Templates {
    templated: Option<Template>,
    distal: Option<Template>,
}

/// `_local_templates(layouts, ctx, c, exclude)`; `exclude` = (layout index, segment index)
/// (python: object identity).
fn local_templates(
    layouts: &[ReadLayout],
    ctx: Option<&SiteContext>,
    c: &StructureCfg,
    exclude: &[(usize, usize)],
) -> Templates {
    let bps: Vec<i64> = ctx.map_or(Vec::new(), |x| [x.left_bp, x.right_bp].into_iter().flatten().collect());
    let mut cands: Vec<TemplateCand> = Vec::new();
    for (k, lay) in layouts.iter().enumerate() {
        if !counts_as_fragment(lay) {
            continue;
        }
        for (si, s) in lay.segments.iter().enumerate() {
            if s.kind != SegKind::Local || s.qlen() < c.templated_min_bp || exclude.contains(&(k, si)) {
                continue;
            }
            if s.identity != 0.0 && s.identity < c.templated_min_identity {
                continue;
            }
            if low_complexity_default(py_slice(&lay.seq, s.q_st, s.q_en)) {
                continue;
            }
            let iv = template_abs(s, ctx);
            let mut d = None;
            if let (Some(iv), Some(&lo), Some(&hi)) = (iv, bps.iter().min(), bps.iter().max()) {
                if iv.0 >= lo - 5 && iv.1 <= hi + 5 {
                    continue; // the TSD itself (a read through the duplication)
                }
                d = bps
                    .iter()
                    .map(|&b| if iv.0 <= b && b <= iv.1 { 0 } else { (iv.0 - b).abs().min((iv.1 - b).abs()) })
                    .min();
            }
            cands.push((iv, d, &lay.frag_key));
        }
    }
    let mut out = Templates::default();
    for (iv, d, _) in &cands {
        let mut fs: Vec<&FragKey> = Vec::new();
        for (v, _, f) in &cands {
            let hit = match (iv, v) {
                (None, None) => true,
                (Some(iv), Some(v)) => v.0 < iv.1 + 20 && iv.0 < v.1 + 20,
                _ => false,
            };
            if hit {
                add_unique(&mut fs, *f);
            }
        }
        let n = fs.len() as i64;
        if n < c.templated_min_fragments {
            continue;
        }
        let slot = if d.is_none_or(|d| d <= c.templated_max_dist) { &mut out.templated } else { &mut out.distal };
        if slot.is_none_or(|o| n > o.1) {
            *slot = Some((*d, n, *iv));
        }
    }
    out
}

/// `_sense_switch(res, j5, lib, cls_, c)` -> (sense end, anti end) on the consensus
fn sense_switch(res: &AssemblyResult, j5: &Segment, lib: &RteLibrary, cls: &str, c: &StructureCfg) -> Option<(i64, i64)> {
    let m = c.switch_min_seg;
    for lay in &res.layouts {
        let el: Vec<&Segment> =
            lay.segments.iter().filter(|s| s.kind == SegKind::Element && class_is(lib, &s.target, cls)).collect();
        for w in el.windows(2) {
            let (a, b) = (w[0], w[1]);
            if !(a.strand > 0 && b.strand < 0) || b.q_st - a.q_en > 10 {
                continue;
            }
            if a.qlen() < m || b.qlen() < m {
                continue;
            }
            if a.t_st >= j5.t_st - 50 && b.t_en > a.t_en + 20 {
                return Some((a.t_en, b.t_en));
            }
        }
    }
    None
}

/// `_inverted_tail_5p(five_layouts, c)`: REF | poly-T with no element segment after it.
fn inverted_tail_5p(five: &[usize], layouts: &[ReadLayout], c: &StructureCfg) -> bool {
    five.iter().any(|&k| {
        let segs = &layouts[k].segments;
        if segs.len() < 2 || segs[0].kind != SegKind::Ref {
            return false;
        }
        let t = &segs[1];
        t.kind == SegKind::PolyA
            && t.strand < 0
            && t.qlen() >= c.inverted_tail_min
            && t.q_st - segs[0].q_en <= c.inverted_tail_max_gap
            && !segs[2..].iter().any(|x| x.kind == SegKind::Element)
    })
}

/// `_flank_masked_frac(lib, seg, flanks)`: soft-masked fraction of the hit's flank interval.
fn flank_masked_frac(flanks: &Ordered<Vec<u8>>, seg: &Segment) -> f64 {
    let fs = flanks.get(&seg.target).map_or(&[][..], |v| v.as_slice());
    let piece = py_slice(fs, seg.t_st.max(0), seg.t_en.max(0));
    if piece.is_empty() {
        return 0.0;
    }
    piece.iter().filter(|c| c.is_ascii_lowercase()).count() as f64 / piece.len() as f64
}

/// `_source_class_ok(lib, sid, cls_)`: `src.get("element_class") or src.get("class") or ""`.
fn source_class_ok(lib: &RteLibrary, sid: &str, cls: &str) -> bool {
    if cls.is_empty() {
        return true;
    }
    let sc = match lib.sources.get(sid) {
        Some(r) => [r.get("element_class"), r.get("class")].into_iter().flatten().find(|v| !v.is_empty()).unwrap_or(""),
        None => "",
    };
    sc.is_empty() || sc == cls
}

/// `_at_rich(piece, max_frac)`
fn at_rich(piece: &[u8], max_frac: f64) -> bool {
    if piece.is_empty() {
        return false;
    }
    let a = piece.iter().filter(|c| c.eq_ignore_ascii_case(&b'A')).count();
    let t = piece.iter().filter(|c| c.eq_ignore_ascii_case(&b'T')).count();
    a.max(t) as f64 / piece.len() as f64 > max_frac
}

/// `_td5_at_junction(five_layouts)`: REF | (POLYA) | sense FLANK5P.
fn td5_at_junction(five: &[usize], layouts: &[ReadLayout]) -> bool {
    five.iter().any(|&k| {
        let segs = &layouts[k].segments;
        segs.len() >= 2
            && segs[0].kind == SegKind::Ref
            && segs[1..]
                .iter()
                .find(|x| x.kind != SegKind::PolyA)
                .is_some_and(|n| n.kind == SegKind::Flank5p && n.strand > 0)
    })
}

/// `_explained_by_consensus(piece, lib, max_frac)`
fn explained_by_consensus(piece: &[u8], lib: &RteLibrary, max_frac: f64) -> bool {
    lib.consensus.iter().any(|(name, seq)| {
        let end = lib.cons_end.get(name).map_or(seq.len() as i64, |&e| e as i64);
        edlib_best(piece, py_slice(seq, 0, end), max_frac, true).is_some()
    })
}

/// `_ref_at_breakpoint(ref, ctx, tol)`: the REF piece ends at an insertion breakpoint (or carries
/// no window coordinates).
fn ref_at_breakpoint(r: &Segment, ctx: Option<&SiteContext>, tol: i64) -> bool {
    let Some(ctx) = ctx else { return true };
    if r.t_st < 0 {
        return true;
    }
    let (bps, a, b): (Vec<i64>, i64, i64) = if ctx.window_seq.is_empty() {
        // genome-free local reference = right_flank + N*30 + left_flank
        let n = ctx.right_flank.len() as i64;
        let mut v = Vec::new();
        if !ctx.right_flank.is_empty() {
            v.push(n);
        }
        if !ctx.left_flank.is_empty() {
            v.push(n + 30);
        }
        (v, r.t_st, r.t_en)
    } else {
        (
            [ctx.left_bp, ctx.right_bp].into_iter().flatten().collect(),
            ctx.window_start + r.t_st,
            ctx.window_start + r.t_en,
        )
    };
    if bps.is_empty() {
        return true;
    }
    bps.iter().any(|&x| (a - x).abs().min((b - x).abs()) <= tol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn py_slice_like_python() {
        let s = b"ABCDEF";
        assert_eq!(py_slice(s, 1, 3), b"BC");
        assert_eq!(py_slice(s, -2, 10), b"EF");
        assert_eq!(py_slice(s, 4, 2), b"");
        assert_eq!(py_slice(s, -1, -1), b"");
        assert!(contains(b"ACGT", b"CG"));
        assert!(!contains(b"ACGT", b"GC"));
    }
}
