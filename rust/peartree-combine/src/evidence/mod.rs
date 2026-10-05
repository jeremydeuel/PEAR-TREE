//! Pooled per-patient junction evidence (TPRT sidecars). OWNER: P5 (this file, store.rs,
//! output.rs); row.rs / dedup.rs / junction.rs: P3; filters.rs: P4.
//!
//! Mirrors `apply_evidence`, `_make_evaluate`, `_judge_chunk`, `EvidencePool.absorb_one_sided`,
//! `_reanchor`, `_member_loci`, `_absorb_target`, `_replace_clips` of
//! src/combine_insertions_evidence.py:766-1243, 1255-1298. SPEC.md §4 (overview §4.0).
//!
//! The pooled independent-fragment gate (`require_independent_fragments`, fail reason
//! `n_independent<k`) is NOT ported (SPEC.md §0); `supported` is still reported. The far-pair
//! split condition is `pside is not None and far_pair_split and prec` (python with the gate off).
//!
//! Memory / parallelism (PLAN.md): the insertion list is split into the store's chunks
//! (consecutive index ranges); each chunk is judged on a rayon pool of `ctx.threads` workers
//! with only its own rows in memory (`Store::chunk_rows`), its insertions mutated in place
//! (far-pair splits, clip replacement), and its records' reads moved to the store
//! (`put_reads` / `finish_chunk`) before the rows are dropped. Records keep counters +
//! consensus only. Chunk results are merged in insertion order, so outputs do not depend on
//! chunking or thread count. `absorb_one_sided` uses a per-(contig, side) sorted position index.

pub mod dedup;
pub mod filters;
pub mod junction;
pub mod output;
pub mod row;
pub mod store;

use crate::context::Ctx;
use crate::evidence::filters::{aligned_part, far_pair_check, slippage_check, to_one_sided, Breakpoints};
use crate::evidence::junction::{evaluate_junction, JunctionRecord, Supported};
use crate::evidence::row::{EvidenceRow, Role};
use crate::evidence::store::Store;
use crate::genome::RefFetch;
use crate::library::{LibraryMatcher, Matcher};
use crate::model::{ContigId, FileId, Insertion, Interner, Member, Side, SideSet, SIDES};
use crate::seq::QualSeq;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::BTreeMap;
use std::path::Path;

/// State kept from apply_evidence for absorb_one_sided and write_evidence_outputs (python
/// `(kept, records, failed, {"pool", "store", ...})`).
pub struct EvidenceState {
    /// python `records`: insertion name -> uid whose record list is current for that name
    /// (setdefault for kept then failed; absorb overwrites; absorbed names popped)
    pub records: FxHashMap<String, u32>,
    /// python `recmap` (`id(ins)` -> [JunctionRecord]), indexed by uid
    pub recmap: Vec<Option<Vec<JunctionRecord>>>,
    /// python `failed` names, sorted (byte order) and unique -- written after the survivors
    pub failed: Vec<String>,
    /// python `by_id` (`_member_loci(ins)` at apply time; absorb extends it), indexed by uid
    pub members: Vec<Vec<Member>>,
    /// sidecar rows + detached reads
    pub store: Store,
    /// basenames with a sidecar (python `have`)
    pub have: Vec<FileId>,
}

/// Per-chunk result (python `_judge_chunk` return value).
struct ChunkOut {
    /// (kept, records) per insertion of the chunk, in order
    out: Vec<(bool, Vec<JunctionRecord>)>,
    reasons: Vec<String>,
    tprt: Vec<String>,
    short: Vec<(&'static str, u32)>,
    replaced: usize,
}

/// Read-only state shared by the chunk workers (python `_CTX`).
struct Judge<'a> {
    ctx: &'a Ctx,
    store: &'a Store,
    /// member loci by position in the apply_evidence input
    members: &'a [Vec<Member>],
    /// have[FileId]
    have: &'a [bool],
    matcher: Option<&'a dyn Matcher>,
    bps: Option<&'a Breakpoints>,
}

/// `apply_evidence(insertions, accepted_files, cfg, breakpoints, threads, shard_dir)`.
///
/// Returns `(insertions, None)` unchanged when no accepted file has a sidecar (python returns
/// None -> legacy path). Otherwise `(kept, Some(state))`: kept insertions in input order (with
/// far-pair splits / clip replacements applied), failed (far_pair / slippage) recorded in
/// `state.failed`. Chunks are evaluated in parallel (rayon, `ctx.threads`), results merged in
/// insertion order -- output identical for any thread count / chunking. SPEC.md §4.0.
pub fn apply_evidence(
    insertions: Vec<Insertion>,
    accepted: &[FileId],
    ctx: &Ctx,
    breakpoints: Option<&Breakpoints>,
    shard_dir: &Path,
) -> Result<(Vec<Insertion>, Option<EvidenceState>), String> {
    let cfg = &ctx.cfg;
    if cfg.require_independent_fragments {
        return Err("require_independent_fragments is not supported by the Rust port (SPEC.md §0)".into());
    }
    // python `have` (computed before creating the shard dir: python returns None -- leaving an
    // empty shard dir -- only when no sidecar exists, which the driver already excludes)
    if !accepted.iter().any(|&f| store::sidecar_path(&ctx.files[f as usize].path).exists()) {
        return Ok((insertions, None));
    }
    let min_ind = cfg.min_independent_fragments;
    let members: Vec<Vec<Member>> = insertions.iter().map(member_loci).collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(ctx.threads.max(1))
        .build()
        .map_err(|e| format!("cannot build thread pool: {e}"))?;
    let (store, chunks) = pool.install(|| Store::build(shard_dir, accepted, &insertions, &members, ctx))?;
    let mut have_mask = vec![false; ctx.files.len()];
    for &f in store.have() {
        have_mask[f as usize] = true;
    }
    let missing: Vec<&str> = accepted
        .iter()
        .filter(|&&f| !have_mask[f as usize])
        .map(|&f| ctx.files[f as usize].basename.as_str())
        .collect();
    if !missing.is_empty() {
        println!(
            "WARNING: {} discovery file(s) have no evidence sidecar; insertions they contribute to have supported=NA: {}{}",
            missing.len(),
            missing[..missing.len().min(5)].join(","),
            if missing.len() > 5 { " ..." } else { "" }
        );
    }
    let slip_on = cfg.slippage_reject;
    let far_on = cfg.far_pair_strict;
    let matcher: Option<LibraryMatcher> =
        if slip_on || far_on { Some(LibraryMatcher::open(&cfg.resolve_rte_library())?) } else { None };

    // split the insertion list into the chunks (consecutive ranges tiling 0..n)
    let mut insertions = insertions;
    let results: Vec<ChunkOut> = {
        let judge = Judge {
            ctx,
            store: &store,
            members: &members,
            have: &have_mask,
            matcher: matcher.as_ref().map(|m| m as &dyn Matcher),
            bps: breakpoints,
        };
        let mut slices: Vec<(usize, usize, &mut [Insertion])> = Vec::with_capacity(chunks.len());
        let mut rest: &mut [Insertion] = &mut insertions;
        for (c, &(lo, hi)) in chunks.iter().enumerate() {
            let (a, b) = std::mem::take(&mut rest).split_at_mut(hi - lo);
            slices.push((c, lo, a));
            rest = b;
        }
        debug_assert!(rest.is_empty());
        pool.install(|| {
            slices
                .into_par_iter()
                .map(|(c, lo, sl)| judge.judge_chunk(c, lo, sl))
                .collect::<Result<Vec<_>, String>>()
        })?
    };

    // merge in insertion order (python apply_evidence after _run_chunks)
    let n_uid = insertions.iter().map(|i| i.uid as usize + 1).max().unwrap_or(0);
    let mut recmap: Vec<Option<Vec<JunctionRecord>>> = (0..n_uid).map(|_| None).collect();
    let mut members_by_uid: Vec<Vec<Member>> = (0..n_uid).map(|_| Vec::new()).collect();
    let mut reasons: BTreeMap<String, u64> = BTreeMap::new();
    let mut tprt: BTreeMap<String, u64> = BTreeMap::new();
    let mut short: BTreeMap<&'static str, u64> = BTreeMap::new();
    let mut n_replaced = 0usize;
    let mut oks: Vec<bool> = Vec::with_capacity(insertions.len());
    for res in results {
        for r in res.reasons {
            *reasons.entry(r).or_default() += 1;
        }
        for r in res.tprt {
            *tprt.entry(r).or_default() += 1;
        }
        for (k, v) in res.short {
            *short.entry(k).or_default() += v as u64;
        }
        n_replaced += res.replaced;
        for (ok, recs) in res.out {
            let k = oks.len();
            recmap[insertions[k].uid as usize] = Some(recs);
            oks.push(ok);
        }
    }
    for (ins, ms) in insertions.iter().zip(members) {
        members_by_uid[ins.uid as usize] = ms;
    }
    let mut kept: Vec<Insertion> = Vec::with_capacity(insertions.len());
    let mut failed_objs: Vec<Insertion> = Vec::new();
    for (ins, ok) in insertions.into_iter().zip(oks) {
        if ok {
            kept.push(ins);
        } else {
            failed_objs.push(ins);
        }
    }
    let mut records: FxHashMap<String, u32> = FxHashMap::default();
    for i in &kept {
        records.entry(i.name(&ctx.contigs)).or_insert(i.uid);
    }
    let kept_names: FxHashSet<String> = records.keys().cloned().collect();
    let mut failed: Vec<String> = Vec::new();
    for i in &failed_objs {
        let name = i.name(&ctx.contigs);
        if !kept_names.contains(&name) {
            records.entry(name.clone()).or_insert(i.uid);
            failed.push(name);
        }
    }
    failed.sort();
    failed.dedup();

    if !tprt.is_empty() {
        let s: Vec<String> = tprt.iter().map(|(k, v)| format!("{k}={v}")).collect();
        println!("TPRT combine filters: {}", s.join(", "));
    }
    println!(
        "evidence sidecars: {}/{} files; {} insertions evaluated in {} chunk(s), {} evidence rows",
        store.have().len(),
        accepted.len(),
        kept.len() + failed_objs.len(),
        chunks.len(),
        store.n_rows()
    );
    let n_sup: usize = kept
        .iter()
        .map(|i| recmap[i.uid as usize].as_ref().map_or(0, |r| r.iter().filter(|r| r.supported == Supported::No).count()))
        .sum();
    println!("junctions below {min_ind} pooled independent fragments (reported as supported=0, not dropped): {n_sup}");
    if cfg.indel_aware_consensus {
        println!("indel-aware consensus replaced {n_replaced} junction clip(s) in combined output");
    }
    if cfg.count_short_overhang {
        let n_used: u64 = records
            .values()
            .filter_map(|&u| recmap[u as usize].as_ref())
            .flat_map(|r| r.iter())
            .map(|r| r.n_short_used as u64)
            .sum();
        let s: Vec<String> = short.iter().map(|(k, v)| format!("{k}={v}")).collect();
        println!(
            "SHORT overhang reads: {n_used} fragment(s) used, {} rejected ({})",
            short.values().sum::<u64>(),
            if s.is_empty() { "-".to_string() } else { s.join(", ") }
        );
    }
    let _ = reasons; // python logs these only with the (dropped) gate
    let have = store.have().to_vec();
    Ok((kept, Some(EvidenceState { records, recmap, failed, members: members_by_uid, store, have })))
}

impl Judge<'_> {
    /// python `_judge_chunk(c)` over `slice` = insertions[lo..lo+len].
    fn judge_chunk(&self, c: usize, lo: usize, slice: &mut [Insertion]) -> Result<ChunkOut, String> {
        let ctx = self.ctx;
        let cfg = &ctx.cfg;
        let min_ind = cfg.min_independent_fragments;
        let rows = self.store.chunk_rows(c, &ctx.contigs)?;
        let mut res = ChunkOut { out: Vec::with_capacity(slice.len()), reasons: vec![], tprt: vec![], short: vec![], replaced: 0 };
        let genome: &dyn RefFetch = &ctx.genome;
        for (j, ins) in slice.iter_mut().enumerate() {
            let ms = &self.members[lo + j];
            let mut ok = true;
            let open = ins.open_side_eff();
            let mut recs: Vec<JunctionRecord> = Vec::with_capacity(2);
            for side in SIDES {
                if Some(side) == open {
                    continue;
                }
                let rec = evaluate(ins, ms, side, ctx, &mut |m, s, out| out.extend_from_slice(rows.get(m, s)));
                res.short.extend(rec.short_reasons.iter().copied());
                recs.push(rec);
            }
            if cfg.far_pair_strict && open.is_none() && recs.len() == 2 {
                let matcher = self.matcher.expect("matcher built when far_pair_strict");
                if let Some((reason, pside)) = far_pair_check(ins, &recs, cfg, matcher, self.bps, Some(genome), &ctx.contigs) {
                    res.tprt.push(format!("far_pair:{reason}"));
                    let prec = pside.and_then(|p| recs.iter().position(|r| r.side == p));
                    match (pside, prec) {
                        (Some(pside), Some(pi)) if cfg.far_pair_split => {
                            let old_name = ins.name(&ctx.contigs);
                            to_one_sided(ins, pside);
                            ins.member_sides = Some(ms.iter().map(|m| (*m, SideSet::only(pside))).collect());
                            let mut rec = recs.swap_remove(pi);
                            rec.insertion_id = ins.name(&ctx.contigs);
                            let mut loci = locus_names(ms, &ctx.contigs);
                            loci.push(old_name);
                            loci.sort();
                            loci.dedup();
                            rec.member_loci = loci.join(",");
                            rec.fail_reason = format!("split_from_far_pair:{reason}");
                            recs = vec![rec];
                            res.tprt.push("far_pair:split_to_one_sided".into());
                        }
                        _ => {
                            for r in recs.iter_mut() {
                                r.fail_reason = format!("far_pair:{reason}");
                            }
                            res.reasons.push(format!("far_pair:{reason}"));
                            ok = false;
                        }
                    }
                }
            }
            if ok && cfg.slippage_reject {
                let matcher = self.matcher.expect("matcher built when slippage_reject");
                let why = slippage_check(ins, &recs, cfg, matcher, genome, &ctx.contigs);
                if !why.is_empty() {
                    res.tprt.push(why.split(':').next().unwrap_or("").to_string());
                    res.reasons.push(why.split('(').next().unwrap_or("").to_string());
                    for r in recs.iter_mut() {
                        r.fail_reason = why.clone();
                    }
                    ok = false;
                }
            }
            if ok && ins.files.iter().all(|&f| self.have[f as usize]) {
                for r in recs.iter_mut() {
                    r.supported = if (r.n_independent as usize) < min_ind { Supported::No } else { Supported::Yes };
                }
            }
            if ok && cfg.indel_aware_consensus {
                res.replaced += replace_clips(ins, &recs);
            }
            // store.detach: reads rendered with rec.insertion_id, rows dropped
            for rec in recs.iter_mut() {
                let text = rec.render_reads(&rec.insertion_id, &ctx.files);
                rec.reads_ref = Some(self.store.put_reads(c, &text));
                rec.rows = None;
            }
            res.out.push((ok, recs));
        }
        drop(rows);
        self.store.finish_chunk(c)?;
        Ok(res)
    }
}

/// python `evaluate(ins, side)` (`_make_evaluate`): pooled rows of every member allowed on
/// `side` (member order, each member's rows in sidecar order), re-anchored, evaluated;
/// `aligned` and `member_loci` set. `lookup(m, side, out)` appends m's rows.
fn evaluate(
    ins: &Insertion,
    ms: &[Member],
    side: Side,
    ctx: &Ctx,
    lookup: &mut dyn FnMut(&Member, Side, &mut Vec<EvidenceRow>),
) -> JunctionRecord {
    let mut pooled = Vec::new();
    let sides_of = member_sides_index(ins);
    for m in ms {
        let allowed = match &sides_of {
            Some(ix) => ix.get(m).copied().unwrap_or(SideSet::BOTH),
            None => ins.allowed_sides(m),
        };
        if allowed.contains(side) {
            lookup(m, side, &mut pooled);
        }
    }
    let pooled = reanchor(pooled, side, ins.junction(side));
    let name = ins.name(&ctx.contigs);
    let genome: &dyn RefFetch = &ctx.genome;
    let mut rec = evaluate_junction(&name, side, pooled, &ctx.cfg, Some(genome), &ctx.contigs, &ctx.files);
    rec.aligned = aligned_part(ins, side);
    let mut loci = locus_names(ms, &ctx.contigs);
    loci.sort();
    loci.dedup();
    if !(loci.len() == 1 && loci[0] == name) {
        rec.member_loci = loci.join(",");
    }
    rec
}

/// hash index of a large `member_sides` (python dict lookup); None = use the linear
/// `Insertion::allowed_sides` (small or absent dicts)
fn member_sides_index(ins: &Insertion) -> Option<FxHashMap<Member, SideSet>> {
    match &ins.member_sides {
        Some(v) if v.len() > 16 => {
            let mut ix = FxHashMap::default();
            for (m, s) in v {
                ix.insert(*m, *s); // dict semantics: last assignment wins (keys are unique anyway)
            }
            Some(ix)
        }
        _ => None,
    }
}

/// locus NAME strings of the members (unsorted, with repeats)
fn locus_names(ms: &[Member], contigs: &Interner) -> Vec<String> {
    ms.iter().map(|(_, l)| l.name(contigs)).collect()
}

/// python `_absorb_target` over a static position index: for each (contig, real side) the
/// insertions whose open side is not that side, sorted by (junction, index). The candidate
/// key `(t is one-sided, |d|, t NAME, input index)` is total, so the result equals python's
/// scan (first minimum in list order) for any scan order.
struct AbsorbIndex {
    by: FxHashMap<(ContigId, Side), Vec<(i64, usize)>>,
}

impl AbsorbIndex {
    fn new(insertions: &[Insertion]) -> AbsorbIndex {
        let mut by: FxHashMap<(ContigId, Side), Vec<(i64, usize)>> = FxHashMap::default();
        for (i, t) in insertions.iter().enumerate() {
            for side in SIDES {
                if t.open_side_eff() == Some(side) {
                    continue;
                }
                // python `_ins_junction(t, side) - pos` would raise on None (poly-A types never
                // reach this stage); such candidates are skipped
                if let Some(p) = t.junction(side) {
                    by.entry((t.contig, side)).or_default().push((p, i));
                }
            }
        }
        for v in by.values_mut() {
            v.sort_unstable();
        }
        AbsorbIndex { by }
    }

    #[allow(clippy::too_many_arguments)]
    fn target(&self, x: usize, contig: ContigId, side: Side, pos: i64, tol: i64, alive: &[bool], one_sided: &[bool], names: &[String]) -> Option<usize> {
        let v = self.by.get(&(contig, side))?;
        let start = v.partition_point(|&(p, _)| p < pos - tol);
        let mut best: Option<(bool, i64, &str, usize)> = None;
        for &(p, i) in &v[start..] {
            if p > pos + tol {
                break;
            }
            if i == x || !alive[i] {
                continue;
            }
            let key = (one_sided[i], (p - pos).abs(), names[i].as_str(), i);
            if best.is_none_or(|b| key < b) {
                best = Some(key);
            }
        }
        best.map(|b| b.3)
    }
}

impl EvidenceState {
    /// `EvidencePool.absorb_one_sided(insertions)` (evidence.py:1193), run by the driver after
    /// the clipped-remap filter. Returns (surviving insertions, n_absorbed). SPEC.md §4.7.
    /// Must not be O(one-sided x all): index candidates per (contig) by junction position; the
    /// selection key `(target is one-sided, |d|, target NAME string)` makes the result
    /// independent of scan order.
    pub fn absorb_one_sided(&mut self, insertions: Vec<Insertion>, ctx: &Ctx) -> (Vec<Insertion>, usize) {
        let cfg = &ctx.cfg;
        let mut tol = cfg.merge_tolerance_bp;
        if cfg.far_pair_strict {
            tol = tol.max(5);
        }
        if tol <= 0 {
            return (insertions, 0);
        }
        let min_ind = cfg.min_independent_fragments;
        let mut insertions = insertions;
        let n = insertions.len();
        let names: Vec<String> = insertions.iter().map(|i| i.name(&ctx.contigs)).collect();
        let one_sided: Vec<bool> = insertions.iter().map(|i| i.open_side_eff().is_some()).collect();
        let mut alive = vec![true; n];
        let mut one: Vec<usize> = (0..n).filter(|&i| one_sided[i]).collect();
        // python sorted(key=(-len(files), name)) -- stable
        one.sort_by(|&a, &b| insertions[b].files.len().cmp(&insertions[a].files.len()).then_with(|| names[a].cmp(&names[b])));
        let index = AbsorbIndex::new(&insertions);
        let mut n_abs = 0usize;
        for xi in one {
            if !alive[xi] {
                continue;
            }
            let x = &insertions[xi];
            let side = x.open_side_eff().expect("one-sided").other();
            let Some(pos) = x.junction(side) else { continue };
            let Some(ti) = index.target(xi, x.contig, side, pos, tol, &alive, &one_sided, &names) else { continue };
            let xc: Option<Box<[u8]>> = x.clipped(side).map(|q| q.seq.clone());
            let tc: Option<Box<[u8]>> = insertions[ti].clipped(side).map(|q| q.seq.clone());
            if let (Some(xc), Some(tc)) = (&xc, &tc) {
                if !crate::seq::clips_agree(&[&tc[..], &xc[..]], 0.6, 8, 6) {
                    continue;
                }
            }
            let x_files = x.files.clone();
            let xu = x.uid as usize;
            let tu = insertions[ti].uid as usize;
            let xm = self.members[xu].clone();
            // members[t] = dict.fromkeys(members[t] + members[x])
            let old: FxHashSet<Member> = self.members[tu].iter().copied().collect();
            let mut seen = old.clone();
            for m in &xm {
                if seen.insert(*m) {
                    self.members[tu].push(*m);
                }
            }
            let t = &mut insertions[ti];
            // ms = dict(member_sides or {}); ms[m] = (side,) for x's members not previously in t
            let mut ms: Vec<(Member, SideSet)> = t.member_sides.take().unwrap_or_default();
            let mut pos_of: FxHashMap<Member, usize> = ms.iter().enumerate().map(|(i, (m, _))| (*m, i)).collect();
            for m in &xm {
                if !old.contains(m) {
                    match pos_of.get(m) {
                        Some(&i) => ms[i].1 = SideSet::only(side),
                        None => {
                            pos_of.insert(*m, ms.len());
                            ms.push((*m, SideSet::only(side)));
                        }
                    }
                }
            }
            t.member_sides = Some(ms);
            // t.files + [f for f in x.files if f not in t.files] (against the ORIGINAL t.files)
            let orig: FxHashSet<FileId> = t.files.iter().copied().collect();
            t.files.extend(x_files.iter().copied().filter(|f| !orig.contains(f)));
            let store = &self.store;
            let new = evaluate(t, &self.members[tu], side, ctx, &mut |m, s, out| out.extend(store.lookup(m, s, &ctx.contigs)));
            let mut recs = self.recmap[tu].take().unwrap_or_default();
            // python replaces the record of `side` (one per side) with `new`; `new` itself
            // still feeds _replace_clips when t had no such record (not reachable: t's `side`
            // is a real junction)
            let mut unplaced = None;
            let slot = recs.iter().position(|r| r.side == side);
            match slot {
                Some(i) => recs[i] = new,
                None => unplaced = Some(new),
            }
            for r in recs.iter_mut() {
                r.supported = if (r.n_independent as usize) >= min_ind { Supported::Yes } else { Supported::No };
            }
            if cfg.indel_aware_consensus {
                let nr = match (slot, &unplaced) {
                    (Some(i), _) => &recs[i],
                    (None, Some(u)) => u,
                    (None, None) => unreachable!(),
                };
                replace_clips(t, std::slice::from_ref(nr));
            }
            self.recmap[tu] = Some(recs);
            self.records.insert(names[ti].clone(), tu as u32);
            alive[xi] = false;
            n_abs += 1;
        }
        let alive_names: FxHashSet<&str> = (0..n).filter(|&i| alive[i]).map(|i| names[i].as_str()).collect();
        for i in 0..n {
            if !alive[i] && !alive_names.contains(names[i].as_str()) {
                self.records.remove(&names[i]);
            }
        }
        let out: Vec<Insertion> = insertions.into_iter().zip(alive).filter_map(|(i, a)| a.then_some(i)).collect();
        (out, n_abs)
    }
}

/// `_member_loci(ins)`: member_loci + (f, ins.name) for f in files, deduplicated in order.
pub fn member_loci(ins: &Insertion) -> Vec<Member> {
    let own = ins.locus();
    let mut seen: FxHashSet<Member> = FxHashSet::default();
    let mut out = Vec::with_capacity(ins.member_loci.len() + ins.files.len());
    for m in ins.member_loci.iter().copied().chain(ins.files.iter().map(|&f| (f, own))) {
        if seen.insert(m) {
            out.push(m);
        }
    }
    out
}

/// `_reanchor(rows, side, junction)` (evidence.py:1257): CLIP/SHORT rows with clip_at >= 0 whose
/// own locus junction (`LocusKey::junction(side)`) differs from `junction` get
/// `clip_at += junction - jm` when the result stays within `0..=len(seq)`; others unchanged.
/// `junction` None -> unchanged.
pub fn reanchor(rows: Vec<EvidenceRow>, side: Side, junction: Option<i64>) -> Vec<EvidenceRow> {
    let Some(junction) = junction else { return rows };
    let mut rows = rows;
    for r in rows.iter_mut() {
        if matches!(r.role, Role::Clip | Role::Short) && r.clip_at >= 0 {
            let jm = r.locus.junction(side);
            if jm != junction {
                let at = r.clip_at + (junction - jm);
                if 0 <= at && at <= r.seq.len() as i64 {
                    r.clip_at = at;
                }
            }
        }
    }
    rows
}

/// python `str.islower()` on ASCII: some cased character and no uppercase one.
fn py_islower(s: &[u8]) -> bool {
    s.iter().any(|b| b.is_ascii_lowercase()) && !s.iter().any(|b| b.is_ascii_uppercase())
}

/// `_replace_clips(ins, recs)` -> number replaced: for each rec with non-empty
/// combined_consensus.seq and an existing clip of rec.side: skip if shorter than the old clip or
/// equal ignoring case; else the clip becomes QualSeq(seq lowercased if the old clip
/// `is_lower()` else uppercased, combined_consensus.score).
pub fn replace_clips(ins: &mut Insertion, recs: &[JunctionRecord]) -> usize {
    let mut n = 0;
    for rec in recs {
        let c = &rec.combined_consensus;
        if c.seq.is_empty() {
            continue;
        }
        let slot = match rec.side {
            Side::Right => &mut ins.right_clipped,
            Side::Left => &mut ins.left_clipped,
        };
        let Some(old) = slot.as_ref() else { continue };
        if c.seq.len() < old.len() || c.seq.eq_ignore_ascii_case(&old.seq) {
            continue;
        }
        let seq = if py_islower(&old.seq) { c.seq.to_ascii_lowercase() } else { c.seq.to_ascii_uppercase() };
        *slot = Some(QualSeq::new(seq, c.score.clone()));
        n += 1;
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::ConsensusResult;
    use crate::model::{InsType, LocusKey, Tok};

    fn ins(contigs: &Interner, uid: u32, name: &str, ty: InsType, open: Option<Side>, files: Vec<FileId>) -> Insertion {
        let k = LocusKey::parse(name, contigs).unwrap();
        Insertion {
            uid,
            contig: k.contig,
            name_start: k.start,
            name_end: k.end,
            ty,
            open_side: open,
            left_clipped: None,
            left_aligned: None,
            left_pos: Some(k.start.pos),
            right_clipped: None,
            right_aligned: None,
            right_pos: Some(k.end.pos),
            member_loci: files.iter().map(|&f| (f, k)).collect(),
            files,
            member_sides: None,
        }
    }

    fn row(contigs: &Interner, locus: &str, role: Role, clip_at: i64, seq: &str) -> EvidenceRow {
        let star = contigs.intern("*");
        EvidenceRow {
            file: 0,
            locus: LocusKey::parse(locus, contigs).unwrap(),
            side: Side::Right,
            role,
            frag: "f".into(),
            r12: 1,
            flag: 0,
            ref_: star,
            pos: -1,
            strand: "*".into(),
            outer: -1,
            mref: star,
            mpos: -1,
            mstrand: "*".into(),
            tlen: 0,
            mapq: 0,
            cigar: "*".into(),
            clip_at,
            seq: seq.as_bytes().into(),
            qual: Box::new([]),
        }
    }

    fn rec(side: Side, seq: &str, n_ind: u32) -> JunctionRecord {
        JunctionRecord {
            insertion_id: "x".into(),
            side,
            rows: Some(vec![]),
            reads_ref: None,
            n_reads: 0,
            n_fragments: 0,
            n_independent: n_ind,
            n_samples: 0,
            n_mates: 0,
            n_duplicates: 0,
            n_cross: 0,
            n_dup_coord: 0,
            n_dup_seq: 0,
            member_loci: String::new(),
            n_short_used: 0,
            n_short_rejected: 0,
            n_short_mate_inside: 0,
            n_independent_no_short: 0,
            short_reasons: vec![],
            supported: Supported::Na,
            consensus: ConsensusResult::default(),
            combined_consensus: ConsensusResult { seq: seq.as_bytes().to_vec(), score: vec![40; seq.len()], ..Default::default() },
            polya_end: false,
            fail_reason: String::new(),
            aligned: vec![],
        }
    }

    #[test]
    fn member_loci_dedup_in_order() {
        let c = Interner::new();
        let mut i = ins(&c, 0, "chr1:100-120", InsType::FullInfo, None, vec![0, 1]);
        let other = LocusKey::parse("chr1:101-120", &c).unwrap();
        i.member_loci = vec![(1, other), (0, i.locus()), (1, other)];
        i.files = vec![0, 1, 0];
        let own = i.locus();
        assert_eq!(member_loci(&i), vec![(1, other), (0, own), (1, own)]);
    }

    #[test]
    fn reanchor_shifts_only_other_loci_clip_rows_in_range() {
        let c = Interner::new();
        let rows = vec![
            row(&c, "chr1:100-120", Role::Clip, 5, "ACGTACGTAC"), // own locus: unchanged
            row(&c, "chr1:100-117", Role::Clip, 5, "ACGTACGTAC"), // jm 117 -> +3 -> 8
            row(&c, "chr1:100-117", Role::Short, 8, "ACGTACGTAC"), // 11 > len 10 -> unchanged
            row(&c, "chr1:100-oneside_123", Role::Short, 5, "ACGTACGTAC"), // -3 -> 2
            row(&c, "chr1:100-117", Role::Mate, 5, "ACGTACGTAC"), // not CLIP/SHORT
            row(&c, "chr1:100-117", Role::Clip, -1, "ACGTACGTAC"), // clip_at < 0
            row(&c, "chr1:100-110", Role::Clip, 0, "ACGTACGTAC"), // +10 -> 10 == len, allowed
        ];
        let got: Vec<i64> = reanchor(rows.clone(), Side::Right, Some(120)).iter().map(|r| r.clip_at).collect();
        assert_eq!(got, vec![5, 8, 8, 2, 5, -1, 10]);
        let same: Vec<i64> = reanchor(rows.clone(), Side::Right, None).iter().map(|r| r.clip_at).collect();
        assert_eq!(same, vec![5, 5, 8, 5, 5, -1, 0]);
        // LEFT uses the start token (100 everywhere) -> unchanged at junction 100
        let left: Vec<i64> = reanchor(rows, Side::Left, Some(100)).iter().map(|r| r.clip_at).collect();
        assert_eq!(left, vec![5, 5, 8, 5, 5, -1, 0]);
    }

    #[test]
    fn replace_clips_rules() {
        let c = Interner::new();
        let mut i = ins(&c, 0, "chr1:100-120", InsType::FullInfo, None, vec![0]);
        i.right_clipped = Some(QualSeq::new(b"acgtac".to_vec(), vec![30; 6]));
        i.left_clipped = Some(QualSeq::new(b"ACGTACGT".to_vec(), vec![30; 8]));
        // right: longer, case-different -> replaced, lowercased (old is lower)
        // left: shorter -> skipped
        let n = replace_clips(&mut i, &[rec(Side::Right, "ACGTACGG", 2), rec(Side::Left, "ACG", 2)]);
        assert_eq!(n, 1);
        assert_eq!(&*i.right_clipped.as_ref().unwrap().seq, b"acgtacgg");
        assert_eq!(&*i.right_clipped.as_ref().unwrap().qual, &[40; 8]);
        // equal ignoring case -> skipped
        assert_eq!(replace_clips(&mut i, &[rec(Side::Left, "acgtacgt", 2)]), 0);
        // longer, old upper -> uppercased
        assert_eq!(replace_clips(&mut i, &[rec(Side::Left, "acgtacgtt", 2)]), 1);
        assert_eq!(&*i.left_clipped.as_ref().unwrap().seq, b"ACGTACGTT");
        // empty consensus / missing clip -> skipped
        i.right_clipped = None;
        assert_eq!(replace_clips(&mut i, &[rec(Side::Right, "ACGTACGTACGT", 2), rec(Side::Left, "", 2)]), 0);
        // python islower: no cased char -> False -> uppercase
        assert!(!py_islower(b"NNNN"));
        assert!(py_islower(b"nnAcn".to_ascii_lowercase().as_slice()));
    }

    /// The indexed `_absorb_target` equals python's literal scan (first minimum of
    /// `(one-sided, d, name)` in list order) on random data incl. duplicate names and dead
    /// candidates.
    #[test]
    fn absorb_target_matches_python_scan() {
        let c = Interner::new();
        let mut seed: u64 = 0x9e3779b97f4a7c15;
        let mut rnd = |m: u64| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed % m
        };
        for round in 0..30 {
            let n = 40 + rnd(80) as usize;
            let mut v = Vec::new();
            for u in 0..n {
                let contig = ["chr1", "chr2"][rnd(2) as usize];
                let l = 1000 + rnd(40) as i64;
                let (name, ty, open) = match rnd(3) {
                    0 => (format!("{contig}:{l}-{}", l + rnd(20) as i64), InsType::FullInfo, None),
                    1 => (format!("{contig}:{l}-oneside_{l}"), InsType::RightDisc, Some(Side::Right)),
                    _ => (format!("{contig}:oneside_{l}-{l}"), InsType::LeftDisc, Some(Side::Left)),
                };
                v.push(ins(&c, u as u32, &name, ty, open, vec![0]));
            }
            let names: Vec<String> = v.iter().map(|i| i.name(&c)).collect();
            let one_sided: Vec<bool> = v.iter().map(|i| i.open_side_eff().is_some()).collect();
            let alive: Vec<bool> = (0..n).map(|_| rnd(5) != 0).collect();
            let index = AbsorbIndex::new(&v);
            for (xi, x) in v.iter().enumerate() {
                let Some(o) = x.open_side_eff() else { continue };
                let side = o.other();
                let pos = x.junction(side).unwrap();
                let tol = 1 + (round % 7) as i64;
                // python literal
                let mut best: Option<((bool, i64, &str), usize)> = None;
                for (ti, t) in v.iter().enumerate() {
                    if ti == xi || !alive[ti] || t.contig != x.contig || t.open_side_eff() == Some(side) {
                        continue;
                    }
                    let d = (t.junction(side).unwrap() - pos).abs();
                    if d <= tol {
                        let key = (t.open_side_eff().is_some(), d, names[ti].as_str());
                        if best.is_none_or(|b| key < b.0) {
                            best = Some((key, ti));
                        }
                    }
                }
                let got = index.target(xi, x.contig, side, pos, tol, &alive, &one_sided, &names);
                assert_eq!(got, best.map(|b| b.1), "round {round} x {}", names[xi]);
            }
        }
        let _ = Tok::pos(0);
    }

    /// apply_evidence / absorb need P3 (row parse, evaluate_junction, render_reads), P4
    /// (aligned_part, filters) and P1 (Config, Genome). Equivalence across thread counts is
    /// covered by `THREAD_CHECK=1 tests/equiv.sh` at integration.
    #[test]
    #[ignore = "needs P1/P3/P4 implementations; covered by tests/equiv.sh (THREAD_CHECK=1) at integration"]
    fn apply_evidence_thread_invariance() {}
}
