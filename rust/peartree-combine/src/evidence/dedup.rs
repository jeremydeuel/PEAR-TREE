//! Fragments and the lenient within-sample dedup (independent_clusters). OWNER: P3.
//!
//! Mirrors `Fragment`, `collapse_fragments`, `DedupParams`, `_hp_compress_cut`, `_raw_cut`,
//! `_semi_close`, `_prefix_close`, `_seq_close_from_start`, `_has_long_run`, `_seq_close`,
//! `_junction_pos`, `_dup_seqs`, `_clip_shift`, `_outer_in_clip`, `_mate_unreliable`, `_is_dup`,
//! `_is_dup_oriented`, `_identical`, `independent_clusters` of
//! src/combine_insertions_evidence.py:178-519. SPEC.md §4.2.
//!
//! NOT a gate: the clusters decide the consensus vote weights (1/|cluster|), the consensus depth
//! unit (cluster id = `group`) and the reported `n_independent` / `supported` columns, all of which
//! reach the outputs (SPEC.md §0 "n_independent decision").

use crate::align::{self, Mode};
use crate::config::Config;
use crate::evidence::row::{needs_revcomp, EvidenceRow, Role};
use crate::model::{ContigId, InputFile, Interner, Side};
use crate::seq::revcomp;
use rustc_hash::FxHashMap;

/// All rows of one template (sample, frag) at one junction. Holds indices into the junction's
/// row slice (`rows[i]`), so a Fragment borrows nothing.
///
/// NB python keys fragments / samples by the sample NAME; here by FileId. They coincide because
/// input basenames are unique (two basenames `x.txt.gz` and `x` would share python's sample
/// name `x` -- not a real input).
#[derive(Clone, Debug)]
pub struct Fragment {
    /// FileId of the sample
    pub file: u32,
    pub frag: Box<str>,
    /// indices of this fragment's rows, in row order
    pub rows: Vec<usize>,
    /// `primary`: among non-MATE rows, min by (role priority, r12) -- first wins on ties
    pub primary: usize,
    /// `mate`: among rows with r12 != primary.r12, min by (0 if MATE else 1, r12) -- first wins
    pub mate: Option<usize>,
}

impl Fragment {
    /// python `Fragment(sample, frag, rows)`; None when it has no primary (MATE rows only).
    fn build(file: u32, frag: Box<str>, idx: Vec<usize>, rows: &[EvidenceRow]) -> Option<Fragment> {
        let mut primary: Option<usize> = None;
        for &i in &idx {
            let r = &rows[i];
            if r.role == Role::Mate {
                continue;
            }
            let better = match primary {
                None => true,
                Some(p) => (r.role.priority(), r.r12) < (rows[p].role.priority(), rows[p].r12),
            };
            if better {
                primary = Some(i);
            }
        }
        let primary = primary?;
        let pr12 = rows[primary].r12;
        let key = |r: &EvidenceRow| (if r.role == Role::Mate { 0u8 } else { 1u8 }, r.r12);
        let mut mate: Option<usize> = None;
        for &i in &idx {
            let r = &rows[i];
            if r.r12 == pr12 {
                continue;
            }
            if mate.map_or(true, |m| key(r) < key(&rows[m])) {
                mate = Some(i);
            }
        }
        Some(Fragment { file, frag, rows: idx, primary, mate })
    }

    /// `Fragment.swapped()`: same rows with primary and mate exchanged (None without a mate).
    pub fn swapped(&self) -> Option<Fragment> {
        let m = self.mate?;
        Some(Fragment { file: self.file, frag: self.frag.clone(), rows: self.rows.clone(), primary: m, mate: Some(self.primary) })
    }
}

/// `collapse_fragments(rows)`: group by (sample, frag) in first-appearance order, build
/// Fragments, drop those without a primary, then STABLE sort by (sample NAME string, primary
/// outer, frag string). `files` gives the sample names.
pub fn collapse_fragments(rows: &[EvidenceRow], files: &[InputFile]) -> Vec<Fragment> {
    let mut index: FxHashMap<(u32, &str), usize> = FxHashMap::default();
    let mut groups: Vec<(u32, &str, Vec<usize>)> = Vec::new();
    for (i, r) in rows.iter().enumerate() {
        let k = (r.file, &*r.frag);
        match index.get(&k) {
            Some(&g) => groups[g].2.push(i),
            None => {
                index.insert(k, groups.len());
                groups.push((r.file, &*r.frag, vec![i]));
            }
        }
    }
    let mut frags: Vec<Fragment> =
        groups.into_iter().filter_map(|(f, fr, idx)| Fragment::build(f, fr.into(), idx, rows)).collect();
    frags.sort_by(|a, b| {
        files[a.file as usize]
            .sample
            .as_bytes()
            .cmp(files[b.file as usize].sample.as_bytes())
            .then(rows[a.primary].outer.cmp(&rows[b.primary].outer))
            .then(a.frag.as_bytes().cmp(b.frag.as_bytes()))
    });
    frags
}

/// python `DedupParams` (`from_cfg`).
#[derive(Clone, Copy, Debug)]
pub struct DedupParams {
    pub tol: i64,
    pub max_edit: i64,
    pub max_edit_frac: f64,
    pub polya_min: usize,
    pub mate_min_mapq: i64,
}

impl DedupParams {
    pub fn from_cfg(cfg: &Config) -> DedupParams {
        DedupParams {
            tol: cfg.dup_coord_tolerance,
            max_edit: cfg.dup_max_edit,
            max_edit_frac: cfg.dup_max_edit_frac,
            polya_min: cfg.polya_min_len,
            mate_min_mapq: cfg.dup_mate_min_mapq,
        }
    }
    /// `budget(n)` -> `crate::pyfmt::dedup_budget`.
    pub fn budget(&self, n: usize) -> i64 {
        crate::pyfmt::dedup_budget(self.max_edit, self.max_edit_frac, n)
    }
}

/// Dedup statistics (`stats` dict: n_dup_coord / n_dup_seq).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DedupStats {
    pub n_dup_coord: u32,
    pub n_dup_seq: u32,
}

/// Result of `independent_clusters`.
#[derive(Clone, Debug, Default)]
pub struct Clusters {
    /// clusters as indices into the `frags` argument; clusters ordered by their smallest member
    /// index, members ascending (python `groups` dict order -- independent of union-find roots)
    pub clusters: Vec<Vec<usize>>,
    pub n_dup: u32,
    pub n_cross: u32,
    pub stats: DedupStats,
}

// ------------------------------------------------------------------ dedup internals
//
// Everything the O(n^2) pair loop needs per row is computed once per call (`RowInfo`), and the
// lenient sequence comparison reuses the buffers of one `Scratch`, so the inner loop does not
// allocate. `V` is a (primary, mate) view of a fragment -- `Fragment.swapped()` without a clone.

/// Per-row values the dedup reads repeatedly.
#[derive(Default)]
struct RowInfo {
    mapped: bool,
    mate_mapped: bool,
    /// `seq.upper()`
    seq_up: Vec<u8>,
    /// `allele_forward_seq(r).upper()`
    af_up: Vec<u8>,
    /// `outward_clip()[0].upper()` (CLIP rows only)
    clip_up: Vec<u8>,
}

fn row_info(r: &EvidenceRow, contigs: &Interner) -> RowInfo {
    let seq_up = r.seq.to_ascii_uppercase();
    let af_up = if needs_revcomp(r) { revcomp(&seq_up) } else { seq_up.clone() };
    let clip_up = if r.role == Role::Clip { r.outward_clip_seq().to_ascii_uppercase() } else { Vec::new() };
    RowInfo { mapped: r.mapped(contigs), mate_mapped: r.mate_mapped(contigs), seq_up, af_up, clip_up }
}

struct Dd<'a> {
    rows: &'a [EvidenceRow],
    info: Vec<RowInfo>,
    p: &'a DedupParams,
}

impl<'a> Dd<'a> {
    /// RowInfo for every row that is some fragment's primary or mate (others stay empty).
    fn new<'f>(rows: &'a [EvidenceRow], frags: impl Iterator<Item = &'f Fragment>, p: &'a DedupParams, contigs: &Interner) -> Dd<'a> {
        let mut info: Vec<RowInfo> = Vec::with_capacity(rows.len());
        info.resize_with(rows.len(), RowInfo::default);
        let mut done = vec![false; rows.len()];
        for f in frags {
            for i in std::iter::once(f.primary).chain(f.mate) {
                if !done[i] {
                    done[i] = true;
                    info[i] = row_info(&rows[i], contigs);
                }
            }
        }
        Dd { rows, info, p }
    }
}

#[derive(Clone, Copy)]
struct V {
    pr: usize,
    mate: Option<usize>,
}

impl V {
    fn of(f: &Fragment) -> V {
        V { pr: f.primary, mate: f.mate }
    }
    fn swapped(self) -> Option<V> {
        self.mate.map(|m| V { pr: m, mate: Some(self.pr) })
    }
}

type Place<'r> = (ContigId, i64, &'r str);

impl<'a> Dd<'a> {
    #[inline]
    fn strand(&self, v: V) -> &'a str {
        &self.rows[v.pr].strand
    }
    #[inline]
    fn outer(&self, v: V) -> i64 {
        self.rows[v.pr].outer
    }
    /// `Fragment.mate_coord()`
    fn mate_coord(&self, v: V) -> Option<Place<'a>> {
        let p = &self.rows[v.pr];
        if !self.info[v.pr].mate_mapped {
            return None;
        }
        Some((p.mref, p.mpos, &p.mstrand))
    }
    /// `Fragment.mate_outer()`
    fn mate_outer(&self, v: V) -> Option<Place<'a>> {
        let m = v.mate?;
        let r = &self.rows[m];
        if self.info[m].mapped && r.outer >= 0 {
            Some((r.ref_, r.outer, &r.strand))
        } else {
            None
        }
    }
    /// `Fragment.clip_seq()`
    fn clip_seq(&self, v: V) -> &[u8] {
        if self.rows[v.pr].role == Role::Clip {
            &self.info[v.pr].clip_up
        } else {
            &self.info[v.pr].seq_up
        }
    }
    /// `Fragment.mate_seq()`
    fn mate_seq(&self, v: V) -> &[u8] {
        match v.mate {
            Some(m) => &self.info[m].af_up,
            None => b"",
        }
    }
    /// `_junction_pos(f)`
    fn junction_pos(&self, v: V) -> Option<i64> {
        let p = &self.rows[v.pr];
        if p.role != Role::Clip || !self.info[v.pr].mapped {
            return None;
        }
        Some(match p.side {
            Side::Left => p.pos,
            Side::Right => p.pos + crate::evidence::row::cigar_ref_len(&p.cigar),
        })
    }
    fn both_clip(&self, a: V, b: V) -> bool {
        self.rows[a.pr].role == Role::Clip && self.rows[b.pr].role == Role::Clip
    }
    /// `_dup_seqs(a, b)`
    fn dup_seqs(&self, a: V, b: V) -> (&[u8], &[u8]) {
        if self.both_clip(a, b) {
            (self.clip_seq(a), self.clip_seq(b))
        } else {
            (&self.info[a.pr].af_up, &self.info[b.pr].af_up)
        }
    }
    /// `_clip_shift(a, b, p)`
    fn clip_shift(&self, a: V, b: V) -> i64 {
        if self.both_clip(a, b) {
            0
        } else {
            self.p.tol
        }
    }
    /// `_outer_in_clip(f)`
    fn outer_in_clip(&self, v: V) -> bool {
        let p = &self.rows[v.pr];
        p.role == Role::Clip && ((p.side == Side::Left && &*p.strand == "+") || (p.side == Side::Right && &*p.strand == "-"))
    }
    /// `_mate_unreliable(f, p)`
    fn mate_unreliable(&self, v: V) -> bool {
        match v.mate {
            Some(m) => self.info[m].mapped && self.rows[m].mapq < self.p.mate_min_mapq,
            None => false,
        }
    }

    /// `_is_dup(a, b, p)`
    fn is_dup(&self, a: V, b: V, s: &mut Scratch) -> &'static str {
        let k = self.is_dup_oriented(a, b, s);
        if !k.is_empty() || self.strand(a) == self.strand(b) {
            return k;
        }
        for (x, y) in [(Some(a), b.swapped()), (a.swapped(), Some(b))] {
            let (Some(x), Some(y)) = (x, y) else { continue };
            if !self.info[x.pr].mapped || !self.info[y.pr].mapped {
                continue;
            }
            let k = self.is_dup_oriented(x, y, s);
            if !k.is_empty() {
                return k;
            }
        }
        ""
    }

    /// `_is_dup_oriented(a, b, p)`
    fn is_dup_oriented(&self, a: V, b: V, s: &mut Scratch) -> &'static str {
        let p = self.p;
        let tol = p.tol;
        if self.strand(a) != self.strand(b) {
            return "";
        }
        if let (Some(ja), Some(jb)) = (self.junction_pos(a), self.junction_pos(b)) {
            if (ja - jb).abs() > tol {
                return "";
            }
        }
        let otol = if self.outer_in_clip(a) && self.outer_in_clip(b) { 2 * tol } else { tol };
        if (self.outer(a) - self.outer(b)).abs() > otol {
            return "";
        }
        let (am, bm) = if self.mate_unreliable(a) || self.mate_unreliable(b) {
            (None, None)
        } else {
            match (self.mate_outer(a), self.mate_outer(b)) {
                (Some(ao), Some(bo)) => (Some(ao), Some(bo)),
                _ => (self.mate_coord(a), self.mate_coord(b)),
            }
        };
        if am.is_none() != bm.is_none() {
            return "";
        }
        if let (Some(am), Some(bm)) = (am, bm) {
            if !(am.0 == bm.0 && am.2 == bm.2 && (am.1 - bm.1).abs() <= tol) {
                return "";
            }
            let (x, y) = self.dup_seqs(a, b);
            return if seq_close_s(x, y, p, self.clip_shift(a, b), s) { "coord" } else { "" };
        }
        let (x, y) = self.dup_seqs(a, b);
        if !seq_close_s(x, y, p, self.clip_shift(a, b), s) {
            return "";
        }
        let (ams, bms) = (self.mate_seq(a), self.mate_seq(b));
        if !ams.is_empty() && !bms.is_empty() && !seq_close_s(ams, bms, p, tol, s) {
            return "";
        }
        "seq"
    }

    /// `_identical(a, b)`
    fn identical(&self, a: V, b: V) -> bool {
        if self.strand(a) != self.strand(b) || self.outer(a) != self.outer(b) {
            return false;
        }
        let am = self.mate_outer(a).or_else(|| self.mate_coord(a));
        let bm = self.mate_outer(b).or_else(|| self.mate_coord(b));
        if am != bm {
            return false;
        }
        if self.info[a.pr].seq_up != self.info[b.pr].seq_up {
            return false;
        }
        self.mate_seq(a) == self.mate_seq(b)
    }
}

/// `independent_clusters(frags, params, stats)` (evidence.py:472): within-sample pairs
/// (samples in first-appearance order of `frags`; each sample's indices STABLE-sorted by primary
/// outer; all pairs i<j in that order, skipping already-connected pairs, `_is_dup` -> union,
/// n_dup += 1, stats by kind), then cross-sample exact identity (`by_key[(strand, outer)]` groups
/// in first-appearance order, pairs in index order, different sample, not connected,
/// `_identical` -> union, n_cross += 1).
pub fn independent_clusters(frags: &[Fragment], rows: &[EvidenceRow], p: &DedupParams, contigs: &Interner) -> Clusters {
    let n = frags.len();
    let dd = Dd::new(rows, frags.iter(), p, contigs);
    let mut s = Scratch::default();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    let v: Vec<V> = frags.iter().map(V::of).collect();
    let mut out = Clusters::default();

    // within-sample (python by_sample, dict order = first appearance)
    let mut sample_ix: Vec<(u32, Vec<usize>)> = Vec::new();
    for (i, f) in frags.iter().enumerate() {
        match sample_ix.iter_mut().find(|(s, _)| *s == f.file) {
            Some((_, idx)) => idx.push(i),
            None => sample_ix.push((f.file, vec![i])),
        }
    }
    for (_, idx) in sample_ix.iter_mut() {
        idx.sort_by_key(|&i| rows[frags[i].primary].outer);
        for a_pos in 0..idx.len() {
            let i = idx[a_pos];
            for &j in &idx[a_pos + 1..] {
                if find(&mut parent, i) == find(&mut parent, j) {
                    continue;
                }
                let kind = dd.is_dup(v[i], v[j], &mut s);
                if !kind.is_empty() {
                    let (rj, ri) = (find(&mut parent, j), find(&mut parent, i));
                    parent[rj] = ri;
                    out.n_dup += 1;
                    if kind == "coord" {
                        out.stats.n_dup_coord += 1;
                    } else {
                        out.stats.n_dup_seq += 1;
                    }
                }
            }
        }
    }
    // cross-sample: exact identity only (python by_key, dict order = first appearance)
    let mut key_ix: FxHashMap<(&str, i64), usize> = FxHashMap::default();
    let mut by_key: Vec<Vec<usize>> = Vec::new();
    for (i, f) in frags.iter().enumerate() {
        let r = &rows[f.primary];
        let k = (&*r.strand, r.outer);
        match key_ix.get(&k) {
            Some(&g) => by_key[g].push(i),
            None => {
                key_ix.insert(k, by_key.len());
                by_key.push(vec![i]);
            }
        }
    }
    for idx in &by_key {
        for a_pos in 0..idx.len() {
            let i = idx[a_pos];
            for &j in &idx[a_pos + 1..] {
                if frags[i].file != frags[j].file
                    && find(&mut parent, i) != find(&mut parent, j)
                    && dd.identical(v[i], v[j])
                {
                    let (rj, ri) = (find(&mut parent, j), find(&mut parent, i));
                    parent[rj] = ri;
                    out.n_cross += 1;
                }
            }
        }
    }
    // groups in order of first member (python dict keyed by root, filled for i = 0..n)
    let mut slot: Vec<usize> = vec![usize::MAX; n];
    for i in 0..n {
        let r = find(&mut parent, i);
        if slot[r] == usize::MAX {
            slot[r] = out.clusters.len();
            out.clusters.push(Vec::new());
        }
        out.clusters[slot[r]].push(i);
    }
    out
}

/// `_is_dup(a, b, p)` -> "" / "coord" / "seq" (as an enum-free &'static str for parity).
pub fn is_dup(a: &Fragment, b: &Fragment, rows: &[EvidenceRow], p: &DedupParams, contigs: &Interner) -> &'static str {
    let dd = Dd::new(rows, [a, b].into_iter(), p, contigs);
    dd.is_dup(V::of(a), V::of(b), &mut Scratch::default())
}

/// `_seq_close(a, b, p, shift)` (exposed for unit tests against the python).
pub fn seq_close(a: &[u8], b: &[u8], p: &DedupParams, shift: i64) -> bool {
    seq_close_s(a, b, p, shift, &mut Scratch::default())
}

// ------------------------------------------------------------------ lenient sequence identity

/// Reusable buffers of the sequence comparison (homopolymer-compressed cuts, reversed copies).
#[derive(Default)]
struct Scratch {
    ha: Vec<u8>,
    hb: Vec<u8>,
    ra: Vec<u8>,
    rb: Vec<u8>,
}

/// `_hp_compress_cut(seq, polya_min)` into `out`.
fn hp_compress_cut(seq: &[u8], polya_min: usize, out: &mut Vec<u8>) {
    out.clear();
    let n = seq.len();
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j < n && seq[j] == seq[i] {
            j += 1;
        }
        out.push(seq[i]);
        if matches!(seq[i], b'A' | b'T') && j - i >= polya_min {
            break;
        }
        i = j;
    }
}

/// `_raw_cut(seq, polya_min)` (a prefix of `seq`).
fn raw_cut(seq: &[u8], polya_min: usize) -> &[u8] {
    let n = seq.len();
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j < n && seq[j] == seq[i] {
            j += 1;
        }
        if matches!(seq[i], b'A' | b'T') && j - i >= polya_min {
            return &seq[..(i + polya_min).min(n)];
        }
        i = j;
    }
    seq
}

/// `_semi_close(a, b, p.budget)`.
fn semi_close(a: &[u8], b: &[u8], p: &DedupParams) -> bool {
    let (a, b) = if a.len() > b.len() { (b, a) } else { (a, b) };
    let m = a.len();
    if m == 0 || b.starts_with(a) {
        return true;
    }
    let k = p.budget(m).min(i32::MAX as i64) as i32;
    align::distance(a, b, Mode::Shw, k, &[]) != -1
}

/// `_prefix_close(a, b, p)`.
fn prefix_close(a: &[u8], b: &[u8], p: &DedupParams, s: &mut Scratch) -> bool {
    hp_compress_cut(a, p.polya_min, &mut s.ha);
    hp_compress_cut(b, p.polya_min, &mut s.hb);
    if semi_close(&s.ha, &s.hb, p) {
        return true;
    }
    semi_close(raw_cut(a, p.polya_min), raw_cut(b, p.polya_min), p)
}

/// `_seq_close_from_start(a, b, p, shift)` (python slices past the end give "").
fn seq_close_from_start(a: &[u8], b: &[u8], p: &DedupParams, shift: i64, s: &mut Scratch) -> bool {
    if shift == 0 {
        return prefix_close(a, b, p, s);
    }
    let mut k: i64 = 0;
    while k <= shift {
        let ku = k as usize;
        if prefix_close(&a[ku.min(a.len())..], b, p, s) || (k != 0 && prefix_close(a, &b[ku.min(b.len())..], p, s)) {
            return true;
        }
        k += 1;
    }
    false
}

/// `_has_long_run(x, n)`: `"A"*n in x or "T"*n in x`.
fn has_long_run(x: &[u8], n: usize) -> bool {
    if n == 0 {
        return true;
    }
    let mut run = 0usize;
    let mut prev = 0u8;
    for &c in x {
        run = if c == prev { run + 1 } else { 1 };
        prev = c;
        if run >= n && matches!(c, b'A' | b'T') {
            return true;
        }
    }
    false
}

/// `_seq_close(a, b, p, shift)`.
fn seq_close_s(a: &[u8], b: &[u8], p: &DedupParams, shift: i64, s: &mut Scratch) -> bool {
    if seq_close_from_start(a, b, p, shift, s) {
        return true;
    }
    if has_long_run(a, p.polya_min) || has_long_run(b, p.polya_min) {
        let mut ra = std::mem::take(&mut s.ra);
        let mut rb = std::mem::take(&mut s.rb);
        ra.clear();
        ra.extend(a.iter().rev());
        rb.clear();
        rb.extend(b.iter().rev());
        let r = seq_close_from_start(&ra, &rb, p, shift, s);
        s.ra = ra;
        s.rb = rb;
        return r;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::row::p3_fixture::{self, Fixture};

    fn check_fixture(fx: &Fixture) -> (usize, u32, u32) {
        let (mut n, mut dups, mut cross) = (0, 0, 0);
        for (j, (locus, side, rows)) in fx.doc["junctions"].as_array().unwrap().iter().zip(&fx.junctions).map(|(j, (l, s, r))| (j, (l, s, r))) {
            let ctx = format!("{locus} {}", side.as_str());
            let frags = collapse_fragments(rows, &fx.files);
            let exp = j["frags"].as_array().unwrap();
            assert_eq!(frags.len(), exp.len(), "{ctx}");
            for (f, e) in frags.iter().zip(exp) {
                assert_eq!(fx.files[f.file as usize].sample, e[0].as_str().unwrap(), "{ctx}");
                assert_eq!(&*f.frag, e[1].as_str().unwrap(), "{ctx}");
                assert_eq!(f.primary as i64, e[2].as_i64().unwrap(), "{ctx}");
                assert_eq!(f.mate.map(|m| m as i64), e[3].as_i64(), "{ctx}");
            }
            for c in ["A", "B", "C"] {
                let p = DedupParams::from_cfg(&p3_fixture::config(&fx.doc["cfgs"][c]));
                let e = &j[c];
                let got = independent_clusters(&frags, rows, &p, &fx.contigs);
                let ec: Vec<Vec<usize>> = serde_json::from_value(e["clusters"].clone()).unwrap();
                assert_eq!(got.clusters, ec, "{ctx} cfg {c}");
                assert_eq!(got.n_dup as i64, e["n_dup"].as_i64().unwrap(), "{ctx} cfg {c}");
                assert_eq!(got.n_cross as i64, e["n_cross"].as_i64().unwrap(), "{ctx} cfg {c}");
                assert_eq!(got.stats.n_dup_coord as i64, e["n_dup_coord"].as_i64().unwrap(), "{ctx} cfg {c}");
                assert_eq!(got.stats.n_dup_seq as i64, e["n_dup_seq"].as_i64().unwrap(), "{ctx} cfg {c}");
                n += 1;
                dups += got.n_dup;
                cross += got.n_cross;
            }
        }
        (n, dups, cross)
    }

    #[test]
    fn fragments_and_clusters_match_python() {
        let (n, dups, cross) = check_fixture(&p3_fixture::checked_in());
        assert!(n >= 120 && dups > 100 && cross > 0, "{n} {dups} {cross}");
    }

    /// All 2872 junctions of the E2E fixture x 3 configs (`P3_FULL_DIR`, see p3_fixture).
    #[test]
    #[ignore]
    fn fragments_and_clusters_match_python_full() {
        let Some(fx) = p3_fixture::full() else { return eprintln!("P3_FULL_DIR not set: skipped") };
        let (n, dups, cross) = check_fixture(&fx);
        eprintln!("full fixture: {n} junction x config evaluations, {dups} dups, {cross} cross");
    }

    #[test]
    fn seq_close_matches_python() {
        let fx = p3_fixture::checked_in();
        let ps = [
            DedupParams { tol: 5, max_edit: 3, max_edit_frac: 0.02, polya_min: 8, mate_min_mapq: 20 },
            DedupParams { tol: 5, max_edit: 2, max_edit_frac: 0.05, polya_min: 6, mate_min_mapq: 20 },
        ];
        let cases = fx.doc["seq_close"].as_array().unwrap();
        let (mut t, mut f) = (0, 0);
        let mut s = Scratch::default();
        for c in cases {
            let (a, b) = (c[0].as_str().unwrap().as_bytes(), c[1].as_str().unwrap().as_bytes());
            let shift = c[2].as_i64().unwrap();
            let p = &ps[c[3].as_u64().unwrap() as usize];
            let want = c[4].as_bool().unwrap();
            assert_eq!(seq_close(a, b, p, shift), want, "{:?} {:?} shift {shift}", c[0], c[1]);
            // reused scratch gives the same answer
            assert_eq!(seq_close_s(a, b, p, shift, &mut s), want);
            if want { t += 1 } else { f += 1 }
        }
        assert!(t > 1000 && f > 150, "{t} {f}");
    }

    #[test]
    fn cut_helpers() {
        let mut o = Vec::new();
        hp_compress_cut(b"GGAAAAAAAAACCT", 8, &mut o);
        assert_eq!(o, b"GA");
        hp_compress_cut(b"GGAAAACCT", 8, &mut o);
        assert_eq!(o, b"GACT");
        assert_eq!(raw_cut(b"GGAAAAAAAAAACC", 8), b"GGAAAAAAAA");
        assert_eq!(raw_cut(b"GGAAAC", 8), b"GGAAAC");
        assert!(has_long_run(b"CCTTTT", 4) && !has_long_run(b"CCTTTGGGG", 4) && has_long_run(b"", 0));
    }
}
