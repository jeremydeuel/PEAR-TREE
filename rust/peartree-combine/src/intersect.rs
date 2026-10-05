//! Pool per-sample discovery records by locus. OWNER: P1.
//!
//! Mirrors src/combine_insertions_intersect_insertions.py (`intersect_insertions` and helpers).
//! SPEC.md §3.3 -- read it for the exact ordering rules: the OUTPUT ORDER of this function is
//! the order of every later output file, and it is python dict insertion order.

use crate::align::{self, Mode};
use crate::config::Config;
use crate::model::{FileId, InsType, Insertion, InputFile, Interner, LocusKey, Member, Side, SideSet};
use crate::seq::{clips_agree, polya_trimmed, sequence_matching_score, ScoreSeq};
use rustc_hash::FxHashMap;
use std::hash::Hash;

/// An insertion-ordered multimap (python `dict` of lists keyed in first-seen order).
struct OrderedBuckets<K: Hash + Eq + Copy> {
    keys: Vec<K>,
    vals: Vec<Vec<Insertion>>,
    index: FxHashMap<K, usize>,
}

impl<K: Hash + Eq + Copy> OrderedBuckets<K> {
    fn new() -> Self {
        OrderedBuckets { keys: Vec::new(), vals: Vec::new(), index: FxHashMap::default() }
    }
    fn push(&mut self, k: K, v: Insertion) {
        match self.index.get(&k) {
            Some(&i) => self.vals[i].push(v),
            None => {
                self.index.insert(k, self.vals.len());
                self.keys.push(k);
                self.vals.push(vec![v]);
            }
        }
    }
    fn len(&self) -> usize {
        self.keys.len()
    }
}

type TKey = (u32, Option<i64>, Option<i64>);

/// `_polya_to_one_sided` (ISX:27): the poly-A end keeps its coordinate (from the name token
/// `polyA_<P>`) but no sequence; the name is unchanged.
fn polya_to_one_sided(mut i: Insertion) -> Insertion {
    if i.ty == InsType::RightPolyA {
        i.right_clipped = None;
        i.right_aligned = None;
        i.right_pos = Some(i.name_end.pos);
        i.ty = InsType::RightDisc;
        i.open_side = Some(Side::Right);
    } else {
        i.left_clipped = None;
        i.left_aligned = None;
        i.left_pos = Some(i.name_start.pos);
        i.ty = InsType::LeftDisc;
        i.open_side = Some(Side::Left);
    }
    i
}

/// the clip of the real side of a one-sided record (`_real_clip`)
fn real_clip(i: &Insertion) -> Option<&crate::seq::QualSeq> {
    if i.ty == InsType::LeftDisc {
        i.right_clipped.as_ref()
    } else {
        i.left_clipped.as_ref()
    }
}

/// `_real_sides(i)[0]`
fn first_real_side(i: &Insertion) -> Side {
    match i.ty {
        InsType::RightDisc => Side::Left,
        InsType::LeftDisc => Side::Right,
        _ => Side::Left,
    }
}

fn clip_bytes(i: &Insertion, side: Side) -> Option<&[u8]> {
    i.clipped(side).map(|q| &q.seq[..])
}

fn pos_of(i: &Insertion, side: Side) -> i64 {
    i.junction(side).expect("one-sided record without a real-side coordinate")
}

/// `intersect_insertions(insertions)` with the three settings taken from `cfg`
/// (`keep_polya_one_sided`, `merge_tolerance_bp`, `polya_aware_clip_agreement`).
///
/// Input: every accepted record, files in command-line order, records in file order.
/// Output: surviving insertions in python `full_insertions.values()` order (SPEC.md §3.3), with
/// `uid` = output index (0..n).
pub fn intersect_insertions(records: Vec<Insertion>, cfg: &Config, contigs: &Interner, files: &[InputFile]) -> Vec<Insertion> {
    let keep_polya_one_sided = cfg.keep_polya_one_sided;
    let tol = cfg.merge_tolerance_bp;
    let polya_aware = cfg.polya_aware_clip_agreement;

    // ---- step 1: bucket
    let mut full: OrderedBuckets<TKey> = OrderedBuckets::new();
    let mut polya: OrderedBuckets<TKey> = OrderedBuckets::new();
    let mut disc: OrderedBuckets<TKey> = OrderedBuckets::new();
    for i in records {
        let k: TKey = (i.contig, i.left_pos, i.right_pos);
        match i.ty {
            InsType::FullInfo => full.push(k, i),
            InsType::RightPolyA | InsType::LeftPolyA => polya.push(k, i),
            InsType::RightDisc | InsType::LeftDisc => disc.push(k, i),
        }
    }
    eprintln!(
        "imported {} unique full-information insertions, {} polyA insertions and {} discordant-anchored insertions",
        full.len(),
        polya.len(),
        disc.len()
    );

    // ---- step 2: per full key
    let full_keys: Vec<TKey> = full.keys;
    let mut full_vals: Vec<Option<Insertion>> = Vec::with_capacity(full_keys.len());
    for hits in full.vals {
        full_vals.push(combine_hits(hits, polya_aware));
    }

    // ---- step 3: fuzzy merge of live full insertions
    if tol > 0 {
        let live: Vec<usize> = (0..full_vals.len()).filter(|&i| full_vals[i].is_some()).collect();
        let keys: Vec<(u32, i64, i64)> = live
            .iter()
            .map(|&i| (full_keys[i].0, full_keys[i].1.expect("full key"), full_keys[i].2.expect("full key")))
            .collect();
        let weight: Vec<usize> = live.iter().map(|&i| full_vals[i].as_ref().unwrap().member_loci.len()).collect();
        for cl in fuzzy_clusters(&keys, &weight, tol, contigs) {
            let rep_i = live[cl[0]];
            for &k in &cl[1..] {
                let o_i = live[k];
                let agree = {
                    let rep = full_vals[rep_i].as_ref().unwrap();
                    let o = full_vals[o_i].as_ref().unwrap();
                    [Side::Left, Side::Right]
                        .iter()
                        .all(|&s| shift_tolerant_agree(clip_bytes(rep, s), clip_bytes(o, s), tol, polya_aware))
                };
                if agree {
                    let o = full_vals[o_i].take().unwrap();
                    absorb(full_vals[rep_i].as_mut().unwrap(), &o, &[Side::Left, Side::Right]);
                }
            }
        }
    }

    // ---- step 4: poly-A records
    if keep_polya_one_sided {
        for (k, hits) in polya.keys.iter().zip(polya.vals) {
            for h in hits {
                disc.push(*k, polya_to_one_sided(h));
            }
        }
    }
    // (otherwise the poly-A records are dropped: python `continue`s past the rest of that loop)

    // ---- step 5: discordant / one-sided representatives, keyed by NAME string
    let mut named_keys: Vec<LocusKey> = Vec::new();
    let mut named_vals: Vec<Option<Insertion>> = Vec::new();
    let mut named_index: FxHashMap<LocusKey, usize> = FxHashMap::default();
    for hits in disc.vals {
        let mut rep_i = 0usize;
        for (hi, h) in hits.iter().enumerate().skip(1) {
            let real_side = real_clip(h);
            let rep_side = real_clip(&hits[rep_i]);
            if let Some(rs) = real_side {
                if rep_side.map_or(true, |p| rs.len() > p.len()) {
                    rep_i = hi;
                }
            }
        }
        let mut hits: Vec<Option<Insertion>> = hits.into_iter().map(Some).collect();
        let mut rep = hits[rep_i].take().unwrap();
        if rep.open_side.is_some() {
            for h in hits.iter().flatten() {
                let new_files: Vec<FileId> = h.files.iter().copied().filter(|f| !rep.files.contains(f)).collect();
                rep.files.extend(new_files);
                if rep.ty == InsType::RightDisc {
                    if let (Some(ha), Some(ra)) = (h.left_aligned.as_ref(), rep.left_aligned.as_ref()) {
                        if ha.len() > ra.len() {
                            rep.left_aligned = h.left_aligned.clone();
                        }
                    }
                }
                if rep.ty == InsType::LeftDisc {
                    if let (Some(ha), Some(ra)) = (h.right_aligned.as_ref(), rep.right_aligned.as_ref()) {
                        if ha.len() > ra.len() {
                            rep.right_aligned = h.right_aligned.clone();
                        }
                    }
                }
            }
        }
        let name = rep.locus();
        match named_index.get(&name) {
            Some(&i) => named_vals[i] = Some(rep),
            None => {
                named_index.insert(name, named_vals.len());
                named_keys.push(name);
                named_vals.push(Some(rep));
            }
        }
    }

    // ---- step 6: fuzzy merge of one-sided loci
    if tol > 0 {
        fuzzy_one_sided(&mut named_vals, tol, polya_aware, contigs);
    }

    // ---- step 7
    let mut out: Vec<Insertion> = full_vals.into_iter().chain(named_vals).flatten().collect();
    for (n, i) in out.iter_mut().enumerate() {
        i.uid = n as u32;
    }
    out
}

/// One full-key bucket -> the surviving record (None = clips/flanks disagree).
fn combine_hits(hits: Vec<Insertion>, polya_aware: bool) -> Option<Insertion> {
    if hits.len() == 1 {
        return hits.into_iter().next();
    }
    let q = |f: &dyn Fn(&Insertion) -> &Option<crate::seq::QualSeq>| -> Vec<ScoreSeq> {
        hits.iter().map(|s| ScoreSeq::Qual(f(s).as_ref().expect("full insertion lacks a side"))).collect()
    };
    if polya_aware {
        let rc: Vec<&[u8]> = hits.iter().map(|s| &s.right_clipped.as_ref().unwrap().seq[..]).collect();
        let lc: Vec<&[u8]> = hits.iter().map(|s| &s.left_clipped.as_ref().unwrap().seq[..]).collect();
        if !(clips_agree(&rc, 0.6, 8, 6) && clips_agree(&lc, 0.6, 8, 6)) {
            return None;
        }
    } else if sequence_matching_score(&q(&|s| &s.right_clipped)) < 0.6 {
        return None;
    }
    if !polya_aware && sequence_matching_score(&q(&|s| &s.left_clipped)) < 0.6 {
        return None;
    }
    if sequence_matching_score(&q(&|s| &s.left_aligned)) < 0.6 {
        return None;
    }
    if sequence_matching_score(&q(&|s| &s.right_aligned)) < 0.6 {
        return None;
    }
    let mut it = hits.into_iter();
    let mut combined = it.next().unwrap();
    for h in it {
        crate::insertion::merge_into(&mut combined, h);
    }
    Some(combined)
}

/// `_fuzzy_one_sided` (ISX:133) over the named (one-sided) values.
fn fuzzy_one_sided(vals: &mut [Option<Insertion>], tol: i64, polya_aware: bool, contigs: &Interner) {
    // groups keyed by (contig, real side); items: (index, name string)
    let mut order: Vec<(u32, Side)> = Vec::new();
    let mut groups: FxHashMap<(u32, Side), Vec<(usize, String)>> = FxHashMap::default();
    for (i, v) in vals.iter().enumerate() {
        if let Some(v) = v {
            if v.open_side.is_some() {
                let g = (v.contig, first_real_side(v));
                groups.entry(g).or_insert_with(|| {
                    order.push(g);
                    Vec::new()
                });
                groups.get_mut(&g).unwrap().push((i, v.name(contigs)));
            }
        }
    }
    for g in order {
        let s = g.1;
        let mut items = groups.remove(&g).unwrap();
        items.sort_by(|a, b| {
            let (va, vb) = (vals[a.0].as_ref().unwrap(), vals[b.0].as_ref().unwrap());
            (std::cmp::Reverse(va.files.len()), pos_of(va, s), &a.1).cmp(&(std::cmp::Reverse(vb.files.len()), pos_of(vb, s), &b.1))
        });
        let mut reps: Vec<usize> = Vec::new();
        for (vi, _) in items {
            let found = {
                let v = vals[vi].as_ref().unwrap();
                reps.iter().copied().find(|&ri| {
                    let r = vals[ri].as_ref().unwrap();
                    (pos_of(r, s) - pos_of(v, s)).abs() <= tol && shift_tolerant_agree(clip_bytes(r, s), clip_bytes(v, s), tol, polya_aware)
                })
            };
            match found {
                None => reps.push(vi),
                Some(ri) => {
                    let v = vals[vi].take().unwrap();
                    absorb(vals[ri].as_mut().unwrap(), &v, &[s]);
                }
            }
        }
    }
}

/// `_shift_tolerant_agree(a, b, tol, polya_aware)` (intersect:57): None on either side -> true;
/// polya_trimmed (polya_min 8) both when polya_aware else uppercase; shorter = x; len(x) < 6 ->
/// true; `edlib HW distance(x[:20], y[:20+tol+4], k=max(1, len(q)//4)) != -1`.
pub fn shift_tolerant_agree(a: Option<&[u8]>, b: Option<&[u8]>, tol: i64, polya_aware: bool) -> bool {
    let (Some(a), Some(b)) = (a, b) else { return true };
    let (mut x, mut y) = if polya_aware {
        (polya_trimmed(a, 8), polya_trimmed(b, 8))
    } else {
        (a.to_ascii_uppercase(), b.to_ascii_uppercase())
    };
    if x.len() > y.len() {
        std::mem::swap(&mut x, &mut y);
    }
    if x.len() < 6 {
        return true;
    }
    let q = &x[..x.len().min(20)];
    let ylim = (20i64 + tol + 4).clamp(0, y.len() as i64) as usize;
    let k = (q.len() / 4).max(1) as i32;
    align::distance(q, &y[..ylim], Mode::Hw, k, &[]) != -1
}

/// `_fuzzy_clusters(keys, weight, tol)` (intersect:92) over full-insertion keys
/// (contig, left_pos, right_pos): order by (-weight, contig NAME string, L, R); bucket width
/// `max(tol, 1)` on L; a key joins the lowest-index cluster whose REPRESENTATIVE has
/// |dL| <= tol and |dR| <= tol among buckets b-1, b, b+1 of the key's L bucket (python floor
/// division -- positions are non-negative), else starts a cluster (registered in its own
/// bucket only). Returns clusters as index lists into `keys`, representative first.
pub fn fuzzy_clusters(keys: &[(u32, i64, i64)], weight: &[usize], tol: i64, contigs: &Interner) -> Vec<Vec<usize>> {
    let names: Vec<&'static str> = keys.iter().map(|k| contigs.name(k.0)).collect();
    let mut order: Vec<usize> = (0..keys.len()).collect();
    order.sort_by(|&a, &b| {
        (std::cmp::Reverse(weight[a]), names[a], keys[a].1, keys[a].2).cmp(&(std::cmp::Reverse(weight[b]), names[b], keys[b].1, keys[b].2))
    });
    let w = tol.max(1);
    let mut buckets: FxHashMap<(u32, i64), Vec<usize>> = FxHashMap::default();
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    for ki in order {
        let k = keys[ki];
        let b0 = k.1.div_euclid(w);
        let mut hit: Option<usize> = None;
        for b in [b0 - 1, b0, b0 + 1] {
            if let Some(list) = buckets.get(&(k.0, b)) {
                for &ci in list {
                    let rep = keys[clusters[ci][0]];
                    if (rep.1 - k.1).abs() <= tol && (rep.2 - k.2).abs() <= tol && hit.map_or(true, |h| ci < h) {
                        hit = Some(ci);
                    }
                }
            }
        }
        match hit {
            None => {
                clusters.push(vec![ki]);
                buckets.entry((k.0, b0)).or_default().push(clusters.len() - 1);
            }
            Some(ci) => clusters[ci].push(ki),
        }
    }
    clusters
}

/// `_absorb(target, other, sides)` (intersect:116): member loci of `other` (its member_loci +
/// (f, other.name) for f in other.files, deduplicated in order) not yet in target.member_loci are
/// appended and get `member_sides[m] = sides` (target.member_sides starts from a copy of the
/// existing dict or empty); files of `other` not in target.files are appended. (Mates are not
/// stored.)
pub fn absorb(target: &mut Insertion, other: &Insertion, sides: &[Side]) {
    let set = SideSet(sides.iter().fold(0u8, |a, &s| a | SideSet::only(s).0));
    let mut ms: Vec<(Member, SideSet)> = target.member_sides.take().unwrap_or_default();
    let olocus = other.locus();
    let mut mem: Vec<Member> = Vec::with_capacity(other.member_loci.len() + other.files.len());
    for m in other.member_loci.iter().copied().chain(other.files.iter().map(|&f| (f, olocus))) {
        if !mem.contains(&m) {
            mem.push(m);
        }
    }
    for m in mem {
        if !target.member_loci.contains(&m) {
            target.member_loci.push(m);
            match ms.iter_mut().find(|(k, _)| *k == m) {
                Some(e) => e.1 = set,
                None => ms.push((m, set)),
            }
        }
    }
    target.member_sides = if ms.is_empty() { None } else { Some(ms) };
    // python: `target.files + [f for f in other.files if f not in target.files]` (the filter sees
    // the ORIGINAL target.files)
    let new_files: Vec<FileId> = other.files.iter().copied().filter(|f| !target.files.contains(f)).collect();
    target.files.extend(new_files);
}

/// helper used by absorb / evidence: `(f, name)` member for each file id of `ins`.
pub fn own_members(ins: &Insertion) -> Vec<(FileId, crate::model::LocusKey)> {
    ins.files.iter().map(|&f| (f, ins.locus())).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::insertion::parse_discovery_file;
    use crate::region_filter::filter_dense_regions;
    use std::path::Path;

    /// expected clusters from python `_fuzzy_clusters(keys, weight, tol)` (as indices into keys)
    fn check_clusters(keys: &[(&str, i64, i64)], weight: &[usize], tol: i64, expect: &[&[usize]]) {
        let contigs = Interner::new();
        let k: Vec<(u32, i64, i64)> = keys.iter().map(|(c, l, r)| (contigs.intern(c), *l, *r)).collect();
        let got = fuzzy_clusters(&k, weight, tol, &contigs);
        let want: Vec<Vec<usize>> = expect.iter().map(|c| c.to_vec()).collect();
        assert_eq!(got, want);
    }

    #[test]
    fn shift_tolerant_agree_matches_python() {
        assert_eq!(shift_tolerant_agree(Some(b"GTGGCGAACAGTATTGACCTGGCCGATGCT"), Some(b"ATGGCGAACAGTATTAACCTGGCCTATGCT"), 5, true), true);
        assert_eq!(shift_tolerant_agree(Some(b"ACGGGAGCAGGTCGCCTCAAGATAAGAGTA"), Some(b"ACTGGAGCATCTCGCCTCATGTTAAGAGGA"), 1, true), true);
        assert_eq!(shift_tolerant_agree(Some(b"CACC"), Some(b"CACC"), 5, false), true);
        assert_eq!(shift_tolerant_agree(Some(b"GCAAGGCAGACG"), Some(b"ACAGCAATGACG"), 10, false), false);
        assert_eq!(shift_tolerant_agree(Some(b"CTCT"), Some(b"CTCT"), 10, false), true);
        assert_eq!(shift_tolerant_agree(Some(b"TTTC"), Some(b"AAAAAAAAAAAATTTC"), 10, false), true);
        assert_eq!(shift_tolerant_agree(Some(b"gatgctagttctaaggtgtcggacctacgtgcttgacccacgacg"), Some(b"GATGCTAGTTCTAAGGTGTCGGACCTACGTGCTTGACCCACGACG"), 5, true), true);
        assert_eq!(shift_tolerant_agree(Some(b"TCGGTAAGCTTAAACTTCTTCAGGCGCACC"), Some(b"GGTAAGCTTAAACTTCTTCAGGCGCACCGAGTGC"), 10, true), true);
        assert_eq!(shift_tolerant_agree(Some(b"acacggtgtatgcggacgcacattcgacca"), Some(b"AAAAAAAAAAAAACACGGTGTATGCGGACTCAGCTTAGGCCA"), 10, false), true);
        assert_eq!(shift_tolerant_agree(Some(b"AGCCTAACAACCGGCCCAGCTTCGTTCGAAAATGACTTTCAGAGT"), Some(b"AAAAAAAAAAAAAGCCTAACAACCGGCCCAGCTTCGTTCGAAAATGACTTTCAGAGT"), 5, true), true);
        assert_eq!(shift_tolerant_agree(Some(b"TACCCAGTAGCC"), Some(b"AGATGGTGTTGTTCTTTCACGTCCAAAATGTGTAT"), 10, true), false);
        assert_eq!(shift_tolerant_agree(Some(b"AGCCGCCCTCAG"), Some(b"AAAAAAAAAAAAGGCCTGCCTCTG"), 1, false), false);
        assert_eq!(shift_tolerant_agree(Some(b"ACGAATTTTTAATTTTTCATTTCACCTAGG"), Some(b"GCGAATGTTTACTTTTTCATTTTACCTAGT"), 1, true), true);
        assert_eq!(shift_tolerant_agree(Some(b"ctgcctttccactaacatcactcgccccat"), Some(b"AAAAAAAAAAAAGCCTTTCCACTAACATCACTCGCCCCATTTCACA"), 1, false), false);
        assert_eq!(shift_tolerant_agree(Some(b"cgacggttcggc"), Some(b"GGTTCGGCACTTAA"), 1, true), false);
        assert_eq!(shift_tolerant_agree(Some(b"ACCC"), Some(b"AGTTAACCCCGCCCCGAATATGAACAGTAGCTTCG"), 10, false), true);
        assert_eq!(shift_tolerant_agree(Some(b"ACGTGAGTAATTTGTCGCAGTTAGGAGCTTCACATCTGGCGCCGT"), Some(b"TGAGTAATTTGTCGCAGTTAGGAGCTTCACATCTGGCGCCGTAACACT"), 5, true), true);
        assert_eq!(shift_tolerant_agree(Some(b"ccacgcgagtgcggtcgttaggtgttgact"), Some(b"AAAAAAAAAAAACGACGCGAGTGCGGTCGTTATGTGCTGACT"), 10, false), true);
        assert_eq!(shift_tolerant_agree(Some(b"ATGCTGAGCCGAGAGAAAGCATCTGATAAT"), Some(b"ATGCTGAGCCGAGAGAAAGCATCTGATAAT"), 10, false), true);
        assert_eq!(shift_tolerant_agree(Some(b"CGTA"), Some(b"AGTCGGACGTTCTCCAACTAAATACAGGTTCACCG"), 10, true), true);
        assert_eq!(shift_tolerant_agree(Some(b"ATCACACAATAT"), Some(b"AACGGACTCTAT"), 10, true), false);
        assert_eq!(shift_tolerant_agree(Some(b"AAAC"), Some(b"AACGGGAGT"), 10, false), true);
        assert_eq!(shift_tolerant_agree(Some(b"atga"), Some(b"GTGA"), 5, false), true);
        assert_eq!(shift_tolerant_agree(Some(b"ctgctgaaatgc"), Some(b"CGGCCGATATAA"), 1, false), false);
        assert!(shift_tolerant_agree(None, Some(b"ACGT"), 5, false));
    }
    #[test]
    fn fuzzy_clusters_match_python() {
        check_clusters(&[("chr10", 149, 164), ("chr10", 141, 152), ("chr10", 124, 149), ("chr2", 147, 148), ("chr2", 142, 170), ("chr10", 144, 151), ("chr1", 138, 152), ("chr1", 128, 134), ("chr1", 118, 132), ("chr10", 153, 159), ("chr2", 147, 157), ("chr10", 104, 116), ("chr1", 143, 148), ("chr1", 123, 138), ("chr10", 104, 119), ("chr2", 132, 159), ("chrX", 143, 149), ("chr10", 112, 138), ("chrX", 112, 121), ("chrX", 117, 124), ("chr2", 102, 115), ("chr10", 121, 134), ("chr1", 136, 147), ("chr10", 115, 141), ("chr1", 109, 128), ("chr2", 138, 152), ("chrX", 135, 152), ("chrX", 108, 116), ("chr10", 135, 138), ("chr2", 126, 130), ("chr10", 133, 137), ("chr2", 156, 180), ("chr1", 110, 117), ("chrX", 110, 112), ("chrX", 150, 163), ("chr2", 156, 174), ("chr10", 155, 159), ("chr2", 145, 158), ("chr1", 103, 116), ("chr1", 101, 129), ("chr2", 104, 113), ("chrX", 104, 120), ("chrX", 154, 163), ("chr1", 128, 135), ("chrX", 142, 158), ("chr2", 157, 173), ("chr10", 127, 129), ("chr2", 136, 148), ("chr10", 155, 177), ("chr2", 141, 148), ("chrX", 123, 139), ("chr2", 143, 169), ("chr1", 144, 167), ("chr1", 139, 160), ("chrX", 113, 134), ("chr2", 151, 180), ("chr1", 128, 143), ("chr2", 143, 167), ("chr10", 129, 139)], &[2, 4, 1, 2, 4, 4, 2, 2, 3, 3, 4, 4, 3, 2, 2, 2, 3, 1, 1, 2, 4, 4, 1, 4, 4, 3, 3, 3, 4, 3, 2, 4, 1, 2, 4, 3, 1, 3, 2, 2, 2, 3, 3, 3, 2, 1, 4, 1, 2, 1, 4, 3, 1, 1, 1, 2, 1, 2, 1], 6, &[&[24], &[11, 14], &[23, 17], &[21, 46], &[28, 30, 58], &[1, 5], &[20, 40], &[4, 51, 57], &[10, 37], &[31, 35, 55], &[50], &[34, 42], &[8, 13], &[43, 7], &[12, 6], &[9, 0, 36], &[29], &[25, 47, 49], &[41, 27], &[26], &[16], &[39], &[38], &[48], &[15], &[3], &[33], &[19, 18], &[44], &[32], &[56], &[22], &[53], &[52], &[2], &[45], &[54]]);
        check_clusters(&[("chr10", 111, 119), ("chr10", 101, 101), ("chr1", 105, 134), ("chr1", 112, 116), ("chrX", 121, 123), ("chr2", 120, 129), ("chrX", 147, 162), ("chr2", 121, 122), ("chr1", 116, 121), ("chr2", 105, 107), ("chr1", 144, 152), ("chr10", 150, 177), ("chr2", 121, 137), ("chrX", 109, 115), ("chr1", 148, 152), ("chrX", 124, 133), ("chr1", 114, 123), ("chr1", 151, 166), ("chr1", 104, 122), ("chr10", 112, 137), ("chrX", 151, 165), ("chr10", 139, 141), ("chrX", 136, 149), ("chr10", 100, 106), ("chr10", 106, 132), ("chrX", 115, 139), ("chr2", 132, 145), ("chr2", 146, 147), ("chr1", 114, 137), ("chr1", 114, 130), ("chr2", 113, 133), ("chrX", 139, 145), ("chr10", 113, 122), ("chr2", 108, 113), ("chr1", 114, 128), ("chr2", 152, 174), ("chr2", 125, 135), ("chr2", 103, 127), ("chr2", 105, 114), ("chr1", 120, 136), ("chr10", 109, 114), ("chr10", 129, 129), ("chr10", 120, 123), ("chr2", 143, 165), ("chrX", 133, 142), ("chr1", 106, 127), ("chr1", 139, 151), ("chrX", 130, 132), ("chr2", 151, 172), ("chr10", 128, 138), ("chrX", 145, 158), ("chr2", 134, 148), ("chr2", 139, 140), ("chr1", 149, 163), ("chr1", 140, 169), ("chr2", 108, 109), ("chr10", 104, 118), ("chr1", 119, 140), ("chr1", 154, 178), ("chr2", 127, 143)], &[1, 2, 4, 1, 1, 1, 3, 2, 1, 1, 3, 2, 4, 2, 2, 2, 4, 4, 3, 3, 1, 2, 4, 1, 1, 3, 4, 4, 2, 2, 3, 4, 4, 2, 2, 2, 4, 1, 3, 2, 1, 3, 4, 2, 3, 3, 2, 3, 2, 4, 1, 1, 2, 3, 1, 3, 1, 1, 3, 2], 12, &[&[2, 16, 18, 45, 34, 29, 28], &[17, 53, 58, 54], &[32, 42, 56, 24, 40, 0], &[49, 41, 21], &[12, 36, 26, 30, 59, 5], &[27, 52, 51], &[22, 31, 44, 50], &[10, 46, 14], &[19], &[38, 55, 33, 9], &[25, 15], &[47, 4], &[6, 20], &[39, 57], &[1, 23], &[11], &[7], &[43, 48, 35], &[13], &[3, 8], &[37]]);
        check_clusters(&[("chr2", 123, 132), ("chr2", 139, 150), ("chrX", 124, 133), ("chr1", 114, 114), ("chrX", 148, 168), ("chr10", 152, 181), ("chr1", 156, 179), ("chr10", 148, 152), ("chr2", 116, 132), ("chr2", 124, 137), ("chr2", 108, 115), ("chr2", 142, 168), ("chr1", 122, 150), ("chr10", 154, 164), ("chr10", 154, 177), ("chr1", 150, 177), ("chrX", 121, 136), ("chrX", 150, 173), ("chr10", 146, 156), ("chr2", 115, 117), ("chr1", 107, 117), ("chr1", 157, 182), ("chr1", 114, 125), ("chr1", 139, 141), ("chrX", 147, 148), ("chr10", 155, 169), ("chrX", 119, 144), ("chrX", 120, 148), ("chr2", 146, 172), ("chr2", 147, 174), ("chr2", 136, 165), ("chr1", 138, 156), ("chr1", 130, 144), ("chrX", 100, 128), ("chr10", 113, 119), ("chr2", 134, 145), ("chr1", 141, 170), ("chr1", 129, 147), ("chrX", 101, 123), ("chr10", 127, 129), ("chr10", 133, 142), ("chr2", 106, 113), ("chrX", 110, 122), ("chr1", 159, 172), ("chr10", 120, 129), ("chr2", 132, 155), ("chr10", 131, 148), ("chr1", 142, 169), ("chr10", 138, 150), ("chr10", 111, 111), ("chr1", 155, 173), ("chr2", 103, 132), ("chr1", 113, 129), ("chr1", 157, 173), ("chr10", 132, 146), ("chr10", 135, 141), ("chr10", 109, 129), ("chrX", 151, 151)], &[4, 2, 3, 3, 2, 4, 2, 4, 1, 1, 1, 3, 2, 2, 3, 2, 2, 2, 2, 2, 1, 4, 2, 3, 4, 1, 4, 1, 4, 1, 1, 4, 2, 3, 4, 2, 2, 3, 4, 2, 2, 2, 2, 4, 3, 4, 3, 3, 2, 2, 4, 1, 2, 2, 3, 1, 3, 2], 1, &[&[31], &[50], &[21], &[43], &[34], &[7], &[5], &[0], &[45], &[28], &[38], &[26], &[24], &[3], &[37], &[23], &[47, 36], &[56], &[44], &[46], &[54], &[14], &[11], &[33], &[2], &[52], &[22], &[12], &[32], &[15], &[6], &[53], &[49], &[39], &[40], &[48], &[18], &[13], &[41], &[19], &[35], &[1], &[42], &[16], &[4], &[17], &[57], &[20], &[55], &[25], &[51], &[10], &[8], &[9], &[30], &[29], &[27]]);
        check_clusters(&[("chrX", 153, 167), ("chrX", 130, 138), ("chrX", 133, 139), ("chrX", 137, 153), ("chr10", 132, 137), ("chr10", 104, 115), ("chrX", 104, 116), ("chr1", 122, 145), ("chrX", 121, 132), ("chrX", 141, 145), ("chrX", 155, 181), ("chr1", 102, 129), ("chrX", 122, 138), ("chrX", 127, 146), ("chr2", 110, 127), ("chr1", 143, 147), ("chr2", 143, 170), ("chrX", 150, 160), ("chr10", 121, 146), ("chr10", 135, 152), ("chrX", 141, 146), ("chr2", 107, 111), ("chr1", 139, 149), ("chrX", 128, 143), ("chr2", 123, 139), ("chr1", 122, 139), ("chr2", 140, 155), ("chr1", 121, 129), ("chrX", 139, 158), ("chr2", 101, 112), ("chrX", 104, 115), ("chr1", 117, 145), ("chr2", 118, 144), ("chrX", 110, 132), ("chrX", 101, 103), ("chr10", 113, 114), ("chr10", 109, 118), ("chr10", 114, 115), ("chrX", 116, 119), ("chr1", 109, 126), ("chr1", 149, 178), ("chr10", 127, 153), ("chr10", 102, 125), ("chrX", 154, 177), ("chrX", 127, 129), ("chr10", 138, 142), ("chr2", 102, 104), ("chr1", 110, 113), ("chr1", 101, 111), ("chr10", 107, 121), ("chr10", 106, 111), ("chr10", 138, 149), ("chr10", 123, 126), ("chrX", 120, 132), ("chrX", 116, 130), ("chr10", 130, 130), ("chr10", 110, 115), ("chr10", 150, 161), ("chr1", 128, 144), ("chr1", 150, 164)], &[1, 4, 4, 1, 3, 4, 2, 1, 2, 4, 2, 4, 2, 1, 1, 3, 4, 2, 4, 4, 3, 4, 2, 3, 4, 2, 3, 2, 1, 3, 3, 3, 3, 2, 4, 3, 1, 4, 1, 2, 4, 1, 4, 3, 4, 1, 1, 2, 1, 4, 3, 1, 4, 4, 3, 1, 4, 3, 1, 1], 1, &[&[11], &[40], &[42], &[5], &[49], &[56], &[37, 35], &[18], &[52], &[19], &[21], &[24], &[16], &[34], &[53, 8], &[44], &[1], &[2], &[9, 20], &[31], &[15], &[50], &[4], &[57], &[29], &[32], &[26], &[30, 6], &[54], &[23], &[43], &[39], &[47], &[27], &[25], &[22], &[33], &[12], &[17], &[10], &[48], &[7], &[58], &[59], &[36], &[41], &[55], &[45], &[51], &[46], &[14], &[38], &[13], &[3], &[28], &[0]]);
        check_clusters(&[("chr2", 113, 115), ("chr2", 117, 142), ("chr2", 113, 142), ("chrX", 149, 167), ("chr2", 129, 149), ("chr2", 125, 146), ("chrX", 107, 108), ("chr10", 151, 172), ("chr2", 103, 122), ("chr10", 122, 142), ("chrX", 154, 161), ("chr2", 152, 168), ("chr1", 128, 143), ("chr1", 105, 107), ("chr1", 113, 127), ("chrX", 156, 178), ("chr1", 146, 155), ("chr2", 153, 182), ("chr10", 108, 128), ("chr1", 141, 146), ("chr2", 121, 126), ("chr10", 158, 187), ("chr10", 130, 157), ("chr10", 116, 124), ("chr1", 114, 119), ("chr2", 149, 151), ("chrX", 134, 153), ("chrX", 113, 116), ("chrX", 158, 173), ("chr2", 143, 144), ("chrX", 114, 134), ("chrX", 130, 156), ("chr10", 159, 167), ("chr10", 133, 154), ("chr1", 135, 145), ("chrX", 156, 161), ("chr10", 157, 172), ("chrX", 131, 160), ("chr2", 136, 147), ("chr1", 135, 150), ("chr2", 110, 120), ("chr1", 123, 135), ("chr1", 108, 123), ("chr2", 121, 133), ("chr10", 120, 144), ("chr1", 120, 126), ("chrX", 107, 116), ("chrX", 140, 151), ("chr2", 130, 159), ("chr10", 134, 161), ("chr10", 123, 129), ("chr10", 119, 128), ("chr10", 145, 163), ("chr1", 126, 126), ("chr10", 135, 137), ("chr10", 132, 148), ("chr1", 148, 174), ("chr10", 142, 145), ("chr2", 159, 162), ("chr10", 143, 161)], &[1, 3, 1, 4, 1, 3, 3, 1, 4, 3, 2, 1, 2, 2, 2, 1, 2, 1, 3, 3, 4, 4, 1, 1, 4, 1, 3, 2, 4, 3, 1, 1, 1, 4, 4, 2, 3, 3, 2, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 2, 3, 1, 4, 4, 4, 1, 1, 2, 4, 2], 6, &[&[24, 42], &[53, 45], &[34, 39, 19], &[33, 22, 55], &[54], &[52, 59], &[21], &[8], &[20], &[58], &[3, 10], &[28, 15], &[41], &[18], &[9, 44], &[50, 51], &[36, 7, 32], &[40, 0], &[1, 2], &[5, 4], &[29], &[6], &[37, 31], &[26, 47], &[13], &[14], &[12], &[16], &[49], &[57], &[43], &[38], &[46, 27], &[35], &[56], &[23], &[48], &[25], &[11], &[17], &[30]]);
    }

    // ------------------------------------------------------------------ fixture (3 colonies)

    fn fixture_dir() -> Option<std::path::PathBuf> {
        let p = std::env::var("P1_FIXTURE_DIR").unwrap_or_else(|_| {
            "/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad/work/e2e".into()
        });
        let p = std::path::PathBuf::from(p);
        if p.join("S1.discovery.txt.gz").exists() {
            Some(p)
        } else {
            eprintln!("skipping: fixture {} not found", p.display());
            None
        }
    }

    fn fnv(s: &str) -> u64 {
        let mut h: u64 = 0xcbf29ce484222325;
        for b in s.bytes() {
            h = (h ^ b as u64).wrapping_mul(0x100000001b3);
        }
        h
    }

    fn qstr(q: &Option<crate::seq::QualSeq>) -> String {
        match q {
            None => "-".into(),
            Some(q) => {
                let quals: Vec<String> = q.qual.iter().map(|x| x.to_string()).collect();
                format!("{}/{}", String::from_utf8_lossy(&q.seq), quals.join(","))
            }
        }
    }

    fn opt<T: std::fmt::Display>(x: Option<T>) -> String {
        x.map(|v| v.to_string()).unwrap_or_else(|| "None".into())
    }

    /// the canonical per-insertion line the python dump (scratch p1dump.py) hashes
    fn canon(i: &Insertion, contigs: &Interner, files: &[InputFile]) -> String {
        let name = i.name(contigs);
        let ms = match &i.member_sides {
            None => String::new(),
            Some(v) => v
                .iter()
                .map(|((f, l), s)| {
                    let mut sides = vec![];
                    if s.contains(Side::Left) {
                        sides.push("LEFT");
                    }
                    if s.contains(Side::Right) {
                        sides.push("RIGHT");
                    }
                    format!("{}@{}={}", files[*f as usize].basename, l.name(contigs), sides.join(","))
                })
                .collect::<Vec<_>>()
                .join(";"),
        };
        [
            name,
            (i.ty as u8).to_string(),
            match i.open_side {
                None => "None".into(),
                Some(s) => s.as_str().into(),
            },
            opt(i.left_pos),
            opt(i.right_pos),
            i.files.iter().map(|&f| files[f as usize].basename.clone()).collect::<Vec<_>>().join(","),
            i.member_loci.iter().map(|(f, l)| format!("{}@{}", files[*f as usize].basename, l.name(contigs))).collect::<Vec<_>>().join(";"),
            ms,
            qstr(&i.left_clipped),
            qstr(&i.left_aligned),
            qstr(&i.right_clipped),
            qstr(&i.right_aligned),
        ]
        .join("|")
    }

    fn test_cfg(tol: i64, polya_aware: bool, keep_polya: bool) -> Config {
        let j = serde_json::json!({
            "combine_insertions": {
                "genome_2bit": "g.2bit", "exclude_files_with_many_insertions": 100000,
                "samtools_executable": "samtools", "bowtie2_executable": "bowtie2",
                "bowtie2_index": "i", "bowtie2_index2": "i2", "bowtie2_index2_lo": "lo",
                "clean_remap_max_insertion": 20,
                "keep_polya_one_sided": keep_polya, "merge_tolerance_bp": tol,
                "polya_aware_clip_agreement": polya_aware
            },
            "genotyping": {"max_bases": 50}
        });
        Config::from_json(&j, None).unwrap()
    }

    /// expected counts and hashes from python (scratch/p1dump.py: parseFile -> contig filter ->
    /// intersect_insertions -> filter_dense_regions(100, cut)) over the same three files.
    #[test]
    fn fixture_matches_python() {
        let Some(dir) = fixture_dir() else { return };
        let paths: Vec<String> = (1..=3).map(|i| dir.join(format!("S{i}.discovery.txt.gz")).to_string_lossy().into_owned()).collect();
        let files: Vec<InputFile> = paths.iter().map(|p| InputFile::new(p)).collect();
        let contigs = Interner::new();
        use rayon::prelude::*;
        let imports: Vec<_> = files
            .par_iter()
            .enumerate()
            .map(|(n, f)| parse_discovery_file(Path::new(&f.path), n as u32, &contigs).unwrap())
            .collect();
        let mut raw = vec![];
        for imp in imports {
            raw.extend(imp.records);
        }
        assert_eq!(raw.len(), 1758);
        // (tol, polya_aware, keep_polya, rf_cutoff, n_intersect, hash, n_rf_kept, removed, hot, rf_hash)
        let cases: &[(i64, bool, bool, usize, usize, u64, usize, usize, usize, u64)] = &[
            (0, false, false, 4, 1449, 0x678bfe0aa2ba52ed, 1429, 20, 4, 0xaa92dc72dc1b166a),
            (0, true, false, 2, 1446, 0xd13159a3b9ab4f0f, 1010, 436, 203, 0xe75aa72293ef5709),
            (5, false, false, 2, 1441, 0x75d42a7aa62d2802, 1018, 423, 197, 0xf02edcc0b0ad8615),
            (5, true, true, 2, 1438, 0xc5dccd10dcfde531, 1019, 419, 194, 0x38e1a60bf9b3e7ea),
            (10, true, true, 4, 1429, 0x5d2fd3a666004c3a, 1413, 16, 3, 0xdad1695c568de7d7),
            (0, false, true, 4, 1451, 0x479950e8afbdc196, 1431, 20, 4, 0xafd5a5a27e0178ad),
            (3, false, true, 2, 1445, 0x854a39faadf91786, 1016, 429, 199, 0x77ff5e4636b812d5),
        ];
        for &(tol, pa, kp, cut, n, h, nk, rem, hot, hk) in cases {
            let cfg = test_cfg(tol, pa, kp);
            let out = intersect_insertions(raw.clone(), &cfg, &contigs, &files);
            let text = out.iter().map(|i| canon(i, &contigs, &files)).collect::<Vec<_>>().join("\n");
            assert_eq!(out.len(), n, "intersect n tol={tol} pa={pa} kp={kp}");
            assert_eq!(fnv(&text), h, "intersect hash tol={tol} pa={pa} kp={kp}");
            assert!(out.iter().enumerate().all(|(i, x)| x.uid as usize == i));
            let (kept, removed, nhot) = filter_dense_regions(out, 100, cut);
            let ktext = kept.iter().map(|i| canon(i, &contigs, &files)).collect::<Vec<_>>().join("\n");
            assert_eq!((kept.len(), removed, nhot), (nk, rem, hot), "rf tol={tol} pa={pa} kp={kp}");
            assert_eq!(fnv(&ktext), hk, "rf hash tol={tol} pa={pa} kp={kp}");
        }
    }
}
