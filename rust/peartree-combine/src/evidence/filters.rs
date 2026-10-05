//! Evidence-level TPRT filters (far pair, slippage) and their helpers. OWNER: P4.
//!
//! Mirrors `discovery_breakpoints`, `_colonies`, `_far_pair_check`, `_inside_mates`,
//! `_to_one_sided`, `_carries_element`, `_slippage_check`, `_outward_clip`, `_aligned_part` of
//! src/combine_insertions_evidence.py:1294-1452 and 358-363. SPEC.md §4.5, §4.6.

use crate::config::Config;
use crate::evidence::junction::JunctionRecord;
use crate::evidence::row::{allele_forward_seq, Role};
use crate::genome::RefFetch;
use crate::library::Matcher;
use crate::model::{ContigId, FileId, InsType, Insertion, Interner, Side, Tok};
use crate::seq::revcomp;
use crate::tprt::{far_geometry, far_pair_verdict, outward_reference, slippage_junction, FarPairInput};
use rustc_hash::FxHashMap;

/// `discovery_breakpoints(records)`: {(contig, side): [(pos, sample)]} over EVERY accepted
/// per-sample discovery record before intersect, in record order. LEFT when type not in (2, 5)
/// and left_pos is set; RIGHT when type not in (1, 4) and right_pos is set. Sample = files[0].
pub type Breakpoints = FxHashMap<(ContigId, Side), Vec<(i64, FileId)>>;

/// Sample id used for a record without files (python `"?"`; never happens for parsed records).
pub const UNKNOWN_SAMPLE: FileId = FileId::MAX;

pub fn discovery_breakpoints(records: &[Insertion]) -> Breakpoints {
    let mut bp: Breakpoints = FxHashMap::default();
    for i in records {
        let s = i.files.first().copied().unwrap_or(UNKNOWN_SAMPLE);
        if !matches!(i.ty, InsType::LeftPolyA | InsType::LeftDisc) {
            if let Some(p) = i.left_pos {
                bp.entry((i.contig, Side::Left)).or_default().push((p, s));
            }
        }
        if !matches!(i.ty, InsType::RightPolyA | InsType::RightDisc) {
            if let Some(p) = i.right_pos {
                bp.entry((i.contig, Side::Right)).or_default().push((p, s));
            }
        }
    }
    bp
}

/// `_aligned_part(ins, side)`: RIGHT -> `str(right_aligned.revcomp()).upper()`; LEFT ->
/// `str(left_aligned).upper()`; "" when absent.
pub fn aligned_part(ins: &Insertion, side: Side) -> Vec<u8> {
    match side {
        Side::Right => ins.right_aligned.as_ref().map(|a| revcomp(&a.seq).to_ascii_uppercase()).unwrap_or_default(),
        Side::Left => ins.left_aligned.as_ref().map(|a| a.seq.to_ascii_uppercase()).unwrap_or_default(),
    }
}

/// `_outward_clip(rec, ins, with_mates=False)`: with_mates -> rec.consensus, else
/// combined_consensus if its seq is non-empty else consensus; non-empty -> its seq uppercased;
/// else the insertion's clipped seq of rec.side uppercased; else "".
pub fn outward_clip(rec: &JunctionRecord, ins: &Insertion, with_mates: bool) -> Vec<u8> {
    let c = if with_mates || rec.combined_consensus.seq.is_empty() { &rec.consensus } else { &rec.combined_consensus };
    if !c.seq.is_empty() {
        return c.seq.to_ascii_uppercase();
    }
    ins.clipped(rec.side).map(|q| q.seq.to_ascii_uppercase()).unwrap_or_default()
}

/// `_inside_mates(rec, min_mapq=20, max_dist=1000)`: needs `rec.rows` (not detached).
pub fn inside_mates(rec: &JunctionRecord, contigs: &Interner) -> Vec<Vec<u8>> {
    const MIN_MAPQ: i64 = 20;
    const MAX_DIST: i64 = 1000;
    let rows = rec.rows.as_ref().expect("_inside_mates needs the record's rows (not detached)");
    // ev: (sample, frag) -> LAST CLIP/DISC row (python dict assignment)
    let mut ev: FxHashMap<(FileId, &str), usize> = FxHashMap::default();
    for (i, r) in rows.iter().enumerate() {
        if matches!(r.role, Role::Clip | Role::Disc) {
            ev.insert((r.file, &r.frag), i);
        }
    }
    let mut seen: FxHashMap<(FileId, &str), ()> = FxHashMap::default();
    let mut out = Vec::new();
    for r in rows {
        let k = (r.file, &*r.frag);
        if r.role != Role::Mate || seen.contains_key(&k) || r.seq.is_empty() || &*r.seq == b"*" {
            continue;
        }
        let Some(&pi) = ev.get(&k) else { continue };
        let p = &rows[pi];
        let inside = !r.mapped(contigs) || r.ref_ != p.ref_ || (r.pos - p.pos).abs() > MAX_DIST || r.mapq < MIN_MAPQ;
        if inside {
            let s = allele_forward_seq(r).to_ascii_uppercase();
            seen.insert(k, ());
            out.push(if rec.side == Side::Left { revcomp(&s) } else { s });
        }
    }
    out
}

/// `_colonies(ins, rec, breakpoints, tol)`: samples of rec's evidence-role rows plus every
/// discovery breakpoint sample within tol of the insertion's junction. Sorted unique FileIds.
pub fn colonies(ins: &Insertion, rec: &JunctionRecord, bps: Option<&Breakpoints>, tol: i64) -> Vec<u32> {
    let pos = ins.junction(rec.side).expect("_colonies: junction of a real side");
    let rows = rec.rows.as_ref().expect("_colonies needs the record's rows (not detached)");
    let mut out: Vec<u32> = rows.iter().filter(|r| r.role.is_evidence()).map(|r| r.file).collect();
    if let Some(v) = bps.and_then(|b| b.get(&(ins.contig, rec.side))) {
        for &(p, s) in v {
            if (p - pos).abs() <= tol {
                out.push(s);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// `_far_pair_check(ins, recs, cfg, matcher, breakpoints, ref_fetch)` -> None (not a far pair /
/// passes) or Some((reason, polya side)). `recs` = [LEFT, RIGHT] records (both present: the
/// caller checks `open_side is None and len(recs) == 2`). Reasons are the python strings
/// ("no_polarity", "no_element_on_complex_clip", "element_antisense", "element_class_conflict",
/// "colony_mismatch", "polya_side_slippage"); "few_fragments" is never produced (gate dropped).
pub fn far_pair_check(
    ins: &Insertion,
    recs: &[JunctionRecord],
    cfg: &Config,
    matcher: &dyn Matcher,
    bps: Option<&Breakpoints>,
    ref_fetch: Option<&dyn RefFetch>,
    contigs: &Interner,
) -> Option<(String, Option<Side>)> {
    let gap = ins.right_pos.expect("far pair: right_pos") - ins.left_pos.expect("far pair: left_pos");
    if !far_geometry(gap, cfg) {
        return None;
    }
    let tol = cfg.merge_tolerance_bp.max(cfg.far_pair_colony_tol);
    // python `by = {r.side: r for r in recs}`: the last record of a side wins
    let mut by: [Option<&JunctionRecord>; 2] = [None, None];
    for r in recs {
        by[r.side as usize] = Some(r);
    }
    let mut inp = FarPairInput { clips: [vec![], vec![]], colonies: [vec![], vec![]], inside_mates: [vec![], vec![]] };
    for (si, rec) in by.iter().enumerate() {
        let Some(rec) = rec else { continue };
        let old = ins.clipped(rec.side).map(|q| q.seq.to_ascii_uppercase()).unwrap_or_default();
        let mut clips: Vec<Vec<u8>> = Vec::with_capacity(2);
        for x in [outward_clip(rec, ins, false), old] {
            if !x.is_empty() && !clips.contains(&x) {
                clips.push(x);
            }
        }
        inp.clips[si] = clips;
        inp.colonies[si] = colonies(ins, rec, bps, tol);
        inp.inside_mates[si] = inside_mates(rec, contigs);
    }
    let (mut reason, pside) = far_pair_verdict(&inp, matcher, cfg);
    if reason.is_empty() {
        if let Some(fetch) = ref_fetch {
            let ps = pside.expect("a passing verdict has a poly-A side");
            let jn = ins.junction(ps).expect("far pair: junction");
            let (line, j) = outward_reference(fetch, contigs.name(ins.contig), jn, ps, 80);
            let rec = by[ps as usize].expect("far pair: poly-A side record");
            if !slippage_junction(&outward_clip(rec, ins, false), &line, j, cfg).is_empty() {
                reason = "polya_side_slippage";
            }
        }
    }
    if reason.is_empty() {
        None
    } else {
        Some((reason.to_string(), pside))
    }
}

/// `_to_one_sided(ins, real_side)`: drop the other side (clipped/aligned None, mates cleared),
/// LEFT real: right_pos = left_pos, type 4, open_side RIGHT, name `c:{L}-oneside_{L}`;
/// RIGHT real: left_pos = right_pos, type 5, open_side LEFT, name `c:oneside_{R}-{R}`.
pub fn to_one_sided(ins: &mut Insertion, real_side: Side) {
    match real_side {
        Side::Left => {
            let l = ins.left_pos.expect("_to_one_sided: left_pos");
            ins.right_clipped = None;
            ins.right_aligned = None;
            ins.right_pos = Some(l);
            ins.ty = InsType::RightDisc;
            ins.open_side = Some(Side::Right);
            ins.name_start = Tok::pos(l);
            ins.name_end = Tok::one_side(l);
        }
        Side::Right => {
            let r = ins.right_pos.expect("_to_one_sided: right_pos");
            ins.left_clipped = None;
            ins.left_aligned = None;
            ins.left_pos = Some(r);
            ins.ty = InsType::LeftDisc;
            ins.open_side = Some(Side::Left);
            ins.name_start = Tok::one_side(r);
            ins.name_end = Tok::pos(r);
        }
    }
}

/// `_carries_element(rec, ins, matcher)`.
pub fn carries_element(rec: &JunctionRecord, ins: &Insertion, matcher: &dyn Matcher, contigs: &Interner) -> bool {
    if matcher.hit(&outward_clip(rec, ins, false), true).is_some() {
        return true;
    }
    let mut hits = 0;
    for m in inside_mates(rec, contigs) {
        if matcher.hit(&m, true).is_some() || matcher.hit(&revcomp(&m), true).is_some() {
            hits += 1;
            if hits >= 2 {
                return true;
            }
        }
    }
    false
}

/// `_slippage_check(ins, recs, cfg, matcher, ref_fetch)` -> "" or `slippage:{SIDE}({why})`.
pub fn slippage_check(
    ins: &Insertion,
    recs: &[JunctionRecord],
    cfg: &Config,
    matcher: &dyn Matcher,
    ref_fetch: &dyn RefFetch,
    contigs: &Interner,
) -> String {
    // python dict `slip` keyed by side: insertion order = first record of a side, value = last
    let mut slip: Vec<(Side, &'static str)> = Vec::with_capacity(2);
    for r in recs {
        let jn = ins.junction(r.side).expect("_slippage_check: junction of a real side");
        let (line, j) = outward_reference(ref_fetch, contigs.name(ins.contig), jn, r.side, 80);
        let oc = outward_clip(r, ins, false);
        let mut why = slippage_junction(&oc, &line, j, cfg);
        if !why.is_empty()
            && !r.consensus.seq.is_empty()
            && r.consensus.seq.len() > oc.len()
            && slippage_junction(&r.consensus.seq.to_ascii_uppercase(), &line, j, cfg).is_empty()
        {
            why = "";
        }
        match slip.iter_mut().find(|(s, _)| *s == r.side) {
            Some(e) => e.1 = why,
            None => slip.push((r.side, why)),
        }
    }
    let slip_of = |s: Side| slip.iter().find(|(x, _)| *x == s).map(|(_, w)| *w).unwrap_or("");
    for &(s, why) in &slip {
        if why.is_empty() {
            continue;
        }
        let carried = recs
            .iter()
            .filter(|o| o.side != s && slip_of(o.side).is_empty())
            .any(|o| carries_element(o, ins, matcher, contigs));
        if !carried {
            return format!("slippage:{}({})", s.as_str(), why);
        }
    }
    String::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{LocusKey, TokKind};
    use crate::seq::QualSeq;

    fn qs(s: &[u8]) -> QualSeq {
        QualSeq { seq: s.to_vec().into_boxed_slice(), qual: vec![30u8; s.len()].into_boxed_slice() }
    }

    fn ins(ty: InsType, l: Option<i64>, r: Option<i64>, files: Vec<FileId>) -> Insertion {
        Insertion {
            uid: 0,
            contig: 0,
            name_start: Tok::pos(l.unwrap_or(0)),
            name_end: Tok::pos(r.unwrap_or(0)),
            ty,
            open_side: None,
            left_clipped: None,
            left_aligned: None,
            left_pos: l,
            right_clipped: None,
            right_aligned: None,
            right_pos: r,
            files,
            member_loci: vec![],
            member_sides: None,
        }
    }

    #[test]
    fn breakpoints_by_type() {
        let recs = vec![
            ins(InsType::FullInfo, Some(10), Some(20), vec![1]),
            ins(InsType::RightPolyA, Some(30), None, vec![2]),
            ins(InsType::LeftPolyA, None, Some(40), vec![3]),
            ins(InsType::RightDisc, Some(50), Some(55), vec![4]),
            ins(InsType::LeftDisc, Some(60), Some(65), vec![]),
        ];
        let bp = discovery_breakpoints(&recs);
        assert_eq!(bp[&(0, Side::Left)], vec![(10, 1), (30, 2), (50, 4)]);
        assert_eq!(bp[&(0, Side::Right)], vec![(20, 1), (40, 3), (65, UNKNOWN_SAMPLE)]);
    }

    #[test]
    fn one_sided_names_and_aligned_part() {
        let c = Interner::new();
        c.intern("chr1");
        let mut i = ins(InsType::FullInfo, Some(100), Some(115), vec![0]);
        i.left_aligned = Some(qs(b"acgT"));
        i.right_aligned = Some(qs(b"aacG"));
        i.right_clipped = Some(qs(b"TTT"));
        assert_eq!(aligned_part(&i, Side::Left), b"ACGT".to_vec());
        assert_eq!(aligned_part(&i, Side::Right), b"CGTT".to_vec());
        let mut a = i.clone();
        to_one_sided(&mut a, Side::Left);
        assert_eq!(a.name(&c), "chr1:100-oneside_100");
        assert_eq!((a.ty, a.open_side, a.right_pos, a.right_clipped.is_none()), (InsType::RightDisc, Some(Side::Right), Some(100), true));
        assert_eq!(a.open_side_eff(), Some(Side::Right));
        let mut b = i.clone();
        to_one_sided(&mut b, Side::Right);
        assert_eq!(b.name(&c), "chr1:oneside_115-115");
        assert_eq!(b.name_start.kind, TokKind::OneSide);
        assert_eq!((b.ty, b.left_pos, b.left_aligned.is_none()), (InsType::LeftDisc, Some(115), true));
        assert_eq!(LocusKey::parse(&b.name(&c), &c), Some(b.locus()));
        assert_eq!(aligned_part(&b, Side::Left), Vec::<u8>::new());
    }

    fn row(file: FileId, side: Side, role: Role, frag: &str) -> crate::evidence::row::EvidenceRow {
        crate::evidence::row::EvidenceRow {
            file,
            locus: LocusKey { contig: 0, start: Tok::pos(100), end: Tok::pos(200) },
            side,
            role,
            frag: frag.into(),
            r12: 1,
            flag: 0,
            ref_: 0,
            pos: 100,
            strand: "+".into(),
            outer: -1,
            mref: 0,
            mpos: -1,
            mstrand: "*".into(),
            tlen: 0,
            mapq: 60,
            cigar: "*".into(),
            clip_at: -1,
            seq: b"ACGT".to_vec().into_boxed_slice(),
            qual: Box::new([]),
        }
    }

    fn rec(side: Side, cons: &[u8], rows: Vec<crate::evidence::row::EvidenceRow>) -> JunctionRecord {
        use crate::consensus::ConsensusResult;
        use crate::evidence::junction::Supported;
        JunctionRecord {
            insertion_id: "chr1:100-200".into(),
            side,
            rows: Some(rows),
            reads_ref: None,
            n_reads: 0,
            n_fragments: 0,
            n_independent: 0,
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
            consensus: ConsensusResult { seq: cons.to_vec(), ..ConsensusResult::default() },
            combined_consensus: ConsensusResult::default(),
            polya_end: false,
            fail_reason: String::new(),
            aligned: vec![],
        }
    }

    #[test]
    fn outward_clip_fallbacks() {
        let mut i = ins(InsType::FullInfo, Some(100), Some(200), vec![0]);
        i.left_clipped = Some(qs(b"ggcc"));
        let mut r = rec(Side::Left, b"acg", vec![]);
        assert_eq!(outward_clip(&r, &i, false), b"ACG".to_vec());
        r.combined_consensus.seq = b"tt".to_vec();
        assert_eq!(outward_clip(&r, &i, false), b"TT".to_vec());
        assert_eq!(outward_clip(&r, &i, true), b"ACG".to_vec());
        r.consensus.seq.clear();
        r.combined_consensus.seq.clear();
        assert_eq!(outward_clip(&r, &i, false), b"GGCC".to_vec());
        let r2 = rec(Side::Right, b"", vec![]);
        assert_eq!(outward_clip(&r2, &i, false), Vec::<u8>::new());
    }

    #[test]
    fn colonies_rows_and_breakpoints() {
        let i = ins(InsType::FullInfo, Some(100), Some(200), vec![0]);
        let r = rec(
            Side::Left,
            b"",
            vec![row(3, Side::Left, Role::Clip, "a"), row(7, Side::Left, Role::Mate, "a"), row(1, Side::Left, Role::Short, "b")],
        );
        let mut bp: Breakpoints = FxHashMap::default();
        bp.insert((0, Side::Left), vec![(95, 5), (94, 6), (105, 3), (100, 2)]);
        bp.insert((0, Side::Right), vec![(100, 9)]);
        assert_eq!(colonies(&i, &r, Some(&bp), 5), vec![1, 2, 3, 5]);
        assert_eq!(colonies(&i, &r, None, 5), vec![1, 3]);
    }

    /// reference with a T tract at [200, 215) (RIGHT junction 200 points into it)
    fn genome() -> (crate::tprt::tests::VecFetch, Vec<u8>) {
        let mut g = Vec::new();
        let mut x: u32 = 12345;
        for _ in 0..400 {
            x = x.wrapping_mul(1103515245).wrapping_add(12345);
            g.push(b"ACG"[((x >> 16) % 3) as usize]);
        }
        for b in &mut g[200..215] {
            *b = b'T';
        }
        (crate::tprt::tests::VecFetch(g.clone()), g)
    }

    #[test]
    fn far_pair_and_slippage() {
        use crate::tprt::tests::{cfg, lh, StubMatcher};
        let c = Interner::new();
        c.intern("chr1");
        let cf = cfg();
        let (f, g) = genome();
        let m = StubMatcher(vec![("GGCCGGCCGGCCGGCCGGCC", lh("L1", b'+', 30))]);
        let i = ins(InsType::FullInfo, Some(100), Some(200), vec![0]);
        let el = b"AGGCCGGCCGGCCGGCCGGCCA".to_vec();
        // RIGHT clip: slipped T tract + the reference continuing after it
        let mut slipped = b"TTTTTTTTTTTT".to_vec();
        slipped.extend_from_slice(&g[215..245]);
        let mut real_tail = b"TTTTTTTTTTTT".to_vec();
        real_tail.extend_from_slice(b"GATCGATCGGATCCATGCAT");
        let recs = |l: &[u8], r: &[u8], fl: FileId, fr: FileId| {
            vec![
                rec(Side::Left, l, vec![row(fl, Side::Left, Role::Clip, "a")]),
                rec(Side::Right, r, vec![row(fr, Side::Right, Role::Clip, "b")]),
            ]
        };
        // passes the verdict; slippage at the poly-A side only with a reference
        let rs = recs(&el, &slipped, 1, 1);
        assert_eq!(far_pair_check(&i, &rs, &cf, &m, None, None, &c), None);
        assert_eq!(
            far_pair_check(&i, &rs, &cf, &m, None, Some(&f), &c),
            Some(("polya_side_slippage".to_string(), Some(Side::Right)))
        );
        // not a far pair
        let near = ins(InsType::FullInfo, Some(190), Some(200), vec![0]);
        assert_eq!(far_pair_check(&near, &rs, &cf, &m, None, Some(&f), &c), None);
        // colonies differ
        let rs2 = recs(&el, &real_tail, 1, 2);
        assert_eq!(
            far_pair_check(&i, &rs2, &cf, &m, None, Some(&f), &c),
            Some(("colony_mismatch".to_string(), Some(Side::Right)))
        );
        let rs3 = recs(&el, &real_tail, 1, 1);
        assert_eq!(far_pair_check(&i, &rs3, &cf, &m, None, Some(&f), &c), None);
        // slippage check: RIGHT slips; LEFT carries the element -> kept; otherwise rejected
        assert_eq!(slippage_check(&i, &rs, &cf, &m, &f, &c), "");
        let rs4 = recs(b"ACGTTGCAACGTAGCTAGCATT", &slipped, 1, 1);
        assert_eq!(slippage_check(&i, &rs4, &cf, &m, &f, &c), "slippage:RIGHT(repeat_shifted_reference)");
        // the mate-extended consensus (longer, not slippage) rescues
        let mut rs5 = recs(b"ACGTTGCAACGTAGCTAGCATT", &slipped, 1, 1);
        rs5[1].combined_consensus.seq = slipped.clone();
        rs5[1].consensus.seq = [real_tail.clone(), b"ACGATCGAGGATCCAGTAC".to_vec()].concat();
        assert_eq!(slippage_check(&i, &rs5, &cf, &m, &f, &c), "");
        assert!(carries_element(&rs[0], &i, &m, &c));
        assert!(!carries_element(&rs4[0], &i, &m, &c));
    }

    // inside_mates with MATE rows needs P3's `EvidenceRow::mapped` / `allele_forward_seq`
    // (todo!() in this worktree); it is exercised by tests/equiv.sh after integration.
}
