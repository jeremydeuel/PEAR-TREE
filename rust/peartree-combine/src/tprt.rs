//! Combine-level TPRT filter primitives. OWNER: P4.
//!
//! Mirrors src/combine_insertions_tprt_filters.py:45-311 (reference repeats, slippage rule,
//! far-pair verdict). The fuzzy-merge helpers (hp_compress / polya_trimmed / clips_agree) live
//! in `crate::seq` (P1). All sequences OUTWARD from the junction. SPEC.md §4.6.

use crate::align;
use crate::config::Config;
use crate::genome::RefFetch;
use crate::library::{LibHit, Matcher};
use crate::model::Side;
use crate::seq::revcomp;

/// `outward_reference(fetch, contig, junction, side, width=80)` -> (line, j).
/// lo = max(0, junction - width); s = fetch(contig, lo, junction + width) uppercased;
/// j = junction - lo; LEFT -> (revcomp(s), len(s) - j) else (s, j). NB python returns
/// ("", 0) only on an exception, which `get_sequence` never raises (it returns ""): an empty
/// fetch therefore yields ("", j) / ("", -j) -- callers only test `not line`.
pub fn outward_reference(fetch: &dyn RefFetch, contig: &str, junction: i64, side: Side, width: i64) -> (Vec<u8>, i64) {
    let lo = (junction - width).max(0);
    let s = fetch.fetch(contig, lo, junction + width).to_ascii_uppercase();
    let j = junction - lo;
    match side {
        Side::Left => {
            let n = s.len() as i64;
            (revcomp(&s), n - j)
        }
        Side::Right => (s, j),
    }
}

/// `repeat_at(line, j, max_period=6, slack=2)` -> (start, end, unit).
pub fn repeat_at(line: &[u8], j: i64, max_period: usize, slack: i64) -> (i64, i64, Vec<u8>) {
    let mut best: (i64, i64, Vec<u8>) = (j, j, Vec::new());
    let n = line.len() as i64;
    let at = |i: i64| line[i as usize];
    for p in 1..=max_period as i64 {
        let lo = (j - slack).max(0);
        let hi = (n - p).min(j + slack); // inclusive (python range(.., hi + 1))
        let mut st = lo;
        while st <= hi {
            let mut a = st;
            while a > 0 && at(a - 1) == at(a - 1 + p) {
                a -= 1;
            }
            let mut b = st;
            while b + p < n && at(b + p) == at(b) {
                b += 1;
            }
            b += p;
            if b - a >= 2 * p && b - a > best.1 - best.0 && (a <= j + slack && b >= j - slack) {
                best = (a, b, line[a as usize..(a + p) as usize].to_vec());
            }
            st += 1;
        }
    }
    best
}

/// `strip_repeat(s, unit, max_mismatch_frac=0.125)`.
pub fn strip_repeat(s: &[u8], unit: &[u8], max_mismatch_frac: f64) -> usize {
    if unit.is_empty() || s.is_empty() {
        return 0;
    }
    let p = unit.len();
    let mut best = 0;
    for ph in 0..p {
        let mut mm: i64 = 0;
        let mut last_ok = 0;
        let mut prev_bad = false;
        for (i, &c) in s.iter().enumerate() {
            if c == unit[(i + ph) % p] {
                last_ok = i + 1;
                prev_bad = false;
            } else {
                mm += 1;
                if prev_bad || mm > 1i64.max((max_mismatch_frac * (i + 1) as f64) as i64) {
                    break;
                }
                prev_bad = true;
            }
        }
        best = best.max(last_ok);
    }
    best
}

/// `structured_len(s, polya_min=8)`.
pub fn structured_len(s: &[u8], polya_min: usize) -> usize {
    let n = s.len();
    let mut i = 0;
    while i < n {
        let mut k = i;
        while k < n && s[k] == s[i] {
            k += 1;
        }
        if k - i >= polya_min {
            return i;
        }
        i = k;
    }
    n
}

/// python `s[j:]` for any integer j.
fn py_slice_from(s: &[u8], j: i64) -> &[u8] {
    let n = s.len() as i64;
    let k = if j < 0 { (n + j).max(0) } else { j.min(n) };
    &s[k as usize..]
}

/// `slippage_junction(clip, line, j, cfg)` -> "" | "repeat_only" | "repeat_junk" |
/// "repeat_shifted_reference".
pub fn slippage_junction(clip: &[u8], line: &[u8], j: i64, cfg: &Config) -> &'static str {
    if clip.is_empty() || line.is_empty() {
        return "";
    }
    let min_run = cfg.slippage_min_ref_run_combine as i64;
    let min_str = cfg.slippage_min_str_len as i64;
    let min_struct = cfg.slippage_min_structured;
    let polya_min = cfg.polya_min_len;
    let (a, b, unit) = repeat_at(line, j, cfg.slippage_max_period, 2);
    if unit.is_empty() {
        return "";
    }
    let need = if unit.len() == 1 { min_run } else { min_str.max(3 * unit.len() as i64) };
    if (b - a) < need {
        return "";
    }
    let clip = clip.to_ascii_uppercase();
    let k = strip_repeat(&clip, &unit, 0.125);
    let rest = &clip[k..];
    if k == 0 && !rest.is_empty() {
        // the clip opens with its OWN long homopolymer of another base (e.g. a poly-A tail on
        // the other strand): not a continuation of this tract
        let mut r = 0;
        while r < rest.len() && rest[r] == rest[0] {
            r += 1;
        }
        if r >= polya_min && !unit.contains(&rest[0]) {
            return "";
        }
    }
    let ref_out = py_slice_from(line, j);
    let ref_rest = &ref_out[strip_repeat(ref_out, &unit, 0.125)..];
    if structured_len(rest, polya_min) < min_struct {
        return "repeat_only";
    }
    if unit.len() == 1 {
        // post-homopolymer SBS phasing junk: the 'rest' is still dominated by the tract base
        let w = &rest[..rest.len().min(40)];
        let cnt = w.iter().filter(|&&c| c == unit[0]).count();
        if (cnt as f64) >= cfg.slippage_junk_frac * w.len() as f64 {
            return "repeat_junk";
        }
    }
    let m = rest.len().min(12);
    let budget = 1usize.max(m / 6);
    if ref_rest.len() >= m {
        let t = &ref_rest[..ref_rest.len().min(m + budget + 3)];
        if align::distance(&rest[..m], t, align::Mode::Shw, budget as i32, &[]) != -1 {
            return "repeat_shifted_reference";
        }
    }
    ""
}

/// `leading_polyt(seq, n=10)`.
pub fn leading_polyt(seq: &[u8], n: usize) -> bool {
    if seq.len() < n {
        return false;
    }
    let t = seq[..n].iter().filter(|c| c.eq_ignore_ascii_case(&b'T')).count();
    t * 5 >= n * 4
}

/// `after_polyt(seq)`.
pub fn after_polyt(seq: &[u8]) -> Vec<u8> {
    let s = seq.to_ascii_uppercase();
    let mut i: usize = 0;
    let mut bad = 0;
    while i < s.len() {
        if s[i] == b'T' {
            bad = 0;
        } else {
            bad += 1;
            if bad >= 2 {
                i -= 1;
                break;
            }
        }
        i += 1;
    }
    let r = &s[i..];
    let k = r.iter().position(|&c| c != b'T').unwrap_or(r.len());
    r[k..].to_vec()
}

/// `far_geometry(gap, cfg)`: gap < -far_pair_max_tsd_deletion or gap > far_pair_tsd_max.
pub fn far_geometry(gap: i64, cfg: &Config) -> bool {
    gap < -cfg.far_pair_max_tsd_deletion || gap > cfg.far_pair_tsd_max
}

/// `colonies_consistent(a, b, frac)`: share a colony and `|a ^ b| <= max(1, int(frac*|a|b|))`.
/// Colonies are sample-name ids (FileId); `a`, `b` sorted and unique (sets).
pub fn colonies_consistent(a: &[u32], b: &[u32], frac: f64) -> bool {
    let inter = a.iter().filter(|x| b.binary_search(x).is_ok()).count();
    if inter == 0 {
        return false;
    }
    let union = a.len() + b.len() - inter;
    let sym = union - inter;
    (sym as i64) <= 1i64.max((frac * union as f64) as i64)
}

/// Per-side inputs of `far_pair_verdict`. Index 0 = LEFT, 1 = RIGHT (python dict keys; a side
/// missing from the python dict is an empty Vec here, which python's `.get(s, ())` equals).
pub struct FarPairInput {
    /// candidate outward clips, best first, deduplicated, non-empty
    pub clips: [Vec<Vec<u8>>; 2],
    /// colonies (sorted, unique) per side
    pub colonies: [Vec<u32>; 2],
    /// `_inside_mates` per side
    pub inside_mates: [Vec<Vec<u8>>; 2],
}

fn sidx(s: Side) -> usize {
    match s {
        Side::Left => 0,
        Side::Right => 1,
    }
}

/// `far_pair_verdict(clips, colonies, matcher, cfg, inside_mates, n_ind=None)` -> (reason, polya
/// side). The `n_ind` / "few_fragments" branch belongs to the dropped pooled gate and is NOT
/// ported (python passes n_ind=None whenever require_independent_fragments is off).
pub fn far_pair_verdict(inp: &FarPairInput, matcher: &dyn Matcher, cfg: &Config) -> (&'static str, Option<Side>) {
    let n_pa = cfg.far_pair_min_polya;
    let pa: Vec<Side> = [Side::Left, Side::Right]
        .into_iter()
        .filter(|&s| inp.clips[sidx(s)].iter().any(|c| leading_polyt(c, n_pa)))
        .collect();
    if pa.len() != 1 {
        return ("no_polarity", None);
    }
    let pside = pa[0];
    let cside = pside.other();
    let anti = cfg.far_pair_allow_antisense;
    let mut h: Option<LibHit> = inp.clips[sidx(cside)].iter().find_map(|c| matcher.hit(c, false));
    if h.is_none() {
        let mh: Vec<LibHit> = inp.inside_mates[sidx(cside)]
            .iter()
            .filter_map(|m| matcher.hit(m, false))
            .filter(|x| x.class != "FLANK" && (anti || x.strand == b'+'))
            .collect();
        if mh.len() < 2 {
            return ("no_element_on_complex_clip", Some(pside));
        }
        // python max(key=mlen): the first maximal element
        let mut best = &mh[0];
        for x in &mh[1..] {
            if x.mlen > best.mlen {
                best = x;
            }
        }
        h = Some(best.clone());
    }
    let h = h.unwrap();
    if h.strand != b'+' && !anti {
        return ("element_antisense", Some(pside));
    }
    let tail = inp.clips[sidx(pside)].iter().find_map(|c| matcher.hit(&after_polyt(c), false));
    if let Some(t) = tail {
        if t.class != h.class {
            return ("element_class_conflict", Some(pside));
        }
    }
    if !colonies_consistent(&inp.colonies[0], &inp.colonies[1], cfg.far_pair_colony_frac) {
        return ("colony_mismatch", Some(pside));
    }
    ("", Some(pside))
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// A Config with the python defaults of every key P4 reads (P1's `from_json` not needed).
    pub(crate) fn cfg() -> Config {
        Config {
            genome_2bit: String::new(),
            exclude_files_with_many_insertions: 0,
            samtools_executable: String::new(),
            bowtie2_executable: String::new(),
            bowtie2_index: String::new(),
            bowtie2_index2: String::new(),
            bowtie2_index2_lo: String::new(),
            genotyping_max_bases: 30,
            clean_remap_max_insertion: 0,
            clean_remap_min_as: -15,
            trim_far_flank_before_remap: false,
            keep_polya_one_sided: false,
            merge_tolerance_bp: 0,
            polya_aware_clip_agreement: false,
            require_independent_fragments: false,
            min_independent_fragments: 2,
            indel_aware_consensus: false,
            dup_coord_tolerance: 5,
            dup_max_edit: 3,
            dup_max_edit_frac: 0.02,
            polya_min_len: 8,
            dup_mate_min_mapq: 20,
            count_short_overhang: false,
            short_overhang_min_bases: 5,
            short_overhang_min_ref_mismatch: 2,
            short_mate_max_dist: 1000,
            short_mate_min_mapq: 20,
            slippage_reject: false,
            far_pair_strict: false,
            far_pair_split: true,
            rte_library: "resources/rte_library".into(),
            slippage_min_ref_run_combine: 8,
            slippage_min_str_len: 12,
            slippage_min_structured: 10,
            slippage_max_period: 6,
            slippage_junk_frac: 0.5,
            far_pair_max_tsd_deletion: 30,
            far_pair_tsd_max: 40,
            far_pair_min_polya: 10,
            far_pair_allow_antisense: false,
            far_pair_colony_frac: 0.2,
            far_pair_colony_tol: 5,
            foldback_filter: false,
            foldback_k: 20,
            foldback_min_short: 12,
            foldback_window: 50,
            foldback_max_mismatch: 1,
            foldback_min_entropy: 1.0,
            config_dir: None,
        }
    }

    /// Stub matcher: a hit when the uppercased query (>= 20 bp) contains a key (first key wins).
    pub(crate) struct StubMatcher(pub Vec<(&'static str, LibHit)>);
    impl Matcher for StubMatcher {
        fn hit(&self, seq: &[u8], _wf: bool) -> Option<LibHit> {
            let s = String::from_utf8_lossy(&seq.to_ascii_uppercase()).into_owned();
            if s.len() < 20 {
                return None;
            }
            self.0.iter().find(|(k, _)| s.contains(k)).map(|(_, h)| h.clone())
        }
    }
    pub(crate) fn lh(class: &'static str, strand: u8, mlen: i64) -> LibHit {
        LibHit { ctg: format!("{class}_X"), class, strand, mlen, q_st: 0, q_en: mlen }
    }

    /// A RefFetch over one in-memory sequence with py2bit clamping semantics.
    pub(crate) struct VecFetch(pub Vec<u8>);
    impl RefFetch for VecFetch {
        fn fetch(&self, _n: &str, s: i64, e: i64) -> Vec<u8> {
            let e = e.min(self.0.len() as i64);
            if s >= e || s < 0 {
                return vec![];
            }
            self.0[s as usize..e as usize].to_ascii_uppercase()
        }
    }

    #[test]
    fn small_primitives() {
        assert_eq!(repeat_at(b"ACGTAAAAAAAAAGCT", 6, 6, 2), (4, 13, b"A".to_vec()));
        assert_eq!(repeat_at(b"ACGTACGT", 4, 6, 2), (0, 8, b"ACGT".to_vec()));
        assert_eq!(repeat_at(b"ACG", 1, 6, 2), (1, 1, vec![]));
        assert_eq!(strip_repeat(b"AAAAAGAAAC", b"A", 0.125), 9);
        assert_eq!(strip_repeat(b"CACACACATT", b"AC", 0.125), 8);
        assert_eq!(structured_len(b"ACGTAAAAAAAAG", 8), 4);
        assert_eq!(structured_len(b"ACGT", 8), 4);
        assert!(leading_polyt(b"TTTTTTTTAA", 10));
        assert!(!leading_polyt(b"TTTTTTTAAA", 10));
        assert!(!leading_polyt(b"TTTTTTTT", 10));
        assert_eq!(after_polyt(b"ttttttttgcaggt"), b"GCAGGT".to_vec());
        assert_eq!(after_polyt(b"TTTTATTTTGCA"), b"GCA".to_vec());
        assert_eq!(after_polyt(b"GC"), b"GC".to_vec());
        assert_eq!(after_polyt(b"TTTT"), Vec::<u8>::new());
        assert!(colonies_consistent(&[1, 2, 3], &[1, 2], 0.2));
        assert!(!colonies_consistent(&[1, 2, 3], &[1], 0.2));
        assert!(!colonies_consistent(&[1], &[2], 0.2));
        let c = cfg();
        assert!(far_geometry(-31, &c) && far_geometry(41, &c) && !far_geometry(-30, &c) && !far_geometry(40, &c));
    }

    #[test]
    fn outward_reference_orientation() {
        let f = VecFetch(b"acgtTTTTGGGG".to_vec());
        assert_eq!(outward_reference(&f, "c", 4, Side::Right, 3), (b"CGTTTT".to_vec(), 3));
        assert_eq!(outward_reference(&f, "c", 4, Side::Left, 3), (b"AAAACG".to_vec(), 3));
        assert_eq!(outward_reference(&f, "c", 1, Side::Right, 3), (b"ACGT".to_vec(), 1));
        assert_eq!(outward_reference(&f, "c", 1, Side::Left, 3), (b"ACGT".to_vec(), 3));
        assert_eq!(outward_reference(&f, "c", 100, Side::Left, 3), (vec![], -3));
    }

    #[test]
    fn verdict_paths() {
        let c = cfg();
        let m = StubMatcher(vec![
            ("GGCCGGCCGGCCGGCCGGCC", lh("L1", b'+', 30)),
            ("CCGGAACCGGAACCGGAACC", lh("L1", b'-', 30)),
            ("ATATCGCGATATCGCGATAT", lh("ALU", b'+', 25)),
        ]);
        let pt = b"TTTTTTTTTTTTACGTAGCTAGCTAGCATCGATCGA".to_vec();
        let el = b"AGGCCGGCCGGCCGGCCGGCCA".to_vec();
        let mk = |l: Vec<Vec<u8>>, r: Vec<Vec<u8>>, cl: Vec<u32>, cr: Vec<u32>, ml: Vec<Vec<u8>>| FarPairInput {
            clips: [l, r],
            colonies: [cl, cr],
            inside_mates: [ml, vec![]],
        };
        assert_eq!(far_pair_verdict(&mk(vec![], vec![], vec![], vec![], vec![]), &m, &c), ("no_polarity", None));
        assert_eq!(
            far_pair_verdict(&mk(vec![pt.clone()], vec![pt.clone()], vec![], vec![], vec![]), &m, &c),
            ("no_polarity", None)
        );
        assert_eq!(
            far_pair_verdict(&mk(vec![el.clone()], vec![pt.clone()], vec![1], vec![1], vec![]), &m, &c),
            ("", Some(Side::Right))
        );
        assert_eq!(
            far_pair_verdict(&mk(vec![el.clone()], vec![pt.clone()], vec![1], vec![2], vec![]), &m, &c),
            ("colony_mismatch", Some(Side::Right))
        );
        let anti = b"ACCGGAACCGGAACCGGAACCA".to_vec();
        assert_eq!(
            far_pair_verdict(&mk(vec![anti], vec![pt.clone()], vec![1], vec![1], vec![]), &m, &c),
            ("element_antisense", Some(Side::Right))
        );
        let mut ptc = b"TTTTTTTTTTTTT".to_vec();
        ptc.extend_from_slice(b"GGATATCGCGATATCGCGATAT");
        assert_eq!(
            far_pair_verdict(&mk(vec![el.clone()], vec![ptc], vec![1], vec![1], vec![]), &m, &c),
            ("element_class_conflict", Some(Side::Right))
        );
        // no clip hit: >= 2 sense non-FLANK inside-mate hits rescue
        let short = b"ACGTAC".to_vec();
        assert_eq!(
            far_pair_verdict(&mk(vec![short.clone()], vec![pt.clone()], vec![1], vec![1], vec![el.clone()]), &m, &c),
            ("no_element_on_complex_clip", Some(Side::Right))
        );
        assert_eq!(
            far_pair_verdict(
                &mk(vec![short], vec![pt.clone()], vec![1], vec![1], vec![el.clone(), el.clone()]),
                &m,
                &c
            ),
            ("", Some(Side::Right))
        );
    }

    /// python reference values (stratified sample of the randomized differential cases)
    const CASES: &[&str] = &[
        "slip\tAAAAAAAAAAAAAAAAAAACCTGCTCGTTATGCCAGAC\tCCATGTTCGGCTCTAGAAATTGTGGATCTAAACCAGGGTG\t20\t\t20,24,TG",
        "slip\tGGAACTTCGATGAAGTGAT\tGATCGTTGTACTACGTACGATCAAGATCAAGATCAAGATCAA\t22\t\t18,42,GATCAA",
        "slip\tAAGTAGACGTCGAAGTTATACCGCAGAAGTGTC\tTCCACTAGCTCAAGTTTGGAAGTCGAAGTCGAAGTCGAAGTCGAAG\t20\t\t18,46,GAAGTC",
        "slip\tGTCACGCCATGCCATGCCATGCCATGCCATCACG\tTGTCCGCCGGCGGGAGGAAGTGTGCGAACTGTATGAACGT\t0\t\t0,0,",
        "slip\tACGACGACGACGACGACGACGCTAGCCCGTTAGGTG\tAGTAAATTCTGCCAGCGATCCCACGACGACGACGACGACGACGACGACGAC\t-3\t\t-3,-3,",
        "slip\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAACACTC\tAGCATCCAGTCAACGCGAGGGATTTTTTTTTTTTTTTTTTTTTTAT\t-3\t\t-3,-3,",
        "slip\tTCTGTCTGACTATCTGTCTGTCTGTCGTGGTATTTATCCGC\tAAAACAATTCCACCGGTTCTATCTGTCTGTCTGTCTGTCTGTCTGTC\t18\t\t16,18,T",
        "slip\tTCACACTCTTTTTTTTTTTTGCGTTTTTTTTTTTTGCG\tGCACGCGTCATAAGCAGGAGTCACACTCACACTCACTCCC\t23\t\t20,36,TCACAC",
        "slip\tAA\tTTGAGCTTTGAGGTCAGATTCAGAACTCCAACCGATGCTATTATTATTATTATTATTATTATCATCGCAGAAATAAAGCG\t38\trepeat_only\t38,62,TAT",
        "slip\tCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCC\tATCGTCCGGCTAACGTTCCCCCCCCCCCCCCTCCCCCCCCCC\t19\trepeat_only\t17,31,C",
        "slip\tAGGGCG\tCCTCTGAATGTCACGGGACCGAAAAAAAAAAAAACCTCTG\t20\trepeat_only\t21,34,A",
        "slip\tCCCCCCCCCCCCCCCCCCAACCCCCCCCCCCCCAA\tGTTGTTACGACCGACATCCCCCCCCCCCCCCCCCCCCCCCCC\t18\trepeat_only\t17,42,C",
        "slip\tGTGTGTGTGCGTACACGT\tAGATCCAACATAACCCGGCAGGTGTGTGTGTGTGTGTGTGTGTGTG\t46\trepeat_only\t21,46,GT",
        "slip\tggag\tTCTATATTGGCGATTCCGCCCGACCTAGAGCGGGATCGTCCCCCCCCCCCCCCCCCCCCCCCCCCCGCGCCGTCGTACTT\t41\trepeat_only\t39,66,C",
        "slip\tt\tTGGGTAGTACGGCAATGGGGCTGTTGTTGTTGTTGTTGTTGTTGTT\t19\trepeat_only\t21,46,TGT",
        "slip\tAGTTGCACCC\tAGCTTCAATCATATGAGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT\t20\trepeat_only\t17,47,GT",
        "slip\tTTTTTTTTTTTTTTTTTTTTTTTTCGTTCGACCCA\tTAATTCGCGTGCCAAGTTGTATTTTTTTTCGTTCGACCCA\t22\trepeat_shifted_reference\t21,29,T",
        "slip\tAAAAAAAAAAAAAAAAAAAAAAAACCTGTCACTACC\tAGTAGGTTGAGCGGGTGCCTCCCTGTCAGAATGTACAAAAAAAAAAAAAAAAAAAAAAAAACCTGTCACTACCTGATATG\t40\trepeat_shifted_reference\t36,61,A",
        "slip\tAAAAAAAAAAAAAAAAAATTCACCCTTTTTTCAGAGCTA\tACCACCCTCCTTGATAGTGTAAGAAACAGAGGGAGAAAAAAAAAAAAAATTCACCCTTTTTTCAGAGCTATACACAGAGT\t40\trepeat_shifted_reference\t35,49,A",
        "slip\tCCTGTAGTTCGGTTCCTCCTGAACAAGTGTC\tAAGGTCTTACACCTGAAAACTAAAATGCGCTAATACCCCCCCCCGCCTGTAGTTCGGTTCCTCCTGAACAAGTGTCGACT\t40\trepeat_shifted_reference\t35,44,C",
        "slip\ttgaatatgaatatgaattaatcgtggtgtgt\tCCGACCTGTATAATGGAGGAAAATTAGAGTATAGGGCCTTATGAATATGAATATGAATATGAATTAATCGTGGTGTGTGG\t40\trepeat_shifted_reference\t39,64,TATGAA",
        "slip\tGGTCAAAAGGTCAAGGTCAAGGGCATAAGGATGTAA\tGCATAGAACCTGAACACATATGTCATCACAAAACACATAAGGGAGGTCAAGGTCAAGGTCAAGGGCATAAGGATGTAACG\t41\trepeat_shifted_reference\t43,64,AGGTCA",
        "slip\tAAAAAAAAAAAAAAAAAAAAAAAATTGTAGATGAC\tGTTCTATAGTAGGCTTGGCAAAAAAAAAATTGTAGATGAC\t22\trepeat_shifted_reference\t19,29,A",
        "slip\tAAAAAAAAAAAAAAGAAGTCCTGTGGAGTCACTACCCCTGG\tGTACCCCCTCTATTTTTAGAACCGAAGCTTTAGCATGTCAAAAAAAAGAAGTCCTGTGGAGTCACTACCCCTGGTAACGA\t40\trepeat_shifted_reference\t39,47,A",
        "slip\tTATTTCTTTTTTCATT\tTTGCTTCCCGGACTCCTATCTTTTTTTTTTTTTTTTTTTTTTTT\t22\trepeat_junk\t20,44,T",
        "slip\tGGGAGAGGGGGATAGC\tACGACTGGTCAGCAGGCGTTCCGGGGGGGGGGGGGGGGTT\t20\trepeat_junk\t22,38,G",
        "slip\tggggaggggggggggcaggggtggcgt\tCCTGTTCTTCGCAACTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\t46\trepeat_junk\t16,46,G",
        "slip\tCCGCCCGCGCCCCCCTCCCGCCTACTCCAGGCGTTC\tGACATCTCTGGCGTCTTCTACCCCCCCCCCCCCCCCCCCCC\t41\trepeat_junk\t20,41,C",
        "slip\tTTTTTTGTGTTTTTTTCGCGTAGTC\tGCGCGCCGAACAGGTCAGTTTTTTTTGGACCTACAGGTTA\t20\trepeat_junk\t18,26,T",
        "slip\tGGGGTGGGGGCGTGAGGGGGC\tTTTACGTGAATGGGTTAATAGCTTGGGGGGGGGTTGGTTC\t23\trepeat_junk\t24,33,G",
        "slip\tACCCTTTTCGCCTCCCGC\tGGTGACAATACCCCAGCCCCCCCCCCCCCCCCCCCCCCCCC\t41\trepeat_junk\t16,41,C",
        "slip\tgggggaggggggtggagggagggcgatactgag\tTGGCTCAGGGCGCGCTGGGGGGGGGGGGGGAGGAGCATAA\t20\trepeat_junk\t16,30,G",
        "strip\tAAAAAAAAAAAAAAAAAAAAGAC\tA\t22",
        "strip\tTGTGTC\tCACT\t1",
        "strip\tGGTCGATCGATCGATCGATCGATCATATC\tTCGA\t24",
        "strip\tTGAGTGAGTGAGTGAGTGAGTGTA\tTGAG\t22",
        "strip\tGGGGGGGGGGGGGGGGGGGGGAT\tG\t21",
        "strip\tCCCTCACCT\tTCACCT\t9",
        "strip\tTTTTTAG\tT\t5",
        "strip\tGTCGACGACGG\tCGA\t10",
        "strip\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTT\tCCC\t30",
        "strip\tCGAAGCGATGCGATGCGATGCGATTCGC\tATGCG\t27",
        "apt\ttcgtgtatgtttttatttttt\tCGTGTATGTTTTTATTTTTT\t0",
        "apt\tTTTTTTTTGCTTTTGACTCTCATAATCT\tGCTTTTGACTCTCATAATCT\t1",
        "apt\tttttttttttttttctgtgttggggtaaa\tGGGGTAAA\t1",
        "apt\tTTTTTTTTTTTTTTTTTTTTGACAT\tGACAT\t1",
        "apt\tTTTTTTTTTTTTTTTTTTTTTTTTTTGTTTA\t\t1",
        "apt\tCTTTTTTTCCACCGC\tCCACCGC\t0",
        "apt\tATTTTTTTTTTTCGGCAT\tCGGCAT\t1",
        "apt\tTATTGTTTTTTTTTGTTTTTTTCCTT\tCCTT\t1",
        "apt\tTTATATTTTATTTCTCTATCTCGAAGATTCTCTAT\tCGAAGATTCTCTAT\t0",
        "apt\tTTCATCACTCACCAGG\tCATCACTCACCAGG\t0",
        "col\t0,2,7,8,9,10,11\t1,6\t0",
        "col\t7\t\t0",
        "col\t0,1,2,3,6,9,10,11\t6,7\t0",
        "col\t0,1,2,3,5,6,8,10,11\t0,1,3,4,5,8,10,11\t0",
        "col\t0,2,4,6,7,8,9\t0,1,3,6,7,8,9,10\t0",
        "col\t\t1,2,5,7,10\t0",
        "col\t9\t0,1,5,8,9\t0",
        "col\t0,1,2,5,6,7,8,9,10,11\t0,1,2,3,5,6,7,9,10,11\t1",
        "col\t0,2,4,7,9,10\t2,6,7,8,9\t0",
        "col\t2,3,4,5,6,7,8,9,10\t1,2,3,4,6,9,10,11\t0",
    ];

    #[test]
    fn python_differential_literals() {
        assert_eq!(check_cases(CASES.iter().copied()), CASES.len());
    }

    /// Differential check against python on a large randomized case file (generator script in
    /// the P4 scratch dir); run with `P4_TPRT_CASES=<file> cargo test ... -- --ignored`.
    #[test]
    #[ignore]
    fn python_differential_file() {
        let p = std::env::var("P4_TPRT_CASES").unwrap();
        let s = std::fs::read_to_string(p).unwrap();
        let n = check_cases(s.lines());
        eprintln!("{n} cases identical");
    }

    fn check_cases<'a>(lines: impl Iterator<Item = &'a str>) -> usize {
        let c = cfg();
        let mut n = 0;
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            match f[0] {
                "slip" => {
                    let j: i64 = f[3].parse().unwrap();
                    let (a, b, u) = repeat_at(f[2].as_bytes(), j, 6, 2);
                    assert_eq!(format!("{a},{b},{}", String::from_utf8(u).unwrap()), f[5], "repeat_at {line}");
                    assert_eq!(slippage_junction(f[1].as_bytes(), f[2].as_bytes(), j, &c), f[4], "slip {line}");
                }
                "strip" => {
                    assert_eq!(strip_repeat(f[1].as_bytes(), f[2].as_bytes(), 0.125).to_string(), f[3], "{line}");
                }
                "apt" => {
                    assert_eq!(String::from_utf8(after_polyt(f[1].as_bytes())).unwrap(), f[2], "{line}");
                    assert_eq!(if leading_polyt(f[1].as_bytes(), 10) { "1" } else { "0" }, f[3], "{line}");
                }
                "col" => {
                    let parse = |s: &str| -> Vec<u32> {
                        if s.is_empty() {
                            vec![]
                        } else {
                            s.split(',').map(|x| x.parse().unwrap()).collect()
                        }
                    };
                    let r = colonies_consistent(&parse(f[1]), &parse(f[2]), 0.2);
                    assert_eq!(if r { "1" } else { "0" }, f[3], "{line}");
                }
                _ => panic!("bad case {line}"),
            }
            n += 1;
        }
        n
    }
}
