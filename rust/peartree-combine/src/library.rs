//! RTE consensus-library matcher (python `LibraryMatcher`, mappy). OWNER: P4.
//!
//! Mirrors src/combine_insertions_tprt_filters.py:171-229. SPEC.md §4.6 "LibraryMatcher".
//! Index options exactly as python: `k=11, w=3, min_cnt=1, min_chain_score=15,
//! min_dp_score=20 (map_opt.min_dp_max), best_n=3`, no preset, alignment on (mappy always sets
//! MM_F_CIGAR). Element index = `<lib>/consensus.fa`; flank indices = `flanks_3p.fa.gz`,
//! `flanks_5p_sva.fa.gz` when present (skipped when missing or empty).
//!
//! Implementation: the minimap2 C API through the `minimap2` crate's FFI (`minimap2::ffi`),
//! replicating mappy's `Aligner.__cinit__` / `Aligner.map` call sequence exactly
//! (`mm_set_opt(NULL)`, flag |= MM_F_CIGAR, uni-part index, `mm_idx_reader_read` of the first
//! part, `mm_mapopt_update`, `mm_idx_index_name`; per query `mm_map(idx, strlen(seq), seq, ..,
//! name=NULL)`). One `mm_tbuf_t` per thread (thread-local), so `hit` can be called from any
//! number of rayon workers concurrently.

use std::path::Path;

/// `element_class(name)`: "L1" / "ALU" / "SVA" by case-insensitive prefix, else "FLANK".
pub fn element_class(name: &str) -> &'static str {
    let u = name.to_ascii_uppercase();
    for k in ["L1", "ALU", "SVA"] {
        if u.starts_with(k) {
            return k;
        }
    }
    "FLANK"
}

/// python hit tuple `(name, class, strand, mlen, q_st, q_en)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LibHit {
    pub ctg: String,
    /// element_class(ctg) for the element index, "FLANK" for flank indices
    pub class: &'static str,
    /// b'+' = the query as given is sense to the consensus
    pub strand: u8,
    pub mlen: i64,
    pub q_st: i64,
    pub q_en: i64,
}

/// Library matcher interface (a trait so unit tests can stub it). `Sync`: shared by rayon
/// workers.
pub trait Matcher: Sync {
    /// `LibraryMatcher.hit(seq, with_flanks)`: uppercase; len < 20 -> None; over the element
    /// index then (with_flanks) each flank index, in that order, over each index's hits in
    /// minimap2 order: keep the first hit with the strictly largest `mlen >= min_match (18)`.
    /// (The python result cache does not change results; cache or not.)
    fn hit(&self, seq: &[u8], with_flanks: bool) -> Option<LibHit>;
}

/// The minimap2-backed matcher (feature `mappy`).
pub struct LibraryMatcher {
    #[cfg(feature = "mappy")]
    elem: mm::Index,
    #[cfg(feature = "mappy")]
    flank: Vec<mm::Index>,
    #[allow(dead_code)]
    min_match: i64,
}

impl LibraryMatcher {
    /// `LibraryMatcher(library_dir, flanks=True, min_match=18)`. Error when consensus.fa cannot be
    /// indexed (python RuntimeError).
    #[cfg(feature = "mappy")]
    pub fn open(library_dir: &Path) -> Result<LibraryMatcher, String> {
        let ep = library_dir.join("consensus.fa");
        let elem = mm::Index::open(&ep).ok_or_else(|| format!("cannot index {}", ep.display()))?;
        let mut flank = Vec::new();
        for f in ["flanks_3p.fa.gz", "flanks_5p_sva.fa.gz"] {
            let p = library_dir.join(f);
            if p.exists() {
                if let Some(a) = mm::Index::open(&p) {
                    flank.push(a);
                }
            }
        }
        Ok(LibraryMatcher { elem, flank, min_match: 18 })
    }

    #[cfg(not(feature = "mappy"))]
    pub fn open(library_dir: &Path) -> Result<LibraryMatcher, String> {
        Err(format!(
            "cannot index {}: peartree-combine was built without the `mappy` feature",
            library_dir.join("consensus.fa").display()
        ))
    }
}

impl Matcher for LibraryMatcher {
    #[cfg(feature = "mappy")]
    fn hit(&self, seq: &[u8], with_flanks: bool) -> Option<LibHit> {
        let s = seq.to_ascii_uppercase();
        if s.len() < 20 {
            return None;
        }
        // mappy passes `strlen(seq)`: a NUL byte would truncate (never occurs in DNA input)
        let s: &[u8] = match s.iter().position(|&c| c == 0) {
            Some(i) => &s[..i],
            None => &s,
        };
        let mut best: Option<LibHit> = None;
        let nflank = if with_flanks { self.flank.len() } else { 0 };
        for ai in 0..=nflank {
            let al = if ai == 0 { &self.elem } else { &self.flank[ai - 1] };
            al.map(s, |h| {
                if h.mlen >= self.min_match && best.as_ref().is_none_or(|b| h.mlen > b.mlen) {
                    let ctg = al.name(h.rid);
                    best = Some(LibHit {
                        class: if ai == 0 { element_class(ctg) } else { "FLANK" },
                        ctg: ctg.to_string(),
                        strand: if h.rev { b'-' } else { b'+' },
                        mlen: h.mlen,
                        q_st: h.q_st,
                        q_en: h.q_en,
                    });
                }
            });
        }
        best
    }

    #[cfg(not(feature = "mappy"))]
    fn hit(&self, _seq: &[u8], _with_flanks: bool) -> Option<LibHit> {
        unreachable!("LibraryMatcher cannot be constructed without the `mappy` feature")
    }
}

/// Thin safe wrapper over the minimap2 C API (exactly the mappy call sequence).
#[cfg(feature = "mappy")]
pub(crate) mod mm {
    use minimap2::ffi;
    use std::cell::RefCell;
    use std::ffi::{CStr, CString};
    use std::os::raw::c_char;
    use std::path::Path;

    /// One raw minimap2 hit.
    #[derive(Clone, Copy, Debug)]
    pub struct RawHit {
        pub rid: i32,
        pub rev: bool,
        pub mlen: i64,
        pub q_st: i64,
        pub q_en: i64,
        pub r_st: i64,
        pub r_en: i64,
        pub blen: i64,
    }

    /// A built uni-part minimap2 index + the mapping options mappy would use with it.
    pub struct Index {
        idx: *mut ffi::mm_idx_t,
        mapopt: ffi::mm_mapopt_t,
        names: Vec<String>,
    }

    // SAFETY: after construction the index and options are only read; mm_map is thread-safe
    // given a distinct mm_tbuf_t per concurrent call (minimap2 API contract), which `with_tbuf`
    // guarantees (thread-local).
    unsafe impl Send for Index {}
    unsafe impl Sync for Index {}

    struct TBuf(*mut ffi::mm_tbuf_t);
    impl Drop for TBuf {
        fn drop(&mut self) {
            // SAFETY: created by mm_tbuf_init, destroyed once.
            unsafe { ffi::mm_tbuf_destroy(self.0) }
        }
    }
    thread_local! {
        static TBUF: RefCell<Option<TBuf>> = const { RefCell::new(None) };
    }

    fn with_tbuf<R>(f: impl FnOnce(*mut ffi::mm_tbuf_t) -> R) -> R {
        TBUF.with(|c| {
            let mut c = c.borrow_mut();
            // SAFETY: plain allocation
            let b = c.get_or_insert_with(|| TBuf(unsafe { ffi::mm_tbuf_init() }));
            f(b.0)
        })
    }

    impl Index {
        /// mappy `Aligner(path, k=11, w=3, min_cnt=1, min_chain_score=15, min_dp_score=20,
        /// best_n=3)`; None when python's `not aligner` (no index could be read).
        pub fn open(path: &Path) -> Option<Index> {
            use std::os::unix::ffi::OsStrExt;
            let cpath = CString::new(path.as_os_str().as_bytes()).ok()?;
            // SAFETY: zeroed POD option structs are immediately initialised by mm_set_opt.
            let mut io: ffi::mm_idxopt_t = unsafe { std::mem::zeroed() };
            let mut mo: ffi::mm_mapopt_t = unsafe { std::mem::zeroed() };
            unsafe { ffi::mm_set_opt(std::ptr::null(), &mut io, &mut mo) };
            mo.flag |= ffi::MM_F_CIGAR as i64; // always perform alignment
            io.batch_size = 0x7fff_ffff_ffff_ffff; // always build a uni-part index
            io.k = 11;
            io.w = 3;
            mo.min_cnt = 1;
            mo.min_chain_score = 15;
            mo.min_dp_max = 20;
            mo.best_n = 3;
            // SAFETY: valid C string and option structs; reader closed below.
            let r = unsafe { ffi::mm_idx_reader_open(cpath.as_ptr(), &io, std::ptr::null()) };
            if r.is_null() {
                return None;
            }
            let idx = unsafe { ffi::mm_idx_reader_read(r, 3) }; // mappy n_threads=3; first part only
            unsafe { ffi::mm_idx_reader_close(r) };
            if idx.is_null() {
                return None;
            }
            unsafe {
                ffi::mm_mapopt_update(&mut mo, idx);
                ffi::mm_idx_index_name(idx);
            }
            let n = unsafe { (*idx).n_seq } as usize;
            let mut names = Vec::with_capacity(n);
            for i in 0..n {
                // SAFETY: i < n_seq; names are NUL-terminated C strings owned by the index.
                let p = unsafe { (*(*idx).seq.add(i)).name };
                names.push(unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned());
            }
            Some(Index { idx, mapopt: mo, names })
        }

        pub fn name(&self, rid: i32) -> &str {
            &self.names[rid as usize]
        }

        /// `for h in aligner.map(seq)`, calling `f` per hit in minimap2 order.
        pub fn map(&self, seq: &[u8], mut f: impl FnMut(&RawHit)) {
            if seq.is_empty() {
                return;
            }
            let mut n_regs: i32 = 0;
            let regs = with_tbuf(|b| unsafe {
                // SAFETY: seq valid for its length; b is this thread's buffer; idx/mapopt read-only
                ffi::mm_map(
                    self.idx,
                    seq.len() as i32,
                    seq.as_ptr() as *const c_char,
                    &mut n_regs,
                    b,
                    &self.mapopt,
                    std::ptr::null(),
                )
            });
            if regs.is_null() {
                return;
            }
            let n = n_regs.max(0) as usize;
            let mut hits = Vec::with_capacity(n);
            for i in 0..n {
                // SAFETY: regs has n_regs elements; each `p` was malloc'd by minimap2.
                let r = unsafe { &*regs.add(i) };
                hits.push(RawHit {
                    rid: r.rid,
                    rev: r.rev() != 0,
                    mlen: r.mlen as i64,
                    q_st: r.qs as i64,
                    q_en: r.qe as i64,
                    r_st: r.rs as i64,
                    r_en: r.re as i64,
                    blen: r.blen as i64,
                });
                unsafe { ffi::free(r.p as *mut _) };
            }
            unsafe { ffi::free(regs as *mut _) };
            for h in &hits {
                f(h);
            }
        }
    }

    impl Drop for Index {
        fn drop(&mut self) {
            // SAFETY: built by mm_idx_reader_read, destroyed once.
            unsafe { ffi::mm_idx_destroy(self.idx) }
        }
    }
}

#[cfg(all(test, feature = "mappy"))]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn lib_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/rte_library")
    }

    #[test]
    fn element_class_prefixes() {
        assert_eq!(element_class("L1HS"), "L1");
        assert_eq!(element_class("AluYa5"), "ALU");
        assert_eq!(element_class("sva_e"), "SVA");
        assert_eq!(element_class("chr1:123"), "FLANK");
    }

    /// Dump every raw hit of every query in `$P4_MAPPY_QUERIES` (one per line) against the three
    /// indices, in the format of the mappy parity script (`n idx ctg strand mlen q_st q_en r_st
    /// r_en blen`), to `$P4_MAPPY_OUT`. Run manually:
    /// `P4_MAPPY_QUERIES=q.txt P4_MAPPY_OUT=rust.tsv cargo test --release dump_raw_hits -- --ignored`
    #[test]
    #[ignore]
    fn dump_raw_hits() {
        use std::io::Write;
        let q = std::env::var("P4_MAPPY_QUERIES").unwrap();
        let o = std::env::var("P4_MAPPY_OUT").unwrap();
        let d = lib_dir();
        let idx: Vec<mm::Index> = ["consensus.fa", "flanks_3p.fa.gz", "flanks_5p_sva.fa.gz"]
            .iter()
            .map(|f| mm::Index::open(&d.join(f)).unwrap())
            .collect();
        let mut out = std::io::BufWriter::new(std::fs::File::create(o).unwrap());
        for (n, line) in std::fs::read_to_string(q).unwrap().lines().enumerate() {
            for (ai, al) in idx.iter().enumerate() {
                al.map(line.trim().as_bytes(), |h| {
                    writeln!(
                        out,
                        "{n}\t{ai}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        al.name(h.rid),
                        if h.rev { -1 } else { 1 },
                        h.mlen,
                        h.q_st,
                        h.q_en,
                        h.r_st,
                        h.r_en,
                        h.blen
                    )
                    .unwrap();
                });
            }
        }
    }

    /// `LibraryMatcher.hit` results generated with mappy 2.31 (reference venv) for fixture clips.
    #[test]
    fn hits_match_mappy() {
        let d = lib_dir();
        let m = LibraryMatcher::open(&d).unwrap();
        for (q, wf, exp) in EXPECTED {
            let got = m.hit(q.as_bytes(), *wf);
            let exp = exp.map(|(c, cl, st, ml, a, b)| LibHit {
                ctg: c.to_string(),
                class: cl,
                strand: st,
                mlen: ml,
                q_st: a,
                q_en: b,
            });
            assert_eq!(got, exp, "query {q} with_flanks={wf}");
        }
        // < 20 bp never hits; lowercase is uppercased
        assert_eq!(m.hit(b"ACGTACGTACGTACGTACG", true), None);
        let (q, wf, _) = EXPECTED.iter().find(|e| e.2.is_some()).unwrap();
        assert_eq!(m.hit(q.to_ascii_lowercase().as_bytes(), *wf), m.hit(q.as_bytes(), *wf));
    }

    #[test]
    fn concurrent_use() {
        let m = LibraryMatcher::open(&lib_dir()).unwrap();
        let serial: Vec<_> = EXPECTED.iter().map(|(q, wf, _)| m.hit(q.as_bytes(), *wf)).collect();
        use rayon::prelude::*;
        for _ in 0..4 {
            let par: Vec<_> = EXPECTED.par_iter().map(|(q, wf, _)| m.hit(q.as_bytes(), *wf)).collect();
            assert_eq!(par, serial);
        }
    }

    type Exp = (&'static str, bool, Option<(&'static str, &'static str, u8, i64, i64, i64)>);
    const EXPECTED: &[Exp] = &[
        ("GGAGGCTGAGGCAGGAGAAT", false, Some(("ALU_YA5", "ALU", b'+', 20, 0, 20))),
        ("GTAGAGATGGGGTCTCTCTA", true, None),
        ("CTGCCTGATTCTCCTGCCTCA", false, Some(("SVA_E", "SVA", b'+', 21, 0, 21))),
        ("AGGTCATTCTGGTCTGCTGAG", true, None),
        ("CAGGCGCGTGCCACCACGCCCGGCT", false, None),
        ("TTTTTTTTTCTTTTTTTTTTTTTTT", true, Some(("SVA_chr10_12530761_r", "FLANK", b'-', 25, 0, 25))),
        ("GGGAGGCTGAGGCAGGAGAATGGCGTG", false, Some(("ALU_YA5", "ALU", b'+', 27, 0, 27))),
        ("GCGTGGTGGCACGTGCCTGTAGTCCCA", true, Some(("L1_chrX_30805015_r", "FLANK", b'-', 27, 0, 27))),
        ("AATAGCTTAGAGAAGTCTTTGGAAGCC", false, None),
        ("AGTGTGAGCCACCACACCCGGGGGCACCC", true, Some(("SVA_chr19_53478125_r", "FLANK", b'-', 20, 1, 21))),
        ("CAGTAGCACAATCTGGGCTCACTGCAACC", false, None),
        ("TTTTTGTCCTTTTTTTTTTTTTTTTTTTTT", true, Some(("L1_chr2_148639164_f", "FLANK", b'-', 24, 0, 25))),
        ("AGTGTCTGTTCAGGTCCTTCGCCCACTTTT", false, Some(("L1PA3", "L1", b'-', 29, 0, 30))),
        ("CTACTTTCTAATTTTGTAACTTAGAAAATA", true, None),
        ("GTCCTCACCGCCCCGACAGCGTCCTCACGG", false, None),
        ("TTAAAGAGAAGACATTACAATTTATTACAGA", true, Some(("L1_chr7_64347657_r", "FLANK", b'-', 29, 2, 31))),
        ("AAAAAAAAAAAAAATAGGCCAGGCACGGTGGCT", false, None),
        ("TGGCTCTGGGAAGCATCAGCCATGGAAGTTTAC", true, None),
        ("AGTGGAGCCACTCCTGCAGCAACTGCACCTGCTT", false, None),
        ("CGAAGCTCAGGAGACTCCGTTCGCACAAAACGCT", true, None),
        ("ACATGGAGGGCCTTCTCTAAGAAGTCGGCCAACGC", false, None),
        ("TTCTTTTTTTTTTTTTTCTCTCTTCTAAATTTTTTT", true, Some(("SVA_hg38_chr13_49377430_f", "FLANK", b'+', 23, 3, 27))),
        ("CCACTATCACAAGCTTCCTAGAGACATGAGGACTTG", false, None),
        ("TAGAAAAGTCAAAAGATAACAAGTGTTGGTGAGAACA", true, Some(("L1_chrX_94043586_r", "FLANK", b'-', 30, 1, 35))),
        ("GGGAATACTGCTGGACCAAGCACATCTGCAGGAGAAG", false, None),
        ("AGAAAATAATAATAATAATAATATTAAAAGGTAAACAG", true, Some(("L1_chr2_130377260_r", "FLANK", b'-', 25, 2, 27))),
        ("GCTGGCAGCTGGCAGCTGGCATGCTGAGACAAGGCATG", false, None),
        ("TTGTTTTTGGTTTTTGTTGTTTTTTTTGAGATATAGTCT", true, Some(("L1_chr14_52901282_u/-", "FLANK", b'+', 29, 3, 33))),
        ("ATATATAAAGAAAAGAATAGAACAAACAAACATAATGAA", false, None),
        ("GCACATGTACCCTAAAGCTTAGAGTATAATAAAAAAAAAA", true, Some(("L1_chr8_92647061_r", "FLANK", b'-', 39, 0, 40))),
        ("AGGGCCCTCTTTACCTTGCAGAGGCGAGCACCGTCAGGAG", false, None),
        ("TGCTCTCACCGTGGCCTGTCCACGGTCCAGGTCCATCTCA", true, None),
        ("CAAGGGTGAGCCACTGCACCCAGCCCCCATCCATTTCCTTT", false, None),
        ("CAACCAAAACATTAAAAAAACAGAAAAAAAAAAAACAAAAA", true, Some(("L1_chr1_164930326_r", "FLANK", b'+', 30, 10, 41))),
        ("TTTCAGCTCTTTCAGGAATTGCCACACTGCATTTCACAATG", false, None),
        ("TGCCTGTAATCCCAGCTACTCGGGAGGCTGAGTCAGCAGAAT", true, Some(("ALU_YB8", "ALU", b'+', 38, 1, 42))),
        ("TGGAGCAAACGATGTTACTCAAGGATGCACTTTACTGCTCTTCC", false, None),
        ("GCCATTTTCCCCACTGCAGGGGTAAGGTTCAGCCAGAGCACCCA", true, None),
        ("AAAATACTTCAAAAAATAATAATAACAATTCAACAATTAAAAATA", false, None),
        ("TGGTTTACAAGACAGTCCTGGAACTCAGCCTCCCTGAGGAGCCTC", true, None),
        ("TAGGTTAACCTCATGTATCAAAAGTTTTTATTCAGGTTTTGAGTGTTTTA", false, None),
        ("AAAAAATAGGGAAAGAAAAAAAACTATACAGTAAAAAAAAAAATAAAAGAAAA", true, Some(("SVA_chr7_74483462_r", "FLANK", b'+', 26, 10, 37))),
        ("AAAAAAAAAAAAAAAAAAAAAGGAAAATAGAAACAGGCCAACCAAAATGTCTGGGT", false, None),
        ("AAATCAAAAAAAAAAACCGAACAGCAAAAAAAAGGCAAGAAACAAGAAAATTAATAA", true, None),
        ("TTTTTTTTTTTTTTTTTTTTATTTTTTTTATTTTTTGTTTTTTTTTTTTTTTTTTTTTTT", false, None),
        ("GTATAGTATATACTGTATATACTATATAGTATAGTATATACTGTATATACTATATAGTATAGTATA", true, Some(("L1_chr4_95304645_f", "FLANK", b'-', 57, 1, 66))),
        ("AAAAAAATACAAAGAGGAAGAAAAAAAAAAAAAAAACGAAGTAAAAAAAAAAGATAAGTAAATGAGA", false, None),
        ("ACTGGCGGGGCGCTCGCAGCGCGGGCCGCCGGGAGGATGAATGTTTTATGGGATCAGGTGGTTTGACGCG", true, None),
        ("TGGTTGTATTGATTTATTTTTTTTTTGTTTTGTTTTTTTTTCATTTGTTTTTTCTTTTGTATACGTTGTATTTC", false, None),
        ("GGGCTAAGTGTATGGAGAGGGTCCCCAAACAAGAGCTGTGGCTGCCAAGTCTCTTCCTTCCTTAAGGAAGCTTGAGGTAGCTGG", true, None),
        ("CAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", false, None),
        ("ACAGTAAAATAAACAAAGCAAGACAAAGAAAAAAAAAAAAACTAACAAAATATAAATAAAAAAATGCATCAAAATAAAACAAAAAAAAAACACAAAAGAAAA", true, Some(("L1_chr2_172801429_f", "FLANK", b'+', 60, 31, 102))),
        ("ATAAACAAAAAACAAAAATCAACTACTTACGAACATAAAAAAAGACAAAAAAGAACTAGAAAAAAAATGCATAAATAAAAAAAAAATAACAAAATAAAAAAAAA", false, None),
        ("ACATAAAGATGAAAAAAAAGAAAGACAAAAAACGACAATTAAACGAAAACAAAAAAAAAAATAAAAAAACAAACTTAAAAAAACTATAGAAAAGAAAAAAAAAAAA", true, Some(("L1_chr2_169724675_f", "FLANK", b'+', 50, 11, 72))),
        ("CAGGAGAATGGTGTGAACCGGGGAGGCGGAGCTTGCAGTGAGCAGAGATAACGCCACTGCACTCTAGCCTGGGCGACAGGGTGAGACTGTCTCAAAAAAAAAAAAAAA", false, Some(("ALU_YA5", "ALU", b'+', 85, 0, 93))),
        ("AAAAAAAATAAAAAAAAAAAATTAGCCGGGTGTGGTGGAGTGCACCTGTAGTCCCAGCTGCTGGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCAGAGCTTGC", true, Some(("L1_chrX_33005657_f", "FLANK", b'-', 98, 8, 111))),
        ("TTCAAAAGAGTGACTCCTTGTTCTGCTGTTTCCCACCTCTCGCTACTGTACTTGACCAATCTTAAAGTGAATCTTATGCTCTGTGAATTATATATATATATATATTTTTTTTTTTTTTTTT", false, None),
        ("ATGCAGGCATGCAATCATGCCTGCAAGCAAGCACATCATGGAGAATGGGGTGTCCACCCCCTCAAGCATTTATCCTTTGAGTTGCAAACAATCCAGTTGCACTCTTAATTTTTTTTTTTTTTTTTTTTT", true, Some(("L1_chrX_25901741_r", "FLANK", b'-', 94, 3, 116))),
        ("TTTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTGATCTCCTGACCTCGTGATCCACCCGCCTTGGCCTCCCAAAGTGTTGGGATTACAGGCGTGAGCCACCGCGCCCAGCCAAAATTGC", false, Some(("ALU_Y", "ALU", b'-', 127, 1, 132))),
        ("GCATTAGGAGATATACCTAATGCTTAATGACGAGTTACTGGGTGCAGCACACCAGCATGGCACGTGTATACATATGTAACTAACCTGCACGTTGTGCACATATACCCTAAAACTTAAAGTATAATAATAATAAAATAAAATAAA", true, Some(("L1_chr3_123292579_f", "FLANK", b'+', 139, 0, 144))),
        ("GCCAAGGCAGACAGATCACGAGGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTTTGCTAAAAATACAAAAAAATTAGCCGGGCATGGTGGCGACACCTGTAATCCCAGCTACTCGGGAGGCTGAGGCGGGAGAATC", false, Some(("ALU_Y", "ALU", b'+', 131, 0, 150))),
        ("TAAGACAAGGACGGGTGCGGTGGCTCATGCCTATAATCCCAGCACTTTGGGAGGCCAAGGCGGGCAGATCACAAGGTCAGGAGATTGAGACCATCCTGGCTAACATGGTGAAACCCCATCTCTACTAAAAATACAAAAAAATTAGCCAGGT", true, Some(("SVA_chr2_110485709_r", "FLANK", b'+', 129, 13, 150))),
        ("AGAGATCAGATTCACCAAATAAAACCTTCAGCTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAATAAATCAAAAAAAACAAAAAAAAATAAATGGAAAAAAAAAAAAAAA", false, None),
        ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAACAAAAAATAATAAGAATCAAAAGAAAAGAAAGAAAAAGATAAAAAACAAAAAACAACATAAAAAAATAAAAAAAAAGGAACAAAAAAAAGAAAAAGAGAAAAAA", true, Some(("L1_chr22_49488441_r", "FLANK", b'-', 67, 66, 144))),
        ("GGAGGATCGCTTGAGCCTAGGAGTTCAAAACCAGCCTGAGCAACATAGTGAGACCCCATCTCAATTCAATAGTTAAAAAAATAAAATTAAAAGGCTATAGGTGAACTGCAGATGAAGAGATATATAGGGTGAGGTCAGGAAGGATCCTGTA", false, None),
        ("TTTTATTATACTCTAAGTTTTAGGGTACATGTGCACATTGTGCAGGTTAGTTACATATGTATACATGTGCCATGCTGGTGCGCTGCACCCACTAATGTGTCATCTAGCATTAGGTATATCTCCCAATGCTATCCCTCCCCCCTCCCCCGAC", true, Some(("L1_chr6_146374016_r", "FLANK", b'-', 151, 0, 151))),
        ("GTACCAGAGCTGTCTGTAGACCATGCAAAGGTAGCATCCTATTGGCAGCAGGTATGGAGGAGAGGAGGGAAGAGTTGCTGTTGGGTGACTGGCTGGCTGGAGGTAGAGCTAGTCCCTGAGTCTGGGAGGCTGAGAAAGGGCCAGGTTGAAT", false, None),
        ("ACGCATGCACCACCACGCCCAGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGATGGTCTCGATCTTCTGACCTTGTGATCCACCCACCTCGGCCTCCCAAAGTGTTGGGATTACAGGCGTGAACCACTGCCC", true, Some(("L1_chr7_142376221_r", "FLANK", b'-', 137, 2, 149))),
        ("CTCTACTAAAAATACAAAATTAGCCAGGCGTGGTGGTGCATGCCTGTAATCCCAGCTACTCAGGAGGCTGAGGTGGGAGAATCATTTGAACCCAGGAGGCAGAGATTGCAGTGAGCCAAGATCACATCACTGCACTCCAGCCTGGATGACA", false, Some(("ALU_Y", "ALU", b'+', 124, 0, 145))),
        ("TTCAGCATCTCCTCCTCCTCTCCTTCTCTGTGCAGCTGCCCTGCCGCACTAGTGAGGTTTCTGTTCCTCAGCACCACCGTTAGTTCCTGTCTCTGCCTCAGGGCCTTTGCACTTGCTATTGCTGCTGCCCAGTGCGCTCTGCACATGGCCA", true, Some(("L1_chr4_81689177_f", "FLANK", b'+', 61, 21, 102))),
        ("GCCTCCATCATTCATTCATTCACTCACTCACTCATTCAACAAACAAACAATGGGCACCTACCCACCACGGGCCAGGCCTGCCTGGGTTCTGGGGTTCAACTGTGAGTGCACAAACACACTAACAGCCCTTGGAGGACTCACAGTCCAGGGA", false, None),
        ("AGGAGACAGGGCAGAATCCCAGACCATCAGGCTCCTTTCTGCTATCTGTGGATCAAAGGCTGTGTGTGGATCCAACTCCCATACTCTGTGGATCCCTTCCTCCTGTCTCCTCTGCCACCCCCTACCTAAGCCCAGAAACTTCTCTCTGCCT", true, Some(("L1_chr10_123570750_f", "FLANK", b'+', 33, 3, 43))),
        ("GACGGGAGGATCCCTTGAGGCCAGTTCGAGGCTGCAGTGAGCTGCGCTTGTTGCACTGCACTCCAGCCTGGGTGACAGAGCAAGAACTTTCTCAAAAAAAAAAAAAAAAAAAAAAGGCCGGGTGCGGTGGTTCCCAGCACTTTGGGAGGCC", false, Some(("ALU_YA5", "ALU", b'+', 30, 53, 85))),
        ("GACCCCGGCTCACATGTTTCTCTCACTTGGCTGCAGCTCTGGCCCTGCCCCATCTCCCACCCACGTGTGTTGTCATGAAAATGAAACAAGGTGGCGCTGGTGAGAAGCAGGTGGAAGGCAGGGCTGCTGGCCACAGGCTGCTGTGAGGATC", true, Some(("L1_chr3_106273373_r", "FLANK", b'-', 56, 14, 84))),
        ("TGTTATATATAATATAATATTATATAATATTATAGAATATAATATTATATAATATAGAATATTATATTATATAATAATATATGATATTATTTTATATAATAATATATGATATTATTTTATATAATAATATATGATATTATATTATATAATA", false, None),
        ("TTTGTTTTTTTTGCTCTTTTTTTTTTTCTTTTTTTTTTGTTGTTTAATCTTCTTTTATCGTATTTTTTTTTCCTTGTTTATTTTCTTGTTTTTTTTATTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTC", true, Some(("SVA_chr6_151707445_f", "FLANK", b'-', 69, 67, 150))),
        ("TTGTTTGTTTATTTAGAGACAGAGTCTCGCTCTGTTGCCCAGGCAGTGGCGCGATCTCGGCTCACTGCAACCTCCACCTCCTGGGCTCAAGCAATTCTCCTGCTTCGGCGTCCCAAGTAGCTAGGACTACAGGTGCGCAGCACCACGCCCG", false, Some(("ALU_Y", "ALU", b'-', 120, 15, 151))),
        ("CTGACCACTGGGATATCATTTTTTGATTTGGAAATGATCTCAATGAGTTTGTGGTTTGGGTTTCCCGGGCTCTGGGCTCTGTTCCCTGGAGACCGCAGGGTGGGGCATGAAAAGCAAGATCCTAAAAAGCTCCTGGGGCGTTTCTCTGTAT", true, Some(("L1_chr20_11676071_f", "FLANK", b'+', 52, 79, 134))),
        ("TCATGGAGATGAAAGGACAGGCGAGCCCAGGAGGGCACACACCATGGCACGGCAGGCAGACAGACATAGACATAATGCCCATATGGCAGCATGATCTCTCACACACACACATGCTCACATGGGCAGCACACTCACACACACGTGCATATGG", false, None),
        ("GATGGAATCCTGTGAGAGAAGATAATGGGGATGGACAGCATCTAGGGATTCTAACTCTAGGCATCCCCACAGGAAGACAGAGCAAATCGGTTAGAGCGAATAATCAAAATTACAATAGGAGAAAAGTTTCCTTAGGTTGCAAAAAGGTCTG", true, Some(("L1_chr4_172875229_r", "FLANK", b'-', 24, 100, 125))),
        ("ATAATATATGACATTGTATTATATAATAATATATAACATTGTATTATATAATAATATATAATATTGTATTATTTAATAATATATAATATATTATATAATAATATATATTATATTATATATTAATATATAATATTATATAATAGAATATTAT", false, None),
        ("AGGCTGAGGAAGGAGAATCGCTTGAACCCGGGAGGCAGAGGTTGCAGTGAGCCGAGAATGCACCACTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCTGTAAGAAGAGAAAGAGAAAGAGAAGGTGAAAGAGAGGAGAAGGAGCAGGAG", true, Some(("SVA_chr10_100975879_f", "FLANK", b'-', 109, 0, 125))),
        ("TTAGCTCATCTTGGAATTGAATTTTTTTTTTTTTTTTTTTTTTTTGAGACAGTGTCGCTCTGTCGCCCAGGCTGAAGTGCAGTGGCATGATCTCGGTTCACTGCAGTCTCCACCTCCCGGGTTCAAGTGATTCTCCTGCCTCAGTCTCCTG", false, Some(("ALU_Y", "ALU", b'-', 92, 45, 149))),
        ("TTACTTTCCCTTCCACTCTCTTCTCCTCTCCTCTTCTCTCCCCTCCTCCCCTCCCCTCCCCTCCCCTCCCCTCCTCTCCCCCCTCCCCTCCTCCCCTCCCCTCCCCTCCTCTCCTCCCCTCCCCTCCCCTCCCCTCCTCTCCTCTCCTCCC", true, Some(("L1_chr7_92421318_f", "FLANK", b'+', 130, 6, 151))),
        ("CTTCCACAGTGTTTGTGTCCCTGGGTACTTAAAGATTAGGGAGTGGTGATGACTCTTAACGAGCATGCTGCCTTCAAGCATCTGTTTAACAAAGCACATCTTGCACCGCCCTTAATCCATTTAACCCTGAGTGGACACAGCACATGTTTCA", false, Some(("SVA_E", "SVA", b'-', 149, 0, 151))),
        ("CCCTATTTAACAAATGGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATCAATTCAAGATGGATTAAAGATTTAAACGTTAGACCTAAAACCATAAAAACCCTAGAAGAAAA", true, Some(("L1_chr8_92647061_r", "FLANK", b'+', 150, 0, 151))),
        ("ATTGTGCCTCTGCACTCCAGCCTGGGCAACGGAGTGAGACTCTGTCAAAAAAAAAAAAAAAAAGAAAATATCACTTATCTACCAGTTGTCTTTCTGTGGCTTGGGTTTGGTGTGCTCTCTGTTCAAAGGACGGGGTAAATTGCTAAGAGAA", false, Some(("ALU_YA5", "ALU", b'+', 36, 5, 46))),
        ("CGACCATCTTGCTCACACACAGACACAGACCACAGTTCACACATCTCTGTCTCCCTTCAACCCGGTACAACCCCCTCCTACCTGTTCACAGAGCACCACATACACTCGGACGCTTGGCAGGACTTAGCATAGCTACACTGTTACATTTAAT", true, Some(("L1_chr6_158789722_r", "FLANK", b'-', 34, 12, 97))),
        ("TTCACTTGTTTATCTGCTGACCTTCCCTCCACTATTGTCCTATGACCCTGCCAAATCCCCCTCTGCGAGAAACACCCAAGAATGATCAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", false, Some(("SVA_F", "SVA", b'+', 97, 0, 101))),
        ("CTTCAAGTGATCCTCCCGCCTCAGCCTCCCAAAATGATGGGATTACAGGCATGAGTTACCCCATCTGGCGCAAGAATATTAGTGATTTCTATCATAATCATTGGTTATTTACTTATTCCAACAGTTAATGAGCACCTGCTATGTGCTGGGA", true, Some(("SVA_chr10_100774933_f", "FLANK", b'+', 81, 7, 113))),
        ("GGCTAACACGGTGATATCCCATCTCTACTAAAATACAAAAAATTAGCTGGGCATGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGCGGCAGGAGAATGGCGTGAACCTGGGAGGCAGAGCTTGCAGTGAGCCTAGACTGCGCC", false, Some(("ALU_Y", "ALU", b'+', 139, 0, 151))),
        ("CCACAACCAATTTTTTAAAAAGAAAGAAAAATAAAATATGCTATTCATAGTGCATAATCAACATAATCCCACTAGCAAGTACAGTATCGTACTAAATGAACCATTTGATTTGCAAATTTTCCAACGTCGTGAATCTCACAACTTTTGTGCC", true, Some(("L1_chr4_134564267_r", "FLANK", b'+', 64, 28, 124))),
        ("CAATTCCCACCTATGAGTGAGAATATGCGGTGTTTGGTTTTTTGTTCTTGCGATAGTTTACTGAGAATGATGGTTTCCAATTTCATCCATGTCCCTACAAAGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGTATA", false, Some(("L1HS", "L1", b'-', 150, 0, 151))),
        ("TGATGGTGATGATAATGATAATGATGGTGATGGAGGTGATAATGATGATGGTGGTGGTGATGGTGATGATGAAGATAATGATGGTGATGGTGGTGATAATGATGATGGTGGTGGTGATGGTGATGATAATGATAATGATAATGATGGTGAT", true, Some(("L1_chr2_174185677_r", "FLANK", b'+', 132, 0, 151))),
        ("TTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAGCACAACTAAATCCACACTGCACACAGACGCAGACAGAAAGCCTTCAAGTGGCTCTGTTTTCTGCTCCCTGCCTTGCCAGGTCCACAAGCAGAG", false, None),
        ("TCATCTCATCATTTCATCATTTCATCTCATTTCTTCTCATCATTTCATCTCATCATTTTATCTCATTTCATCTCATCTCATTTCAATTTCATTTCATTATTTCATTTCATTTCACTTCATTTCATCTCATCATTTTATCTCATCTCATTTC", true, Some(("L1_chr14_8517112_f", "FLANK", b'+', 145, 0, 151))),
        ("AGGGTTGTGCTGTTGTTATCCCCATTTTATGGAGGAAGAAACAGAGGCTCCGTGGGTCCAGAAATCCACCCAGACTCCAACGCAGGGCCGATGACACTCAGCGGGCATCAGCAGAATGGCTGAGAGCCCCAGGCATCCAGTGAGAGGCCCC", false, None),
        ("TTGGCTTAGCAGCTCCTATTTTCTCCATCTAACAGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCGGGATCTCGGCTCACTGCAAGCTCCGCCTCCCGGGTTCAC", true, Some(("L1_chr13_30525317_f", "FLANK", b'-', 113, 38, 151))),
        ("ATTTTTTTCCAAGTAGATATTGGAAATGATAGACTGACAGTCAGGAGAATTTGGGAATGGATTATAAAACAGGTGGCTTGGGAATACTTACGAAAGGAAGTGGTGACCTCTATGGATCCCCTGTATTCAGTTAGAATGTGTCATGCCAAAT", false, None),
        ("GCTCTTAGGCAAGAGGCCTGGGGAGCCTGCTCACCCCTTCCACCACGTGAGGACAAGGCTATAAGGCATCTATGAGAAAATGGCCCTCATCAGACACCGAATCTGCAGGCACATTCATCTTGGATTTCCCAGCCTCCAGAACTGTGAACAA", true, Some(("L1_chr12_44066939_r", "FLANK", b'+', 62, 74, 147))),
        ("GGGTTCAGGTGATTCTCCTGCCTCAGCCTCTTGAATAAATGGGATTACAGGCGCCCGCCACCATGCCTGGCTAATGTTTGTATTTTTAGTAGAGATGGGATTTCACCATGTTGGTCAGGCTGATCTCAAACTCCTGACGTCAAATGATTGG", false, Some(("ALU_YB8", "ALU", b'-', 119, 0, 142))),
        ("ACACAGAAAATAGCAGAGTCTTTTAGCTCCTGCCAAAAGTCTGCACACAGGCCCTTTATTTACGGCTAAGTGTAAATGACTGATACAGCCAACTAATAATCAACATGATAATCATTCTTATATAGCTCTGTTTTTCAGAGGAGGAAACGAT", true, Some(("L1_chr14_57324334_r", "FLANK", b'-', 20, 80, 100))),
        ("ATGTTTTTTTTTTTGCTTATTTTTTTTTTTTTTCTTTTTTTTTTTTCGTTTTGTTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTG", false, None),
        ("CTCACTGCAACCTCTGCCTCCTGGGTTCAAGTGATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGACTACAGGCATGCGCCACCACGCCAGGCTAATTTTAGTACTTTTGTAGAGACCTGGTTTCGCCATGTTGGCCAGGCTGGTCTTGAA", true, Some(("L1_chr6_24682781_r", "FLANK", b'+', 140, 0, 151))),
        ("GATGATAGGGGCCGGGCGCGGTGGCTCACGCCTGTAACCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCATGAGGTCAGGAGATCGAGACCATCCTGGCTAACAAGGTGAAACCCCGTCTCTACTAAAAATACAAACAATTAGCCGGGC", false, Some(("ALU_YB8", "ALU", b'+', 140, 9, 151))),
        ("TTTCCGCTCAGGACGCTGCGGGTGATGGGGAAGACTCTGAGGGAGTCACTGCTATACCCCTGACTCGGCTGCCCATTTCCTCATCTGCTTGGCCGACGCTTCTTCTCTGGCACCAGCTGCATCCTGGTTTTTTTCTGGCTTTCTTCAAAAA", true, Some(("L1_chr4_192496713_f", "FLANK", b'+', 20, 28, 49))),
        ("CCCTTTCTCTTTCCCTTTTCCTTTTCCTTTCCCCTTTCCCCTTTCTTTTTTCCTTTCCTTTCCTTTCCCTTCTGCTTTCCTTTCCCTTCCGCTTTCCTTTCCTTTCCCTTCCGCTTTCCTTTCCTTTCCCTTCCACTCTCTTCTCCTCTCC", false, None),
        ("TGGCCAGGCTAGTCTCGAATTCCTGGCCTCAAGTGATCTGCCTGCTTTGGCCTCCTAAAGTGCTGGAATTACTGGCATGAGCCATTGCGCCCAGCTCAACAACCCTTATCTCCCAGGTTCTCACACCGACTTCTTGAGTTACGTAACGAAG", true, Some(("SVA_chr16_28977873_r", "FLANK", b'-', 80, 0, 92))),
        ("GGCCTCCCAAAGTGCTGGGATTACAGGTGTGAGCCACTGCACCTGGCCCATCATCTCTTAAATATCAAGGTAGATGGTTTATATTAGTATTAGCCTACATTTTTTTTAACATATAAAGCATCATATTAATGTTCAATTTCCTGGTTAAATA", false, Some(("ALU_YB8", "ALU", b'-', 44, 0, 48))),
        ("CGAGTAGCTGGGACCACAGGCGTGTGCCACCACACCTGGCTAATTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCACCCAGAGTGCACCGGGCTAGAGTGCAATGGCGCGATCTTGTCTCACTGCCACCTCCACTTCCAGGT", true, Some(("SVA_hg38_chr9_81709537_f", "FLANK", b'+', 117, 2, 150))),
        ("GTGGCCACAAATACTTCTTTCCTCAATGAACTGTTAAAAGGACTCACTGAGAAAAGTATGAGAAAATGACCTCCAGGAGTTTGAGACCAGCCTGACCAACATGGTGAAACCCCGTTTCCACTAAAAATACAAAATAATTAGCTGGGCGTGG", false, Some(("ALU_Y", "ALU", b'+', 68, 73, 151))),
        ("GCCCTCGCTATGCTCTCCGGGTCTGTGCTGAGGGGAACGCAGCTCCGCCCTCGCAAAGGCACACAGCGCCGGCGTGGCGGAGAGGCGGACAGCGGCGGAGAGGCGGACAGCGGCGGCGCGGCGGAGAGGCGGACAGCGGCGGAGAGGCGGA", true, None),
        ("ATCACAAAGGCGACCTCTTTCGTCACATCTGTGTGAGAGGAGCTGGACTTTCTCAGATATGGATGAGCCTCACGACGCACATCCTTTCGTCCCTTGCTTGTCTACCAAGTCCCTATGTGCTACAGGAGGCAGTGGTACAGGCCAGCCGCTC", false, None),
        ("CAGAATGCCATGTTTGTGAGTGTGATCATTCATTCATCATATGGCTGAGGGCAACTGTTAGTAATTGAAGAATCTACGCAACCTCTAAGCTACTTCAATGAGGAAGAGACTATTGAGCTGTGAGGATCTTGGCTGTAGCCAGACTCAAGGT", true, None),
        ("TAGAAGGTTCTATTTCACCTCATGGGTTTACACCCTGGCATCATTTGCTTTCTGAGACATGATTGTCTCCAGGCCAGGGTTGATCCAAGCTCCTTCTTGCTCATGGAGGCAGAGAGCTCTTATTTGTTTTCCTCCCTTTGTAGTGTAATCT", false, None),
        ("TTCTGCACCAATGCTTTAAGATATAACCTGGTAAGTAATAACTCCTCAAAATGCAAAGCCTAATTTACTGCAAAATCTTGTTTCCTTATATTCTAAAAGCCTATAGCTTTTAGAATTCACCAGTTCTAAGGCTGATGAAATGAGAAAGCAG", true, None),
        ("ATTAAGTGAGTTCATGAAAGCCCTCGGTGTGTGGCATGCGATAACCACTACCTAAGCATAAGGTATGTTTACATGCACTGGTTGTTTACCACGACGCACCAGCATTCATTCAAGTTATTTTCTTCCAGAGTCTGCCTGGAATGTAGAACCG", false, None),
        ("TCTTTTTTCCTCGTTTCGTGCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", true, None),
        ("AGCACAAATCTAAAATTGTAACACAGGCGACTATCAATCTTTTCCTCTATGATTTGAGCATTTTGTGTTTTGCTTAAAAACTCAGCCCTTGTCCCAAGGTCATAGAGATATTCTCCCGTATTTTCTACTGATCCATCCCCCATATTTAGGT", false, None),
        ("CGGATACACGACTCTAGCCACCAGTGTGACCCTGTTAAAAGCCTCGGAAGTGGAAGAGATTCTGGATGGCAACGATGAGAAGTACAAGGCTGTGTCCATCAGCACAGAGCCCCCCACCTACCTCAGGTAATGCGTTCCTGGCCAGGGCATC", true, None),
    ];
}
