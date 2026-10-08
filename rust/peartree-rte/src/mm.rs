//! mappy (minimap2 python binding) equivalent. FOUNDATION (implemented).
//!
//! Replicates `mappy.Aligner.__cinit__` / `Aligner.map` of mappy 2.31 (minimap2/python/mappy.pyx
//! + cmappy.h) on the minimap2 2.30 C library bundled by the `minimap2` crate:
//!
//!   * `mm_set_opt(NULL)`, then `mm_set_opt(preset)` when a preset is given,
//!     `map_opt.flag |= MM_F_CIGAR` (always align), `idx_opt.batch_size = 0x7fff_ffff_ffff_ffff`
//!     (uni-part index), then only the options python passed (k, w, min_cnt, min_chain_score,
//!     min_dp_score -> `min_dp_max`, best_n);
//!   * from a path (`Aligner(fn)`): `mm_idx_reader_open` / `mm_idx_reader_read(r, 3)` (first part
//!     only) / `mm_mapopt_update` / `mm_idx_index_name`;
//!   * from a sequence (`Aligner(seq=...)`): `mappy_idx_seq` = `mm_idx_str(w, k, flag & 1,
//!     bucket_bits, 1, [seq], ["N/A"])`, `mm_mapopt_update`, then `map_opt.mid_occ = 1000`;
//!   * `map(seq)`: `mm_map(idx, strlen(seq), seq, .., name=NULL)`; every region becomes one
//!     [`Hit`] in minimap2 order (primary and secondary alike, exactly like mappy).
//!
//! Python call sites and their options (keep these presets identical when porting):
//!   library.py   `RteLibrary.aligner(kind)`        [`MapOpts::LIBRARY`] from a FASTA path
//!                (consensus / intact (concatenated temp FASTA) / flanks3 / flanks5)
//!   assembly.py  `SiteContext.aligner()`           [`MapOpts::LOCAL`] from `seq=` (window)
//!   assembly.py  `SiteContext.wide_aligner()`      [`MapOpts::SR`] from `seq=` (wide window)
//!   annotator.py `_premrna_fn` (gene model window) [`MapOpts::SR`] from `seq=`
//!   transduction.py `MappyLocator`                 [`MapOpts::SR`] from a path (.mmi / FASTA)
//!
//! Thread safety: an [`Aligner`] is immutable after construction; each thread maps with its own
//! `mm_tbuf_t` (thread-local), so `map` may be called from any number of rayon workers.

use minimap2::ffi;
use std::cell::RefCell;
use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::path::Path;
use std::sync::Arc;

/// The mappy keyword arguments used by tools/rte (None = not passed = minimap2 default/preset).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MapOpts {
    pub preset: Option<&'static str>,
    pub k: Option<i16>,
    pub w: Option<i16>,
    pub min_cnt: Option<i32>,
    pub min_chain_score: Option<i32>,
    pub min_dp_score: Option<i32>,
    pub best_n: Option<i32>,
}

impl MapOpts {
    /// library.py `RteLibrary.MAP_OPTS = dict(k=11, w=3, min_chain_score=18, min_dp_score=25, best_n=8)`
    pub const LIBRARY: MapOpts = MapOpts {
        preset: None,
        k: Some(11),
        w: Some(3),
        min_cnt: None,
        min_chain_score: Some(18),
        min_dp_score: Some(25),
        best_n: Some(8),
    };
    /// assembly.py `SiteContext.aligner()`: `k=11, w=3, min_chain_score=18, min_dp_score=25, best_n=6`
    pub const LOCAL: MapOpts = MapOpts { best_n: Some(6), ..MapOpts::LIBRARY };
    /// `preset="sr"` (wide window, pre-mRNA window, MappyLocator)
    pub const SR: MapOpts = MapOpts {
        preset: Some("sr"),
        k: None,
        w: None,
        min_cnt: None,
        min_chain_score: None,
        min_dp_score: None,
        best_n: None,
    };
}

/// One mappy `Alignment` (all fields tools/rte reads).
#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    /// `h.ctg` (reference sequence name; "N/A" for a `seq=` aligner)
    pub ctg: Arc<str>,
    pub ctg_len: i64,
    pub r_st: i64,
    pub r_en: i64,
    /// +1 / -1
    pub strand: i32,
    pub q_st: i64,
    pub q_en: i64,
    pub mapq: u32,
    /// `[[len, op], ...]`, op 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7='=' 8=X
    pub cigar: Vec<(u32, u32)>,
    pub is_primary: bool,
    pub mlen: i64,
    pub blen: i64,
    pub nm: i64,
}

impl Hit {
    /// `h.mlen / h.blen if h.blen else 0.0` (assembly._hit_identity)
    pub fn identity(&self) -> f64 {
        if self.blen != 0 {
            self.mlen as f64 / self.blen as f64
        } else {
            0.0
        }
    }
}

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

/// A built minimap2 index + the mapping options mappy would use with it.
pub struct Aligner {
    idx: *mut ffi::mm_idx_t,
    mapopt: ffi::mm_mapopt_t,
    names: Vec<Arc<str>>,
    lens: Vec<i64>,
}

// SAFETY: after construction the index and options are only read; mm_map is thread-safe given a
// distinct mm_tbuf_t per concurrent call (minimap2 API contract), which `with_tbuf` guarantees.
unsafe impl Send for Aligner {}
unsafe impl Sync for Aligner {}

fn init_opts(o: &MapOpts) -> (ffi::mm_idxopt_t, ffi::mm_mapopt_t) {
    // SAFETY: zeroed POD option structs are immediately initialised by mm_set_opt.
    let mut io: ffi::mm_idxopt_t = unsafe { std::mem::zeroed() };
    let mut mo: ffi::mm_mapopt_t = unsafe { std::mem::zeroed() };
    unsafe { ffi::mm_set_opt(std::ptr::null(), &mut io, &mut mo) };
    if let Some(p) = o.preset {
        let c = CString::new(p).expect("preset name");
        let rc = unsafe { ffi::mm_set_opt(c.as_ptr(), &mut io, &mut mo) };
        assert!(rc == 0, "minimap2 rejected preset {p:?}");
    }
    mo.flag |= ffi::MM_F_CIGAR as i64; // always perform alignment
    io.batch_size = 0x7fff_ffff_ffff_ffff; // always build a uni-part index
    if let Some(k) = o.k {
        io.k = k;
    }
    if let Some(w) = o.w {
        io.w = w;
    }
    if let Some(v) = o.min_cnt {
        mo.min_cnt = v;
    }
    if let Some(v) = o.min_chain_score {
        mo.min_chain_score = v;
    }
    if let Some(v) = o.min_dp_score {
        mo.min_dp_max = v;
    }
    if let Some(v) = o.best_n {
        mo.best_n = v;
    }
    (io, mo)
}

impl Aligner {
    fn finish(idx: *mut ffi::mm_idx_t, mapopt: ffi::mm_mapopt_t) -> Aligner {
        // SAFETY: idx is a valid index; names are NUL-terminated C strings owned by it.
        let n = unsafe { (*idx).n_seq } as usize;
        let mut names = Vec::with_capacity(n);
        let mut lens = Vec::with_capacity(n);
        for i in 0..n {
            let s = unsafe { &*(*idx).seq.add(i) };
            let name = if s.name.is_null() {
                String::new()
            } else {
                unsafe { CStr::from_ptr(s.name) }.to_string_lossy().into_owned()
            };
            names.push(Arc::from(name.as_str()));
            lens.push(s.len as i64);
        }
        Aligner { idx, mapopt, names, lens }
    }

    /// `mappy.Aligner(path, **opts)`; None when python's `not aligner` (no index could be read:
    /// missing / empty file). Reads FASTA (plain or gzipped) or a prebuilt `.mmi`.
    pub fn from_path(path: &Path, opts: MapOpts) -> Option<Aligner> {
        use std::os::unix::ffi::OsStrExt;
        let cpath = CString::new(path.as_os_str().as_bytes()).ok()?;
        let (io, mut mo) = init_opts(&opts);
        // SAFETY: valid C string and option structs; the reader is closed below.
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
        Some(Aligner::finish(idx, mo))
    }

    /// `mappy.Aligner(seq=seq, **opts)` (one reference sequence named "N/A").
    pub fn from_seq(seq: &[u8], opts: MapOpts) -> Option<Aligner> {
        let (io, mut mo) = init_opts(&opts);
        // mappy_idx_seq: calloc(len + 1) copy, so an embedded NUL truncates exactly like C
        let mut s = seq.to_vec();
        s.push(0);
        let name = CString::new("N/A").unwrap();
        let mut seqs: [*const c_char; 1] = [s.as_ptr() as *const c_char];
        let mut names: [*const c_char; 1] = [name.as_ptr()];
        // SAFETY: one NUL-terminated sequence and name, valid for the call (mm_idx_str copies).
        let idx = unsafe {
            ffi::mm_idx_str(
                io.w as i32,
                io.k as i32,
                (io.flag & 1) as i32,
                io.bucket_bits as i32,
                1,
                seqs.as_mut_ptr(),
                names.as_mut_ptr(),
            )
        };
        if idx.is_null() {
            return None;
        }
        unsafe { ffi::mm_mapopt_update(&mut mo, idx) };
        mo.mid_occ = 1000; // mappy: don't filter high-occ seeds
        Some(Aligner::finish(idx, mo))
    }

    /// Reference sequence names in index order.
    pub fn names(&self) -> &[Arc<str>] {
        &self.names
    }

    /// `for h in aligner.map(seq)`: all regions in minimap2 order.
    pub fn map(&self, seq: &[u8]) -> Vec<Hit> {
        // mappy: `if ((map_opt.flag & 4) and (idx.flag & 2)): return` (index without sequences)
        // SAFETY: idx valid for the aligner's lifetime
        if (self.mapopt.flag & 4) != 0 && (unsafe { (*self.idx).flag } & 2) != 0 {
            return Vec::new();
        }
        // mm_map(.., strlen(seq), ..): an embedded NUL ends the query like in C
        let qlen = seq.iter().position(|&b| b == 0).unwrap_or(seq.len());
        if qlen == 0 {
            return Vec::new();
        }
        let mut n_regs: i32 = 0;
        let regs = with_tbuf(|b| unsafe {
            // SAFETY: seq valid for qlen bytes; b is this thread's buffer; idx/mapopt read-only
            ffi::mm_map(self.idx, qlen as i32, seq.as_ptr() as *const c_char, &mut n_regs, b, &self.mapopt, std::ptr::null())
        });
        if regs.is_null() {
            return Vec::new();
        }
        let n = n_regs.max(0) as usize;
        let mut hits = Vec::with_capacity(n);
        for i in 0..n {
            // SAFETY: regs has n_regs elements; each `p` was malloc'd by minimap2.
            let r = unsafe { &*regs.add(i) };
            let mut cigar = Vec::new();
            let mut n_ambi = 0i64;
            if !r.p.is_null() {
                let p = unsafe { &*r.p };
                n_ambi = p.n_ambi() as i64;
                let cg = unsafe { p.cigar.as_slice(p.n_cigar as usize) };
                cigar.extend(cg.iter().map(|&c| (c >> 4, c & 0xf)));
            }
            let rid = r.rid as usize;
            hits.push(Hit {
                ctg: self.names[rid].clone(),
                ctg_len: self.lens[rid],
                r_st: r.rs as i64,
                r_en: r.re as i64,
                strand: if r.rev() != 0 { -1 } else { 1 },
                q_st: r.qs as i64,
                q_en: r.qe as i64,
                mapq: r.mapq(),
                cigar,
                is_primary: r.id == r.parent,
                mlen: r.mlen as i64,
                blen: r.blen as i64,
                nm: r.blen as i64 - r.mlen as i64 + n_ambi,
            });
            unsafe { ffi::free(r.p as *mut _) };
        }
        unsafe { ffi::free(regs as *mut _) };
        hits
    }
}

impl std::fmt::Debug for Aligner {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Aligner({} seqs)", self.names.len())
    }
}

impl Drop for Aligner {
    fn drop(&mut self) {
        // SAFETY: built by mm_idx_reader_read / mm_idx_str, destroyed once.
        unsafe { ffi::mm_idx_destroy(self.idx) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seq_aligner_finds_itself() {
        // a 300 bp pseudo-random reference; a 100 bp slice maps back exactly
        let mut x: u64 = 12345;
        let reference: Vec<u8> = (0..300)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect();
        let al = Aligner::from_seq(&reference, MapOpts::LOCAL).unwrap();
        assert_eq!(&*al.names()[0], "N/A");
        let hits = al.map(&reference[100..200]);
        assert!(!hits.is_empty());
        let h = &hits[0];
        assert_eq!((h.r_st, h.r_en, h.q_st, h.q_en, h.strand), (100, 200, 0, 100, 1));
        assert_eq!(h.mlen, 100);
        assert!(h.is_primary);
        let rc = crate::sequtil::rc(&reference[150..250]);
        let hits = al.map(&rc);
        assert_eq!(hits[0].strand, -1);
        assert_eq!((hits[0].r_st, hits[0].r_en), (150, 250));
    }
}
