//! Sidecar row routing + detached reads storage (python `_ShardStore`). OWNER: P5.
//!
//! Mirrors src/combine_insertions_evidence.py:838-980 in PURPOSE, not in file format: the
//! on-disk layout is private (the python shard files are not an output). Contract:
//!   * one streaming pass per sidecar (files in parallel), keeping only rows whose
//!     (file, parsed locus) is a member locus of some insertion, routed to every chunk that
//!     needs them; within a (member, side) the row order is the sidecar line order;
//!   * `chunk_rows(c)` gives a chunk's rows keyed by (member, side);
//!   * `lookup(member, side)` serves the parent-side re-evaluation in absorb_one_sided;
//!   * evaluated chunks move each record's rendered reads (FASTA text) to disk
//!     (`put_reads` -> `ReadsRef`), so peak memory is bounded by the chunks in flight;
//!   * the scratch directory `<stem>.evidence_shards` is created fresh (rm -rf first) and
//!     removed by `cleanup()` after the evidence outputs are written.
//! Sidecar path rule (python `sidecar_path`): `<txt.gz>.evidence.tsv.gz` (Rust discovery), else
//! `<stem>.evidence.tsv.gz` when only that one exists. SPEC.md §4.0, §7.
//!
//! Layout as built (PLAN.md "Memory / parallelism design"):
//!   * `route`: (file, locus) -> chunk(s) listing it (one u32 per member; a side table for the
//!     rare member shared by several chunks). Held for the whole run (absorb lookups).
//!   * rows: one scratch file per sidecar (`rows.<FileId>.bin`). The routing task of a sidecar
//!     appends wanted raw lines to per-chunk buffers; when the buffered total reaches
//!     `FLUSH_BYTES` every non-empty buffer is written as one raw-deflate (level 1) block,
//!     recording `blocks[chunk] = [(offset, len)]`. Peak per routing task: ~FLUSH_BYTES. Blocks
//!     of a chunk are in flush order, so a file's lines of a chunk stay in sidecar order.
//!   * reads: `put_reads` appends to the chunk's in-memory buffer (offset relative to the
//!     chunk); `finish_chunk` writes the whole buffer to `reads.bin` (one contiguous run per
//!     chunk -> sequential reads in output order) and frees it.
//!   * python routing quirk kept: the locus field is taken from the raw line (with its newline),
//!     so a `locus` column that is LAST in the header only matches on an unterminated final line
//!     (python `line.split("\t", li + 1)` keeps the "\n"). Matching is exact string equality
//!     with the member id (tokens must be canonical: no sign / leading zero that `Tok` would
//!     normalise).

use crate::context::Ctx;
use crate::evidence::junction::ReadsRef;
use crate::evidence::row::{EvidenceRow, SidecarHeader};
use crate::model::{ContigId, FileId, Insertion, Interner, LocusKey, Member, Side, Tok, TokKind};
use flate2::read::DeflateDecoder;
use flate2::write::DeflateEncoder;
use flate2::Compression;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::VecDeque;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

/// Buffered wanted-row bytes per routing task before its chunk buffers are flushed as blocks.
pub const FLUSH_BYTES: usize = 16 << 20;
/// python `_rows_cache` size (parent-side lookups).
const LRU_CHUNKS: usize = 4;
/// high bit of a `route` value: index into `multi` instead of a chunk id
const MULTI: u32 = 0x8000_0000;

/// python `sidecar_path(txt_gz)`.
pub fn sidecar_path(txt_gz: &str) -> PathBuf {
    let rust = format!("{txt_gz}.evidence.tsv.gz");
    let stem = txt_gz.strip_suffix(".txt.gz").unwrap_or(txt_gz);
    let alt = format!("{stem}.evidence.tsv.gz");
    if !Path::new(&rust).exists() && Path::new(&alt).exists() {
        PathBuf::from(alt)
    } else {
        PathBuf::from(rust)
    }
}

/// python `_ShardStore.build` chunking: `max(1, min(256, max(4*threads, n//2000)))` chunks of
/// `ceil(n / n_chunks)` consecutive insertions; `[(0, 0)]` when n = 0.
pub fn chunk_ranges(n: usize, threads: usize) -> Vec<(usize, usize)> {
    let n_chunks = 1.max(256.min((4 * threads.max(1)).max(n / 2000)));
    let size = 1.max(n.div_ceil(n_chunks));
    let v: Vec<(usize, usize)> = (0..n).step_by(size).map(|lo| (lo, n.min(lo + size))).collect();
    if v.is_empty() {
        vec![(0, 0)]
    } else {
        v
    }
}

/// Member loci per chunk above which a python chunk is split (memory bound: a chunk's parsed
/// rows are held while it is judged, `threads` chunks at a time; rows scale with members).
pub const MAX_CHUNK_MEMBERS: usize = 2000;

/// `chunk_ranges`, then every chunk holding more than `max_members` member loci split into
/// consecutive pieces of at most `max_members` (a single insertion is never split). Outputs
/// never depend on the chunking; this only bounds the rows in flight as sidecars grow.
pub fn chunk_ranges_bounded(members: &[Vec<Member>], threads: usize, max_members: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for (lo, hi) in chunk_ranges(members.len(), threads) {
        let mut start = lo;
        let mut acc = 0usize;
        for (k, ms) in members[lo..hi].iter().enumerate() {
            let k = lo + k;
            if acc > 0 && acc + ms.len() > max_members {
                out.push((start, k));
                start = k;
                acc = 0;
            }
            acc += ms.len();
        }
        out.push((start, hi));
    }
    out
}

/// Rows of one chunk, keyed by (member, side), each list in sidecar order.
#[derive(Default)]
pub struct ChunkRows {
    map: FxHashMap<(Member, Side), Vec<EvidenceRow>>,
}

impl ChunkRows {
    /// python `lookup(m, side)` -> rows (empty when none).
    pub fn get(&self, m: &Member, side: Side) -> &[EvidenceRow] {
        self.map.get(&(*m, side)).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// append a row under ((row.file, row.locus), row.side)
    pub fn push(&mut self, r: EvidenceRow) {
        self.map.entry(((r.file, r.locus), r.side)).or_default().push(r);
    }

    pub fn n_rows(&self) -> usize {
        self.map.values().map(|v| v.len()).sum()
    }
}

/// member -> chunk(s) whose insertions list it (python `route`; first = `locus_chunk`).
#[derive(Default)]
struct Route {
    map: FxHashMap<Member, u32>,
    multi: Vec<Vec<u32>>,
}

impl Route {
    /// chunks must be added in non-decreasing order (python `if not lst or lst[-1] != c`)
    fn add(&mut self, m: Member, c: u32) {
        match self.map.get_mut(&m) {
            None => {
                self.map.insert(m, c);
            }
            Some(v) if *v & MULTI != 0 => {
                let l = &mut self.multi[(*v & !MULTI) as usize];
                if *l.last().unwrap() != c {
                    l.push(c);
                }
            }
            Some(v) => {
                if *v != c {
                    self.multi.push(vec![*v, c]);
                    *v = MULTI | (self.multi.len() - 1) as u32;
                }
            }
        }
    }

    fn chunks(&self, m: &Member) -> &[u32] {
        match self.map.get(m) {
            None => &[],
            Some(v) if *v & MULTI != 0 => &self.multi[(*v & !MULTI) as usize],
            Some(v) => std::slice::from_ref(v),
        }
    }
}

/// One sidecar's scratch file.
struct ScratchFile {
    file: FileId,
    path: PathBuf,
    /// raw header line (parsed into a `SidecarHeader` on first use)
    header: String,
    /// blocks[chunk] = [(offset, compressed len)]
    blocks: Vec<Vec<(u64, u32)>>,
}

/// Rendered reads of evaluated records.
struct Reads {
    path: PathBuf,
    /// (append handle, end offset)
    writer: Mutex<(Option<File>, u64)>,
    reader: OnceLock<File>,
    /// in-memory text of chunks not yet finished
    pending: Vec<Mutex<Vec<u8>>>,
    /// file offset of each finished chunk's text
    base: Vec<OnceLock<u64>>,
}

pub struct Store {
    dir: PathBuf,
    have: Vec<FileId>,
    n_rows: u64,
    route: Route,
    n_chunks: usize,
    files: Vec<ScratchFile>,
    headers: OnceLock<Result<Vec<SidecarHeader>, String>>,
    reads: Reads,
    lru: Mutex<VecDeque<(usize, Arc<ChunkRows>)>>,
    /// per-member row index for absorb lookups (`prefetch_members`); None = chunk LRU only
    member_ix: Option<MemberIx>,
}

/// `prefetch_members` result: every wanted member's routed lines (raw, "\n"-terminated, in
/// chunk-line order) as one raw-deflate blob in `lookup.bin`. `len == 0` = member has no rows.
struct MemberIx {
    fh: File,
    map: FxHashMap<Member, (u64, u32, u32)>,
}

impl Store {
    /// Route every wanted sidecar row (python `_ShardStore.build`). `members[k]` = member loci of
    /// `insertions[k]`. Returns the store and the chunks as half-open insertion index ranges
    /// (any chunking is allowed: outputs never depend on it). Routing runs on the current rayon
    /// pool (one task per sidecar).
    pub fn build(dir: &Path, accepted: &[FileId], _insertions: &[Insertion], members: &[Vec<Member>], ctx: &Ctx) -> Result<(Store, Vec<(usize, usize)>), String> {
        let sidecars: Vec<(FileId, PathBuf)> = accepted
            .iter()
            .map(|&f| (f, sidecar_path(&ctx.files[f as usize].path)))
            .filter(|(_, p)| p.exists())
            .collect();
        let chunks = chunk_ranges_bounded(members, ctx.threads, MAX_CHUNK_MEMBERS);
        let st = Store::build_from(dir, &sidecars, members, &chunks, &ctx.contigs, FLUSH_BYTES)?;
        Ok((st, chunks))
    }

    /// `build` without a Ctx: `sidecars` = (file, existing sidecar path) in input order;
    /// `chunks` must tile `0..members.len()` in order.
    pub fn build_from(
        dir: &Path,
        sidecars: &[(FileId, PathBuf)],
        members: &[Vec<Member>],
        chunks: &[(usize, usize)],
        contigs: &Interner,
        flush_bytes: usize,
    ) -> Result<Store, String> {
        let _ = std::fs::remove_dir_all(dir);
        std::fs::create_dir_all(dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
        let mut route = Route::default();
        for (c, &(lo, hi)) in chunks.iter().enumerate() {
            for ms in &members[lo..hi] {
                for m in ms {
                    route.add(*m, c as u32);
                }
            }
        }
        let n_chunks = chunks.len();
        let routed: Vec<(ScratchFile, u64)> = sidecars
            .par_iter()
            .map(|(f, sp)| {
                let path = dir.join(format!("rows.{f}.bin"));
                route_file(*f, sp, &path, &route, n_chunks, contigs, flush_bytes)
            })
            .collect::<Result<Vec<_>, String>>()?;
        let n_rows = routed.iter().map(|(_, n)| n).sum();
        let files: Vec<ScratchFile> = routed.into_iter().map(|(s, _)| s).collect();
        let reads_path = dir.join("reads.bin");
        let wf = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(&reads_path)
            .map_err(|e| format!("cannot create {}: {e}", reads_path.display()))?;
        println!(
            "evidence shards: {n_rows} rows of {} member loci routed into {n_chunks} chunk(s) under {}",
            route.map.len(),
            dir.display()
        );
        Ok(Store {
            dir: dir.to_path_buf(),
            have: sidecars.iter().map(|(f, _)| *f).collect(),
            n_rows,
            route,
            n_chunks,
            files,
            headers: OnceLock::new(),
            reads: Reads {
                path: reads_path,
                writer: Mutex::new((Some(wf), 0)),
                reader: OnceLock::new(),
                pending: (0..n_chunks).map(|_| Mutex::new(Vec::new())).collect(),
                base: (0..n_chunks).map(|_| OnceLock::new()).collect(),
            },
            lru: Mutex::new(VecDeque::new()),
            member_ix: None,
        })
    }

    /// basenames with an existing sidecar (python `have`), in input order
    pub fn have(&self) -> &[FileId] {
        &self.have
    }

    /// number of routed rows (log only)
    pub fn n_rows(&self) -> u64 {
        self.n_rows
    }

    pub fn n_chunks(&self) -> usize {
        self.n_chunks
    }

    /// Raw routed lines of chunk `c` (without "\n"), sidecar by sidecar in input order, each
    /// file's lines in sidecar order. `f(index into have(), line)`.
    pub fn for_each_chunk_line(&self, c: usize, mut f: impl FnMut(usize, &[u8]) -> Result<(), String>) -> Result<(), String> {
        let mut comp = Vec::new();
        let mut raw = Vec::new();
        for (fi, sf) in self.files.iter().enumerate() {
            let blocks = &sf.blocks[c];
            if blocks.is_empty() {
                continue;
            }
            let fh = File::open(&sf.path).map_err(|e| format!("cannot open {}: {e}", sf.path.display()))?;
            for &(off, len) in blocks {
                comp.resize(len as usize, 0);
                fh.read_exact_at(&mut comp, off).map_err(|e| format!("read {}: {e}", sf.path.display()))?;
                raw.clear();
                DeflateDecoder::new(&comp[..])
                    .read_to_end(&mut raw)
                    .map_err(|e| format!("inflate {}: {e}", sf.path.display()))?;
                // every routed line is non-empty (it holds a locus) and "\n"-terminated
                for line in raw.split(|&b| b == b'\n').filter(|l| !l.is_empty()) {
                    f(fi, line)?;
                }
            }
        }
        Ok(())
    }

    fn headers(&self) -> Result<&Vec<SidecarHeader>, String> {
        self.headers
            .get_or_init(|| self.files.iter().map(|sf| SidecarHeader::parse(&sf.header)).collect())
            .as_ref()
            .map_err(|e| e.clone())
    }

    /// Rows of chunk `chunk`, parsed (python `_load_rows`: rows whose field count differs from
    /// the header's are skipped by `EvidenceRow::parse`). NB `contigs` added to the skeleton
    /// signature: rows intern `ref`/`mref` into the run's interner (`Ctx::contigs`).
    pub fn chunk_rows(&self, chunk: usize, contigs: &Interner) -> Result<ChunkRows, String> {
        let headers = self.headers()?;
        let mut out = ChunkRows::default();
        self.for_each_chunk_line(chunk, |fi, line| {
            let s = std::str::from_utf8(line).map_err(|e| format!("sidecar row not UTF-8: {e}"))?;
            if let Some(r) = EvidenceRow::parse(s, &headers[fi], self.files[fi].file, contigs)? {
                out.push(r);
            }
            Ok(())
        })?;
        Ok(out)
    }

    /// parent-side lookup (absorb_one_sided re-evaluation): the rows of the first chunk listing
    /// `m` (every chunk listing m holds all of m's rows); 4-chunk LRU of parsed chunks like
    /// python `_rows_cache`. Panics on scratch I/O failure (the signature has no error path).
    /// NB `contigs` added to the skeleton signature (see `chunk_rows`).
    pub fn lookup(&self, m: &Member, side: Side, contigs: &Interner) -> Vec<EvidenceRow> {
        let Some(&c) = self.route.chunks(m).first() else {
            return Vec::new();
        };
        if let Some(ix) = &self.member_ix {
            if let Some(&(off, len, fi)) = ix.map.get(m) {
                return self.lookup_indexed(ix, off, len, fi as usize, m, side, contigs).unwrap_or_else(|e| panic!("evidence store lookup: {e}"));
            }
        }
        let c = c as usize;
        let rows = {
            let mut lru = self.lru.lock().unwrap();
            if let Some(p) = lru.iter().position(|(k, _)| *k == c) {
                lru[p].1.clone()
            } else {
                let r = Arc::new(self.chunk_rows(c, contigs).unwrap_or_else(|e| panic!("evidence store lookup: {e}")));
                if lru.len() >= LRU_CHUNKS {
                    lru.pop_front();
                }
                lru.push_back((c, r.clone()));
                r
            }
        };
        rows.get(m, side).to_vec()
    }

    /// Index the rows of `wanted` members for `lookup` (absorb_one_sided): one parallel pass
    /// over all chunks (on the current rayon pool) writes each wanted member's lines -- those
    /// of the first chunk listing it, exactly what the chunk LRU path would parse -- as one
    /// compressed blob to `lookup.bin`. A lookup then reads and parses only that member's
    /// rows instead of a whole chunk (the LRU thrashes: absorb visits members in name order,
    /// not chunk order). Rows and their order are those of `chunk_rows(c).get(m, side)`.
    pub fn prefetch_members(&mut self, wanted: &rustc_hash::FxHashSet<Member>, contigs: &Interner) -> Result<(), String> {
        let headers = self.headers()?.clone();
        let path = self.dir.join("lookup.bin");
        let wf = OpenOptions::new()
            .create(true)
            .truncate(true)
            .read(true)
            .write(true)
            .open(&path)
            .map_err(|e| format!("cannot create {}: {e}", path.display()))?;
        let w = Mutex::new((wf, 0u64));
        let this = &*self;
        let parts: Vec<Vec<(Member, u64, u32, u32)>> = (0..this.n_chunks)
            .into_par_iter()
            .map(|c| -> Result<Vec<(Member, u64, u32, u32)>, String> {
                let mut per: FxHashMap<Member, (u32, Vec<u8>)> = FxHashMap::default();
                let mut order: Vec<Member> = Vec::new();
                this.for_each_chunk_line(c, |fi, line| {
                    let s = std::str::from_utf8(line).map_err(|e| format!("sidecar row not UTF-8: {e}"))?;
                    if let Some(r) = EvidenceRow::parse(s, &headers[fi], this.files[fi].file, contigs)? {
                        let m = (r.file, r.locus);
                        if wanted.contains(&m) && this.route.chunks(&m).first() == Some(&(c as u32)) {
                            let e = per.entry(m).or_insert_with(|| {
                                order.push(m);
                                (fi as u32, Vec::new())
                            });
                            e.1.extend_from_slice(line);
                            e.1.push(b'\n');
                        }
                    }
                    Ok(())
                })?;
                let mut blobs = Vec::with_capacity(order.len());
                for m in order {
                    let (fi, raw) = per.remove(&m).unwrap();
                    let mut enc = DeflateEncoder::new(Vec::new(), Compression::fast());
                    enc.write_all(&raw).map_err(|e| e.to_string())?;
                    blobs.push((m, fi, enc.finish().map_err(|e| e.to_string())?));
                }
                let mut g = w.lock().unwrap();
                let mut out = Vec::with_capacity(blobs.len());
                for (m, fi, b) in blobs {
                    g.0.write_all(&b).map_err(|e| format!("write {}: {e}", path.display()))?;
                    out.push((m, g.1, b.len() as u32, fi));
                    g.1 += b.len() as u64;
                }
                Ok(out)
            })
            .collect::<Result<_, String>>()?;
        let (mut fh, _) = w.into_inner().unwrap();
        fh.flush().map_err(|e| format!("write {}: {e}", path.display()))?;
        let mut map: FxHashMap<Member, (u64, u32, u32)> = wanted.iter().map(|m| (*m, (0, 0, 0))).collect();
        for (m, off, len, fi) in parts.into_iter().flatten() {
            map.insert(m, (off, len, fi));
        }
        self.member_ix = Some(MemberIx { fh, map });
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn lookup_indexed(&self, ix: &MemberIx, off: u64, len: u32, fi: usize, m: &Member, side: Side, contigs: &Interner) -> Result<Vec<EvidenceRow>, String> {
        let mut out = Vec::new();
        if len == 0 {
            return Ok(out);
        }
        let headers = self.headers()?;
        let mut comp = vec![0u8; len as usize];
        ix.fh.read_exact_at(&mut comp, off).map_err(|e| format!("read lookup.bin: {e}"))?;
        let mut raw = Vec::new();
        DeflateDecoder::new(&comp[..]).read_to_end(&mut raw).map_err(|e| format!("inflate lookup.bin: {e}"))?;
        for line in raw.split(|&b| b == b'\n').filter(|l| !l.is_empty()) {
            let s = std::str::from_utf8(line).map_err(|e| format!("sidecar row not UTF-8: {e}"))?;
            if let Some(r) = EvidenceRow::parse(s, &headers[fi], self.files[fi].file, contigs)? {
                if (r.file, r.locus) == *m && r.side == side {
                    out.push(r);
                }
            }
        }
        Ok(out)
    }

    /// store a record's rendered reads text; thread-safe
    pub fn put_reads(&self, chunk: usize, text: &str) -> ReadsRef {
        let mut buf = self.reads.pending[chunk].lock().unwrap();
        assert!(self.reads.base[chunk].get().is_none(), "put_reads after finish_chunk({chunk})");
        let offset = buf.len() as u64;
        buf.extend_from_slice(text.as_bytes());
        ReadsRef { chunk: chunk as u32, offset, len: text.len() as u32 }
    }

    /// Move chunk `chunk`'s pending reads text to the reads scratch file (one contiguous run)
    /// and free it. Call once, after the chunk's last `put_reads` (added to the skeleton API).
    pub fn finish_chunk(&self, chunk: usize) -> Result<(), String> {
        let buf = std::mem::take(&mut *self.reads.pending[chunk].lock().unwrap());
        let mut w = self.reads.writer.lock().unwrap();
        let base = w.1;
        let fh = w.0.as_mut().ok_or("evidence store already cleaned up")?;
        fh.write_all(&buf).map_err(|e| format!("write {}: {e}", self.reads.path.display()))?;
        w.1 += buf.len() as u64;
        self.reads.base[chunk].set(base).map_err(|_| format!("finish_chunk({chunk}) called twice"))?;
        Ok(())
    }

    pub fn read_reads(&self, r: ReadsRef) -> String {
        let c = r.chunk as usize;
        let bytes = match self.reads.base[c].get() {
            Some(&base) => {
                let fh = self.reads.reader.get_or_init(|| {
                    File::open(&self.reads.path).unwrap_or_else(|e| panic!("open {}: {e}", self.reads.path.display()))
                });
                let mut b = vec![0u8; r.len as usize];
                fh.read_exact_at(&mut b, base + r.offset)
                    .unwrap_or_else(|e| panic!("read {}: {e}", self.reads.path.display()));
                b
            }
            None => {
                let buf = self.reads.pending[c].lock().unwrap();
                buf[r.offset as usize..r.offset as usize + r.len as usize].to_vec()
            }
        };
        String::from_utf8(bytes).expect("reads text is UTF-8")
    }

    pub fn cleanup(&mut self) {
        self.lru.lock().unwrap().clear();
        self.member_ix = None;
        self.reads.writer.lock().unwrap().0 = None;
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// Exact-string locus match: parse `contig:start-end` (contig = before the LAST ':') accepting
/// only canonical tokens -- the member id is `LocusKey::name`, so a non-canonical token string
/// could never equal it in python. Contigs never interned cannot be members.
fn parse_locus_exact(s: &[u8], cache: &mut FxHashMap<Vec<u8>, Option<ContigId>>, contigs: &Interner) -> Option<LocusKey> {
    let colon = s.iter().rposition(|&b| b == b':')?;
    let (c, rest) = (&s[..colon], &s[colon + 1..]);
    let dash = rest.iter().position(|&b| b == b'-')?;
    let start = canonical_tok(&rest[..dash])?;
    let end = canonical_tok(&rest[dash + 1..])?;
    let contig = match cache.get(c) {
        Some(v) => (*v)?,
        None => {
            let v = std::str::from_utf8(c).ok().and_then(|s| contigs.get(s));
            cache.insert(c.to_vec(), v);
            v?
        }
    };
    Some(LocusKey { contig, start, end })
}

fn canonical_tok(s: &[u8]) -> Option<Tok> {
    let t = Tok::parse(s)?;
    let pre = match t.kind {
        TokKind::Pos => 0,
        TokKind::PolyA => 6,
        TokKind::Disc => 5,
        TokKind::OneSide => 8,
    };
    let num = &s[pre..];
    let digits = num.strip_prefix(b"-").unwrap_or(num);
    let ok = !digits.is_empty()
        && digits.iter().all(u8::is_ascii_digit)
        && (digits.len() == 1 || digits[0] != b'0')
        && !(num[0] == b'-' && digits == b"0");
    ok.then_some(t)
}

/// One routing task: stream sidecar `sp`, keep rows of routed members, write blocks to `out`.
fn route_file(
    file: FileId,
    sp: &Path,
    out: &Path,
    route: &Route,
    n_chunks: usize,
    contigs: &Interner,
    flush_bytes: usize,
) -> Result<(ScratchFile, u64), String> {
    let fh = File::open(sp).map_err(|e| format!("cannot open {}: {e}", sp.display()))?;
    let mut rd = BufReader::with_capacity(1 << 18, flate2::read::MultiGzDecoder::new(BufReader::with_capacity(1 << 18, fh)));
    let mut line = Vec::with_capacity(1024);
    rd.read_until(b'\n', &mut line).map_err(|e| format!("read {}: {e}", sp.display()))?;
    let header = String::from_utf8_lossy(line.strip_suffix(b"\n").unwrap_or(&line)).into_owned();
    // python `header.index("locus")` (ValueError when absent)
    let li = header
        .split('\t')
        .position(|c| c == "locus")
        .ok_or_else(|| format!("{}: sidecar header has no 'locus' column", sp.display()))?;
    let mut w = File::create(out).map_err(|e| format!("cannot create {}: {e}", out.display()))?;
    let mut bufs: Vec<Vec<u8>> = vec![Vec::new(); n_chunks];
    let mut blocks: Vec<Vec<(u64, u32)>> = vec![Vec::new(); n_chunks];
    let mut buffered = 0usize;
    let mut off = 0u64;
    let mut n_rows = 0u64;
    let mut cache = FxHashMap::default();
    let flush = |bufs: &mut Vec<Vec<u8>>, blocks: &mut Vec<Vec<(u64, u32)>>, off: &mut u64, w: &mut File| -> Result<(), String> {
        for (c, b) in bufs.iter_mut().enumerate() {
            if b.is_empty() {
                continue;
            }
            let mut enc = DeflateEncoder::new(Vec::with_capacity(b.len() / 3 + 64), Compression::fast());
            enc.write_all(b).map_err(|e| e.to_string())?;
            let z = enc.finish().map_err(|e| e.to_string())?;
            w.write_all(&z).map_err(|e| format!("write {}: {e}", out.display()))?;
            blocks[c].push((*off, z.len() as u32));
            *off += z.len() as u64;
            b.clear();
        }
        Ok(())
    };
    loop {
        line.clear();
        let n = rd.read_until(b'\n', &mut line).map_err(|e| format!("read {}: {e}", sp.display()))?;
        if n == 0 {
            break;
        }
        // python `p = line.split("\t", li + 1)`; `p[li]` keeps a trailing "\n" when last
        let Some(field) = line.split(|&b| b == b'\t').nth(li) else {
            continue;
        };
        let Some(key) = parse_locus_exact(field, &mut cache, contigs) else {
            continue;
        };
        let cs = route.chunks(&(file, key));
        if cs.is_empty() {
            continue;
        }
        n_rows += 1;
        let body = line.strip_suffix(b"\n").unwrap_or(&line);
        for &c in cs {
            let b = &mut bufs[c as usize];
            b.extend_from_slice(body);
            b.push(b'\n');
            buffered += body.len() + 1;
        }
        if buffered >= flush_bytes {
            flush(&mut bufs, &mut blocks, &mut off, &mut w)?;
            buffered = 0;
        }
    }
    flush(&mut bufs, &mut blocks, &mut off, &mut w)?;
    Ok((ScratchFile { file, path: out.to_path_buf(), header, blocks }, n_rows))
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;

    fn tmpdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("peartree-p5-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn write_gz(path: &Path, text: &str) {
        let mut e = GzEncoder::new(File::create(path).unwrap(), Compression::default());
        e.write_all(text.as_bytes()).unwrap();
        e.finish().unwrap();
    }

    fn route_of(members: &[Vec<Member>], chunks: &[(usize, usize)], m: Member) -> Vec<u32> {
        let mut r = Route::default();
        for (c, &(lo, hi)) in chunks.iter().enumerate() {
            for ms in &members[lo..hi] {
                for m in ms {
                    r.add(*m, c as u32);
                }
            }
        }
        r.chunks(&m).to_vec()
    }

    #[test]
    fn chunking_matches_python() {
        // python: n_chunks = max(1, min(256, max(4*threads, n//2000))); size = ceil(n/n_chunks)
        assert_eq!(chunk_ranges(0, 4), vec![(0, 0)]);
        assert_eq!(chunk_ranges(3, 4), vec![(0, 1), (1, 2), (2, 3)]);
        assert_eq!(chunk_ranges(10, 1), vec![(0, 3), (3, 6), (6, 9), (9, 10)]);
        let c = chunk_ranges(1_000_000, 4);
        assert_eq!(c.len(), 256);
        assert_eq!(c.last().unwrap().1, 1_000_000);
        let c = chunk_ranges(300_000, 8);
        assert_eq!(c.len(), 150);
        assert!(c.windows(2).all(|w| w[0].1 == w[1].0));
    }

    #[test]
    fn bounded_chunking_splits_large_chunks() {
        let m = |n: usize| vec![(0u32, LocusKey::parse("chr1:1-2", &Interner::new()).unwrap()); n];
        // 10 insertions, threads 1 -> python chunks [0,3) [3,6) [6,9) [9,10)
        let members: Vec<Vec<Member>> = vec![m(2), m(2), m(2), m(5), m(1), m(1), m(1), m(1), m(1), m(9)];
        assert_eq!(chunk_ranges_bounded(&members, 1, 100), chunk_ranges(10, 1));
        assert_eq!(chunk_ranges_bounded(&members, 1, 4), vec![(0, 2), (2, 3), (3, 4), (4, 6), (6, 9), (9, 10)]);
    }

    #[test]
    fn canonical_tokens() {
        for ok in [&b"120"[..], b"0", b"oneside_5", b"polyA_77", b"disc_9", b"-4"] {
            assert!(canonical_tok(ok).is_some(), "{}", String::from_utf8_lossy(ok));
        }
        for bad in [&b"0120"[..], b"+120", b"-0", b"120\n", b"oneside_05", b"", b"x"] {
            assert!(canonical_tok(bad).is_none(), "{}", String::from_utf8_lossy(bad));
        }
    }

    /// Synthetic sidecars -> routed lines identical for any thread count / flush size, in
    /// sidecar order per file, routed to every chunk listing the member; unwanted rows,
    /// non-canonical loci and other files' loci dropped.
    #[test]
    fn routing_is_deterministic_and_ordered() {
        let d = tmpdir("route");
        let contigs = Interner::new();
        let key = |s: &str| LocusKey::parse(s, &contigs).unwrap();
        let m0: Member = (0, key("chr1:100-120"));
        let m1: Member = (0, key("chr1:500-oneside_500"));
        let m2: Member = (1, key("chr2:7-9"));
        let m3: Member = (1, key("chr1:100-120"));
        let mut a = String::from("locus\tside\trole\tfrag\tseq\n");
        let mut b = String::from("frag\tlocus\tside\n");
        for i in 0..400 {
            let l = ["chr1:100-120", "chr1:500-oneside_500", "chr2:7-9", "chr1:0100-120", "chrX:1-2"][i % 5];
            a.push_str(&format!("{l}\t{}\tCLIP\tf{i}\tACGT\n", if i % 2 == 0 { "LEFT" } else { "RIGHT" }));
            b.push_str(&format!("g{i}\t{l}\tRIGHT\n"));
        }
        b.push_str("short_row\n"); // fewer fields than li+1 -> skipped
        let pa = d.join("A.txt.gz.evidence.tsv.gz");
        let pb = d.join("B.evidence.tsv.gz");
        write_gz(&pa, &a);
        write_gz(&pb, &b);
        // insertion 0: m0, m3; 1: m1; 2: m0 (shared with insertion 0); 3: m2
        let members = vec![vec![m0, m3], vec![m1], vec![m0], vec![m2]];
        let sidecars = vec![(0u32, pa.clone()), (1u32, pb.clone())];
        let chunks = vec![(0, 1), (1, 2), (2, 3), (3, 4)];
        let collect = |threads: usize, flush: usize| {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            let sd = d.join(format!("shards-{threads}-{flush}"));
            let st = pool.install(|| Store::build_from(&sd, &sidecars, &members, &chunks, &contigs, flush)).unwrap();
            let mut per_chunk = Vec::new();
            for c in 0..chunks.len() {
                let mut v = Vec::new();
                st.for_each_chunk_line(c, |fi, l| {
                    v.push((fi, String::from_utf8(l.to_vec()).unwrap()));
                    Ok(())
                })
                .unwrap();
                per_chunk.push(v);
            }
            (st.n_rows(), per_chunk)
        };
        let base = collect(1, FLUSH_BYTES);
        for (t, f) in [(1, 1), (4, 37), (8, 1000), (3, FLUSH_BYTES)] {
            assert_eq!(collect(t, f), base, "threads={t} flush={f}");
        }
        let (n_rows, per_chunk) = base;
        // A: 80 rows m0 + 80 rows m1; B: 80 rows m2 + 80 rows m3
        assert_eq!(n_rows, 320);
        let c0 = &per_chunk[0];
        assert_eq!(c0.len(), 160);
        assert!(c0[..80].iter().all(|(fi, l)| *fi == 0 && l.starts_with("chr1:100-120\t")));
        assert!(c0[80..].iter().all(|(fi, l)| *fi == 1 && l.contains("\tchr1:100-120\t")));
        let frags: Vec<usize> = c0[..80].iter().map(|(_, l)| l.split('\t').nth(3).unwrap()[1..].parse().unwrap()).collect();
        assert!(frags.windows(2).all(|w| w[0] < w[1]), "sidecar order kept");
        assert_eq!(per_chunk[1].len(), 80);
        assert!(per_chunk[1].iter().all(|(_, l)| l.starts_with("chr1:500-oneside_500\t")));
        assert_eq!(per_chunk[2], c0[..80].to_vec(), "shared member routed to both chunks");
        assert_eq!(per_chunk[3].len(), 80);
        assert_eq!(route_of(&members, &chunks, m0), vec![0, 2]);
        assert_eq!(route_of(&members, &chunks, m2), vec![3]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn locus_last_column_matches_only_unterminated_line_like_python() {
        let d = tmpdir("lastcol");
        let contigs = Interner::new();
        let k = LocusKey::parse("chr1:100-120", &contigs).unwrap();
        let p = d.join("A.evidence.tsv.gz");
        write_gz(&p, "side\tlocus\nLEFT\tchr1:100-120\nRIGHT\tchr1:100-120");
        let st = Store::build_from(&d.join("s"), &[(0, p)], &[vec![(0, k)]], &[(0, 1)], &contigs, 10).unwrap();
        assert_eq!(st.n_rows(), 1);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_locus_column_is_an_error() {
        let d = tmpdir("nolocus");
        let contigs = Interner::new();
        let p = d.join("A.evidence.tsv.gz");
        write_gz(&p, "side\trole\nLEFT\tCLIP\n");
        assert!(Store::build_from(&d.join("s"), &[(0, p)], &[], &[(0, 0)], &contigs, 10).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn reads_roundtrip_pending_and_finished() {
        let d = tmpdir("reads");
        let contigs = Interner::new();
        let mut st = Store::build_from(&d.join("s"), &[], &[vec![], vec![]], &[(0, 1), (1, 2)], &contigs, 10).unwrap();
        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        let refs: Vec<Vec<(ReadsRef, String)>> = pool.install(|| {
            (0..2usize)
                .into_par_iter()
                .map(|c| {
                    (0..50)
                        .map(|i| {
                            let t = format!(">x{c}_{i}|LEFT|CLIP|s|f|1\n{}\n", "ACGT".repeat(i));
                            (st.put_reads(c, &t), t)
                        })
                        .collect::<Vec<_>>()
                })
                .collect()
        });
        for (r, t) in &refs[1] {
            assert_eq!(&st.read_reads(*r), t, "pending chunk served from memory");
        }
        st.finish_chunk(1).unwrap();
        st.finish_chunk(0).unwrap();
        for v in &refs {
            for (r, t) in v {
                assert_eq!(&st.read_reads(*r), t);
            }
        }
        assert!(st.finish_chunk(0).is_err());
        st.cleanup();
        assert!(!d.join("s").exists());
        let _ = std::fs::remove_dir_all(&d);
    }

    /// chunk_rows / lookup parse rows with P3's `EvidenceRow::parse`.
    #[test]
    fn chunk_rows_and_lookup_group_by_member_side() {
        let d = tmpdir("rows");
        let contigs = Interner::new();
        let k = LocusKey::parse("chr1:100-120", &contigs).unwrap();
        let p = d.join("A.evidence.tsv.gz");
        write_gz(&p, "locus\tside\trole\tfrag\nchr1:100-120\tLEFT\tCLIP\ta\nchr1:100-120\tRIGHT\tCLIP\tb\nchr1:100-120\tLEFT\tMATE\tc\n");
        let st = Store::build_from(&d.join("s"), &[(0, p)], &[vec![(0, k)]], &[(0, 1)], &contigs, 10).unwrap();
        let cr = st.chunk_rows(0, &contigs).unwrap();
        let l = cr.get(&(0, k), Side::Left);
        assert_eq!(l.len(), 2);
        assert_eq!(&*l[0].frag, "a");
        assert_eq!(&*l[1].frag, "c");
        assert_eq!(cr.get(&(0, k), Side::Right).len(), 1);
        assert_eq!(st.lookup(&(0, k), Side::Left, &contigs).len(), 2);
        let _ = std::fs::remove_dir_all(&d);
    }
}
