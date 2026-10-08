//! Memory-bounded access to the pooled reads FASTAs. FOUNDATION (implemented, tested).
//!
//! Replaces `inputs.read_reads_fa` + `RteAnnotator.load_evidence`, which held every read of every
//! insertion as Python objects (PD49229, 722 colonies: killed at 32 GB).
//!
//! # Why an index + spill file, not a grouped stream
//! combine (python `write_evidence_outputs`, rust `evidence/output.rs`) writes reads per name in
//! `names` order, but the files are NOT grouped by locus in practice: PD51635's
//! `insertions.reads.fa.gz` has 5,337 of 116,637 loci in 2-17 separate runs (one-sided
//! `...-oneside_N` loci re-emitted with identical reads up to millions of records apart; python
//! appends every copy, so the duplicates are part of the reference behaviour). The genotype2
//! extra-reads file (`tools/genotype_extra_reads.py merge`) is written colony by colony, so a
//! locus is scattered over up to one run per colony by construction.
//!
//! # Design
//! [`ReadStore::build`] makes ONE pass over the (gzipped) FASTA. Every record of a wanted locus
//! is appended, compactly encoded (4-bit bases, see `packed.rs`), to an anonymous spill file
//! (created in `tmp_dir`, unlinked immediately: it vanishes when the process exits, however it
//! exits); the in-memory index maps locus -> byte runs in that file (adjacent records of one
//! locus merge into one run) + record count. [`ReadStore::reads`] then `pread`s one locus's runs
//! in file order. RAM = the index (~100 B per locus) + whatever loci are in flight; disk = about
//! 0.9x the uncompressed sequence bytes (PD51635: 1.39 GB spill for 1.57 Gbases; peak RSS of a
//! full index + read-back pass 120 MB). Reads of one locus come back in file order,
//! duplicates included -- exactly the list python's `read_reads_fa` built.
//!
//! [`cap_reads`] is `RteAnnotator._cap_reads` (rte_max_reads, default 400), applied to the
//! complete per-locus list exactly where python applies it (in `annotate`, after pooling the
//! genotype reads), so the kept reads and their order are identical.
//!
//! [`LocusLoader`] bundles the evidence table + the two stores and yields one locus's
//! [`LocusData`]; [`map_loci`] runs a closure over loci with rayon in bounded chunks, results in
//! input order (deterministic regardless of scheduling).

use crate::inputs::{
    open_text, parse_read_header, read_evidence_tsv, EvidenceRead, EvidenceTable, InsertionEvidence, GT_ROLE_PREFIX,
};
use crate::packed::PackedSeq;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static SPILL_SEQ: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Debug, Default)]
struct LocusRuns {
    /// (byte offset, byte length) in the spill file, in file order
    runs: Vec<(u64, u32)>,
    n: u32,
}

/// Build-time statistics (also printed by `peartree-rte scan`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StoreStats {
    pub records_seen: u64,
    pub records_kept: u64,
    pub loci: u64,
    pub runs: u64,
    pub multi_run_loci: u64,
    pub max_reads_per_locus: u64,
    pub spill_bytes: u64,
}

/// Index + spill file of one reads FASTA.
pub struct ReadStore {
    file: Option<File>,
    index: FxHashMap<Box<str>, LocusRuns>,
    pub stats: StoreStats,
}

/// Options of [`ReadStore::build`].
#[derive(Clone, Debug)]
pub struct StoreOptions {
    /// keep only records whose (upper-cased) role starts with this (genotype reads: "GT_")
    pub role_prefix: Option<String>,
    /// directory of the (immediately unlinked) spill file
    pub tmp_dir: PathBuf,
}

impl ReadStore {
    /// An empty store (file absent: python skips a missing sidecar).
    pub fn empty() -> ReadStore {
        ReadStore { file: None, index: FxHashMap::default(), stats: StoreStats::default() }
    }

    /// One pass over `path` (plain or gzip): spill the records of `wanted` loci (all when None).
    /// Parsing is python `read_reads_fa`: a '>' line starts a record named by its first
    /// whitespace-delimited token; other non-empty lines are stripped and concatenated.
    pub fn build(path: &Path, wanted: Option<&FxHashSet<String>>, opts: &StoreOptions) -> io::Result<ReadStore> {
        let mut rdr = open_text(path)?;
        std::fs::create_dir_all(&opts.tmp_dir)?;
        let spill_path = opts.tmp_dir.join(format!(
            ".peartree-rte.{}.{}.spill",
            std::process::id(),
            SPILL_SEQ.fetch_add(1, Ordering::Relaxed)
        ));
        let file = File::options().read(true).write(true).create_new(true).open(&spill_path)?;
        // unlink now: the open handle keeps the data; nothing is left behind on any exit
        std::fs::remove_file(&spill_path)?;
        let mut w = BufWriter::with_capacity(1 << 20, file.try_clone()?);
        let mut store = ReadStore { file: None, index: FxHashMap::default(), stats: StoreStats::default() };
        let mut off: u64 = 0;
        let mut name: Option<String> = None;
        let mut seq: Vec<u8> = Vec::new();
        let mut buf: Vec<u8> = Vec::new();
        let mut rec: Vec<u8> = Vec::new();
        let mut line: Vec<u8> = Vec::new();
        let flush = |name: &Option<String>,
                         seq: &[u8],
                         store: &mut ReadStore,
                         w: &mut BufWriter<File>,
                         off: &mut u64,
                         rec: &mut Vec<u8>|
         -> io::Result<()> {
            let Some(name) = name else { return Ok(()) };
            store.stats.records_seen += 1;
            let h = parse_read_header(name);
            if let Some(wt) = wanted {
                if !wt.contains(&h.insertion_id) {
                    return Ok(());
                }
            }
            if let Some(p) = &opts.role_prefix {
                if !h.role.starts_with(p.as_str()) {
                    return Ok(());
                }
            }
            rec.clear();
            encode_record(rec, &h.side, &h.role, &h.sample, &h.frag, &h.r12, seq)?;
            w.write_all(rec)?;
            let len = rec.len() as u64;
            let e = store.index.entry(h.insertion_id.into_boxed_str()).or_default();
            match e.runs.last_mut() {
                Some(last) if last.0 + last.1 as u64 == *off && (last.1 as u64 + len) <= u32::MAX as u64 => {
                    last.1 += len as u32
                }
                _ => e.runs.push((*off, len as u32)),
            }
            e.n += 1;
            *off += len;
            store.stats.records_kept += 1;
            Ok(())
        };
        loop {
            line.clear();
            if rdr.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            // python text mode: universal newlines, then rstrip("\n")
            while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
                line.pop();
            }
            if line.first() == Some(&b'>') {
                flush(&name, &seq, &mut store, &mut w, &mut off, &mut rec)?;
                let rest = String::from_utf8_lossy(&line[1..]);
                name = Some(rest.split_whitespace().next().unwrap_or("").to_string());
                seq.clear();
            } else if !line.is_empty() {
                buf.clear();
                buf.extend_from_slice(line.trim_ascii());
                seq.extend_from_slice(&buf);
            }
        }
        flush(&name, &seq, &mut store, &mut w, &mut off, &mut rec)?;
        w.flush()?;
        drop(w);
        store.stats.spill_bytes = off;
        store.stats.loci = store.index.len() as u64;
        for v in store.index.values() {
            store.stats.runs += v.runs.len() as u64;
            if v.runs.len() > 1 {
                store.stats.multi_run_loci += 1;
            }
            store.stats.max_reads_per_locus = store.stats.max_reads_per_locus.max(v.n as u64);
        }
        store.file = Some(file);
        Ok(store)
    }

    /// Number of stored reads of `locus` (0 when absent).
    pub fn n_reads(&self, locus: &str) -> usize {
        self.index.get(locus).map_or(0, |v| v.n as usize)
    }

    pub fn contains(&self, locus: &str) -> bool {
        self.index.contains_key(locus)
    }

    /// Every stored locus (hash order -- sort before using it for output).
    pub fn loci(&self) -> Vec<String> {
        self.index.keys().map(|k| k.to_string()).collect()
    }

    /// Every stored read of `locus`, in file order (duplicates included).
    pub fn reads(&self, locus: &str) -> io::Result<Vec<EvidenceRead>> {
        let (Some(file), Some(e)) = (&self.file, self.index.get(locus)) else { return Ok(Vec::new()) };
        let mut out = Vec::with_capacity(e.n as usize);
        let mut buf = Vec::new();
        for &(off, len) in &e.runs {
            buf.resize(len as usize, 0);
            file.read_exact_at(&mut buf, off)?;
            let mut p = 0;
            while p < buf.len() {
                let (r, used) = decode_record(&buf[p..]).ok_or_else(|| io::Error::other("corrupt spill record"))?;
                out.push(r);
                p += used;
            }
        }
        Ok(out)
    }
}

fn put_str(rec: &mut Vec<u8>, s: &str) -> io::Result<()> {
    let n = u16::try_from(s.len()).map_err(|_| io::Error::other(format!("header field too long: {}", s.len())))?;
    rec.extend_from_slice(&n.to_le_bytes());
    rec.extend_from_slice(s.as_bytes());
    Ok(())
}

fn encode_record(rec: &mut Vec<u8>, side: &str, role: &str, sample: &str, frag: &str, r12: &str, seq: &[u8]) -> io::Result<()> {
    rec.extend_from_slice(&[0u8; 4]); // total length, patched below
    for s in [side, role, sample, frag, r12] {
        put_str(rec, s)?;
    }
    PackedSeq::encode(seq).write_to(rec);
    let n = u32::try_from(rec.len()).map_err(|_| io::Error::other("read record too long"))?;
    rec[..4].copy_from_slice(&n.to_le_bytes());
    Ok(())
}

fn decode_record(buf: &[u8]) -> Option<(EvidenceRead, usize)> {
    let total = u32::from_le_bytes(buf.get(..4)?.try_into().ok()?) as usize;
    let b = buf.get(..total)?;
    let mut p = 4;
    let mut strs: [Box<str>; 5] = Default::default();
    for s in strs.iter_mut() {
        let n = u16::from_le_bytes(b.get(p..p + 2)?.try_into().ok()?) as usize;
        *s = std::str::from_utf8(b.get(p + 2..p + 2 + n)?).ok()?.into();
        p += 2 + n;
    }
    let (seq, _) = PackedSeq::read_from(&b[p..])?;
    let [side, role, sample, frag, r12] = strs;
    Some((EvidenceRead { side, role, sample, frag, r12, seq }, total))
}

/// `RteAnnotator._cap_reads(reads, cap)`: keep every junction read (CLIP / POLYA / SPAN, the
/// first `cap` of them), then thin the rest evenly (`rest[int(i * len(rest) / room)]`). NOTE the
/// order: junction reads first, then the sampled rest (python builds `keep` that way).
pub fn cap_reads(reads: Vec<EvidenceRead>, cap: usize) -> Vec<EvidenceRead> {
    if reads.len() <= cap {
        return reads;
    }
    let is_junction = |r: &EvidenceRead| matches!(&*r.role, "CLIP" | "POLYA" | "SPAN");
    let mut keep: Vec<EvidenceRead> = Vec::new();
    let mut rest: Vec<EvidenceRead> = Vec::new();
    for r in reads {
        if is_junction(&r) {
            if keep.len() < cap {
                keep.push(r);
            }
        } else {
            rest.push(r);
        }
    }
    let room = cap - keep.len();
    if room > 0 && !rest.is_empty() {
        let step = rest.len() as f64 / room as f64;
        let picks: Vec<usize> = (0..room).map(|i| (i as f64 * step) as usize).collect();
        for i in picks {
            keep.push(rest[i].clone());
        }
    }
    keep
}

/// What one locus contributes (python: `self.evidence.get(key)` + `self.gt_reads.get(key)`).
#[derive(Clone, Debug, Default)]
pub struct LocusData {
    /// junction rows + combine reads (empty when the locus is in neither sidecar)
    pub evidence: InsertionEvidence,
    /// GT_* reads of the genotype_reads file (kept apart; pooled by the annotator)
    pub gt_reads: Vec<EvidenceRead>,
}

/// Evidence table + combine reads store + genotype reads store.
pub struct LocusLoader {
    pub evidence: EvidenceTable,
    pub reads: ReadStore,
    pub gt: ReadStore,
    /// a genotype_reads file was loaded (python `has_gt_reads`: adds the gt_* columns)
    pub has_gt_reads: bool,
}

impl LocusLoader {
    /// `load_evidence(evidence_path, reads_path, wanted, gt_reads_path)`: each path may be None
    /// or missing on disk (skipped, like python).
    pub fn open(
        evidence: Option<&Path>,
        reads: Option<&Path>,
        gt_reads: Option<&Path>,
        wanted: &FxHashSet<String>,
        tmp_dir: &Path,
    ) -> io::Result<LocusLoader> {
        fn exists(p: Option<&Path>) -> Option<&Path> {
            p.filter(|p| p.exists())
        }
        let ev = match exists(evidence) {
            Some(p) => read_evidence_tsv(p, Some(wanted))?,
            None => EvidenceTable::default(),
        };
        let opts = StoreOptions { role_prefix: None, tmp_dir: tmp_dir.to_path_buf() };
        let rd = match exists(reads) {
            Some(p) => ReadStore::build(p, Some(wanted), &opts)?,
            None => ReadStore::empty(),
        };
        let gopts = StoreOptions { role_prefix: Some(GT_ROLE_PREFIX.to_string()), ..opts };
        let (gt, has_gt) = match exists(gt_reads) {
            Some(p) => (ReadStore::build(p, Some(wanted), &gopts)?, true),
            None => (ReadStore::empty(), false),
        };
        Ok(LocusLoader { evidence: ev, reads: rd, gt, has_gt_reads: has_gt })
    }

    pub fn load(&self, locus: &str) -> io::Result<LocusData> {
        let mut evidence = InsertionEvidence::new(locus);
        if let Some(js) = self.evidence.get(locus) {
            evidence.junctions = js.clone();
        }
        evidence.reads = self.reads.reads(locus)?;
        Ok(LocusData { evidence, gt_reads: self.gt.reads(locus)? })
    }
}

/// Run `f` over `loci` with rayon, `chunk` loci in flight at a time; results in `loci` order.
pub fn map_loci<T, R, F>(loci: &[T], chunk: usize, f: F) -> Vec<R>
where
    T: Sync,
    R: Send,
    F: Fn(&T) -> R + Sync + Send,
{
    let mut out = Vec::with_capacity(loci.len());
    for c in loci.chunks(chunk.max(1)) {
        let part: Vec<R> = c.par_iter().map(&f).collect();
        out.extend(part);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;

    fn tmpdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("rte_stream_{tag}_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn write_gz(path: &Path, text: &str) {
        let mut e = GzEncoder::new(File::create(path).unwrap(), Compression::fast());
        e.write_all(text.as_bytes()).unwrap();
        e.finish().unwrap();
    }

    fn names(rs: &[EvidenceRead]) -> Vec<String> {
        rs.iter().map(|r| format!("{}|{}|{}|{}|{}", r.side, r.role, r.sample, r.frag, r.r12)).collect()
    }

    const FA: &str = "\
>chrA:10-20|LEFT|CLIP|S1|f1|1 some comment
ACGT
acgt
>chrB:5-oneside_5|x|RIGHT|MATE|S1|f2|2
NNNN
>chrA:10-20|RIGHT|clip|S2|f3|1
GGGG

>chrC:1-2|LEFT|CLIP|S3|f9|1
TTTT
>chrA:10-20|LEFT|GT_MATE|S4|f4|2
CC*C
>chrB:5-oneside_5|x|RIGHT|MATE|S1|f2|2
NNNN
";

    #[test]
    fn split_locus_pipes_and_file_order() {
        let d = tmpdir("split");
        let p = d.join("r.fa.gz");
        write_gz(&p, FA);
        let opts = StoreOptions { role_prefix: None, tmp_dir: d.clone() };
        let s = ReadStore::build(&p, None, &opts).unwrap();
        // chrA's records are non-contiguous (interrupted by chrB / chrC): 2 runs, 3 reads
        let a = s.reads("chrA:10-20").unwrap();
        assert_eq!(names(&a), vec!["LEFT|CLIP|S1|f1|1", "RIGHT|CLIP|S2|f3|1", "LEFT|GT_MATE|S4|f4|2"]);
        assert_eq!(a[0].seq(), b"ACGTACGT");
        assert_eq!(a[2].seq(), b"CC*C");
        // a locus name containing '|' and duplicated records (python keeps both copies)
        let b = s.reads("chrB:5-oneside_5|x").unwrap();
        assert_eq!(names(&b), vec!["RIGHT|MATE|S1|f2|2", "RIGHT|MATE|S1|f2|2"]);
        assert_eq!(s.stats.records_seen, 6);
        assert_eq!(s.stats.loci, 3);
        assert_eq!(s.stats.multi_run_loci, 2);
        assert_eq!(s.reads("nope").unwrap().len(), 0);
        // no spill file left behind in tmp_dir
        assert_eq!(std::fs::read_dir(&d).unwrap().filter(|e| e.as_ref().unwrap().file_name().to_string_lossy().ends_with(".spill")).count(), 0);
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn wanted_and_gt_role_filter() {
        let d = tmpdir("gt");
        let p = d.join("g.fa");
        std::fs::write(&p, FA).unwrap();
        let w: FxHashSet<String> = ["chrA:10-20".to_string(), "chrC:1-2".to_string()].into_iter().collect();
        let opts = StoreOptions { role_prefix: Some("GT_".into()), tmp_dir: d.clone() };
        let s = ReadStore::build(&p, Some(&w), &opts).unwrap();
        assert_eq!(names(&s.reads("chrA:10-20").unwrap()), vec!["LEFT|GT_MATE|S4|f4|2"]);
        assert!(!s.contains("chrC:1-2")); // wanted but no GT_ read
        assert!(!s.contains("chrB:5-oneside_5|x")); // not wanted
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn loader_pools_like_python() {
        let d = tmpdir("loader");
        let fa = d.join("r.fa.gz");
        write_gz(&fa, FA);
        let ev = d.join("e.tsv");
        std::fs::write(&ev, "insertion_id\tside\tn_independent\nchrA:10-20\tLEFT\t2\n").unwrap();
        let w: FxHashSet<String> = ["chrA:10-20".to_string(), "chrZ".to_string()].into_iter().collect();
        let l = LocusLoader::open(Some(&ev), Some(&fa), Some(&fa), &w, &d).unwrap();
        assert!(l.has_gt_reads);
        let x = l.load("chrA:10-20").unwrap();
        assert_eq!(x.evidence.junctions.len(), 1);
        assert_eq!(x.evidence.reads.len(), 3);
        assert_eq!(names(&x.gt_reads), vec!["LEFT|GT_MATE|S4|f4|2"]);
        let z = l.load("chrZ").unwrap();
        assert!(z.evidence.reads.is_empty() && z.evidence.junctions.is_empty() && z.gt_reads.is_empty());
        let l2 = LocusLoader::open(None, None, Some(&d.join("missing.fa")), &w, &d).unwrap();
        assert!(!l2.has_gt_reads);
        std::fs::remove_dir_all(&d).ok();
    }

    fn mk(role: &str, i: usize) -> EvidenceRead {
        EvidenceRead {
            side: "LEFT".into(),
            role: role.into(),
            sample: "S".into(),
            frag: i.to_string().into(),
            r12: "1".into(),
            seq: PackedSeq::encode(b"A"),
        }
    }

    #[test]
    fn cap_reads_like_python() {
        // python: reads = [MATE 0..9, CLIP 10, MATE 11..14, POLYA 15, GT_MATE 16..19]; _cap_reads(reads, 6)
        //   -> ["10", "15", "0", "4", "9", "14"] (rest has 18, room 4, step 4.5)
        let mut reads = Vec::new();
        for i in 0..20 {
            let role = match i {
                10 => "CLIP",
                15 => "POLYA",
                16..=19 => "GT_MATE",
                _ => "MATE",
            };
            reads.push(mk(role, i));
        }
        let got: Vec<String> = cap_reads(reads.clone(), 6).iter().map(|r| r.frag.to_string()).collect();
        assert_eq!(got, vec!["10", "15", "0", "4", "9", "14"]);
        assert_eq!(cap_reads(reads.clone(), 20).len(), 20);
        // more junction reads than the cap: the first `cap` junction reads only
        let j: Vec<EvidenceRead> = (0..5).map(|i| mk("SPAN", i)).chain((5..8).map(|i| mk("MATE", i))).collect();
        let got: Vec<String> = cap_reads(j, 3).iter().map(|r| r.frag.to_string()).collect();
        assert_eq!(got, vec!["0", "1", "2"]);
    }

    #[test]
    fn map_loci_keeps_order() {
        let v: Vec<usize> = (0..1000).collect();
        let r = map_loci(&v, 7, |x| x * 2);
        assert_eq!(r, (0..1000).map(|x| x * 2).collect::<Vec<_>>());
    }
}
