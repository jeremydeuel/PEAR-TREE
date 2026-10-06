//! 2bit genome access. OWNER: P1.
//!
//! Mirrors src/combine_insertions_get_sequence.py (`get_sequence`, py2bit semantics).
//! SPEC.md §2.3.

use rustc_hash::FxHashMap;
use std::fs::File;
use std::os::unix::fs::FileExt;
use std::path::Path;

/// Reference access used by the genotyping output, the slippage / far-pair filters and the
/// SHORT-overhang check (python `ref_fetch`). Must be `Sync`: called from rayon workers.
pub trait RefFetch: Sync {
    /// python `get_sequence(seqname, start, end)` -> UPPERCASE bases, or empty on any failure.
    fn fetch(&self, seqname: &str, start: i64, end: i64) -> Vec<u8>;
}

/// A .2bit file (UCSC format v0, both endiannesses, 32-bit offsets; 64-bit offset variant
/// optional). Reads with positional reads (`FileExt::read_at`) so one handle is shared by all
/// threads; only the header/index (+ per-sequence N-block and length tables, loaded lazily or at
/// open) is held in memory. Soft-mask blocks are ignored (py2bit default `storeMasked=False`
/// returns uppercase); N blocks yield 'N'.
pub struct Genome {
    file: File,
    /// in file index order
    seqs: Vec<SeqInfo>,
    by_name: FxHashMap<String, usize>,
}

struct SeqInfo {
    name: String,
    len: u64,
    /// file offset of the first packed-DNA byte
    dna_off: u64,
    /// (start, end) half-open, sorted, non-overlapping
    n_blocks: Vec<(u64, u64)>,
}

struct Rd<'a> {
    f: &'a File,
    swap: bool,
}

impl<'a> Rd<'a> {
    fn bytes(&self, off: u64, n: usize) -> Result<Vec<u8>, String> {
        let mut b = vec![0u8; n];
        self.f.read_exact_at(&mut b, off).map_err(|e| format!("2bit read error at {off}: {e}"))?;
        Ok(b)
    }
    fn u32_at(&self, off: u64) -> Result<u32, String> {
        let b = self.bytes(off, 4)?;
        let a = [b[0], b[1], b[2], b[3]];
        Ok(if self.swap { u32::from_be_bytes(a) } else { u32::from_le_bytes(a) })
    }
    fn u64_at(&self, off: u64) -> Result<u64, String> {
        let b = self.bytes(off, 8)?;
        let mut a = [0u8; 8];
        a.copy_from_slice(&b);
        Ok(if self.swap { u64::from_be_bytes(a) } else { u64::from_le_bytes(a) })
    }
    fn u32s(&self, off: u64, n: usize) -> Result<Vec<u32>, String> {
        let b = self.bytes(off, n * 4)?;
        Ok(b.chunks_exact(4)
            .map(|c| {
                let a = [c[0], c[1], c[2], c[3]];
                if self.swap { u32::from_be_bytes(a) } else { u32::from_le_bytes(a) }
            })
            .collect())
    }
}

const SIG: u32 = 0x1A41_2743;

impl Genome {
    pub fn open(path: &Path) -> Result<Genome, String> {
        let file = File::open(path).map_err(|e| format!("cannot open 2bit {}: {e}", path.display()))?;
        let ctx = |e: String| format!("{}: {e}", path.display());
        let mut sig = [0u8; 4];
        file.read_exact_at(&mut sig, 0).map_err(|e| ctx(format!("cannot read header: {e}")))?;
        let swap = if u32::from_le_bytes(sig) == SIG {
            false
        } else if u32::from_be_bytes(sig) == SIG {
            true
        } else {
            return Err(ctx("not a 2bit file (bad signature)".into()));
        };
        let rd = Rd { f: &file, swap };
        let version = rd.u32_at(4).map_err(ctx)?;
        if version > 1 {
            return Err(ctx(format!("unsupported 2bit version {version}")));
        }
        let count = rd.u32_at(8).map_err(ctx)? as usize;
        // index
        let mut pos: u64 = 16;
        let mut index: Vec<(String, u64)> = Vec::with_capacity(count);
        for _ in 0..count {
            let nl = rd.bytes(pos, 1).map_err(ctx)?[0] as usize;
            let name = String::from_utf8_lossy(&rd.bytes(pos + 1, nl).map_err(ctx)?).into_owned();
            pos += 1 + nl as u64;
            let off = if version == 1 {
                let o = rd.u64_at(pos).map_err(ctx)?;
                pos += 8;
                o
            } else {
                let o = rd.u32_at(pos).map_err(ctx)? as u64;
                pos += 4;
                o
            };
            index.push((name, off));
        }
        let mut seqs = Vec::with_capacity(count);
        let mut by_name = FxHashMap::default();
        for (name, off) in index {
            let len = rd.u32_at(off).map_err(ctx)? as u64;
            let nb = rd.u32_at(off + 4).map_err(ctx)? as usize;
            let starts = rd.u32s(off + 8, nb).map_err(ctx)?;
            let sizes = rd.u32s(off + 8 + 4 * nb as u64, nb).map_err(ctx)?;
            let mb_off = off + 8 + 8 * nb as u64;
            let mb = rd.u32_at(mb_off).map_err(ctx)? as u64;
            // mask block table, then the reserved word, then the packed DNA
            let dna_off = mb_off + 4 + 8 * mb + 4;
            let n_blocks = starts.iter().zip(&sizes).map(|(&s, &z)| (s as u64, s as u64 + z as u64)).collect();
            by_name.entry(name.clone()).or_insert(seqs.len());
            seqs.push(SeqInfo { name, len, dna_off, n_blocks });
        }
        Ok(Genome { file, seqs, by_name })
    }

    /// Sequence names in 2bit index order.
    pub fn seqnames(&self) -> Vec<&str> {
        self.seqs.iter().map(|s| s.name.as_str()).collect()
    }

    /// Length of a sequence by exact name.
    pub fn seq_len(&self, name: &str) -> Option<u64> {
        self.by_name.get(name).map(|&i| self.seqs[i].len)
    }

    /// Resolve a seqname like `get_sequence` (index into `seqs`), None if not found.
    fn resolve(&self, seqname: &str) -> Option<usize> {
        if let Some(&i) = self.by_name.get(seqname) {
            return Some(i);
        }
        let mut name: String = seqname.to_string();
        if seqname == "MT" {
            name = "chrM".to_string();
        } else {
            let new_seqname = format!("chr{seqname}");
            if let Some(&i) = self.by_name.get(new_seqname.as_str()) {
                return Some(i);
            }
            if seqname.chars().count() > 2 {
                if let Some(stripped) = seqname.strip_suffix(".1") {
                    name = stripped.to_string();
                }
                if let Some(s) = self.seqs.iter().find(|s| s.name.contains(name.as_str())) {
                    name = s.name.clone();
                }
            }
        }
        self.by_name.get(name.as_str()).copied()
    }
}

impl RefFetch for Genome {
    /// Name mapping (get_sequence.py:37-55): if `seqname` not present: 'MT' -> 'chrM'; else if
    /// `chr{seqname}` present use it; else if len(seqname) > 2: strip a trailing '.1', then the
    /// FIRST name (index order) that contains it as a substring. Still absent -> "".
    /// `start >= end` -> "". py2bit then: `end` is clamped to the sequence length; `start < 0`
    /// or `start >= clamped end` -> error -> "" (verified against py2bit, SPEC.md §2.3).
    fn fetch(&self, seqname: &str, start: i64, end: i64) -> Vec<u8> {
        let Some(si) = self.resolve(seqname) else { return Vec::new() };
        if start >= end {
            return Vec::new();
        }
        let sq = &self.seqs[si];
        let end = (end as i128).min(sq.len as i128);
        if start < 0 || (start as i128) >= end {
            return Vec::new();
        }
        let (start, end) = (start as u64, end as u64);
        let b0 = start / 4;
        let b1 = (end - 1) / 4;
        let mut buf = vec![0u8; (b1 - b0 + 1) as usize];
        if self.file.read_exact_at(&mut buf, sq.dna_off + b0).is_err() {
            return Vec::new();
        }
        const BASES: [u8; 4] = *b"TCAG";
        let mut out = Vec::with_capacity((end - start) as usize);
        for p in start..end {
            let byte = buf[(p / 4 - b0) as usize];
            let shift = 6 - 2 * (p % 4);
            out.push(BASES[((byte >> shift) & 3) as usize]);
        }
        // N blocks (sorted): overwrite the overlap
        let first = sq.n_blocks.partition_point(|&(_, e)| e <= start);
        for &(bs, be) in &sq.n_blocks[first..] {
            if bs >= end {
                break;
            }
            let a = bs.max(start);
            let z = be.min(end);
            for p in a..z {
                out[(p - start) as usize] = b'N';
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the hg38 2bit of the equivalence-harness fixture; tests are skipped when it is absent
    fn hg38() -> Option<Genome> {
        let p = std::env::var("P1_GENOME_2BIT").unwrap_or_else(|_| {
            "/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad/genomes/hg38.2bit".into()
        });
        if !Path::new(&p).exists() {
            eprintln!("skipping: {p} not found");
            return None;
        }
        Some(Genome::open(Path::new(&p)).unwrap())
    }

    /// expected values from `combine_insertions_get_sequence.get_sequence` (py2bit)
    #[test]
    fn fetch_matches_py2bit() {
        let Some(g) = hg38() else { return };
        assert_eq!(g.fetch("chr1", 9990, 10010), b"NNNNNNNNNNTAACCCTAAC".to_vec(), "chr1:9990-10010");
        assert_eq!(g.fetch("chr1", 10000, 10040), b"TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC".to_vec(), "chr1:10000-10040");
        assert_eq!(g.fetch("chr1", 100000, 100050), b"ACTAAGCACACAGAGAATAATGTCTAGAATCTGAGTGCCATGTTATCAAA".to_vec(), "chr1:100000-100050");
        assert_eq!(g.fetch("chr22", 10510000, 10510100), b"GAATTCTTGTGTTTATATAATAAGATGTCCTATAATTTCTGTTTGGAATATAAAATCAGCAACTAATATGTATTTTCAAAGCATTATAAATACAGAGTGC".to_vec(), "chr22:10510000-10510100");
        assert_eq!(g.fetch("22", 20000000, 20000037), b"AGCCTCAGGAAGGAGGCAGTGCTGCCAGCCCTTGGGG".to_vec(), "22:20000000-20000037");
        assert_eq!(g.fetch("MT", 0, 30), b"GATCACAGGTCTATCACCCTATTAACCACT".to_vec(), "MT:0-30");
        assert_eq!(g.fetch("chrM", 16560, 16580), b"ATCACGATG".to_vec(), "chrM:16560-16580");
        assert_eq!(g.fetch("chrM", 16560, 16600), b"ATCACGATG".to_vec(), "chrM:16560-16600");
        assert_eq!(g.fetch("chrM", 16571, 16590), b"".to_vec(), "chrM:16571-16590");
        assert_eq!(g.fetch("chrM", -3, 5), b"".to_vec(), "chrM:-3-5");
        assert_eq!(g.fetch("1", 500000, 500003), b"AGG".to_vec(), "1:500000-500003");
        assert_eq!(g.fetch("chr1", 500000, 500000), b"".to_vec(), "chr1:500000-500000");
        assert_eq!(g.fetch("chr1", 500010, 500000), b"".to_vec(), "chr1:500010-500000");
        assert_eq!(g.fetch("chrUn_KI270742v1", 100, 160), b"ACCCCAACATTTGGACATTAAACAAAATACTTCAGAAAAACCCATTGTCCAAGAAAAGTT".to_vec(), "chrUn_KI270742v1:100-160");
        assert_eq!(g.fetch("KI270742.1", 100, 130), b"ACCCCAACATTTGGACATTAAACAAAATAC".to_vec(), "KI270742.1:100-130");
        assert_eq!(g.fetch("chrEBV", 1, 3), b"".to_vec(), "chrEBV:1-3");
        assert_eq!(g.fetch("nonexistent", 1, 10), b"".to_vec(), "nonexistent:1-10");
        assert_eq!(g.fetch("chr2", 90000000, 90000300), b"CTCACTGACATTTGTGCTTATGTGATTTTTTCAAAAAAATTCAGATGTCAATGAGAATATTGTGCCGCCTCAGTTTTATTTATTTTTATTTTTTTAACTTTTGTTTTAGGTTCAGGGATATATGTGAAGTTTTGTTACATAACTGAACTTGTGCCATGGGGGTTCCTTGTACAGATTACTTTGTCACCCAGGTATTATTCCCAGTGCCCAATAGTTATCTTTTCTGCTCCTTTCCTTTCTTCCACCCTCCACCCTCAGGTAGACCCCAGTGTGTATTGTTCCCTTATTTGTGTTCATGAG".to_vec(), "chr2:90000000-90000300");
        assert_eq!(g.fetch("chr5", 70000000, 70000020), b"CGAATGCATGCACATATAGC".to_vec(), "chr5:70000000-70000020");
        assert_eq!(g.fetch("chr1", 248956400, 248956430), b"NNNNNNNNNNNNNNNNNNNNNN".to_vec(), "chr1:248956400-248956430");
        assert_eq!(g.fetch("chr1", 248956420, 248956430), b"NN".to_vec(), "chr1:248956420-248956430");
        assert_eq!(g.fetch("chr1", 143184570, 143184610), b"NNNNNNNNNNNNNNNNNGAATTCAATGCAATCATCGAATG".to_vec(), "chr1:143184570-143184610");
        assert_eq!(g.fetch("chr1", 143184587, 143184600), b"GAATTCAATGCAA".to_vec(), "chr1:143184587-143184600");
        assert_eq!(g.fetch("chr22", 10509990, 10510010), b"NNNNNNNNNNGAATTCTTGT".to_vec(), "chr22:10509990-10510010");
        assert_eq!(g.fetch("chr1", 12345, 12346), b"C".to_vec(), "chr1:12345-12346");
        assert_eq!(g.fetch("chr1", 12346, 12349), b"AGA".to_vec(), "chr1:12346-12349");
        assert_eq!(g.fetch("chr1", 12347, 12351), b"GACC".to_vec(), "chr1:12347-12351");
        assert_eq!(g.fetch("chr1", 139999990, 140000010), b"NNNNNNNNNNNNNNNNNNNN".to_vec(), "chr1:139999990-140000010");
        assert_eq!(g.fetch("chrX", 1000000, 1000101), b"TGTAGAAACATTAGCCTGGCTAACAAGGTGAAACCCCATCTCTACTAACAATACAAAATATTGGTTGGGCGTGGTGGCGGGTGCTTGTAATCCCAGCTACT".to_vec(), "chrX:1000000-1000101");
        assert_eq!(g.fetch("chr3", 5000001, 5000002), b"C".to_vec(), "chr3:5000001-5000002");
        assert_eq!(g.fetch("chr3", 5000003, 5000004), b"G".to_vec(), "chr3:5000003-5000004");
        assert_eq!(g.fetch("chr3", 5000002, 5000005), b"TGA".to_vec(), "chr3:5000002-5000005");
        let mut names = g.seqnames().into_iter().take(5).collect::<Vec<_>>();
        assert_eq!(names.remove(0), "chr1");
    }
}
