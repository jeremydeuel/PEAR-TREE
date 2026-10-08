//! tools/rte/genome.py -- genome sequence access. FOUNDATION (implemented).
//!
//! One interface ([`Genome::fetch`]: 0-based half-open, upper-case, empty when out of range or
//! the contig is unknown), two back-ends:
//! * [`TwoBit`] -- py2bit on a .2bit (`genome_2bit` / `remap_2bit`). Reader copied from
//!   rust/peartree-combine/src/genome.rs (positional reads, `Sync`, only the index + N-block
//!   tables in memory), with tools/rte's semantics on top: contig lookup exact, else the
//!   chr-prefix alternative (chr8 <-> 8); `start, end = max(0, start), min(len, end)`;
//!   `end <= start` -> "".
//! * [`FastaGenome`] -- small FASTA files (tests, fixtures, the e2e `reduced.fa`), fully in
//!   memory. A record named `contig:start-end` is a REGION of `contig` starting at `start`;
//!   `fetch` returns the slice of the FIRST region containing `start` (python semantics,
//!   including clipping at the region end).
//!
//! [`open_genome`] picks by suffix (`.2bit` -> TwoBit, else FASTA), like python `open_genome`.

use crate::inputs::open_text;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::BufRead;
use std::os::unix::fs::FileExt;
use std::path::Path;

/// python genome object (`fetch`, `length`).
pub trait Genome: Send + Sync {
    /// `genome.fetch(contig, start, end)`
    fn fetch(&self, contig: &str, start: i64, end: i64) -> Vec<u8>;
    /// `genome.length(contig)` (0 when unknown)
    fn length(&self, contig: &str) -> i64;
}

fn alt(contig: &str) -> String {
    match contig.strip_prefix("chr") {
        Some(s) => s.to_string(),
        None => format!("chr{contig}"),
    }
}

/// `open_genome(spec)`: None / "" -> None; `.2bit` -> TwoBit; anything else -> FastaGenome.
pub fn open_genome(spec: Option<&str>) -> Result<Option<Box<dyn Genome>>, String> {
    let Some(s) = spec.filter(|s| !s.is_empty()) else { return Ok(None) };
    if s.ends_with(".2bit") {
        Ok(Some(Box::new(TwoBit::open(Path::new(s))?)))
    } else {
        Ok(Some(Box::new(FastaGenome::open(&[Path::new(s)])?)))
    }
}

// ------------------------------------------------------------------------------------- 2bit

pub struct TwoBit {
    file: File,
    seqs: Vec<SeqInfo>,
    by_name: FxHashMap<String, usize>,
}

struct SeqInfo {
    len: u64,
    dna_off: u64,
    n_blocks: Vec<(u64, u64)>,
}

struct Rd<'a> {
    f: &'a File,
    swap: bool,
}

impl Rd<'_> {
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
        let a: [u8; 8] = b[..8].try_into().unwrap();
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

impl TwoBit {
    pub fn open(path: &Path) -> Result<TwoBit, String> {
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
            let dna_off = mb_off + 4 + 8 * mb + 4;
            let n_blocks = starts.iter().zip(&sizes).map(|(&s, &z)| (s as u64, s as u64 + z as u64)).collect();
            by_name.entry(name).or_insert(seqs.len());
            seqs.push(SeqInfo { len, dna_off, n_blocks });
        }
        Ok(TwoBit { file, seqs, by_name })
    }

    fn resolve(&self, contig: &str) -> Option<usize> {
        self.by_name.get(contig).or_else(|| self.by_name.get(&alt(contig))).copied()
    }

    /// py2bit `sequence(c, start, end)` for 0 <= start < end <= len (upper-case, N blocks 'N').
    fn raw(&self, sq: &SeqInfo, start: u64, end: u64) -> Vec<u8> {
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
        let first = sq.n_blocks.partition_point(|&(_, e)| e <= start);
        for &(bs, be) in &sq.n_blocks[first..] {
            if bs >= end {
                break;
            }
            for p in bs.max(start)..be.min(end) {
                out[(p - start) as usize] = b'N';
            }
        }
        out
    }
}

impl Genome for TwoBit {
    fn fetch(&self, contig: &str, start: i64, end: i64) -> Vec<u8> {
        let Some(i) = self.resolve(contig) else { return Vec::new() };
        let sq = &self.seqs[i];
        let (s, e) = (start.max(0), end.min(sq.len as i64));
        if e <= s {
            return Vec::new();
        }
        self.raw(sq, s as u64, e as u64)
    }
    fn length(&self, contig: &str) -> i64 {
        self.resolve(contig).map_or(0, |i| self.seqs[i].len as i64)
    }
}

// ------------------------------------------------------------------------------------ FASTA

/// python `FastaGenome`: contig -> [(offset, upper-case sequence)] in file order.
#[derive(Default)]
pub struct FastaGenome {
    regions: FxHashMap<String, Vec<(i64, Vec<u8>)>>,
}

impl FastaGenome {
    /// Read FASTA files (plain or gzipped); record name = first word of the header.
    pub fn open(paths: &[&Path]) -> Result<FastaGenome, String> {
        let mut g = FastaGenome::default();
        for p in paths {
            for (name, seq) in read_fasta(p)? {
                g.add(&name, &seq);
            }
        }
        Ok(g)
    }

    /// `FastaGenome(records={...})` / `add(name, seq)`.
    pub fn add(&mut self, name: &str, seq: &[u8]) {
        let (contig, off) = match parse_region(name) {
            Some((c, s, _)) => (c, s),
            None => (name.to_string(), 0),
        };
        self.regions.entry(contig).or_default().push((off, seq.to_ascii_uppercase()));
    }

    fn resolve(&self, contig: &str) -> Option<&Vec<(i64, Vec<u8>)>> {
        self.regions.get(contig).or_else(|| self.regions.get(&alt(contig)))
    }
}

/// `^(.+):(\d+)-(\d+)$` -> (contig, start, end)
pub fn parse_region(name: &str) -> Option<(String, i64, i64)> {
    let (c, rest) = name.rsplit_once(':')?;
    let (a, b) = rest.split_once('-')?;
    if c.is_empty() || a.is_empty() || b.is_empty() || !a.bytes().all(|x| x.is_ascii_digit()) || !b.bytes().all(|x| x.is_ascii_digit()) {
        return None;
    }
    Some((c.to_string(), a.parse().ok()?, b.parse().ok()?))
}

impl Genome for FastaGenome {
    fn fetch(&self, contig: &str, start: i64, end: i64) -> Vec<u8> {
        let Some(rs) = self.resolve(contig) else { return Vec::new() };
        for (off, seq) in rs {
            let n = seq.len() as i64;
            if start >= *off && start < off + n {
                let a = (start - off).max(0);
                let b = (end - off).min(n).max(0);
                return if b > a { seq[a as usize..b as usize].to_vec() } else { Vec::new() };
            }
        }
        Vec::new()
    }
    fn length(&self, contig: &str) -> i64 {
        self.resolve(contig).map_or(0, |rs| rs.iter().map(|(o, s)| o + s.len() as i64).max().unwrap_or(0))
    }
}

/// `mappy.fastx_read(path)` for FASTA: (name = first header word, sequence) in file order.
pub fn read_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>, String> {
    let rdr = open_text(path).map_err(|e| format!("cannot open {}: {e}", path.display()))?;
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in rdr.split(b'\n') {
        let mut line = line.map_err(|e| format!("{}: {e}", path.display()))?;
        while matches!(line.last(), Some(b'\r')) {
            line.pop();
        }
        if line.first() == Some(&b'>') {
            let h = String::from_utf8_lossy(&line[1..]);
            out.push((h.split_whitespace().next().unwrap_or("").to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.iter().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_regions_like_python() {
        let mut g = FastaGenome::default();
        g.add("chr1:1000-1010", b"acgtacgtac");
        g.add("chr2", b"GGGGCCCC");
        // python FastaGenome.fetch
        assert_eq!(g.fetch("chr1", 1002, 1006), b"GTAC");
        assert_eq!(g.fetch("1", 1008, 1020), b"AC"); // chr-alias, clipped at the region end
        assert_eq!(g.fetch("chr1", 999, 1005), b""); // start outside every region
        assert_eq!(g.fetch("chr2", 2, 4), b"GG");
        assert_eq!(g.fetch("chr3", 0, 4), b"");
        assert_eq!(g.length("chr1"), 1010);
        assert_eq!(parse_region("chr1:5-9"), Some(("chr1".into(), 5, 9)));
        assert_eq!(parse_region("HLA:A:5-9"), Some(("HLA:A".into(), 5, 9)));
        assert_eq!(parse_region("chr1"), None);
    }

    #[test]
    fn twobit_like_py2bit_when_available() {
        let p = "/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad/genomes/hg38.2bit";
        if !Path::new(p).exists() {
            eprintln!("skipping: {p} not found");
            return;
        }
        let g = TwoBit::open(Path::new(p)).unwrap();
        // python TwoBitGenome(hg38).fetch(...)
        assert_eq!(g.fetch("chr1", 9990, 10010), b"NNNNNNNNNNTAACCCTAAC");
        assert_eq!(g.fetch("22", 20000000, 20000037), b"AGCCTCAGGAAGGAGGCAGTGCTGCCAGCCCTTGGGG");
        assert_eq!(g.fetch("chrM", 16560, 16600), b"ATCACGATG");
        assert_eq!(g.fetch("chrM", -3, 5), b"GATCA"); // clamped (combine's get_sequence returns "")
        assert_eq!(g.fetch("nonexistent", 1, 10), b"");
        assert_eq!(g.length("chrM"), 16569);
    }
}
