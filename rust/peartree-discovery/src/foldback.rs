//! SPEC-9: cruciform fold-back junction gate (`foldback_filter`).
//!
//! Low-input libraries made by enzymatic fragmentation turn cruciform DNA at inverted
//! repeats into hairpin fragments: the resolved hairpin carries a short stretch of the
//! *opposite strand* of the sequence beside it (Ellis et al. 2021 Nat Protoc 16:841,
//! Fig. 3g). Aligned, such a read is soft-clipped at the fold, and the clip is the reverse
//! complement of reference a few bases away. Many fragments fold at near-identical points,
//! so duplicate marking misses them and the clip cluster passes every fragment floor
//! (PD45886b_lo0019, 12:81477172: 84 reads in 170 bp, mates sharing a start with mirror
//! clips `97M54S`/`54S97M`). A real insertion's clip is element or transduced sequence, not
//! an inverted copy of its own flank — measured 2/2573 germline MEIs within ±50 bp at k=20.
//!
//! The reference is read from a UCSC .2bit by random access (a few bytes per breakpoint);
//! nothing is held in memory but the sequence index.

use std::cell::RefCell;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};

use rustc_hash::FxHashMap;

use crate::filters::{clip_entropy, longest_homopolymer_run};
use crate::config::CLIP_LEFT;

/// A probe below this entropy (bits) or with a homopolymer run >= `PROBE_MAX_RUN` is not
/// judged: poly-A/T and STR clips match inverted reference tracts by chance.
const PROBE_MIN_ENTROPY: f64 = 1.5;
const PROBE_MAX_RUN: usize = 8;

/// Random-access reader for a UCSC .2bit file (versions 0 and 1, either byte order).
pub struct TwoBit {
    file: RefCell<BufReader<File>>,
    /// contig -> (offset of the packed DNA, length)
    seqs: FxHashMap<String, (u64, u64)>,
}

impl TwoBit {
    pub fn open(path: &str) -> io::Result<TwoBit> {
        let mut f = BufReader::new(File::open(path)?);
        let mut hdr = [0u8; 16];
        f.read_exact(&mut hdr)?;
        let le = match u32::from_le_bytes(hdr[0..4].try_into().unwrap()) {
            0x1A41_2743 => true,
            0x4327_411A => false,
            _ => return Err(io::Error::new(io::ErrorKind::InvalidData, format!("{path}: not a .2bit file"))),
        };
        let u32_at = |b: &[u8]| {
            let a: [u8; 4] = b.try_into().unwrap();
            if le { u32::from_le_bytes(a) } else { u32::from_be_bytes(a) }
        };
        let version = u32_at(&hdr[4..8]);
        let count = u32_at(&hdr[8..12]) as usize;
        let mut index = Vec::with_capacity(count);
        for _ in 0..count {
            let mut n = [0u8; 1];
            f.read_exact(&mut n)?;
            let mut name = vec![0u8; n[0] as usize];
            f.read_exact(&mut name)?;
            let off = if version == 1 {
                let mut b = [0u8; 8];
                f.read_exact(&mut b)?;
                if le { u64::from_le_bytes(b) } else { u64::from_be_bytes(b) }
            } else {
                let mut b = [0u8; 4];
                f.read_exact(&mut b)?;
                u32_at(&b) as u64
            };
            index.push((String::from_utf8_lossy(&name).into_owned(), off));
        }
        let mut seqs = FxHashMap::default();
        for (name, off) in index {
            // record: dnaSize, nBlockCount, nBlockStarts[], nBlockSizes[], maskBlockCount,
            // maskBlockStarts[], maskBlockSizes[], reserved, packedDna
            f.seek(SeekFrom::Start(off))?;
            let mut b = [0u8; 8];
            f.read_exact(&mut b)?;
            let (size, n_blocks) = (u32_at(&b[0..4]) as u64, u32_at(&b[4..8]) as u64);
            f.seek(SeekFrom::Current((8 * n_blocks) as i64))?;
            let mut m = [0u8; 4];
            f.read_exact(&mut m)?;
            let mask_blocks = u32_at(&m) as u64;
            let dna = off + 8 + 8 * n_blocks + 4 + 8 * mask_blocks + 4;
            seqs.insert(name, (dna, size));
        }
        Ok(TwoBit { file: RefCell::new(f), seqs })
    }

    /// The .2bit name for a BAM contig: as is, else with `chr` added or removed
    /// (hs37d5 `12` -> hg19 `chr12`); `MT`/`chrM` never match (rCRS vs hg19's chrM).
    fn resolve(&self, contig: &str) -> Option<(u64, u64)> {
        if matches!(contig, "MT" | "M" | "chrM" | "chrMT") {
            return None;
        }
        if let Some(&v) = self.seqs.get(contig) {
            return Some(v);
        }
        let alt = match contig.strip_prefix("chr") {
            Some(s) => s.to_string(),
            None => format!("chr{contig}"),
        };
        self.seqs.get(&alt).copied()
    }

    pub fn len(&self, contig: &str) -> Option<u64> {
        self.resolve(contig).map(|(_, l)| l)
    }

    /// Uppercase reference bases [start, end) (0-based, clamped to the contig). N blocks
    /// read as the packed placeholder (T); harmless here, a T run is never a judged probe.
    pub fn fetch(&self, contig: &str, start: i64, end: i64) -> Option<Vec<u8>> {
        let (dna, len) = self.resolve(contig)?;
        let s = start.max(0) as u64;
        let e = (end.max(0) as u64).min(len);
        if s >= e {
            return Some(Vec::new());
        }
        let (b0, b1) = (s / 4, e.div_ceil(4));
        let mut packed = vec![0u8; (b1 - b0) as usize];
        let mut f = self.file.borrow_mut();
        f.seek(SeekFrom::Start(dna + b0)).ok()?;
        f.read_exact(&mut packed).ok()?;
        const BASES: [u8; 4] = *b"TCAG";
        Some((s..e).map(|p| {
            let byte = packed[(p / 4 - b0) as usize];
            BASES[((byte >> (6 - 2 * (p % 4))) & 3) as usize]
        }).collect())
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// The `k` clip bases next to the junction, uppercased (`CLIP_LEFT` stores
/// `[clipped][unclipped]`, so the end of the clip; `CLIP_RIGHT` the start). `None` when the
/// clip is shorter than `k`, holds a non-ACGT base, or is low-complexity.
pub fn probe(side: i32, clipped: &[u8], k: usize) -> Option<Vec<u8>> {
    if k == 0 || clipped.len() < k {
        return None;
    }
    let p: Vec<u8> = if side == CLIP_LEFT { &clipped[clipped.len() - k..] } else { &clipped[..k] }
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    if !p.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
        return None;
    }
    if clip_entropy(&p) < PROBE_MIN_ENTROPY || longest_homopolymer_run(&p).0 >= PROBE_MAX_RUN {
        return None;
    }
    Some(p)
}

/// True if the probe's reverse complement occurs in `reference` (the window around the
/// breakpoint).
pub fn is_foldback(probe: &[u8], reference: &[u8]) -> bool {
    let rc = revcomp(probe);
    reference.windows(rc.len()).any(|w| w == rc.as_slice())
}

/// Gate for one breakpoint: true = fold-back (drop). Contigs the .2bit lacks are kept.
pub fn breakpoint_is_foldback(
    tb: &TwoBit,
    contig: &str,
    side: i32,
    clipped: &[u8],
    breakpoint: i64,
    k: usize,
    window: i64,
) -> bool {
    let Some(p) = probe(side, clipped, k) else { return false };
    match tb.fetch(contig, breakpoint - window, breakpoint + window + 1) {
        Some(r) => is_foldback(&p, &r),
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::CLIP_RIGHT;

    /// Pack `seqs` into a minimal version-0 little-endian .2bit (no N or mask blocks).
    fn write_2bit(path: &std::path::Path, seqs: &[(&str, &[u8])]) {
        let mut out = Vec::new();
        out.extend_from_slice(&0x1A41_2743u32.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        out.extend_from_slice(&(seqs.len() as u32).to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        let idx_len: usize = seqs.iter().map(|(n, _)| 1 + n.len() + 4).sum();
        let mut off = 16 + idx_len;
        let mut recs = Vec::new();
        for (name, s) in seqs {
            out.push(name.len() as u8);
            out.extend_from_slice(name.as_bytes());
            out.extend_from_slice(&(off as u32).to_le_bytes());
            let mut r = Vec::new();
            r.extend_from_slice(&(s.len() as u32).to_le_bytes());
            r.extend_from_slice(&0u32.to_le_bytes()); // nBlockCount
            r.extend_from_slice(&0u32.to_le_bytes()); // maskBlockCount
            r.extend_from_slice(&0u32.to_le_bytes()); // reserved
            for chunk in s.chunks(4) {
                let mut byte = 0u8;
                for (i, b) in chunk.iter().enumerate() {
                    let code = match b { b'T' => 0, b'C' => 1, b'A' => 2, _ => 3 };
                    byte |= code << (6 - 2 * i);
                }
                r.push(byte);
            }
            off += r.len();
            recs.push(r);
        }
        for r in recs {
            out.extend_from_slice(&r);
        }
        std::fs::write(path, out).unwrap();
    }

    // hg19 chr12:[81477142, 81477242) (0-based), the PD45886 lo0019 fold-back site; the
    // "TSD" GTCACCCAGCCTGGAGTGCAATG starts at 81477172
    const SITE: &[u8] = b"TTTTTTTCTTTTGAGACGGAGTCTCACTCTGTCACCCAGCCTGGAGTGCAATGGCATGATCTCTGCTCACTGCAACCTCCATCTCCCCGCTTCAACCATT";

    fn tb() -> (tempdir::Dir, TwoBit) {
        let d = tempdir::Dir::new();
        let p = d.0.join("t.2bit");
        let mut chr = vec![b'G'; 81477142];
        chr.extend_from_slice(SITE);
        chr.extend(std::iter::repeat_n(b'C', 300));
        write_2bit(&p, &[("chr12", &chr), ("chr1", b"ACGTACGTAC")]);
        let t = TwoBit::open(p.to_str().unwrap()).unwrap();
        (d, t)
    }

    mod tempdir {
        pub struct Dir(pub std::path::PathBuf);
        impl Dir {
            pub fn new() -> Dir {
                let p = std::env::temp_dir().join(format!("foldback-test-{}-{:?}", std::process::id(), std::thread::current().id()));
                std::fs::create_dir_all(&p).unwrap();
                Dir(p)
            }
        }
        impl Drop for Dir {
            fn drop(&mut self) {
                let _ = std::fs::remove_dir_all(&self.0);
            }
        }
    }

    #[test]
    fn twobit_fetch_and_chr_toggle() {
        let (_d, t) = tb();
        assert_eq!(t.fetch("chr1", 0, 10).unwrap(), b"ACGTACGTAC");
        assert_eq!(t.fetch("1", 2, 6).unwrap(), b"GTAC", "hs37d5 name resolves to chr1");
        assert_eq!(t.fetch("chr1", 8, 50).unwrap(), b"AC", "clamped to the contig end");
        assert_eq!(t.fetch("12", 81477142, 81477142 + 10).unwrap(), &SITE[..10]);
        assert!(t.fetch("MT", 0, 5).is_none());
        assert!(t.fetch("GL000192.1", 0, 5).is_none());
        assert_eq!(t.len("12"), Some(81477142 + SITE.len() as u64 + 300));
    }

    #[test]
    fn lo0019_left_clip_is_a_foldback() {
        let (_d, t) = tb();
        // combine's L junction: clip `gggtgacagagtgagactcc` | aligned GTCACCCAGCC... at 81477172
        let clipped = b"TGGGTGACAGAGTGAGACTCC";
        assert!(breakpoint_is_foldback(&t, "12", CLIP_LEFT, clipped, 81477172, 20, 50));
        // and the R junction: aligned ...GTGCAATG | clip `agcagagatcatgccatt` (18 bp)
        assert!(breakpoint_is_foldback(&t, "12", CLIP_RIGHT, b"AGCAGAGATCATGCCATT", 81477195, 18, 50));
    }

    /// The real UCSC file, when available (PEARTREE_TEST_HG19_2BIT=/path/hg19.2bit); skipped
    /// otherwise. Pins the reader to UCSC's packing, not just to `write_2bit` above.
    #[test]
    fn real_hg19_lo0019_site() {
        let Ok(p) = std::env::var("PEARTREE_TEST_HG19_2BIT") else { return };
        let t = TwoBit::open(&p).unwrap();
        assert_eq!(t.len("12"), Some(133_851_895));
        assert_eq!(t.fetch("12", 81477142, 81477242).unwrap(), SITE);
        assert!(breakpoint_is_foldback(&t, "12", CLIP_LEFT, b"TGGGTGACAGAGTGAGACTCC", 81477172, 20, 50));
    }

    #[test]
    fn element_clip_and_low_complexity_are_kept() {
        let (_d, t) = tb();
        // an Alu 5' end clip (GGCCGGGCGCGGTGGCTCACGCCTG) is not an inverted copy of this flank
        assert!(!breakpoint_is_foldback(&t, "12", CLIP_RIGHT, b"GGCCGGGCGCGGTGGCTCACGCCTG", 81477195, 20, 50));
        // a poly-A clip is never judged
        assert!(probe(CLIP_LEFT, b"AAAAAAAAAAAAAAAAAAAAAAAA", 20).is_none());
        // too short to judge
        assert!(probe(CLIP_LEFT, b"GGGTGACAGAG", 20).is_none());
        // out of the window
        assert!(!breakpoint_is_foldback(&t, "12", CLIP_LEFT, b"TGGGTGACAGAGTGAGACTCC", 81477172 + 5000, 20, 50));
    }
}
