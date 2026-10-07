//! Reference genome access: FASTA (+ `.fai`) or UCSC 2bit, behind one trait. Owner: A.
//! `fetch` returns UPPERCASE bases, `N` for coordinates outside the contig.
//!
//! Both readers are random-access and keep memory small: a FASTA is read through its `.fai`
//! (seek + one `read_exact` per fetch), a 2bit through its index plus, per contig touched, the
//! N-block table (the packed bases are read only for the requested range; mask blocks are
//! ignored because the output is uppercase anyway).

// until the driver is wired in, parts of this module are unused in the binary

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};

pub trait RefSeq: Send {
    /// `genome[start, end)` of `chr`, 0-based half-open, uppercase; positions outside the contig
    /// are filled with `N` so callers never need to clamp. Unknown contig -> Err.
    fn fetch(&mut self, chr: &str, start: i64, end: i64) -> io::Result<Vec<u8>>;
    fn contig_len(&self, chr: &str) -> Option<i64>;
}

/// Open `path` as FASTA (needs `<path>.fai`) or 2bit (by extension `.2bit`).
pub fn open_reference(path: &str) -> io::Result<Box<dyn RefSeq>> {
    if path.to_ascii_lowercase().ends_with(".2bit") {
        Ok(Box::new(TwoBit::open(path)?))
    } else {
        Ok(Box::new(Fasta::open(path)?))
    }
}

/// Uppercase a base; anything that is not A/C/G/T becomes `N`.
#[inline]
pub fn norm_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        c @ (b'A' | b'C' | b'G' | b'T') => c,
        _ => b'N',
    }
}

fn bad(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// The reference's own name for `chr`: itself, else the same contig with the `chr` prefix
/// toggled. hs37d5 BAMs (and so their contract loci) say `1`, the UCSC hg19.2bit that serves
/// them says `chr1` (identical sequence for 1..22, X, Y); combine and annotate toggle the same way.
fn resolve<V>(map: &HashMap<String, V>, chr: &str) -> Option<String> {
    if map.contains_key(chr) {
        return Some(chr.to_string());
    }
    let alt = match chr.strip_prefix("chr") {
        Some(rest) => rest.to_string(),
        None => format!("chr{chr}"),
    };
    map.contains_key(&alt).then_some(alt)
}

fn unknown_contig(chr: &str) -> io::Error {
    io::Error::new(io::ErrorKind::NotFound, format!("contig '{chr}' not in the reference"))
}

/// Result buffer for `[start, end)` and the clamped in-contig range `(cs, ce)` (cs >= ce when the
/// request does not touch the contig at all). The buffer is pre-filled with `N`.
fn padded_buffer(start: i64, end: i64, len: i64) -> (Vec<u8>, i64, i64) {
    let n = (end - start).max(0) as usize;
    (vec![b'N'; n], start.max(0), end.min(len))
}

// ---------------------------------------------------------------------------------------------
// FASTA + .fai
// ---------------------------------------------------------------------------------------------

struct FaiEntry {
    len: i64,
    offset: u64,
    line_bases: i64,
    line_width: i64,
}

pub struct Fasta {
    file: File,
    index: HashMap<String, FaiEntry>,
}

impl Fasta {
    pub fn open(path: &str) -> io::Result<Fasta> {
        let fai_path = format!("{path}.fai");
        let fai = File::open(&fai_path).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("cannot open FASTA index {fai_path}: {e} (create it with `samtools faidx {path}`)"),
            )
        })?;
        let mut index = HashMap::new();
        for (i, line) in BufReader::new(fai).lines().enumerate() {
            let line = line?;
            let line = line.trim_end();
            if line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                return Err(bad(format!("{fai_path}:{}: expected 5 tab-separated fields", i + 1)));
            }
            let num = |s: &str| -> io::Result<i64> {
                s.parse().map_err(|_| bad(format!("{fai_path}:{}: bad number '{s}'", i + 1)))
            };
            let entry = FaiEntry {
                len: num(f[1])?,
                offset: num(f[2])? as u64,
                line_bases: num(f[3])?,
                line_width: num(f[4])?,
            };
            if entry.len > 0 && (entry.line_bases <= 0 || entry.line_width < entry.line_bases) {
                return Err(bad(format!("{fai_path}:{}: bad line geometry for '{}'", i + 1, f[0])));
            }
            index.insert(f[0].to_string(), entry);
        }
        let file = File::open(path).map_err(|e| io::Error::new(e.kind(), format!("cannot open reference {path}: {e}")))?;
        Ok(Fasta { file, index })
    }
}

impl RefSeq for Fasta {
    fn fetch(&mut self, chr: &str, start: i64, end: i64) -> io::Result<Vec<u8>> {
        let key = resolve(&self.index, chr).ok_or_else(|| unknown_contig(chr))?;
        let e = &self.index[&key];
        let (mut out, cs, ce) = padded_buffer(start, end, e.len);
        if cs >= ce {
            return Ok(out);
        }
        let byte_of = |p: i64| e.offset + ((p / e.line_bases) * e.line_width + p % e.line_bases) as u64;
        let from = byte_of(cs);
        let to = byte_of(ce - 1) + 1;
        let mut raw = vec![0u8; (to - from) as usize];
        self.file.seek(SeekFrom::Start(from))?;
        self.file.read_exact(&mut raw)?;
        let mut k = (cs - start) as usize;
        let mut got = 0usize;
        for &b in &raw {
            if b == b'\n' || b == b'\r' {
                continue;
            }
            out[k] = norm_base(b);
            k += 1;
            got += 1;
        }
        if got as i64 != ce - cs {
            return Err(bad(format!("FASTA '{chr}' {cs}-{ce}: .fai line geometry does not match the file")));
        }
        Ok(out)
    }

    fn contig_len(&self, chr: &str) -> Option<i64> {
        resolve(&self.index, chr).map(|k| self.index[&k].len)
    }
}

// ---------------------------------------------------------------------------------------------
// UCSC 2bit
// ---------------------------------------------------------------------------------------------

const TWOBIT_SIGNATURE: u32 = 0x1A41_2743;

struct TbDetail {
    /// half-open `(start, end)` N runs, sorted as in the file
    n_blocks: Vec<(i64, i64)>,
    /// file offset of the first packed byte
    dna_off: u64,
}

struct TbContig {
    /// file offset of the sequence record
    offset: u64,
    len: i64,
    detail: Option<TbDetail>,
}

pub struct TwoBit {
    file: File,
    big_endian: bool,
    contigs: HashMap<String, TbContig>,
}

impl TwoBit {
    pub fn open(path: &str) -> io::Result<TwoBit> {
        let file = File::open(path).map_err(|e| io::Error::new(e.kind(), format!("cannot open reference {path}: {e}")))?;
        let mut rd = BufReader::new(&file);
        let mut hdr = [0u8; 16];
        rd.read_exact(&mut hdr).map_err(|e| bad(format!("{path}: truncated 2bit header: {e}")))?;
        let big_endian = match u32::from_le_bytes([hdr[0], hdr[1], hdr[2], hdr[3]]) {
            TWOBIT_SIGNATURE => false,
            s if s.swap_bytes() == TWOBIT_SIGNATURE => true,
            _ => return Err(bad(format!("{path}: not a 2bit file (bad signature)"))),
        };
        let u32_of = |b: [u8; 4]| if big_endian { u32::from_be_bytes(b) } else { u32::from_le_bytes(b) };
        let version = u32_of([hdr[4], hdr[5], hdr[6], hdr[7]]);
        if version > 1 {
            return Err(bad(format!("{path}: unsupported 2bit version {version}")));
        }
        let count = u32_of([hdr[8], hdr[9], hdr[10], hdr[11]]) as usize;
        let mut raw_index: Vec<(String, u64)> = Vec::with_capacity(count);
        for _ in 0..count {
            let mut nl = [0u8; 1];
            rd.read_exact(&mut nl).map_err(|e| bad(format!("{path}: truncated 2bit index: {e}")))?;
            let mut name = vec![0u8; nl[0] as usize];
            rd.read_exact(&mut name).map_err(|e| bad(format!("{path}: truncated 2bit index: {e}")))?;
            let offset = if version == 1 {
                let mut b = [0u8; 8];
                rd.read_exact(&mut b).map_err(|e| bad(format!("{path}: truncated 2bit index: {e}")))?;
                if big_endian { u64::from_be_bytes(b) } else { u64::from_le_bytes(b) }
            } else {
                let mut b = [0u8; 4];
                rd.read_exact(&mut b).map_err(|e| bad(format!("{path}: truncated 2bit index: {e}")))?;
                u32_of(b) as u64
            };
            raw_index.push((String::from_utf8_lossy(&name).into_owned(), offset));
        }
        drop(rd);
        let mut tb = TwoBit { file, big_endian, contigs: HashMap::with_capacity(count) };
        // contig lengths (dnaSize, the first field of each record) are read eagerly so that
        // `contig_len(&self)` works; the block tables are loaded lazily on first `fetch`.
        for (name, offset) in raw_index {
            tb.file.seek(SeekFrom::Start(offset))?;
            let len = tb.read_u32().map_err(|e| bad(format!("{path}: truncated record of '{name}': {e}")))? as i64;
            tb.contigs.insert(name, TbContig { offset, len, detail: None });
        }
        Ok(tb)
    }

    fn read_u32(&mut self) -> io::Result<u32> {
        let mut b = [0u8; 4];
        self.file.read_exact(&mut b)?;
        Ok(if self.big_endian { u32::from_be_bytes(b) } else { u32::from_le_bytes(b) })
    }

    fn load_detail(&mut self, offset: u64) -> io::Result<TbDetail> {
        self.file.seek(SeekFrom::Start(offset + 4))?;
        let n_count = self.read_u32()? as usize;
        let mut starts = Vec::with_capacity(n_count);
        for _ in 0..n_count {
            starts.push(self.read_u32()? as i64);
        }
        let mut n_blocks = Vec::with_capacity(n_count);
        for s in starts {
            let sz = self.read_u32()? as i64;
            n_blocks.push((s, s + sz));
        }
        let mask_count = self.read_u32()? as u64;
        // skip mask starts + sizes (8 bytes per block) and the reserved word
        let dna_off = offset + 4 + 4 + 8 * n_count as u64 + 4 + 8 * mask_count + 4;
        Ok(TbDetail { n_blocks, dna_off })
    }
}

impl RefSeq for TwoBit {
    fn fetch(&mut self, chr: &str, start: i64, end: i64) -> io::Result<Vec<u8>> {
        let key = resolve(&self.contigs, chr).ok_or_else(|| unknown_contig(chr))?;
        let (offset, len, have_detail) = {
            let c = &self.contigs[&key];
            (c.offset, c.len, c.detail.is_some())
        };
        let (mut out, cs, ce) = padded_buffer(start, end, len);
        if cs >= ce {
            return Ok(out);
        }
        if !have_detail {
            let d = self.load_detail(offset)?;
            if let Some(c) = self.contigs.get_mut(&key) {
                c.detail = Some(d);
            }
        }
        let (dna_off, n_blocks) = {
            let d = self.contigs.get(&key).and_then(|c| c.detail.as_ref()).ok_or_else(|| unknown_contig(chr))?;
            (d.dna_off, d.n_blocks.clone())
        };
        let first_byte = (cs / 4) as u64;
        let last_byte = ((ce - 1) / 4) as u64;
        let mut packed = vec![0u8; (last_byte - first_byte + 1) as usize];
        self.file.seek(SeekFrom::Start(dna_off + first_byte))?;
        self.file.read_exact(&mut packed)?;
        const BASES: [u8; 4] = *b"TCAG";
        for p in cs..ce {
            let byte = packed[(p / 4 - cs / 4) as usize];
            let code = (byte >> (6 - 2 * (p % 4))) & 3;
            out[(p - start) as usize] = BASES[code as usize];
        }
        for (bs, be) in n_blocks {
            let (s, e) = (bs.max(cs), be.min(ce));
            for p in s..e {
                out[(p - start) as usize] = b'N';
            }
        }
        Ok(out)
    }

    fn contig_len(&self, chr: &str) -> Option<i64> {
        resolve(&self.contigs, chr).map(|k| self.contigs[&k].len)
    }
}

// ---------------------------------------------------------------------------------------------
// in-memory reference for tests (also usable by other modules' tests)
// ---------------------------------------------------------------------------------------------

#[cfg(test)]
pub struct MemRef(pub HashMap<String, Vec<u8>>);

#[cfg(test)]
impl MemRef {
    pub fn new(contigs: &[(&str, &[u8])]) -> MemRef {
        MemRef(contigs.iter().map(|(n, s)| (n.to_string(), s.iter().map(|&b| norm_base(b)).collect())).collect())
    }
}

#[cfg(test)]
impl RefSeq for MemRef {
    fn fetch(&mut self, chr: &str, start: i64, end: i64) -> io::Result<Vec<u8>> {
        let s = self.0.get(chr).ok_or_else(|| unknown_contig(chr))?;
        let (mut out, cs, ce) = padded_buffer(start, end, s.len() as i64);
        for p in cs..ce.max(cs) {
            out[(p - start) as usize] = s[p as usize];
        }
        Ok(out)
    }
    fn contig_len(&self, chr: &str) -> Option<i64> {
        self.0.get(chr).map(|s| s.len() as i64)
    }
}

/// deterministic pseudo-random ACGT (xorshift) of length n
#[cfg(test)]
pub fn rand_seq(n: usize, mut seed: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed >> 33) as usize % 4]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;

    fn tmpdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("pt_g2_{tag}_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn truth() -> Vec<(String, Vec<u8>)> {
        let mut c1 = rand_seq(1003, 7);
        for b in &mut c1[100..130] {
            *b = b'N';
        }
        let c2 = rand_seq(61, 99);
        vec![("chr1".to_string(), c1), ("HLA-A*01:01:01:01".to_string(), c2)]
    }

    fn check_reader(r: &mut dyn RefSeq, truth: &[(String, Vec<u8>)]) {
        for (name, s) in truth {
            assert_eq!(r.contig_len(name), Some(s.len() as i64));
            let len = s.len() as i64;
            // whole contig, interior ranges, byte/line boundaries, padded ends
            for &(a, b) in &[(0, len), (3, 4), (0, 1), (len - 1, len), (50, 131), (97, 133), (-5, 20), (len - 7, len + 9), (-10, len + 10), (4, 4), (len + 5, len + 9), (-9, -3)] {
                let got = r.fetch(name, a, b).unwrap();
                let want: Vec<u8> = (a..b.max(a)).map(|p| if p < 0 || p >= len { b'N' } else { s[p as usize] }).collect();
                assert_eq!(got, want, "{name}:{a}-{b}");
            }
        }
        assert!(r.fetch("nope", 0, 5).is_err());
        assert_eq!(r.contig_len("nope"), None);
    }

    #[test]
    fn fasta_with_fai_and_padding() {
        let d = tmpdir("fa");
        let t = truth();
        let fa = d.join("mini.fa");
        let mut f = File::create(&fa).unwrap();
        let mut fai = String::new();
        let mut pos = 0u64;
        for (i, (name, s)) in t.iter().enumerate() {
            // lowercase part + CRLF-free 60-col lines; second contig exactly one line wide
            let hdr = format!(">{name} some description\n");
            f.write_all(hdr.as_bytes()).unwrap();
            pos += hdr.len() as u64;
            let lw = 60usize;
            let mut body = Vec::new();
            for ch in s.chunks(lw) {
                let mut line = ch.to_vec();
                if i == 0 {
                    line.make_ascii_lowercase();
                }
                body.extend_from_slice(&line);
                body.push(b'\n');
            }
            f.write_all(&body).unwrap();
            fai.push_str(&format!("{name}\t{}\t{pos}\t{lw}\t{}\n", s.len(), lw + 1));
            pos += body.len() as u64;
        }
        std::fs::write(d.join("mini.fa.fai"), fai).unwrap();
        let mut r = open_reference(fa.to_str().unwrap()).unwrap();
        check_reader(r.as_mut(), &t);
        // missing .fai is a clear error
        std::fs::copy(&fa, d.join("nofai.fa")).unwrap();
        assert!(open_reference(d.join("nofai.fa").to_str().unwrap()).is_err());
    }

    /// Build a 2bit file (version 0) with the given byte order. Lowercase input bases become mask
    /// blocks (ignored by the reader), N runs become N blocks.
    fn make_2bit(seqs: &[(String, Vec<u8>)], be: bool, masked: bool) -> Vec<u8> {
        let w32 = |v: &mut Vec<u8>, x: u32| v.extend_from_slice(&if be { x.to_be_bytes() } else { x.to_le_bytes() });
        let mut records: Vec<Vec<u8>> = Vec::new();
        for (_, s) in seqs {
            let mut r = Vec::new();
            w32(&mut r, s.len() as u32);
            let mut blocks = Vec::new();
            let mut i = 0;
            while i < s.len() {
                if s[i] == b'N' {
                    let j = s[i..].iter().position(|&b| b != b'N').map_or(s.len(), |k| i + k);
                    blocks.push((i as u32, (j - i) as u32));
                    i = j;
                } else {
                    i += 1;
                }
            }
            w32(&mut r, blocks.len() as u32);
            for b in &blocks {
                w32(&mut r, b.0);
            }
            for b in &blocks {
                w32(&mut r, b.1);
            }
            if masked {
                w32(&mut r, 2);
                w32(&mut r, 10);
                w32(&mut r, 40);
                w32(&mut r, 5);
                w32(&mut r, 7);
            } else {
                w32(&mut r, 0);
            }
            w32(&mut r, 0);
            let mut packed = vec![0u8; s.len().div_ceil(4)];
            for (p, &b) in s.iter().enumerate() {
                let code = match b {
                    b'T' => 0u8,
                    b'C' => 1,
                    b'A' => 2,
                    b'G' => 3,
                    _ => 0, // N: packed as T, covered by the N block
                };
                packed[p / 4] |= code << (6 - 2 * (p % 4));
            }
            r.extend_from_slice(&packed);
            records.push(r);
        }
        let mut out = Vec::new();
        w32(&mut out, TWOBIT_SIGNATURE);
        w32(&mut out, 0);
        w32(&mut out, seqs.len() as u32);
        w32(&mut out, 0);
        let index_len: usize = seqs.iter().map(|(n, _)| 1 + n.len() + 4).sum();
        let mut off = (16 + index_len) as u32;
        for ((n, _), rec) in seqs.iter().zip(&records) {
            out.push(n.len() as u8);
            out.extend_from_slice(n.as_bytes());
            w32(&mut out, off);
            off += rec.len() as u32;
        }
        for rec in records {
            out.extend_from_slice(&rec);
        }
        out
    }

    #[test]
    fn twobit_both_byte_orders_and_padding() {
        let d = tmpdir("tb");
        let t = truth();
        for (be, masked) in [(false, false), (true, false), (false, true), (true, true)] {
            let p = d.join(format!("mini_{be}_{masked}.2bit"));
            std::fs::write(&p, make_2bit(&t, be, masked)).unwrap();
            let mut r = open_reference(p.to_str().unwrap()).unwrap();
            check_reader(r.as_mut(), &t);
        }
        let p = d.join("junk.2bit");
        std::fs::write(&p, b"this is not a 2bit file at all").unwrap();
        assert!(open_reference(p.to_str().unwrap()).is_err());
    }

    #[test]
    fn chr_prefix_is_toggled_for_hs37d5_names() {
        // hs37d5 contract loci say `1`; the hg19.2bit serving them says `chr1` (and vice versa)
        let d = tmpdir("alias");
        let t = truth();
        let p = d.join("alias.2bit");
        std::fs::write(&p, make_2bit(&t, false, false)).unwrap();
        let mut r = open_reference(p.to_str().unwrap()).unwrap();
        assert_eq!(r.contig_len("1"), Some(t[0].1.len() as i64));
        assert_eq!(r.fetch("1", 200, 260).unwrap(), r.fetch("chr1", 200, 260).unwrap());
        assert!(r.fetch("chrHLA-A*01:01:01:01", 0, 5).is_ok());
        assert!(r.fetch("2", 0, 5).is_err());
    }

    #[test]
    fn norm_base_maps_iupac_to_n() {
        assert_eq!(norm_base(b'a'), b'A');
        assert_eq!(norm_base(b'R'), b'N');
        assert_eq!(norm_base(b'n'), b'N');
        assert_eq!(norm_base(b'-'), b'N');
    }
}
