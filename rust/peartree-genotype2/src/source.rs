//! Indexed region reading over BAM or CRAM, behind one trait.
//!
//! The genotyping driver only needs: the header, and an iterator of the records
//! overlapping a region. BAM and CRAM expose that through different concrete
//! reader/record types (`bam::Record` vs `sam::alignment::RecordBuf`), both of
//! which implement `sam::alignment::Record` — so `AnyRecord::as_dyn()` hands the
//! decoder a uniform `&dyn Record` and the rest of the pipeline is format-blind.
//! CRAM needs a reference FASTA (bases are stored as diffs against it).
//!
//! I/O (SPEC "Driver / I/O"): the legacy genotyper read the file unbuffered — every BGZF
//! block cost two `read` syscalls, and on Lustre each of those is a network round trip.
//! Here the `File` sits behind a [`SeekBufReader`] of `io_buffer_bytes` (4 MiB) and — unlike
//! `std::io::BufReader`, whose `seek` always discards the buffer — a seek that lands inside the
//! buffered range is served from memory. Loci are processed in coordinate order, so consecutive
//! index queries usually seek a short way forward, i.e. into data that is already buffered.
//!
//! Fill size is adaptive (`io_fill_bytes`, 256 KiB): after a seek that missed the buffer the
//! first read fetches only `io_fill_bytes`, and every further sequential fill doubles up to the
//! buffer capacity. A 17k-locus contract over a 100 GB BAM puts consecutive loci several MiB
//! apart, so with a fixed 4 MiB fill almost every window paid a 4 MiB Lustre read for the
//! ~50 KB it needed (PD37590: 7.7k reads, 31 GB, 12 min per colony — the legacy speed); the
//! small first fill makes the sparse case cheap and the doubling keeps dense regions streaming.
//!
//! BAM queries go through [`BamQuery`], which stops at the first record starting past the region
//! instead of draining every index chunk (noodles' `Query` does the latter: ~4x the BGZF blocks
//! inflated per query, and with a large buffer every extra chunk seek would refill it).

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::alignment::Record as AlignmentRecord;
use noodles_sam::Header;
use std::fs::File;
use std::io::{self, BufRead, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

/// Process-wide physical I/O counters of every [`SeekBufReader`] (file `read` calls, bytes
/// read, seeks that missed the buffer) — the numbers that matter on Lustre. Relaxed atomics,
/// one update per buffer fill / miss.
pub static IO_READS: AtomicU64 = AtomicU64::new(0);
pub static IO_BYTES: AtomicU64 = AtomicU64::new(0);
pub static IO_SEEKS: AtomicU64 = AtomicU64::new(0);

/// Snapshot (reads, bytes, seeks) of the physical I/O counters.
pub fn io_counters() -> (u64, u64, u64) {
    (IO_READS.load(Ordering::Relaxed), IO_BYTES.load(Ordering::Relaxed), IO_SEEKS.load(Ordering::Relaxed))
}

fn count_read(n: usize) {
    IO_READS.fetch_add(1, Ordering::Relaxed);
    IO_BYTES.fetch_add(n as u64, Ordering::Relaxed);
}

/// An owned record from either backend, viewable as a `&dyn Record`.
pub enum AnyRecord {
    Bam(bam::Record),
    Cram(noodles_sam::alignment::RecordBuf),
}

impl AnyRecord {
    #[inline]
    pub fn as_dyn(&self) -> &dyn AlignmentRecord {
        match self {
            AnyRecord::Bam(r) => r,
            AnyRecord::Cram(r) => r,
        }
    }
}

/// A coordinate-indexed alignment source (BAM or CRAM).
pub trait RegionSource {
    fn header(&self) -> &Header;
    /// Records overlapping `region`, in coordinate order (matching pysam fetch).
    fn query<'a>(
        &'a mut self,
        region: &Region,
    ) -> io::Result<Box<dyn Iterator<Item = io::Result<AnyRecord>> + 'a>>;
}

/// A buffered reader whose `seek` keeps the buffer when the target lies inside it.
///
/// `std::io::BufReader::seek(SeekFrom::Start(..))` always drops its buffer, so with a large
/// buffer every BGZF seek (one or more per index query) would re-read `capacity` bytes. This
/// reader tracks the file offset of its buffer and only goes to the file on a miss.
pub struct SeekBufReader<R> {
    inner: R,
    buf: Box<[u8]>,
    /// next unread byte in `buf`
    pos: usize,
    /// valid bytes in `buf`
    filled: usize,
    /// file offset of `buf[0]`; the inner reader's cursor is at `buf_start + filled`
    buf_start: u64,
    /// bytes the first fill after a buffer miss asks for (<= capacity)
    fill_min: usize,
    /// bytes the next fill asks for: `fill_min` after a miss, doubling per sequential fill
    next_fill: usize,
}

impl<R: Read + Seek> SeekBufReader<R> {
    /// Fixed fills: every fill asks for `capacity` bytes.
    #[allow(dead_code)]
    pub fn with_capacity(capacity: usize, inner: R) -> io::Result<Self> {
        Self::with_fill(capacity, capacity, inner)
    }

    /// Adaptive fills: `fill_min` bytes right after a buffer miss, doubling for each further
    /// sequential fill, never more than `capacity`.
    pub fn with_fill(capacity: usize, fill_min: usize, mut inner: R) -> io::Result<Self> {
        let buf_start = inner.stream_position()?;
        let capacity = capacity.max(1);
        let fill_min = fill_min.clamp(1, capacity);
        Ok(SeekBufReader {
            inner,
            buf: vec![0u8; capacity].into_boxed_slice(),
            pos: 0,
            filled: 0,
            buf_start,
            fill_min,
            next_fill: fill_min,
        })
    }

    #[inline]
    fn logical_position(&self) -> u64 {
        self.buf_start + self.pos as u64
    }
}

impl<R: Read> Read for SeekBufReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        // large read with an empty buffer: bypass (same policy as std's BufReader)
        if self.pos == self.filled && out.len() >= self.buf.len() {
            self.buf_start += self.filled as u64;
            self.pos = 0;
            self.filled = 0;
            let n = self.inner.read(out)?;
            count_read(n);
            self.buf_start += n as u64;
            return Ok(n);
        }
        let avail = self.fill_buf()?;
        let n = avail.len().min(out.len());
        out[..n].copy_from_slice(&avail[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl<R: Read> BufRead for SeekBufReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.pos >= self.filled {
            // the inner cursor sits at buf_start + filled: continue sequentially from there
            self.buf_start += self.filled as u64;
            self.pos = 0;
            self.filled = 0;
            let want = self.next_fill.min(self.buf.len());
            self.filled = loop {
                match self.inner.read(&mut self.buf[..want]) {
                    Ok(n) => {
                        count_read(n);
                        break n;
                    }
                    Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                    Err(e) => return Err(e),
                }
            };
            // sequential consumption: the next fill may be larger (dense region)
            self.next_fill = want.saturating_mul(2).min(self.buf.len());
        }
        Ok(&self.buf[self.pos..self.filled])
    }

    fn consume(&mut self, amt: usize) {
        self.pos = (self.pos + amt).min(self.filled);
    }
}

impl<R: Read + Seek> Seek for SeekBufReader<R> {
    fn seek(&mut self, target: SeekFrom) -> io::Result<u64> {
        let abs = match target {
            SeekFrom::Start(p) => Some(p),
            SeekFrom::Current(off) => {
                let cur = self.logical_position() as i128 + off as i128;
                if cur < 0 {
                    return Err(io::Error::new(io::ErrorKind::InvalidInput, "seek before start of file"));
                }
                Some(cur as u64)
            }
            SeekFrom::End(_) => None,
        };
        match abs {
            Some(p) if p >= self.buf_start && p <= self.buf_start + self.filled as u64 => {
                // hit: stay in memory
                self.pos = (p - self.buf_start) as usize;
                Ok(p)
            }
            Some(p) => {
                IO_SEEKS.fetch_add(1, Ordering::Relaxed);
                let p = self.inner.seek(SeekFrom::Start(p))?;
                self.buf_start = p;
                self.pos = 0;
                self.filled = 0;
                self.next_fill = self.fill_min;
                Ok(p)
            }
            None => {
                IO_SEEKS.fetch_add(1, Ordering::Relaxed);
                let p = self.inner.seek(target)?;
                self.buf_start = p;
                self.pos = 0;
                self.filled = 0;
                self.next_fill = self.fill_min;
                Ok(p)
            }
        }
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        Ok(self.logical_position())
    }
}

/// Generic over the byte reader so the buffered (`.bai`) and the plain (`.csi`-only, which
/// needs `noodles-csi`, not a direct dependency) paths share one implementation.
pub struct BamSource<R> {
    reader: bam::io::IndexedReader<bgzf::io::Reader<R>>,
    header: Header,
}

impl<R: Read + Seek> RegionSource for BamSource<R> {
    fn header(&self) -> &Header {
        &self.header
    }
    fn query<'a>(
        &'a mut self,
        region: &Region,
    ) -> io::Result<Box<dyn Iterator<Item = io::Result<AnyRecord>> + 'a>> {
        let reference_sequence_id = self.header.reference_sequences().get_index_of(region.name()).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, format!("region reference sequence does not exist in the header: {:?}", region.name()))
        })?;
        let interval = region.interval();
        // the index's own chunk list (bins + linear-index min offset), as noodles' query uses
        let chunks: Vec<(bgzf::VirtualPosition, bgzf::VirtualPosition)> =
            self.reader.index().query(reference_sequence_id, interval)?.iter().map(|c| (c.start(), c.end())).collect();
        let start = interval.start().map(usize::from).unwrap_or(1);
        let end = interval.end().map(usize::from).unwrap_or(usize::MAX);
        Ok(Box::new(BamQuery {
            reader: &mut self.reader,
            chunks: chunks.into_iter(),
            chunk_end: None,
            reference_sequence_id,
            start,
            end,
            done: false,
        }))
    }
}

/// Indexed BAM query that, unlike noodles' `Query`, STOPS at the first record starting past the
/// region end. noodles reads every chunk the index returns to its end (the leaf 16 kb bin plus the
/// parent bins' chunks), i.e. tens of KB of BGZF blocks inflated per query for a 1 kb window; in a
/// coordinate-sorted BAM every record after the first one starting beyond the region (or on a later
/// reference) cannot intersect it, so the rest of the chunk list is skipped.
pub struct BamQuery<'a, R> {
    reader: &'a mut bam::io::IndexedReader<bgzf::io::Reader<R>>,
    chunks: std::vec::IntoIter<(bgzf::VirtualPosition, bgzf::VirtualPosition)>,
    chunk_end: Option<bgzf::VirtualPosition>,
    reference_sequence_id: usize,
    /// 1-based inclusive region bounds
    start: usize,
    end: usize,
    done: bool,
}

impl<R: Read + Seek> BamQuery<'_, R> {
    fn next_record(&mut self) -> io::Result<Option<bam::Record>> {
        let mut record = bam::Record::default();
        loop {
            if self.done {
                return Ok(None);
            }
            let chunk_end = match self.chunk_end {
                Some(e) => e,
                None => match self.chunks.next() {
                    Some((s, e)) => {
                        self.reader.get_mut().seek(s)?;
                        self.chunk_end = Some(e);
                        e
                    }
                    None => {
                        self.done = true;
                        return Ok(None);
                    }
                },
            };
            if self.reader.get_ref().virtual_position() >= chunk_end {
                self.chunk_end = None;
                continue;
            }
            if self.reader.read_record(&mut record)? == 0 {
                self.done = true;
                return Ok(None);
            }
            match record.reference_sequence_id().transpose()? {
                Some(id) if id == self.reference_sequence_id => {}
                Some(id) if id < self.reference_sequence_id => continue,
                // a later reference, or the unplaced-unmapped tail: nothing more can intersect
                _ => {
                    self.done = true;
                    return Ok(None);
                }
            }
            let Some(start) = record.alignment_start().transpose()?.map(usize::from) else { continue };
            if start > self.end {
                self.done = true;
                return Ok(None);
            }
            let end = record.alignment_end().transpose()?.map(usize::from).unwrap_or(start);
            if end >= self.start {
                return Ok(Some(record));
            }
        }
    }
}

impl<R: Read + Seek> Iterator for BamQuery<'_, R> {
    type Item = io::Result<AnyRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(Some(r)) => Some(Ok(AnyRecord::Bam(r))),
            Ok(None) => None,
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

pub struct CramSource<R> {
    reader: cram::io::IndexedReader<R>,
    header: Header,
}

impl<R: Read + Seek> RegionSource for CramSource<R> {
    fn header(&self) -> &Header {
        &self.header
    }
    fn query<'a>(
        &'a mut self,
        region: &Region,
    ) -> io::Result<Box<dyn Iterator<Item = io::Result<AnyRecord>> + 'a>> {
        let it = self.reader.query(&self.header, region)?.map(|r| r.map(AnyRecord::Cram));
        Ok(Box::new(it))
    }
}

fn with_ext(path: &str, ext: &str) -> String {
    format!("{path}.{ext}")
}

pub fn is_cram(path: &str) -> bool {
    path.ends_with(".cram")
}

/// Open `path` as an indexed source with the default 4 MiB read buffer / 256 KiB first fill.
#[allow(dead_code)]
pub fn open_source(path: &str, reference: Option<&str>) -> io::Result<Box<dyn RegionSource>> {
    open_source_buffered(path, reference, 4 << 20, 256 << 10)
}

/// Open `path` as an indexed source whose file reads go through a `buffer_bytes` buffer with
/// adaptive fills starting at `fill_bytes` after every buffer miss (see the module docs).
/// CRAM (`.cram`) requires a reference FASTA (`.fai`-indexed; a `.2bit` cannot decode CRAM)
/// whose @SQ names match the CRAM header. BAM needs no reference.
pub fn open_source_buffered(
    path: &str,
    reference: Option<&str>,
    buffer_bytes: usize,
    fill_bytes: usize,
) -> io::Result<Box<dyn RegionSource>> {
    if is_cram(path) {
        let ref_path = reference.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "CRAM input requires a reference FASTA (--reference <ref.fa>)")
        })?;
        if ref_path.ends_with(".2bit") {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "CRAM input {path} needs a FASTA reference (with .fai) to decode its bases; \
                     got the 2bit file {ref_path}. Pass --reference <the BAM's assembly>.fa"
                ),
            ));
        }
        if !Path::new(&with_ext(ref_path, "fai")).exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("CRAM decoding needs a FASTA index: {} not found (samtools faidx {ref_path})", with_ext(ref_path, "fai")),
            ));
        }
        let fa = fasta::io::indexed_reader::Builder::default().build_from_path(ref_path)?;
        let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fa));
        let index = cram::crai::fs::read(with_ext(path, "crai"))?;
        let file = SeekBufReader::with_fill(buffer_bytes, fill_bytes, File::open(path)?)?;
        let mut reader = cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .set_index(index)
            .build_from_reader(file)?;
        let header = reader.read_header()?;
        Ok(Box::new(CramSource { reader, header }))
    } else {
        let bai = with_ext(path, "bai");
        if Path::new(&bai).exists() {
            let index = bam::bai::fs::read(&bai)?;
            let file = SeekBufReader::with_fill(buffer_bytes, fill_bytes, File::open(path)?)?;
            let mut reader = bam::io::indexed_reader::Builder::default().set_index(index).build_from_reader(file)?;
            let header = reader.read_header()?;
            Ok(Box::new(BamSource { reader, header }))
        } else {
            // .csi only: noodles' own path-based builder (reads the CSI) with an unbuffered file
            let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(path)?;
            let header = reader.read_header()?;
            Ok(Box::new(BamSource { reader, header }))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use std::io::Cursor;

    /// (name, flags, start) of every record noodles' own query and our early-stopping
    /// `BamQuery` return for each region — they must be identical, in order.
    fn assert_query_matches_noodles(path: &str, contig: &str, regions: &[(usize, usize)], buffer: usize) -> usize {
        let mut ours = open_source_buffered(path, None, buffer, (buffer / 4).max(1)).unwrap();
        let mut theirs = bam::io::indexed_reader::Builder::default().build_from_path(path).unwrap();
        let header = theirs.read_header().unwrap();
        let key = |r: &dyn AlignmentRecord| {
            (r.name().map(|n| n.to_vec()), r.flags().unwrap().bits(), r.alignment_start().map(|p| usize::from(p.unwrap())))
        };
        let mut total = 0;
        for &(s, e) in regions {
            let region = Region::new(contig, Position::new(s).unwrap()..=Position::new(e).unwrap());
            let a: Vec<_> = ours.query(&region).unwrap().map(|r| key(r.unwrap().as_dyn())).collect();
            let b: Vec<_> = theirs.query(&header, &region).unwrap().map(|r| key(&r.unwrap())).collect();
            assert_eq!(a, b, "region {contig}:{s}-{e}");
            total += a.len();
        }
        total
    }

    #[test]
    fn bam_query_matches_noodles_on_test_bam() {
        let bam = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test_data/test.bam");
        let mut regions = vec![(1, 100), (32_992_170, 32_992_178), (32_993_172, 32_993_172), (40_000_000, 40_001_000)];
        for s in (32_990_000..32_995_000).step_by(97) {
            for len in [1, 10, 150, 700] {
                regions.push((s, s + len));
            }
        }
        // small buffer too: exercises fills across BGZF block boundaries
        for buffer in [64, 4 << 20] {
            assert!(assert_query_matches_noodles(bam, "13", &regions, buffer) > 1000);
        }
        let mut ours = open_source_buffered(bam, None, 4 << 20, 4 << 20).unwrap();
        let region = Region::new("13", Position::new(32_992_170).unwrap()..=Position::new(32_992_178).unwrap());
        assert!(ours.query(&region).unwrap().count() > 0);
    }

    /// Same check over a larger BAM: `PEARTREE_SMOKE_BAM=<bam> PEARTREE_SMOKE_CONTIG=chr22`.
    #[test]
    #[ignore]
    fn bam_query_matches_noodles_on_smoke_bam() {
        let Ok(bam) = std::env::var("PEARTREE_SMOKE_BAM") else { return };
        let contig = std::env::var("PEARTREE_SMOKE_CONTIG").unwrap_or_else(|_| "chr22".into());
        let mut regions = Vec::new();
        let mut x: u64 = 12345;
        for _ in 0..2000 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let s = 16_000_000 + (x >> 33) as usize % 34_000_000;
            regions.push((s, s + 1 + (x >> 20) as usize % 1200));
        }
        let n = assert_query_matches_noodles(&bam, &contig, &regions, 4 << 20);
        println!("{} regions, {n} records identical", regions.len());
    }

    /// Counts inner seeks so the test can prove in-buffer seeks stay in memory.
    struct Counting {
        c: Cursor<Vec<u8>>,
        seeks: usize,
        reads: usize,
    }
    impl Read for Counting {
        fn read(&mut self, b: &mut [u8]) -> io::Result<usize> {
            self.reads += 1;
            self.c.read(b)
        }
    }
    impl Seek for Counting {
        fn seek(&mut self, p: SeekFrom) -> io::Result<u64> {
            self.seeks += 1;
            self.c.seek(p)
        }
    }

    #[test]
    fn seekbuf_reads_and_seeks_like_a_file() {
        let data: Vec<u8> = (0..1000u32).map(|i| (i % 251) as u8).collect();
        let inner = Counting { c: Cursor::new(data.clone()), seeks: 0, reads: 0 };
        let mut r = SeekBufReader::with_capacity(64, inner).unwrap();
        let seeks0 = r.inner.seeks; // stream_position() in the constructor
        let mut b = [0u8; 10];
        r.read_exact(&mut b).unwrap();
        assert_eq!(&b[..], &data[0..10]);
        // seek inside the buffered 64 bytes: no inner seek, no new read
        let reads = r.inner.reads;
        r.seek(SeekFrom::Start(40)).unwrap();
        r.read_exact(&mut b).unwrap();
        assert_eq!(&b[..], &data[40..50]);
        r.seek(SeekFrom::Start(5)).unwrap();
        r.read_exact(&mut b).unwrap();
        assert_eq!(&b[..], &data[5..15]);
        assert_eq!(r.inner.seeks, seeks0);
        assert_eq!(r.inner.reads, reads);
        assert_eq!(r.stream_position().unwrap(), 15);
        // relative seek inside the buffer
        r.seek(SeekFrom::Current(10)).unwrap();
        r.read_exact(&mut b).unwrap();
        assert_eq!(&b[..], &data[25..35]);
        // miss: one inner seek, data correct
        r.seek(SeekFrom::Start(700)).unwrap();
        assert_eq!(r.inner.seeks, seeks0 + 1);
        r.read_exact(&mut b).unwrap();
        assert_eq!(&b[..], &data[700..710]);
        // reading across the buffer end continues sequentially (no seek)
        let mut big = vec![0u8; 200];
        r.read_exact(&mut big).unwrap();
        assert_eq!(&big[..], &data[710..910]);
        assert_eq!(r.inner.seeks, seeks0 + 1);
        // read to EOF
        let mut rest = Vec::new();
        r.read_to_end(&mut rest).unwrap();
        assert_eq!(&rest[..], &data[910..]);
        // seek from end
        r.seek(SeekFrom::End(-3)).unwrap();
        let mut t = Vec::new();
        r.read_to_end(&mut t).unwrap();
        assert_eq!(&t[..], &data[997..]);
    }

    #[test]
    fn seekbuf_adaptive_fill_grows_and_resets() {
        let data: Vec<u8> = (0..4096u32).map(|i| (i % 253) as u8).collect();
        let inner = Counting { c: Cursor::new(data.clone()), seeks: 0, reads: 0 };
        let mut r = SeekBufReader::with_fill(1024, 64, inner).unwrap();
        // first fill after the (implicit) miss: 64 bytes, then 128, 256, ... up to the capacity
        let mut b = [0u8; 10];
        r.read_exact(&mut b).unwrap();
        assert_eq!(r.filled, 64);
        r.seek(SeekFrom::Start(64)).unwrap(); // end of the buffer: a hit, next fill is sequential
        r.read_exact(&mut b).unwrap();
        assert_eq!((r.buf_start, r.filled), (64, 128));
        assert_eq!(&b[..], &data[64..74]);
        r.seek(SeekFrom::Start(192)).unwrap();
        r.read_exact(&mut b).unwrap();
        assert_eq!((r.buf_start, r.filled), (192, 256));
        r.seek(SeekFrom::Start(448)).unwrap();
        r.read_exact(&mut b).unwrap();
        assert_eq!((r.buf_start, r.filled), (448, 512));
        r.seek(SeekFrom::Start(960)).unwrap();
        r.read_exact(&mut b).unwrap();
        assert_eq!((r.buf_start, r.filled), (960, 1024)); // capped at the capacity
        // a miss resets the fill to fill_min, data still correct
        let seeks = r.inner.seeks;
        r.seek(SeekFrom::Start(3000)).unwrap();
        r.read_exact(&mut b).unwrap();
        assert_eq!(r.inner.seeks, seeks + 1);
        assert_eq!((r.buf_start, r.filled), (3000, 64));
        assert_eq!(&b[..], &data[3000..3010]);
        // fixed-fill constructor keeps the old behaviour: every fill asks for the capacity
        let inner = Counting { c: Cursor::new(data.clone()), seeks: 0, reads: 0 };
        let mut r = SeekBufReader::with_capacity(256, inner).unwrap();
        r.read_exact(&mut b).unwrap();
        assert_eq!(r.filled, 256);
    }
}
