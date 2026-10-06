//! Indexed region reading over BAM or CRAM, behind one trait.
//!
//! The genotyping driver only needs: the header, and an iterator of the records
//! overlapping a region. BAM and CRAM expose that through different concrete
//! reader/record types (`bam::Record` vs `sam::alignment::RecordBuf`), both of
//! which implement `sam::alignment::Record` — so `AnyRecord::as_dyn()` hands the
//! decoder a uniform `&dyn Record` and the rest of the pipeline is format-blind.
//! CRAM needs a reference FASTA (bases are stored as diffs against it).

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::alignment::Record as AlignmentRecord;
use noodles_sam::Header;
use std::fs::File;
use std::io;

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

pub struct BamSource {
    reader: bam::io::IndexedReader<bgzf::io::Reader<File>>,
    header: Header,
}

impl RegionSource for BamSource {
    fn header(&self) -> &Header {
        &self.header
    }
    fn query<'a>(
        &'a mut self,
        region: &Region,
    ) -> io::Result<Box<dyn Iterator<Item = io::Result<AnyRecord>> + 'a>> {
        // disjoint field borrows: reader mutable, header shared, same lifetime.
        let it = self.reader.query(&self.header, region)?.map(|r| r.map(AnyRecord::Bam));
        Ok(Box::new(it))
    }
}

pub struct CramSource {
    reader: cram::io::IndexedReader<File>,
    header: Header,
}

impl RegionSource for CramSource {
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

/// Open `path` as an indexed source. CRAM (`.cram`) requires a reference FASTA
/// (`.fai`-indexed) whose @SQ names match the CRAM header.
pub fn open_source(path: &str, reference: Option<&str>) -> io::Result<Box<dyn RegionSource>> {
    if path.ends_with(".cram") {
        let ref_path = reference.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "CRAM input requires a reference FASTA (--reference <ref.fa>)")
        })?;
        let fa = fasta::io::indexed_reader::Builder::default().build_from_path(ref_path)?;
        let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fa));
        let mut reader = cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(path)?;
        let header = reader.read_header()?;
        Ok(Box::new(CramSource { reader, header }))
    } else {
        let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(path)?;
        let header = reader.read_header()?;
        Ok(Box::new(BamSource { reader, header }))
    }
}
