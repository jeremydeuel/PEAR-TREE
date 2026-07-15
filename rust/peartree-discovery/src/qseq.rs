//! Port of src/quality_seq.py + src/revcomp.py
//!
//! A DNA sequence (ASCII bytes) paired with a per-base quality/score track
//! (i32, because consensus scores can exceed the 0-255 phred range).

#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'G' => b'C', b'C' => b'G',
        b'a' => b't', b't' => b'a', b'g' => b'c', b'c' => b'g',
        other => other,
    }
}

pub fn revcomp_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[derive(Clone, Debug)]
pub struct QualitySeq {
    pub seq: Vec<u8>,
    pub qual: Vec<i32>,
}

impl QualitySeq {
    pub fn new(seq: Vec<u8>, qual: Vec<i32>) -> Self {
        assert_eq!(seq.len(), qual.len());
        QualitySeq { seq, qual }
    }

    pub fn empty() -> Self {
        QualitySeq { seq: Vec::new(), qual: Vec::new() }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Python slice `self[start:stop]` with step 1: negative indices count from
    /// the end, out-of-range indices clamp, start>stop yields empty. Never panics.
    pub fn pyslice(&self, start: Option<isize>, stop: Option<isize>) -> QualitySeq {
        let len = self.seq.len() as isize;
        let mut s = start.unwrap_or(0);
        let mut e = stop.unwrap_or(len);
        if s < 0 { s += len; }
        if e < 0 { e += len; }
        s = s.clamp(0, len);
        e = e.clamp(0, len);
        if s > e { s = e; }
        let (s, e) = (s as usize, e as usize);
        QualitySeq { seq: self.seq[s..e].to_vec(), qual: self.qual[s..e].to_vec() }
    }

    pub fn revcomp(&self) -> QualitySeq {
        let seq = revcomp_bytes(&self.seq);
        let mut qual = self.qual.clone();
        qual.reverse();
        QualitySeq { seq, qual }
    }

    pub fn upper(&self) -> QualitySeq {
        QualitySeq { seq: self.seq.iter().map(|b| b.to_ascii_uppercase()).collect(), qual: self.qual.clone() }
    }

    pub fn concat(&self, other: &QualitySeq) -> QualitySeq {
        let mut seq = self.seq.clone();
        seq.extend_from_slice(&other.seq);
        let mut qual = self.qual.clone();
        qual.extend_from_slice(&other.qual);
        QualitySeq { seq, qual }
    }

    /// base at position (like `.seq(i)`)
    pub fn base(&self, i: usize) -> u8 {
        self.seq[i]
    }

    /// `str(self) == other`
    pub fn eq_bytes(&self, other: &[u8]) -> bool {
        self.seq == other
    }

    /// FASTQ record, matching QualitySeq.fastq(): phred base 33, clamped at 126.
    pub fn fastq(&self, title: &str) -> String {
        let mut out = String::with_capacity(self.seq.len() * 2 + title.len() + 6);
        out.push('@');
        out.push_str(title);
        out.push('\n');
        out.push_str(std::str::from_utf8(&self.seq).unwrap());
        out.push_str("\n+\n");
        for &q in &self.qual {
            let mut c = q + 33;
            if c > 126 { c = 126; }
            out.push(c as u8 as char);
        }
        out.push('\n');
        out
    }
}
