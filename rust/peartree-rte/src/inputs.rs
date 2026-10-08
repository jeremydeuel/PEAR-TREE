//! tools/rte/inputs.py -- the combine -> annotate sidecars. FOUNDATION (implemented).
//!
//! * `<P>.insertions.evidence.tsv.gz`: one row per insertion junction, read by header name
//!   (`read_evidence_tsv`). Small (one row per junction) -> held in memory, wanted loci only.
//! * `<P>.insertions.reads.fa.gz` / `<P>.insertions.genotype_reads.fa.gz`: NEVER held whole --
//!   see `stream.rs` (`ReadStore`). This module only defines the record types and the header
//!   parser shared with it.
//!
//! Parity notes (Python behaviour kept on purpose):
//! * `cross_sample_identical` is read from a column of THAT name; combine writes
//!   `n_cross_sample_identical`, so in practice it is always 0 (inputs.py:63). Kept identical.
//! * `_int(v)` = `int(float(v))`, anything unparsable (or NaN/inf) -> 0; `_float(v)` = `float(v)`
//!   else 0.0. `clip_consensus` / `beyond_polya` are `.strip(".")`-ed, `side` upper-cased.
//! * A repeated (locus, side) row replaces the earlier one in place (Python dict semantics).

use crate::packed::PackedSeq;
use flate2::read::MultiGzDecoder;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// One evidence.tsv row (inputs.JunctionEvidence).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct JunctionEvidence {
    pub side: String,
    pub n_reads: i64,
    pub n_fragments: i64,
    pub n_independent: i64,
    pub n_samples: i64,
    pub n_mates: i64,
    pub supported: i64,
    pub clip_consensus: String,
    pub polya_len_median: f64,
    pub polya_len_range: String,
    pub beyond_polya: String,
    pub beyond_polya_support: i64,
    pub cross_sample_identical: i64,
    pub n_short_used: i64,
    pub n_short_mate_inside: i64,
}

/// One pooled read (inputs.EvidenceRead). `side` / `role` upper-cased like Python; `seq` packed.
#[derive(Clone, Debug, PartialEq)]
pub struct EvidenceRead {
    pub side: Box<str>,
    pub role: Box<str>,
    pub sample: Box<str>,
    pub frag: Box<str>,
    pub r12: Box<str>,
    pub seq: PackedSeq,
}

impl EvidenceRead {
    /// `fragment_key` = (sample, frag)
    pub fn fragment_key(&self) -> (&str, &str) {
        (&self.sample, &self.frag)
    }
    /// The sequence (upper case; Python keeps the file's case but every consumer upper-cases).
    pub fn seq(&self) -> Vec<u8> {
        self.seq.decode()
    }
    /// genotype2 extra-pass read (`GT_*` role): classification evidence only.
    pub fn is_gt(&self) -> bool {
        self.role.starts_with(GT_ROLE_PREFIX)
    }
}

pub const GT_ROLE_PREFIX: &str = "GT_";

/// inputs.InsertionEvidence: junction rows (side -> row, Python dict order) + reads.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InsertionEvidence {
    pub insertion_id: String,
    pub junctions: Vec<JunctionEvidence>,
    pub reads: Vec<EvidenceRead>,
}

impl InsertionEvidence {
    pub fn new(id: &str) -> InsertionEvidence {
        InsertionEvidence { insertion_id: id.to_string(), ..Default::default() }
    }
    /// `ev.junctions.get(side)`
    pub fn junction(&self, side: &str) -> Option<&JunctionEvidence> {
        self.junctions.iter().find(|j| j.side == side)
    }
    /// `ev.junctions[j.side] = j` (replace in place, else append)
    pub fn set_junction(&mut self, j: JunctionEvidence) {
        match self.junctions.iter_mut().find(|x| x.side == j.side) {
            Some(x) => *x = j,
            None => self.junctions.push(j),
        }
    }
}

/// A parsed reads-FASTA header `insertion_id|side|role|sample|frag|r12` (inputs.read_reads_fa
/// `flush`): a locus name may itself contain '|' (more than 6 fields -> the last five are the
/// tail); fewer than 6 fields are padded with ''.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadHeader {
    pub insertion_id: String,
    pub side: String,
    pub role: String,
    pub sample: String,
    pub frag: String,
    pub r12: String,
}

/// `name` is the header after '>' up to the first whitespace (python `line[1:].split()[0]`).
pub fn parse_read_header(name: &str) -> ReadHeader {
    let mut parts: Vec<&str> = name.split('|').collect();
    while parts.len() < 6 {
        parts.push("");
    }
    let n = parts.len();
    let (iid, tail) = if n > 6 { (parts[..n - 5].join("|"), &parts[n - 5..]) } else { (parts[0].to_string(), &parts[1..6]) };
    ReadHeader {
        insertion_id: iid,
        side: tail[0].to_uppercase(),
        role: tail[1].to_uppercase(),
        sample: tail[2].to_string(),
        frag: tail[3].to_string(),
        r12: tail[4].to_string(),
    }
}

/// Open a plain or gzipped (any number of members, like python's gzip) text file.
pub fn open_text(path: &Path) -> std::io::Result<Box<dyn BufRead + Send>> {
    let f = File::open(path)?;
    let r: Box<dyn Read + Send> = if path.to_string_lossy().ends_with(".gz") {
        Box::new(MultiGzDecoder::new(BufReader::with_capacity(1 << 20, f)))
    } else {
        Box::new(f)
    };
    Ok(Box::new(BufReader::with_capacity(1 << 20, r)))
}

/// python `_int(v)`: `int(float(v))`, else 0.
pub fn py_int_field(v: Option<&str>) -> i64 {
    match v.and_then(|s| s.trim().parse::<f64>().ok()) {
        Some(x) if x.is_finite() => x.trunc() as i64,
        _ => 0,
    }
}

/// python `_float(v)`: `float(v)`, else 0.0.
pub fn py_float_field(v: Option<&str>) -> f64 {
    v.and_then(|s| s.trim().parse::<f64>().ok()).unwrap_or(0.0)
}

fn strip_dots(s: &str) -> String {
    s.trim_matches('.').to_string()
}

impl JunctionEvidence {
    /// `JunctionEvidence.from_row(r)` with `get(col)` = the row's value for a header name.
    pub fn from_row<'a>(get: impl Fn(&str) -> Option<&'a str>) -> JunctionEvidence {
        let s = |k: &str| get(k).unwrap_or("");
        JunctionEvidence {
            side: s("side").to_uppercase(),
            n_reads: py_int_field(get("n_reads")),
            n_fragments: py_int_field(get("n_fragments")),
            n_independent: py_int_field(get("n_independent")),
            n_samples: py_int_field(get("n_samples")),
            n_mates: py_int_field(get("n_mates")),
            supported: py_int_field(get("supported")),
            clip_consensus: strip_dots(s("clip_consensus")),
            polya_len_median: py_float_field(get("polya_len_median")),
            polya_len_range: s("polya_len_range").to_string(),
            beyond_polya: strip_dots(s("beyond_polya")),
            beyond_polya_support: py_int_field(get("beyond_polya_support")),
            cross_sample_identical: py_int_field(get("cross_sample_identical")),
            n_short_used: py_int_field(get("n_short_used")),
            n_short_mate_inside: py_int_field(get("n_short_mate_inside")),
        }
    }
}

/// Junction rows per insertion (only `wanted` loci when given), in memory.
pub type EvidenceTable = FxHashMap<String, Vec<JunctionEvidence>>;

/// `read_evidence_tsv(path)` restricted to `wanted`. Blank lines skipped; the first non-blank
/// line is the header; rows without `insertion_id` skipped.
pub fn read_evidence_tsv(path: &Path, wanted: Option<&FxHashSet<String>>) -> std::io::Result<EvidenceTable> {
    let mut rdr = open_text(path)?;
    let mut out: EvidenceTable = FxHashMap::default();
    let mut header: Option<FxHashMap<String, usize>> = None;
    let mut line = String::new();
    loop {
        line.clear();
        if rdr.read_line(&mut line)? == 0 {
            break;
        }
        if line.trim().is_empty() {
            continue;
        }
        let l = line.trim_end_matches(['\n', '\r']);
        let fields: Vec<&str> = l.split('\t').collect();
        let Some(h) = &header else {
            // csv.DictReader: a duplicated column name -> the LAST one wins
            header = Some(fields.iter().enumerate().map(|(i, f)| (f.to_string(), i)).collect());
            continue;
        };
        let get = |k: &str| h.get(k).and_then(|&i| fields.get(i).copied());
        let Some(iid) = get("insertion_id").filter(|s| !s.is_empty()) else { continue };
        if let Some(w) = wanted {
            if !w.contains(iid) {
                continue;
            }
        }
        let j = JunctionEvidence::from_row(get);
        let v = out.entry(iid.to_string()).or_default();
        match v.iter_mut().find(|x| x.side == j.side) {
            Some(x) => *x = j,
            None => v.push(j),
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_with_pipes_in_locus() {
        let h = parse_read_header("chr1:5-oneside_5|x|y|LEFT|clip|S1|f7|1");
        assert_eq!(h.insertion_id, "chr1:5-oneside_5|x|y");
        assert_eq!((h.side.as_str(), h.role.as_str(), h.sample.as_str(), h.frag.as_str(), h.r12.as_str()),
                   ("LEFT", "CLIP", "S1", "f7", "1"));
        let h = parse_read_header("chr1:5-9|right|gt_mate|S1|f7|2");
        assert_eq!(h.insertion_id, "chr1:5-9");
        assert_eq!((h.side.as_str(), h.role.as_str()), ("RIGHT", "GT_MATE"));
        // fewer than 6 fields: padded
        let h = parse_read_header("loc|LEFT|CLIP");
        assert_eq!((h.insertion_id.as_str(), h.side.as_str(), h.role.as_str(), h.sample.as_str()), ("loc", "LEFT", "CLIP", ""));
        let h = parse_read_header("");
        assert_eq!(h.insertion_id, "");
    }

    #[test]
    fn evidence_tsv_by_header_name() {
        let dir = std::env::temp_dir().join(format!("rte_ev_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("ev.tsv");
        std::fs::write(&p, "\ninsertion_id\tside\tn_independent\tclip_consensus\tpolya_len_median\tn_samples\n\
            a|b\tleft\t3.0\t..acgtACGT.\t12.5\t2\n\
            a|b\tRIGHT\tnan\t\tx\n\
            c\tLEFT\t1\tAC\t\t1\n\
            \tLEFT\t1\tAC\t\t1\n\
            a|b\tLEFT\t4\tGG\t\t1\n").unwrap();
        let t = read_evidence_tsv(&p, None).unwrap();
        let a = &t["a|b"];
        assert_eq!(a.len(), 2);
        assert_eq!(a[0].side, "LEFT");
        assert_eq!(a[0].n_independent, 4); // replaced in place by the later LEFT row
        assert_eq!(a[0].clip_consensus, "GG");
        assert_eq!(a[1].side, "RIGHT");
        assert_eq!(a[1].n_independent, 0);
        assert_eq!(a[1].polya_len_median, 0.0);
        assert_eq!(a[1].n_samples, 0); // missing trailing field -> default
        let w: FxHashSet<String> = ["c".to_string()].into_iter().collect();
        let t = read_evidence_tsv(&p, Some(&w)).unwrap();
        assert_eq!(t.len(), 1);
        assert_eq!(t["c"][0].clip_consensus, "AC");
        std::fs::remove_dir_all(&dir).ok();
    }
}
