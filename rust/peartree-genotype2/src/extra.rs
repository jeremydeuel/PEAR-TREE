//! Extra evidence at loci a colony carries but did not discover (`gt_extra_reads`, default off).
//!
//! Discovery reports a colony for a locus only when its clip evidence clears the discovery
//! gates; the joint step later calls many more colonies carriers. Their BAMs hold further
//! junction reads and discordant pairs whose inside mates show more of the inserted element
//! (5' truncation point, inversion, transduction). For every locus where this colony
//!   (1) is NOT a discovery member (`--members`: the samples combine kept reads from),
//!   (1b) is NOT germline: discovered in at most `gt_extra_max_member_frac` (0.5) of the patient's
//!       colonies (`#colonies N` of the members table = discovery files of the run); a locus
//!       most colonies discovered is already well sampled, and its carriers would only spend
//!       the mate budget, and
//!   (2) passes the ALT gate: >= `gt_extra_min_alt` evidence reads realigned as `Alt` (the
//!       genotyper's own read assignment: LLR >= `llr_informative`, explained >=
//!       `min_explained_frac`),
//! the driver appends to `<out>.extra_reads.fa.gz`, in the insertions.reads.fa.gz conventions
//! (`>locus|SIDE|ROLE|sample|frag|r12`, allele/site-forward sequence):
//!   GT_CLIP / GT_POLYA  the Alt reads as stored (aligned at the site = site-forward); GT_POLYA
//!                       when the junction-facing soft clip is a pure poly-A/T
//!   GT_DISC             discordant anchors in the flank pointing at a junction (last base within
//!                       `gt_extra_disc_span` before R / first base within it after L, MAPQ >=
//!                       `gt_extra_anchor_mapq`, primary, not 0x400)
//!   GT_MATE             their inside mates, fetched by random access through the index, oriented
//!                       like discovery's `mate_site_forward` / combine's `allele_forward_seq`
//! The GT_ prefix keeps these out of every evidence rule that reads the combine roles
//! (cluster/somatic_table.py's junction_support). The genotype rows never change: the extra
//! pass runs after the locus is called, with its own query.

use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

pub const ROLE_CLIP: &str = "GT_CLIP";
pub const ROLE_POLYA: &str = "GT_POLYA";
pub const ROLE_DISC: &str = "GT_DISC";
pub const ROLE_MATE: &str = "GT_MATE";

/// Sidecar of a per-colony genotype file.
pub fn sidecar_path(out: &str) -> String {
    format!("{out}.extra_reads.fa.gz")
}

/// Colony name of an output path (`genotypes/<SAMPLE>.txt.gz` -> `SAMPLE`), the same stem the
/// joint step uses for the matrix columns.
pub fn sample_from_out(out: &str) -> String {
    let base = std::path::Path::new(out).file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
    let s = base.strip_suffix(".gz").unwrap_or(&base);
    s.strip_suffix(".txt").or_else(|| s.strip_suffix(".tsv")).unwrap_or(s).to_string()
}

/// Discovery memberships as the extra pass needs them.
#[derive(Clone, Debug, Default)]
pub struct Members {
    /// loci this colony is a discovery member of (its reads are in insertions.reads.fa.gz)
    pub mine: FxHashSet<String>,
    /// loci discovered in more than `gt_extra_max_member_frac` of the patient's colonies
    /// (germline / near-clonal: their structure is already well sampled; empty without
    /// `#colonies`)
    pub germline: FxHashSet<String>,
    /// the patient's colony count (`#colonies N` header of the members table)
    pub n_colonies: Option<usize>,
}

/// Why a locus gets (or does not get) the extra pass.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Candidate {
    Yes,
    /// this colony discovered it
    Member,
    /// discovered in > max_member_frac of all colonies
    Germline,
}

impl Members {
    pub fn candidate(&self, locus: &str) -> Candidate {
        if self.mine.contains(locus) {
            Candidate::Member
        } else if self.germline.contains(locus) {
            Candidate::Germline
        } else {
            Candidate::Yes
        }
    }
}

/// `n_members / n_colonies > max_frac` (strict).
pub fn is_germline(n_members: usize, n_colonies: usize, max_frac: f64) -> bool {
    n_colonies > 0 && n_members as f64 > max_frac * n_colonies as f64
}

/// The discovery memberships of `sample`. `path` is either the members table built once per
/// patient by tools/genotype_extra_reads.py (`#colonies N` = colonies in the run, then
/// `locus<TAB>sample,sample,...`), or combine's `insertions.reads.fa.gz` itself (headers
/// `locus|SIDE|ROLE|sample|frag|r12`; a locus name containing `|` keeps the last five fields as
/// the tail; no colony count there, so no germline filter). Gzip or plain.
pub fn load_members(path: &str, sample: &str, max_member_frac: f64) -> io::Result<Members> {
    let f = File::open(path).map_err(|e| io::Error::new(e.kind(), format!("cannot open {path}: {e}")))?;
    let r: Box<dyn Read> = if path.ends_with(".gz") { Box::new(MultiGzDecoder::new(f)) } else { Box::new(f) };
    members_from(BufReader::with_capacity(1 << 20, r), sample, max_member_frac)
}

pub(crate) fn members_from<R: BufRead>(reader: R, sample: &str, max_member_frac: f64) -> io::Result<Members> {
    let mut out = Members::default();
    let mut counts: Vec<(String, usize)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            let parts: Vec<&str> = h.trim_end().split('|').collect();
            if parts.len() >= 6 && parts[parts.len() - 3] == sample {
                let locus = parts[..parts.len() - 5].join("|");
                if !out.mine.contains(&locus) {
                    out.mine.insert(locus);
                }
            }
            continue;
        }
        if let Some(n) = line.strip_prefix("#colonies") {
            out.n_colonies = n.trim().parse().ok();
            continue;
        }
        if line.is_empty() || line.starts_with('#') || line.starts_with("locus\t") {
            continue;
        }
        let Some((locus, samples)) = line.split_once('\t') else { continue }; // FASTA sequence lines
        let mut n = 0usize;
        let mut me = false;
        for s in samples.trim_end().split(',').filter(|s| !s.is_empty()) {
            n += 1;
            me |= s == sample;
        }
        if me {
            out.mine.insert(locus.to_string());
        }
        counts.push((locus.to_string(), n));
    }
    if let Some(nc) = out.n_colonies {
        out.germline = counts.into_iter().filter(|&(_, n)| is_germline(n, nc, max_member_frac)).map(|(l, _)| l).collect();
    }
    Ok(out)
}

/// Fragment id of a read pair: FNV-1a 64 of the qname, 16 hex digits (insertions.reads.fa.gz
/// uses an opaque 16-hex fragment id too; GT_DISC and its GT_MATE share it).
pub fn frag_id(qname: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in qname {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            x => x,
        })
        .collect()
}

/// A mate's sequence in site-forward (allele) orientation: in an FR pair it is opposite to its
/// anchor, so the stored sequence is kept when the strands differ, else reverse-complemented
/// (discovery `model::mate_site_forward`, combine `allele_forward_seq`). An unmapped mate is
/// stored as sequenced with 0x10 clear, which this handles the same way.
pub fn mate_site_forward(stored: &[u8], mate_reverse: bool, anchor_reverse: bool) -> Vec<u8> {
    if mate_reverse != anchor_reverse {
        stored.to_vec()
    } else {
        revcomp(stored)
    }
}

/// A junction-facing clip that is a poly-A/T tail: >= 10 bp and >= 90% A or >= 90% T.
pub fn is_pure_polya(clip: &[u8]) -> bool {
    if clip.len() < 10 {
        return false;
    }
    let a = clip.iter().filter(|&&b| b == b'A' || b == b'a').count();
    let t = clip.iter().filter(|&&b| b == b'T' || b == b't').count();
    a.max(t) * 10 >= clip.len() * 9
}

/// One FASTA record, insertions.reads.fa.gz layout.
#[allow(clippy::too_many_arguments)]
pub fn push_record(out: &mut String, locus: &str, side: &str, role: &str, sample: &str, frag: &str, r12: u8, seq: &[u8]) {
    use std::fmt::Write;
    let _ = writeln!(out, ">{locus}|{side}|{role}|{sample}|{frag}|{r12}");
    out.push_str(&String::from_utf8_lossy(seq));
    out.push('\n');
}

/// Counters of the extra pass (summed per chunk, printed once).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ExtraStats {
    /// loci of which this colony is not a discovery member and that were called (ok rows)
    pub candidates: usize,
    /// loci this colony did not discover that were skipped as germline (members > max frac)
    pub germline_skipped: usize,
    /// ... of which passed the ALT gate (extra pass run)
    pub gated: usize,
    pub clip: usize,
    pub polya: usize,
    pub disc: usize,
    pub mates: usize,
    /// anchors whose mate was not found at its recorded position (or had no position)
    pub mates_missing: usize,
    /// anchors whose mate was not fetched: per-locus or per-colony cap
    pub mates_capped: usize,
    pub micros: u128,
}

impl ExtraStats {
    pub fn add(&mut self, o: &ExtraStats) {
        self.candidates += o.candidates;
        self.germline_skipped += o.germline_skipped;
        self.gated += o.gated;
        self.clip += o.clip;
        self.polya += o.polya;
        self.disc += o.disc;
        self.mates += o.mates;
        self.mates_missing += o.mates_missing;
        self.mates_capped += o.mates_capped;
        self.micros += o.micros;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn members_from_table_and_reads_fa() {
        let tsv = "locus\tsamples\nchr1:100-115\tS1,S2\nchr2:5-6\tS22\nchr3:1-oneside_1\tS2\n";
        let m = members_from(io::Cursor::new(tsv), "S2", 0.5).unwrap();
        assert!(m.mine.contains("chr1:100-115") && m.mine.contains("chr3:1-oneside_1"));
        assert!(!m.mine.contains("chr2:5-6"), "S22 is not S2");
        assert!(m.n_colonies.is_none() && m.germline.is_empty(), "no #colonies: no germline filter");
        let fa = ">chr1:100-115|LEFT|CLIP|S1|abc|1\nACGT\n>chr1:100-115|LEFT|MATE|S2|abc|2\nACGT\n\
                  >odd|name:7-9|RIGHT|DISC|S2|f|1\nAAAA\n>chr9:1-2|LEFT|CLIP|S3|x|1\nA\n";
        let m = members_from(io::Cursor::new(fa), "S2", 0.5).unwrap();
        let mut v: Vec<&String> = m.mine.iter().collect();
        v.sort();
        assert_eq!(v, vec!["chr1:100-115", "odd|name:7-9"]);
        assert!(m.germline.is_empty());
    }

    #[test]
    fn germline_loci_by_member_fraction() {
        // 4 colonies in the run: 3/4 > 0.5 germline, 2/4 = 0.5 NOT (strict), 1/4 not
        let tsv = "#colonies 4\nlocus\tsamples\nA:1-2\tS1,S2,S3\nB:1-2\tS1,S2\nC:1-2\tS3\nD:1-2\tS1,S2,S3,S4\n";
        let m = members_from(io::Cursor::new(tsv), "S4", 0.5).unwrap();
        assert_eq!(m.n_colonies, Some(4));
        assert_eq!(m.candidate("A:1-2"), Candidate::Germline);
        assert_eq!(m.candidate("B:1-2"), Candidate::Yes);
        assert_eq!(m.candidate("C:1-2"), Candidate::Yes);
        assert_eq!(m.candidate("D:1-2"), Candidate::Member, "membership is reported first");
        assert_eq!(m.candidate("E:1-2"), Candidate::Yes, "a locus nobody discovered");
        // the fraction is the config key; colonies without reads still count in the denominator
        let m = members_from(io::Cursor::new(tsv), "S4", 0.25).unwrap();
        assert_eq!(m.candidate("B:1-2"), Candidate::Germline);
        assert_eq!(m.candidate("C:1-2"), Candidate::Yes);
        assert!(is_germline(51, 100, 0.5) && !is_germline(50, 100, 0.5) && !is_germline(3, 0, 0.5));
    }

    #[test]
    fn orientation_and_ids() {
        // mirrors peartree-discovery model::mate_site_forward_follows_the_fr_pair
        let s = b"AACCGT".to_vec();
        assert_eq!(mate_site_forward(&s, true, false), s);
        assert_eq!(mate_site_forward(&s, false, true), s);
        assert_eq!(mate_site_forward(&s, true, true), b"ACGGTT".to_vec());
        assert_eq!(mate_site_forward(&s, false, false), b"ACGGTT".to_vec());
        assert_eq!(frag_id(b"read1"), frag_id(b"read1"));
        assert_ne!(frag_id(b"read1"), frag_id(b"read2"));
        assert_eq!(frag_id(b"").len(), 16);
        assert!(is_pure_polya(b"AAAAAAAAAAAAAAAAAAAC"));
        assert!(is_pure_polya(b"TTTTTTTTTT"));
        assert!(!is_pure_polya(b"AAAAAAAAA"), "too short");
        assert!(!is_pure_polya(b"AAAAAAAACCAAAAAAAAGG"));
        assert_eq!(sample_from_out("/x/genotypes/PD1_lo0001.txt.gz"), "PD1_lo0001");
        assert_eq!(sample_from_out("S1.txt"), "S1");
        assert_eq!(sidecar_path("g/S1.txt.gz"), "g/S1.txt.gz.extra_reads.fa.gz");
        let mut o = String::new();
        push_record(&mut o, "chr1:1-2", "LEFT", ROLE_MATE, "S1", "00ff", 2, b"ACGT");
        assert_eq!(o, ">chr1:1-2|LEFT|GT_MATE|S1|00ff|2\nACGT\n");
    }
}
