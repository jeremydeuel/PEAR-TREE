//! Contract + combined-consensus parsing. Owner: A.
//!
//! * `load_loci`: `<patient>.genotyping[.tprt].txt.gz` -> ordered loci (+ the 12-bp fallback
//!   consensus per side).
//! * `load_combined`: `<patient>.combined.txt.gz` -> per locus, per side `JunctionConsensus`
//!   (case-split: uppercase = reference flank, lowercase = inserted part; qualities phred+33).
//!
//! ## ONE-SIDED loci in the combined file (checked on e2e_phylo/combine/P1.combined.txt.gz)
//!
//! combine_insertions writes a one-sided locus under its OWN one-sided name with only the real
//! side's record: `contig:L-oneside_L` has just `@contig:L-oneside_L:L`, and
//! `contig:oneside_R-R` has just `@contig:oneside_R-R:R` (e2e_phylo: 19 + 22 such records, none
//! with the other side; a handful of one-sided names are in the combined file but were dropped
//! from the contract, e.g. by combine's per-side exclusions). The contract
//! (`genotyping.tprt.txt.gz`, built by src/genotyping_contract_oneside.py) uses exactly the same
//! names (`>chr22:22776772-oneside_22776772`, `>chr22:oneside_24003149-24003149`), so the lookup
//! `combined[locus.name]` hits directly and returns a `ContractSides` with only the real side
//! `Some` (`left` for `L-oneside_L`, `right` for `oneside_R-R`). No re-keying is needed; the
//! record name minus its trailing `:L`/`:R` is the key.
//!
//! Contract-only flank sequences: `@LEFT_REFERENCE` = `genome[L - n, L)` and `@RIGHT_REFERENCE` =
//! `genome[R, R + n)` (n = 12), i.e. the reference sequence just OUTSIDE each breakpoint, not
//! the reference-aligned part of the combined record. The haplotype builder never uses
//! `flank_*` (it takes the flanks from the genome), only `ins_*`.

// until the driver is wired in, parts of this module are unused in the binary
#![allow(dead_code)]

use crate::refseq::norm_base;
use crate::types::{JunctionConsensus, Locus};
use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// Token discovery writes on the missing end of a one-sided locus (`oneside_<pos>`).
pub const ONESIDE_TOKEN: &str = "oneside_";

/// Quality assigned to contract-only (12-bp) consensus bases.
pub const CONTRACT_Q: u8 = 30;

/// A contract entry: the locus plus whatever consensus the contract itself carries (12 bp).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ContractSides {
    pub left: Option<JunctionConsensus>,
    pub right: Option<JunctionConsensus>,
}

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

fn open_gz(path: &str) -> io::Result<BufReader<MultiGzDecoder<File>>> {
    let f = File::open(path).map_err(|e| io::Error::new(e.kind(), format!("cannot open {path}: {e}")))?;
    Ok(BufReader::new(MultiGzDecoder::new(f)))
}

fn clean(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| norm_base(b)).collect()
}

fn contract_side(ins: &Option<Vec<u8>>, flank: &Option<Vec<u8>>) -> Option<JunctionConsensus> {
    let ins = ins.as_ref()?;
    let ins_seq = clean(ins);
    let flank_seq = flank.as_ref().map(|f| clean(f)).unwrap_or_default();
    Some(JunctionConsensus {
        ins_qual: vec![CONTRACT_Q; ins_seq.len()],
        flank_qual: vec![CONTRACT_Q; flank_seq.len()],
        ins_seq,
        flank_seq,
        from_contract_only: true,
    })
}

#[derive(Default)]
struct RawEntry {
    left_ins: Option<Vec<u8>>,
    left_ref: Option<Vec<u8>>,
    right_ins: Option<Vec<u8>>,
    right_ref: Option<Vec<u8>>,
}

pub fn load_loci(path: &str) -> io::Result<Vec<(Locus, ContractSides)>> {
    let reader = open_gz(path)?;
    let mut raw: Vec<(Locus, RawEntry)> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut dup_skip = false; // true while inside a duplicate (ignored) entry
    let mut n_dup = 0usize;
    let mut status: Option<String> = None;
    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| io::Error::new(e.kind(), format!("{path}: {e}")))?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        match line.as_bytes()[0] {
            b'>' => {
                let name = line[1..].trim();
                let locus = parse_locus_name(name).map_err(|e| invalid(format!("{path}:{}: {e}", i + 1)))?;
                status = None;
                if seen.insert(name.to_string()) {
                    dup_skip = false;
                    raw.push((locus, RawEntry::default()));
                } else {
                    dup_skip = true;
                    n_dup += 1;
                }
            }
            b'@' => status = Some(line[1..].trim().to_string()),
            _ => {
                if dup_skip {
                    continue;
                }
                let Some((_, e)) = raw.last_mut() else { continue };
                let seq = Some(line.as_bytes().to_vec());
                match status.as_deref() {
                    Some("RIGHT_INSERTION") => e.right_ins = seq,
                    Some("RIGHT_REFERENCE") => e.right_ref = seq,
                    Some("LEFT_INSERTION") => e.left_ins = seq,
                    Some("LEFT_REFERENCE") => e.left_ref = seq,
                    _ => {}
                }
            }
        }
    }
    if n_dup > 0 {
        eprintln!("warning: {path}: {n_dup} duplicate locus entries ignored (first one wins)");
    }
    Ok(raw
        .into_iter()
        .map(|(locus, e)| {
            let mut sides = ContractSides {
                left: contract_side(&e.left_ins, &e.left_ref),
                right: contract_side(&e.right_ins, &e.right_ref),
            };
            // an open end carries no junction whatever the file says
            if locus.left_open {
                sides.left = None;
            }
            if locus.right_open {
                sides.right = None;
            }
            (locus, sides)
        })
        .collect())
}

/// Per locus name: (left, right) consensus from the combined file. Missing side -> None.
/// Keys are the record names without the trailing `:L` / `:R` (one-sided names included, see the
/// module comment). Case-split: lowercase = inserted part, uppercase = reference flank; both
/// are returned uppercase with the qualities split alongside.
pub fn load_combined(path: &str) -> io::Result<HashMap<String, ContractSides>> {
    let reader = open_gz(path)?;
    let mut lines = reader.lines();
    let mut out: HashMap<String, ContractSides> = HashMap::new();
    let mut lineno = 0usize;
    loop {
        // header: skip blank lines between records
        let header = loop {
            match next_line(&mut lines, &mut lineno)? {
                None => return Ok(out),
                Some(l) if l.trim().is_empty() => continue,
                Some(l) => break l,
            }
        };
        let at = lineno;
        let trunc = || invalid(format!("{path}: truncated record starting at line {at}"));
        let seq = next_line(&mut lines, &mut lineno)?.ok_or_else(trunc)?;
        let plus = next_line(&mut lines, &mut lineno)?.ok_or_else(trunc)?;
        let qual = next_line(&mut lines, &mut lineno)?.ok_or_else(trunc)?;
        if !header.starts_with('@') || !plus.starts_with('+') {
            return Err(invalid(format!("{path}: malformed FASTQ record at line {at}")));
        }
        let title = header[1..].split_whitespace().next().unwrap_or("");
        let (name, side) = title
            .rsplit_once(':')
            .ok_or_else(|| invalid(format!("{path}: line {at}: record name '{title}' has no :L/:R suffix")))?;
        if seq.len() != qual.len() {
            return Err(invalid(format!("{path}: line {at}: sequence/quality length mismatch for {title}")));
        }
        let jc = split_case(seq.as_bytes(), qual.as_bytes());
        let entry = out.entry(name.to_string()).or_default();
        match side {
            "L" => entry.left = Some(jc),
            "R" => entry.right = Some(jc),
            _ => return Err(invalid(format!("{path}: line {at}: record '{title}' must end in :L or :R"))),
        }
    }
}

/// Split a combined-file sequence by case. Uppercase bases (+ their qualities) -> flank,
/// lowercase -> inserted part. Quality bytes are phred+33.
type GzLines = std::io::Lines<BufReader<MultiGzDecoder<File>>>;

fn next_line(lines: &mut GzLines, lineno: &mut usize) -> io::Result<Option<String>> {
    match lines.next() {
        None => Ok(None),
        Some(l) => {
            *lineno += 1;
            Ok(Some(l?.trim_end_matches(['\r', '\n']).to_string()))
        }
    }
}

fn split_case(seq: &[u8], qual: &[u8]) -> JunctionConsensus {
    let mut jc = JunctionConsensus { from_contract_only: false, ..Default::default() };
    for (&b, &q) in seq.iter().zip(qual) {
        let phred = q.saturating_sub(33);
        if b.is_ascii_lowercase() {
            jc.ins_seq.push(norm_base(b));
            jc.ins_qual.push(phred);
        } else {
            jc.flank_seq.push(norm_base(b));
            jc.flank_qual.push(phred);
        }
    }
    jc
}

fn parse_end(tok: &str, name: &str, which: &str) -> io::Result<(i64, bool)> {
    let (num, open) = match tok.strip_prefix(ONESIDE_TOKEN) {
        Some(rest) => (rest, true),
        None => (tok, false),
    };
    if num.is_empty() || !num.bytes().all(|b| b.is_ascii_digit()) {
        return Err(invalid(format!("bad {which} position '{tok}' in locus name {name}")));
    }
    let pos = num.parse().map_err(|_| invalid(format!("bad {which} position '{tok}' in locus name {name}")))?;
    Ok((pos, open))
}

/// Parse a locus name `<contig>:<L>-<R>` (either end may carry `oneside_`); split from the right
/// so contig names containing ':' or '-' work. For a one-sided locus both positions are set to
/// the real breakpoint.
pub fn parse_locus_name(name: &str) -> io::Result<Locus> {
    let (chr, pos) = name.rsplit_once(':').ok_or_else(|| invalid(format!("bad locus name: {name}")))?;
    let (left, right) = pos.rsplit_once('-').ok_or_else(|| invalid(format!("bad locus name: {name}")))?;
    if chr.is_empty() {
        return Err(invalid(format!("empty contig in locus name: {name}")));
    }
    let (mut left_pos, left_open) = parse_end(left, name, "left")?;
    let (mut right_pos, right_open) = parse_end(right, name, "right")?;
    if left_open && right_open {
        return Err(invalid(format!("both ends open in locus name {name}")));
    }
    if right_open {
        right_pos = left_pos;
    } else if left_open {
        left_pos = right_pos;
    }
    Ok(Locus { name: name.to_string(), chr: chr.to_string(), left_pos, right_pos, left_open, right_open })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::LocusKind;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    fn gz_file(tag: &str, text: &str) -> String {
        let p = std::env::temp_dir().join(format!("pt_g2_contract_{tag}_{}.gz", std::process::id()));
        let mut e = GzEncoder::new(File::create(&p).unwrap(), Compression::fast());
        e.write_all(text.as_bytes()).unwrap();
        e.finish().unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn parse_names() {
        let l = parse_locus_name("chr22:18731099-18731120").unwrap();
        assert_eq!((l.chr.as_str(), l.left_pos, l.right_pos, l.left_open, l.right_open), ("chr22", 18731099, 18731120, false, false));
        assert_eq!(l.kind(), LocusKind::Tsd);
        // HLA-style contig with ':' and '-'
        let l = parse_locus_name("HLA-A*01:01:01:01:100-200").unwrap();
        assert_eq!((l.chr.as_str(), l.left_pos, l.right_pos), ("HLA-A*01:01:01:01", 100, 200));
        let l = parse_locus_name("chrUn_KI270742v1:oneside_500-500").unwrap();
        assert_eq!(l.chr, "chrUn_KI270742v1");
        // one-sided, both orientations
        let l = parse_locus_name("chr22:22776772-oneside_22776772").unwrap();
        assert!(l.right_open && !l.left_open && l.is_one_sided());
        assert_eq!((l.left_pos, l.right_pos), (22776772, 22776772));
        assert_eq!(l.kind(), LocusKind::OneSided);
        let l = parse_locus_name("HLA-DRB1*15:01:01:01:oneside_321-321").unwrap();
        assert!(l.left_open && !l.right_open);
        assert_eq!((l.chr.as_str(), l.left_pos, l.right_pos), ("HLA-DRB1*15:01:01:01", 321, 321));
        // rejections
        for bad in ["chr1", "chr1:100", "chr1:oneside_1-oneside_2", "chr1:abc-5", "chr1:5-x7", "chr1:-5", "chr1:5-", ":5-6", "chr1:oneside_-5", "chr1:1-oneside_"] {
            assert!(parse_locus_name(bad).is_err(), "{bad} should be rejected");
        }
    }

    #[test]
    fn kinds() {
        let k = |n: &str| parse_locus_name(n).unwrap().kind();
        assert_eq!(k("c:1000-1015"), LocusKind::Tsd);
        assert_eq!(k("c:1000-1040"), LocusKind::Tsd);
        assert_eq!(k("c:1000-1041"), LocusKind::FarDuplication);
        assert_eq!(k("c:1000-1000"), LocusKind::Blunt);
        assert_eq!(k("c:1000-1001"), LocusKind::Blunt);
        assert_eq!(k("c:1000-995"), LocusKind::TsdDeletion);
        assert_eq!(k("c:1000-970"), LocusKind::TsdDeletion);
        assert_eq!(k("c:1000-969"), LocusKind::FarDeletion);
        assert_eq!(k("c:10000-5000"), LocusKind::FarDeletion);
        assert_eq!(k("c:1000-1600"), LocusKind::FarDuplication);
    }

    #[test]
    fn contract_loading() {
        let text = "\
>chr22:18731099-18731120
@RIGHT_INSERTION
ATTATATGACAC
@RIGHT_REFERENCE
GTAATATAACAT
@LEFT_INSERTION
TATGACACATAA
@LEFT_REFERENCE
TATAATATTTAC

>chr22:22776772-oneside_22776772
@LEFT_INSERTION
AAAAAAAAAAAA
@LEFT_REFERENCE
CCTATCAGATAT
>chr22:oneside_24003149-24003149
@RIGHT_INSERTION
gcactccagcaR
@RIGHT_REFERENCE
GTAATATAACAT
>chr22:18731099-18731120
@LEFT_INSERTION
CCCCCCCCCCCC
";
        let p = gz_file("loci", text);
        let loci = load_loci(&p).unwrap();
        assert_eq!(loci.len(), 3, "duplicate entry ignored");
        let (l0, s0) = &loci[0];
        assert_eq!(l0.name, "chr22:18731099-18731120");
        let (r, l) = (s0.right.as_ref().unwrap(), s0.left.as_ref().unwrap());
        assert_eq!(r.ins_seq, b"ATTATATGACAC");
        assert_eq!(r.flank_seq, b"GTAATATAACAT");
        assert_eq!(l.ins_seq, b"TATGACACATAA", "first duplicate wins");
        assert!(r.from_contract_only && l.from_contract_only);
        assert!(r.ins_qual.iter().all(|&q| q == 30) && r.ins_qual.len() == 12);
        let (l1, s1) = &loci[1];
        assert!(l1.right_open && s1.right.is_none() && s1.left.as_ref().unwrap().ins_seq == b"AAAAAAAAAAAA");
        let (l2, s2) = &loci[2];
        assert!(l2.left_open && s2.left.is_none());
        assert_eq!(s2.right.as_ref().unwrap().ins_seq, b"GCACTCCAGCAN", "uppercased, non-ACGT -> N");
        assert!(load_loci("/nonexistent/contract.txt.gz").is_err());
    }

    #[test]
    fn combined_parsing() {
        // :L = lowercase ins then UPPER flank; :R = UPPER flank then lowercase ins; qualities
        // are phred+33 ('k' = 74, 'F' = 37, '~' = 93); '@' as first quality char must not be
        // taken for a header; a blank line between records; a contig name containing ':'.
        let text = "\
@HLA-A*01:01:01:01:100-120:L
ttgcaACGTAC
+
kkkkk~~~~~~

@HLA-A*01:01:01:01:100-120:R
ACGTACgggaaa
+
@@@@@@FFFFFF
@chr22:22776772-oneside_22776772:L
aaaaaaaCCGG
+
IIIIIII####
";
        let p = gz_file("comb", text);
        let m = load_combined(&p).unwrap();
        assert_eq!(m.len(), 2);
        let e = &m["HLA-A*01:01:01:01:100-120"];
        let l = e.left.as_ref().unwrap();
        assert_eq!(l.ins_seq, b"TTGCA");
        assert_eq!(l.flank_seq, b"ACGTAC");
        assert_eq!(l.ins_qual, vec![74; 5]);
        assert_eq!(l.flank_qual, vec![93; 6]);
        assert!(!l.from_contract_only);
        let r = e.right.as_ref().unwrap();
        assert_eq!(r.flank_seq, b"ACGTAC");
        assert_eq!(r.flank_qual, vec![31; 6]);
        assert_eq!(r.ins_seq, b"GGGAAA");
        assert_eq!(r.ins_qual, vec![37; 6]);
        // one-sided: keyed by the one-sided name, only the real side present
        let o = &m["chr22:22776772-oneside_22776772"];
        assert!(o.right.is_none());
        assert_eq!(o.left.as_ref().unwrap().ins_seq, b"AAAAAAA");
        // malformed files error instead of panicking
        let p = gz_file("trunc", "@a:1-2:L\nACGT\n+\n");
        assert!(load_combined(&p).is_err());
        let p = gz_file("badq", "@a:1-2:L\nACGT\n+\nIII\n");
        assert!(load_combined(&p).is_err());
        let p = gz_file("badside", "@a:1-2:X\nACGT\n+\nIIII\n");
        assert!(load_combined(&p).is_err());
    }
}
