//! Insertion contract parsing and locus genotype summarisation — port of
//! src/genotyping_insertion.py (`Insertion.import_file`, `_side_call`,
//! `summarise_evidence`).

use crate::config::*;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// A candidate insertion locus with its four breakpoint consensus sequences.
/// `left_ref`/`right_ref` are the reference-flank consensuses, `left_clipped`/
/// `right_clipped` the inserted-junction (alt) consensuses. Any may be absent.
#[derive(Clone, Debug)]
pub struct Insertion {
    pub name: String,
    pub chr: String,
    pub left_pos: i64,
    pub right_pos: i64,
    pub left_clipped: Option<Vec<u8>>,
    pub right_clipped: Option<Vec<u8>>,
    pub left_ref: Option<Vec<u8>>,
    pub right_ref: Option<Vec<u8>>,
}

impl Insertion {
    /// Parse a locus name "<contig>:<left>-<right>". Split from the right so that
    /// contigs whose own names contain ':' or '-' parse correctly.
    pub fn new(name: &str) -> io::Result<Insertion> {
        let (chr, pos) = name
            .rsplit_once(':')
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, format!("bad locus name: {name}")))?;
        let (left, right) = pos
            .rsplit_once('-')
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, format!("bad locus name: {name}")))?;
        let left_pos = left
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("bad left pos in {name}")))?;
        let right_pos = right
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("bad right pos in {name}")))?;
        Ok(Insertion {
            name: name.to_string(),
            chr: chr.to_string(),
            left_pos,
            right_pos,
            left_clipped: None,
            right_clipped: None,
            left_ref: None,
            right_ref: None,
        })
    }

    /// Parse a gzipped genotyping contract (the `.genotyping.txt.gz` from
    /// combine_insertions): `>name` starts a locus, `@STATUS` selects which
    /// consensus the following sequence line fills.
    pub fn import_file(path: &str) -> io::Result<Vec<Insertion>> {
        let file = File::open(path)?;
        let reader = BufReader::new(MultiGzDecoder::new(file));
        let mut out: Vec<Insertion> = Vec::new();
        let mut status: Option<String> = None;
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let first = line.as_bytes()[0];
            if first == b'>' {
                out.push(Insertion::new(&line[1..])?);
                status = None;
                continue;
            }
            if first == b'@' {
                status = Some(line[1..].to_string());
                continue;
            }
            let Some(ins) = out.last_mut() else { continue };
            match status.as_deref() {
                Some("RIGHT_INSERTION") => ins.right_clipped = Some(line.as_bytes().to_vec()),
                Some("LEFT_INSERTION") => ins.left_clipped = Some(line.as_bytes().to_vec()),
                Some("RIGHT_REFERENCE") => ins.right_ref = Some(line.as_bytes().to_vec()),
                Some("LEFT_REFERENCE") => ins.left_ref = Some(line.as_bytes().to_vec()),
                _ => {}
            }
        }
        Ok(out)
    }
}

// per-read side classification kinds
const KIND_REF: u8 = 0;
const KIND_ALT: u8 = 1;
const KIND_ART: u8 = 2;

/// Classify one breakpoint side of one read as ref/alt/art. Returns
/// `Some((kind, margin))` or None (side not covered, or an exact tie). Port of
/// `Insertion._side_call`.
fn side_call(side: (i64, i64, i64), art_min: i64) -> Option<(u8, i64)> {
    let (r, a, art) = side;
    if r == 0 && a == 0 && art == 0 {
        return None;
    }
    if art > r && art > a && art >= art_min {
        return Some((KIND_ART, art));
    }
    if a > r {
        return Some((KIND_ALT, a - r));
    }
    if r > a {
        return Some((KIND_REF, r - a));
    }
    None
}

/// A finished locus call plus the raw vote counts, ready for the output row.
pub struct Call {
    pub genotype: &'static str,
    pub score_gt: i64,
    pub score_other: i64,
    pub n_ref: i64,
    pub n_alt: i64,
    pub n_art: i64,
}

impl Call {
    fn simple(genotype: &'static str, score_gt: i64, score_other: i64) -> Call {
        Call { genotype, score_gt, score_other, n_ref: 0, n_alt: 0, n_art: 0 }
    }
}

/// Genotype a locus from its per-read (left_side, right_side) evidence tuples.
/// Port of `Insertion.summarise_evidence`: each spanning read casts one vote and
/// the call is driven by the variant allele fraction, with the same low-coverage
/// recovery and zygosity-confidence gates.
pub fn summarise_evidence(reads: &[((i64, i64, i64), (i64, i64, i64))], cfg: &GenotypingConfig) -> Call {
    let art_min = cfg.art_min_score;

    if reads.is_empty() {
        return Call::simple(GT_NO_COVERAGE, 0, 0);
    }

    let (mut n_ref, mut n_alt, mut n_art) = (0i64, 0i64, 0i64);
    let (mut ref_score, mut alt_score, mut art_score) = (0i64, 0i64, 0i64);

    for &(left, right) in reads {
        let mut covered: Vec<(u8, i64)> = Vec::with_capacity(2);
        if let Some(s) = side_call(left, art_min) {
            covered.push(s);
        }
        if let Some(s) = side_call(right, art_min) {
            covered.push(s);
        }
        if covered.is_empty() {
            continue;
        }
        let has_art = covered.iter().any(|&(k, _)| k == KIND_ART);
        let has_alt = covered.iter().any(|&(k, _)| k == KIND_ALT);
        if has_art {
            n_art += 1;
            art_score += covered.iter().filter(|&&(k, _)| k == KIND_ART).map(|&(_, m)| m).sum::<i64>();
            continue;
        }
        // both junctions of a single short read match the inserted element: chimeric.
        let all_alt = covered.iter().all(|&(k, _)| k == KIND_ALT);
        if cfg.double_alt_is_artefact && all_alt && covered.len() == 2 {
            n_art += 1;
            art_score += covered.iter().map(|&(_, m)| m).sum::<i64>();
            continue;
        }
        if has_alt {
            n_alt += 1;
            alt_score += covered.iter().filter(|&&(k, _)| k == KIND_ALT).map(|&(_, m)| m).sum::<i64>();
        } else {
            n_ref += 1;
            ref_score += covered.iter().filter(|&&(k, _)| k == KIND_REF).map(|&(_, m)| m).sum::<i64>();
        }
    }

    let informative = n_ref + n_alt;
    let total = informative + n_art;

    if total == 0 {
        return Call::simple(GT_NO_COVERAGE, 0, 0);
    }
    // artefact-dominated locus.
    if n_art >= cfg.min_artefact_reads && (n_art as f64) >= cfg.artefact_read_fraction * (total as f64) {
        return Call { genotype: GT_ARTEFACT, score_gt: art_score, score_other: ref_score.max(alt_score), n_ref, n_alt, n_art };
    }
    if informative == 0 {
        return Call { genotype: GT_NO_COVERAGE, score_gt: 0, score_other: 0, n_ref, n_alt, n_art };
    }

    let vaf = n_alt as f64 / informative as f64;
    let confident_ins = n_alt >= cfg.min_supporting_reads && alt_score >= cfg.min_score_for_call;
    let recovered_ins = cfg.recover_low_coverage_presence
        && !confident_ins
        && n_alt >= 1
        && alt_score >= cfg.min_score_for_call
        && n_ref < cfg.min_supporting_reads;

    let mk = |g: &'static str, a: i64, b: i64| Call { genotype: g, score_gt: a, score_other: b, n_ref, n_alt, n_art };

    if vaf >= cfg.vaf_hom_min {
        if confident_ins && informative >= cfg.min_reads_for_zygosity {
            return mk(GT_HOMOZYGOUS, alt_score, ref_score);
        }
        if confident_ins || recovered_ins {
            return mk(GT_INSERTION, alt_score, ref_score);
        }
        return mk(GT_INSERTION_UNCERTAIN, alt_score, ref_score);
    }
    if vaf >= cfg.vaf_het_min {
        if confident_ins {
            return mk(GT_HETEROZYGOUS, alt_score, ref_score);
        }
        if recovered_ins {
            return mk(GT_INSERTION, alt_score, ref_score);
        }
        return mk(GT_INSERTION_UNCERTAIN, alt_score, ref_score);
    }
    if vaf > cfg.vaf_wildtype_max {
        return mk(GT_WILDTYPE_UNCERTAIN, ref_score, alt_score);
    }
    let confident_wt = n_ref >= cfg.min_supporting_reads && ref_score >= cfg.min_score_for_call;
    let gt = if confident_wt { GT_WILDTYPE } else { GT_WILDTYPE_UNCERTAIN };
    mk(gt, ref_score, alt_score)
}

#[cfg(test)]
mod tests {
    use super::*;

    const REF_SIDE: (i64, i64, i64) = (300, -100, 0);
    const ALT_SIDE: (i64, i64, i64) = (-100, 300, 0);
    const ART_SIDE: (i64, i64, i64) = (-50, -50, 90);
    const NONE_SIDE: (i64, i64, i64) = (0, 0, 0);

    fn cfg() -> GenotypingConfig {
        GenotypingConfig::default()
    }

    #[test]
    fn heterozygous_balanced() {
        // 3 insertion reads (alt/ref) + 3 wild-type reads -> het band, confident.
        let mut reads = vec![(ALT_SIDE, REF_SIDE); 3];
        reads.extend(vec![(REF_SIDE, REF_SIDE); 3]);
        let c = summarise_evidence(&reads, &cfg());
        assert_eq!(c.genotype, GT_HETEROZYGOUS);
        assert_eq!(c.n_alt, 3);
        assert_eq!(c.n_ref, 3);
    }

    #[test]
    fn insertion_present_zygosity_unclear() {
        // 3 alt-only reads: VAF 1.0 but < min_reads_for_zygosity informative -> insertion.
        let reads = vec![(ALT_SIDE, NONE_SIDE); 3];
        let c = summarise_evidence(&reads, &cfg());
        assert_eq!(c.genotype, GT_INSERTION);
    }

    #[test]
    fn homozygous_needs_enough_reads() {
        let reads = vec![(ALT_SIDE, NONE_SIDE); 6];
        let c = summarise_evidence(&reads, &cfg());
        assert_eq!(c.genotype, GT_HOMOZYGOUS);
    }

    #[test]
    fn single_alt_read_recovered() {
        let reads = vec![(ALT_SIDE, NONE_SIDE); 1];
        let c = summarise_evidence(&reads, &cfg());
        assert_eq!(c.genotype, GT_INSERTION);
    }

    #[test]
    fn both_ends_alt_is_artefact() {
        // a single read whose BOTH junctions match the element -> chimeric artefact.
        let reads = vec![(ALT_SIDE, ALT_SIDE); 3];
        let c = summarise_evidence(&reads, &cfg());
        assert_eq!(c.genotype, GT_ARTEFACT);
    }

    #[test]
    fn artefact_dominated() {
        let reads = vec![(ART_SIDE, ART_SIDE); 4];
        let c = summarise_evidence(&reads, &cfg());
        assert_eq!(c.genotype, GT_ARTEFACT);
        assert_eq!(c.n_art, 4);
    }

    #[test]
    fn clean_wildtype() {
        let reads = vec![(REF_SIDE, REF_SIDE); 5];
        let c = summarise_evidence(&reads, &cfg());
        assert_eq!(c.genotype, GT_WILDTYPE);
    }

    #[test]
    fn no_coverage_empty() {
        let c = summarise_evidence(&[], &cfg());
        assert_eq!(c.genotype, GT_NO_COVERAGE);
    }

    #[test]
    fn parse_locus_name_with_colons() {
        let i = Insertion::new("HLA-A*01:01:01:32992169-32992177").unwrap();
        assert_eq!(i.chr, "HLA-A*01:01:01");
        assert_eq!(i.left_pos, 32992169);
        assert_eq!(i.right_pos, 32992177);
    }
}
