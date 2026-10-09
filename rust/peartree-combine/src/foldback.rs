//! Combine-stage cruciform fold-back filter (`foldback_filter`), the pooled counterpart of
//! discovery's SPEC-9 gate.
//!
//! Low-input enzymatic-fragmentation libraries turn cruciform DNA at inverted repeats into
//! hairpin fragments whose soft clip is an inverted copy of the adjacent reference (Ellis et al.
//! 2021 Nat Protoc, Fig. 3). Discovery drops the clear cases per colony; what survives there has
//! clips too short in one colony, or a hairpin with a mismatch. Combine judges the POOLED clip of
//! each end (the one written to `combined.txt.gz`), allowing `foldback_max_mismatch` mismatches.
//!
//! Per end: the junction-proximal `k` clip bases (a clip of `min_short`..`k`-1 bases whole) are
//! looked up, reverse-complemented, within ±`window` bp of the junction. `Clear` = a full `k`-mer
//! matched; `Maybe` = a shorter clip matched. A locus is dropped when either end is `Clear` or
//! both ends are `Clear`/`Maybe` (Jeremy 2026-10-09: a short possible fold-back only goes with a
//! fold-back partner). Low-complexity probes (homopolymer >= 8, entropy < `min_entropy`) are never
//! judged. PD45886, ±50 bp, <= 1 mismatch: 17 of the 93 private calls that pass the discovery gate,
//! 6/2573 germline MEIs.

use crate::genome::RefFetch;
use crate::model::{Insertion, Interner, TokKind};

const PROBE_MAX_RUN: usize = 8;

#[derive(Clone, Copy, Debug)]
pub struct Params {
    pub k: usize,
    pub min_short: usize,
    pub window: i64,
    pub max_mismatch: usize,
    pub min_entropy: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fold {
    Clear,
    Maybe,
    No,
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

fn entropy(s: &[u8]) -> f64 {
    let mut c = [0usize; 4];
    for b in s {
        match b {
            b'A' => c[0] += 1,
            b'C' => c[1] += 1,
            b'G' => c[2] += 1,
            b'T' => c[3] += 1,
            _ => {}
        }
    }
    let n: usize = c.iter().sum();
    if n == 0 {
        return 0.0;
    }
    c.iter().filter(|&&x| x > 0).map(|&x| x as f64 / n as f64).map(|p| -p * p.log2()).sum()
}

fn longest_run(s: &[u8]) -> usize {
    let (mut best, mut cur, mut prev) = (0, 0, 0u8);
    for &b in s {
        cur = if b == prev { cur + 1 } else { 1 };
        prev = b;
        best = best.max(cur);
    }
    best
}

/// The junction-proximal probe of a REFERENCE-orientation clip (`left` = LEFT end: the clip
/// precedes the junction, so its last bases; RIGHT: its first), uppercased. None when shorter
/// than `min_short`, non-ACGT or low-complexity.
pub fn probe(clip: &[u8], left: bool, fp: &Params) -> Option<Vec<u8>> {
    let n = clip.len().min(fp.k);
    if n == 0 || n < fp.min_short.min(fp.k) {
        return None;
    }
    let s = if left { &clip[clip.len() - n..] } else { &clip[..n] };
    let p: Vec<u8> = s.iter().map(|b| b.to_ascii_uppercase()).collect();
    if !p.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
        return None;
    }
    if entropy(&p) < fp.min_entropy || longest_run(&p) >= PROBE_MAX_RUN {
        return None;
    }
    Some(p)
}

/// Fewest mismatches of `q` against any gapless window of `r` (None when `r` is shorter).
fn best_mismatch(q: &[u8], r: &[u8]) -> Option<usize> {
    if r.len() < q.len() {
        return None;
    }
    r.windows(q.len()).map(|w| w.iter().zip(q).filter(|(a, b)| a != b).count()).min()
}

/// Verdict for one end.
pub fn classify(genome: &dyn RefFetch, contig: &str, clip: &[u8], left: bool, junction: i64, fp: &Params) -> Fold {
    let Some(p) = probe(clip, left, fp) else { return Fold::No };
    let r = genome.fetch(contig, (junction - fp.window).max(0), junction + fp.window + 1);
    match best_mismatch(&revcomp(&p), &r) {
        Some(mm) if mm <= fp.max_mismatch => {
            if p.len() >= fp.k { Fold::Clear } else { Fold::Maybe }
        }
        _ => Fold::No,
    }
}

/// Drop rule on the two end verdicts.
pub fn drops(l: Fold, r: Fold) -> bool {
    l == Fold::Clear || r == Fold::Clear || (l != Fold::No && r != Fold::No)
}

/// True if the insertion is a fold-back artefact. Ends without a clip, or whose name token is not a
/// junction coordinate (`polyA_` / `disc_` / `oneside_`), count as `No`.
pub fn insertion_is_foldback(ins: &Insertion, contigs: &Interner, genome: &dyn RefFetch, fp: &Params) -> bool {
    let contig = contigs.name(ins.contig);
    let end = |clip: Option<&crate::seq::QualSeq>, left: bool| {
        let tok = if left { &ins.name_start } else { &ins.name_end };
        match clip {
            Some(c) if tok.kind == TokKind::Pos => classify(genome, contig, &c.seq, left, tok.pos, fp),
            _ => Fold::No,
        }
    };
    drops(end(ins.left_clipped.as_ref(), true), end(ins.right_clipped.as_ref(), false))
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Fake(Vec<u8>, i64); // sequence starting at offset .1 on contig "12"
    impl RefFetch for Fake {
        fn fetch(&self, seqname: &str, start: i64, end: i64) -> Vec<u8> {
            if seqname != "12" {
                return Vec::new();
            }
            let (s, e) = ((start - self.1).max(0) as usize, ((end - self.1).max(0) as usize).min(self.0.len()));
            if s >= e { Vec::new() } else { self.0[s..e].to_vec() }
        }
    }

    // hg19 chr12:[81477142, 81477242), PD45886 lo0019 (combine L/R junctions at 81477172 / 81477195)
    const SITE: &[u8] = b"TTTTTTTCTTTTGAGACGGAGTCTCACTCTGTCACCCAGCCTGGAGTGCAATGGCATGATCTCTGCTCACTGCAACCTCCATCTCCCCGCTTCAACCATT";
    const FP: Params = Params { k: 20, min_short: 12, window: 50, max_mismatch: 1, min_entropy: 1.0 };

    fn g() -> Fake {
        Fake(SITE.to_vec(), 81477142)
    }

    #[test]
    fn lo0019_ends_in_reference_orientation() {
        // combined.txt.gz: L `tgggtgacagagtgagactcc|GTCACCC...`, R `...GTGCAATG|agcagagatcatgccatt`
        assert_eq!(classify(&g(), "12", b"tgggtgacagagtgagactcc", true, 81477172, &FP), Fold::Clear);
        assert_eq!(classify(&g(), "12", b"agcagagatcatgccatt", false, 81477195, &FP), Fold::Maybe);
        assert!(drops(Fold::Clear, Fold::Maybe));
    }

    #[test]
    fn one_mismatch_hairpin_is_caught_two_are_not() {
        // flip one base of the L clip (an imperfect inverted repeat)
        assert_eq!(classify(&g(), "12", b"tgggtgacagaCtgagactcc", true, 81477172, &FP), Fold::Clear);
        assert_eq!(classify(&g(), "12", b"tgggtgacagaCtgagaTtcc", true, 81477172, &FP), Fold::No);
        let strict = Params { max_mismatch: 0, ..FP };
        assert_eq!(classify(&g(), "12", b"tgggtgacagaCtgagactcc", true, 81477172, &strict), Fold::No);
    }

    #[test]
    fn strand_and_complexity() {
        // the L clip's DIRECT copy of the flank (not inverted) is no fold-back
        assert_eq!(classify(&g(), "12", b"GGAGTCTCACTCTGTCACCC", true, 81477172, &FP), Fold::No);
        assert!(probe(b"aaaaaaaaaaaaaaaaaaaaaaaa", true, &FP).is_none());
        assert!(probe(b"gagactcc", true, &FP).is_none(), "below min_short");
        assert_eq!(classify(&g(), "7", b"tgggtgacagagtgagactcc", true, 81477172, &FP), Fold::No, "contig absent");
    }

    /// Real data, opt-in: PEARTREE_TEST_COMBINED=<P>.combined.txt.gz PEARTREE_TEST_2BIT=<genome>.2bit
    /// PEARTREE_TEST_OUT=<file> -- writes the loci this filter drops (one per line).
    #[test]
    fn real_combined_scan() {
        use std::io::{BufRead, Write};
        let (Ok(c), Ok(tb), Ok(out)) = (
            std::env::var("PEARTREE_TEST_COMBINED"),
            std::env::var("PEARTREE_TEST_2BIT"),
            std::env::var("PEARTREE_TEST_OUT"),
        ) else {
            return;
        };
        let genome = crate::genome::Genome::open(std::path::Path::new(&tb)).unwrap();
        let rd = std::io::BufReader::new(flate2::read::MultiGzDecoder::new(std::fs::File::open(c).unwrap()));
        let mut ends: std::collections::BTreeMap<String, (Fold, Fold)> = Default::default();
        let mut lines = rd.lines();
        while let Some(Ok(h)) = lines.next() {
            let seq = lines.next().unwrap().unwrap();
            lines.next();
            lines.next();
            let name = &h[1..];
            let (locus, side) = name.rsplit_once(':').unwrap();
            let (contig, span) = locus.rsplit_once(':').unwrap();
            let (a, b) = span.split_once('-').unwrap();
            let left = side == "L";
            let Ok(pos) = (if left { a } else { b }).parse::<i64>() else { continue };
            let clip: Vec<u8> = seq.bytes().filter(|b| b.is_ascii_lowercase()).collect();
            let f = classify(&genome, contig, &clip, left, pos, &FP);
            let e = ends.entry(locus.to_string()).or_insert((Fold::No, Fold::No));
            if left { e.0 = f } else { e.1 = f }
        }
        let mut w = std::fs::File::create(out).unwrap();
        for (locus, (l, r)) in &ends {
            if drops(*l, *r) {
                writeln!(w, "{locus}").unwrap();
            }
        }
    }

    #[test]
    fn short_possible_needs_a_foldback_partner() {
        assert!(drops(Fold::Maybe, Fold::Maybe));
        assert!(drops(Fold::No, Fold::Clear));
        assert!(!drops(Fold::Maybe, Fold::No));
        assert!(!drops(Fold::No, Fold::No));
    }
}
