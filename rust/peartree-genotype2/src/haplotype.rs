//! Per-locus haplotype model construction (SPEC "Per-locus haplotype model"). Owner: A.

// until the driver is wired in, parts of this module are unused in the binary
#![allow(dead_code)]

use crate::config::Config;
use crate::contract::ContractSides;
use crate::refseq::{norm_base, RefSeq};
use crate::types::{Anchor, Hyp, JunctionConsensus, Locus, LocusKind, LocusModel, Segment, GENOME_Q};
use std::io;

/// Quality for consensus bases whose quality vector is missing / the wrong length.
const FALLBACK_Q: u8 = 30;

/// Accumulates one segment (sequence + qualities + anchors + junction columns).
#[derive(Default)]
struct SegBuilder {
    seq: Vec<u8>,
    qual: Vec<u8>,
    anchors: Vec<Anchor>,
    cols: Vec<usize>,
}

impl SegBuilder {
    /// Append `genome[ref_pos, ref_pos + seq.len())` and register its anchor.
    fn genome(&mut self, ref_pos: i64, seq: &[u8]) {
        self.anchors.push(Anchor { ref_pos, idx: self.seq.len(), len: seq.len() });
        self.seq.extend_from_slice(seq);
        self.qual.extend(std::iter::repeat(GENOME_Q).take(seq.len()));
    }

    /// Append consensus (inserted) bases with their own qualities.
    fn consensus(&mut self, seq: &[u8], qual: &[u8]) {
        self.seq.extend_from_slice(seq);
        self.qual.extend_from_slice(qual);
    }

    /// Mark the current end of the segment as a junction column.
    fn junction(&mut self) {
        self.cols.push(self.seq.len());
    }

    fn finish(self, hyp: Hyp, label: &'static str) -> Segment {
        debug_assert_eq!(self.seq.len(), self.qual.len());
        Segment { hyp, label, seq: self.seq, qual: self.qual, anchors: self.anchors, junction_cols: self.cols }
    }
}

/// A reference-only segment `genome[lo, hi)` with junction columns at the given reference
/// coordinates (those inside the segment).
fn ref_segment(reference: &mut dyn RefSeq, chr: &str, lo: i64, hi: i64, junctions: &[i64], label: &'static str) -> io::Result<Segment> {
    let g = reference.fetch(chr, lo, hi)?;
    let mut b = SegBuilder::default();
    b.genome(lo, &g);
    for &j in junctions {
        if j >= lo && j <= hi {
            let c = (j - lo) as usize;
            if !b.cols.contains(&c) {
                b.cols.push(c);
            }
        }
    }
    b.cols.sort_unstable();
    Ok(b.finish(Hyp::Ref, label))
}

/// Uppercased consensus insertion (non-ACGT -> N) with a quality vector of matching length.
fn insertion_of(side: &Option<JunctionConsensus>, which: &str) -> Result<(Vec<u8>, Vec<u8>), String> {
    let jc = side.as_ref().ok_or_else(|| format!("missing {which} junction consensus"))?;
    if jc.ins_seq.is_empty() {
        return Err(format!("empty {which} insertion consensus"));
    }
    let seq: Vec<u8> = jc.ins_seq.iter().map(|&b| norm_base(b)).collect();
    let qual = if jc.ins_qual.len() == seq.len() { jc.ins_qual.clone() } else { vec![FALLBACK_Q; seq.len()] };
    Ok((seq, qual))
}

/// True when `overlap` is acceptable evidence of a real shared sequence: at least 3 distinct
/// bases and no single base above 80% (homopolymers / poly-A tails / dinucleotide repeats do
/// not prove that the two junction consensuses are the same molecule).
fn complex_enough(overlap: &[u8]) -> bool {
    let mut counts = [0usize; 4];
    let mut n = 0usize;
    for &b in overlap {
        let i = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => continue,
        };
        counts[i] += 1;
        n += 1;
    }
    if n == 0 {
        return false;
    }
    let distinct = counts.iter().filter(|&&c| c > 0).count();
    let max = *counts.iter().max().unwrap_or(&0);
    distinct >= 3 && (max as f64) <= 0.8 * n as f64
}

/// Longest suffix(`ins_r`) == prefix(`ins_l`) overlap of at least `merge_overlap_min` bases,
/// with at most `merge_overlap_max_mismatch_frac` mismatches (N is a wildcard), that passes the
/// complexity gate. Returns the overlap length.
fn find_overlap(ins_r: &[u8], ins_l: &[u8], cfg: &Config) -> Option<usize> {
    let kmin = cfg.merge_overlap_min.max(1);
    let kmax = ins_r.len().min(ins_l.len());
    if kmax < kmin {
        return None;
    }
    for k in (kmin..=kmax).rev() {
        let suffix = &ins_r[ins_r.len() - k..];
        let prefix = &ins_l[..k];
        let mismatches = suffix.iter().zip(prefix).filter(|(&a, &b)| a != b && a != b'N' && b != b'N').count();
        if mismatches as f64 > cfg.merge_overlap_max_mismatch_frac * k as f64 {
            continue;
        }
        // judge complexity on the bases both consensuses agree on (prefer a non-N base)
        let merged: Vec<u8> = suffix.iter().zip(prefix).map(|(&a, &b)| if a == b'N' { b } else { a }).collect();
        if complex_enough(&merged) {
            return Some(k);
        }
    }
    None
}

/// Merge `ins_r` (5' part, first) and `ins_l` (3' part) over an overlap of `k` bases: the
/// overlap takes, per base, the higher-quality call (equal bases keep the larger quality).
fn merge_insertions(ins_r: &(Vec<u8>, Vec<u8>), ins_l: &(Vec<u8>, Vec<u8>), k: usize) -> (Vec<u8>, Vec<u8>) {
    let (rs, rq) = ins_r;
    let (ls, lq) = ins_l;
    let rstart = rs.len() - k;
    let mut seq = rs.clone();
    let mut qual = rq.clone();
    for i in 0..k {
        let (a, qa) = (rs[rstart + i], rq[rstart + i]);
        let (b, qb) = (ls[i], lq[i]);
        let (base, q) = if a == b {
            (a, qa.max(qb))
        } else if a == b'N' {
            (b, qb)
        } else if b == b'N' {
            (a, qa)
        } else if qa >= qb {
            (a, (qa - qb).max(2))
        } else {
            (b, (qb - qa).max(2))
        };
        seq[rstart + i] = base;
        qual[rstart + i] = q;
    }
    seq.extend_from_slice(&ls[k..]);
    qual.extend_from_slice(&lq[k..]);
    (seq, qual)
}

/// Build the segments / windows / breakpoints for one locus. `sides` = the best available
/// consensus (combined if present, else contract 12-bp). Never panics on bad input: returns
/// `LocusModel::error(..)` for a locus that cannot be modelled (missing consensus, contig not in
/// the reference, ...). I/O errors from the reference propagate.
pub fn build_model(locus: &Locus, sides: &ContractSides, reference: &mut dyn RefSeq, cfg: &Config) -> io::Result<LocusModel> {
    let kind = locus.kind();
    let f = cfg.flank;
    let (l, r) = (locus.left_pos, locus.right_pos);
    let chr = locus.chr.as_str();

    let Some(clen) = reference.contig_len(chr) else {
        return Ok(LocusModel::error(locus.clone(), format!("contig '{chr}' not in the reference")));
    };
    if l.min(r) < 0 || l.max(r) > clen {
        return Ok(LocusModel::error(locus.clone(), format!("locus position outside contig '{chr}' (length {clen})")));
    }

    // ---- one-sided: only the real breakpoint ----
    if locus.is_one_sided() {
        let (bp, which, label, side) = if locus.right_open {
            (l, "left", "ALT_L", &sides.left)
        } else {
            (r, "right", "ALT_R", &sides.right)
        };
        let (iseq, iqual) = match insertion_of(side, which) {
            Ok(x) => x,
            Err(e) => return Ok(LocusModel::error(locus.clone(), e)),
        };
        let ref_seg = ref_segment(reference, chr, bp - f, bp + f, &[bp], "REF")?;
        let mut b = SegBuilder::default();
        if locus.right_open {
            // ALT_L = insL ++ genome[L, L+F)
            b.consensus(&iseq, &iqual);
            b.junction();
            let g = reference.fetch(chr, bp, bp + f)?;
            b.genome(bp, &g);
        } else {
            // ALT_R = genome[R-F, R) ++ insR
            let g = reference.fetch(chr, bp - f, bp)?;
            b.genome(bp - f, &g);
            b.junction();
            b.consensus(&iseq, &iqual);
        }
        return Ok(LocusModel {
            locus: locus.clone(),
            kind,
            windows: vec![(bp, bp)],
            breakpoints: vec![bp],
            segments: vec![ref_seg, b.finish(Hyp::Alt, label)],
            alt_ref_junction_fraction: 0.0,
            alt_full: false,
            error: None,
        });
    }

    // ---- two-sided ----
    let ins_r = match insertion_of(&sides.right, "right") {
        Ok(x) => x,
        Err(e) => return Ok(LocusModel::error(locus.clone(), e)),
    };
    let ins_l = match insertion_of(&sides.left, "left") {
        Ok(x) => x,
        Err(e) => return Ok(LocusModel::error(locus.clone(), e)),
    };

    let gap = r - l;
    let near = gap.abs() <= f;
    let (lo, hi) = (l.min(r), l.max(r));

    let mut segments: Vec<Segment> = Vec::new();
    if near {
        segments.push(ref_segment(reference, chr, lo - f, hi + f, &[l, r], "REF")?);
    } else {
        segments.push(ref_segment(reference, chr, l - f, l + f, &[l], "REF_L")?);
        segments.push(ref_segment(reference, chr, r - f, r + f, &[r], "REF_R")?);
    }

    let flank_r = reference.fetch(chr, r - f, r)?; // genome[R-F, R)
    let flank_l = reference.fetch(chr, l, l + f)?; // genome[L, L+F)

    let overlap = find_overlap(&ins_r.0, &ins_l.0, cfg);
    let alt_full = overlap.is_some();
    if let Some(k) = overlap {
        // ALT_FULL = genome[R-F, R) ++ merge(insR, insL) ++ genome[L, L+F)
        let (mseq, mqual) = merge_insertions(&ins_r, &ins_l, k);
        let mut b = SegBuilder::default();
        b.genome(r - f, &flank_r);
        b.junction();
        b.consensus(&mseq, &mqual);
        b.junction();
        b.genome(l, &flank_l);
        segments.push(b.finish(Hyp::Alt, "ALT_FULL"));
    } else {
        // ALT_R = genome[R-F, R) ++ insR
        let mut b = SegBuilder::default();
        b.genome(r - f, &flank_r);
        b.junction();
        b.consensus(&ins_r.0, &ins_r.1);
        segments.push(b.finish(Hyp::Alt, "ALT_R"));
        // ALT_L = insL ++ genome[L, L+F)
        let mut b = SegBuilder::default();
        b.consensus(&ins_l.0, &ins_l.1);
        b.junction();
        b.genome(l, &flank_l);
        segments.push(b.finish(Hyp::Alt, "ALT_L"));
    }

    let windows = if near { vec![(lo, hi)] } else { vec![(l, l), (r, r)] };
    let breakpoints = if l == r { vec![l] } else { vec![l, r] };
    let alt_ref_junction_fraction = if kind == LocusKind::FarDuplication && gap >= cfg.dup_retained_min_span { 0.5 } else { 0.0 };

    Ok(LocusModel {
        locus: locus.clone(),
        kind,
        windows,
        breakpoints,
        segments,
        alt_ref_junction_fraction,
        alt_full,
        error: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::parse_locus_name;
    use crate::refseq::{rand_seq, MemRef};

    const CHR: &str = "chrT";
    const CLEN: usize = 20_000;

    fn genome() -> MemRef {
        MemRef::new(&[(CHR, &rand_seq(CLEN, 12345))])
    }

    fn jc(ins: &[u8], q: u8) -> JunctionConsensus {
        JunctionConsensus {
            ins_seq: ins.to_vec(),
            ins_qual: vec![q; ins.len()],
            flank_seq: Vec::new(),
            flank_qual: Vec::new(),
            from_contract_only: false,
        }
    }

    fn sides(ins_l: &[u8], ins_r: &[u8]) -> ContractSides {
        ContractSides { left: Some(jc(ins_l, 40)), right: Some(jc(ins_r, 35)) }
    }

    fn model(name: &str, s: &ContractSides, cfg: &Config) -> (LocusModel, MemRef) {
        let mut g = genome();
        let locus = parse_locus_name(name).unwrap();
        let m = build_model(&locus, s, &mut g, cfg).unwrap();
        (m, g)
    }

    fn seg<'a>(m: &'a LocusModel, label: &str) -> &'a Segment {
        m.segments.iter().find(|s| s.label == label).unwrap_or_else(|| panic!("no segment {label} in {:?}", m.segments.iter().map(|s| s.label).collect::<Vec<_>>()))
    }

    /// Segment index of query base 0 of a read the aligner placed at ref position `p` with a
    /// leading soft clip `s` (the formula documented on `Segment::anchors`).
    fn anchor_index(sg: &Segment, p: i64, s: i64) -> Option<i64> {
        sg.anchors.iter().find(|a| a.ref_pos <= p && p < a.ref_pos + a.len as i64).map(|a| a.idx as i64 + (p - a.ref_pos) - s)
    }

    /// every anchor really is genome[ref_pos, ref_pos+len) at segment idx, with GENOME_Q
    fn check_anchors(sg: &Segment, g: &mut MemRef) {
        assert_eq!(sg.seq.len(), sg.qual.len());
        for a in &sg.anchors {
            let want = g.fetch(CHR, a.ref_pos, a.ref_pos + a.len as i64).unwrap();
            assert_eq!(&sg.seq[a.idx..a.idx + a.len], &want[..], "{} anchor {:?}", sg.label, a);
            assert!(sg.qual[a.idx..a.idx + a.len].iter().all(|&q| q == GENOME_Q));
        }
    }

    const INS_R: &[u8] = b"GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTT";
    const INS_L: &[u8] = b"CCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACAAAAAAAAAAAA";

    #[test]
    fn tsd_gap15() {
        let cfg = Config::default();
        let (l, r) = (5000i64, 5015i64);
        let (m, mut g) = model("chrT:5000-5015", &sides(INS_L, INS_R), &cfg);
        assert!(m.error.is_none());
        assert_eq!(m.kind, LocusKind::Tsd);
        assert_eq!(m.windows, vec![(l, r)]);
        assert_eq!(m.breakpoints, vec![l, r]);
        assert_eq!(m.alt_ref_junction_fraction, 0.0);
        assert!(!m.alt_full);
        let labels: Vec<_> = m.segments.iter().map(|s| s.label).collect();
        assert_eq!(labels, vec!["REF", "ALT_R", "ALT_L"]);
        let f = cfg.flank;
        let rf = seg(&m, "REF");
        assert_eq!(rf.hyp, Hyp::Ref);
        assert_eq!(rf.seq.len() as i64, (r - l) + 2 * f);
        assert_eq!(rf.junction_cols, vec![f as usize, (f + 15) as usize]);
        check_anchors(rf, &mut g);
        // ALT_R = genome[R-F, R) ++ insR
        let ar = seg(&m, "ALT_R");
        assert_eq!(ar.hyp, Hyp::Alt);
        assert_eq!(ar.seq.len(), f as usize + INS_R.len());
        assert_eq!(&ar.seq[f as usize..], INS_R);
        assert_eq!(ar.junction_cols, vec![f as usize]);
        assert_eq!(ar.anchors, vec![Anchor { ref_pos: r - f, idx: 0, len: f as usize }]);
        assert!(ar.qual[f as usize..].iter().all(|&q| q == 35));
        check_anchors(ar, &mut g);
        // ALT_L = insL ++ genome[L, L+F)
        let al = seg(&m, "ALT_L");
        assert_eq!(&al.seq[..INS_L.len()], INS_L);
        assert_eq!(al.seq.len(), INS_L.len() + f as usize);
        assert_eq!(al.junction_cols, vec![INS_L.len()]);
        assert_eq!(al.anchors, vec![Anchor { ref_pos: l, idx: INS_L.len(), len: f as usize }]);
        assert!(al.qual[..INS_L.len()].iter().all(|&q| q == 40));
        check_anchors(al, &mut g);
        // a read starting at ref position p (no clip) maps to the expected index
        assert_eq!(anchor_index(rf, l - 10, 0), Some(f - 10));
        assert_eq!(anchor_index(rf, r + 20, 0), Some(f + 15 + 20));
        assert_eq!(anchor_index(ar, r - 50, 0), Some(f - 50)); // flank part
        assert_eq!(anchor_index(ar, r + 1, 0), None); // not in the genome part of ALT_R
        assert_eq!(anchor_index(al, l + 7, 3), Some((INS_L.len() + 7 - 3) as i64));
        assert_eq!(anchor_index(al, l - 1, 0), None);
    }

    #[test]
    fn blunt_and_tsd_deletion() {
        let cfg = Config::default();
        let f = cfg.flank;
        let (m, mut g) = model("chrT:7000-7000", &sides(INS_L, INS_R), &cfg);
        assert_eq!(m.kind, LocusKind::Blunt);
        assert_eq!(m.windows, vec![(7000, 7000)]);
        assert_eq!(m.breakpoints, vec![7000]);
        assert_eq!(seg(&m, "REF").junction_cols, vec![f as usize]);
        assert_eq!(seg(&m, "REF").seq.len() as i64, 2 * f);
        for s in &m.segments {
            check_anchors(s, &mut g);
        }
        let (m, mut g) = model("chrT:7000-7001", &sides(INS_L, INS_R), &cfg);
        assert_eq!(m.kind, LocusKind::Blunt);
        assert_eq!(m.windows, vec![(7000, 7001)]);
        for s in &m.segments {
            check_anchors(s, &mut g);
        }

        // TSD_DELETION: R < L, genome[R, L) is deleted
        let (m, mut g) = model("chrT:8020-8000", &sides(INS_L, INS_R), &cfg);
        assert_eq!(m.kind, LocusKind::TsdDeletion);
        assert_eq!(m.windows, vec![(8000, 8020)]);
        assert_eq!(m.breakpoints, vec![8020, 8000]);
        let rf = seg(&m, "REF");
        assert_eq!(rf.seq.len() as i64, 20 + 2 * f);
        assert_eq!(rf.junction_cols, vec![f as usize, (f + 20) as usize]);
        check_anchors(rf, &mut g);
        let ar = seg(&m, "ALT_R");
        assert_eq!(ar.anchors[0].ref_pos, 8000 - f);
        check_anchors(ar, &mut g);
        let al = seg(&m, "ALT_L");
        assert_eq!(al.anchors[0].ref_pos, 8020);
        check_anchors(al, &mut g);
    }

    #[test]
    fn far_deletion_two_windows() {
        let cfg = Config::default();
        let f = cfg.flank;
        let (m, mut g) = model("chrT:12000-7000", &sides(INS_L, INS_R), &cfg);
        assert_eq!(m.kind, LocusKind::FarDeletion);
        assert_eq!(m.windows, vec![(12000, 12000), (7000, 7000)]);
        assert_eq!(m.breakpoints, vec![12000, 7000]);
        assert_eq!(m.alt_ref_junction_fraction, 0.0);
        let labels: Vec<_> = m.segments.iter().map(|s| s.label).collect();
        assert_eq!(labels, vec!["REF_L", "REF_R", "ALT_R", "ALT_L"]);
        let rl = seg(&m, "REF_L");
        assert_eq!(rl.anchors, vec![Anchor { ref_pos: 12000 - f, idx: 0, len: 2 * f as usize }]);
        assert_eq!(rl.junction_cols, vec![f as usize]);
        let rr = seg(&m, "REF_R");
        assert_eq!(rr.anchors[0].ref_pos, 7000 - f);
        assert_eq!(rr.junction_cols, vec![f as usize]);
        // reads at either breakpoint land in the right REF segment
        assert_eq!(anchor_index(rl, 12000 + 5, 0), Some(f + 5));
        assert_eq!(anchor_index(rr, 7000 - 5, 0), Some(f - 5));
        assert_eq!(anchor_index(rl, 7000, 0), None);
        for s in &m.segments {
            check_anchors(s, &mut g);
        }
    }

    #[test]
    fn far_duplication_fraction() {
        let cfg = Config::default();
        let (m, mut g) = model("chrT:3000-3600", &sides(INS_L, INS_R), &cfg);
        assert_eq!(m.kind, LocusKind::FarDuplication);
        assert_eq!(m.alt_ref_junction_fraction, 0.5);
        assert_eq!(m.windows, vec![(3000, 3000), (3600, 3600)]);
        for s in &m.segments {
            check_anchors(s, &mut g);
        }
        // below dup_retained_min_span (150) but still a far duplication (> F would need a bigger
        // gap; use a lowered flank so gap 100 is far): fraction stays 0
        let mut c2 = Config::default();
        c2.flank = 60;
        let (m, _) = model("chrT:3000-3100", &sides(INS_L, INS_R), &c2);
        assert_eq!(m.kind, LocusKind::FarDuplication);
        assert_eq!(m.alt_ref_junction_fraction, 0.0);
        assert_eq!(m.windows.len(), 2);
        // near duplication/TSD never gets the fraction
        let (m, _) = model("chrT:3000-3030", &sides(INS_L, INS_R), &cfg);
        assert_eq!(m.alt_ref_junction_fraction, 0.0);
    }

    #[test]
    fn one_sided_both_orientations() {
        let cfg = Config::default();
        let f = cfg.flank;
        // L-oneside_L: real junction at L, only the LEFT consensus
        let s = ContractSides { left: Some(jc(INS_L, 40)), right: None };
        let (m, mut g) = model("chrT:9000-oneside_9000", &s, &cfg);
        assert!(m.error.is_none());
        assert_eq!(m.kind, LocusKind::OneSided);
        assert_eq!(m.windows, vec![(9000, 9000)]);
        assert_eq!(m.breakpoints, vec![9000]);
        let labels: Vec<_> = m.segments.iter().map(|s| s.label).collect();
        assert_eq!(labels, vec!["REF", "ALT_L"]);
        let rf = seg(&m, "REF");
        assert_eq!(rf.seq.len() as i64, 2 * f);
        assert_eq!(rf.junction_cols, vec![f as usize]);
        let al = seg(&m, "ALT_L");
        assert_eq!(&al.seq[..INS_L.len()], INS_L);
        assert_eq!(al.anchors, vec![Anchor { ref_pos: 9000, idx: INS_L.len(), len: f as usize }]);
        for s in &m.segments {
            check_anchors(s, &mut g);
        }
        // oneside_R-R: real junction at R, only the RIGHT consensus
        let s = ContractSides { left: None, right: Some(jc(INS_R, 35)) };
        let (m, mut g) = model("chrT:oneside_9500-9500", &s, &cfg);
        assert!(m.error.is_none());
        assert_eq!(m.windows, vec![(9500, 9500)]);
        let labels: Vec<_> = m.segments.iter().map(|s| s.label).collect();
        assert_eq!(labels, vec!["REF", "ALT_R"]);
        let ar = seg(&m, "ALT_R");
        assert_eq!(&ar.seq[f as usize..], INS_R);
        assert_eq!(ar.anchors, vec![Anchor { ref_pos: 9500 - f, idx: 0, len: f as usize }]);
        assert_eq!(ar.junction_cols, vec![f as usize]);
        for s in &m.segments {
            check_anchors(s, &mut g);
        }
        // the wrong side present -> error row, not a panic
        let s = ContractSides { left: None, right: Some(jc(INS_R, 35)) };
        let (m, _) = model("chrT:9000-oneside_9000", &s, &cfg);
        assert!(m.error.is_some());
    }

    #[test]
    fn alt_full_merge() {
        let cfg = Config::default();
        let f = cfg.flank;
        // a 90-bp insertion: insR = its first 60 bp, insL = its last 55 bp, 25 bp shared
        let ins = rand_seq(90, 4242);
        let ins_r = &ins[..60];
        let ins_l = &ins[35..];
        let (m, mut g) = model("chrT:6000-6012", &sides(ins_l, ins_r), &cfg);
        assert!(m.alt_full);
        let labels: Vec<_> = m.segments.iter().map(|s| s.label).collect();
        assert_eq!(labels, vec!["REF", "ALT_FULL"]);
        let af = seg(&m, "ALT_FULL");
        assert_eq!(af.hyp, Hyp::Alt);
        assert_eq!(af.seq.len(), f as usize + 90 + f as usize);
        assert_eq!(&af.seq[f as usize..f as usize + 90], &ins[..]);
        assert_eq!(af.junction_cols, vec![f as usize, f as usize + 90]);
        assert_eq!(
            af.anchors,
            vec![Anchor { ref_pos: 6012 - f, idx: 0, len: f as usize }, Anchor { ref_pos: 6000, idx: f as usize + 90, len: f as usize }]
        );
        check_anchors(af, &mut g);
        // reads placed right of R map into the second anchor, left of L into the first (the two
        // genome anchors of a TSD overlap on [L, R): a read there has BOTH candidate indices)
        assert_eq!(anchor_index(af, 6012 + 4, 0), Some(f + 90 + 12 + 4));
        assert_eq!(anchor_index(af, 6000 - 4, 0), Some(f - 16));
        // overlap with one mismatch (within 10% of 25) still merges; higher quality base wins
        let mut ins_l_mm = ins[35..].to_vec();
        ins_l_mm[3] = if ins_l_mm[3] == b'A' { b'C' } else { b'A' };
        let (m, _) = model("chrT:6000-6012", &sides(&ins_l_mm, ins_r), &cfg);
        assert!(m.alt_full);
        let af = seg(&m, "ALT_FULL");
        // insL has quality 40 > insR 35, so the merged base is insL's, with quality 40-35 = 5
        assert_eq!(af.seq[f as usize + 38], ins_l_mm[3]);
        assert_eq!(af.qual[f as usize + 38], 5);
    }

    #[test]
    fn no_merge_cases() {
        let cfg = Config::default();
        // no overlap at all
        let a = rand_seq(60, 1);
        let b = rand_seq(60, 2);
        let (m, _) = model("chrT:6000-6012", &sides(&b, &a), &cfg);
        assert!(!m.alt_full);
        assert_eq!(m.segments.len(), 3);
        // overlap shorter than merge_overlap_min (20): 15 shared bases
        let ins = rand_seq(100, 9);
        let (m, _) = model("chrT:6000-6012", &sides(&ins[45..], &ins[..60]), &cfg);
        assert!(!m.alt_full);
        // 30-bp poly-A overlap (the classic solo poly-A tail on both sides) must not merge
        let mut r = rand_seq(40, 11);
        r.extend(std::iter::repeat(b'A').take(30));
        let mut l = vec![b'A'; 30];
        l.extend(rand_seq(40, 12));
        let (m, _) = model("chrT:6000-6012", &sides(&l, &r), &cfg);
        assert!(!m.alt_full, "poly-A overlap must not trigger ALT_FULL");
        assert_eq!(m.segments.iter().filter(|s| s.hyp == Hyp::Alt).count(), 2);
        // low complexity: AT repeat (2 distinct bases) overlap
        let at: Vec<u8> = (0..30).map(|i| if i % 2 == 0 { b'A' } else { b'T' }).collect();
        let mut r = rand_seq(40, 21);
        r.extend(&at);
        let mut l = at.clone();
        l.extend(rand_seq(40, 22));
        let (m, _) = model("chrT:6000-6012", &sides(&l, &r), &cfg);
        assert!(!m.alt_full);
        // too many mismatches in an otherwise long overlap
        let ins = rand_seq(100, 33);
        let mut l = ins[40..].to_vec();
        for i in (0..l.len().min(20)).step_by(2) {
            l[i] = if l[i] == b'G' { b'T' } else { b'G' };
        }
        let (m, _) = model("chrT:6000-6012", &sides(&l, &ins[..60]), &cfg);
        assert!(!m.alt_full);
    }

    #[test]
    fn errors_not_panics() {
        let cfg = Config::default();
        let s = sides(INS_L, INS_R);
        let mut g = genome();
        // unknown contig
        let m = build_model(&parse_locus_name("chrX:100-110").unwrap(), &s, &mut g, &cfg).unwrap();
        assert!(m.error.as_deref().unwrap().contains("not in the reference"));
        assert!(m.segments.is_empty());
        // missing consensus on a two-sided locus
        let s1 = ContractSides { left: Some(jc(INS_L, 40)), right: None };
        let m = build_model(&parse_locus_name("chrT:5000-5010").unwrap(), &s1, &mut g, &cfg).unwrap();
        assert!(m.error.as_deref().unwrap().contains("right"));
        // empty insertion
        let s2 = ContractSides { left: Some(jc(b"", 40)), right: Some(jc(INS_R, 35)) };
        assert!(build_model(&parse_locus_name("chrT:5000-5010").unwrap(), &s2, &mut g, &cfg).unwrap().error.is_some());
        // outside the contig
        let m = build_model(&parse_locus_name("chrT:19990-20100").unwrap(), &s, &mut g, &cfg).unwrap();
        assert!(m.error.is_some());
        // near the contig start: flank padded with N, no error
        let m = build_model(&parse_locus_name("chrT:100-110").unwrap(), &s, &mut g, &cfg).unwrap();
        assert!(m.error.is_none());
        let ar = seg(&m, "ALT_R");
        assert_eq!(ar.seq[0], b'N');
        // lowercase / IUPAC consensus is normalised
        let sl = ContractSides { left: Some(jc(b"acgtRyn", 40)), right: Some(jc(b"ttgcaaa", 35)) };
        let m = build_model(&parse_locus_name("chrT:5000-5010").unwrap(), &sl, &mut g, &cfg).unwrap();
        let ar = seg(&m, "ALT_R");
        assert_eq!(&ar.seq[300..], b"TTGCAAA");
        let al = seg(&m, "ALT_L");
        assert_eq!(&al.seq[..7], b"ACGTNNN");
    }

    /// Real-data smoke test (ignored by default):
    ///   PT_E2E=<e2e_phylo dir> cargo test --release smoke_real_data -- --ignored --nocapture
    #[test]
    #[ignore]
    fn smoke_real_data() {
        use std::collections::BTreeMap;
        let Ok(dir) = std::env::var("PT_E2E") else { return };
        let loci = crate::contract::load_loci(&format!("{dir}/genotype/P1.genotyping.tprt.txt.gz")).unwrap();
        let comb = crate::contract::load_combined(&format!("{dir}/combine/P1.combined.txt.gz")).unwrap();
        let mut reference = crate::refseq::open_reference(&format!("{dir}/ref/reduced.fa")).unwrap();
        let cfg = Config::default();
        let hit = loci.iter().filter(|(l, _)| comb.contains_key(&l.name)).count();
        println!("contract loci {}, combined records {}, contract loci covered by combined {hit}", loci.len(), comb.len());
        let mut per_kind: BTreeMap<&str, (usize, usize, usize)> = BTreeMap::new(); // n, errors, merges
        let (mut lens_l, mut lens_r) = (Vec::new(), Vec::new());
        let mut errors: BTreeMap<String, usize> = BTreeMap::new();
        for (locus, contract_sides) in &loci {
            let sides = comb.get(&locus.name).cloned().unwrap_or_else(|| contract_sides.clone());
            let m = build_model(locus, &sides, reference.as_mut(), &cfg).unwrap();
            let e = per_kind.entry(m.kind.as_str()).or_default();
            e.0 += 1;
            if let Some(err) = &m.error {
                e.1 += 1;
                *errors.entry(err.clone()).or_default() += 1;
                continue;
            }
            if m.alt_full {
                e.2 += 1;
            }
            if let Some(x) = &sides.left { lens_l.push(x.ins_seq.len()); }
            if let Some(x) = &sides.right { lens_r.push(x.ins_seq.len()); }
            for sg in &m.segments {
                assert_eq!(sg.seq.len(), sg.qual.len());
                for a in &sg.anchors {
                    let want = reference.fetch(&locus.chr, a.ref_pos, a.ref_pos + a.len as i64).unwrap();
                    assert_eq!(&sg.seq[a.idx..a.idx + a.len], &want[..]);
                }
            }
        }
        for (k, (n, e, mg)) in &per_kind {
            println!("{k:20} n={n:4} errors={e:3} alt_full={mg:3}");
        }
        for (e, n) in &errors {
            println!("error x{n}: {e}");
        }
        let med = |v: &mut Vec<usize>| { v.sort_unstable(); v.get(v.len() / 2).copied() };
        println!("median insL {:?} (n={}), median insR {:?} (n={})", med(&mut lens_l), lens_l.len(), med(&mut lens_r), lens_r.len());
    }
}
