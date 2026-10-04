//! Genotyping driver — port of src/genotype.py + Insertion.genotype.
//!
//! For each locus: count spanning depth (the high-coverage gate), then fetch the
//! spanning reads, apply the same read gates and per-fragment dedup as the Python
//! genotyper, score each read on both junctions, and summarise to a call. Loci are
//! independent, so the work is split across threads (each with its own indexed
//! reader) and the rows are re-assembled in contract order — the output matches
//! the single-threaded Python run byte-for-byte, with ONE deliberate exception:
//! at high-coverage loci the depth count early-exits (see `genotype_one`), so the
//! reported `coverage` is capped at the threshold rather than the full pileup total.

use crate::config::*;
use crate::evidence::{qleft, qright};
use crate::insertion::{summarise_evidence_opts, Insertion, RefAdjust};
use crate::read::{decode_light, name, qual, reference_name, seq};
use crate::source::{open_source, RegionSource};

use noodles_core::{Position, Region};
use noodles_sam::Header;
use rustc_hash::FxHashSet;
use std::io::{self, Write};
use std::time::Instant;

/// Loci between stderr progress heartbeats (and the flush cadence for the streamed,
/// single-threaded output). Small enough that a hang localises to a narrow window.
const HEARTBEAT_EVERY: usize = 1000;

/// on-disk header + row format, identical to src/genotype.py (OUTPUT_HEADER / format_row).
pub const OUTPUT_HEADER: &str =
    "insertion\tgenotype\tscore_genotype\tscore_alternative\tcoverage\tn_alt\tn_ref\tn_art\n";

struct Row {
    genotype: &'static str,
    score_gt: i64,
    score_other: i64,
    coverage: i64,
    n_alt: i64,
    n_ref: i64,
    n_art: i64,
}

fn format_row(name: &str, r: &Row) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        name, r.genotype, r.score_gt, r.score_other, r.coverage, r.n_alt, r.n_ref, r.n_art
    )
}

/// Build a noodles region for the pysam-style 0-based half-open interval
/// [start0, end0): 1-based inclusive [start0+1, end0].
fn region(chr: &[u8], start0: i64, end0: i64) -> Region {
    let s = Position::new((start0 + 1).max(1) as usize).unwrap();
    let e = Position::new(end0.max(1) as usize).unwrap();
    Region::new(chr.to_vec(), s..=e)
}

#[inline]
fn overlaps(ref_start: i64, ref_end: i64, start0: i64, end0: i64) -> bool {
    ref_start < end0 && ref_end > start0
}

/// Genotype `insertions` against `path` in order, invoking `emit` for each formatted
/// output row. Opens its own indexed reader (thread-local); `reference` is required for
/// CRAM. Emits a stderr heartbeat every `HEARTBEAT_EVERY` loci (and at the end) naming the
/// last locus reached, so a live `tail -f` of the job's stderr tracks progress and a hang
/// localises to the <=HEARTBEAT_EVERY loci after the last heartbeat. `tag` labels the
/// worker in multi-thread runs (empty for the single-threaded stream).
fn process_chunk_with<F: FnMut(String) -> io::Result<()>>(
    path: &str,
    reference: Option<&str>,
    insertions: &[Insertion],
    cfg: &GenotypingConfig,
    tag: &str,
    mut emit: F,
) -> io::Result<()> {
    let mut src = open_source(path, reference)?;
    // Clone the header once so decode can borrow it while `src` is borrowed mutably
    // by the region queries (the two would otherwise conflict).
    let header = src.header().clone();
    let total = insertions.len();
    let t0 = Instant::now();

    for (k, ins) in insertions.iter().enumerate() {
        // Blanket per-locus error handling mirrors the Python worker's try/except:
        // one bad locus becomes an `error` row, never a crash.
        let row = match genotype_one(src.as_mut(), &header, ins, cfg) {
            Ok(row) => row,
            Err(e) => {
                eprintln!("genotyping failed for {}: {e}", ins.name);
                Row { genotype: GT_ERROR, score_gt: 0, score_other: 0, coverage: 0, n_alt: 0, n_ref: 0, n_art: 0 }
            }
        };
        emit(format_row(&ins.name, &row))?;
        if (k + 1) % HEARTBEAT_EVERY == 0 || k + 1 == total {
            eprintln!("  {tag}{}/{} loci, {}s (last {})", k + 1, total, t0.elapsed().as_secs(), ins.name);
        }
    }
    Ok(())
}

/// Collecting wrapper for the multi-thread path: genotype a chunk and return its rows in
/// contract order.
fn process_chunk(
    path: &str,
    reference: Option<&str>,
    insertions: &[Insertion],
    cfg: &GenotypingConfig,
    tag: &str,
) -> io::Result<Vec<String>> {
    let mut rows = Vec::with_capacity(insertions.len());
    process_chunk_with(path, reference, insertions, cfg, tag, |row| {
        rows.push(row);
        Ok(())
    })?;
    Ok(rows)
}

/// Breakpoint windows `(lo, hi)` (0-based, inclusive breakpoint coordinates) a locus is
/// depth-gated and fetched over. Legacy: ONE window `[min(L,R), max(L,R)]` (the TSD, a
/// target-site deletion's deleted bases, or a blunt junction). With `split_breakpoint_span`
/// set, a far pair (`|R - L|` above it: L1-mediated deletion / duplication) gets one 1-bp
/// window per breakpoint, in coordinate order. A one-sided locus names one breakpoint twice,
/// so its single window is that breakpoint either way.
pub fn breakpoint_windows(ins: &Insertion, cfg: &GenotypingConfig) -> Vec<(i64, i64)> {
    let lo = ins.left_pos.min(ins.right_pos);
    let hi = ins.left_pos.max(ins.right_pos);
    if cfg.split_breakpoint_span > 0 && hi - lo > cfg.split_breakpoint_span {
        vec![(lo, lo), (hi, hi)]
    } else {
        vec![(lo, hi)]
    }
}

/// One side's score triple (ref, alt, art).
pub type SideScore = (i64, i64, i64);

/// True when a read at a one-sided locus is a junction read of the MISSING end: its
/// alignment ends, on the open side, in a soft clip of >= `one_sided_open_min_clip` bases
/// within `one_sided_open_window` bp of the real breakpoint. For `L-oneside_L` (right end
/// open) that is a trailing clip (`..M nS`) whose clip point (reference_end) lies near L;
/// for `oneside_R-R` (left end open) a leading clip (`nS M..`) starting near R. Such a read's
/// aligned part crosses the real breakpoint in reference configuration, so scoring only the
/// real side would turn an insertion read into a ref vote.
pub fn is_open_side_junction_read(
    ins: &Insertion,
    cfg: &GenotypingConfig,
    cigar: &[(u8, usize)],
    span: (i64, i64),
) -> bool {
    if !cfg.one_sided_loci || cfg.one_sided_open_window <= 0 || !ins.is_one_sided() {
        return false;
    }
    let w = cfg.one_sided_open_window;
    let min_clip = cfg.one_sided_open_min_clip.max(1) as usize;
    if ins.right_open {
        // real LEFT breakpoint at left_pos; the missing RIGHT junction clips a read's 3' end
        match cigar.last() {
            Some(&(4, len)) if len >= min_clip => (span.1 - ins.left_pos).abs() <= w,
            _ => false,
        }
    } else {
        match cigar.first() {
            Some(&(4, len)) if len >= min_clip => (span.0 - ins.right_pos).abs() <= w,
            _ => false,
        }
    }
}

/// Per-read (left_side, right_side) junction scores, `Ok(None)` for a read that must not
/// vote (a one-sided locus's open-side junction read). A one-sided locus (with
/// `one_sided_loci`) never scores its open end: it carries no consensus, and its position
/// repeats the real breakpoint. `Err(())` = a needed consensus is missing (-> locus `error`).
pub fn score_read(
    ins: &Insertion,
    cfg: &GenotypingConfig,
    cigar: &[(u8, usize)],
    span: (i64, i64),
    s: &[u8],
    q: &[u8],
) -> Result<Option<(SideScore, SideScore)>, ()> {
    if is_open_side_junction_read(ins, cfg, cigar, span) {
        return Ok(None);
    }
    let (reference_start, reference_end) = span;
    let left = if cfg.one_sided_loci && ins.left_open {
        Ok((0, 0, 0))
    } else {
        qleft(cigar, reference_start, reference_end, s, q, ins.left_pos, ins.left_ref.as_deref(), ins.left_clipped.as_deref())
    };
    let right = if cfg.one_sided_loci && ins.right_open {
        Ok((0, 0, 0))
    } else {
        qright(cigar, reference_start, reference_end, s, q, ins.right_pos, ins.right_ref.as_deref(), ins.right_clipped.as_deref())
    };
    Ok(Some((left?, right?)))
}

/// Reference-tally corrections for this locus kind: halve for single-junction evidence
/// (a one-sided locus genotyped on its real side, or a far pair split into two windows),
/// discount the duplicated junction of a far duplication. Legacy loci get none.
pub fn ref_adjust(ins: &Insertion, cfg: &GenotypingConfig) -> RefAdjust {
    let one_sided = cfg.one_sided_loci && ins.is_one_sided();
    let split = breakpoint_windows(ins, cfg).len() > 1;
    RefAdjust {
        halve: cfg.halve_single_junction_ref && (one_sided || split),
        dup_discount: cfg.dup_ref_discount_min_span > 0
            && !ins.is_one_sided()
            && ins.right_pos - ins.left_pos >= cfg.dup_ref_discount_min_span,
    }
}

fn genotype_one(
    src: &mut dyn RegionSource,
    header: &Header,
    ins: &Insertion,
    cfg: &GenotypingConfig,
) -> io::Result<Row> {
    let chr = ins.chr.as_bytes().to_vec();
    let windows = breakpoint_windows(ins, cfg);

    // --- coverage: pysam bam.count(chr, max(0,lo), hi+1) with read_callback='all' ---
    // One window (legacy) = the whole locus. Split far pairs: each breakpoint is gated on
    // its own depth (high-coverage if EITHER exceeds the threshold); `coverage` reports the
    // sum of the breakpoint depths.
    let mut coverage: i64 = 0;
    let mut high = false;
    for &(lo, hi) in &windows {
        let cov_start = lo.max(0);
        let cov_end = hi + 1;
        let mut wcov: i64 = 0;
        let reg = region(&chr, cov_start, cov_end);
        let query = src.query(&reg)?;
        for rec in query {
            let rec = rec?;
            let r = decode_light(rec.as_dyn(), header)?;
            // pysam bam.count() counts EVERY overlapping alignment (verified: default
            // 'all' == 'nofilter' here, and the tally includes secondary/supplementary/
            // qcfail/duplicate records). No flag filtering — just precise overlap.
            if overlaps(r.reference_start, r.reference_end, cov_start, cov_end) {
                wcov += 1;
                // Early-exit once the high-coverage call is already decided. A pileup locus
                // (satellite / rDNA / mismapping stack with millions of reads) is flagged
                // GT_HIGH_COVERAGE no matter its exact depth, so draining the whole pileup
                // just to finish counting is pure I/O with zero effect on the output — and
                // it is what wedges genotyping on such loci. This DIVERGES from pysam's full
                // count on purpose: for high-coverage loci `coverage` is now capped at
                // reads_for_high_coverage + 1 (a ">threshold" sentinel), not the true total.
                // Loci at or below the threshold are still counted in full and unchanged.
                if wcov > cfg.reads_for_high_coverage {
                    break;
                }
            }
        }
        coverage += wcov;
        if wcov > cfg.reads_for_high_coverage {
            high = true;
            break;
        }
    }

    if high {
        return Ok(Row { genotype: GT_HIGH_COVERAGE, score_gt: coverage, score_other: 0, coverage, n_alt: 0, n_ref: 0, n_art: 0 });
    }

    // --- evidence: pysam bam.fetch(chr, max(0,lo-1), hi+1) + read gates + dedup ---
    // (per window; the qname dedup set is shared, so a fragment seen at one breakpoint is
    // not counted again at the other)
    let mut reads: Vec<(SideScore, SideScore)> = Vec::new();
    let mut seen: FxHashSet<Vec<u8>> = FxHashSet::default();
    for &(lo, hi) in &windows {
        let g_start = (lo - 1).max(0);
        let g_end = hi + 1;
        let reg = region(&chr, g_start, g_end);
        let query = src.query(&reg)?;
        for rec in query {
            let rec = rec?;
            let dyn_rec = rec.as_dyn();
            let r = decode_light(dyn_rec, header)?;
            if !overlaps(r.reference_start, r.reference_end, g_start, g_end) {
                continue; // drop index over-fetch (fetch yields only overlapping reads)
            }
            // gates, in the Python order (genotyping_insertion.py::genotype)
            if r.mapq < cfg.min_mapq {
                continue;
            }
            if r.is_secondary || r.is_supplementary || r.is_unmapped || r.is_qcfail || r.is_duplicate {
                continue;
            }
            match r.reference_sequence_id.and_then(|id| reference_name(header, id)) {
                Some(n) if n == chr => {}
                _ => continue,
            }
            let s = seq(dyn_rec);
            let q = qual(dyn_rec);
            if s.is_empty() || q.is_empty() {
                continue;
            }
            // dedup by fragment: first qname wins (overlapping mates / PCR dups count once)
            let qn = name(dyn_rec);
            if !seen.insert(qn) {
                continue;
            }
            match score_read(ins, cfg, &r.cigar, (r.reference_start, r.reference_end), &s, &q) {
                Ok(Some(lr)) => reads.push(lr),
                Ok(None) => {} // open-side junction read of a one-sided locus: no vote
                // a spanning read reached scoring but a consensus was missing: the Python
                // code raises here and the locus is reported as `error` (coverage known).
                Err(()) => {
                    return Ok(Row { genotype: GT_ERROR, score_gt: 0, score_other: 0, coverage, n_alt: 0, n_ref: 0, n_art: 0 });
                }
            }
        }
    }

    let call = summarise_evidence_opts(&reads, cfg, ref_adjust(ins, cfg));
    Ok(Row {
        genotype: call.genotype,
        score_gt: call.score_gt,
        score_other: call.score_other,
        coverage,
        n_alt: call.n_alt,
        n_ref: call.n_ref,
        n_art: call.n_art,
    })
}

/// Genotype `insertions` (already parsed) against one alignment file, writing the
/// 8-column gzip output to `writer`. `reference` is required for CRAM input.
///
/// Single-threaded (the cluster default): rows are STREAMED to `writer` and flushed on the
/// heartbeat cadence, so the output file grows during the run — a hung or killed task
/// leaves a partial, decodable `.txt.gz` whose last row is exactly the last locus done.
/// Multi-threaded: loci are partitioned across `threads` workers (each with its own indexed
/// reader) and re-assembled in contract order; that path can only write once all workers
/// finish, but each still heartbeats its progress to stderr.
pub fn run<W: Write>(
    insertions: &[Insertion],
    input: &str,
    reference: Option<&str>,
    cfg: &GenotypingConfig,
    threads: usize,
    writer: &mut W,
) -> io::Result<()> {
    let n = insertions.len();
    let threads = threads.max(1).min(n.max(1));

    writer.write_all(OUTPUT_HEADER.as_bytes())?;

    if threads == 1 {
        let mut done = 0usize;
        process_chunk_with(input, reference, insertions, cfg, "", |row| {
            writer.write_all(row.as_bytes())?;
            done += 1;
            // A sync-flush of the gzip stream every HEARTBEAT_EVERY rows makes the partial
            // output on disk decodable up to that point without ending the stream.
            if done % HEARTBEAT_EVERY == 0 {
                writer.flush()?;
            }
            Ok(())
        })?;
        writer.flush()?;
        return Ok(());
    }

    // contiguous chunks preserve contract order after concatenation.
    let base = n / threads;
    let rem = n % threads;
    let mut bounds = Vec::with_capacity(threads + 1);
    bounds.push(0usize);
    let mut acc = 0usize;
    for t in 0..threads {
        acc += base + if t < rem { 1 } else { 0 };
        bounds.push(acc);
    }

    let mut chunk_rows: Vec<io::Result<Vec<String>>> = Vec::with_capacity(threads);
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(threads);
        for t in 0..threads {
            let (lo, hi) = (bounds[t], bounds[t + 1]);
            let slice = &insertions[lo..hi];
            let tag = format!("[t{t}] ");
            handles.push(scope.spawn(move || process_chunk(input, reference, slice, cfg, &tag)));
        }
        for h in handles {
            chunk_rows.push(h.join().expect("genotyping worker thread panicked"));
        }
    });

    for chunk in chunk_rows {
        for row in chunk? {
            writer.write_all(row.as_bytes())?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Parse a gzipped genotyping contract into its loci (shared by single and batch runs).
pub fn load_contract(insertion_file: &str) -> io::Result<Vec<Insertion>> {
    Insertion::import_file(insertion_file)
}

#[cfg(test)]
mod tests {
    //! Per-read scoring of the TPRT locus kinds on synthetic reads. Filler bases are `N`
    //! (never scored), so only the 12 bases next to a breakpoint decide a side.
    use super::*;

    const L_ALT: &[u8] = b"ACGTACGTACGT"; // left clip consensus (next to the junction)
    const L_REF: &[u8] = b"TTGGCCAATTGG"; // genome[L-12, L)
    const R_REF: &[u8] = b"GATCGATCGATC"; // genome[R, R+12)
    const R_ALT: &[u8] = b"TGCATGCATGCA"; // right clip consensus

    fn locus(name: &str, left: bool, right: bool) -> Insertion {
        let mut i = Insertion::new(name).unwrap();
        if left {
            i.left_clipped = Some(L_ALT.to_vec());
            i.left_ref = Some(L_REF.to_vec());
        }
        if right {
            i.right_clipped = Some(R_ALT.to_vec());
            i.right_ref = Some(R_REF.to_vec());
        }
        i
    }

    /// A read: `n` bases of N with `motif` placed at query offset `at`.
    fn read(n: usize, motif: &[u8], at: usize) -> (Vec<u8>, Vec<u8>) {
        let mut s = vec![b'N'; n];
        s[at..at + motif.len()].copy_from_slice(motif);
        (s, vec![30u8; n])
    }

    fn tprt_cfg() -> GenotypingConfig {
        GenotypingConfig {
            min_mapq: 60,
            one_sided_loci: true,
            split_breakpoint_span: 40,
            one_sided_open_window: 50,
            dup_ref_discount_min_span: 150,
            halve_single_junction_ref: true,
            ..GenotypingConfig::default()
        }
    }

    fn kind(side: SideScore) -> &'static str {
        match side {
            (0, 0, 0) => "-",
            (r, a, _) if a > r => "alt",
            _ => "ref",
        }
    }

    fn call(
        ins: &Insertion,
        cfg: &GenotypingConfig,
        cigar: &[(u8, usize)],
        start: i64,
        sq: &(Vec<u8>, Vec<u8>),
    ) -> Option<(&'static str, &'static str)> {
        let rlen: usize = cigar.iter().filter(|(op, _)| matches!(op, 0 | 2 | 3 | 7 | 8)).map(|&(_, l)| l).sum();
        score_read(ins, cfg, cigar, (start, start + rlen as i64), &sq.0, &sq.1)
            .unwrap()
            .map(|(l, r)| (kind(l), kind(r)))
    }

    #[test]
    fn windows_legacy_and_split() {
        let legacy = GenotypingConfig::default();
        let cfg = tprt_cfg();
        let tsd = locus("chr1:1000-1015", true, true);
        let tsd_del = locus("chr1:1000-980", true, true);
        let l1del = locus("chr1:5000-1000", true, true);
        let l1dup = locus("chr1:1000-1300", true, true);
        let one = locus("chr1:1000-oneside_1000", true, false);
        assert_eq!(breakpoint_windows(&l1del, &legacy), vec![(1000, 5000)]);
        assert_eq!(breakpoint_windows(&tsd, &cfg), vec![(1000, 1015)]);
        assert_eq!(breakpoint_windows(&tsd_del, &cfg), vec![(980, 1000)]);
        assert_eq!(breakpoint_windows(&l1del, &cfg), vec![(1000, 1000), (5000, 5000)]);
        assert_eq!(breakpoint_windows(&l1dup, &cfg), vec![(1000, 1000), (1300, 1300)]);
        assert_eq!(breakpoint_windows(&one, &cfg), vec![(1000, 1000)]);
        // reference corrections: legacy none; split far pairs halve; far dup also discounts
        assert_eq!(ref_adjust(&tsd, &cfg), RefAdjust::default());
        assert_eq!(ref_adjust(&tsd_del, &cfg), RefAdjust::default());
        assert_eq!(ref_adjust(&l1del, &cfg), RefAdjust { halve: true, dup_discount: false });
        assert_eq!(ref_adjust(&l1dup, &cfg), RefAdjust { halve: true, dup_discount: true });
        assert_eq!(ref_adjust(&one, &cfg), RefAdjust { halve: true, dup_discount: false });
        assert_eq!(ref_adjust(&l1dup, &legacy), RefAdjust::default());
    }

    #[test]
    fn l1_deletion_reads_score_their_own_breakpoint() {
        // L1-mediated deletion: genome[R, L) deleted, L=5000 > R=1000.
        let ins = locus("chr1:5000-1000", true, true);
        let cfg = tprt_cfg();
        // alt at L: [20 clip][80 aligned from L]
        assert_eq!(call(&ins, &cfg, &[(4, 20), (0, 80)], 5000, &read(100, L_ALT, 8)), Some(("alt", "-")));
        // ref at L: 100M from 4950, the 12 bases before L carry genome[L-12, L)
        assert_eq!(call(&ins, &cfg, &[(0, 100)], 4950, &read(100, L_REF, 38)), Some(("ref", "-")));
        // alt at R: [80 aligned up to R-1][20 clip]
        assert_eq!(call(&ins, &cfg, &[(0, 80), (4, 20)], 920, &read(100, R_ALT, 80)), Some(("-", "alt")));
        // ref at R: 100M from 950, genome[R, R+12) at offset 50
        assert_eq!(call(&ins, &cfg, &[(0, 100)], 950, &read(100, R_REF, 50)), Some(("-", "ref")));
        // per-colony het: 1 alt + 1 ref per breakpoint, x3 -> halved ref -> heterozygous
        let mut reads = Vec::new();
        for _ in 0..3 {
            reads.push(((-1, 360, 0), (0, 0, 0)));
            reads.push(((360, -1, 0), (0, 0, 0)));
            reads.push(((0, 0, 0), (-1, 360, 0)));
            reads.push(((0, 0, 0), (360, -1, 0)));
        }
        let c = summarise_evidence_opts(&reads, &cfg, ref_adjust(&ins, &cfg));
        assert_eq!((c.genotype, c.n_alt, c.n_ref), (GT_HETEROZYGOUS, 6, 3));
    }

    #[test]
    fn target_site_deletion_alt_reads_never_double_alt() {
        // L=1000 > R=980: genome[980, 1000) deleted. A clipped alt read covers only its own
        // junction (the other breakpoint lies outside its aligned span) -> never `artefact`.
        let ins = locus("chr1:1000-980", true, true);
        let cfg = GenotypingConfig::default();
        assert_eq!(call(&ins, &cfg, &[(4, 20), (0, 80)], 1000, &read(100, L_ALT, 8)), Some(("alt", "-")));
        assert_eq!(call(&ins, &cfg, &[(0, 80), (4, 20)], 900, &read(100, R_ALT, 80)), Some(("-", "alt")));
        // a reference read spans both breakpoints: genome[980, 992) at offset 50, genome[988,
        // 1000) at offset 58 (consistent: both are the reference) -> ref on both sides
        let mut s = vec![b'N'; 150];
        s[50..62].copy_from_slice(R_REF);
        s[58..70].copy_from_slice(L_REF);
        let sq = (s, vec![30u8; 150]);
        assert_eq!(call(&ins, &cfg, &[(0, 150)], 930, &sq).unwrap().0, "ref");
    }

    #[test]
    fn blunt_junction() {
        let ins = locus("chr1:1000-1000", true, true);
        let cfg = GenotypingConfig::default();
        assert_eq!(call(&ins, &cfg, &[(4, 20), (0, 80)], 1000, &read(100, L_ALT, 8)), Some(("alt", "-")));
        assert_eq!(call(&ins, &cfg, &[(0, 80), (4, 20)], 920, &read(100, R_ALT, 80)), Some(("-", "alt")));
        let mut s = vec![b'N'; 100];
        s[38..50].copy_from_slice(L_REF);
        s[50..62].copy_from_slice(R_REF);
        assert_eq!(call(&ins, &cfg, &[(0, 100)], 950, &(s, vec![30u8; 100])), Some(("ref", "ref")));
    }

    #[test]
    fn one_sided_scores_real_side_only() {
        let left_real = locus("chr1:1000-oneside_1000", true, false);
        let right_real = locus("chr1:oneside_2000-2000", false, true);
        let cfg = tprt_cfg();
        // real-side alt / ref reads
        assert_eq!(call(&left_real, &cfg, &[(4, 20), (0, 80)], 1000, &read(100, L_ALT, 8)), Some(("alt", "-")));
        assert_eq!(call(&left_real, &cfg, &[(0, 100)], 950, &read(100, L_REF, 38)), Some(("ref", "-")));
        assert_eq!(call(&right_real, &cfg, &[(0, 80), (4, 20)], 1920, &read(100, R_ALT, 80)), Some(("-", "alt")));
        // the missing end's junction read (3' clip 15 bp after L): its aligned part crosses L
        // in reference configuration -> no vote instead of a false ref vote
        assert_eq!(call(&left_real, &cfg, &[(0, 65), (4, 35)], 950, &read(100, L_REF, 38)), None);
        // ... and the mirror for a right-real locus (5' clip 15 bp before R)
        assert_eq!(call(&right_real, &cfg, &[(4, 30), (0, 70)], 1985, &read(100, R_REF, 45)), None);
        // a clip far from the breakpoint, or too short, is an ordinary read
        assert_eq!(call(&left_real, &cfg, &[(0, 140), (4, 10)], 950, &read(150, L_REF, 38)), Some(("ref", "-")));
        assert_eq!(call(&left_real, &cfg, &[(0, 96), (4, 4)], 950, &read(100, L_REF, 38)), Some(("ref", "-")));
        // without one_sided_open_window the junction read votes ref (the bias it removes)
        let no_win = GenotypingConfig { one_sided_open_window: 0, ..tprt_cfg() };
        assert_eq!(call(&left_real, &no_win, &[(0, 65), (4, 35)], 950, &read(100, L_REF, 38)), Some(("ref", "-")));
        // with one_sided_loci off the open side has no consensus -> locus error (legacy)
        let legacy = GenotypingConfig::default();
        let sq = read(100, L_REF, 38);
        assert!(score_read(&left_real, &legacy, &[(0, 100)], (950, 1050), &sq.0, &sq.1).is_err());
    }
}
