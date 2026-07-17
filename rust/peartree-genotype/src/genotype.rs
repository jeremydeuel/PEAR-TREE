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
use crate::insertion::{summarise_evidence, Insertion};
use crate::read::{decode_light, name, qual, reference_name, seq};
use crate::source::{open_source, RegionSource};

use noodles_core::{Position, Region};
use noodles_sam::Header;
use rustc_hash::FxHashSet;
use std::io::{self, Write};

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

/// Genotype every insertion in `insertions` against `path`, returning the
/// formatted output rows in order. Opens its own indexed reader (thread-local).
/// `reference` is required for CRAM.
fn process_chunk(
    path: &str,
    reference: Option<&str>,
    insertions: &[Insertion],
    cfg: &GenotypingConfig,
) -> io::Result<Vec<String>> {
    let mut src = open_source(path, reference)?;
    // Clone the header once so decode can borrow it while `src` is borrowed mutably
    // by the region queries (the two would otherwise conflict).
    let header = src.header().clone();

    let mut rows = Vec::with_capacity(insertions.len());
    for ins in insertions {
        // Blanket per-locus error handling mirrors the Python worker's try/except:
        // one bad locus becomes an `error` row, never a crash.
        let row = match genotype_one(src.as_mut(), &header, ins, cfg) {
            Ok(row) => row,
            Err(e) => {
                eprintln!("genotyping failed for {}: {e}", ins.name);
                Row { genotype: GT_ERROR, score_gt: 0, score_other: 0, coverage: 0, n_alt: 0, n_ref: 0, n_art: 0 }
            }
        };
        rows.push(format_row(&ins.name, &row));
    }
    Ok(rows)
}

fn genotype_one(
    src: &mut dyn RegionSource,
    header: &Header,
    ins: &Insertion,
    cfg: &GenotypingConfig,
) -> io::Result<Row> {
    let chr = ins.chr.as_bytes().to_vec();
    let lo = ins.left_pos.min(ins.right_pos);
    let hi = ins.left_pos.max(ins.right_pos);

    // --- coverage: pysam bam.count(chr, max(0,lo), hi+1) with read_callback='all' ---
    let cov_start = lo.max(0);
    let cov_end = hi + 1;
    let mut coverage: i64 = 0;
    {
        let reg = region(&chr, cov_start, cov_end);
        let query = src.query(&reg)?;
        for rec in query {
            let rec = rec?;
            let r = decode_light(rec.as_dyn(), header)?;
            // pysam bam.count() counts EVERY overlapping alignment (verified: default
            // 'all' == 'nofilter' here, and the tally includes secondary/supplementary/
            // qcfail/duplicate records). No flag filtering — just precise overlap.
            if overlaps(r.reference_start, r.reference_end, cov_start, cov_end) {
                coverage += 1;
                // Early-exit once the high-coverage call is already decided. A pileup locus
                // (satellite / rDNA / mismapping stack with millions of reads) is flagged
                // GT_HIGH_COVERAGE no matter its exact depth, so draining the whole pileup
                // just to finish counting is pure I/O with zero effect on the output — and
                // it is what wedges genotyping on such loci. This DIVERGES from pysam's full
                // count on purpose: for high-coverage loci `coverage` is now capped at
                // reads_for_high_coverage + 1 (a ">threshold" sentinel), not the true total.
                // Loci at or below the threshold are still counted in full and unchanged.
                if coverage > cfg.reads_for_high_coverage {
                    break;
                }
            }
        }
    }

    if coverage > cfg.reads_for_high_coverage {
        return Ok(Row { genotype: GT_HIGH_COVERAGE, score_gt: coverage, score_other: 0, coverage, n_alt: 0, n_ref: 0, n_art: 0 });
    }

    // --- evidence: pysam bam.fetch(chr, max(0,lo-1), hi+1) + read gates + dedup ---
    let g_start = (lo - 1).max(0);
    let g_end = hi + 1;
    let mut reads: Vec<((i64, i64, i64), (i64, i64, i64))> = Vec::new();
    let mut seen: FxHashSet<Vec<u8>> = FxHashSet::default();
    {
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
            let left = qleft(
                &r.cigar, r.reference_start, r.reference_end, &s, &q,
                ins.left_pos, ins.left_ref.as_deref(), ins.left_clipped.as_deref(),
            );
            let right = qright(
                &r.cigar, r.reference_start, r.reference_end, &s, &q,
                ins.right_pos, ins.right_ref.as_deref(), ins.right_clipped.as_deref(),
            );
            match (left, right) {
                (Ok(l), Ok(rr)) => reads.push((l, rr)),
                // a spanning read reached scoring but a consensus was missing: the Python
                // code raises here and the locus is reported as `error` (coverage known).
                _ => {
                    return Ok(Row { genotype: GT_ERROR, score_gt: 0, score_other: 0, coverage, n_alt: 0, n_ref: 0, n_art: 0 });
                }
            }
        }
    }

    let call = summarise_evidence(&reads, cfg);
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
/// 8-column gzip output to `writer`. Loci are partitioned across `threads` workers
/// (each with its own indexed reader) and re-assembled in contract order.
/// `reference` is required for CRAM input.
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
            handles.push(scope.spawn(move || process_chunk(input, reference, slice, cfg)));
        }
        for h in handles {
            chunk_rows.push(h.join().expect("genotyping worker thread panicked"));
        }
    });

    writer.write_all(OUTPUT_HEADER.as_bytes())?;
    for chunk in chunk_rows {
        for row in chunk? {
            writer.write_all(row.as_bytes())?;
        }
    }
    Ok(())
}

/// Parse a gzipped genotyping contract into its loci (shared by single and batch runs).
pub fn load_contract(insertion_file: &str) -> io::Result<Vec<Insertion>> {
    Insertion::import_file(insertion_file)
}
