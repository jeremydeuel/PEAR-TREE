//! Per-colony genotyping driver (SPEC "Driver / I/O"). Owner: D. Owns source.rs / read.rs too.
//!
//! Pipeline per alignment file:
//!   1. every `LocusModel` is built up front from the reference (`haplotype::build_model`);
//!   2. loci are sorted by (contig order in the BAM header, min breakpoint) and processed in
//!      that order — consecutive index queries then hit nearby, usually already-buffered data;
//!   3. `--threads N` splits the sorted list into N contiguous chunks, one reader each;
//!   4. per locus, ONE indexed query per fetch window feeds both the depth gate (early exit
//!      above `reads_for_high_coverage`) and the evidence (gated, deduped reads -> `score_read`)
//!      and the discordant-anchor count;
//!   5. `model::call_locus` -> `output::format_row`; rows are streamed in processing order
//!      with a sync flush + stderr heartbeat every `heartbeat_every` loci.
//!
//! The record-level logic (`gate`, `passes_mapq`, `is_discordant_anchor`, `Dedup`,
//! `collect_locus`) works on decoded [`LightRec`]s behind the small [`LocusReads`] /
//! [`RecordView`] traits, so it is unit-tested without a BAM and without the scorer.

use crate::config::Config;
use crate::contract::ContractSides;
use crate::read::{decode_light_into, name, qual_into, seq_into, LightRec};
use crate::source::{io_counters, is_cram, open_source_buffered, AnyRecord, RegionSource};
use crate::types::{output_header, AltSide, Locus, LocusModel, ReadClass, ReadInput, ReadObs, Status};
use crate::{extra, haplotype, model, output, readlik, refseq};

use noodles_core::{Position, Region};
use noodles_sam::Header;
use rustc_hash::FxHashSet;
use std::collections::HashMap;
use std::io::{self, Write};
use std::ops::Range;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;
use std::time::Instant;

/// A trailing (leading) soft clip "faces" the R (L) junction when the clip point lies within
/// this many bp of it (aligners extend a few bases into a poly-A / microhomology, or clip a
/// mismatching base or two next to the junction).
const CLIP_FACING_TOL: i64 = 10;

/// The widened query (breakpoint window ± `disc_span`, needed for the discordant anchors) is
/// abandoned after this many records times (`reads_for_high_coverage` + 1) and the narrow
/// window `[lo-1, hi+1)` is re-queried instead, so a pileup NEXT TO a breakpoint (satellite,
/// rDNA) can never wedge a locus — the cost of the extension stays bounded.
const WIDE_QUERY_CAP_FACTOR: usize = 10;
const WIDE_QUERY_CAP_MIN: usize = 2000;

// ------------------------------------------------------------------------------------------
// public entry points
// ------------------------------------------------------------------------------------------

/// Genotype `loci` (contract order; the driver sorts internally) against one alignment file,
/// writing rows to `writer` (header included). `reference_path` is the genome (FASTA/2bit) used
/// for haplotypes AND for CRAM decoding. `combined` = the per-locus consensus from `--combined`
/// (None -> contract fallback with one warning). `extra` (only with `gt_extra_reads`): this
/// colony's discovery memberships + the sidecar writer of the extra pass (extra.rs).
#[allow(clippy::too_many_arguments)]
pub fn run<W: Write>(
    loci: &[(Locus, ContractSides)],
    combined: Option<&std::collections::HashMap<String, ContractSides>>,
    input: &str,
    reference_path: &str,
    cfg: &Config,
    threads: usize,
    writer: &mut W,
    extra: Option<ExtraArgs<'_>>,
) -> io::Result<()> {
    // open the alignment first: a bad input (e.g. CRAM with a .2bit reference) fails before the
    // model build
    let mut src0 = open_source_buffered(input, Some(reference_path), cfg.io_buffer_bytes, cfg.io_fill_bytes)?;
    let header = src0.header().clone();
    let t0 = Instant::now();
    let models = build_models(loci, combined, reference_path, cfg)?;
    let t_models = t0.elapsed();
    let order = sort_order(&models, |chr| contig_index(&header, chr));
    let chunks = partition(order.len(), threads);
    eprintln!(
        "models: {} loci built in {:.1}s; processing in header/coordinate order, {} chunk(s), io buffer {} KiB (first fill {} KiB)",
        models.len(),
        t_models.as_secs_f64(),
        chunks.len(),
        cfg.io_buffer_bytes >> 10,
        cfg.io_fill_bytes >> 10
    );

    if cfg.ref_bias != 1.0 || !cfg.ref_bias_kind.is_empty() || !cfg.ref_bias_grid.is_empty() {
        let kinds: Vec<String> = cfg.ref_bias_kind.iter().map(|(k, b)| format!("{k}:{b}")).collect();
        eprintln!(
            "reference bias: het alt fraction corrected with b = {} (per kind: {}); pl_het_b profile at b in {:?}",
            cfg.ref_bias,
            if kinds.is_empty() { "none".to_string() } else { kinds.join(",") },
            cfg.ref_bias_grid
        );
    }
    writer.write_all(output_header(&cfg.noise_frac_grid, &cfg.ref_bias_grid).as_bytes())?;
    let t_geno = Instant::now();
    let io0 = io_counters();
    // the extra pass: per chunk, the members set + an even share of the per-colony mate budget
    let n_chunks = chunks.len().max(1);
    let meta = extra.as_ref().map(|e| (e.members, e.sample));
    let mut xout: Option<&mut dyn Write> = extra.map(|e| e.out);
    let new_state = || {
        meta.map(|(members, sample)| ExtraState {
            members,
            sample,
            mate_budget: cfg.gt_extra_max_mate_fetches / n_chunks,
            stats: extra::ExtraStats::default(),
            buf: String::new(),
        })
    };
    let mut extra_stats = extra::ExtraStats::default();

    if chunks.len() <= 1 {
        let mut xs = new_state();
        process_chunk(
            src0.as_mut(),
            &header,
            &models,
            &order,
            cfg,
            "",
            xs.as_mut(),
            &mut |text| match xout.as_mut() {
                Some(w) => w.write_all(text.as_bytes()),
                None => Ok(()),
            },
            &mut |row, heartbeat| {
                writer.write_all(row.as_bytes())?;
                if heartbeat {
                    // sync-flush: the partial .txt.gz on disk decodes up to the last heartbeat
                    writer.flush()?;
                }
                Ok(())
            },
        )?;
        if let Some(x) = xs {
            extra_stats.add(&x.stats);
        }
    } else {
        // Chunk 0 runs on this thread and STREAMS to `writer` (so a killed job still leaves
        // the first chunk's rows); chunks 1.. run on worker threads, each with its own reader,
        // and are appended in chunk order after the join (their extra reads likewise).
        type ChunkOut = (String, String, Option<extra::ExtraStats>);
        let mut collected: Vec<io::Result<ChunkOut>> = Vec::new();
        let mut xs0 = new_state();
        std::thread::scope(|scope| -> io::Result<()> {
            let mut handles = Vec::with_capacity(chunks.len() - 1);
            for (t, range) in chunks.iter().enumerate().skip(1) {
                let idx = &order[range.clone()];
                let models = &models;
                let tag = format!("[t{t}] ");
                let mut xs = new_state();
                handles.push(scope.spawn(move || -> io::Result<ChunkOut> {
                    let mut src = open_source_buffered(input, Some(reference_path), cfg.io_buffer_bytes, cfg.io_fill_bytes)?;
                    let header = src.header().clone();
                    let (mut out, mut xtext) = (String::new(), String::new());
                    process_chunk(
                        src.as_mut(),
                        &header,
                        models,
                        idx,
                        cfg,
                        &tag,
                        xs.as_mut(),
                        &mut |text| {
                            xtext.push_str(text);
                            Ok(())
                        },
                        &mut |row, _| {
                            out.push_str(row);
                            Ok(())
                        },
                    )?;
                    Ok((out, xtext, xs.map(|x| x.stats)))
                }));
            }
            let first = process_chunk(
                src0.as_mut(),
                &header,
                &models,
                &order[chunks[0].clone()],
                cfg,
                "[t0] ",
                xs0.as_mut(),
                &mut |text| match xout.as_mut() {
                    Some(w) => w.write_all(text.as_bytes()),
                    None => Ok(()),
                },
                &mut |row, hb| {
                    writer.write_all(row.as_bytes())?;
                    if hb {
                        writer.flush()?;
                    }
                    Ok(())
                },
            );
            for h in handles {
                collected.push(h.join().unwrap_or_else(|_| Err(io::Error::other("genotyping worker thread panicked"))));
            }
            first
        })?;
        if let Some(x) = xs0 {
            extra_stats.add(&x.stats);
        }
        for chunk in collected {
            let (rows, xtext, xstats) = chunk?;
            writer.write_all(rows.as_bytes())?;
            if let Some(w) = xout.as_mut() {
                w.write_all(xtext.as_bytes())?;
            }
            if let Some(s) = xstats {
                extra_stats.add(&s);
            }
        }
    }
    writer.flush()?;
    if let (Some(w), Some((_, sample))) = (xout.as_mut(), meta) {
        w.flush()?;
        let s = &extra_stats;
        eprintln!(
            "extra reads ({sample}): {} undiscovered loci called, {} passed the ALT gate (>= {} Alt read(s)); wrote {} GT_CLIP, \
             {} GT_POLYA, {} GT_DISC, {} GT_MATE ({} mates not found, {} over the fetch cap); extra pass {:.1}s",
            s.candidates,
            s.gated,
            cfg.gt_extra_min_alt,
            s.clip,
            s.polya,
            s.disc,
            s.mates,
            s.mates_missing,
            s.mates_capped,
            s.micros as f64 / 1e6
        );
    }
    let secs = t_geno.elapsed().as_secs_f64();
    let io1 = io_counters();
    eprintln!(
        "genotyped {} loci in {:.1}s ({:.2} ms/locus, {} thread(s)); file I/O: {} reads, {:.1} MiB, {} buffer-miss seeks",
        models.len(),
        secs,
        if models.is_empty() { 0.0 } else { 1000.0 * secs / models.len() as f64 },
        chunks.len().max(1),
        io1.0 - io0.0,
        (io1.1 - io0.1) as f64 / (1u64 << 20) as f64,
        io1.2 - io0.2
    );
    Ok(())
}

/// Validate an input alignment file up front (exists, indexed).
pub fn check_input(path: &str) -> Result<(), String> {
    if !Path::new(path).is_file() {
        return Err(format!("input file {path} does not exist"));
    }
    let candidates: Vec<String> =
        if is_cram(path) { vec![format!("{path}.crai")] } else { vec![format!("{path}.bai"), format!("{path}.csi")] };
    if candidates.iter().any(|c| Path::new(c).is_file()) {
        Ok(())
    } else {
        Err(format!("input {path} is not indexed: none of {} found (samtools index {path})", candidates.join(", ")))
    }
}

// ------------------------------------------------------------------------------------------
// models, ordering, chunking
// ------------------------------------------------------------------------------------------

/// The consensus to model a locus with: the combined file's sides, each missing side filled
/// from the contract's 12-bp consensus. Returns (sides, whether a NEEDED side fell back).
fn merged_sides(locus: &Locus, contract: &ContractSides, combined: Option<&HashMap<String, ContractSides>>) -> (ContractSides, bool) {
    let Some(map) = combined else {
        return (contract.clone(), false); // main.rs already warned about the missing --combined
    };
    let c = map.get(&locus.name);
    let left = c.and_then(|c| c.left.clone());
    let right = c.and_then(|c| c.right.clone());
    let fell_back = (left.is_none() && !locus.left_open) || (right.is_none() && !locus.right_open);
    let sides = ContractSides {
        left: left.or_else(|| contract.left.clone()),
        right: right.or_else(|| contract.right.clone()),
    };
    (sides, fell_back)
}

fn build_models(
    loci: &[(Locus, ContractSides)],
    combined: Option<&HashMap<String, ContractSides>>,
    reference_path: &str,
    cfg: &Config,
) -> io::Result<Vec<LocusModel>> {
    let mut reference = refseq::open_reference(reference_path)
        .map_err(|e| io::Error::new(e.kind(), format!("cannot open reference {reference_path}: {e}")))?;
    let mut models = Vec::with_capacity(loci.len());
    let mut n_fallback = 0usize;
    let mut first_fallback: Option<&str> = None;
    let mut n_error = 0usize;
    for (locus, contract) in loci {
        let (sides, fell_back) = merged_sides(locus, contract, combined);
        if fell_back {
            n_fallback += 1;
            first_fallback.get_or_insert(&locus.name);
        }
        let built = catch_unwind(AssertUnwindSafe(|| haplotype::build_model(locus, &sides, reference.as_mut(), cfg)));
        let m = match built {
            Ok(Ok(m)) => m,
            Ok(Err(e)) => LocusModel::error(locus.clone(), format!("reference: {e}")),
            Err(_) => LocusModel::error(locus.clone(), "haplotype model construction panicked"),
        };
        if let Some(reason) = &m.error {
            n_error += 1;
            if n_error <= 10 {
                eprintln!("locus {}: cannot be modelled ({reason}) -> error row", locus.name);
            }
        }
        models.push(m);
    }
    if n_error > 10 {
        eprintln!("... {n_error} loci in total cannot be modelled (error rows)");
    }
    if n_fallback > 0 {
        eprintln!(
            "WARNING: {n_fallback}/{} loci lack a junction consensus in --combined; their missing side(s) use the \
             contract's 12-bp consensus (first: {})",
            loci.len(),
            first_fallback.unwrap_or("?")
        );
    }
    Ok(models)
}

fn contig_index(header: &Header, chr: &str) -> Option<usize> {
    header.reference_sequences().get_index_of(chr.as_bytes())
}

/// Processing order: (contig index in the alignment header, min(L, R), contract index).
/// Contigs missing from the header go last (they become `error` rows).
pub(crate) fn sort_order(models: &[LocusModel], contig_index: impl Fn(&str) -> Option<usize>) -> Vec<usize> {
    let mut keyed: Vec<(usize, i64, usize)> = models
        .iter()
        .enumerate()
        .map(|(i, m)| {
            let l = &m.locus;
            (contig_index(&l.chr).unwrap_or(usize::MAX), l.left_pos.min(l.right_pos), i)
        })
        .collect();
    keyed.sort_unstable();
    keyed.into_iter().map(|(_, _, i)| i).collect()
}

/// `threads` contiguous, near-equal ranges over `n` items (never more ranges than items; at
/// least one range, possibly empty, so a run with no loci still writes its header).
pub(crate) fn partition(n: usize, threads: usize) -> Vec<Range<usize>> {
    let k = threads.max(1).min(n.max(1));
    let (base, rem) = (n / k, n % k);
    let mut out = Vec::with_capacity(k);
    let mut lo = 0;
    for t in 0..k {
        let hi = lo + base + usize::from(t < rem);
        out.push(lo..hi);
        lo = hi;
    }
    out
}

/// Heartbeat (stderr line + sync flush) after locus `done` (1-based) of `total`.
pub(crate) fn is_heartbeat(done: usize, total: usize, every: usize) -> bool {
    done == total || (every > 0 && done % every == 0)
}

// ------------------------------------------------------------------------------------------
// chunk / locus processing
// ------------------------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn process_chunk(
    src: &mut dyn RegionSource,
    header: &Header,
    models: &[LocusModel],
    idx: &[usize],
    cfg: &Config,
    tag: &str,
    mut extra: Option<&mut ExtraState<'_>>,
    emit_extra: &mut dyn FnMut(&str) -> io::Result<()>,
    emit: &mut dyn FnMut(&str, bool) -> io::Result<()>,
) -> io::Result<()> {
    let total = idx.len();
    let t0 = Instant::now();
    let mut sum = LocusStats::default();
    for (k, &i) in idx.iter().enumerate() {
        let m = &models[i];
        let (row, stats) = genotype_locus(src, header, m, cfg, extra.as_deref_mut());
        sum.fetched += stats.fetched;
        sum.gated += stats.gated;
        sum.requeried += stats.requeried;
        if let Some(x) = extra.as_deref_mut() {
            if !x.buf.is_empty() {
                emit_extra(&x.buf)?;
                x.buf.clear();
            }
        }
        let hb = is_heartbeat(k + 1, total, cfg.heartbeat_every);
        emit(&row, hb)?;
        if hb {
            eprintln!("  {tag}{}/{} loci, {}s (last {})", k + 1, total, t0.elapsed().as_secs(), m.locus.name);
        }
    }
    if total > 0 {
        eprintln!(
            "  {tag}records streamed {} ({:.0}/locus), evidence reads gated {}, pileup re-queries {}",
            sum.fetched,
            sum.fetched as f64 / total as f64,
            sum.gated,
            sum.requeried
        );
    }
    Ok(())
}

/// One locus -> one output row (+ stream diagnostics). Never fails: I/O errors and panics in
/// the scorer / model become an `error` row (with the coverage when it is already known).
///
/// With `extra` (gt_extra_reads), a locus this colony did not discover keeps its scored reads
/// and, once the row is made, gets the extra pass if it passes the ALT gate. The row never
/// depends on it; a failing extra pass only warns and writes nothing for the locus.
fn genotype_locus(
    src: &mut dyn RegionSource,
    header: &Header,
    m: &LocusModel,
    cfg: &Config,
    extra: Option<&mut ExtraState<'_>>,
) -> (String, LocusStats) {
    let name = m.locus.name.as_str();
    let n_prof = cfg.noise_frac_grid.len() + cfg.ref_bias_grid.len();
    if m.error.is_some() {
        return (output::format_simple_row(name, m.kind, Status::Error, 0, n_prof), LocusStats::default());
    }
    let extra = extra.filter(|x| is_extra_candidate(x.members, name));
    let mut kept: Vec<KeptRead> = Vec::new();
    let collected = catch_unwind(AssertUnwindSafe(|| -> io::Result<Outcome> {
        let contig = contig_index(header, &m.locus.chr).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, format!("contig {} not in the alignment header", m.locus.chr))
        })?;
        let mut reads = SourceReads::new(src, header, &m.locus.chr, contig);
        let keep = if extra.is_some() { Some(&mut kept) } else { None };
        collect_locus_keep(m, contig, cfg, &mut reads, &mut |r| readlik::score_read(m, r, cfg), keep)
    }));
    let (coverage, obs, n_disc, stats) = match collected {
        Ok(Ok(Outcome::HighCoverage)) => {
            let row = output::format_simple_row(name, m.kind, Status::HighCoverage, cfg.reads_for_high_coverage + 1, n_prof);
            return (row, LocusStats::default());
        }
        Ok(Ok(Outcome::Collected { coverage, obs, n_disc, stats })) => (coverage, obs, n_disc, stats),
        Ok(Err(e)) => {
            eprintln!("genotyping failed for {name}: {e}");
            return (output::format_simple_row(name, m.kind, Status::Error, 0, n_prof), LocusStats::default());
        }
        Err(_) => {
            eprintln!("genotyping failed for {name}: panic while collecting / scoring reads");
            return (output::format_simple_row(name, m.kind, Status::Error, 0, n_prof), LocusStats::default());
        }
    };
    let row = match catch_unwind(AssertUnwindSafe(|| {
        if obs.is_empty() && n_disc == 0 {
            return output::format_simple_row(name, m.kind, Status::NoReads, coverage, n_prof);
        }
        let call = model::call_locus(&obs, n_disc, cfg, cfg.ref_bias_for(m.kind.as_str()));
        output::format_row(name, m.kind, &call, coverage, n_disc)
    })) {
        Ok(row) => row,
        Err(_) => {
            eprintln!("genotyping failed for {name}: panic in the genotype model");
            output::format_simple_row(name, m.kind, Status::Error, coverage, n_prof)
        }
    };
    if let Some(x) = extra {
        x.stats.candidates += 1;
        if passes_alt_gate(&obs, cfg) {
            let t = Instant::now();
            let mut text = String::new();
            let mut budget = x.mate_budget;
            let sample = x.sample;
            let done = catch_unwind(AssertUnwindSafe(|| -> io::Result<extra::ExtraStats> {
                let contig = contig_index(header, &m.locus.chr).ok_or_else(|| io::Error::other("contig not in the header"))?;
                let mut reads = SourceReads::new(&mut *src, header, &m.locus.chr, contig);
                let (anchors, mut st) = collect_extra(m, contig, cfg, sample, &obs, &kept, &mut reads, &mut text)?;
                let mut mates = SourceMates { src: &mut *src, header, light: LightRec::default() };
                fetch_extra_mates(m, sample, &anchors, &mut mates, &mut budget, &mut text, &mut st)?;
                Ok(st)
            }));
            match done {
                Ok(Ok(st)) => {
                    x.stats.add(&st);
                    x.stats.gated += 1;
                    x.mate_budget = budget;
                    x.buf.push_str(&text);
                }
                Ok(Err(e)) => eprintln!("extra reads failed for {name}: {e}"),
                Err(_) => eprintln!("extra reads failed for {name}: panic"),
            }
            x.stats.micros += t.elapsed().as_micros();
        }
    }
    (row, stats)
}

// ------------------------------------------------------------------------------------------
// record access behind small traits (BAM/CRAM in production, vectors in tests)
// ------------------------------------------------------------------------------------------

/// One streamed record: the cheap decoded fields, its name, and on-demand seq/qual.
pub(crate) trait RecordView {
    fn light(&self) -> &LightRec;
    fn qname(&self) -> &[u8];
    /// Decode sequence and phred qualities into the (reused) buffers.
    fn seq_qual(&self, seq: &mut Vec<u8>, qual: &mut Vec<u8>) -> io::Result<()>;
}

/// Streams the records overlapping a 0-based half-open interval of the locus contig, in
/// coordinate order, to a visitor that returns `Ok(false)` to stop early.
pub(crate) trait LocusReads {
    fn visit(&mut self, start: i64, end: i64, f: &mut dyn FnMut(&dyn RecordView) -> io::Result<bool>) -> io::Result<()>;
}

struct SourceReads<'s> {
    src: &'s mut dyn RegionSource,
    header: &'s Header,
    chr: &'s str,
    contig_len: Option<i64>,
    light: LightRec,
}

impl<'s> SourceReads<'s> {
    fn new(src: &'s mut dyn RegionSource, header: &'s Header, chr: &'s str, contig: usize) -> Self {
        let contig_len = header.reference_sequences().get_index(contig).map(|(_, rs)| usize::from(rs.length()) as i64);
        SourceReads { src, header, chr, contig_len, light: LightRec::default() }
    }
}

struct BoundRecord<'a> {
    light: &'a LightRec,
    rec: &'a AnyRecord,
}

impl RecordView for BoundRecord<'_> {
    fn light(&self) -> &LightRec {
        self.light
    }
    fn qname(&self) -> &[u8] {
        name(self.rec.as_dyn())
    }
    fn seq_qual(&self, seq: &mut Vec<u8>, qual: &mut Vec<u8>) -> io::Result<()> {
        seq_into(self.rec.as_dyn(), seq);
        qual_into(self.rec.as_dyn(), qual)
    }
}

impl LocusReads for SourceReads<'_> {
    fn visit(&mut self, start: i64, end: i64, f: &mut dyn FnMut(&dyn RecordView) -> io::Result<bool>) -> io::Result<()> {
        let start = start.max(0);
        let end = match self.contig_len {
            Some(len) => end.min(len),
            None => end,
        };
        if end <= start {
            return Ok(());
        }
        // pysam-style [start, end) -> noodles 1-based inclusive [start+1, end]
        let s = Position::new(start as usize + 1).expect("start + 1 >= 1");
        let e = Position::new(end as usize).expect("end >= 1");
        let region = Region::new(self.chr, s..=e);
        let SourceReads { src, header, light, .. } = self;
        for rec in src.query(&region)? {
            let rec = rec?;
            decode_light_into(rec.as_dyn(), header, light)?;
            if !overlaps(light.reference_start, light.reference_end, start, end) {
                continue; // index over-fetch / unmapped placed mates (empty aligned span)
            }
            if !f(&BoundRecord { light, rec: &rec })? {
                break;
            }
        }
        Ok(())
    }
}

// ------------------------------------------------------------------------------------------
// record-level logic
// ------------------------------------------------------------------------------------------

#[inline]
fn overlaps(ref_start: i64, ref_end: i64, start0: i64, end0: i64) -> bool {
    ref_start < end0 && ref_end > start0
}

/// Junction geometry the gates need. `r_junction` = R when the right junction
/// (`genome[..R) ++ INS`) is real, `l_junction` = L when the left one (`INS ++ genome[L..)`) is.
pub(crate) struct Geometry<'m> {
    contig: usize,
    breakpoints: &'m [i64],
    r_junction: Option<i64>,
    l_junction: Option<i64>,
}

impl<'m> Geometry<'m> {
    pub(crate) fn new(model: &'m LocusModel, contig: usize) -> Self {
        let l = &model.locus;
        // `L-oneside_L` (right_open): the real junction is the LEFT one at L;
        // `oneside_R-R` (left_open): the real junction is the RIGHT one at R.
        Geometry {
            contig,
            breakpoints: &model.breakpoints,
            r_junction: (!l.right_open).then_some(l.right_pos),
            l_junction: (!l.left_open).then_some(l.left_pos),
        }
    }
}

/// Counted towards `coverage` / the high-coverage gate: primary, mapped, not 0x400.
#[inline]
pub(crate) fn depth_eligible(r: &LightRec) -> bool {
    !(r.is_unmapped || r.is_secondary || r.is_supplementary || r.is_duplicate)
}

/// MAPQ rule: `min_mapq`, or `min_mapq_clipped` for a read whose soft clip of at least
/// `min_clip_for_lowmapq` bases sits at a junction-facing end (a trailing clip at R, a leading
/// clip at L): such a read's MAPQ reflects only its short aligned flank.
pub(crate) fn passes_mapq(r: &LightRec, geo: &Geometry, cfg: &Config) -> bool {
    if r.mapq >= cfg.min_mapq {
        return true;
    }
    if r.mapq < cfg.min_mapq_clipped {
        return false;
    }
    let min_clip = cfg.min_clip_for_lowmapq.max(1);
    let right = geo
        .r_junction
        .is_some_and(|rp| r.trailing_softclip() >= min_clip && (r.reference_end - rp).abs() <= CLIP_FACING_TOL);
    let left = geo
        .l_junction
        .is_some_and(|lp| r.leading_softclip() >= min_clip && (r.reference_start - lp).abs() <= CLIP_FACING_TOL);
    right || left
}

/// Discordant anchor (SPEC): primary, MAPQ >= `min_mapq`, mate unmapped / other contig /
/// `|tlen| > disc_max_tlen` / same strand, and a FORWARD anchor whose last aligned base lies in
/// `[R - disc_span, R + 5]` or a REVERSE anchor starting in `[L - 5, L + disc_span]` (its mate
/// points into the insertion). Flag/contig gates are applied by the caller (`gate`).
pub(crate) fn is_discordant_anchor(r: &LightRec, geo: &Geometry, cfg: &Config) -> bool {
    anchor_with(r, geo, cfg.min_mapq, cfg.disc_span, cfg.disc_max_tlen)
}

/// `is_discordant_anchor` with explicit MAPQ floor / span (the extra pass uses its own).
pub(crate) fn anchor_with(r: &LightRec, geo: &Geometry, min_mapq: u8, span: i64, max_tlen: i64) -> bool {
    if !r.is_paired || r.mapq < min_mapq {
        return false;
    }
    let abnormal = r.mate_unmapped
        || r.mate_reference_sequence_id != r.reference_sequence_id
        || r.tlen.abs() > max_tlen
        || r.is_reverse == r.mate_reverse;
    if !abnormal {
        return false;
    }
    if r.is_reverse {
        geo.l_junction.is_some_and(|lp| r.reference_start >= lp - 5 && r.reference_start <= lp + span)
    } else {
        let last = r.reference_end - 1;
        geo.r_junction.is_some_and(|rp| last >= rp - span && last <= rp + 5)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Gate {
    /// passes the read gates and overlaps a breakpoint: decode + score
    Evidence,
    /// a discordant anchor (counted in `n_disc` unless its qname is an evidence read)
    Discordant,
    Skip,
}

/// Per-record classification against one locus (pure; dedup happens afterwards).
pub(crate) fn gate(r: &LightRec, geo: &Geometry, cfg: &Config) -> Gate {
    if r.is_unmapped || r.is_secondary || r.is_supplementary || r.is_qcfail || r.is_duplicate {
        return Gate::Skip;
    }
    if r.reference_sequence_id != Some(geo.contig) {
        return Gate::Skip;
    }
    let on_bp = geo.breakpoints.iter().any(|&bp| overlaps(r.reference_start, r.reference_end, bp - 1, bp + 1));
    if on_bp && passes_mapq(r, geo, cfg) {
        Gate::Evidence
    } else if is_discordant_anchor(r, geo, cfg) {
        Gate::Discordant
    } else {
        Gate::Skip
    }
}

/// Lenient-dedup key: the fragment's (start, end, strand). For a pair with both mates on the
/// contig: [min(pos, mate_pos), +|tlen|) and the strand of read 1; otherwise the read's own
/// unclipped span and strand.
pub(crate) fn fragment_key(r: &LightRec) -> (i64, i64, u8) {
    if r.is_paired && !r.mate_unmapped && r.mate_reference_sequence_id == r.reference_sequence_id && r.tlen != 0 && r.mate_start >= 0 {
        let s = r.reference_start.min(r.mate_start);
        let strand = if r.is_first { r.is_reverse } else { !r.is_reverse };
        (s, s + r.tlen.abs(), u8::from(strand))
    } else {
        let s = r.reference_start - r.leading_softclip() as i64;
        let e = r.reference_end + r.trailing_softclip() as i64;
        (s, e, 2 + u8::from(r.is_reverse))
    }
}

/// Per-locus read dedup: by qname across all the locus's windows (first wins); with `lenient`
/// also by fragment (start, end, strand) — catches duplicates the marker missed.
pub(crate) struct Dedup {
    lenient: bool,
    names: FxHashSet<Vec<u8>>,
    frags: FxHashSet<(i64, i64, u8)>,
}

impl Dedup {
    pub(crate) fn new(lenient: bool) -> Self {
        Dedup { lenient, names: FxHashSet::default(), frags: FxHashSet::default() }
    }
    /// True when this read is the first of its fragment. The qname is remembered either way,
    /// so the mate of a lenient-rejected duplicate is rejected too.
    pub(crate) fn admit(&mut self, qname: &[u8], r: &LightRec) -> bool {
        if self.names.contains(qname) {
            return false;
        }
        self.names.insert(qname.to_vec());
        !self.lenient || self.frags.insert(fragment_key(r))
    }
    pub(crate) fn has_name(&self, qname: &[u8]) -> bool {
        self.names.contains(qname)
    }
}

/// Diagnostics of one locus's record stream.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct LocusStats {
    /// records streamed (all queries)
    pub fetched: usize,
    /// records passing the evidence gates with non-empty seq/qual (before dedup)
    pub gated: usize,
    /// windows whose widened query hit the record cap and were re-queried narrow
    pub requeried: usize,
}

#[derive(Debug)]
pub(crate) enum Outcome {
    HighCoverage,
    Collected { coverage: i64, obs: Vec<ReadObs>, n_disc: i64, stats: LocusStats },
}

/// Stream every fetch window of `model` once and collect: the depth (with the high-coverage
/// early exit, which also ends the evidence collection), the scored evidence reads (gated,
/// deduped), and the discordant anchors.
///
/// Query per window `(lo, hi)`: `[lo - 1 - ext, hi + 1 + ext)` with `ext = max(disc_span, 5)` —
/// the breakpoint window of the legacy fetch, widened so the discordant anchors (which end up
/// to `disc_span` bp before R / start up to `disc_span` after L, without overlapping the
/// breakpoint) come out of the SAME record stream. Depth counts only records overlapping
/// `[lo, hi + 1)` (the legacy count window); evidence only records overlapping a breakpoint.
#[cfg_attr(not(test), allow(dead_code))] // production calls collect_locus_keep
pub(crate) fn collect_locus(
    model: &LocusModel,
    contig: usize,
    cfg: &Config,
    reads: &mut dyn LocusReads,
    score: &mut dyn FnMut(&ReadInput) -> ReadObs,
) -> io::Result<Outcome> {
    collect_locus_keep(model, contig, cfg, reads, score, None)
}

/// `collect_locus` that also keeps, index-aligned with `obs`, what the extra pass needs of each
/// scored read (`keep`, only for loci that may get an extra pass; None = the plain path).
pub(crate) fn collect_locus_keep(
    model: &LocusModel,
    contig: usize,
    cfg: &Config,
    reads: &mut dyn LocusReads,
    score: &mut dyn FnMut(&ReadInput) -> ReadObs,
    mut keep: Option<&mut Vec<KeptRead>>,
) -> io::Result<Outcome> {
    let geo = Geometry::new(model, contig);
    let thr = cfg.reads_for_high_coverage;
    let ext = cfg.disc_span.max(5);
    let cap = (thr.max(0) as usize + 1).saturating_mul(WIDE_QUERY_CAP_FACTOR).max(WIDE_QUERY_CAP_MIN);

    let mut dedup = Dedup::new(cfg.lenient_dedup);
    let mut disc: FxHashSet<Vec<u8>> = FxHashSet::default();
    let mut obs: Vec<ReadObs> = Vec::new();
    let mut stats = LocusStats::default();
    let mut coverage: i64 = 0;
    let (mut seq, mut qual) = (Vec::new(), Vec::new());

    for &(lo, hi) in &model.windows {
        let (d0, d1) = (lo, hi + 1);
        let mut wcov: i64 = 0;
        for (qs, qe, wide) in [(lo - 1 - ext, hi + 1 + ext, true), (lo - 1, hi + 1, false)] {
            wcov = 0;
            let mut streamed = 0usize;
            let mut tripped = false;
            let mut capped = false;
            reads.visit(qs, qe, &mut |rv| {
                streamed += 1;
                if wide && streamed > cap {
                    capped = true;
                    return Ok(false);
                }
                stats.fetched += 1;
                let r = rv.light();
                if depth_eligible(r) && overlaps(r.reference_start, r.reference_end, d0, d1) {
                    wcov += 1;
                    if wcov > thr {
                        tripped = true;
                        return Ok(false);
                    }
                }
                match gate(r, &geo, cfg) {
                    Gate::Evidence => {
                        rv.seq_qual(&mut seq, &mut qual)?;
                        if seq.is_empty() || qual.len() != seq.len() {
                            return Ok(true);
                        }
                        stats.gated += 1;
                        if dedup.admit(rv.qname(), r) {
                            let input = ReadInput {
                                seq: &seq,
                                qual: &qual,
                                cigar: &r.cigar,
                                ref_start: r.reference_start,
                                ref_end: r.reference_end,
                                reverse: r.is_reverse,
                            };
                            obs.push(score(&input));
                            if let Some(k) = keep.as_deref_mut() {
                                k.push(KeptRead {
                                    qname: rv.qname().to_vec(),
                                    seq: seq.clone(),
                                    lead: r.leading_softclip(),
                                    trail: r.trailing_softclip(),
                                    is_first: r.is_first,
                                });
                            }
                        }
                    }
                    Gate::Discordant => {
                        let q = rv.qname();
                        if !disc.contains(q) {
                            disc.insert(q.to_vec());
                        }
                    }
                    Gate::Skip => {}
                }
                Ok(true)
            })?;
            if tripped {
                return Ok(Outcome::HighCoverage);
            }
            if !capped {
                break; // the widened query completed: no narrow re-query
            }
            // pileup next to the breakpoint: re-stream only the breakpoint window. Reads
            // already scored are rejected by the qname dedup; depth is recounted from zero.
            stats.requeried += 1;
        }
        coverage += wcov;
    }
    let n_disc = disc.iter().filter(|q| !dedup.has_name(q)).count() as i64;
    Ok(Outcome::Collected { coverage, obs, n_disc, stats })
}

// ------------------------------------------------------------------------------------------
// extra evidence for undiscovered carriers (`gt_extra_reads`, extra.rs)
// ------------------------------------------------------------------------------------------

/// What the extra pass needs of a scored read (kept index-aligned with `obs`).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct KeptRead {
    pub qname: Vec<u8>,
    /// as stored (reference-forward = site-forward)
    pub seq: Vec<u8>,
    pub lead: usize,
    pub trail: usize,
    pub is_first: bool,
}

/// A fetched mate: its stored sequence and strand bit.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct MateRec {
    pub seq: Vec<u8>,
    pub is_reverse: bool,
}

/// Random access to one read of a pair by (contig, 0-based position, qname, which read).
pub(crate) trait MateSource {
    fn fetch_mate(&mut self, contig: usize, pos: i64, qname: &[u8], want_first: bool) -> io::Result<Option<MateRec>>;
}

/// Records scanned at a mate's position before giving up (a pileup there cannot wedge a locus).
const MATE_SCAN_CAP: usize = 5000;

struct SourceMates<'s> {
    src: &'s mut dyn RegionSource,
    header: &'s Header,
    light: LightRec,
}

impl MateSource for SourceMates<'_> {
    fn fetch_mate(&mut self, contig: usize, pos: i64, qname: &[u8], want_first: bool) -> io::Result<Option<MateRec>> {
        let Some((chr, _)) = self.header.reference_sequences().get_index(contig) else { return Ok(None) };
        let chr = String::from_utf8_lossy(chr.as_ref()).into_owned();
        let p = Position::new(pos.max(0) as usize + 1).expect("pos + 1 >= 1");
        let region = Region::new(chr.as_str(), p..=p);
        let SourceMates { src, header, light } = self;
        // raw query, NOT `SourceReads::visit`: an unmapped mate is placed at its anchor's
        // position with an empty aligned span, which the overlap filter there would drop
        for (n, rec) in src.query(&region)?.enumerate() {
            if n >= MATE_SCAN_CAP {
                break;
            }
            let rec = rec?;
            if name(rec.as_dyn()) != qname {
                continue;
            }
            decode_light_into(rec.as_dyn(), header, light)?;
            if light.is_secondary || light.is_supplementary || (light.is_paired && light.is_first != want_first) {
                continue;
            }
            let mut seq = Vec::new();
            seq_into(rec.as_dyn(), &mut seq);
            return Ok(Some(MateRec { seq, is_reverse: light.is_reverse }));
        }
        Ok(None)
    }
}

/// What `run` needs for the extra pass: this colony's discovery memberships (loci), its name
/// in the sidecar headers, and the sidecar writer.
pub struct ExtraArgs<'a> {
    pub members: &'a FxHashSet<String>,
    pub sample: &'a str,
    pub out: &'a mut dyn Write,
}

/// Per-chunk state of the extra pass: who we are, which loci we discovered, the per-colony
/// mate-fetch budget, counters and the pending sidecar text (drained after every locus).
pub(crate) struct ExtraState<'a> {
    pub members: &'a FxHashSet<String>,
    pub sample: &'a str,
    pub mate_budget: usize,
    pub stats: extra::ExtraStats,
    pub buf: String,
}

/// Side of an Alt read: the junction of its best alt segment; for ALT_FULL (both) the end that
/// carries the longer soft clip (a leading clip faces L, a trailing one R).
pub(crate) fn alt_read_side(side: AltSide, lead: usize, trail: usize) -> &'static str {
    match side {
        AltSide::Left => "LEFT",
        AltSide::Right => "RIGHT",
        _ if trail > lead => "RIGHT",
        _ => "LEFT",
    }
}

/// The membership check: the extra pass is for loci this colony did NOT discover (its own
/// discovery reads are already in insertions.reads.fa.gz).
pub(crate) fn is_extra_candidate(members: &FxHashSet<String>, locus: &str) -> bool {
    !members.contains(locus)
}

/// The ALT gate: enough reads the realigner assigned to the alt haplotype.
pub(crate) fn passes_alt_gate(obs: &[ReadObs], cfg: &Config) -> bool {
    obs.iter().filter(|o| o.class == ReadClass::Alt).count() >= cfg.gt_extra_min_alt
}

/// A discordant anchor of the extra pass, waiting for its mate fetch.
#[derive(Clone, Debug)]
pub(crate) struct ExtraAnchor {
    qname: Vec<u8>,
    rec: LightRec,
    side: &'static str,
    seq: Vec<u8>,
    /// the read is already written as GT_CLIP / GT_POLYA (no GT_DISC record; the mate is fetched)
    written: bool,
}

/// The extra pass of one gated locus, part 1: the kept Alt reads (GT_CLIP / GT_POLYA), then ONE
/// query per fetch window widened by `gt_extra_disc_span` for the discordant anchors (at most
/// `gt_extra_max_mates`, returned for `fetch_extra_mates`). Records go to `out`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn collect_extra(
    model: &LocusModel,
    contig: usize,
    cfg: &Config,
    sample: &str,
    obs: &[ReadObs],
    kept: &[KeptRead],
    reads: &mut dyn LocusReads,
    out: &mut String,
) -> io::Result<(Vec<ExtraAnchor>, extra::ExtraStats)> {
    let mut st = extra::ExtraStats::default();
    let locus = model.locus.name.as_str();
    let mut written: FxHashSet<Vec<u8>> = FxHashSet::default();
    for (o, k) in obs.iter().zip(kept) {
        if o.class != ReadClass::Alt || st.clip + st.polya >= cfg.gt_extra_max_reads {
            continue;
        }
        let side = alt_read_side(o.alt_side, k.lead, k.trail);
        let n = k.seq.len();
        let clip: &[u8] = if side == "LEFT" { &k.seq[..k.lead.min(n)] } else { &k.seq[n - k.trail.min(n)..] };
        let role = if extra::is_pure_polya(clip) {
            st.polya += 1;
            extra::ROLE_POLYA
        } else {
            st.clip += 1;
            extra::ROLE_CLIP
        };
        extra::push_record(out, locus, side, role, sample, &extra::frag_id(&k.qname), if k.is_first { 1 } else { 2 }, &k.seq);
        written.insert(k.qname.clone());
    }

    // discordant anchors: the evidence reads' flag gates, own MAPQ floor and span
    let geo = Geometry::new(model, contig);
    let span = cfg.gt_extra_disc_span.max(5);
    let cap = (cfg.reads_for_high_coverage.max(0) as usize + 1).saturating_mul(WIDE_QUERY_CAP_FACTOR).max(WIDE_QUERY_CAP_MIN);
    let mut dedup = Dedup::new(cfg.lenient_dedup);
    let mut anchors: Vec<ExtraAnchor> = Vec::new();
    let (mut seq, mut qual) = (Vec::new(), Vec::new());
    for &(lo, hi) in &model.windows {
        let mut streamed = 0usize;
        reads.visit(lo - 1 - span, hi + 1 + span, &mut |rv| {
            streamed += 1;
            if streamed > cap || anchors.len() >= cfg.gt_extra_max_mates {
                return Ok(false);
            }
            let r = rv.light();
            if r.is_unmapped || r.is_secondary || r.is_supplementary || r.is_qcfail || r.is_duplicate {
                return Ok(true);
            }
            if r.reference_sequence_id != Some(contig) || !anchor_with(r, &geo, cfg.gt_extra_anchor_mapq, span, cfg.disc_max_tlen) {
                return Ok(true);
            }
            if !dedup.admit(rv.qname(), r) {
                return Ok(true);
            }
            rv.seq_qual(&mut seq, &mut qual)?;
            let q = rv.qname();
            anchors.push(ExtraAnchor {
                qname: q.to_vec(),
                rec: r.clone(),
                // a reverse anchor right of L points into the 3' end (LEFT), a forward one left of R
                side: if r.is_reverse { "LEFT" } else { "RIGHT" },
                seq: seq.clone(),
                written: written.contains(q),
            });
            Ok(true)
        })?;
    }
    Ok((anchors, st))
}

/// The extra pass, part 2 (after the anchor stream, so it may reuse the reader): GT_DISC for each
/// anchor not already written, and its inside mate fetched by random access (GT_MATE), while the
/// per-colony `budget` lasts.
pub(crate) fn fetch_extra_mates(
    model: &LocusModel,
    sample: &str,
    anchors: &[ExtraAnchor],
    mates: &mut dyn MateSource,
    budget: &mut usize,
    out: &mut String,
    st: &mut extra::ExtraStats,
) -> io::Result<()> {
    let locus = model.locus.name.as_str();
    for a in anchors {
        let frag = extra::frag_id(&a.qname);
        let r12 = if a.rec.is_first { 1 } else { 2 };
        if !a.written && !a.seq.is_empty() {
            extra::push_record(out, locus, a.side, extra::ROLE_DISC, sample, &frag, r12, &a.seq);
            st.disc += 1;
        }
        // an unmapped mate sits at its anchor's position (its mate fields = the anchor's own)
        let (Some(mc), true) = (a.rec.mate_reference_sequence_id, a.rec.mate_start >= 0) else {
            st.mates_missing += 1;
            continue;
        };
        if *budget == 0 {
            st.mates_capped += 1;
            continue;
        }
        *budget -= 1;
        match mates.fetch_mate(mc, a.rec.mate_start, &a.qname, !a.rec.is_first)? {
            Some(m) if !m.seq.is_empty() => {
                let s = extra::mate_site_forward(&m.seq, m.is_reverse, a.rec.is_reverse);
                extra::push_record(out, locus, a.side, extra::ROLE_MATE, sample, &frag, 3 - r12, &s);
                st.mates += 1;
            }
            _ => st.mates_missing += 1,
        }
    }
    Ok(())
}

// ------------------------------------------------------------------------------------------
// tests
// ------------------------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{LocusKind, ReadClass};

    /// Locus name -> Locus (test-local parser; contract.rs is another package's stub).
    fn parse(name: &str) -> Locus {
        let (chr, pos) = name.rsplit_once(':').unwrap();
        let (l, r) = pos.rsplit_once('-').unwrap();
        let end = |s: &str| match s.strip_prefix("oneside_") {
            Some(v) => (v.parse::<i64>().unwrap(), true),
            None => (s.parse::<i64>().unwrap(), false),
        };
        let ((lp, lo), (rp, ro)) = (end(l), end(r));
        // `L-oneside_L`: the RIGHT end is open; `oneside_R-R`: the LEFT end is open
        Locus { name: name.into(), chr: chr.into(), left_pos: lp, right_pos: rp, left_open: lo, right_open: ro }
    }

    /// Model with the SPEC windows / breakpoints (what haplotype.rs computes), no segments.
    fn model_for(name: &str, flank: i64) -> LocusModel {
        let locus = parse(name);
        let mut m = LocusModel::error(locus.clone(), "x");
        m.error = None;
        let (l, r) = (locus.left_pos, locus.right_pos);
        let (lo, hi) = (l.min(r), l.max(r));
        if locus.is_one_sided() {
            m.windows = vec![(l, l)];
            m.breakpoints = vec![l];
        } else if hi - lo <= flank {
            m.windows = vec![(lo, hi)];
            m.breakpoints = if l == r { vec![l] } else { vec![lo, hi] };
        } else {
            m.windows = vec![(lo, lo), (hi, hi)];
            m.breakpoints = vec![lo, hi];
        }
        m
    }

    fn dummy_obs() -> ReadObs {
        ReadObs { ll_ref: 0.0, ll_alt: 0.0, class: ReadClass::Uninformative, explained_frac: 1.0, crosses_junction: false, alt_side: crate::types::AltSide::None }
    }

    /// A mapped, primary, properly paired forward read on contig 0.
    fn rec(start: i64, cigar: &[(u8, usize)], mapq: u8) -> LightRec {
        let rl: i64 = cigar.iter().filter(|(op, _)| matches!(op, 0 | 2 | 3 | 7 | 8)).map(|&(_, l)| l as i64).sum();
        LightRec {
            reference_sequence_id: Some(0),
            reference_start: start,
            reference_end: start + rl,
            mapq,
            is_paired: true,
            is_first: true,
            mate_reverse: true,
            mate_reference_sequence_id: Some(0),
            mate_start: start + 200,
            tlen: 350,
            cigar: cigar.to_vec(),
            ..Default::default()
        }
    }

    struct TRec {
        light: LightRec,
        name: Vec<u8>,
        seq: Vec<u8>,
    }
    impl RecordView for TRec {
        fn light(&self) -> &LightRec {
            &self.light
        }
        fn qname(&self) -> &[u8] {
            &self.name
        }
        fn seq_qual(&self, seq: &mut Vec<u8>, qual: &mut Vec<u8>) -> io::Result<()> {
            seq.clear();
            seq.extend_from_slice(&self.seq);
            qual.clear();
            qual.extend(std::iter::repeat(30).take(self.seq.len()));
            Ok(())
        }
    }
    /// In-memory "BAM": records sorted by start; `visit` yields the overlapping ones.
    struct VecReads {
        recs: Vec<TRec>,
        queries: Vec<(i64, i64)>,
        yielded: usize,
    }
    impl VecReads {
        fn new(mut v: Vec<(&str, LightRec)>) -> Self {
            v.sort_by_key(|(_, r)| r.reference_start);
            let recs = v
                .into_iter()
                .map(|(n, light)| {
                    let len: usize = light.cigar.iter().filter(|(op, _)| matches!(op, 0 | 1 | 4 | 7 | 8)).map(|&(_, l)| l).sum();
                    TRec { light, name: n.as_bytes().to_vec(), seq: vec![b'A'; len] }
                })
                .collect();
            VecReads { recs, queries: Vec::new(), yielded: 0 }
        }
    }
    impl LocusReads for VecReads {
        fn visit(&mut self, start: i64, end: i64, f: &mut dyn FnMut(&dyn RecordView) -> io::Result<bool>) -> io::Result<()> {
            self.queries.push((start, end));
            for r in &self.recs {
                if overlaps(r.light.reference_start, r.light.reference_end, start, end) {
                    self.yielded += 1;
                    if !f(r)? {
                        break;
                    }
                }
            }
            Ok(())
        }
    }

    fn run_collect(m: &LocusModel, cfg: &Config, reads: &mut VecReads) -> (Outcome, usize) {
        let mut n_scored = 0usize;
        let out = collect_locus(m, 0, cfg, reads, &mut |_| {
            n_scored += 1;
            dummy_obs()
        })
        .unwrap();
        (out, n_scored)
    }

    #[test]
    fn mapq_rule_rescues_only_junction_facing_clips() {
        let cfg = Config::default(); // min_mapq 60, clipped 20, min clip 20
        let m = model_for("chr1:1000-1015", 300); // TSD: L=1000, R=1015
        let geo = Geometry::new(&m, 0);
        // high MAPQ: always
        assert!(passes_mapq(&rec(900, &[(0, 150)], 60), &geo, &cfg));
        // low MAPQ, no clip: no
        assert!(!passes_mapq(&rec(900, &[(0, 150)], 30), &geo, &cfg));
        // trailing 30-bp clip ending at R (alt read of the right junction): yes
        assert!(passes_mapq(&rec(895, &[(0, 120), (4, 30)], 30), &geo, &cfg));
        // ... within the tolerance (aligner ran 4 bp into the insertion): yes
        assert!(passes_mapq(&rec(899, &[(0, 120), (4, 30)], 30), &geo, &cfg));
        // ... but MAPQ below min_mapq_clipped: no
        assert!(!passes_mapq(&rec(895, &[(0, 120), (4, 30)], 19), &geo, &cfg));
        // ... clip too short: no
        assert!(!passes_mapq(&rec(895, &[(0, 135), (4, 15)], 30), &geo, &cfg));
        // leading clip starting at L (alt read of the left junction): yes
        assert!(passes_mapq(&rec(1000, &[(4, 40), (0, 110)], 25), &geo, &cfg));
        // leading clip at R / trailing clip at L face the wrong way: no
        assert!(!passes_mapq(&rec(1015, &[(4, 40), (0, 110)], 25), &geo, &cfg));
        assert!(!passes_mapq(&rec(880, &[(0, 120), (4, 30)], 25), &geo, &cfg));
        // a clip far from any junction: no
        assert!(!passes_mapq(&rec(700, &[(0, 120), (4, 30)], 30), &geo, &cfg));
        // one-sided `L-oneside_L` (real LEFT junction): leading clip at L yes, trailing at L no
        let one = model_for("chr1:1000-oneside_1000", 300);
        let g1 = Geometry::new(&one, 0);
        assert!(passes_mapq(&rec(1000, &[(4, 40), (0, 110)], 25), &g1, &cfg));
        assert!(!passes_mapq(&rec(880, &[(0, 120), (4, 30)], 25), &g1, &cfg));
        // `oneside_R-R` (real RIGHT junction): trailing clip at R yes
        let oner = model_for("chr1:oneside_1000-1000", 300);
        let g2 = Geometry::new(&oner, 0);
        assert!(passes_mapq(&rec(880, &[(0, 120), (4, 30)], 25), &g2, &cfg));
        assert!(!passes_mapq(&rec(1000, &[(4, 40), (0, 110)], 25), &g2, &cfg));
    }

    #[test]
    fn discordant_anchor_classification() {
        let cfg = Config::default(); // disc_span 300, disc_max_tlen 1000
        let m = model_for("chr1:5000-5010", 300); // L=5000, R=5010
        let geo = Geometry::new(&m, 0);
        let fwd = |end_excl: i64| rec(end_excl - 150, &[(0, 150)], 60);
        let rev = |start: i64| {
            let mut r = rec(start, &[(0, 150)], 60);
            r.is_reverse = true;
            r.mate_reverse = false;
            r
        };
        // a normal pair is never discordant
        assert!(!is_discordant_anchor(&fwd(4900), &geo, &cfg));
        // forward anchor ending near R, mate unmapped / other contig / far / same strand
        let mut a = fwd(4900);
        a.mate_unmapped = true;
        assert!(is_discordant_anchor(&a, &geo, &cfg));
        let mut b = fwd(4900);
        b.mate_reference_sequence_id = Some(3);
        b.tlen = 0;
        assert!(is_discordant_anchor(&b, &geo, &cfg));
        let mut c = fwd(4900);
        c.tlen = 5000;
        assert!(is_discordant_anchor(&c, &geo, &cfg));
        let mut d = fwd(4900);
        d.mate_reverse = false; // same strand
        assert!(is_discordant_anchor(&d, &geo, &cfg));
        // forward window [R - 300, R + 5] on the LAST aligned base (end_excl - 1)
        let mut e = fwd(5010 - 300 + 1);
        e.mate_unmapped = true;
        assert!(is_discordant_anchor(&e, &geo, &cfg));
        let mut e2 = fwd(5010 - 300);
        e2.mate_unmapped = true;
        assert!(!is_discordant_anchor(&e2, &geo, &cfg));
        let mut e3 = fwd(5010 + 6);
        e3.mate_unmapped = true;
        assert!(is_discordant_anchor(&e3, &geo, &cfg));
        let mut e4 = fwd(5010 + 7);
        e4.mate_unmapped = true;
        assert!(!is_discordant_anchor(&e4, &geo, &cfg));
        // reverse anchor starting in [L - 5, L + 300]
        let mut f = rev(5100);
        f.mate_unmapped = true;
        assert!(is_discordant_anchor(&f, &geo, &cfg));
        let mut f2 = rev(4995);
        f2.mate_unmapped = true;
        assert!(is_discordant_anchor(&f2, &geo, &cfg));
        let mut f3 = rev(4994);
        f3.mate_unmapped = true;
        assert!(!is_discordant_anchor(&f3, &geo, &cfg));
        let mut f4 = rev(5301);
        f4.mate_unmapped = true;
        assert!(!is_discordant_anchor(&f4, &geo, &cfg));
        // wrong orientation: a forward read right of L / reverse read left of R
        let mut g = rec(5100, &[(0, 150)], 60);
        g.mate_unmapped = true;
        assert!(!is_discordant_anchor(&g, &geo, &cfg));
        let mut h = rev(4700);
        h.mate_unmapped = true;
        assert!(!is_discordant_anchor(&h, &geo, &cfg));
        // MAPQ below min_mapq, unpaired: no
        let mut i = fwd(4900);
        i.mate_unmapped = true;
        i.mapq = 59;
        assert!(!is_discordant_anchor(&i, &geo, &cfg));
        let mut j = fwd(4900);
        j.mate_unmapped = true;
        j.is_paired = false;
        assert!(!is_discordant_anchor(&j, &geo, &cfg));
        // one-sided `L-oneside_L`: only the LEFT junction anchors (reverse, right of L)
        let one = model_for("chr1:5000-oneside_5000", 300);
        let g1 = Geometry::new(&one, 0);
        assert!(is_discordant_anchor(&f, &g1, &cfg));
        let mut k = fwd(4950);
        k.mate_unmapped = true;
        assert!(!is_discordant_anchor(&k, &g1, &cfg));
        // gate(): duplicates and secondary anchors are skipped
        let mut dup = a.clone();
        dup.is_duplicate = true;
        assert_eq!(gate(&a, &geo, &cfg), Gate::Discordant);
        assert_eq!(gate(&dup, &geo, &cfg), Gate::Skip);
    }

    #[test]
    fn gate_flags_contig_overlap() {
        let cfg = Config::default();
        let m = model_for("chr1:1000-1015", 300);
        let geo = Geometry::new(&m, 0);
        let ok = rec(950, &[(0, 150)], 60);
        assert_eq!(gate(&ok, &geo, &cfg), Gate::Evidence);
        for f in [
            |r: &mut LightRec| r.is_secondary = true,
            |r: &mut LightRec| r.is_supplementary = true,
            |r: &mut LightRec| r.is_unmapped = true,
            |r: &mut LightRec| r.is_qcfail = true,
            |r: &mut LightRec| r.is_duplicate = true,
            |r: &mut LightRec| r.reference_sequence_id = Some(1),
            |r: &mut LightRec| r.mapq = 59,
        ] {
            let mut r = ok.clone();
            f(&mut r);
            assert_eq!(gate(&r, &geo, &cfg), Gate::Skip);
        }
        // between the breakpoints of a TSD but touching neither [L-1, L+1) nor [R-1, R+1)
        let inside = rec(1002, &[(0, 10)], 60);
        assert_eq!(gate(&inside, &geo, &cfg), Gate::Skip);
        // ends exactly at L-1 (covers L-1): overlaps [L-1, L+1)
        let edge = rec(850, &[(0, 150)], 60);
        assert_eq!(gate(&edge, &geo, &cfg), Gate::Evidence);
        let before = rec(849, &[(0, 150)], 60);
        assert_eq!(gate(&before, &geo, &cfg), Gate::Skip);
    }

    #[test]
    fn dedup_by_qname_and_lenient_fragment() {
        let a = rec(1000, &[(0, 150)], 60); // pair [1000, 1350), read1 forward
        let mut mate = rec(1200, &[(0, 150)], 60);
        mate.is_first = false;
        mate.is_reverse = true;
        mate.mate_start = 1000;
        mate.tlen = -350;
        let dup_other_name = a.clone();
        let mut strict = Dedup::new(false);
        assert!(strict.admit(b"r1", &a));
        assert!(!strict.admit(b"r1", &mate)); // overlapping mate counts once
        assert!(strict.admit(b"r2", &dup_other_name)); // unmarked PCR dup survives strict
        assert_eq!(fragment_key(&a), fragment_key(&mate)); // both mates name one fragment
        let mut lenient = Dedup::new(true);
        assert!(lenient.admit(b"r1", &a));
        assert!(!lenient.admit(b"r2", &dup_other_name)); // same start/end/strand: dup
        let mut r2mate = mate.clone();
        r2mate.mate_start = 1000;
        assert!(!lenient.admit(b"r2", &r2mate)); // its mate: rejected by qname
        assert!(lenient.has_name(b"r2"));
        // same start, other strand of the fragment: distinct
        let mut flipped = a.clone();
        flipped.is_reverse = true;
        flipped.mate_reverse = false;
        assert!(lenient.admit(b"r3", &flipped));
        // different fragment end: distinct
        let mut longer = a.clone();
        longer.tlen = 400;
        assert!(lenient.admit(b"r4", &longer));
        // mate unmapped: the read's own unclipped span
        let mut u1 = rec(2000, &[(4, 20), (0, 130)], 60);
        u1.mate_unmapped = true;
        let mut u2 = rec(2000, &[(4, 20), (0, 130)], 60);
        u2.mate_unmapped = true;
        assert_eq!(fragment_key(&u1), (1980, 2130, 2));
        assert!(lenient.admit(b"u1", &u1));
        assert!(!lenient.admit(b"u2", &u2));
    }

    #[test]
    fn collect_counts_depth_evidence_disc_and_dedups() {
        let cfg = Config::default();
        let m = model_for("chr1:1000-1015", 300);
        let mut disc_mate_unmapped = rec(700, &[(0, 150)], 60); // ends at 849 (>= R-300=715)
        disc_mate_unmapped.mate_unmapped = true;
        let mut dup = rec(960, &[(0, 150)], 60);
        dup.is_duplicate = true;
        let mut mate = rec(990, &[(0, 150)], 60);
        mate.is_reverse = true;
        let mut lowq = rec(960, &[(0, 150)], 10);
        lowq.mate_start = 0;
        let mut disc_also_evidence = rec(900, &[(0, 115), (4, 35)], 60); // ends at R: evidence
        disc_also_evidence.mate_unmapped = true;
        let mut reads = VecReads::new(vec![
            ("d1", disc_mate_unmapped),
            ("a", rec(950, &[(0, 150)], 60)),
            ("dup", dup),
            ("a", mate), // mate of "a": deduped
            ("low", lowq),
            ("j", disc_also_evidence),
            ("far", rec(3000, &[(0, 150)], 60)), // outside the widened query
        ]);
        let (out, n_scored) = run_collect(&m, &cfg, &mut reads);
        let Outcome::Collected { coverage, obs, n_disc, stats } = out else { panic!("high coverage") };
        // depth window [1000, 1016): a, a-mate, low, j (dup and d1 excluded)
        assert_eq!(coverage, 4);
        // evidence: a, j (a-mate deduped, low MAPQ dropped, dup dropped)
        assert_eq!((obs.len(), n_scored), (2, 2));
        assert_eq!(stats.gated, 3);
        // disc: d1 only ("j" is an evidence read)
        assert_eq!(n_disc, 1);
        // ONE query for the single window, widened by disc_span
        assert_eq!(reads.queries, vec![(999 - 300, 1016 + 300)]);
        assert_eq!(stats.fetched, 6);
    }

    #[test]
    fn depth_early_exit_stops_evidence_too() {
        let cfg = Config { reads_for_high_coverage: 5, ..Config::default() };
        let m = model_for("chr1:1000-1015", 300);
        let v: Vec<(String, LightRec)> = (0..50).map(|i| (format!("r{i}"), rec(900 + i, &[(0, 150)], 60))).collect();
        let mut reads = VecReads::new(v.iter().map(|(n, r)| (n.as_str(), r.clone())).collect());
        let (out, n_scored) = run_collect(&m, &cfg, &mut reads);
        assert!(matches!(out, Outcome::HighCoverage));
        assert_eq!(reads.yielded, 6, "stream stops at threshold + 1");
        assert_eq!(n_scored, 5, "no read scored after the gate trips");
        // exactly at the threshold: no trip, all counted
        let v: Vec<(String, LightRec)> = (0..5).map(|i| (format!("r{i}"), rec(900 + i, &[(0, 150)], 60))).collect();
        let mut reads = VecReads::new(v.iter().map(|(n, r)| (n.as_str(), r.clone())).collect());
        let (out, _) = run_collect(&m, &cfg, &mut reads);
        assert!(matches!(out, Outcome::Collected { coverage: 5, .. }));
        // duplicates do not count towards the gate
        let v: Vec<(String, LightRec)> = (0..50)
            .map(|i| {
                let mut r = rec(900 + i, &[(0, 150)], 60);
                r.is_duplicate = i >= 3;
                (format!("r{i}"), r)
            })
            .collect();
        let mut reads = VecReads::new(v.iter().map(|(n, r)| (n.as_str(), r.clone())).collect());
        let (out, _) = run_collect(&m, &cfg, &mut reads);
        assert!(matches!(out, Outcome::Collected { coverage: 3, .. }));
    }

    #[test]
    fn pileup_beside_the_breakpoint_falls_back_to_the_narrow_window() {
        let cfg = Config { reads_for_high_coverage: 5, ..Config::default() }; // cap = max(2000, 60)
        let m = model_for("chr1:1000-1015", 300);
        // 3000 reads stacked at 750..850 (inside the widened query, not over the window)
        let mut v: Vec<(String, LightRec)> = (0..3000).map(|i| (format!("p{i}"), rec(700, &[(0, 100)], 60))).collect();
        v.push(("x".into(), rec(950, &[(0, 150)], 60)));
        v.push(("y".into(), rec(990, &[(0, 150)], 60)));
        let mut reads = VecReads::new(v.iter().map(|(n, r)| (n.as_str(), r.clone())).collect());
        let (out, n_scored) = run_collect(&m, &cfg, &mut reads);
        let Outcome::Collected { coverage, stats, .. } = out else { panic!() };
        assert_eq!(stats.requeried, 1);
        assert_eq!(reads.queries, vec![(699, 1316), (999, 1016)]);
        assert_eq!((coverage, n_scored), (2, 2));
    }

    #[test]
    fn split_windows_share_dedup_and_sum_coverage() {
        let cfg = Config::default();
        let m = model_for("chr1:1000-1500", 300); // far duplication: two windows
        assert_eq!(m.windows, vec![(1000, 1000), (1500, 1500)]);
        let mut reads = VecReads::new(vec![
            ("a", rec(950, &[(0, 150)], 60)),
            ("b", rec(1450, &[(0, 150)], 60)),
            ("c", rec(1100, &[(0, 150)], 60)), // in both widened queries, overlaps no breakpoint
        ]);
        let (out, n) = run_collect(&m, &cfg, &mut reads);
        let Outcome::Collected { coverage, .. } = out else { panic!() };
        assert_eq!((coverage, n), (2, 2));
        assert_eq!(reads.queries.len(), 2);
    }

    #[test]
    fn locus_sorting_by_header_order() {
        let names = ["chr2:500-510", "chr1:900-910", "chrUn:5-6", "chr1:100-oneside_100", "chr2:50-30", "chr1:5000-1000"];
        let models: Vec<LocusModel> = names.iter().map(|n| model_for(n, 300)).collect();
        // header order chr2 before chr1
        let idx = |c: &str| match c {
            "chr2" => Some(0),
            "chr1" => Some(1),
            _ => None,
        };
        let order = sort_order(&models, idx);
        let sorted: Vec<&str> = order.iter().map(|&i| names[i]).collect();
        // min(L, R): chr1:5000-1000 sorts at 1000, after chr1:900-910
        assert_eq!(sorted, vec!["chr2:50-30", "chr2:500-510", "chr1:100-oneside_100", "chr1:900-910", "chr1:5000-1000", "chrUn:5-6"]);
        assert_eq!(models[5].kind, LocusKind::FarDeletion);
    }

    #[test]
    fn chunk_partitioning() {
        assert_eq!(partition(10, 3), vec![0..4, 4..7, 7..10]);
        assert_eq!(partition(2, 8), vec![0..1, 1..2]);
        assert_eq!(partition(0, 4), vec![0..0]);
        assert_eq!(partition(5, 1), vec![0..5]);
        assert_eq!(partition(5, 0), vec![0..5]);
        let p = partition(31_000, 7);
        assert_eq!(p.len(), 7);
        assert_eq!(p.last().unwrap().end, 31_000);
        assert!(p.windows(2).all(|w| w[0].end == w[1].start));
    }

    #[test]
    fn heartbeat_cadence() {
        let beats: Vec<usize> = (1..=2500).filter(|&d| is_heartbeat(d, 2500, 1000)).collect();
        assert_eq!(beats, vec![1000, 2000, 2500]);
        assert_eq!((1..=7).filter(|&d| is_heartbeat(d, 7, 0)).collect::<Vec<_>>(), vec![7]);
    }

    #[test]
    fn merged_sides_fill_missing_side_from_contract() {
        use crate::types::JunctionConsensus;
        let jc = |s: &[u8], contract: bool| JunctionConsensus { ins_seq: s.to_vec(), from_contract_only: contract, ..Default::default() };
        let locus = parse("chr1:1000-1015");
        let contract = ContractSides { left: Some(jc(b"AAAA", true)), right: Some(jc(b"CCCC", true)) };
        let mut map = HashMap::new();
        map.insert(locus.name.clone(), ContractSides { left: Some(jc(b"GGGGGGGG", false)), right: None });
        let (s, fb) = merged_sides(&locus, &contract, Some(&map));
        assert!(fb);
        assert_eq!(s.left.unwrap().ins_seq, b"GGGGGGGG");
        assert_eq!(s.right.unwrap().ins_seq, b"CCCC");
        let one = parse("chr1:1000-oneside_1000"); // right open: only the left side is needed
        map.insert(one.name.clone(), ContractSides { left: Some(jc(b"TTTT", false)), right: None });
        assert!(!merged_sides(&one, &ContractSides::default(), Some(&map)).1);
        assert!(!merged_sides(&locus, &contract, None).1);
    }

    #[test]
    fn check_input_requires_index() {
        let dir = std::env::temp_dir().join(format!("pt2_check_input_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("x.bam");
        std::fs::write(&bam, b"").unwrap();
        let p = bam.to_str().unwrap();
        assert!(check_input(p).unwrap_err().contains("not indexed"));
        std::fs::write(dir.join("x.bam.csi"), b"").unwrap();
        assert!(check_input(p).is_ok());
        assert!(check_input(dir.join("nope.bam").to_str().unwrap()).unwrap_err().contains("does not exist"));
        let cram = dir.join("y.cram");
        std::fs::write(&cram, b"").unwrap();
        assert!(check_input(cram.to_str().unwrap()).is_err());
        std::fs::write(dir.join("y.cram.crai"), b"").unwrap();
        assert!(check_input(cram.to_str().unwrap()).is_ok());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // ---------------------------------------------------------------------------------------
    // extra pass (gt_extra_reads)
    // ---------------------------------------------------------------------------------------

    fn alt_obs(side: AltSide) -> ReadObs {
        ReadObs { ll_ref: -50.0, ll_alt: -5.0, class: ReadClass::Alt, explained_frac: 1.0, crosses_junction: true, alt_side: side }
    }

    fn kept(name: &str, seq: &[u8], lead: usize, trail: usize) -> KeptRead {
        KeptRead { qname: name.as_bytes().to_vec(), seq: seq.to_vec(), lead, trail, is_first: true }
    }

    /// Mates by (contig, pos, qname, first?) -> (stored seq, reverse); records every request.
    #[derive(Default)]
    struct FakeMates {
        mates: Vec<((usize, i64, Vec<u8>, bool), MateRec)>,
        asked: Vec<(usize, i64, Vec<u8>, bool)>,
    }
    impl MateSource for FakeMates {
        fn fetch_mate(&mut self, contig: usize, pos: i64, qname: &[u8], want_first: bool) -> io::Result<Option<MateRec>> {
            let k = (contig, pos, qname.to_vec(), want_first);
            self.asked.push(k.clone());
            Ok(self.mates.iter().find(|(m, _)| *m == k).map(|(_, r)| r.clone()))
        }
    }

    fn records(text: &str) -> Vec<(String, String)> {
        let l: Vec<&str> = text.lines().collect();
        l.chunks(2).map(|c| (c[0].to_string(), c[1].to_string())).collect()
    }

    #[test]
    fn extra_membership_and_alt_gate() {
        let mut members = FxHashSet::default();
        members.insert("chr1:1000-1015".to_string());
        assert!(!is_extra_candidate(&members, "chr1:1000-1015"), "a discovery member gets no extra pass");
        assert!(is_extra_candidate(&members, "chr1:2000-2015"));
        assert!(is_extra_candidate(&members, "chr1:1000-1016"), "exact locus names only");
        let cfg = Config::default(); // gt_extra_min_alt 1
        let mut obs = vec![dummy_obs(), dummy_obs()];
        assert!(!passes_alt_gate(&obs, &cfg));
        obs.push(alt_obs(AltSide::Right));
        assert!(passes_alt_gate(&obs, &cfg));
        let two = Config { gt_extra_min_alt: 2, ..Config::default() };
        assert!(!passes_alt_gate(&obs, &two));
        // an Unexplained / Ref read never opens the gate
        obs.push(ReadObs { class: ReadClass::Unexplained, ..alt_obs(AltSide::None) });
        assert!(!passes_alt_gate(&obs, &two));
        assert!(!Config::default().gt_extra_reads, "the key defaults OFF");
    }

    #[test]
    fn collect_keep_is_index_aligned_and_off_by_default() {
        let cfg = Config::default();
        let m = model_for("chr1:1000-1015", 300);
        let mk = || {
            let mut mate = rec(990, &[(0, 150)], 60);
            mate.is_reverse = true;
            VecReads::new(vec![
                ("a", rec(950, &[(0, 150)], 60)),
                ("b", rec(900, &[(0, 115), (4, 35)], 60)),
                ("a", mate), // deduped: not kept either
            ])
        };
        let mut kept = Vec::new();
        let mut reads = mk();
        let out = collect_locus_keep(&m, 0, &cfg, &mut reads, &mut |_| dummy_obs(), Some(&mut kept)).unwrap();
        let Outcome::Collected { obs, .. } = out else { panic!() };
        assert_eq!(obs.len(), kept.len());
        assert_eq!(kept.iter().map(|k| k.qname.as_slice()).collect::<Vec<_>>(), vec![b"b".as_slice(), b"a"]);
        assert_eq!((kept[0].lead, kept[0].trail, kept[0].seq.len()), (0, 35, 150));
        // without `keep` the outcome is the same
        let mut reads = mk();
        let plain = collect_locus(&m, 0, &cfg, &mut reads, &mut |_| dummy_obs()).unwrap();
        let mut reads = mk();
        let with = collect_locus_keep(&m, 0, &cfg, &mut reads, &mut |_| dummy_obs(), Some(&mut Vec::new())).unwrap();
        assert_eq!(format!("{plain:?}"), format!("{with:?}"));
    }

    #[test]
    fn extra_reads_sides_roles_and_anchor_gates() {
        let cfg = Config { gt_extra_max_mates: 3, ..Config::default() }; // span 500, anchor MAPQ 20
        let m = model_for("chr1:5000-5010", 300); // L = 5000, R = 5010
        let polya = [b"ACGTACGTAC".as_slice(), &[b'A'; 20]].concat();
        let obs = vec![alt_obs(AltSide::Right), dummy_obs(), alt_obs(AltSide::Left), alt_obs(AltSide::Both)];
        let kept = vec![
            kept("r1", b"CCCCCCCCCCGGGGG", 0, 5), // trailing clip at R -> RIGHT, GT_CLIP
            kept("ref", b"ACGT", 0, 0),          // not Alt: skipped
            kept("l1", &polya, 20, 0),           // leading clip half A: LEFT, GT_CLIP (not a pure tail)
            kept("f1", b"TTTTTTTTTTTTGGG", 12, 3), // ALT_FULL, longer leading clip -> LEFT, poly-T -> GT_POLYA
        ];
        // anchors: forward ending 450 bp before R (inside span 500, outside the genotyper's 300),
        // MAPQ 25 (>= 20, < min_mapq 60), mate unmapped; a duplicate; a low-MAPQ one; a reverse
        // one right of L with its mate on chr3; a normal pair; one beyond the span
        let fwd = |end_excl: i64, mapq: u8| {
            let mut r = rec(end_excl - 100, &[(0, 100)], mapq);
            r.mate_unmapped = true;
            r.mate_start = end_excl - 100; // placed at the anchor
            r
        };
        let mut rev = rec(5200, &[(0, 100)], 60);
        rev.is_reverse = true;
        rev.mate_reverse = false;
        rev.mate_reference_sequence_id = Some(3);
        rev.mate_start = 777;
        rev.tlen = 0;
        rev.is_first = false;
        let mut dup = fwd(4700, 60);
        dup.is_duplicate = true;
        let mut reads = VecReads::new(vec![
            ("d1", fwd(5010 - 450, 25)),
            ("dup", dup),
            ("low", fwd(4800, 19)),
            ("rv", rev),
            ("normal", rec(4800, &[(0, 100)], 60)),
            ("far", fwd(5010 - 520, 60)),
            ("r1", fwd(5000, 60)), // also written as GT_CLIP: no GT_DISC, mate still fetched
        ]);
        let mut out = String::new();
        let (anchors, mut st) = collect_extra(&m, 0, &cfg, "S9", &obs, &kept, &mut reads, &mut out).unwrap();
        assert_eq!(reads.queries, vec![(5000 - 1 - 500, 5010 + 1 + 500)], "one widened query per window");
        let names: Vec<&[u8]> = anchors.iter().map(|a| a.qname.as_slice()).collect();
        assert_eq!(names, vec![b"d1".as_slice(), b"r1", b"rv"], "dup / low MAPQ / normal / beyond-span skipped");
        assert_eq!(anchors.iter().map(|a| a.side).collect::<Vec<_>>(), vec!["RIGHT", "RIGHT", "LEFT"]);
        let mut mates = FakeMates::default();
        mates.mates.push(((0, 4460, b"d1".to_vec(), false), MateRec { seq: b"AACCGT".to_vec(), is_reverse: false }));
        mates.mates.push(((3, 777, b"rv".to_vec(), true), MateRec { seq: b"AACCGT".to_vec(), is_reverse: false }));
        let mut budget = 100;
        fetch_extra_mates(&m, "S9", &anchors, &mut mates, &mut budget, &mut out, &mut st).unwrap();
        assert_eq!(budget, 97);
        let recs = records(&out);
        let heads: Vec<&str> = recs.iter().map(|(h, _)| h.as_str()).collect();
        let f = |q: &[u8]| extra::frag_id(q);
        assert_eq!(
            heads,
            vec![
                format!(">chr1:5000-5010|RIGHT|GT_CLIP|S9|{}|1", f(b"r1")),
                format!(">chr1:5000-5010|LEFT|GT_CLIP|S9|{}|1", f(b"l1")),
                format!(">chr1:5000-5010|LEFT|GT_POLYA|S9|{}|1", f(b"f1")),
                format!(">chr1:5000-5010|RIGHT|GT_DISC|S9|{}|1", f(b"d1")),
                format!(">chr1:5000-5010|RIGHT|GT_MATE|S9|{}|2", f(b"d1")),
                format!(">chr1:5000-5010|LEFT|GT_DISC|S9|{}|2", f(b"rv")),
                format!(">chr1:5000-5010|LEFT|GT_MATE|S9|{}|1", f(b"rv")),
            ]
        );
        // orientation: forward anchor + unmapped (0x10 clear) mate -> reverse-complemented;
        // reverse anchor + forward-stored mate -> kept as stored
        assert_eq!(recs[4].1, "ACGGTT");
        assert_eq!(recs[6].1, "AACCGT");
        assert_eq!(recs[0].1, "CCCCCCCCCCGGGGG", "junction reads as stored (site-forward)");
        // r1's mate was requested (no record: not in the fake) and counted missing
        assert!(mates.asked.contains(&(0, 4900, b"r1".to_vec(), false)));
        assert_eq!((st.clip, st.polya, st.disc, st.mates, st.mates_missing), (2, 1, 2, 2, 1));
    }

    #[test]
    fn extra_caps_per_locus_and_per_colony() {
        let cfg = Config { gt_extra_max_mates: 2, gt_extra_max_reads: 1, ..Config::default() };
        let m = model_for("chr1:5000-5010", 300);
        let v: Vec<(String, LightRec)> = (0..6)
            .map(|i| {
                let mut r = rec(4800 + i, &[(0, 100)], 60);
                r.mate_unmapped = true;
                r.mate_start = 4800 + i;
                (format!("a{i}"), r)
            })
            .collect();
        let mut reads = VecReads::new(v.iter().map(|(n, r)| (n.as_str(), r.clone())).collect());
        let obs = vec![alt_obs(AltSide::Right), alt_obs(AltSide::Right)];
        let kept = vec![kept("x", b"ACGTACGTAC", 0, 4), kept("y", b"ACGTACGTAC", 0, 4)];
        let mut out = String::new();
        let (anchors, mut st) = collect_extra(&m, 0, &cfg, "S", &obs, &kept, &mut reads, &mut out).unwrap();
        assert_eq!(anchors.len(), 2, "gt_extra_max_mates anchors per locus");
        assert_eq!(st.clip, 1, "gt_extra_max_reads junction reads per locus");
        let mut mates = FakeMates::default();
        let mut budget = 1; // the per-colony fetch budget runs out after one mate
        fetch_extra_mates(&m, "S", &anchors, &mut mates, &mut budget, &mut out, &mut st).unwrap();
        assert_eq!(mates.asked.len(), 1);
        assert_eq!((budget, st.mates_capped, st.disc), (0, 1, 2));
    }

    // ---------------------------------------------------------------------------------------
    // real-BAM smoke test of the I/O path (fetch + gates + dedup, no scoring / no model).
    //   cargo test --release -- --ignored --nocapture smoke_
    // Defaults to test_data/test.bam + its contract; PEARTREE_SMOKE_BAM / PEARTREE_SMOKE_CONTRACT
    // point it at other data. Prints one TSV line per locus (name, coverage or HIGH, fetched,
    // gated, kept, n_disc) and the ms/locus of the whole pass.
    // ---------------------------------------------------------------------------------------

    fn contract_names(path: &str) -> Vec<String> {
        use std::io::{BufRead, BufReader};
        let f = std::fs::File::open(path).unwrap();
        BufReader::new(flate2::read::MultiGzDecoder::new(f))
            .lines()
            .map(|l| l.unwrap())
            .filter_map(|l| l.strip_prefix('>').map(|s| s.trim().to_string()))
            .collect()
    }

    #[test]
    #[ignore]
    fn smoke_io_path_real_bam() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
        let bam = std::env::var("PEARTREE_SMOKE_BAM").unwrap_or_else(|_| format!("{root}/test_data/test.bam"));
        let contract = std::env::var("PEARTREE_SMOKE_CONTRACT")
            .unwrap_or_else(|_| format!("{root}/test_data/test_step2.txt.genotyping.txt.gz"));
        let cfg = Config::load(std::env::var("PEARTREE_SMOKE_CONFIG").ok().as_deref()).unwrap();
        check_input(&bam).unwrap();
        let models: Vec<LocusModel> = contract_names(&contract).iter().map(|n| model_for(n, cfg.flank)).collect();
        let t0 = Instant::now();
        let io0 = io_counters();
        let mut src = open_source_buffered(&bam, None, cfg.io_buffer_bytes, cfg.io_fill_bytes).unwrap();
        let header = src.header().clone();
        let order = sort_order(&models, |c| contig_index(&header, c));
        let mut lines = Vec::new();
        let mut decoded = 0usize;
        for &i in &order {
            let m = &models[i];
            let contig = contig_index(&header, &m.locus.chr).expect("contig in header");
            let mut reads = SourceReads::new(src.as_mut(), &header, &m.locus.chr, contig);
            // "scoring" = touch the decoded read so the decode cost is part of the timing
            let out = collect_locus(m, contig, &cfg, &mut reads, &mut |r| {
                decoded += r.seq.len() + r.qual.len() + r.cigar.len();
                dummy_obs()
            })
            .unwrap();
            lines.push(match out {
                Outcome::HighCoverage => format!("{}\tHIGH\t-\t-\t-\t-", m.locus.name),
                Outcome::Collected { coverage, obs, n_disc, stats } => {
                    format!("{}\t{coverage}\t{}\t{}\t{}\t{n_disc}", m.locus.name, stats.fetched, stats.gated, obs.len())
                }
            });
        }
        let secs = t0.elapsed().as_secs_f64();
        let io1 = io_counters();
        println!("#locus\tcoverage\tfetched\tgated\tkept\tn_disc");
        for l in &lines {
            println!("{l}");
        }
        println!(
            "# {} loci in {:.3}s = {:.3} ms/locus (io_buffer_bytes {}, decoded {} bytes)",
            models.len(),
            secs,
            1000.0 * secs / models.len().max(1) as f64,
            cfg.io_buffer_bytes,
            decoded
        );
        println!(
            "# file I/O: {} reads, {:.1} MiB, {} buffer-miss seeks",
            io1.0 - io0.0,
            (io1.1 - io0.1) as f64 / (1u64 << 20) as f64,
            io1.2 - io0.2
        );
        assert!(!lines.is_empty());
    }
}
