//! Port of src/discovery.py (Discovery).

use std::fs::File;
use std::io::{self, Read, Write};
use std::num::NonZeroUsize;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::alignment::Record as AlignmentRecord;
use noodles_sam::Header;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::config::*;
use crate::coverage::Coverage;
use crate::evidence::{
    frag_hash, select_lowest, short_candidate, DiscLite, EvExtra, EvRec, MateReq, PaLite, ReqKind, Role, ShortCollector,
    ShortLite, Sidecar, Target, Emitted,
    FLAG_PAIRED, FLAG_REVERSE, FLAG_SUPPLEMENTARY,
};
use crate::exons::GeneModel;
use crate::filters::{both_clips_slippage, clean_clipped_seq, is_adapter, is_low_complexity, is_slippage_clip, longest_homopolymer_run, mean_kmer_diversity};
use crate::intervals::IntervalIndex;
use crate::model::{join, Breakpoint};
use crate::polya::PolyABreakpoint;
use crate::qseq::{revcomp_bytes, QualitySeq};
use crate::read::BamRead;
use crate::stats::Stats;

/// Post-TSD pairing modes, in greedy priority order (lower first).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum PairMode {
    TsdDeletion = 0,
    Blunt = 1,
    L1Del = 2,
    L1Dup = 3,
}

/// One buffered emission of `output` (see there).
struct Out<'a> {
    left: Emit<'a>,
    right: Emit<'a>,
    pal: Vec<&'a PolyABreakpoint>,
    par: Vec<&'a PolyABreakpoint>,
    dead: bool,
}

/// True when >= 80 % of the first `n` bases of a stored (junction-outward) clip are T —
/// the poly-A tail seen from the junction (both sides store it so; tolerant of the
/// interruptions/sequencing errors real and simulated poly-A tails carry).
fn leading_polyt(seq: &[u8], n: usize) -> bool {
    leading_run_of(seq, n, b'T')
}

fn leading_polya(seq: &[u8], n: usize) -> bool {
    leading_run_of(seq, n, b'A')
}

fn leading_run_of(seq: &[u8], n: usize, base: u8) -> bool {
    n > 0 && seq.len() >= n && seq[..n].iter().filter(|&&b| b.to_ascii_uppercase() == base).count() * 5 >= n * 4
}

/// SPEC-8b junction spare (`clip_slippage_junction_spare` = k > 0): the first k bases of a
/// stored clip (junction-proximal) are structured sequence — no homopolymer >= 8, >= 3
/// distinct bases. A short insert (orphan transduction tag, short element) followed by a
/// long poly-A lowers the WHOLE-clip entropy below the SPEC-8b threshold although the
/// junction itself is not a tract continuation; reference-tract slippage is homopolymer
/// right at the junction and is never spared.
fn junction_structured(seq: &[u8], k: usize) -> bool {
    if k == 0 || seq.len() < k {
        return false;
    }
    let pre = &seq[..k];
    let mut seen = [false; 4];
    for b in pre {
        match b.to_ascii_uppercase() {
            b'A' => seen[0] = true,
            b'C' => seen[1] = true,
            b'G' => seen[2] = true,
            b'T' => seen[3] = true,
            _ => {}
        }
    }
    longest_homopolymer_run(pre).0 < 8 && seen.iter().filter(|&&x| x).count() >= 3
}

/// A clip that looks like element body: long enough, no homopolymer start in either
/// orientation, not low-complexity.
fn complex_clip(seq: &[u8], n: usize) -> bool {
    seq.len() >= MIN_CLIP_LEN && !leading_polyt(seq, n) && !leading_polya(seq, n) && !is_low_complexity(seq, 0.8)
}

/// TPRT polarity of a LEFT/RIGHT pair: exactly one side is the poly-A tail, the other a
/// complex element clip (used to admit far-apart L1-mediated pairs).
fn polarised_pair(lb: &Breakpoint, rb: &Breakpoint, n: usize) -> bool {
    let (lp, rp) = (leading_polyt(&lb.clipped.seq, n), leading_polyt(&rb.clipped.seq, n));
    (lp && complex_clip(&rb.clipped.seq, n)) || (rp && complex_clip(&lb.clipped.seq, n))
}

/// Pairing mode of a LEFT/RIGHT pair outside the normal TSD window (None = not admitted).
/// gap = r - l: [-max_target_site_deletion, -1] target-site deletion; [0, tsd_min) blunt;
/// polarised and [-max_l1_mediated_span, -max_target_site_deletion-1] L1-mediated
/// deletion; polarised and (tsd_max, max_l1_mediated_span] L1-mediated duplication.
fn pair_mode(c: &DiscoveryConfig, lb: &Breakpoint, rb: &Breakpoint) -> Option<PairMode> {
    let gap = rb.breakpoint - lb.breakpoint;
    if gap >= c.tsd_min && gap <= c.tsd_max {
        return None; // normal TSD window: the legacy loop's business
    }
    if c.max_target_site_deletion > 0 && gap < 0 && gap >= -c.max_target_site_deletion {
        return Some(PairMode::TsdDeletion);
    }
    if c.allow_blunt_pairs && gap >= 0 && gap < c.tsd_min {
        return Some(PairMode::Blunt);
    }
    let span = c.max_l1_mediated_span;
    if span > 0 && polarised_pair(lb, rb, c.l1_mediated_min_polya) {
        if gap < -c.max_target_site_deletion.max(0) && gap >= -span {
            return Some(PairMode::L1Del);
        }
        if gap > c.tsd_max && gap <= span {
            return Some(PairMode::L1Dup);
        }
    }
    None
}

/// Greedy one-to-one matching of the not-yet-paired (`used != 2`) LEFT/RIGHT breakpoints
/// (both sorted by position) under `pair_mode`, ordered by mode rank, |gap|, LEFT position,
/// indices — deterministic. `reject` is the SPEC-8b both-clips gate. Returns (li, ri).
fn greedy_extra_pairs<F: Fn(&Breakpoint, &Breakpoint) -> bool>(
    c: &DiscoveryConfig,
    l: &[&Breakpoint],
    r: &[&Breakpoint],
    used_l: &[u8],
    used_r: &[u8],
    reject: F,
) -> Vec<(usize, usize)> {
    let reach = c.max_target_site_deletion.max(0).max(c.tsd_max).max(c.max_l1_mediated_span.max(0));
    let mut edges: Vec<(u8, i64, i64, usize, usize)> = Vec::new();
    for (li, lb) in l.iter().enumerate() {
        if used_l[li] == 2 {
            continue;
        }
        let lo = r.partition_point(|b| b.breakpoint < lb.breakpoint - reach);
        for (ri, rb) in r.iter().enumerate().skip(lo) {
            if rb.breakpoint > lb.breakpoint + reach {
                break;
            }
            if used_r[ri] == 2 {
                continue;
            }
            if let Some(m) = pair_mode(c, lb, rb) {
                edges.push((m as u8, (rb.breakpoint - lb.breakpoint).abs(), lb.breakpoint, li, ri));
            }
        }
    }
    edges.sort_unstable();
    let (mut ul, mut ur) = (vec![false; l.len()], vec![false; r.len()]);
    let mut out = Vec::new();
    for (_, _, _, li, ri) in edges {
        if ul[li] || ur[ri] || reject(l[li], r[ri]) {
            continue;
        }
        ul[li] = true;
        ur[ri] = true;
        out.push((li, ri));
    }
    out
}

/// One-sided locus gate: fragment floor, optional poly-A/T-tail requirement, and the
/// SPEC-8 reference-tract slippage test, applied here even when `slippage_filter` is off
/// (a lone poly-A clip at a reference A-tract is the main FP source), plus a
/// low-complexity-flank veto (satellite/tract mismap).
fn one_sided_ok(c: &DiscoveryConfig, b: &Breakpoint) -> bool {
    if b.n_frags < c.one_sided_min_fragments {
        return false;
    }
    if c.one_sided_require_polya && !leading_polyt(&b.clipped.seq, c.one_sided_min_polya) {
        return false;
    }
    if is_slippage_clip(b.side, &b.clipped.seq, &b.unclipped.seq, c.slippage_min_ref_run, c.slippage_min_clip_frac, c.slippage_max_period.max(1)) {
        return false;
    }
    !is_low_complexity(&b.unclipped.seq, 0.8)
}

#[derive(Clone, Copy)]
enum BpRef {
    Left(usize),
    Right(usize),
    PolyA(usize),
}

#[derive(Clone, Copy)]
enum Emit<'a> {
    Bp(&'a Breakpoint),
    /// One-sided locus: the missing end (no reads), named `oneside_<pos>` with the real
    /// breakpoint's coordinate.
    Open(i64),
    Pa(&'a PolyABreakpoint),
    /// Feature A: a discordant-cluster end (coordinate only, no reads). Dormant since the
    /// both-sided rescue rework — Feature A now requires a direct clip on the missing side
    /// (see `discordant_rescue`), so no call is emitted with a bare discordant end. Retained
    /// with the disc-cluster machinery for a possible future one-sided mode.
    #[allow(dead_code)]
    Disc(&'a DiscordantCluster),
}

/// Feature A: one discordant read-pair observation — a confidently-placed anchor
/// read whose mate is discordant (different contig, or far away on the same one).
/// Its innermost coordinate votes for an insertion breakpoint on the given `role`
/// side (a forward anchor bounds the insertion on its right → fills the LEFT/smaller
/// coordinate; a reverse anchor → the RIGHT/larger coordinate). `mate_ref_id`/
/// `mate_pos` record where the mate landed, for the D3 RTE-origin check.
#[derive(Clone)]
struct DiscordantObs {
    contig: String,
    role: i32, // CLIP_LEFT or CLIP_RIGHT
    pos: i64,
    mate_ref_id: Option<usize>,
    mate_pos: i64,
}

/// Feature A: a cluster of >= `discordant_min_reads` discordant observations sharing
/// a contig, role and (approximate) breakpoint. It can stand in for a missing
/// reciprocal breakpoint during output pairing. `mate_dests` are the mate landing
/// sites, used by the D3 RTE-origin check.
// Dormant since the both-sided rescue rework (fields written by cluster_discordant but no
// longer read for emission); retained for a possible future one-sided discordant mode.
#[allow(dead_code)]
#[derive(Clone)]
struct DiscordantCluster {
    contig: String,
    role: i32,
    pos: i64,
    n_reads: usize,
    /// mate landing sites; consumed by the D3 RTE-origin check.
    mate_dests: Vec<(Option<usize>, i64)>,
    /// D3: fraction of `mate_dests` landing in the RTE track (0.0 until computed).
    rte_origin: f64,
}

/// Open the BAM as a record reader. With `bam_threads > 1` the BGZF blocks are
/// decoded on a worker pool (SPD-2); otherwise the single-threaded decoder is used,
/// preserving the documented single-core footprint. Both yield identical records in
/// identical order, so output stays byte-identical.
fn open_bam(path: &str, bam_threads: usize) -> io::Result<bam::io::Reader<Box<dyn Read>>> {
    let file = File::open(path)?;
    let inner: Box<dyn Read> = if bam_threads > 1 {
        let wc = NonZeroUsize::new(bam_threads).unwrap_or(NonZeroUsize::MIN);
        Box::new(bgzf::io::MultithreadedReader::with_worker_count(wc, file))
    } else {
        Box::new(bgzf::io::Reader::new(file))
    };
    Ok(bam::io::Reader::from(inner))
}

/// True when the input is a CRAM file (by extension).
fn is_cram(path: &str) -> bool {
    path.ends_with(".cram")
}

/// True when a coordinate index (`.bai` or `.csi`) sits next to the BAM. Contig-level
/// parallelism (SPD-3/SPD-4c) needs an index to fetch each contig independently; without
/// one the single-threaded scan is used so a default `--threads > 1` stays safe on
/// un-indexed input.
fn has_bam_index(path: &str) -> bool {
    std::path::Path::new(&format!("{path}.bai")).exists()
        || std::path::Path::new(&format!("{path}.csi")).exists()
}

/// SPEC-3/4 coverage pre-pass per-record accumulation, shared by the BAM and CRAM
/// scan loops. Bins read starts per contig; flushes the previous contig's bins on a
/// contig change (input is coordinate-sorted).
fn cov_accumulate(
    cov: &mut Coverage,
    cur_id: &mut Option<usize>,
    cur_name: &mut String,
    bins: &mut Vec<u32>,
    read: &BamRead,
    bin_size: i64,
) {
    let Some(id) = read.reference_sequence_id else { return };
    if read.reference_start < 0 {
        return;
    }
    if *cur_id != Some(id) {
        if cur_id.is_some() {
            cov.set_contig(std::mem::take(cur_name), std::mem::take(bins));
        }
        *cur_id = Some(id);
        *cur_name = read.reference_name().unwrap_or_default();
        *bins = Vec::new();
    }
    let b = (read.reference_start / bin_size) as usize;
    if bins.len() <= b {
        bins.resize(b + 1, 0);
    }
    bins[b] += 1;
}

/// Open a CRAM reader backed by a reference FASTA repository (built from the indexed
/// reference). CRAM stores read bases as differences against the reference, so the
/// FASTA — with names matching the CRAM header @SQ — is required to decode sequences.
fn open_cram(path: &str, reference: Option<&str>) -> io::Result<cram::io::Reader<File>> {
    let ref_path = reference.ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidInput, "CRAM input requires a reference FASTA (--reference <ref.fa>)")
    })?;
    let fa = fasta::io::indexed_reader::Builder::default().build_from_path(ref_path)?;
    let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fa));
    cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(path)
}

/// True if an alternative-alignment CIGAR (XA or SA tag) covers essentially the
/// whole read end-to-end (clips fewer than MIN_CLIP_LEN bases). Mirrors the
/// Python `_alt_is_full_length`.
fn alt_is_full_length(cigar: &str) -> bool {
    let mut aligned: usize = 0;
    let mut clip: usize = 0;
    let mut num: usize = 0;
    let mut saw_digit = false;
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch - b'0') as usize;
            saw_digit = true;
        } else {
            if saw_digit {
                match ch {
                    b'S' | b'H' => clip += num,
                    b'M' | b'=' | b'X' => aligned += num,
                    _ => {}
                }
            }
            num = 0;
            saw_digit = false;
        }
    }
    aligned > 0 && clip < MIN_CLIP_LEN
}

/// True if the read has any XA/SA alternative alignment spanning the whole read.
/// Such reads map contiguously elsewhere and are not genuine junctions.
fn maps_fully_elsewhere(read: &BamRead<'_>) -> bool {
    // (tag string, index of the CIGAR field within a comma-separated entry)
    for (tag, cigar_idx) in [(read.xa(), 2usize), (read.sa(), 3usize)] {
        if let Some(s) = tag {
            for entry in s.split(';') {
                if entry.is_empty() {
                    continue;
                }
                let fields: Vec<&str> = entry.split(',').collect();
                if fields.len() > cigar_idx && alt_is_full_length(fields[cigar_idx]) {
                    return true;
                }
            }
        }
    }
    false
}

pub struct Discovery {
    temporary_breakpoints: Vec<Breakpoint>,
    final_left_breakpoints: Vec<Breakpoint>,
    final_right_breakpoints: Vec<Breakpoint>,
    /// Feature A both-sided rescue: clip clusters that passed all quality gates but fell
    /// below the normal evidence floor (>= discordant_partner_min_reads, < evidence_floor).
    /// Not part of normal output; only `discordant_rescue` pairs one of these (the missing
    /// junction) with a solid anchor breakpoint. Empty unless `discordant_anchor` is on.
    subfloor_left_breakpoints: Vec<Breakpoint>,
    subfloor_right_breakpoints: Vec<Breakpoint>,
    polya: Vec<PolyABreakpoint>,
    reference_name: Option<String>,
    filepath: String,
    config: DiscoveryConfig,
    exclude: Option<IntervalIndex>,
    rm_mask: Option<IntervalIndex>,
    coverage: Coverage,
    /// SPD-3: when set, extract_chimeric processes only this contig (indexed fetch).
    /// Used by the parallel per-contig workers; None = the full single-threaded scan.
    only_contig: Option<usize>,
    /// Feature A: discordant read-pair observations collected during the scan (only
    /// when `discordant_anchor` is on). Persist across contigs like `polya`.
    discordant_obs: Vec<DiscordantObs>,
    /// Feature A: clusters built from `discordant_obs` in a prepass, consumed by the
    /// output rescue. Empty unless `discordant_anchor` is on.
    discordant_clusters: Vec<DiscordantCluster>,
    /// D3: RTE track for the mate-origin check (from `discordant_rte_track`). Set on
    /// the parent only (clustering runs after the SPD-3 merge); workers leave it None.
    discordant_rte: Option<IntervalIndex>,
    /// D5: exon model for the splice / processed-pseudogene annotation (from
    /// `exon_annotation`). Parent only; consumed by `splice_annotate`.
    exon_model: Option<GeneModel>,
    /// Reference FASTA path for CRAM decoding (`--reference`); None for BAM input.
    reference_path: Option<String>,
    stats: Stats,
    bam_threads: usize, // BGZF decode workers (SPD-2); 1 = single-threaded
    /// TPRT sidecar: the current contig's discordant-anchor observations (compact, no
    /// sequence), attached to nearby final breakpoints in `cleanup`. Empty unless
    /// `evidence_sidecar` is on.
    sc_disc_tmp: Vec<DiscLite>,
    /// TPRT sidecar: short-overhang candidate collector for the current contig (None
    /// unless `evidence_sidecar` && `short_overhang_evidence`).
    sc_short: Option<ShortCollector>,
}

impl Discovery {
    pub fn new(
        filepath: String,
        bam_threads: usize,
        config: DiscoveryConfig,
        exclude: Option<IntervalIndex>,
        rm_mask: Option<IntervalIndex>,
    ) -> Self {
        let coverage = Coverage::new(config.coverage_bin_size);
        Discovery {
            temporary_breakpoints: Vec::new(),
            final_left_breakpoints: Vec::new(),
            final_right_breakpoints: Vec::new(),
            subfloor_left_breakpoints: Vec::new(),
            subfloor_right_breakpoints: Vec::new(),
            polya: Vec::new(),
            reference_name: None,
            filepath,
            config,
            exclude,
            rm_mask,
            coverage,
            only_contig: None,
            discordant_obs: Vec::new(),
            discordant_clusters: Vec::new(),
            discordant_rte: None,
            exon_model: None,
            reference_path: None,
            stats: Stats::default(),
            bam_threads: bam_threads.max(1),
            sc_disc_tmp: Vec::new(),
            sc_short: None,
        }
    }

    /// True when either coverage-based gate (SPEC-3/SPEC-4) is enabled, or the
    /// discordant one-sided-call coverage ceiling needs the median populated.
    fn coverage_enabled(&self) -> bool {
        self.config.coverage_mask
            || self.config.adaptive_evidence
            || (self.config.discordant_anchor && self.config.discordant_coverage_max_mult.is_some())
    }

    /// Stricter local-coverage ceiling for a one-sided discordant call. True (pass) when
    /// no ceiling is set, coverage was not estimated, or the real breakpoint's local depth
    /// is within `discordant_coverage_max_mult` x the genome median.
    fn disc_coverage_ok(&self, rn: &str, pos: i64) -> bool {
        let Some(mult) = self.config.discordant_coverage_max_mult else { return true };
        let med = self.coverage.median();
        if med <= 0.0 {
            return true;
        }
        (self.coverage.local(rn, pos) as f64) <= mult * med
    }

    /// Reject a one-sided discordant call whose real breakpoint's mate reads are a
    /// low-diversity satellite array. True (pass) when no threshold is set or the mates
    /// are too few/short to judge (`None` diversity) — the gate only fires on positive
    /// evidence of low complexity. Dormant since the both-sided rescue rework (no one-sided
    /// disc-cluster calls to gate); retained for a possible future one-sided mode.
    #[allow(dead_code)]
    fn disc_mates_ok(&self, bp: &Breakpoint) -> bool {
        let Some(thr) = self.config.discordant_mate_min_kmer_div else { return true };
        match mean_kmer_diversity(&bp.mate_seqs, 4) {
            Some(div) => div >= thr,
            None => true,
        }
    }

    /// SPEC-3/4 pre-pass: bin read starts per contig and compute the median depth,
    /// so the global median is known before the per-contig cleanups run. Contig name
    /// is resolved once per contig (coordinate-sorted), so no per-read allocation.
    fn estimate_coverage(&mut self) -> io::Result<()> {
        let bin_size = self.config.coverage_bin_size.max(1);
        let mut cur_id: Option<usize> = None;
        let mut cur_name = String::new();
        let mut bins: Vec<u32> = Vec::new();
        if is_cram(&self.filepath) {
            let ref_path = self.reference_path.clone();
            let mut reader = open_cram(&self.filepath, ref_path.as_deref())?;
            let header = reader.read_header()?;
            for result in reader.records(&header) {
                let rec = result?;
                let read = BamRead::from_record(&rec, &header)?;
                cov_accumulate(&mut self.coverage, &mut cur_id, &mut cur_name, &mut bins, &read, bin_size);
            }
        } else {
            let mut reader = open_bam(&self.filepath, self.bam_threads)?;
            let header = reader.read_header()?;
            let mut record = bam::Record::default();
            while reader.read_record(&mut record)? != 0 {
                let read = BamRead::from_record(&record, &header)?;
                cov_accumulate(&mut self.coverage, &mut cur_id, &mut cur_name, &mut bins, &read, bin_size);
            }
        }
        if cur_id.is_some() {
            self.coverage.set_contig(cur_name, bins);
        }
        self.coverage.finalize(self.config.coverage_sample_size);
        // Region-slice runs pin the genome median from the full BAM: local bins (near the
        // truth loci) stay exact, but the slice-estimated median would be inflated.
        if let Some(m) = self.config.coverage_median_override {
            self.coverage.set_median(m);
        }
        Ok(())
    }

    /// SPEC-4: evidence floor for a breakpoint at `pos`, scaled by local/median
    /// coverage but never below the base `min_evidence_reads_per_breakpoint`.
    fn evidence_floor(&self, pos: i64) -> usize {
        // TPRT fragment mode replaces the read floor with a fragment floor (same scaling).
        let base = self
            .config
            .min_evidence_fragments_per_sample
            .unwrap_or(self.config.min_evidence_reads_per_breakpoint);
        if !self.config.adaptive_evidence {
            return base;
        }
        let med = self.coverage.median();
        if med <= 0.0 {
            return base;
        }
        let name = self.reference_name.as_deref().unwrap_or_default();
        let local = self.coverage.local(name, pos) as f64;
        let scaled = (base as f64 * (local / med)).round() as i64;
        scaled.max(base as i64) as usize
    }

    /// OBS-1 reject-counter sidecar, serialised as JSON.
    pub fn stats_json(&self) -> String {
        self.stats.to_json()
    }

    /// D3: install the RTE track used for the discordant mate-origin check.
    pub fn set_discordant_rte(&mut self, ix: Option<IntervalIndex>) {
        self.discordant_rte = ix;
    }

    /// D5: install the exon model used for the splice / pseudogene annotation.
    pub fn set_exon_model(&mut self, m: Option<GeneModel>) {
        self.exon_model = m;
    }

    /// Install the reference FASTA path used to decode CRAM input.
    pub fn set_reference_path(&mut self, p: Option<String>) {
        self.reference_path = p;
    }

    /// Run only the coverage pre-pass and return the estimated genome median. Used by the
    /// `coverage-median` step to precompute the true full-BAM median that region-slice runs
    /// then pin via `PEARTREE_COVERAGE_MEDIAN`. Ignores `coverage_median_override` so it
    /// always reports the freshly estimated value.
    pub fn compute_coverage_median(&mut self) -> io::Result<f64> {
        let saved = self.config.coverage_median_override.take();
        self.estimate_coverage()?;
        self.config.coverage_median_override = saved;
        Ok(self.coverage.median())
    }

    /// True if a breakpoint at (rn, pos) would survive the output-time masks
    /// (contig-OK + SPEC-5/3/7), i.e. it is actually emitted.
    fn bp_visible(&self, rn: &str, pos: i64) -> bool {
        if !self.contig_ok_output(rn) || self.excluded(rn, pos) {
            return false;
        }
        if self.config.coverage_mask && self.coverage.median() > 0.0 {
            let thr = self.config.coverage_mask_multiplier * self.coverage.median();
            if self.coverage.local(rn, pos) as f64 > thr {
                return false;
            }
        }
        if let Some(rm) = &self.rm_mask {
            if rm.contains(rn, pos) {
                return false;
            }
        }
        true
    }

    /// D5 (Feature B): scan output-visible breakpoints for the processed-pseudogene
    /// signature — mates hitting >= `splice_min_exons` distinct exons of one gene while
    /// skipping the introns between them — and return the `<out>.splice.tsv` body
    /// (empty unless `splice_hallmark` is on with an exon model). Non-gating: the main
    /// output is untouched.
    pub fn splice_annotate(&self) -> io::Result<Vec<u8>> {
        let mut buf: Vec<u8> = Vec::new();
        if !self.config.splice_hallmark {
            return Ok(buf);
        }
        let Some(model) = &self.exon_model else { return Ok(buf) };
        let names = self.reference_names().unwrap_or_default();
        writeln!(buf, "contig\tbreakpoint\tside\tgene\tn_exons\tintron_bp\tspan_bp")?;
        for (side, bps) in [("LEFT", &self.final_left_breakpoints), ("RIGHT", &self.final_right_breakpoints)] {
            for bp in bps.iter() {
                if !self.bp_visible(&bp.reference_name, bp.breakpoint) {
                    continue;
                }
                self.splice_row(&mut buf, model, &names, side, bp)?;
            }
        }
        Ok(buf)
    }

    /// Emit a splice row for one breakpoint if its mate destinations span multiple
    /// exons of a single gene with an intron skipped.
    fn splice_row(&self, buf: &mut Vec<u8>, model: &GeneModel, names: &[String], side: &str, bp: &Breakpoint) -> io::Result<()> {
        // per gene: distinct exon ranks hit, and the min-begin / max-end / summed
        // exon length over those exons.
        struct GeneHit {
            ranks: FxHashSet<usize>,
            min_begin: i64,
            max_end: i64,
            exon_len: FxHashMap<usize, i64>,
        }
        let mut genes: FxHashMap<String, GeneHit> = FxHashMap::default();
        for &(mref, mpos) in &bp.mate_dests {
            let Some(id) = mref else { continue };
            let Some(name) = names.get(id) else { continue };
            if let Some(exon) = model.lookup(name, mpos) {
                let g = genes.entry(exon.gene.clone()).or_insert_with(|| GeneHit {
                    ranks: FxHashSet::default(),
                    min_begin: i64::MAX,
                    max_end: i64::MIN,
                    exon_len: FxHashMap::default(),
                });
                g.ranks.insert(exon.rank);
                g.min_begin = g.min_begin.min(exon.begin);
                g.max_end = g.max_end.max(exon.end);
                g.exon_len.insert(exon.rank, exon.end - exon.begin);
            }
        }
        for (gene, h) in genes {
            if h.ranks.len() < self.config.splice_min_exons {
                continue;
            }
            let span = h.max_end - h.min_begin;
            let summed: i64 = h.exon_len.values().sum();
            // intron skipped: the genomic span exceeds the summed exon lengths.
            if span > summed {
                writeln!(
                    buf,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    bp.reference_name,
                    bp.breakpoint,
                    side,
                    gene,
                    h.ranks.len(),
                    span - summed,
                    span
                )?;
            }
        }
        Ok(())
    }

    /// Contig filter for the read-scan (extract_chimeric). With an allowlist set the
    /// name must be in it; otherwise the legacy `len(name) <= 5` + not-MT heuristic.
    fn contig_ok_extract(&self, name: &str) -> bool {
        match &self.config.contig_allowlist {
            Some(set) => set.contains(name),
            None => name.len() <= 5 && name != "MT" && name != "chrM",
        }
    }

    /// Contig filter at output. Legacy path checks only `len(name) <= 5` (MT/chrM
    /// clipped breakpoints were already dropped in the scan, but MT polyA mates were
    /// historically emitted), so keep that exact behaviour when no allowlist is set.
    fn contig_ok_output(&self, name: &str) -> bool {
        match &self.config.contig_allowlist {
            Some(set) => set.contains(name),
            None => name.len() <= 5,
        }
    }

    /// True if a breakpoint at (name, pos) falls in an exclude-BED region (SPEC-5).
    fn excluded(&self, name: &str, pos: i64) -> bool {
        self.exclude.as_ref().is_some_and(|ix| ix.contains(name, pos))
    }

    #[allow(clippy::too_many_arguments)]
    fn add_breakpoint(
        &mut self,
        side: i32,
        breakpoint: i64,
        qname: &str,
        clipped: QualitySeq,
        unclipped: QualitySeq,
        is_read1: bool,
        is_forward: bool,
        exclude: bool,
        mapq: u8,
        flag: u16,
    ) {
        let mut bp = Breakpoint::new(
            side,
            self.reference_name.clone().unwrap_or_default(),
            breakpoint,
            Some(qname.to_string()),
            clipped.upper(),
            unclipped.upper(),
            Some(is_read1),
            Some(is_forward),
            exclude,
            mapq,
        );
        bp.flag = flag;
        self.temporary_breakpoints.push(bp);
    }

    /// Debug-gated (`PEARTREE_MEM_DEBUG`) phase-boundary line attributing resident
    /// bytes to each accumulator. `tmp` is the peak `temporary_breakpoints` slice
    /// for this tag (empty when it has already been drained).
    fn mem_report(&self, tag: &str, tmp: &[Breakpoint]) {
        let disc_bytes: usize = self
            .discordant_obs
            .iter()
            .map(|o| std::mem::size_of::<DiscordantObs>() + o.contig.capacity())
            .sum();
        crate::mem::report(
            tag,
            tmp,
            &self.final_left_breakpoints,
            &self.final_right_breakpoints,
            &self.polya,
            self.discordant_obs.len(),
            disc_bytes,
        );
    }

    fn cleanup(&mut self) {
        // short-overhang candidates of the just-finished contig (empty when off)
        let shorts: Vec<ShortLite> = self.sc_short.as_mut().map(|c| c.finish()).unwrap_or_default();
        if self.reference_name.is_none() {
            self.sc_disc_tmp.clear();
            return;
        }
        if self.temporary_breakpoints.is_empty() {
            self.sc_disc_tmp.clear();
            return;
        }
        let (l0, r0) = (self.final_left_breakpoints.len(), self.final_right_breakpoints.len());
        // temporary_breakpoints holds the just-finished contig's full clipped-read
        // set here, right before it is drained — its true per-contig peak.
        if crate::mem::enabled() {
            let tag = format!("cleanup {}", self.reference_name.as_deref().unwrap_or("?"));
            self.mem_report(&tag, &self.temporary_breakpoints);
        }
        let mut left_bps: Vec<Breakpoint> = Vec::new();
        let mut right_bps: Vec<Breakpoint> = Vec::new();
        for bp in self.temporary_breakpoints.drain(..) {
            if bp.side == CLIP_LEFT {
                left_bps.push(bp);
            } else {
                right_bps.push(bp);
            }
        }
        // stable sort by breakpoint position
        left_bps.sort_by_key(|b| b.breakpoint);
        right_bps.sort_by_key(|b| b.breakpoint);

        for (side_bps, out) in [
            (left_bps, true),  // true -> final_left
            (right_bps, false),
        ] {
            let mut groups: Vec<Vec<Breakpoint>> = Vec::new();
            let mut current: Vec<Breakpoint> = Vec::new();
            for bp in side_bps {
                if let Some(last) = current.last() {
                    if (bp.breakpoint - last.breakpoint).abs() < self.config.cluster_window {
                        current.push(bp);
                        continue;
                    } else {
                        groups.push(std::mem::take(&mut current));
                    }
                }
                current = vec![bp];
            }
            if !current.is_empty() {
                groups.push(current);
            }
            for g in groups {
                let floor = self.evidence_floor(g[0].breakpoint);
                // Feature A: also build the consensus at the lower partner floor, so a
                // quality-passing clip cluster below the normal evidence floor is kept as a
                // candidate missing junction for discordant_rescue. join() only differs by
                // floor in its low-support branch, so at n_reads >= floor this is identical
                // to the final consensus (which we take instead, below). Uses throwaway stats
                // to avoid double-counting the OBS-1 reject counters.
                let sub = if self.config.discordant_anchor {
                    let mut sink = Stats::default();
                    join(g.clone(), &self.config, self.config.discordant_partner_min_reads, &mut sink)
                } else {
                    None
                };
                match join(g, &self.config, floor, &mut self.stats) {
                    Some(joined) => {
                        if out {
                            self.final_left_breakpoints.push(joined);
                        } else {
                            self.final_right_breakpoints.push(joined);
                        }
                    }
                    // below the normal floor but a valid clip cluster >= partner floor:
                    // keep only as a discordant missing-junction candidate.
                    None => {
                        if let Some(s) = sub {
                            if out {
                                self.subfloor_left_breakpoints.push(s);
                            } else {
                                self.subfloor_right_breakpoints.push(s);
                            }
                        }
                    }
                }
            }
        }
        // TPRT sidecar: attach this contig's discordant anchors to its new breakpoints.
        let disc = std::mem::take(&mut self.sc_disc_tmp);
        if self.config.evidence_sidecar && !disc.is_empty() {
            let span = self.config.sidecar_disc_span;
            let cap = self.config.max_evidence_reads_per_breakpoint;
            attach_disc(&mut self.final_left_breakpoints[l0..], &disc, true, span, cap);
            attach_disc(&mut self.final_right_breakpoints[r0..], &disc, false, span, cap);
        }
        if !shorts.is_empty() {
            let c = &self.config;
            let (w, m, cap) = (c.short_overhang_window, c.short_overhang_max, c.max_short_per_breakpoint);
            attach_short(&mut self.final_left_breakpoints[l0..], &shorts, true, w, m, cap);
            attach_short(&mut self.final_right_breakpoints[r0..], &shorts, false, w, m, cap);
        }
    }

    pub fn extract_chimeric(&mut self) -> io::Result<()> {
        let reject_fullmap = self.config.reject_fully_mapping_reads;
        let min_mapq = self.config.min_mapq;
        self.reference_name = None;
        let mut current_ref_id: Option<usize> = None;
        if let Some(cid) = self.only_contig {
            // SPD-3 worker path: fetch just this contig via the index.
            let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(&self.filepath)?;
            let header = reader.read_header()?;
            let name: Vec<u8> = header
                .reference_sequences()
                .get_index(cid)
                .map(|(k, _)| {
                    let b: &[u8] = k.as_ref();
                    b.to_vec()
                })
                .unwrap_or_default();
            let region = Region::new(name, ..);
            for rec in reader.query(&header, &region)? {
                let rec = rec?;
                self.handle_record(&rec, &header, min_mapq, reject_fullmap, &mut current_ref_id)?;
            }
        } else if is_cram(&self.filepath) {
            // native CRAM full scan (reference required to decode read sequences).
            let ref_path = self.reference_path.clone();
            let mut reader = open_cram(&self.filepath, ref_path.as_deref())?;
            let header = reader.read_header()?;
            for result in reader.records(&header) {
                let rec = result?;
                self.handle_record(&rec, &header, min_mapq, reject_fullmap, &mut current_ref_id)?;
            }
        } else {
            let mut reader = open_bam(&self.filepath, self.bam_threads)?;
            let header = reader.read_header()?;
            let mut record = bam::Record::default();
            while reader.read_record(&mut record)? != 0 {
                self.handle_record(&record, &header, min_mapq, reject_fullmap, &mut current_ref_id)?;
            }
        }
        self.cleanup();
        Ok(())
    }

    /// Per-read body shared by the full-scan and the SPD-3 per-contig indexed paths.
    /// A read-level `continue` in the original single loop is a `return Ok(())` here.
    fn handle_record(
        &mut self,
        record: &dyn AlignmentRecord,
        header: &Header,
        min_mapq: u8,
        reject_fullmap: bool,
        current_ref_id: &mut Option<usize>,
    ) -> io::Result<()> {
        let read = BamRead::from_record(record, header)?;

        if read.mapq < min_mapq {
            // Mate-anchored rescue: keep a soft-clipped read below the MAPQ floor if its
            // mate maps uniquely (MQ >= min_mapq). Mirrors Python src/discovery.py for the
            // proper-pair case. Feature A extends it to DISCORDANT (non-proper) pairs when
            // `discordant_anchor` is on: a read whose body maps poorly INTO a repetitive RTE
            // (low MAPQ) but whose clip pins the junction is exactly the missing-side
            // evidence we want — a few genomic bases in the clip localise the breakpoint, and
            // the uniquely-mapped mate is the required high-MAPQ read of the pair. Otherwise
            // the read only survives as a polyA mate.
            let has_clip = (read.left_is_soft && read.left_len >= MIN_CLIP_LEN)
                || (read.right_is_soft && read.right_len >= MIN_CLIP_LEN);
            let mate_unique = read.mate_is_mapped && read.mq().map_or(false, |mq| mq >= min_mapq);
            let rescued = self.config.mate_anchor_rescue
                && has_clip
                && mate_unique
                && (read.is_proper_pair || self.config.discordant_anchor);
            if !rescued {
                // NB: this path never checked the 0x400 flag (legacy, kept for
                // byte-identity); with `ignore_dup_flag` dups are kept everywhere anyway.
                // `drop_dup_in_polya_path` applies the clip-path drop rule here too.
                if self.config.drop_dup_in_polya_path && drop_read(&read, self.config.ignore_dup_flag) {
                    return Ok(());
                }
                if let Some(mut b) = PolyABreakpoint::find_polya(&read) {
                    if self.config.evidence_sidecar {
                        let mate_ok = read.mate_is_mapped && read.mate_pos >= 0;
                        b.sc = Some(PaLite {
                            flag: read.flag,
                            ref_id: if read.mapped { read.reference_sequence_id.map_or(-1, |i| i as i32) } else { -1 },
                            pos: if read.mapped { read.reference_start } else { -1 },
                            mref: if mate_ok { read.mate_ref_id.map_or(-1, |i| i as i32) } else { -1 },
                            mpos: if mate_ok { read.mate_pos } else { -1 },
                        });
                    }
                    self.polya.push(b);
                }
                return Ok(());
            }
            // else: fall through and process this clip as a breakpoint
        }
        if drop_read(&read, self.config.ignore_dup_flag) {
            return Ok(());
        }
        // SPD-5a: detect contig change by integer reference id; resolve the name
        // string only when it changes, not for every read.
        let ref_id = match read.reference_sequence_id {
            Some(id) => id,
            None => return Ok(()), // high-mapq unmapped: skip (would crash pysam)
        };
        if *current_ref_id != Some(ref_id) {
            self.cleanup();
            *current_ref_id = Some(ref_id);
            self.reference_name = read.reference_name();
        }
        let ref_name = self.reference_name.as_deref().unwrap_or_default();
        if !self.contig_ok_extract(ref_name) {
            return Ok(());
        }

        // Feature A: collect a discordant read-pair observation from this primary,
        // high-MAPQ, contig-OK anchor. Discordant = paired-but-not-proper with a
        // mapped mate on a different contig, or far away on the same one (mate strand
        // is not in the record, so the non-FR criterion is intentionally omitted).
        // Independent of the clip logic below, so a read may also become a breakpoint.
        if self.config.discordant_anchor
            && !read.is_proper_pair
            && read.mate_is_mapped
            && !read.is_supplementary
        {
            if let Some(mref) = read.mate_ref_id {
                let same_contig = read.reference_sequence_id == read.mate_ref_id;
                let discordant = !same_contig
                    || (read.mate_pos - read.reference_start).abs() > self.config.discordant_max_tlen;
                // Track-free RTE-origin proxy: keep only mates that map ambiguously
                // (MQ <= threshold), i.e. into a repeat family, when the gate is set.
                let mate_ambiguous = match self.config.discordant_mate_max_mapq {
                    Some(maxq) => read.mq().map_or(false, |q| q <= maxq),
                    None => true,
                };
                if discordant && mate_ambiguous && read.mate_pos >= 0 {
                    let (role, pos) = if read.is_reverse {
                        (CLIP_RIGHT, read.reference_start)
                    } else {
                        (CLIP_LEFT, read.reference_end)
                    };
                    self.discordant_obs.push(DiscordantObs {
                        contig: ref_name.to_string(),
                        role,
                        pos,
                        mate_ref_id: Some(mref),
                        mate_pos: read.mate_pos,
                    });
                    self.stats.disc_obs += 1;
                }
            }
        }

        // TPRT sidecar: compact discordant-anchor observation (mate unmapped, on another
        // contig, far away, or same strand) from a primary high-MAPQ read. Attached to a
        // nearby final breakpoint in cleanup(); sequence fetched in the mate pass.
        if self.config.evidence_sidecar
            && read.flag & FLAG_PAIRED != 0
            && !read.is_proper_pair
            && !read.is_supplementary
            && read.mapped
            && read.mapq >= min_mapq
            && read.reference_start >= 0
        {
            let mate_ok = read.mate_is_mapped && read.mate_pos >= 0 && read.mate_ref_id.is_some();
            let discordant = !mate_ok
                || read.mate_ref_id != read.reference_sequence_id
                || (read.mate_pos - read.reference_start).abs() > self.config.discordant_max_tlen
                || read.mate_is_reverse == read.is_reverse;
            if discordant {
                self.sc_disc_tmp.push(DiscLite {
                    frag: frag_hash(read.name_bytes()),
                    flag: read.flag,
                    ref_id: ref_id as i32,
                    start: read.reference_start,
                    end: read.reference_end,
                });
            }
        }

        if !read.has_cigar {
            return Ok(());
        }

        // TPRT sidecar: short-overhang candidates (primary, MAPQ >= min_mapq; the 0x400
        // rule was applied by drop_read above).
        if self.config.evidence_sidecar && self.config.short_overhang_evidence {
            if read.mapped && !read.is_supplementary && read.mapq >= min_mapq && read.reference_start >= 0 {
                let (w, m) = (self.config.short_overhang_window, self.config.short_overhang_max);
                self.sc_short.get_or_insert_with(|| ShortCollector::new(w, m)).push(ShortLite {
                    frag: frag_hash(read.name_bytes()),
                    flag: read.flag,
                    ref_id: ref_id as i32,
                    start: read.reference_start,
                    end: read.reference_end,
                    lead_soft: read.lead_soft as u32,
                    trail_soft: read.trail_soft as u32,
                });
            }
        }

        let left_soft = read.left_is_soft;
        let right_soft = read.right_is_soft;
        let left_len = read.left_len;
        let right_len = read.right_len;

        let clip: Option<i32> = if left_soft && !right_soft {
            Some(CLIP_LEFT)
        } else if right_soft && !left_soft {
            Some(CLIP_RIGHT)
        } else if right_soft && left_soft {
            if left_len > right_len {
                Some(CLIP_LEFT)
            } else if right_len > left_len {
                Some(CLIP_RIGHT)
            } else {
                None
            }
        } else {
            None
        };

        let Some(clip) = clip else { return Ok(()) };

        // drop reads that map contiguously elsewhere in the reference
        // (XA/SA full-length alt): not real chimeric junctions.
        if reject_fullmap && maps_fully_elsewhere(&read) {
            return Ok(());
        }

        // cruciform / short-indel exclusion via SA tag
        let mut exclude_flag = false;
        if read.is_supplementary {
            if let Some(sa) = read.sa() {
                for part in sa.split(';') {
                    if part.is_empty() {
                        continue; // inner loop over SA parts
                    }
                    let fields: Vec<&str> = part.split(',').collect();
                    if fields.len() >= 2 && fields[0] == ref_name {
                        if let Ok(start) = fields[1].parse::<i64>() {
                            if (read.reference_start - start).abs() < self.config.exclude_same_contig_supplementary {
                                exclude_flag = true;
                            }
                        }
                    }
                }
            }
        }

        // surviving clip candidate: now decode the heavy fields (SPD-1).
        let seq = read.seq();
        let qname = read.query_name();
        let full = QualitySeq::new(seq.clone(), read.qual());
        let n = seq.len();

        if clip == CLIP_LEFT {
            // adapter check on the last MIN_CLIP_LEN of the clipped part, revcomped
            let clip_part = &seq[..left_len];
            let sub = &clip_part[clip_part.len().saturating_sub(MIN_CLIP_LEN)..];
            if is_adapter(&revcomp_bytes(sub)) {
                return Ok(());
            }
            let clipped = full.pyslice(None, Some(left_len as isize));
            let unclipped = if right_soft {
                full.pyslice(Some(left_len as isize), Some(-(right_len as isize)))
            } else {
                full.pyslice(Some(left_len as isize), None)
            };
            self.add_breakpoint(clip, read.reference_start, &qname, clipped, unclipped, read.is_read1, read.is_forward(), exclude_flag, read.mapq, read.flag);
            if let Some(c) = self.sc_short.as_mut() {
                c.hot(read.reference_start, true);
            }
        } else {
            // CLIP_RIGHT: adapter check on the first MIN_CLIP_LEN of the clipped part
            let clip_part = &seq[n - right_len..];
            let sub = &clip_part[..clip_part.len().min(MIN_CLIP_LEN)];
            if is_adapter(sub) {
                return Ok(());
            }
            let clipped = full.pyslice(Some(-(right_len as isize)), None);
            let unclipped = if left_soft {
                full.pyslice(Some(left_len as isize), Some(-(right_len as isize)))
            } else {
                full.pyslice(None, Some(-(right_len as isize)))
            };
            self.add_breakpoint(clip, read.reference_end, &qname, clipped, unclipped, read.is_read1, read.is_forward(), exclude_flag, read.mapq, read.flag);
            if let Some(c) = self.sc_short.as_mut() {
                c.hot(read.reference_end, false);
            }
        }
        Ok(())
    }

    fn get_mates(&self) -> (FxHashSet<String>, FxHashSet<String>, FxHashMap<String, BpRef>) {
        let mut read1_bp: FxHashSet<String> = FxHashSet::default();
        let mut read2_bp: FxHashSet<String> = FxHashSet::default();
        let mut qmap: FxHashMap<String, BpRef> = FxHashMap::default();

        for (idx, bp) in self.final_left_breakpoints.iter().enumerate() {
            for (read_1, qname) in &bp.mates {
                if *read_1 {
                    read1_bp.insert(qname.clone());
                } else {
                    read2_bp.insert(qname.clone());
                }
                qmap.insert(qname.clone(), BpRef::Left(idx));
            }
        }
        for (idx, bp) in self.final_right_breakpoints.iter().enumerate() {
            for (read_1, qname) in &bp.mates {
                if *read_1 {
                    read1_bp.insert(qname.clone());
                } else {
                    read2_bp.insert(qname.clone());
                }
                qmap.insert(qname.clone(), BpRef::Right(idx));
            }
        }
        for (idx, pa) in self.polya.iter().enumerate() {
            if qmap.contains_key(&pa.qname) {
                continue; // real breakpoints have priority
            }
            qmap.insert(pa.qname.clone(), BpRef::PolyA(idx));
            if pa.is_read1 {
                read2_bp.insert(pa.qname.clone());
            } else {
                read1_bp.insert(pa.qname.clone());
            }
        }
        let union: FxHashSet<String> = read1_bp.intersection(&read2_bp).cloned().collect();
        for q in &union {
            read1_bp.remove(q);
            read2_bp.remove(q);
        }
        (read1_bp, read2_bp, qmap)
    }

    pub fn find_mates(&mut self) -> io::Result<()> {
        // The second (mate) pass is a single linear decode of the file. For BAM the
        // decode is already parallelised by `bam_threads` (N-way BGZF), which beats a
        // contig-partitioned scan (big chromosomes dominate and unbalance the workers);
        // for CRAM the scan is single-threaded. Two alternatives were tried and dropped:
        // per-mate indexed point-fetch (one query per mate, ~1e4-1e5 genome-wide — far
        // slower than one scan) and a contig-parallel scan (slower than BGZF-serial on
        // BAM, and the CRAM index query returned a different record set than the linear
        // scan — not byte-identical). The linear scan stays the validated path.
        self.find_mates_scan()
    }

    /// Validated path: a second linear pass matching mates by qname.
    fn find_mates_scan(&mut self) -> io::Result<()> {
        let (read1_mates, read2_mates, qmap) = self.get_mates();
        // The qname hashmaps live for the whole mate scan and are dropped on return,
        // so this is the only point that sees them — report their footprint here.
        if crate::mem::enabled() {
            let map_bytes: usize = qmap.keys().map(|k| k.capacity()).sum::<usize>()
                + qmap.len() * std::mem::size_of::<(String, BpRef)>()
                + read1_mates.iter().map(|k| k.capacity()).sum::<usize>()
                + read2_mates.iter().map(|k| k.capacity()).sum::<usize>()
                + (read1_mates.len() + read2_mates.len()) * std::mem::size_of::<String>();
            let rss = crate::mem::rss_bytes()
                .map(|b| format!("rss={:.0}MiB ", b as f64 / 1048576.0))
                .unwrap_or_default();
            eprintln!(
                "[mem] find_mates maps: {rss}qmap={} read1={} read2={} (~{:.1}MiB of qname keys)",
                qmap.len(),
                read1_mates.len(),
                read2_mates.len(),
                map_bytes as f64 / 1048576.0,
            );
        }
        // Feature B: only pay for mate-destination capture when a consumer is enabled.
        let capture_dests = self.config.splice_hallmark || self.config.discordant_anchor;
        let min_mapq = self.config.min_mapq;
        if is_cram(&self.filepath) {
            let ref_path = self.reference_path.clone();
            let mut reader = open_cram(&self.filepath, ref_path.as_deref())?;
            let header = reader.read_header()?;
            for result in reader.records(&header) {
                let rec = result?;
                let read = BamRead::from_record(&rec, &header)?;
                self.apply_mate_record(&read, &read1_mates, &read2_mates, &qmap, capture_dests, min_mapq);
            }
        } else {
            let mut reader = open_bam(&self.filepath, self.bam_threads)?;
            let header = reader.read_header()?;
            let mut record = bam::Record::default();
            while reader.read_record(&mut record)? != 0 {
                let read = BamRead::from_record(&record, &header)?;
                self.apply_mate_record(&read, &read1_mates, &read2_mates, &qmap, capture_dests, min_mapq);
            }
        }
        Ok(())
    }

    /// TPRT sidecar, before the mate pass: (1) drop the payload of every final
    /// breakpoint that cannot be emitted (not output-visible, and no opposite-side
    /// breakpoint in the TSD window nor a poly-A read whose mate lands within poly-A
    /// pairing range) — a conservative superset of `output`'s pairing, so memory is spent
    /// only on plausible loci; (2) build the sorted capture requests: mates (the
    /// `has_mate` orientation, or every CLIP/DISC read with `fetch_all_mates`), capped at
    /// `max_mates_per_breakpoint` by lowest fragment hash, plus each DISC anchor itself.
    fn sidecar_prepare(&mut self, emitted: Option<&Emitted>) -> io::Result<Vec<MateReq>> {
        let names = self.reference_names()?;
        let name_id: FxHashMap<&str, i32> = names.iter().enumerate().map(|(i, n)| (n.as_str(), i as i32)).collect();
        let slack: i64 = 1000; // mate start vs end + read length, generous
        // predicted poly-A end positions (from the anchoring mate), per contig
        let mut pa_pos: FxHashMap<&str, Vec<i64>> = FxHashMap::default();
        for pa in &self.polya {
            if let Some(l) = &pa.sc {
                if l.mref >= 0 {
                    if let Some(n) = names.get(l.mref as usize) {
                        pa_pos.entry(n.as_str()).or_default().push(l.mpos);
                    }
                }
            }
        }
        for v in pa_pos.values_mut() {
            v.sort_unstable();
        }
        let positions = |bps: &[Breakpoint], keep: Option<&[bool]>| {
            let mut m: FxHashMap<String, Vec<i64>> = FxHashMap::default();
            for (i, b) in bps.iter().enumerate() {
                if keep.map_or(true, |k| k[i]) {
                    m.entry(b.reference_name.clone()).or_default().push(b.breakpoint);
                }
            }
            for v in m.values_mut() {
                v.sort_unstable();
            }
            m
        };
        let lpos = positions(&self.final_left_breakpoints, None);
        let rpos = positions(&self.final_right_breakpoints, None);
        let any_in = |v: Option<&Vec<i64>>, lo: i64, hi: i64| -> bool {
            v.is_some_and(|v| {
                let i = v.partition_point(|&x| x < lo);
                i < v.len() && v[i] <= hi
            })
        };
        let (mut tmin, tmax, far) = (self.config.tsd_min, self.config.tsd_max, self.config.polya_far_dist);
        // post-TSD pairing modes widen the partner window (a conservative superset of
        // `extra_pairs`); one-sided loci need no partner at all.
        let c = &self.config;
        if c.max_target_site_deletion > 0 {
            tmin = tmin.min(-c.max_target_site_deletion);
        }
        if c.allow_blunt_pairs {
            tmin = tmin.min(0);
        }
        // L1-mediated far pairs (`pair_mode`) need TPRT polarity: one side a poly-A/T
        // tail, the other a complex clip. Widening the plain partner window to the full
        // span instead kept practically every breakpoint at 30x (a partner of either
        // kind lies within 50 kb almost everywhere) and blew the mate pass up to 20-37M
        // requests / 17-30 GB. Test polarity per side: tail <-> complex partners only.
        let span = c.max_l1_mediated_span;
        let n_pa = c.l1_mediated_min_polya;
        let polar = |bps: &[Breakpoint]| {
            let (mut tail, mut cx): (FxHashMap<String, Vec<i64>>, FxHashMap<String, Vec<i64>>) = Default::default();
            if span > 0 {
                for b in bps {
                    if leading_polyt(&b.clipped.seq, n_pa) {
                        tail.entry(b.reference_name.clone()).or_default().push(b.breakpoint);
                    } else if complex_clip(&b.clipped.seq, n_pa) {
                        cx.entry(b.reference_name.clone()).or_default().push(b.breakpoint);
                    }
                }
            }
            for v in tail.values_mut().chain(cx.values_mut()) {
                v.sort_unstable();
            }
            (tail, cx)
        };
        let (l_tail, l_cx) = polar(&self.final_left_breakpoints);
        let (r_tail, r_cx) = polar(&self.final_right_breakpoints);
        let far_partner = |b: &Breakpoint, tail: &FxHashMap<String, Vec<i64>>, cx: &FxHashMap<String, Vec<i64>>| {
            if span <= 0 {
                return false;
            }
            let (rn, p) = (b.reference_name.as_str(), b.breakpoint);
            (leading_polyt(&b.clipped.seq, n_pa) && any_in(cx.get(rn), p - span, p + span))
                || (complex_clip(&b.clipped.seq, n_pa) && any_in(tail.get(rn), p - span, p + span))
        };
        let lone_ok = |b: &Breakpoint| c.one_sided_loci && one_sided_ok(c, b);
        let mut keep_l = vec![false; self.final_left_breakpoints.len()];
        let mut keep_r = vec![false; self.final_right_breakpoints.len()];
        if let Some(e) = emitted {
            // exact: the loci output emits (dry run)
            for (i, b) in self.final_left_breakpoints.iter().enumerate() {
                keep_l[i] = b.ev.is_some() && e.left[i];
            }
            for (i, b) in self.final_right_breakpoints.iter().enumerate() {
                keep_r[i] = b.ev.is_some() && e.right[i];
            }
        } else {
        for (i, b) in self.final_left_breakpoints.iter().enumerate() {
            let (rn, p) = (b.reference_name.as_str(), b.breakpoint);
            keep_l[i] = b.ev.is_some()
                && self.bp_visible(rn, p)
                && (any_in(rpos.get(rn), p + tmin, p + tmax)
                    || any_in(pa_pos.get(rn), p - slack, p + far + slack)
                    || far_partner(b, &r_tail, &r_cx)
                    || lone_ok(b));
        }
        for (i, b) in self.final_right_breakpoints.iter().enumerate() {
            let (rn, p) = (b.reference_name.as_str(), b.breakpoint);
            keep_r[i] = b.ev.is_some()
                && self.bp_visible(rn, p)
                && (any_in(lpos.get(rn), p - tmax, p - tmin)
                    || any_in(pa_pos.get(rn), p - far - slack, p + slack)
                    || far_partner(b, &l_tail, &l_cx)
                    || lone_ok(b));
        }
        }
        if self.config.evidence_sidecar {
            let kl = keep_l.iter().filter(|&&k| k).count();
            let kr = keep_r.iter().filter(|&&k| k).count();
            eprintln!(
                "evidence sidecar: kept {kl}/{} LEFT and {kr}/{} RIGHT breakpoints for the mate pass",
                keep_l.len(),
                keep_r.len()
            );
        }
        // poly-A reads that may pair with a kept breakpoint (on the mate's contig)
        let klpos = positions(&self.final_left_breakpoints, Some(&keep_l));
        let krpos = positions(&self.final_right_breakpoints, Some(&keep_r));
        let mut reqs: Vec<MateReq> = Vec::new();
        for (i, pa) in self.polya.iter_mut().enumerate() {
            let Some(l) = pa.sc else { continue };
            let want = match emitted {
                Some(e) => e.polya[i],
                None => {
                    let Some(n) = (l.mref >= 0).then(|| names.get(l.mref as usize)).flatten() else { continue };
                    let (lo, hi) = (l.mpos - far - slack, l.mpos + far + slack);
                    any_in(klpos.get(n.as_str()), lo, hi) || any_in(krpos.get(n.as_str()), lo, hi)
                }
            };
            if want {
                pa.ev = Some(Box::default());
                let hash = frag_hash(pa.qname.as_bytes());
                if pa.sc_mate_routed {
                    // the anchoring primary mate (legacy mate-pass routing, recorded)
                    reqs.push(MateReq {
                        hash,
                        idx: i as u32,
                        ref_id: -1,
                        pos: -1,
                        r12: if pa.is_read1 { 2 } else { 1 },
                        kind: ReqKind::Mate,
                        target: Target::PolyA,
                        supp: false,
                    });
                }
                reqs.push(MateReq {
                    hash,
                    idx: i as u32,
                    ref_id: l.ref_id,
                    pos: l.pos,
                    r12: crate::evidence::r12_of(l.flag),
                    kind: ReqKind::PolySelf,
                    target: Target::PolyA,
                    supp: false,
                });
            }
        }
        let (all, cap) = (self.config.fetch_all_mates, self.config.max_mates_per_breakpoint);
        let short_cap = self.config.max_short_per_breakpoint;
        for (target, bps, keep) in [
            (Target::Left, &mut self.final_left_breakpoints, &keep_l),
            (Target::Right, &mut self.final_right_breakpoints, &keep_r),
        ] {
            for (i, bp) in bps.iter_mut().enumerate() {
                if !keep[i] {
                    bp.ev = None;
                    continue;
                }
                let ref_id = name_id.get(bp.reference_name.as_str()).copied().unwrap_or(-1);
                if let Some(ev) = bp.ev.as_deref() {
                    build_requests(ev, target, i as u32, ref_id, all, cap, &mut reqs);
                    build_short_requests(ev, target, i as u32, all, short_cap, &mut reqs);
                }
            }
        }
        reqs.sort_by_key(|r| (r.hash, r.target as u8, r.idx));
        Ok(reqs)
    }

    /// TPRT sidecar pass, after the legacy mate pass: a dry run of `output` marks the
    /// breakpoints and poly-A reads of every locus it will emit, capture requests are
    /// built for exactly those, and one more linear BAM pass fetches the records. (The
    /// old capture inside the mate pass had to guess the emitted set before poly-A reads
    /// were placed; at 30x its superset kept ~87% of all breakpoints: 32M requests,
    /// 28 GB.) With `discordant_anchor` the rescue's loci are not covered by the dry
    /// run, so the conservative superset rule is used instead.
    pub fn sidecar_pass(&mut self) -> io::Result<()> {
        let emitted = if self.config.discordant_anchor {
            None
        } else {
            let mut sink = io::sink();
            let mut out_sink = io::sink();
            let mut hm: Vec<u8> = Vec::new();
            let mut sc = Sidecar {
                w: &mut sink,
                names: Vec::new(),
                stats: Default::default(),
                collect: Some(Emitted {
                    left: vec![false; self.final_left_breakpoints.len()],
                    right: vec![false; self.final_right_breakpoints.len()],
                    polya: vec![false; self.polya.len()],
                }),
            };
            self.output(&mut out_sink, &mut hm, Some(&mut sc))?;
            sc.collect
        };
        let reqs = self.sidecar_prepare(emitted.as_ref())?;
        eprintln!("evidence sidecar: {} capture requests (sidecar pass)", reqs.len());
        if reqs.is_empty() {
            return Ok(());
        }
        if is_cram(&self.filepath) {
            let ref_path = self.reference_path.clone();
            let mut reader = open_cram(&self.filepath, ref_path.as_deref())?;
            let header = reader.read_header()?;
            for result in reader.records(&header) {
                let rec = result?;
                let read = BamRead::from_record(&rec, &header)?;
                self.sidecar_capture(&read, &reqs);
            }
        } else {
            let mut reader = open_bam(&self.filepath, self.bam_threads)?;
            let header = reader.read_header()?;
            let mut record = bam::Record::default();
            while reader.read_record(&mut record)? != 0 {
                let read = BamRead::from_record(&record, &header)?;
                self.sidecar_capture(&read, &reqs);
            }
        }
        Ok(())
    }

    /// TPRT sidecar, per mate-pass record: capture a primary record matching a request
    /// (by qname hash, r12 and expected placement) as a MATE or DISC row.
    fn sidecar_capture(&mut self, read: &BamRead, reqs: &[MateReq]) {
        if read.is_secondary || read.is_qcfail {
            return;
        }
        if read.is_duplicate && !self.config.ignore_dup_flag {
            return;
        }
        let h = frag_hash(read.name_bytes());
        let mut k = reqs.partition_point(|r| r.hash < h);
        while k < reqs.len() && reqs[k].hash == h {
            let q = reqs[k];
            k += 1;
            if !q.matches(read) {
                continue;
            }
            let (role, clip_at) = match q.kind {
                ReqKind::Mate => (Role::Mate, -1),
                ReqKind::DiscSelf => (Role::Disc, -1),
                ReqKind::PolySelf => (Role::PolyA, -1),
                ReqKind::ShortSelf => {
                    // read offset of the junction, from the breakpoint coordinate + CIGAR
                    let b = if q.target == Target::Left {
                        self.final_left_breakpoints[q.idx as usize].breakpoint
                    } else {
                        self.final_right_breakpoints[q.idx as usize].breakpoint
                    };
                    (Role::Short, read.query_offset_at(b) as i32)
                }
                ReqKind::ClipSelf => {
                    // read offset of the junction in the stored (reference-forward) sequence
                    let at = if q.target == Target::Left { read.left_len } else { read.record_len().saturating_sub(read.right_len) };
                    (Role::Clip, at as i32)
                }
            };
            let r = EvRec::from_read(read, role, clip_at);
            match q.target {
                Target::PolyA => {
                    if let Some(e) = self.polya[q.idx as usize].ev.as_mut() {
                        if role == Role::Mate {
                            e.mate = Some(r);
                        } else {
                            e.read = Some(r);
                        }
                    }
                }
                Target::Left | Target::Right => {
                    let bp = if q.target == Target::Left {
                        &mut self.final_left_breakpoints[q.idx as usize]
                    } else {
                        &mut self.final_right_breakpoints[q.idx as usize]
                    };
                    if let Some(ev) = bp.ev.as_mut() {
                        match role {
                            Role::Mate => ev.mates.push(r),
                            Role::Disc => ev.disc.push(r),
                            Role::Short => ev.short.push(r),
                            _ => ev.clip.push(r),
                        }
                    }
                }
            }
        }
    }

    /// TPRT sidecar: write every read record of one emitted locus. `pa_left` /
    /// `pa_right` are the poly-A reads pooled at a poly-A end (the paired one first, then
    /// same-clip poly-A reads within `cluster_window`); empty for a breakpoint end.
    fn sidecar_locus(
        &self,
        sc: &mut Sidecar,
        left: &Emit,
        right: &Emit,
        pa_left: &[&PolyABreakpoint],
        pa_right: &[&PolyABreakpoint],
    ) -> io::Result<()> {
        if let Some(c) = sc.collect.as_mut() {
            let fl = self.final_left_breakpoints.as_ptr() as usize;
            let fr = self.final_right_breakpoints.as_ptr() as usize;
            let fp = self.polya.as_ptr() as usize;
            let bsz = std::mem::size_of::<Breakpoint>();
            let psz = std::mem::size_of::<PolyABreakpoint>();
            for (is_left, e, pas) in [(true, left, pa_left), (false, right, pa_right)] {
                match e {
                    Emit::Bp(b) => {
                        let (base, v) = if is_left { (fl, &mut c.left) } else { (fr, &mut c.right) };
                        v[(*b as *const Breakpoint as usize - base) / bsz] = true;
                    }
                    Emit::Pa(_) => {
                        for pa in pas.iter().take(self.config.max_evidence_reads_per_breakpoint) {
                            c.polya[(*pa as *const PolyABreakpoint as usize - fp) / psz] = true;
                        }
                    }
                    Emit::Disc(_) | Emit::Open(_) => {}
                }
            }
            return Ok(());
        }
        let locus = locus_name(left, right);
        sc.stats.loci += 1;
        for (side, e, pas) in [("LEFT", left, pa_left), ("RIGHT", right, pa_right)] {
            match e {
                Emit::Bp(b) => {
                    if let Some(ev) = b.ev.as_deref() {
                        for r in ev.clip.iter().chain(&ev.disc).chain(&ev.short).chain(&ev.mates) {
                            sc.row(&locus, side, r)?;
                        }
                        sc.stats.mates_per_bp.push(ev.mates.len() as u32);
                        if self.config.short_overhang_evidence {
                            sc.stats.short_per_bp.push(ev.short.len() as u32);
                        }
                    }
                }
                Emit::Pa(_) => {
                    // same caps as a breakpoint side; `pas` is in a fixed order (paired
                    // read first, then by position), so the truncation is deterministic.
                    let (mut n_mates, max_m) = (0usize, self.config.max_mates_per_breakpoint);
                    for pa in pas.iter().take(self.config.max_evidence_reads_per_breakpoint) {
                        if let Some(pe) = pa.ev.as_deref() {
                            if let Some(r) = &pe.read {
                                sc.row(&locus, side, r)?;
                            }
                            if let Some(m) = &pe.mate {
                                if n_mates < max_m {
                                    sc.row(&locus, side, m)?;
                                    n_mates += 1;
                                }
                            }
                        }
                    }
                }
                Emit::Disc(_) | Emit::Open(_) => {}
            }
        }
        Ok(())
    }

    /// Per-record body of the mate scan, shared by the BAM and CRAM loops: match a
    /// mate by qname and attach its sequence / landing site to the breakpoint (or set
    /// the polyA mate).
    #[allow(clippy::too_many_arguments)]
    fn apply_mate_record(
        &mut self,
        read: &BamRead,
        read1_mates: &FxHashSet<String>,
        read2_mates: &FxHashSet<String>,
        qmap: &FxHashMap<String, BpRef>,
        capture_dests: bool,
        min_mapq: u8,
    ) {
        if drop_read(read, self.config.ignore_dup_flag) {
            return;
        }
        // query_name is needed for every read (membership check), so decode it here.
        let qname = read.query_name();
        if read.is_read1 {
            if !read1_mates.contains(&qname) {
                return;
            }
        } else if !read2_mates.contains(&qname) {
            return;
        }
        let Some(&bpref) = qmap.get(&qname) else { return };
        match bpref {
            BpRef::PolyA(i) => {
                self.polya[i].set_mate(read, min_mapq);
                // TPRT sidecar: the anchoring PRIMARY mate of the poly-A read is captured
                // in the sidecar pass (only if this read's locus is emitted).
                if !read.is_supplementary && self.config.evidence_sidecar {
                    self.polya[i].sc_mate_routed = true;
                }
            }
            BpRef::Left(i) | BpRef::Right(i) => {
                let seq = clean_clipped_seq(&QualitySeq::new(read.seq(), read.qual()));
                let seq = if read.is_forward() { seq } else { seq.revcomp() };
                // Feature B: capture the mate's landing site for the splice check.
                let dest = (read.reference_sequence_id, read.reference_start);
                let bp = match bpref {
                    BpRef::Left(_) => &mut self.final_left_breakpoints[i],
                    BpRef::Right(_) => &mut self.final_right_breakpoints[i],
                    _ => unreachable!(),
                };
                bp.mate_seqs.push(seq);
                if capture_dests {
                    bp.mate_dests.push(dest);
                }
            }
        }
    }

    pub fn discovery(&mut self) -> io::Result<()> {
        // SPD-3 per-contig parallelism uses indexed BAM fetch; CRAM and un-indexed BAM
        // fall back to the single-threaded scan (which still gets N-way BGZF decode via
        // bam_threads). Both paths are validated byte-identical on the real WGS BAM.
        if self.config.contig_threads > 1 && !is_cram(&self.filepath) && has_bam_index(&self.filepath) {
            return self.discovery_parallel();
        }
        if self.coverage_enabled() {
            self.estimate_coverage()?;
        }
        self.extract_chimeric()?;
        if crate::mem::enabled() {
            self.mem_report("after extract (single-threaded)", &[]);
        }
        self.find_mates()?;
        // extend_mates() is a no-op in the Python (operates on the already-emptied
        // temporary_breakpoints); intentionally omitted.
        self.cluster_discordant();
        Ok(())
    }

    /// Feature A: fold `discordant_obs` into `discordant_clusters` (single-linkage
    /// within `discordant_window`, per contig+role), keeping clusters with at least
    /// `discordant_min_reads` distinct reads. No-op unless `discordant_anchor` is on.
    fn cluster_discordant(&mut self) {
        if !self.config.discordant_anchor || self.discordant_obs.is_empty() {
            return;
        }
        let window = self.config.discordant_window.max(0);
        let min_reads = self.config.discordant_min_reads.max(1);
        // group key (contig, role); within each, sort by pos and single-link.
        let mut obs = std::mem::take(&mut self.discordant_obs);
        obs.sort_by(|a, b| {
            a.contig.cmp(&b.contig).then(a.role.cmp(&b.role)).then(a.pos.cmp(&b.pos))
        });
        let mut clusters: Vec<DiscordantCluster> = Vec::new();
        let mut i = 0;
        while i < obs.len() {
            let mut j = i + 1;
            while j < obs.len()
                && obs[j].contig == obs[i].contig
                && obs[j].role == obs[i].role
                && obs[j].pos - obs[j - 1].pos <= window
            {
                j += 1;
            }
            let group = &obs[i..j];
            if group.len() >= min_reads {
                // representative position = median of the (sorted) group
                let pos = group[group.len() / 2].pos;
                let mate_dests = group.iter().map(|o| (o.mate_ref_id, o.mate_pos)).collect();
                clusters.push(DiscordantCluster {
                    contig: obs[i].contig.clone(),
                    role: obs[i].role,
                    pos,
                    n_reads: group.len(),
                    mate_dests,
                    rte_origin: 0.0,
                });
            }
            i = j;
        }
        self.stats.disc_clusters = clusters.len() as u64;

        // D3: score each cluster's mate-origin against the RTE track (reuse the SPEC-7
        // RepeatMasker loader, applied to the mate landing site). Non-gating unless
        // `discordant_rte_only`; here we only compute the fraction + a diagnostic count.
        if let Some(rte) = &self.discordant_rte {
            let names = self.reference_names().unwrap_or_default();
            let mut rejected = 0u64;
            for c in clusters.iter_mut() {
                let (mut hits, mut total) = (0usize, 0usize);
                for &(mref, mpos) in &c.mate_dests {
                    let Some(id) = mref else { continue };
                    let Some(name) = names.get(id) else { continue };
                    total += 1;
                    if rte.contains(name, mpos) {
                        hits += 1;
                    }
                }
                c.rte_origin = if total > 0 { hits as f64 / total as f64 } else { 0.0 };
                if c.rte_origin < self.config.discordant_rte_min {
                    rejected += 1;
                }
            }
            self.stats.disc_rejected_rte = rejected;
        }

        self.discordant_clusters = clusters;
    }

    /// Resolve the header's reference sequence names in id order (for the D3 mate
    /// contig lookup). Opened once per clustering pass, only when the RTE track is set.
    pub fn reference_names(&self) -> io::Result<Vec<String>> {
        let header = if is_cram(&self.filepath) {
            open_cram(&self.filepath, self.reference_path.as_deref())?.read_header()?
        } else {
            open_bam(&self.filepath, 1)?.read_header()?
        };
        Ok(header
            .reference_sequences()
            .keys()
            .map(|k| {
                let b: &[u8] = k.as_ref();
                String::from_utf8_lossy(b).into_owned()
            })
            .collect())
    }

    /// Feature A: find a discordant cluster on `contig` with the given `role` whose
    /// representative position lies in [lo, hi]. Returns the one with the most reads.
    /// Dormant since the both-sided rescue rework; retained for a possible future one-sided
    /// mode.
    #[allow(dead_code)]
    fn find_disc_cluster(&self, contig: &str, role: i32, lo: i64, hi: i64) -> Option<&DiscordantCluster> {
        self.discordant_clusters
            .iter()
            .filter(|c| c.contig == contig && c.role == role && c.pos >= lo && c.pos <= hi)
            .filter(|c| !self.config.discordant_rte_only || c.rte_origin >= self.config.discordant_rte_min)
            .max_by_key(|c| c.n_reads)
    }

    /// SPD-3: process contigs in parallel, each in its own worker over an indexed
    /// fetch, then merge in contig-id order (which equals the coordinate-sorted
    /// single-threaded order) so the breakpoint indices — and thus output — match.
    /// ⚠ Reproduces the single-threaded output on the local multi-contig BAM, but
    /// the merge/order guarantees need the real-WGS differential before this is
    /// trusted. Default (contig_threads = 1) is the validated single-threaded path.
    fn discovery_parallel(&mut self) -> io::Result<()> {
        if self.coverage_enabled() {
            self.estimate_coverage()?;
        }
        let n_contigs = {
            let mut r = bam::io::indexed_reader::Builder::default().build_from_path(&self.filepath)?;
            let h = r.read_header()?;
            h.reference_sequences().len()
        };
        let filepath = self.filepath.clone();
        let config = self.config.clone();
        let coverage = self.coverage.clone();
        let nthreads = self.config.contig_threads.min(n_contigs.max(1)).max(1);
        let next = std::sync::atomic::AtomicUsize::new(0);
        let results: std::sync::Mutex<Vec<(usize, Discovery)>> = std::sync::Mutex::new(Vec::new());
        let first_err: std::sync::Mutex<Option<io::Error>> = std::sync::Mutex::new(None);
        std::thread::scope(|s| {
            for _ in 0..nthreads {
                s.spawn(|| loop {
                    let cid = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                    if cid >= n_contigs || first_err.lock().unwrap().is_some() {
                        break;
                    }
                    match Discovery::process_contig(&filepath, cid, &config, &coverage) {
                        Ok(w) => results.lock().unwrap().push((cid, w)),
                        Err(e) => {
                            *first_err.lock().unwrap() = Some(e);
                            break;
                        }
                    }
                });
            }
        });
        if let Some(e) = first_err.into_inner().unwrap() {
            return Err(e);
        }
        let mut results = results.into_inner().unwrap();
        results.sort_by_key(|(cid, _)| *cid);
        for (_, mut w) in results {
            self.final_left_breakpoints.append(&mut w.final_left_breakpoints);
            self.final_right_breakpoints.append(&mut w.final_right_breakpoints);
            self.subfloor_left_breakpoints.append(&mut w.subfloor_left_breakpoints);
            self.subfloor_right_breakpoints.append(&mut w.subfloor_right_breakpoints);
            self.polya.append(&mut w.polya);
            self.stats.merge(&w.stats);
            self.discordant_obs.append(&mut w.discordant_obs);
        }
        if crate::mem::enabled() {
            self.mem_report("parent after merge", &[]);
        }
        self.find_mates()?;
        if crate::mem::enabled() {
            self.mem_report("parent after find_mates", &[]);
        }
        self.cluster_discordant();
        Ok(())
    }

    /// One SPD-3 worker: run extract_chimeric over a single contig and return the
    /// worker (owning that contig's breakpoints / polyA / stats).
    fn process_contig(filepath: &str, cid: usize, config: &DiscoveryConfig, coverage: &Coverage) -> io::Result<Discovery> {
        let mut w = Discovery::new(filepath.to_string(), 1, config.clone(), None, None);
        w.only_contig = Some(cid);
        w.coverage = coverage.clone();
        w.extract_chimeric()?;
        if crate::mem::enabled() {
            w.mem_report(&format!("worker cid={cid} done"), &[]);
        }
        Ok(w)
    }

    pub fn output<W: Write>(&self, writer: &mut W, hallmarks: &mut Vec<u8>, mut sc: Option<&mut Sidecar>) -> io::Result<()> {
        let cw = self.config.cluster_window;
        // SENS-5: hallmark annotation is non-gating — the FASTQ output below is
        // identical whether or not it is enabled; only this sidecar is added.
        let hm = self.config.hallmark_score;
        if hm {
            hallmarks.extend_from_slice(b"contig\tleft\tright\ttsd\tpolya_purity\ten_motif\tscore\n");
        }
        // group by reference, preserving the Python dict-key ordering of polyA
        let mut polya_order: Vec<String> = Vec::new();
        let mut polya_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        for (i, pa) in self.polya.iter().enumerate() {
            let (rn, bp) = (pa.reference_name.as_ref(), pa.breakpoint);
            let (Some(rn), Some(_)) = (rn, bp) else { continue };
            if !polya_map.contains_key(rn) {
                polya_order.push(rn.clone());
                polya_map.insert(rn.clone(), Vec::new());
            }
            polya_map.get_mut(rn).unwrap().push(i);
        }
        let mut left_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        for (i, bp) in self.final_left_breakpoints.iter().enumerate() {
            left_map.entry(bp.reference_name.clone()).or_default().push(i);
        }
        let mut right_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        for (i, bp) in self.final_right_breakpoints.iter().enumerate() {
            right_map.entry(bp.reference_name.clone()).or_default().push(i);
        }

        // union of all reference names, sorted (lexicographic, matching Python sorted())
        let mut union: Vec<String> = polya_map
            .keys()
            .chain(left_map.keys())
            .chain(right_map.keys())
            .cloned()
            .collect();
        union.sort();
        union.dedup();
        for rn in union {
            if !self.contig_ok_output(&rn) {
                continue;
            }
            if !polya_map.contains_key(&rn) {
                polya_order.push(rn.clone());
                polya_map.insert(rn, Vec::new());
            }
        }

        let empty: Vec<usize> = Vec::new();
        for rn in &polya_order {
            if !self.contig_ok_output(rn) {
                continue;
            }
            let mut l: Vec<&Breakpoint> = left_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.final_left_breakpoints[i]).collect();
            let mut r: Vec<&Breakpoint> = right_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.final_right_breakpoints[i]).collect();
            let mut p: Vec<&PolyABreakpoint> = polya_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.polya[i]).collect();
            l.sort_by_key(|b| b.breakpoint);
            r.sort_by_key(|b| b.breakpoint);
            p.sort_by_key(|b| b.breakpoint.unwrap());

            // SPEC-5: drop breakpoints inside exclude-BED regions (no-op if unset)
            if self.exclude.is_some() {
                l.retain(|b| !self.excluded(rn, b.breakpoint));
                r.retain(|b| !self.excluded(rn, b.breakpoint));
                p.retain(|b| !self.excluded(rn, b.breakpoint.unwrap()));
            }

            // SPEC-3: drop breakpoints in local-coverage pileups (> multiplier x median)
            if self.config.coverage_mask && self.coverage.median() > 0.0 {
                let thr = self.config.coverage_mask_multiplier * self.coverage.median();
                l.retain(|b| (self.coverage.local(rn, b.breakpoint) as f64) <= thr);
                r.retain(|b| (self.coverage.local(rn, b.breakpoint) as f64) <= thr);
                p.retain(|b| (self.coverage.local(rn, b.breakpoint.unwrap()) as f64) <= thr);
            }

            // SPEC-7: drop breakpoints inside a young RepeatMasker element (no-op if unset)
            if let Some(rm) = &self.rm_mask {
                l.retain(|b| !rm.contains(rn, b.breakpoint));
                r.retain(|b| !rm.contains(rn, b.breakpoint));
                p.retain(|b| !rm.contains(rn, b.breakpoint.unwrap()));
            }

            // SPEC-8: drop clip clusters that are homopolymer slippage against a reference
            // tract rather than an insertion junction. Poly-A-mate breakpoints (`p`) carry no
            // aligned side to judge, so they are untouched.
            if self.config.slippage_filter {
                let keep = |b: &&Breakpoint| {
                    !is_slippage_clip(
                        b.side,
                        &b.clipped.seq,
                        &b.unclipped.seq,
                        self.config.slippage_min_ref_run,
                        self.config.slippage_min_clip_frac,
                        self.config.slippage_max_period,
                    )
                };
                l.retain(keep);
                r.retain(keep);
            }

            // Emissions are buffered per contig (left, right, poly-A pools) and flushed in
            // order below, so the optional post-TSD pairing modes can add pairs and upgrade a
            // Bp+polyA emission. With every mode off the flushed sequence equals the legacy
            // streamed one (byte-identical).
            let want_pools = sc.is_some();
            let mut outs: Vec<Out> = Vec::new();
            // per-breakpoint use: 0 = free, 1 = in a Bp+polyA emission, 2 = in a Bp+Bp pair
            let mut used_l: Vec<u8> = vec![0; l.len()];
            let mut used_r: Vec<u8> = vec![0; r.len()];
            let mut pa_out_l: Vec<usize> = vec![usize::MAX; l.len()];
            let mut pa_out_r: Vec<usize> = vec![usize::MAX; r.len()];
            let (mut il, mut ir, mut ip) = (0usize, 0usize, 0usize);
            while il < l.len() && ir < r.len() {
                let tsd = r[ir].breakpoint - l[il].breakpoint;
                if tsd < self.config.tsd_min {
                    while ip < p.len() && p[ip].breakpoint.unwrap() - l[il].breakpoint < self.config.polya_near_dist {
                        ip += 1;
                    }
                    if ip != p.len() && p[ip].breakpoint.unwrap() - l[il].breakpoint < self.config.polya_far_dist && p[ip].clip == CLIP_RIGHT {
                        used_l[il] = 1;
                        pa_out_l[il] = outs.len();
                        outs.push(Out {
                            left: Emit::Bp(l[il]),
                            right: Emit::Pa(p[ip]),
                            pal: Vec::new(),
                            par: if want_pools { pa_pool(&p, ip, cw) } else { Vec::new() },
                            dead: false,
                        });
                        il += 1;
                        continue;
                    }
                    ir += 1;
                    continue;
                } else if tsd > self.config.tsd_max {
                    while ip < p.len() && r[ir].breakpoint - p[ip].breakpoint.unwrap() > self.config.polya_far_dist {
                        ip += 1;
                    }
                    if ip != p.len() && r[ir].breakpoint - p[ip].breakpoint.unwrap() > self.config.polya_near_dist && p[ip].clip == CLIP_LEFT {
                        used_r[ir] = 1;
                        pa_out_r[ir] = outs.len();
                        outs.push(Out {
                            left: Emit::Pa(p[ip]),
                            right: Emit::Bp(r[ir]),
                            pal: if want_pools { pa_pool(&p, ip, cw) } else { Vec::new() },
                            par: Vec::new(),
                            dead: false,
                        });
                        ir += 1;
                        continue;
                    }
                    il += 1;
                    continue;
                } else {
                    // SPEC-8b: reject a Bp+Bp insertion whose BOTH clip consensuses are
                    // homopolymer/low-complexity poly-A/T — double-sided reference-tract
                    // slippage, not a real junction. Clips are passed in the same orientation
                    // print_output emits them (left CLIPPED is revcomp'd, right CLIPPED plain),
                    // so the gate matches the emitted consensus exactly. A real MEI keeps a
                    // structured element body on one side, so it survives (one-sided spare).
                    // Both breakpoints are consumed (emitted or rejected) — advance BOTH
                    // pointers. (Previously only `il` advanced, leaving `ir` stuck on the
                    // just-paired right breakpoint; the next left breakpoint then saw a
                    // stale, too-far-left right → negative tsd → it was mis-routed into the
                    // poly-A rescue instead of pairing with its true right partner. Those
                    // mis-paired calls became poly-A-type and were dropped by combine.)
                    used_l[il] = 2;
                    used_r[ir] = 2;
                    if !self.spec8b_rejects(l[il], r[ir]) {
                        outs.push(Out { left: Emit::Bp(l[il]), right: Emit::Bp(r[ir]), pal: Vec::new(), par: Vec::new(), dead: false });
                    }
                    il += 1;
                    ir += 1;
                }
            }
            if self.config.extra_pairing() {
                self.extra_pairs(&l, &r, &mut used_l, &mut used_r, &pa_out_l, &pa_out_r, &mut outs);
            }
            for o in &outs {
                if o.dead {
                    continue;
                }
                if hm {
                    write_hallmark(hallmarks, rn, &o.left, &o.right)?;
                }
                if let Some(s) = sc.as_deref_mut() {
                    self.sidecar_locus(s, &o.left, &o.right, &o.pal, &o.par)?;
                }
                print_output(writer, o.left, o.right)?;
            }
        }
        Ok(())
    }

    /// SPEC-8b gate on a Bp+Bp pair (see `output`).
    fn spec8b_rejects(&self, lb: &Breakpoint, rb: &Breakpoint) -> bool {
        let k = self.config.clip_slippage_junction_spare;
        if k > 0 && (junction_structured(&lb.clipped.seq, k) || junction_structured(&rb.clipped.seq, k)) {
            return false;
        }
        self.config.clip_slippage_filter
            && both_clips_slippage(
                &lb.clipped.revcomp().seq,
                &rb.clipped.seq,
                self.config.clip_slippage_min_run,
                self.config.clip_slippage_max_entropy,
                self.config.clip_slippage_require_same_base,
                self.config.clip_slippage_any_base,
            )
    }

    /// Post-TSD pairing (TPRT modes, all off by default). Runs on the breakpoints the
    /// legacy loop left unpaired (free, or only in a Bp+polyA emission, which a Bp+Bp pair
    /// upgrades): (1) target-site deletion / blunt / L1-mediated pairs, greedily by mode
    /// rank, then |gap|, then position (deterministic); (2) one-sided loci for the rest.
    #[allow(clippy::too_many_arguments)]
    fn extra_pairs<'a>(
        &self,
        l: &[&'a Breakpoint],
        r: &[&'a Breakpoint],
        used_l: &mut [u8],
        used_r: &mut [u8],
        pa_out_l: &[usize],
        pa_out_r: &[usize],
        outs: &mut Vec<Out<'a>>,
    ) {
        let c = &self.config;
        for (li, ri) in greedy_extra_pairs(c, l, r, used_l, used_r, |lb, rb| self.spec8b_rejects(lb, rb)) {
            for (u, idx) in [(used_l[li], pa_out_l[li]), (used_r[ri], pa_out_r[ri])] {
                if u == 1 && idx != usize::MAX {
                    outs[idx].dead = true; // upgraded: Bp+Bp supersedes Bp+polyA
                }
            }
            used_l[li] = 2;
            used_r[ri] = 2;
            outs.push(Out { left: Emit::Bp(l[li]), right: Emit::Bp(r[ri]), pal: Vec::new(), par: Vec::new(), dead: false });
        }
        if !c.one_sided_loci {
            return;
        }
        for (side_l, bps, used) in [(true, l, &*used_l), (false, r, &*used_r)] {
            for (i, b) in bps.iter().enumerate() {
                if used[i] != 0 || !one_sided_ok(c, b) {
                    continue;
                }
                let open = Emit::Open(b.breakpoint);
                let (left, right) = if side_l { (Emit::Bp(b), open) } else { (open, Emit::Bp(b)) };
                outs.push(Out { left, right, pal: Vec::new(), par: Vec::new(), dead: false });
            }
        }
    }

    /// Feature A (both-sided rescue): append discordant-anchored calls after the normal
    /// output. Both junctions of an insertion must carry DIRECT clip evidence: a solid
    /// anchor breakpoint (>= `discordant_min_anchor_reads`, i.e. one that clustered
    /// normally) with no reciprocal real partner in its TSD window is paired with a
    /// *sub-floor* clip on the missing side — a clip cluster that fell below the normal
    /// evidence floor (>= `discordant_partner_min_reads`), typically a single read whose
    /// body maps poorly into a repetitive RTE but whose clip pins the junction, kept alive
    /// by the low-MAPQ discordant-mate rescue in `handle_record`. A discordant mate cluster
    /// alone never completes a call. Only runs when `discordant_anchor` is on; the default
    /// output above is untouched, and the "no real partner in window" test is mutually
    /// exclusive with the main loop's pairing so a real pair is never duplicated.
    pub fn discordant_rescue<W: Write>(&mut self, writer: &mut W, mut sc: Option<&mut Sidecar>) -> io::Result<()> {
        if !self.config.discordant_anchor {
            return Ok(());
        }
        if self.subfloor_left_breakpoints.is_empty() && self.subfloor_right_breakpoints.is_empty() {
            return Ok(());
        }
        let (tsd_min, tsd_max) = (self.config.tsd_min, self.config.tsd_max);
        let in_window = |gap: i64| gap >= tsd_min && gap <= tsd_max;
        let anchor_min = self.config.discordant_min_anchor_reads;
        let partner_min = self.config.discordant_partner_min_reads;

        let by_contig = |bps: &[Breakpoint]| {
            let mut m: FxHashMap<String, Vec<usize>> = FxHashMap::default();
            for (i, bp) in bps.iter().enumerate() {
                m.entry(bp.reference_name.clone()).or_default().push(i);
            }
            m
        };
        let left_map = by_contig(&self.final_left_breakpoints);
        let right_map = by_contig(&self.final_right_breakpoints);
        let sub_left_map = by_contig(&self.subfloor_left_breakpoints);
        let sub_right_map = by_contig(&self.subfloor_right_breakpoints);

        let mut contigs: Vec<String> = left_map.keys().chain(right_map.keys()).cloned().collect();
        contigs.sort();
        contigs.dedup();

        // strongest sub-floor clip of the given side within [lo, hi] (sorted-by-pos input).
        fn pick_partner<'a>(sorted: &[&'a Breakpoint], lo: i64, hi: i64, min_reads: usize) -> Option<&'a Breakpoint> {
            let start = sorted.partition_point(|b| b.breakpoint < lo);
            let mut best: Option<&'a Breakpoint> = None;
            let mut i = start;
            while i < sorted.len() && sorted[i].breakpoint <= hi {
                let b = sorted[i];
                if b.n_reads >= min_reads && best.map_or(true, |cur| b.n_reads > cur.n_reads) {
                    best = Some(b);
                }
                i += 1;
            }
            best
        }

        let empty: Vec<usize> = Vec::new();
        let mut paired: u64 = 0;
        for rn in &contigs {
            if !self.contig_ok_output(rn) {
                continue;
            }
            // solid anchors (finals), masked exactly as normal output masks them.
            let mut l: Vec<&Breakpoint> = left_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.final_left_breakpoints[i]).collect();
            let mut r: Vec<&Breakpoint> = right_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.final_right_breakpoints[i]).collect();
            l.sort_by_key(|b| b.breakpoint);
            r.sort_by_key(|b| b.breakpoint);
            self.retain_visible(rn, &mut l);
            self.retain_visible(rn, &mut r);
            // sub-floor missing-junction candidates, masked the same way.
            let mut sl: Vec<&Breakpoint> = sub_left_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.subfloor_left_breakpoints[i]).collect();
            let mut sr: Vec<&Breakpoint> = sub_right_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.subfloor_right_breakpoints[i]).collect();
            sl.sort_by_key(|b| b.breakpoint);
            sr.sort_by_key(|b| b.breakpoint);
            self.retain_visible(rn, &mut sl);
            self.retain_visible(rn, &mut sr);

            // a LEFT anchor with no real RIGHT partner in window -> a RIGHT sub-floor clip
            for lb in &l {
                if lb.n_reads < anchor_min {
                    continue;
                }
                if r.iter().any(|rb| in_window(rb.breakpoint - lb.breakpoint)) {
                    continue;
                }
                if !self.disc_coverage_ok(rn, lb.breakpoint) {
                    continue;
                }
                // Reject an anchor whose genomic *flank* (aligned side) is a low-complexity
                // satellite array: a real insertion has a unique/complex flank, whereas a
                // pericentromeric/subtelomeric mismap is satellite on both sides. (The clip
                // itself may be legitimately low-complexity — an Alu poly-A tail or SVA VNTR
                // — so the flank, not the clip, is gated.)
                if is_low_complexity(&lb.unclipped.seq, 0.8) {
                    continue;
                }
                if let Some(mb) = pick_partner(&sr, lb.breakpoint + tsd_min, lb.breakpoint + tsd_max, partner_min) {
                    if is_low_complexity(&mb.unclipped.seq, 0.8) {
                        continue;
                    }
                    if let Some(s) = sc.as_deref_mut() {
                        self.sidecar_locus(s, &Emit::Bp(lb), &Emit::Bp(mb), &[], &[])?;
                    }
                    print_output(writer, Emit::Bp(lb), Emit::Bp(mb))?;
                    paired += 1;
                }
            }
            // a RIGHT anchor with no real LEFT partner in window -> a LEFT sub-floor clip
            for rb in &r {
                if rb.n_reads < anchor_min {
                    continue;
                }
                if l.iter().any(|lb| in_window(rb.breakpoint - lb.breakpoint)) {
                    continue;
                }
                if !self.disc_coverage_ok(rn, rb.breakpoint) {
                    continue;
                }
                if is_low_complexity(&rb.unclipped.seq, 0.8) {
                    continue;
                }
                if let Some(mb) = pick_partner(&sl, rb.breakpoint - tsd_max, rb.breakpoint - tsd_min, partner_min) {
                    if is_low_complexity(&mb.unclipped.seq, 0.8) {
                        continue;
                    }
                    if let Some(s) = sc.as_deref_mut() {
                        self.sidecar_locus(s, &Emit::Bp(mb), &Emit::Bp(rb), &[], &[])?;
                    }
                    print_output(writer, Emit::Bp(mb), Emit::Bp(rb))?;
                    paired += 1;
                }
            }
        }
        self.stats.disc_paired = paired;
        Ok(())
    }

    /// Apply the SPEC-5 (exclude-BED), SPEC-3 (coverage), SPEC-7 (RM) and SPEC-8 (slippage)
    /// retains to a breakpoint list, matching `output`'s masking so a rescue never
    /// resurrects a masked locus.
    fn retain_visible(&self, rn: &str, v: &mut Vec<&Breakpoint>) {
        if self.exclude.is_some() {
            v.retain(|b| !self.excluded(rn, b.breakpoint));
        }
        if self.config.coverage_mask && self.coverage.median() > 0.0 {
            let thr = self.config.coverage_mask_multiplier * self.coverage.median();
            v.retain(|b| (self.coverage.local(rn, b.breakpoint) as f64) <= thr);
        }
        if let Some(rm) = &self.rm_mask {
            v.retain(|b| !rm.contains(rn, b.breakpoint));
        }
        if self.config.slippage_filter {
            v.retain(|b| {
                !is_slippage_clip(
                    b.side,
                    &b.clipped.seq,
                    &b.unclipped.seq,
                    self.config.slippage_min_ref_run,
                    self.config.slippage_min_clip_frac,
                    self.config.slippage_max_period,
                )
            });
        }
    }
}

// --- TPRT sidecar helpers ---

/// Read-level drop rule shared by the clip path and the mate pass: secondary and
/// QC-fail always, 0x400 duplicates unless `ignore_dup_flag`.
#[inline]
fn drop_read(read: &BamRead, ignore_dup_flag: bool) -> bool {
    read.is_secondary || read.is_qcfail || (read.is_duplicate && !ignore_dup_flag)
}

/// Attach discordant anchors (one contig, start-sorted) to that contig's new final
/// breakpoints. LEFT-clipped reads lie downstream of the insertion, so a LEFT junction at
/// B takes a REVERSE anchor starting at/after B (-5 bp) and ending within `span` (it looks
/// back into the element, like the LEFT reverse `has_mate` reads); a RIGHT junction takes
/// a FORWARD anchor ending at/before B (+5) and starting within `span`. Fragments already present as CLIP evidence are skipped; the
/// total CLIP + DISC per side is capped at `cap` (lowest fragment hash).
fn attach_disc(bps: &mut [Breakpoint], disc: &[DiscLite], left: bool, span: i64, cap: usize) {
    for bp in bps.iter_mut() {
        let b = bp.breakpoint;
        let Some(ev) = bp.ev.as_mut() else { continue };
        let lo = disc.partition_point(|d| d.start < b - span);
        let mut sel: Vec<DiscLite> = Vec::new();
        for d in &disc[lo..] {
            if d.start > b + span {
                break;
            }
            let ok = if left {
                d.is_reverse() && d.start >= b - 5 && d.end - b <= span
            } else {
                !d.is_reverse() && d.end <= b + 5 && b - d.start <= span
            };
            if ok && !ev.clip_lite.iter().any(|c| c.frag == d.frag) {
                sel.push(*d);
            }
        }
        let room = cap.saturating_sub(ev.clip_lite.len());
        ev.disc_lite = select_lowest(sel, room, |d| (d.frag, d.flag));
    }
}

/// Attach short-overhang candidates (one contig, start-sorted) to that contig's new final
/// breakpoints: `short_candidate` against the consensus breakpoint, excluding fragments
/// already present as CLIP evidence with the same r12, capped at `cap` (lowest (frag, flag)).
fn attach_short(bps: &mut [Breakpoint], shorts: &[ShortLite], left: bool, window: i64, max: i64, cap: usize) {
    let reach = window.max(max) + 2000; // the longest aligned span considered
    for bp in bps.iter_mut() {
        let b = bp.breakpoint;
        let Some(ev) = bp.ev.as_mut() else { continue };
        let lo = shorts.partition_point(|r| r.start < b - reach);
        let mut sel: Vec<ShortLite> = Vec::new();
        for r in &shorts[lo..] {
            if r.start > b + window {
                break;
            }
            if short_candidate(r, b, left, window, max, MIN_CLIP_LEN as u32)
                && !ev.clip_lite.iter().any(|c| c.frag == r.frag && c.r12() == r.r12())
            {
                sel.push(*r);
            }
        }
        sel.dedup();
        ev.short_lite = select_lowest(sel, cap, |r| (r.frag, r.flag));
    }
}

/// SHORT capture requests for one breakpoint side: each short-overhang read itself
/// (ShortSelf, placement checked) and — with `fetch_all_mates` — its primary mate, capped
/// separately at `cap` (lowest (frag, r12)) and skipping mates the CLIP/DISC requests
/// already cover, so the legacy MATE selection is unchanged.
fn build_short_requests(ev: &EvExtra, target: Target, idx: u32, fetch_all: bool, cap: usize, out: &mut Vec<MateReq>) {
    if ev.short_lite.is_empty() {
        return;
    }
    let req = |hash, r12, kind, ref_id, pos| MateReq { hash, idx, ref_id, pos, r12, kind, target, supp: false };
    let mut cands: Vec<(u64, u8)> = Vec::new();
    for s in &ev.short_lite {
        out.push(req(s.frag, s.r12(), ReqKind::ShortSelf, s.ref_id, s.start));
        if fetch_all && s.flag & FLAG_PAIRED != 0 && s.mate_r12() != 0 {
            cands.push((s.frag, s.mate_r12()));
        }
    }
    cands.sort_unstable();
    cands.dedup();
    let covered = |f: u64, r: u8| {
        ev.clip_lite.iter().any(|c| c.frag == f && (c.mate_r12() == r || (c.r12() == r && c.is_primary())))
            || ev.disc_lite.iter().any(|d| d.frag == f)
    };
    cands.retain(|&(f, r)| !covered(f, r));
    for (hash, r12) in select_lowest(cands, cap, |&c| c) {
        out.push(req(hash, r12, ReqKind::Mate, -1, -1));
    }
}

/// Capture requests for one breakpoint side. CLIP reads: each recorded clipped read
/// itself (ClipSelf). Mates: the legacy `has_mate` orientation (LEFT reverse / RIGHT
/// forward) of CLIP reads, or — with `fetch_all_mates` — the mate of every CLIP (incl.
/// supplementary) and DISC read; deduplicated per (frag, r12), skipping a mate already
/// present as a primary CLIP record, capped at `cap` by lowest (frag, r12). DISC
/// anchors: each anchor itself (DiscSelf).
fn build_requests(ev: &EvExtra, target: Target, idx: u32, ref_id: i32, fetch_all: bool, cap: usize, out: &mut Vec<MateReq>) {
    let left = target == Target::Left;
    let req = |hash, r12, kind, ref_id, pos, supp| MateReq { hash, idx, ref_id, pos, r12, kind, target, supp };
    // (frag, wanted r12)
    let mut cands: Vec<(u64, u8)> = Vec::new();
    for c in &ev.clip_lite {
        out.push(req(c.frag, c.r12(), ReqKind::ClipSelf, ref_id, c.pos, c.flag & FLAG_SUPPLEMENTARY != 0));
        if c.flag & FLAG_PAIRED == 0 || c.mate_r12() == 0 {
            continue;
        }
        let reverse = c.flag & FLAG_REVERSE != 0;
        let has_mate = if left { reverse } else { !reverse };
        if fetch_all || has_mate {
            cands.push((c.frag, c.mate_r12()));
        }
    }
    for d in &ev.disc_lite {
        out.push(req(d.frag, d.r12(), ReqKind::DiscSelf, d.ref_id, d.start, false));
        if fetch_all {
            match d.r12() {
                1 => cands.push((d.frag, 2)),
                2 => cands.push((d.frag, 1)),
                _ => {}
            }
        }
    }
    cands.sort_unstable();
    cands.dedup();
    cands.retain(|&(f, r)| !ev.clip_lite.iter().any(|c| c.frag == f && c.r12() == r && c.is_primary()));
    for (hash, r12) in select_lowest(cands, cap, |&c| c) {
        out.push(req(hash, r12, ReqKind::Mate, -1, -1, false));
    }
}

/// Poly-A reads pooled at an emitted poly-A end: `p[ip]` first, then every other poly-A
/// read of the same clip side within `window` bp (sorted input).
fn pa_pool<'a>(p: &[&'a PolyABreakpoint], ip: usize, window: i64) -> Vec<&'a PolyABreakpoint> {
    let x = p[ip].breakpoint.unwrap_or(0);
    let mut v = vec![p[ip]];
    for (k, q) in p.iter().enumerate() {
        if k != ip && q.clip == p[ip].clip && (q.breakpoint.unwrap_or(i64::MIN / 2) - x).abs() <= window {
            v.push(q);
        }
    }
    v
}

// --- SENS-5 hallmark features (non-gating annotation) ---

fn emit_clip<'a>(e: &'a Emit) -> Option<&'a QualitySeq> {
    match e {
        Emit::Bp(b) => Some(&b.clipped),
        Emit::Pa(p) => p.clipped.as_ref(),
        Emit::Disc(_) | Emit::Open(_) => None,
    }
}

fn emit_unclip<'a>(e: &'a Emit) -> Option<&'a QualitySeq> {
    match e {
        Emit::Bp(b) => Some(&b.unclipped),
        Emit::Pa(_) | Emit::Disc(_) | Emit::Open(_) => None,
    }
}

/// Poly-A/T purity over the terminal (<=12 bp) window of a clip: max(fracA, fracT).
fn terminal_purity(qs: &QualitySeq) -> f64 {
    let n = qs.seq.len();
    if n == 0 {
        return 0.0;
    }
    let w = n.min(12);
    let tail = &qs.seq[n - w..];
    let a = tail.iter().filter(|&&b| b == b'A' || b == b'a').count();
    let t = tail.iter().filter(|&&b| b == b'T' || b == b't').count();
    a.max(t) as f64 / w as f64
}

fn contains_motif(qs: &QualitySeq, motif: &[u8]) -> bool {
    let up: Vec<u8> = qs.seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    up.windows(motif.len()).any(|w| w == motif)
}

/// Compute and write one hallmark TSV line for an emitted pair. Element-class-aware
/// and never gating: poly-A is scored, never required (so ERV is not penalised); the
/// EN motif is only a small tie-breaker.
fn write_hallmark(buf: &mut Vec<u8>, contig: &str, left: &Emit, right: &Emit) -> io::Result<()> {
    let epos = |e: &Emit| -> i64 {
        match e {
            Emit::Bp(b) => b.breakpoint,
            Emit::Pa(p) => p.breakpoint.unwrap_or(0),
            Emit::Disc(c) => c.pos,
            Emit::Open(p) => *p,
        }
    };
    let (l_pos, r_pos): (i64, i64) = (epos(left), epos(right));
    let tsd = r_pos - l_pos;
    let mut purity = 0.0f64;
    for c in [emit_clip(left), emit_clip(right)].into_iter().flatten() {
        purity = purity.max(terminal_purity(c));
    }
    // one fixed EN motif (L1 endonuclease consensus), checked in the genomic flank
    let en = [emit_unclip(left), emit_unclip(right)]
        .into_iter()
        .flatten()
        .any(|u| contains_motif(u, b"TTAAAA"));
    let tsd_reward = if (2..=20).contains(&tsd) { 1.0 } else { 0.0 };
    let score = purity + tsd_reward + if en { 0.1 } else { 0.0 };
    writeln!(buf, "{contig}\t{l_pos}\t{r_pos}\t{tsd}\t{purity:.3}\t{en}\t{score:.3}")
}

/// Locus id `contig:L-R` exactly as in the `.txt.gz` record names.
fn locus_name(left: &Emit, right: &Emit) -> String {
    let (reference_name, left_str, right_str) = match (left, right) {
        (Emit::Pa(pa), Emit::Bp(rb)) => (
            pa.reference_name.clone().unwrap(),
            format!("polyA_{}", pa.breakpoint.unwrap()),
            format!("{}", rb.breakpoint),
        ),
        (Emit::Bp(lb), Emit::Pa(pa)) => (
            pa.reference_name.clone().unwrap(),
            format!("{}", lb.breakpoint),
            format!("polyA_{}", pa.breakpoint.unwrap()),
        ),
        (Emit::Bp(lb), Emit::Bp(rb)) => (
            lb.reference_name.clone(),
            format!("{}", lb.breakpoint),
            format!("{}", rb.breakpoint),
        ),
        // Feature A: a discordant end carries a `disc_<pos>` token and no reads.
        (Emit::Bp(lb), Emit::Disc(c)) => (
            lb.reference_name.clone(),
            format!("{}", lb.breakpoint),
            format!("disc_{}", c.pos),
        ),
        (Emit::Disc(c), Emit::Bp(rb)) => (
            rb.reference_name.clone(),
            format!("disc_{}", c.pos),
            format!("{}", rb.breakpoint),
        ),
        // TPRT one-sided locus: the missing end repeats the real coordinate.
        (Emit::Bp(lb), Emit::Open(p)) => (lb.reference_name.clone(), format!("{}", lb.breakpoint), format!("oneside_{p}")),
        (Emit::Open(p), Emit::Bp(rb)) => (rb.reference_name.clone(), format!("oneside_{p}"), format!("{}", rb.breakpoint)),
        _ => unreachable!("invalid pairing: polyA/disc ends are only paired with a real breakpoint"),
    };
    format!("{}:{}-{}", reference_name, left_str, right_str)
}

fn print_output<W: Write>(writer: &mut W, left: Emit, right: Emit) -> io::Result<()> {
    let bp_name = locus_name(&left, &right);

    match &left {
        Emit::Pa(pa) => {
            let c = pa.clipped.as_ref().unwrap().revcomp();
            writer.write_all(c.fastq(&format!("{}:LEFT:CLIPPED_POLYA", bp_name)).as_bytes())?;
        }
        Emit::Bp(lb) => {
            writer.write_all(lb.clipped.revcomp().fastq(&format!("{}:LEFT:CLIPPED", bp_name)).as_bytes())?;
            writer.write_all(lb.unclipped.fastq(&format!("{}:LEFT:ALIGNED", bp_name)).as_bytes())?;
            for (i, mate) in lb.mate_seqs.iter().enumerate() {
                writer.write_all(mate.revcomp().fastq(&format!("{}:LEFT:MATE{}", bp_name, i)).as_bytes())?;
            }
        }
        // Feature A: discordant left end emits no reads (coordinate only).
        Emit::Disc(_) | Emit::Open(_) => {}
    }
    match &right {
        Emit::Pa(pa) => {
            let c = pa.clipped.as_ref().unwrap();
            writer.write_all(c.fastq(&format!("{}:RIGHT:CLIPPED_POLYA", bp_name)).as_bytes())?;
        }
        Emit::Bp(rb) => {
            writer.write_all(rb.unclipped.revcomp().fastq(&format!("{}:RIGHT:ALIGNED", bp_name)).as_bytes())?;
            writer.write_all(rb.clipped.fastq(&format!("{}:RIGHT:CLIPPED", bp_name)).as_bytes())?;
            for (i, mate) in rb.mate_seqs.iter().enumerate() {
                writer.write_all(mate.revcomp().fastq(&format!("{}:RIGHT:MATE{}", bp_name, i)).as_bytes())?;
            }
        }
        // Feature A: discordant right end emits no reads (coordinate only).
        Emit::Disc(_) | Emit::Open(_) => {}
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::ClipLite;
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::RecordBuf;

    fn read_with_flags(flags: u16) -> RecordBuf {
        RecordBuf::builder()
            .set_flags(Flags::from_bits_truncate(flags))
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(100).unwrap())
            .set_cigar(vec![Op::new(Kind::Match, 4)].into())
            .set_sequence(b"ACGT".to_vec().into())
            .build()
    }

    #[test]
    fn dup_flag_handling_both_settings() {
        let h = Header::default();
        let dup = read_with_flags(0x1 | 0x40 | 0x400);
        let r = BamRead::from_record(&dup, &h).unwrap();
        assert!(drop_read(&r, false), "legacy: 0x400 dropped");
        assert!(!drop_read(&r, true), "ignore_dup_flag: 0x400 kept");
        let plain = read_with_flags(0x1 | 0x40);
        let r = BamRead::from_record(&plain, &h).unwrap();
        assert!(!drop_read(&r, false) && !drop_read(&r, true));
        // secondary / QC-fail are dropped regardless of the dup setting
        for f in [0x100u16, 0x200] {
            let rec = read_with_flags(0x1 | f | 0x400);
            let r = BamRead::from_record(&rec, &h).unwrap();
            assert!(drop_read(&r, false) && drop_read(&r, true));
        }
    }

    fn clip_lite(i: u64, flag: u16) -> ClipLite {
        ClipLite { frag: frag_hash(format!("q{i}").as_bytes()), flag, pos: 1000 }
    }

    fn mates(v: &[MateReq]) -> Vec<MateReq> {
        v.iter().filter(|r| r.kind == ReqKind::Mate).copied().collect()
    }

    #[test]
    fn mate_cap_is_deterministic_and_respects_orientation() {
        // 120 forward RIGHT clip reads (has_mate orientation) + 40 reverse ones
        let mut recs: Vec<ClipLite> = (0..120).map(|i| clip_lite(i, 0x1 | 0x40)).collect();
        recs.extend((120..160).map(|i| clip_lite(i, 0x1 | 0x40 | 0x10)));
        let ev = EvExtra { clip_lite: recs.clone(), ..Default::default() };
        let mut all_reqs = Vec::new();
        build_requests(&ev, Target::Right, 7, 0, false, 50, &mut all_reqs);
        // every clipped read is itself requested (ClipSelf, breakpoint coordinate kept)
        assert_eq!(all_reqs.iter().filter(|r| r.kind == ReqKind::ClipSelf && r.pos == 1000).count(), 160);
        let a = mates(&all_reqs);
        assert_eq!(a.len(), 50);
        assert!(a.iter().all(|r| r.r12 == 2 && r.idx == 7 && r.target == Target::Right));
        // has_mate-only: none of the reverse reads' mates were requested
        let rev: Vec<u64> = recs[120..].iter().map(|r| r.frag).collect();
        assert!(a.iter().all(|r| !rev.contains(&r.hash)));
        // same subset for a permuted input
        recs.reverse();
        recs.rotate_left(33);
        let ev2 = EvExtra { clip_lite: recs.clone(), ..Default::default() };
        let mut b = Vec::new();
        build_requests(&ev2, Target::Right, 7, 0, false, 50, &mut b);
        let key = |v: &[MateReq]| {
            let mut k: Vec<(u64, u8)> = v.iter().map(|r| (r.hash, r.r12)).collect();
            k.sort();
            k
        };
        assert_eq!(key(&a), key(&mates(&b)));
        // the subset is the 50 lowest fragment hashes of the eligible reads
        let mut elig: Vec<u64> = recs.iter().filter(|r| r.flag & 0x10 == 0).map(|r| r.frag).collect();
        elig.sort();
        let mut got: Vec<u64> = a.iter().map(|r| r.hash).collect();
        got.sort();
        assert_eq!(got, elig[..50].to_vec());
        // fetch_all_mates: reverse reads' mates become eligible too; cap still holds
        let mut c = Vec::new();
        build_requests(&ev2, Target::Right, 7, 0, true, 50, &mut c);
        let c = mates(&c);
        assert_eq!(c.len(), 50);
        let mut all: Vec<u64> = recs.iter().map(|r| r.frag).collect();
        all.sort();
        let mut got: Vec<u64> = c.iter().map(|r| r.hash).collect();
        got.sort();
        assert_eq!(got, all[..50].to_vec());
    }

    #[test]
    fn mate_already_present_as_clip_is_not_requested() {
        // both mates of fragment q0 are CLIP records: no MATE request for it, but the
        // supplementary of q1 still asks for q1's mate (the supplementary is not a mate)
        let ev = EvExtra {
            clip_lite: vec![clip_lite(0, 0x1 | 0x40), clip_lite(0, 0x1 | 0x80 | 0x10), clip_lite(1, 0x1 | 0x40 | 0x800)],
            ..Default::default()
        };
        let mut v = Vec::new();
        build_requests(&ev, Target::Right, 0, 0, true, 50, &mut v);
        let m = mates(&v);
        assert_eq!(m.len(), 1);
        assert_eq!((m[0].hash, m[0].r12), (frag_hash(b"q1"), 2));
        // the supplementary self-request demands a supplementary record
        assert!(v.iter().any(|r| r.kind == ReqKind::ClipSelf && r.supp && r.hash == frag_hash(b"q1")));
    }

    #[test]
    fn disc_anchor_sides() {
        // LEFT junction at 1000 takes a reverse anchor downstream; RIGHT a forward upstream
        let mk_bp = |side: i32| {
            let mut b = Breakpoint::new(side, "1".into(), 1000, None, QualitySeq::empty(), QualitySeq::empty(), None, None, false, 0);
            b.ev = Some(Box::new(EvExtra::default()));
            b
        };
        let fwd_up = DiscLite { frag: 1, flag: 0x1 | 0x40, ref_id: 0, start: 800, end: 950 };
        let rev_down = DiscLite { frag: 2, flag: 0x1 | 0x10 | 0x80, ref_id: 0, start: 1050, end: 1200 };
        let far = DiscLite { frag: 3, flag: 0x1 | 0x10 | 0x80, ref_id: 0, start: 1600, end: 1750 };
        let disc = vec![fwd_up, rev_down, far];
        let mut l = vec![mk_bp(CLIP_LEFT)];
        let mut r = vec![mk_bp(CLIP_RIGHT)];
        attach_disc(&mut l, &disc, true, 500, 200);
        attach_disc(&mut r, &disc, false, 500, 200);
        let lf: Vec<u64> = l[0].ev.as_ref().unwrap().disc_lite.iter().map(|d| d.frag).collect();
        let rf: Vec<u64> = r[0].ev.as_ref().unwrap().disc_lite.iter().map(|d| d.frag).collect();
        assert_eq!(lf, vec![2]);
        assert_eq!(rf, vec![1]);
    }

    // ---- TPRT pairing modes ----

    /// deterministic pseudo-random DNA (LCG) — structured, high-entropy clip/flank
    fn dna(seed: u64, n: usize) -> Vec<u8> {
        let mut x = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    /// final (consensus) breakpoint with a STORED clip (junction-outward) and flank
    fn fbp(side: i32, pos: i64, clip: &[u8], frags: usize) -> Breakpoint {
        let flank = dna(pos as u64 + 7, 80);
        let mut b = Breakpoint::new(
            side,
            "1".into(),
            pos,
            None,
            QualitySeq::new(clip.to_vec(), vec![30; clip.len()]),
            QualitySeq::new(flank, vec![30; 80]),
            None,
            None,
            false,
            0,
        );
        b.n_frags = frags;
        b.n_reads = frags;
        b
    }

    fn polyt(n: usize, tail_seed: u64) -> Vec<u8> {
        let mut v = vec![b'T'; n];
        v.extend(dna(tail_seed, 20));
        v
    }

    fn cfg_modes() -> DiscoveryConfig {
        DiscoveryConfig {
            max_target_site_deletion: 30,
            allow_blunt_pairs: true,
            max_l1_mediated_span: 50_000,
            one_sided_loci: true,
            ..DiscoveryConfig::default()
        }
    }

    #[test]
    fn pair_mode_windows() {
        let c = cfg_modes();
        let body = dna(1, 40);
        let tail = polyt(25, 2);
        let l = |p| fbp(CLIP_LEFT, p, &tail, 3);
        let r = |p| fbp(CLIP_RIGHT, p, &body, 3);
        // normal TSD window is the legacy loop's business
        assert_eq!(pair_mode(&c, &l(1000), &r(1015)), None);
        // target-site deletion: RIGHT up to 30 bp LEFT of LEFT
        assert_eq!(pair_mode(&c, &l(1000), &r(985)), Some(PairMode::TsdDeletion));
        assert_eq!(pair_mode(&c, &l(1000), &r(970)), Some(PairMode::TsdDeletion));
        // blunt (0..tsd_min-1)
        assert_eq!(pair_mode(&c, &l(1000), &r(1000)), Some(PairMode::Blunt));
        assert_eq!(pair_mode(&c, &l(1000), &r(1001)), Some(PairMode::Blunt));
        // polarised far pairs: L1-mediated deletion / duplication
        assert_eq!(pair_mode(&c, &l(5000), &r(1000)), Some(PairMode::L1Del));
        assert_eq!(pair_mode(&c, &l(1000), &r(4000)), Some(PairMode::L1Dup));
        assert_eq!(pair_mode(&c, &l(1000), &r(60_000)), None); // beyond the span
        // the same far pairs WITHOUT polarity (both complex, or both poly-A) are refused
        let lc = fbp(CLIP_LEFT, 5000, &dna(3, 40), 3);
        assert_eq!(pair_mode(&c, &lc, &r(1000)), None);
        let rt = fbp(CLIP_RIGHT, 4000, &polyt(25, 4), 3);
        assert_eq!(pair_mode(&c, &l(1000), &rt), None);
        // modes off -> nothing
        let off = DiscoveryConfig::default();
        assert!(!off.extra_pairing());
        for (a, b) in [(1000, 985), (1000, 1000), (5000, 1000), (1000, 4000)] {
            assert_eq!(pair_mode(&off, &l(a), &r(b)), None);
        }
    }

    #[test]
    fn greedy_pairs_prefer_small_gaps_and_are_one_to_one() {
        let c = cfg_modes();
        let tail = polyt(25, 2);
        let body = dna(1, 40);
        let l1 = fbp(CLIP_LEFT, 1000, &tail, 3);
        let l2 = fbp(CLIP_LEFT, 3000, &tail, 3);
        let r1 = fbp(CLIP_RIGHT, 990, &body, 3); // TSD deletion partner of l1 (gap -10)
        let r2 = fbp(CLIP_RIGHT, 2000, &body, 3); // L1 dup for l1 (+1000) / L1 del for l2 (-1000)
        let l = vec![&l1, &l2];
        let r = vec![&r1, &r2];
        let pairs = greedy_extra_pairs(&c, &l, &r, &[0, 0], &[0, 0], |_, _| false);
        assert_eq!(pairs, vec![(0, 0), (1, 1)]);
        // an already-paired (used == 2) breakpoint is never re-used
        let pairs = greedy_extra_pairs(&c, &l, &r, &[2, 0], &[0, 0], |_, _| false);
        assert_eq!(pairs, vec![(1, 1)]);
        // a breakpoint only in a Bp+polyA emission (used == 1) can be upgraded
        let pairs = greedy_extra_pairs(&c, &l, &r, &[1, 0], &[0, 0], |_, _| false);
        assert_eq!(pairs[0], (0, 0));
        // the SPEC-8b gate vetoes a pair
        let pairs = greedy_extra_pairs(&c, &l, &r, &[0, 0], &[0, 0], |a, _| a.breakpoint == 1000);
        assert_eq!(pairs, vec![(1, 1)]);
    }

    #[test]
    fn one_sided_gate() {
        let c = cfg_modes();
        // poly-A tail junction with 2 fragments passes
        assert!(one_sided_ok(&c, &fbp(CLIP_LEFT, 1000, &polyt(20, 5), 2)));
        // fragment floor
        assert!(!one_sided_ok(&c, &fbp(CLIP_LEFT, 1000, &polyt(20, 5), 1)));
        // complex clip needs one_sided_require_polya = false
        let cx = fbp(CLIP_RIGHT, 1000, &dna(9, 40), 3);
        assert!(!one_sided_ok(&c, &cx));
        assert!(one_sided_ok(&DiscoveryConfig { one_sided_require_polya: false, ..cfg_modes() }, &cx));
        // SPEC-8 reference-tract slippage is always applied: LEFT clip T-run continuing a
        // reference T-run at the junction (aligned side starts with T x 12)
        let mut slip = fbp(CLIP_LEFT, 1000, &polyt(20, 5), 4);
        let mut flank = vec![b'T'; 12];
        flank.extend(dna(11, 60));
        slip.unclipped = QualitySeq::new(flank, vec![30; 72]);
        assert!(!one_sided_ok(&c, &slip));
    }

    #[test]
    fn junction_spare_and_polyt_detection() {
        // tolerant poly-T start (2 interruptions in 10)
        assert!(leading_polyt(b"TTTTATTTCTTTTTTT", 10));
        assert!(!leading_polyt(b"TTTAATTCCTTTTTTT", 10));
        assert!(!leading_polyt(b"TTTT", 10));
        // structured junction-proximal prefix vs homopolymer at the junction
        let orphan = [b"GGGGGCTGCGCTAGTCGCATCAAAACTAAG".as_slice(), &[b'A'; 40]].concat();
        assert!(junction_structured(&orphan, 20));
        assert!(!junction_structured(&[b'T'; 40], 20));
        assert!(!junction_structured(&orphan, 0)); // off
        // SPEC-8b spare only with the key on
        let mut d = Discovery::new(String::new(), 1, DiscoveryConfig::default(), None, None);
        let lb = fbp(CLIP_LEFT, 1000, &orphan, 3);
        let rb = fbp(CLIP_RIGHT, 1013, &[b'T'; 60], 3);
        let legacy = d.spec8b_rejects(&lb, &rb);
        d.config.clip_slippage_junction_spare = 20;
        assert!(!d.spec8b_rejects(&lb, &rb));
        // a genuinely double-homopolymer pair is still rejected with the spare on
        let lh = fbp(CLIP_LEFT, 1000, &[b'T'; 60], 3);
        assert!(d.spec8b_rejects(&lh, &rb));
        let _ = legacy;
    }

    #[test]
    fn oneside_locus_names() {
        let b = fbp(CLIP_LEFT, 5000, &polyt(20, 1), 2);
        assert_eq!(locus_name(&Emit::Bp(&b), &Emit::Open(5000)), "1:5000-oneside_5000");
        let r = fbp(CLIP_RIGHT, 7000, &polyt(20, 1), 2);
        assert_eq!(locus_name(&Emit::Open(7000), &Emit::Bp(&r)), "1:oneside_7000-7000");
        // non-TSD Bp+Bp pairs keep plain numeric names (deletion: left > right)
        let l = fbp(CLIP_LEFT, 5015, &polyt(20, 1), 2);
        let rr = fbp(CLIP_RIGHT, 5000, &dna(2, 30), 2);
        assert_eq!(locus_name(&Emit::Bp(&l), &Emit::Bp(&rr)), "1:5015-5000");
        let mut out = Vec::new();
        print_output(&mut out, Emit::Bp(&b), Emit::Open(5000)).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("@1:5000-oneside_5000:LEFT:CLIPPED") && s.contains(":LEFT:ALIGNED") && !s.contains("RIGHT"));
    }

    // ---- SHORT overhang ----

    fn short_rec(flags: u16, start: usize, cigar: Vec<Op>, len: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_flags(Flags::from_bits_truncate(flags))
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(start + 1).unwrap())
            .set_cigar(cigar.into())
            .set_sequence(dna(start as u64, len).into())
            .build()
    }

    #[test]
    fn short_clip_at_from_cigar_both_sides_both_strands() {
        let h = Header::default();
        for strand in [0u16, 0x10] {
            // (a) LEFT: 5S95M at 1000 -> junction at b=1000 is offset 5
            let rec = short_rec(0x1 | 0x40 | strand, 1000, vec![Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 95)], 100);
            let r = BamRead::from_record(&rec, &h).unwrap();
            assert_eq!(r.query_offset_at(1000), 5);
            assert_eq!(r.query_offset_at(1002), 7); // window +-3
            assert_eq!(r.query_offset_at(998), 3);
            // (b) LEFT: unclipped 100M at 990 crossing b=1000 by 10 bp -> offset 10
            let rec = short_rec(0x1 | 0x80 | strand, 990, vec![Op::new(Kind::Match, 100)], 100);
            let r = BamRead::from_record(&rec, &h).unwrap();
            assert_eq!(r.query_offset_at(1000), 10);
            // (a) RIGHT: 93M7S ending at 1000 -> first clipped base offset 93
            let rec = short_rec(0x1 | 0x40 | strand, 907, vec![Op::new(Kind::Match, 93), Op::new(Kind::SoftClip, 7)], 100);
            let r = BamRead::from_record(&rec, &h).unwrap();
            assert_eq!(r.reference_end, 1000);
            assert_eq!(r.query_offset_at(1000), 93);
            // (b) RIGHT: 50M2I48M from 915 ends at 1013; b=1000 -> 85 + 2 (insertion) = 87
            let rec = short_rec(0x1 | 0x80 | strand, 915, vec![Op::new(Kind::Match, 50), Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 48)], 100);
            let r = BamRead::from_record(&rec, &h).unwrap();
            assert_eq!(r.reference_end, 1013);
            assert_eq!(r.query_offset_at(1000), 87);
            // deletion before b shifts the offset back
            let rec = short_rec(0x1 | 0x80 | strand, 900, vec![Op::new(Kind::Match, 50), Op::new(Kind::Deletion, 4), Op::new(Kind::Match, 50)], 100);
            let r = BamRead::from_record(&rec, &h).unwrap();
            assert_eq!(r.query_offset_at(1000), 96);
        }
    }

    fn slite(i: u64, start: i64, end: i64, lead: u32, trail: u32) -> ShortLite {
        ShortLite { frag: frag_hash(format!("s{i}").as_bytes()), flag: 0x1 | 0x40, ref_id: 0, start, end, lead_soft: lead, trail_soft: trail }
    }

    #[test]
    fn short_attach_cap_is_deterministic_and_skips_clip_frags() {
        let mk = |side: i32| {
            let mut b = Breakpoint::new(side, "1".into(), 1000, None, QualitySeq::empty(), QualitySeq::empty(), None, None, false, 0);
            b.ev = Some(Box::new(EvExtra::default()));
            b
        };
        // 150 LEFT (b) candidates + one that is a CLIP read of the cluster (same frag/r12)
        let mut shorts: Vec<ShortLite> = (0..150).map(|i| slite(i, 985 + (i % 10) as i64, 1135, 0, 0)).collect();
        shorts.push(slite(999, 1000, 1140, 4, 0));
        shorts.sort_by_key(|r| (r.start, r.frag, r.flag));
        let mut l = vec![mk(CLIP_LEFT)];
        l[0].ev.as_mut().unwrap().clip_lite.push(ClipLite { frag: frag_hash(b"s999"), flag: 0x1 | 0x40, pos: 1000 });
        attach_short(&mut l, &shorts, true, 3, 20, 100);
        let got: Vec<u64> = l[0].ev.as_ref().unwrap().short_lite.iter().map(|r| r.frag).collect();
        assert_eq!(got.len(), 100);
        assert!(!got.contains(&frag_hash(b"s999")));
        let mut want: Vec<u64> = (0..150).map(|i| frag_hash(format!("s{i}").as_bytes())).collect();
        want.sort();
        let mut g = got.clone();
        g.sort();
        assert_eq!(g, want[..100].to_vec()); // the 100 lowest fragment hashes
        // permuted input -> same subset
        let mut perm = shorts.clone();
        perm.reverse();
        perm.sort_by_key(|r| (r.start, r.frag, r.flag));
        let mut l2 = vec![mk(CLIP_LEFT)];
        attach_short(&mut l2, &perm, true, 3, 20, 100);
        let mut g2: Vec<u64> = l2[0].ev.as_ref().unwrap().short_lite.iter().map(|r| r.frag).filter(|f| *f != frag_hash(b"s999")).collect();
        g2.sort();
        assert_eq!(g2[..99], g[..99]);
        // RIGHT side: none of the LEFT candidates qualify
        let mut r = vec![mk(CLIP_RIGHT)];
        attach_short(&mut r, &shorts, false, 3, 20, 100);
        assert!(r[0].ev.as_ref().unwrap().short_lite.is_empty());
        // requests: ShortSelf per read + capped mates (fetch_all_mates)
        let mut reqs = Vec::new();
        build_short_requests(l[0].ev.as_ref().unwrap(), Target::Left, 3, true, 40, &mut reqs);
        assert_eq!(reqs.iter().filter(|q| q.kind == ReqKind::ShortSelf).count(), 100);
        assert_eq!(reqs.iter().filter(|q| q.kind == ReqKind::Mate).count(), 40);
        let mut reqs = Vec::new();
        build_short_requests(l[0].ev.as_ref().unwrap(), Target::Left, 3, false, 40, &mut reqs);
        assert_eq!(reqs.iter().filter(|q| q.kind == ReqKind::Mate).count(), 0);
    }
}
