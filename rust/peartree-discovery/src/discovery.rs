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
use crate::exons::GeneModel;
use crate::filters::{clean_clipped_seq, is_adapter, is_low_complexity, is_slippage_clip, mean_kmer_diversity};
use crate::intervals::IntervalIndex;
use crate::model::{join, Breakpoint};
use crate::polya::PolyABreakpoint;
use crate::qseq::{revcomp_bytes, QualitySeq};
use crate::read::BamRead;
use crate::stats::Stats;

#[derive(Clone, Copy)]
enum BpRef {
    Left(usize),
    Right(usize),
    PolyA(usize),
}

enum Emit<'a> {
    Bp(&'a Breakpoint),
    Pa(&'a PolyABreakpoint),
    /// Feature A: a discordant-cluster end (coordinate only, no reads).
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
    /// evidence of low complexity.
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
        Ok(())
    }

    /// SPEC-4: evidence floor for a breakpoint at `pos`, scaled by local/median
    /// coverage but never below the base `min_evidence_reads_per_breakpoint`.
    fn evidence_floor(&self, pos: i64) -> usize {
        let base = self.config.min_evidence_reads_per_breakpoint;
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
    ) {
        let bp = Breakpoint::new(
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
        if self.reference_name.is_none() {
            return;
        }
        if self.temporary_breakpoints.is_empty() {
            return;
        }
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
                if let Some(joined) = join(g, &self.config, floor, &mut self.stats) {
                    if out {
                        self.final_left_breakpoints.push(joined);
                    } else {
                        self.final_right_breakpoints.push(joined);
                    }
                }
            }
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
            // mate maps uniquely (proper pair, MQ >= min_mapq). Mirrors Python
            // src/discovery.py. Otherwise the read only survives as a polyA mate.
            let has_clip = (read.left_is_soft && read.left_len >= MIN_CLIP_LEN)
                || (read.right_is_soft && read.right_len >= MIN_CLIP_LEN);
            let rescued = self.config.mate_anchor_rescue
                && has_clip
                && read.is_proper_pair
                && read.mate_is_mapped
                && read.mq().map_or(false, |mq| mq >= min_mapq);
            if !rescued {
                if let Some(b) = PolyABreakpoint::find_polya(&read) {
                    self.polya.push(b);
                }
                return Ok(());
            }
            // else: fall through and process this clip as a breakpoint
        }
        if read.is_secondary || read.is_qcfail || read.is_duplicate {
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

        if !read.has_cigar {
            return Ok(());
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
            self.add_breakpoint(clip, read.reference_start, &qname, clipped, unclipped, read.is_read1, read.is_forward(), exclude_flag, read.mapq);
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
            self.add_breakpoint(clip, read.reference_end, &qname, clipped, unclipped, read.is_read1, read.is_forward(), exclude_flag, read.mapq);
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
        if read.is_secondary || read.is_qcfail || read.is_duplicate {
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
    fn reference_names(&self) -> io::Result<Vec<String>> {
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

    pub fn output<W: Write>(&self, writer: &mut W, hallmarks: &mut Vec<u8>) -> io::Result<()> {
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

            let (mut il, mut ir, mut ip) = (0usize, 0usize, 0usize);
            while il < l.len() && ir < r.len() {
                let tsd = r[ir].breakpoint - l[il].breakpoint;
                if tsd < self.config.tsd_min {
                    while ip < p.len() && p[ip].breakpoint.unwrap() - l[il].breakpoint < self.config.polya_near_dist {
                        ip += 1;
                    }
                    if ip != p.len() && p[ip].breakpoint.unwrap() - l[il].breakpoint < self.config.polya_far_dist && p[ip].clip == CLIP_RIGHT {
                        if hm {
                            write_hallmark(hallmarks, rn, &Emit::Bp(l[il]), &Emit::Pa(p[ip]))?;
                        }
                        print_output(writer, Emit::Bp(l[il]), Emit::Pa(p[ip]))?;
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
                        if hm {
                            write_hallmark(hallmarks, rn, &Emit::Pa(p[ip]), &Emit::Bp(r[ir]))?;
                        }
                        print_output(writer, Emit::Pa(p[ip]), Emit::Bp(r[ir]))?;
                        ir += 1;
                        continue;
                    }
                    il += 1;
                    continue;
                } else {
                    if hm {
                        write_hallmark(hallmarks, rn, &Emit::Bp(l[il]), &Emit::Bp(r[ir]))?;
                    }
                    print_output(writer, Emit::Bp(l[il]), Emit::Bp(r[ir]))?;
                    // Both breakpoints are consumed by this insertion — advance BOTH
                    // pointers. (Previously only `il` advanced, leaving `ir` stuck on the
                    // just-paired right breakpoint; the next left breakpoint then saw a
                    // stale, too-far-left right → negative tsd → it was mis-routed into the
                    // poly-A rescue instead of pairing with its true right partner. Those
                    // mis-paired calls became poly-A-type and were dropped by combine.)
                    il += 1;
                    ir += 1;
                }
            }
        }
        Ok(())
    }

    /// Feature A: append discordant-anchored calls after the normal output. A real
    /// breakpoint with no reciprocal real partner in its TSD window is paired with a
    /// discordant cluster on the missing side (`disc_<pos>` token, no reads for that
    /// end). Only runs when `discordant_anchor` is on, so the default output above is
    /// untouched. The "no real partner in window" test is mutually exclusive with the
    /// main loop's pairing, so a real pair is never duplicated.
    pub fn discordant_rescue<W: Write>(&mut self, writer: &mut W) -> io::Result<()> {
        if !self.config.discordant_anchor || self.discordant_clusters.is_empty() {
            return Ok(());
        }
        let (tsd_min, tsd_max) = (self.config.tsd_min, self.config.tsd_max);
        let in_window = |gap: i64| gap >= tsd_min && gap <= tsd_max;
        // Discordant anchors sit up to ~a fragment length from the junction, so the
        // cluster search reaches `discordant_rescue_span` (falling back to tsd_max),
        // decoupled from the TSD bound used to detect a *real* reciprocal partner.
        let span = self.config.discordant_rescue_span.unwrap_or(tsd_max);

        // group breakpoints per contig, applying the same SPEC-5/3/7 retains as output.
        let mut left_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        for (i, bp) in self.final_left_breakpoints.iter().enumerate() {
            left_map.entry(bp.reference_name.clone()).or_default().push(i);
        }
        let mut right_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        for (i, bp) in self.final_right_breakpoints.iter().enumerate() {
            right_map.entry(bp.reference_name.clone()).or_default().push(i);
        }
        let mut contigs: Vec<String> = left_map.keys().chain(right_map.keys()).cloned().collect();
        contigs.sort();
        contigs.dedup();

        let empty: Vec<usize> = Vec::new();
        let mut paired: u64 = 0;
        for rn in &contigs {
            if !self.contig_ok_output(rn) {
                continue;
            }
            let mut l: Vec<&Breakpoint> = left_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.final_left_breakpoints[i]).collect();
            let mut r: Vec<&Breakpoint> = right_map.get(rn).unwrap_or(&empty).iter().map(|&i| &self.final_right_breakpoints[i]).collect();
            l.sort_by_key(|b| b.breakpoint);
            r.sort_by_key(|b| b.breakpoint);
            self.retain_visible(rn, &mut l);
            self.retain_visible(rn, &mut r);

            // a LEFT breakpoint with no real RIGHT partner in window -> RIGHT-role cluster
            for lb in &l {
                if r.iter().any(|rb| in_window(rb.breakpoint - lb.breakpoint)) {
                    continue;
                }
                if !self.disc_coverage_ok(rn, lb.breakpoint) {
                    continue;
                }
                // Reject a lone breakpoint whose genomic *flank* (aligned side) is a
                // low-complexity satellite array: a real insertion has a unique/complex
                // flank, whereas a pericentromeric/subtelomeric mismap is satellite on
                // both sides. (The element clip itself may be legitimately low-complexity
                // — an Alu poly-A tail or SVA VNTR — so the flank, not the clip, is gated.)
                if is_low_complexity(&lb.unclipped.seq, 0.8) || !self.disc_mates_ok(lb) {
                    continue;
                }
                if let Some(c) = self.find_disc_cluster(rn, CLIP_RIGHT, lb.breakpoint + tsd_min, lb.breakpoint + span) {
                    print_output(writer, Emit::Bp(lb), Emit::Disc(c))?;
                    paired += 1;
                }
            }
            // a RIGHT breakpoint with no real LEFT partner in window -> LEFT-role cluster
            for rb in &r {
                if l.iter().any(|lb| in_window(rb.breakpoint - lb.breakpoint)) {
                    continue;
                }
                if !self.disc_coverage_ok(rn, rb.breakpoint) {
                    continue;
                }
                if is_low_complexity(&rb.unclipped.seq, 0.8) || !self.disc_mates_ok(rb) {
                    continue;
                }
                if let Some(c) = self.find_disc_cluster(rn, CLIP_LEFT, rb.breakpoint - span, rb.breakpoint - tsd_min) {
                    print_output(writer, Emit::Disc(c), Emit::Bp(rb))?;
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

// --- SENS-5 hallmark features (non-gating annotation) ---

fn emit_clip<'a>(e: &'a Emit) -> Option<&'a QualitySeq> {
    match e {
        Emit::Bp(b) => Some(&b.clipped),
        Emit::Pa(p) => p.clipped.as_ref(),
        Emit::Disc(_) => None,
    }
}

fn emit_unclip<'a>(e: &'a Emit) -> Option<&'a QualitySeq> {
    match e {
        Emit::Bp(b) => Some(&b.unclipped),
        Emit::Pa(_) | Emit::Disc(_) => None,
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

fn print_output<W: Write>(writer: &mut W, left: Emit, right: Emit) -> io::Result<()> {
    let (reference_name, left_str, right_str) = match (&left, &right) {
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
        _ => unreachable!("invalid pairing: polyA/disc ends are only paired with a real breakpoint"),
    };
    let bp_name = format!("{}:{}-{}", reference_name, left_str, right_str);

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
        Emit::Disc(_) => {}
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
        Emit::Disc(_) => {}
    }
    Ok(())
}
