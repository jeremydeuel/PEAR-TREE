//! Discovery-relevant configuration, mirroring the `discovery` and `adapters`
//! sections of src/config.py. (Stage 2 will make these load from a file.)

use rustc_hash::FxHashSet;

pub const MIN_MAPQ: u8 = 40;
pub const MIN_CLIP_LEN: usize = 12;
pub const MIN_EVIDENCE_READS_PER_BREAKPOINT: usize = 2;
/// Feature A both-sided rescue: minimum clip support on the SOLID anchor junction (the
/// side that clustered normally). Defaults to the global evidence floor (2); can be raised
/// to demand a stronger anchor.
pub const DISCORDANT_MIN_ANCHOR_READS: usize = 2;
/// Feature A both-sided rescue: minimum clip support on the MISSING junction, i.e. the
/// reciprocal clip that the normal >=2 floor would drop. 1 = a single (possibly rescued
/// low-MAPQ) clip pinning the junction is enough, since it is corroborated by both the
/// TSD-window pairing with the solid anchor and (for low-MAPQ clips) a uniquely-mapped mate.
pub const DISCORDANT_PARTNER_MIN_READS: usize = 1;
pub const MIN_ADAPTERLEN_FOR_CLIP: usize = 4;
pub const MIN_GOOD_BASES: usize = 10;
pub const EXCLUDE_SAME_CONTIG_SUPPLEMENTARY: i64 = 1000;

pub const POLYA_CUTOFF: usize = 12;

// Clustering + TSD-pairing windows used in cleanup()/output(). Extracted verbatim
// from discovery.py (OBS-2); values unchanged, so output stays byte-identical.
// OBS-3 will revisit CLUSTER_WINDOW vs TSD_MAX jointly with SENS-1.
//
// single-linkage cluster window: two breakpoints closer than this join one group
pub const CLUSTER_WINDOW: i64 = 6;
// TSD size bounds when pairing a left with a right breakpoint (right - left):
// below TSD_MIN it is not a TSD; above TSD_MAX the gap is too large to pair.
pub const TSD_MIN: i64 = 2;
pub const TSD_MAX: i64 = 40;
// polyA-rescue proximity to the anchor breakpoint: a polyA within POLYA_NEAR_DIST
// is too close (skipped); it rescues only when strictly inside POLYA_FAR_DIST.
pub const POLYA_NEAR_DIST: i64 = 12;
pub const POLYA_FAR_DIST: i64 = 120;

// drop clipped reads whose XA/SA shows the whole read maps contiguously
// elsewhere (not a real junction). On by default; set PEARTREE_KEEP_FULLMAP=1
// to disable (e.g. to reproduce pre-filter output).
pub fn reject_fully_mapping_reads() -> bool {
    !matches!(std::env::var("PEARTREE_KEEP_FULLMAP").as_deref(), Ok("1") | Ok("true") | Ok("True"))
}

// clip side constants (match the Python ints)
pub const CLIP_RIGHT: i32 = 1;
pub const CLIP_LEFT: i32 = 2;

pub const ADAPTERS: [&[u8]; 4] = [
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    b"AGATCGGAAAGCACACGTCTGAACTCCAGTCA",
    b"AGATCGGAAAGCGTCGTGTAGGGAAAGAGTGT",
];

/// Runtime discovery configuration (OBS-4). Mirrors the tunable parameters of the
/// `discovery` block in src/config.py. Every field defaults to the module constant
/// above, so a run with no `--config` file and no env override is byte-identical to
/// the pre-OBS-4 port.
///
/// ⚠ The default `min_mapq` is 40, matching the *generic* src/config.py. Production
/// configs `config_hs.py` / `config_mm.py` use **60**. Pin `min_mapq` (via `--config`
/// or `PEARTREE_MIN_MAPQ`) on any real-WGS run, or the port is only faithful to the
/// generic config.
#[derive(Clone, Debug)]
pub struct DiscoveryConfig {
    pub min_mapq: u8,
    pub min_evidence_reads_per_breakpoint: usize,
    pub min_good_bases: usize,
    pub exclude_same_contig_supplementary: i64,
    pub cluster_window: i64,
    pub tsd_min: i64,
    pub tsd_max: i64,
    pub polya_near_dist: i64,
    pub polya_far_dist: i64,
    pub reject_fully_mapping_reads: bool,
    /// Mate-anchored rescue: keep a soft-clipped read below `min_mapq` if it is a proper
    /// pair whose mate maps uniquely (mate MAPQ `MQ` tag >= `min_mapq`). The unique mate
    /// anchors the breakpoint, recovering insertions into low-mapability-but-mate-unique
    /// flanks without lowering the global MAPQ floor (centromeric artefacts have BOTH
    /// mates low-MAPQ). Needs the `MQ` tag from `samtools fixmate`. Off by default (a
    /// behaviour change; enable via `--config` / `PEARTREE_MATE_RESCUE` after a real-WGS
    /// check). Byte-identical to Python `src/discovery.py`.
    pub mate_anchor_rescue: bool,
    /// SPEC-5/SENS-4: explicit primary-assembly allowlist. When `Some`, a contig is
    /// processed iff its name is in the set (replacing the `len(name) > 5` + MT/chrM
    /// heuristic). When `None`, the legacy heuristic applies, so output is unchanged.
    pub contig_allowlist: Option<FxHashSet<String>>,
    /// SPEC-5: path to a BED file of regions to drop breakpoints in. `None` = no-op.
    pub exclude_bed: Option<String>,
    /// SPEC-3: drop breakpoints whose local coverage exceeds `coverage_mask_multiplier`
    /// times the genome-wide median (pileup mask). Off by default.
    pub coverage_mask: bool,
    pub coverage_mask_multiplier: f64,
    /// SPEC-4: scale the evidence floor by local/median coverage (never below the
    /// base `min_evidence_reads_per_breakpoint`). Off by default.
    pub adaptive_evidence: bool,
    /// shared coverage-estimator parameters (SPEC-3/4)
    pub coverage_bin_size: i64,
    pub coverage_sample_size: usize,
    /// SENS-1/OBS-3: count evidence reads within +/- this jitter of the modal
    /// breakpoint, not only at the exact mode. 0 = exact (legacy). Tie to the
    /// cluster window (OBS-3). Gate with SPEC-3/4 before defaulting on.
    pub evidence_window: i64,
    /// SENS-7: extend consensus while the best base strictly beats the second-best
    /// (rather than beating the sum of all others). Off = legacy.
    pub consensus_tolerant: bool,
    /// SENS-2: per-locus guard shipped with a lowered `min_mapq`. Drop a breakpoint
    /// when the fraction of its supporting clipped reads with MAPQ below
    /// `lowq_mapq_threshold` exceeds this. `None` = off (no guard).
    pub max_lowq_clip_ratio: Option<f64>,
    pub lowq_mapq_threshold: u8,
    /// SENS-5: write a per-insertion hallmark annotation to `<out>.hallmarks.tsv`
    /// (poly-A purity, TSD length, EN motif). NON-GATING — the main output is
    /// unchanged, so this is byte-identical either way. Off by default.
    pub hallmark_score: bool,
    /// SENS-8: allow clips down to `short_polya_min_clip` (instead of min_clip_len)
    /// when the clipped consensus is a pure poly-A/T terminus. Off by default.
    pub short_polya_clip: bool,
    pub short_polya_min_clip: usize,
    /// SPEC-7: drop breakpoints inside a *young* RepeatMasker element (percent
    /// divergence <= `rm_divergence_max`) from `rm_track` — divergence-gated, not
    /// family membership. `rm_track` must match the BAM's assembly (hs1 tracks are
    /// in-repo; supply a GRCh38 / GRCm39 track for those). Off by default.
    pub rm_self_mask: bool,
    pub rm_track: Option<String>,
    pub rm_divergence_max: f64,
    /// SPEC-8: drop breakpoints minted by homopolymer slippage — the aligned side of the
    /// junction ends in a run of base X >= `slippage_min_ref_run`, and the clip is
    /// >= `slippage_min_clip_frac` base X. That is *more of a tract already in the
    /// reference*, not an insertion. Paired by design: a real MEI's 3' poly-A clip is
    /// spared because the reference at its junction is ordinary sequence. Needs no
    /// reference file — the aligned part of the read is the reference. Off by default.
    ///
    /// Measured 2026-07-16 on real PD44579 discovery output: drops 41% of the known-artefact
    /// 3-5 carrier band and 30% of the contract, costing **0.6% of real germline MEIs**
    /// (7/1144, all rare, max AF 0.316; 0/40 common MEIs touched). `min_ref_run` is nearly
    /// inert (6 vs 12 -> 42% vs 39%); `min_clip_frac` is the whole lever (0.6 -> 41%,
    /// 0.9 -> 9%). Enrichment is only ~1.3-1.4x in every mode: contract-wide attrition at
    /// small recall cost, NOT an FP-specific classifier. It cannot replace min_dispersion.
    ///
    /// ⚠ `min_clip_frac = 0` degrades this to a reference-only "junction is in a tandem
    /// tract" test. Tempting (band 80% vs 41%) but it costs **5x the real recall**: 3.0%
    /// of germline MEIs (34/1144), rising to 4.5% at max_period<=6. The pairing is what
    /// makes the gate safe — keep min_clip_frac > 0.
    ///
    /// ⚠ Recall MUST be scored against `testdata/mei/1kg.sv.vcf.gz` (1000G phase-3 MEIs on
    /// hs37d5 — 1,144 contract loci are known ALU/LINE1/SVA sites), NOT against fp10k, whose
    /// simulated implants understate A/T-context recall cost ~10x (it scored the
    /// reference-only variant at 0.30% vs the truth set's 3.0%), and NOT against the
    /// min_dispersion survivors, which are unverified. Use fp10k for mechanism + regression.
    pub slippage_filter: bool,
    pub slippage_min_ref_run: usize,
    pub slippage_min_clip_frac: f64,
    /// Longest tandem period the gate will recognise at the junction. 1 = poly-A/poly-T
    /// homopolymers only; ~6 also catches (CA)n / (TG)n / (TAAAA)n microsatellites, which
    /// slip by the same mechanism. Periods above ~6 (SVA VNTRs, 19-48bp GC-rich VNTRs) are
    /// deliberately out of scope — see the sweep in the SPEC-8 notes.
    pub slippage_max_period: usize,
    /// SPEC-8b: NEW clip-level slippage gate, separate from `slippage_filter` (which inspects
    /// the *aligned* side, and from `max_homopolymer_len` which inspects the mapped part
    /// adjacent to the junction). This one inspects the two breakpoint CLIP consensuses of a
    /// paired insertion and rejects it only when BOTH clips are homopolymer/low-complexity
    /// poly-A/T — the double-sided signature of bwa soft-clipping a reference poly-A/T tract.
    /// A real MEI has a poly-A/T tail on ONE clip and a structured element body on the other,
    /// so one-sidedness spares it (the entropy term is the per-clip carve-out for a body). Only
    /// gates real breakpoint pairs (Bp+Bp); poly-A-mate and discordant ends carry no second
    /// clip consensus and are untouched. **ON by default** (shipped 2026-07 at min_run=11,
    /// max_entropy=1.95, any_base=true); set `clip_slippage_filter = false` to disable (the
    /// config.discovery.grch38.noslip A/B arm and the unvalidated grch37 config do this).
    ///
    /// Re-validated on analysis/mei9x10 against the phylogeny-breaking FP anchor (884 germline
    /// TP / 4645 tree-break FP). At min_run=11, max_entropy=1.95, clip_slippage_any_base=true it
    /// captures 64.3% of phylogeny-breaking FP at 0.45% germline-TP loss (4 TP). The any-base
    /// homopolymer detector (vs the old A/T-only) adds +243 poly-C/poly-G FP over A/T at ZERO
    /// extra TP loss. A junction-microsatellite term was evaluated and REJECTED: it costs more
    /// germline TP than the FP it adds at every threshold (real MEIs carry tandem repeats in
    /// their poly-A/TSD region on both clips). An HMM edge-Alu rescue was likewise dominated by
    /// simply raising min_run. The <1%-TP-loss frontier is min_run 10-11 / max_entropy 1.95.
    pub clip_slippage_filter: bool,
    pub clip_slippage_min_run: usize,
    pub clip_slippage_max_entropy: f64,
    pub clip_slippage_require_same_base: bool,
    /// Count a homopolymer run of ANY base (poly-A/C/G/T) toward the slippage gate, not just
    /// A/T. On by default — captures poly-C/poly-G junction slippage the A/T-only detector is
    /// blind to (+243 FP at 0 extra germline-TP loss on mei9x10). Set false to restore the
    /// legacy A/T-only behaviour for A/B comparison.
    pub clip_slippage_any_base: bool,
    /// SPD-3: process contigs in parallel across this many worker threads for the extract
    /// pass (BAM only, requires a `.bai`/`.csi` index). 1 = the validated single-threaded
    /// path. The mate pass stays a single linear scan (N-way BGZF decode for BAM).
    pub contig_threads: usize,
    // --- Feature A: discordant-mate anchoring (all off/neutral by default) ---
    /// Master switch: collect discordant read pairs as an evidence source so a
    /// one-sided junction (e.g. a lone poly-A clip) can be paired with a cluster of
    /// discordant mates supplying the missing reciprocal side. Off = no change.
    pub discordant_anchor: bool,
    /// Same-contig template length beyond which a mapped pair counts as discordant
    /// (different-contig and non-FR pairs are always discordant).
    pub discordant_max_tlen: i64,
    /// Minimum distinct discordant pairs to form an anchoring cluster.
    pub discordant_min_reads: usize,
    /// Cluster width for discordant observations (defaults to `cluster_window`).
    pub discordant_window: i64,
    /// RTE track (RepeatMasker `.out[.gz]`) for the mate-origin check. `None` = the
    /// origin fraction is reported as 0 (label only, never blocks).
    pub discordant_rte_track: Option<String>,
    /// Keep RTE copies with percent divergence <= this for the mate-origin test.
    pub discordant_rte_divergence_max: f64,
    /// Gate: require the RTE-origin fraction >= `discordant_rte_min` before a
    /// discordant cluster may act as a partner. Off = label only.
    pub discordant_rte_only: bool,
    pub discordant_rte_min: f64,
    /// Track-free RTE-origin proxy: only collect a discordant observation whose mate
    /// maps *ambiguously* (mate MAPQ `MQ` <= this). A mate that originates in an
    /// inserted young RTE maps to many reference paralogs → low MAPQ; a mate placed
    /// uniquely (high MAPQ) reflects structural/artefactual discordance, not an RTE.
    /// `None` = collect regardless of mate MAPQ (legacy). Needs the `MQ` tag.
    pub discordant_mate_max_mapq: Option<u8>,
    /// Half-width of the search window (bp) from a lone real breakpoint to a discordant
    /// cluster on the missing side, during the output rescue. Discordant anchors sit up
    /// to ~a fragment length from the junction, so the TSD bound (`tsd_max`, ~40 bp) is
    /// far too tight. `None` = use `tsd_max` (legacy behaviour).
    pub discordant_rescue_span: Option<i64>,
    /// Reject a one-sided discordant call whose mate reads are low-diversity (a satellite
    /// array): require the mean distinct-4-mer fraction of the real breakpoint's mate reads
    /// to be >= this. A real MEI's flank-anchored mates are complex genomic/element sequence
    /// (>= ~0.45); pericentromeric/subtelomeric satellite mates fall to ~0.3. Gates the
    /// mates, not the clip, so real poly-A/VNTR element clips are spared. `None` = off.
    pub discordant_mate_min_kmer_div: Option<f64>,
    /// Stricter local-coverage ceiling for a *one-sided* discordant call: reject it if
    /// the real breakpoint's local depth exceeds this multiple of the genome median.
    /// One-sided calls are weaker evidence than a reciprocal breakpoint pair, so they
    /// warrant a tighter pileup gate than the global `coverage_mask_multiplier` (organic
    /// assembly-discordance FPs cluster in ~3-5x pileups that pass the 5x mask). `None` =
    /// fall back to `coverage_mask_multiplier`. Populates coverage even if the mask is off.
    pub discordant_coverage_max_mult: Option<f64>,
    /// Feature A both-sided rescue: minimum clip support required on the SOLID anchor
    /// junction. Both junctions of a discordant-rescued insertion must carry direct clip
    /// evidence — a discordant mate cluster alone never completes a call. This is the floor
    /// for the anchor side (defaults to the global evidence floor of 2).
    pub discordant_min_anchor_reads: usize,
    /// Feature A both-sided rescue: minimum clip support on the MISSING junction — the
    /// reciprocal clip the normal >=2 floor would drop. 1 admits a single (possibly rescued
    /// low-MAPQ) clip that pins the junction within the TSD window opposite the anchor.
    pub discordant_partner_min_reads: usize,
    /// Pin the genome-wide coverage median to this value instead of estimating it from the
    /// BAM (`None` = estimate as usual). Set per-BAM (via `PEARTREE_COVERAGE_MEDIAN` or the
    /// config key) when running on a *region slice* of a BAM: the local bin counts near the
    /// truth loci are still exact, but the genome median would be inflated by the slice, so
    /// the `local/median` ratio gates (SPEC-3, `discordant_coverage_max_mult`) must be given
    /// the true full-BAM median. Emit it with `--step coverage-median`. No effect unless a
    /// coverage-ratio gate is on.
    pub coverage_median_override: Option<f64>,
    // --- Feature B: processed-pseudogene (splice) annotation ---
    /// Write a non-gating `<out>.splice.tsv` flagging candidates whose mate reads span
    /// >= `splice_min_exons` exons of a single reference gene (intron skipped). The
    /// main breakpoint output is unchanged. Off by default.
    pub splice_hallmark: bool,
    /// Exon annotation (BED/GTF with gene_id) for the BAM's assembly. Required when
    /// `splice_hallmark` is on; must match the assembly (external input, not shipped).
    pub exon_annotation: Option<String>,
    /// Distinct same-gene exons a candidate's mates must hit to be flagged.
    pub splice_min_exons: usize,
    // --- TPRT-hallmark overhaul (plans/tprt_hallmarks/SPEC.md); all default-off ---
    /// Write `<out>.evidence.tsv.gz`: one row per read record (CLIP/POLYA/DISC/MATE) of
    /// every emitted locus. Non-gating: the `.txt.gz` output is unchanged by this key.
    pub evidence_sidecar: bool,
    /// Sidecar only: request the mate of EVERY evidence read (CLIP incl. supplementary,
    /// DISC), not only the `has_mate` orientation (LEFT reverse / RIGHT forward) that feeds
    /// the FASTQ :MATE records. The extra mates go to the sidecar only, so the FASTQ is
    /// unchanged. No effect unless `evidence_sidecar` is on.
    pub fetch_all_mates: bool,
    /// Keep 0x400 (duplicate-flagged) reads everywhere in discovery: clip path, mate pass,
    /// sidecar capture. Independence is decided downstream from the reads, never from the
    /// flag. Off = legacy: dups dropped in the clip path and mate pass, while the low-MAPQ
    /// poly-A path never checked the flag (kept as-is for byte-identity).
    pub ignore_dup_flag: bool,
    /// When set, the per-breakpoint evidence floor counts DISTINCT FRAGMENTS (qname hash:
    /// both mates and primary+supplementary of one template count once) instead of reads,
    /// and a single-fragment cluster goes through the normal consensus path when the floor
    /// is <= 1. `adaptive_evidence` scales this base like the read floor. `None` = the
    /// legacy read floor `min_evidence_reads_per_breakpoint`.
    pub min_evidence_fragments_per_sample: Option<usize>,
    /// Sidecar caps per breakpoint side; deterministic subset (lowest fragment hash).
    pub max_mates_per_breakpoint: usize,
    pub max_evidence_reads_per_breakpoint: usize,
    /// Sidecar DISC rows: a non-proper high-MAPQ anchor (mate unmapped / other contig /
    /// far / same strand) pointing at the junction from within this many bp.
    pub sidecar_disc_span: i64,
    // --- TPRT pairing modes (all off by default => output byte-identical) ---
    /// Drop 0x400 duplicates also in the low-MAPQ poly-A path (`find_polya`), which
    /// historically never checked the flag. Off = legacy (byte-identical).
    pub drop_dup_in_polya_path: bool,
    /// Target-site deletion: after the normal TSD pairing, pair a leftover LEFT/RIGHT
    /// breakpoint whose RIGHT lies up to this many bp LEFT of the LEFT (negative TSD).
    /// 0 = off.
    pub max_target_site_deletion: i64,
    /// Blunt / EN-independent joins: pair leftover breakpoints with gap 0..tsd_min-1.
    pub allow_blunt_pairs: bool,
    /// L1-mediated deletion / duplication: pair leftover breakpoints up to this many bp
    /// apart (gap < -max_target_site_deletion, or gap > tsd_max) when the pair is
    /// TPRT-polarised (one clip a poly-A/T tail, the other a complex element clip).
    /// 0 = off.
    pub max_l1_mediated_span: i64,
    /// Poly-T run (stored, junction-outward orientation) that marks the poly-A side of a
    /// polarised L1-mediated pair.
    pub l1_mediated_min_polya: usize,
    /// One-sided loci: emit a still-unpaired breakpoint as `contig:L-oneside_L` /
    /// `contig:oneside_R-R` when it has >= `one_sided_min_fragments` fragments and passes
    /// the SPEC-8 reference-tract slippage test (always applied here).
    pub one_sided_loci: bool,
    pub one_sided_min_fragments: usize,
    /// Only poly-A/T-tail junctions (stored clip starts with >= `one_sided_min_polya` T).
    pub one_sided_require_polya: bool,
    pub one_sided_min_polya: usize,
    /// SPEC-8b junction spare: a Bp+Bp pair is never SPEC-8b-rejected when either stored
    /// clip's first k junction-proximal bases are structured (no homopolymer >= 8, >= 3
    /// distinct bases). Rescues short-tag orphan transductions / short inserts whose long
    /// poly-A drags the whole-clip entropy under the threshold. 0 = off (legacy).
    pub clip_slippage_junction_spare: usize,
    /// Sidecar SHORT rows: primary, non-dup, MAPQ >= min_mapq reads crossing an emitted
    /// junction by only a few bases — (a) a junction-side soft clip of 1..MIN_CLIP_LEN-1
    /// bases starting within ±`short_overhang_window` of the breakpoint, or (b) aligned
    /// through it by 1..`short_overhang_max` bases. Collection only (combine decides).
    /// Needs `evidence_sidecar`. Off = sidecar byte-identical.
    pub short_overhang_evidence: bool,
    pub short_overhang_window: i64,
    pub short_overhang_max: i64,
    pub max_short_per_breakpoint: usize,
    /// How the sidecar pass fetches its records: false = one linear BAM scan (default),
    /// true = index seeks to the requested records' known positions (own POS, or the
    /// requesting read's RNEXT/PNEXT for mates), forward-only, so the records are met
    /// in file order and the sidecar is byte-identical. Falls back to the scan for CRAM,
    /// an unindexed BAM, or any request without a position.
    pub sidecar_indexed_fetch: bool,
    /// Gap (bp) below which neighbouring fetch windows are read through, not re-seeked.
    pub sidecar_fetch_gap: i64,
}

impl Default for DiscoveryConfig {
    fn default() -> Self {
        DiscoveryConfig {
            min_mapq: MIN_MAPQ,
            min_evidence_reads_per_breakpoint: MIN_EVIDENCE_READS_PER_BREAKPOINT,
            min_good_bases: MIN_GOOD_BASES,
            exclude_same_contig_supplementary: EXCLUDE_SAME_CONTIG_SUPPLEMENTARY,
            cluster_window: CLUSTER_WINDOW,
            tsd_min: TSD_MIN,
            tsd_max: TSD_MAX,
            polya_near_dist: POLYA_NEAR_DIST,
            polya_far_dist: POLYA_FAR_DIST,
            reject_fully_mapping_reads: reject_fully_mapping_reads(),
            mate_anchor_rescue: false,
            contig_allowlist: None,
            exclude_bed: None,
            coverage_mask: false,
            coverage_mask_multiplier: 5.0,
            adaptive_evidence: false,
            coverage_bin_size: 500,
            coverage_sample_size: 3000,
            evidence_window: 0,
            consensus_tolerant: false,
            max_lowq_clip_ratio: None,
            lowq_mapq_threshold: 40,
            hallmark_score: false,
            short_polya_clip: false,
            short_polya_min_clip: 7,
            rm_self_mask: false,
            rm_track: None,
            rm_divergence_max: 5.0,
            slippage_filter: false,
            slippage_min_ref_run: 8,
            slippage_min_clip_frac: 0.6,
            slippage_max_period: 1,
            clip_slippage_filter: true, // SHIPPED default (was false) — see doc comment
            clip_slippage_min_run: 11, // was 18 (validated operating point, analysis/mei9x10)
            clip_slippage_max_entropy: 1.95, // was 1.88
            clip_slippage_require_same_base: false,
            clip_slippage_any_base: true,
            contig_threads: 1,
            discordant_anchor: false,
            discordant_max_tlen: 1000,
            discordant_min_reads: 3,
            discordant_window: CLUSTER_WINDOW,
            discordant_rte_track: None,
            discordant_rte_divergence_max: 20.0,
            discordant_rte_only: false,
            discordant_rte_min: 0.5,
            discordant_mate_max_mapq: None,
            discordant_rescue_span: None,
            discordant_mate_min_kmer_div: None,
            discordant_coverage_max_mult: None,
            discordant_min_anchor_reads: DISCORDANT_MIN_ANCHOR_READS,
            discordant_partner_min_reads: DISCORDANT_PARTNER_MIN_READS,
            coverage_median_override: None,
            splice_hallmark: false,
            exon_annotation: None,
            splice_min_exons: 2,
            evidence_sidecar: false,
            fetch_all_mates: false,
            ignore_dup_flag: false,
            min_evidence_fragments_per_sample: None,
            max_mates_per_breakpoint: 50,
            max_evidence_reads_per_breakpoint: 200,
            sidecar_disc_span: 500,
            drop_dup_in_polya_path: false,
            max_target_site_deletion: 0,
            allow_blunt_pairs: false,
            max_l1_mediated_span: 0,
            l1_mediated_min_polya: 10,
            one_sided_loci: false,
            one_sided_min_fragments: 2,
            one_sided_require_polya: true,
            one_sided_min_polya: 10,
            clip_slippage_junction_spare: 0,
            short_overhang_evidence: false,
            short_overhang_window: 3,
            short_overhang_max: 20,
            max_short_per_breakpoint: 100,
            sidecar_indexed_fetch: false,
            sidecar_fetch_gap: 16384,
        }
    }
}

fn parse_num<T: std::str::FromStr>(v: &str) -> Result<T, String> {
    v.parse().map_err(|_| format!("not a valid number: '{v}'"))
}

fn parse_bool(v: &str) -> Result<bool, String> {
    match v {
        "true" | "True" | "1" => Ok(true),
        "false" | "False" | "0" => Ok(false),
        _ => Err(format!("not a valid bool: '{v}'")),
    }
}

impl DiscoveryConfig {
    /// True when any post-TSD pairing mode (deletion / blunt / L1-mediated / one-sided)
    /// is enabled. Off = the legacy pairing only.
    pub fn extra_pairing(&self) -> bool {
        self.max_target_site_deletion > 0 || self.allow_blunt_pairs || self.max_l1_mediated_span > 0 || self.one_sided_loci
    }

    /// Build config: start from defaults, overlay a `key = value` file (if given),
    /// then overlay environment overrides (which always win). `#` starts a comment.
    pub fn load(path: Option<&str>) -> Result<DiscoveryConfig, String> {
        let mut cfg = DiscoveryConfig::default();
        if let Some(p) = path {
            let text = std::fs::read_to_string(p).map_err(|e| format!("cannot read config {p}: {e}"))?;
            for (i, raw) in text.lines().enumerate() {
                let line = raw.split('#').next().unwrap_or("").trim();
                if line.is_empty() {
                    continue;
                }
                let (key, val) = line
                    .split_once('=')
                    .ok_or_else(|| format!("{p}:{}: expected 'key = value'", i + 1))?;
                cfg.set(key.trim(), val.trim()).map_err(|e| format!("{p}:{}: {e}", i + 1))?;
            }
        }
        cfg.apply_env()?;
        Ok(cfg)
    }

    fn set(&mut self, key: &str, val: &str) -> Result<(), String> {
        match key {
            "min_mapq" => self.min_mapq = parse_num(val)?,
            "min_evidence_reads_per_breakpoint" => self.min_evidence_reads_per_breakpoint = parse_num(val)?,
            "min_good_bases" => self.min_good_bases = parse_num(val)?,
            "exclude_same_contig_supplementary" => self.exclude_same_contig_supplementary = parse_num(val)?,
            "cluster_window" => self.cluster_window = parse_num(val)?,
            "tsd_min" => self.tsd_min = parse_num(val)?,
            // accept the Python name `max_bp_window` as an alias for tsd_max
            "tsd_max" | "max_bp_window" => self.tsd_max = parse_num(val)?,
            "polya_near_dist" => self.polya_near_dist = parse_num(val)?,
            "polya_far_dist" => self.polya_far_dist = parse_num(val)?,
            "reject_fully_mapping_reads" => self.reject_fully_mapping_reads = parse_bool(val)?,
            "mate_anchor_rescue" => self.mate_anchor_rescue = parse_bool(val)?,
            // comma-separated inline allowlist, e.g. "1,2,...,X,Y,MT"
            "contig_allowlist" => {
                self.contig_allowlist = Some(val.split(',').map(|s| s.trim()).filter(|s| !s.is_empty()).map(String::from).collect())
            }
            // allowlist from a file, one contig name per line (# comments allowed)
            "contig_allowlist_file" => {
                let text = std::fs::read_to_string(val).map_err(|e| format!("cannot read contig_allowlist_file {val}: {e}"))?;
                self.contig_allowlist = Some(
                    text.lines()
                        .map(|l| l.split('#').next().unwrap_or("").trim())
                        .filter(|s| !s.is_empty())
                        .map(String::from)
                        .collect(),
                )
            }
            "exclude_bed" => self.exclude_bed = Some(val.to_string()),
            "coverage_mask" => self.coverage_mask = parse_bool(val)?,
            "coverage_mask_multiplier" => self.coverage_mask_multiplier = parse_num(val)?,
            "adaptive_evidence" => self.adaptive_evidence = parse_bool(val)?,
            "coverage_bin_size" => self.coverage_bin_size = parse_num(val)?,
            "coverage_sample_size" => self.coverage_sample_size = parse_num(val)?,
            "evidence_window" => self.evidence_window = parse_num(val)?,
            "consensus_tolerant" => self.consensus_tolerant = parse_bool(val)?,
            "max_lowq_clip_ratio" => self.max_lowq_clip_ratio = Some(parse_num(val)?),
            "lowq_mapq_threshold" => self.lowq_mapq_threshold = parse_num(val)?,
            "hallmark_score" => self.hallmark_score = parse_bool(val)?,
            "short_polya_clip" => self.short_polya_clip = parse_bool(val)?,
            "short_polya_min_clip" => self.short_polya_min_clip = parse_num(val)?,
            "rm_self_mask" => self.rm_self_mask = parse_bool(val)?,
            "rm_track" => self.rm_track = Some(val.to_string()),
            "rm_divergence_max" => self.rm_divergence_max = parse_num(val)?,
            "slippage_filter" => self.slippage_filter = parse_bool(val)?,
            "slippage_min_ref_run" => self.slippage_min_ref_run = parse_num(val)?,
            "slippage_min_clip_frac" => self.slippage_min_clip_frac = parse_num(val)?,
            "slippage_max_period" => self.slippage_max_period = parse_num(val)?,
            "clip_slippage_filter" => self.clip_slippage_filter = parse_bool(val)?,
            "clip_slippage_min_run" => self.clip_slippage_min_run = parse_num(val)?,
            "clip_slippage_max_entropy" => self.clip_slippage_max_entropy = parse_num(val)?,
            "clip_slippage_require_same_base" => self.clip_slippage_require_same_base = parse_bool(val)?,
            "clip_slippage_any_base" => self.clip_slippage_any_base = parse_bool(val)?,
            "contig_threads" => self.contig_threads = parse_num(val)?,
            "discordant_anchor" => self.discordant_anchor = parse_bool(val)?,
            "discordant_max_tlen" => self.discordant_max_tlen = parse_num(val)?,
            "discordant_min_reads" => self.discordant_min_reads = parse_num(val)?,
            "discordant_window" => self.discordant_window = parse_num(val)?,
            "discordant_rte_track" => self.discordant_rte_track = Some(val.to_string()),
            "discordant_rte_divergence_max" => self.discordant_rte_divergence_max = parse_num(val)?,
            "discordant_rte_only" => self.discordant_rte_only = parse_bool(val)?,
            "discordant_rte_min" => self.discordant_rte_min = parse_num(val)?,
            "discordant_mate_max_mapq" => self.discordant_mate_max_mapq = Some(parse_num(val)?),
            "discordant_rescue_span" => self.discordant_rescue_span = Some(parse_num(val)?),
            "discordant_mate_min_kmer_div" => self.discordant_mate_min_kmer_div = Some(parse_num(val)?),
            "discordant_coverage_max_mult" => self.discordant_coverage_max_mult = Some(parse_num(val)?),
            "discordant_min_anchor_reads" => self.discordant_min_anchor_reads = parse_num(val)?,
            "discordant_partner_min_reads" => self.discordant_partner_min_reads = parse_num(val)?,
            "coverage_median_override" => self.coverage_median_override = Some(parse_num(val)?),
            "splice_hallmark" => self.splice_hallmark = parse_bool(val)?,
            "exon_annotation" => self.exon_annotation = Some(val.to_string()),
            "splice_min_exons" => self.splice_min_exons = parse_num(val)?,
            "evidence_sidecar" => self.evidence_sidecar = parse_bool(val)?,
            "fetch_all_mates" => self.fetch_all_mates = parse_bool(val)?,
            "ignore_dup_flag" => self.ignore_dup_flag = parse_bool(val)?,
            "min_evidence_fragments_per_sample" => self.min_evidence_fragments_per_sample = Some(parse_num(val)?),
            "max_mates_per_breakpoint" => self.max_mates_per_breakpoint = parse_num(val)?,
            "max_evidence_reads_per_breakpoint" => self.max_evidence_reads_per_breakpoint = parse_num(val)?,
            "sidecar_disc_span" => self.sidecar_disc_span = parse_num(val)?,
            "drop_dup_in_polya_path" => self.drop_dup_in_polya_path = parse_bool(val)?,
            "max_target_site_deletion" => self.max_target_site_deletion = parse_num(val)?,
            "allow_blunt_pairs" => self.allow_blunt_pairs = parse_bool(val)?,
            "max_l1_mediated_span" => self.max_l1_mediated_span = parse_num(val)?,
            "l1_mediated_min_polya" => self.l1_mediated_min_polya = parse_num(val)?,
            "one_sided_loci" => self.one_sided_loci = parse_bool(val)?,
            "one_sided_min_fragments" => self.one_sided_min_fragments = parse_num(val)?,
            "one_sided_require_polya" => self.one_sided_require_polya = parse_bool(val)?,
            "one_sided_min_polya" => self.one_sided_min_polya = parse_num(val)?,
            "clip_slippage_junction_spare" => self.clip_slippage_junction_spare = parse_num(val)?,
            "short_overhang_evidence" => self.short_overhang_evidence = parse_bool(val)?,
            "short_overhang_window" => self.short_overhang_window = parse_num(val)?,
            "short_overhang_max" => self.short_overhang_max = parse_num(val)?,
            "max_short_per_breakpoint" => self.max_short_per_breakpoint = parse_num(val)?,
            "sidecar_indexed_fetch" => self.sidecar_indexed_fetch = parse_bool(val)?,
            "sidecar_fetch_gap" => self.sidecar_fetch_gap = parse_num(val)?,
            other => eprintln!("warning: ignoring unknown config key '{other}'"),
        }
        Ok(())
    }

    fn apply_env(&mut self) -> Result<(), String> {
        if let Ok(v) = std::env::var("PEARTREE_MIN_MAPQ") {
            self.min_mapq = v.trim().parse().map_err(|_| format!("PEARTREE_MIN_MAPQ: not a valid u8: '{v}'"))?;
        }
        // PEARTREE_COVERAGE_MEDIAN pins the genome median (per-BAM, for region-slice runs).
        if let Ok(v) = std::env::var("PEARTREE_COVERAGE_MEDIAN") {
            let m: f64 = v.trim().parse().map_err(|_| format!("PEARTREE_COVERAGE_MEDIAN: not a valid f64: '{v}'"))?;
            self.coverage_median_override = Some(m);
        }
        // PEARTREE_KEEP_FULLMAP=1 keeps full-mapping reads (env wins over the file).
        if matches!(std::env::var("PEARTREE_KEEP_FULLMAP").as_deref(), Ok("1") | Ok("true") | Ok("True")) {
            self.reject_fully_mapping_reads = false;
        }
        // PEARTREE_MATE_RESCUE toggles mate-anchored rescue (env wins over the file).
        match std::env::var("PEARTREE_MATE_RESCUE").as_deref() {
            Ok("1") | Ok("true") | Ok("True") => self.mate_anchor_rescue = true,
            Ok("0") | Ok("false") | Ok("False") => self.mate_anchor_rescue = false,
            _ => {}
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_constants() {
        let c = DiscoveryConfig::default();
        assert_eq!(c.min_mapq, MIN_MAPQ);
        assert_eq!(c.min_evidence_reads_per_breakpoint, MIN_EVIDENCE_READS_PER_BREAKPOINT);
        assert_eq!(c.cluster_window, CLUSTER_WINDOW);
        assert_eq!(c.tsd_max, TSD_MAX);
        assert_eq!(c.polya_far_dist, POLYA_FAR_DIST);
    }

    #[test]
    fn spec8b_ships_on_at_validated_operating_point() {
        // Pins the 2026-07 shipped default so a refactor can't silently revert it.
        let c = DiscoveryConfig::default();
        assert!(c.clip_slippage_filter, "SPEC-8b must be ON by default");
        assert_eq!(c.clip_slippage_min_run, 11);
        assert_eq!(c.clip_slippage_max_entropy, 1.95);
        assert!(c.clip_slippage_any_base, "any-base homopolymer must be ON");
        assert!(!c.clip_slippage_require_same_base);
    }

    #[test]
    fn set_overrides_and_alias() {
        let mut c = DiscoveryConfig::default();
        c.set("min_mapq", "60").unwrap();
        c.set("max_bp_window", "50").unwrap(); // Python-name alias for tsd_max
        c.set("reject_fully_mapping_reads", "false").unwrap();
        assert_eq!(c.min_mapq, 60);
        assert_eq!(c.tsd_max, 50);
        assert!(!c.reject_fully_mapping_reads);
    }

    #[test]
    fn unknown_key_is_ignored_not_errored() {
        let mut c = DiscoveryConfig::default();
        assert!(c.set("no_such_key", "3").is_ok());
    }

    #[test]
    fn bad_value_errors() {
        let mut c = DiscoveryConfig::default();
        assert!(c.set("min_mapq", "notnum").is_err());
        assert!(c.set("reject_fully_mapping_reads", "maybe").is_err());
    }

    #[test]
    fn feature_ab_defaults_are_off() {
        let c = DiscoveryConfig::default();
        assert!(!c.discordant_anchor);
        assert_eq!(c.discordant_max_tlen, 1000);
        assert_eq!(c.discordant_min_reads, 3);
        assert_eq!(c.discordant_min_anchor_reads, DISCORDANT_MIN_ANCHOR_READS);
        assert_eq!(c.discordant_partner_min_reads, DISCORDANT_PARTNER_MIN_READS);
        assert_eq!(c.discordant_window, CLUSTER_WINDOW);
        assert!(c.discordant_rte_track.is_none());
        assert!(!c.discordant_rte_only);
        assert!(!c.splice_hallmark);
        assert!(c.exon_annotation.is_none());
        assert_eq!(c.splice_min_exons, 2);
    }

    #[test]
    fn feature_ab_keys_set() {
        let mut c = DiscoveryConfig::default();
        c.set("discordant_anchor", "true").unwrap();
        c.set("discordant_max_tlen", "2500").unwrap();
        c.set("discordant_min_reads", "5").unwrap();
        c.set("discordant_min_anchor_reads", "6").unwrap();
        c.set("discordant_partner_min_reads", "2").unwrap();
        c.set("discordant_rte_track", "/data/rmsk.out.gz").unwrap();
        c.set("discordant_rte_only", "1").unwrap();
        c.set("discordant_rte_min", "0.75").unwrap();
        c.set("splice_hallmark", "true").unwrap();
        c.set("exon_annotation", "/data/exons.bed").unwrap();
        c.set("splice_min_exons", "3").unwrap();
        assert!(c.discordant_anchor);
        assert_eq!(c.discordant_max_tlen, 2500);
        assert_eq!(c.discordant_min_reads, 5);
        assert_eq!(c.discordant_min_anchor_reads, 6);
        assert_eq!(c.discordant_partner_min_reads, 2);
        assert_eq!(c.discordant_rte_track.as_deref(), Some("/data/rmsk.out.gz"));
        assert!(c.discordant_rte_only);
        assert!((c.discordant_rte_min - 0.75).abs() < 1e-9);
        assert!(c.splice_hallmark);
        assert_eq!(c.exon_annotation.as_deref(), Some("/data/exons.bed"));
        assert_eq!(c.splice_min_exons, 3);
    }

    #[test]
    fn tprt_keys_default_off_and_parse() {
        let c = DiscoveryConfig::default();
        assert!(!c.evidence_sidecar && !c.fetch_all_mates && !c.ignore_dup_flag);
        assert!(c.min_evidence_fragments_per_sample.is_none());
        assert_eq!((c.max_mates_per_breakpoint, c.max_evidence_reads_per_breakpoint), (50, 200));
        let mut c = c;
        c.set("evidence_sidecar", "true").unwrap();
        c.set("fetch_all_mates", "true").unwrap();
        c.set("ignore_dup_flag", "true").unwrap();
        c.set("min_evidence_fragments_per_sample", "1").unwrap();
        assert!(c.evidence_sidecar && c.fetch_all_mates && c.ignore_dup_flag);
        assert_eq!(c.min_evidence_fragments_per_sample, Some(1));
    }

    #[test]
    fn pairing_mode_keys_default_off_and_parse() {
        let c = DiscoveryConfig::default();
        assert!(!c.extra_pairing() && !c.drop_dup_in_polya_path);
        let mut c = c;
        c.set("max_target_site_deletion", "30").unwrap();
        c.set("allow_blunt_pairs", "true").unwrap();
        c.set("max_l1_mediated_span", "50000").unwrap();
        c.set("one_sided_loci", "true").unwrap();
        c.set("one_sided_min_fragments", "3").unwrap();
        c.set("drop_dup_in_polya_path", "true").unwrap();
        assert!(c.extra_pairing() && c.drop_dup_in_polya_path && c.allow_blunt_pairs);
        assert_eq!((c.max_target_site_deletion, c.max_l1_mediated_span, c.one_sided_min_fragments), (30, 50000, 3));
    }

    #[test]
    fn inline_allowlist_parses() {
        let mut c = DiscoveryConfig::default();
        assert!(c.contig_allowlist.is_none());
        c.set("contig_allowlist", "1, 2 ,X,, NC_000014.9").unwrap();
        let set = c.contig_allowlist.unwrap();
        assert_eq!(set.len(), 4); // empty entry between the commas is dropped
        assert!(set.contains("1"));
        assert!(set.contains("X"));
        assert!(set.contains("NC_000014.9"));
        assert!(!set.contains("3"));
    }
}
