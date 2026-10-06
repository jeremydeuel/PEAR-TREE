# peartree-genotype2 — realignment-based, phylogeny-aware SV/MEI genotyper (Rust)

Branch `genotype-v2` (off `tprt-hallmarks`). New crate `rust/peartree-genotype2`, binary
`peartree-genotype2`. The legacy `peartree-genotype` stays untouched for A/B. This file is the
**contract between parallel work packages**: every `pub` type/signature in `src/types.rs`,
`src/config.rs` and the stub modules is fixed. If you must change one, change it in the same
commit and say so in your report.

## Why (measured on the legacy genotyper)

* Scores only the 12 junction-adjacent bases per side, ±1 bp register shift, no indels
  (poly-A length jitter = mismatches). Hard votes → VAF bands → per-geometry fudge factors.
* `min_mapq 60` on reads whose MAPQ reflects only the short aligned flank → reference bias
  (het VAF floor 0.46 on 9x10, 7% deficit); the most informative alt reads are dropped.
* Each colony genotyped alone; clade members with few reads become `wild-type?` / `artefact`
  (PD37590 chr12:81091797, chr7:131206236 under-called; E2E 10-colony: 234/1610 truth-present
  colony calls are no-calls).
* Two index queries per locus (count + fetch), unbuffered reader → Lustre-latency bound
  (~35 ms/locus, 18 min/colony for 31k loci).

## Inputs

| flag | what |
|---|---|
| `--bam` | BAM/CRAM, indexed |
| `--insertions` | `<patient>.genotyping[.tprt].txt.gz` — the locus LIST (names; 12-bp consensus only as fallback) |
| `--combined` | `<patient>.combined.txt.gz` — FULL junction consensus per side (see format), optional but strongly recommended |
| `--reference` | `.fa` (+`.fai`) **or** `.2bit` of the BAM's assembly (farm: `/lustre/scratch126/casm/teams/team273/users/jd43/hg38.2bit`) — required |
| `--out`, `--config`, `--threads`, `--step genotype|genotype_batch|joint`, `--manifest` | as legacy |

### Locus name geometry (0-based coordinates L = left_pos, R = right_pos)

Alt haplotype = `genome[.., R) ++ INS ++ genome[L, ..)`.

| kind | gap = R − L | geometry |
|---|---|---|
| TSD | 2..40 | `genome[L,R)` duplicated (target-site duplication) |
| BLUNT | 0..1 | |
| TSD_DELETION | −30..−1 | `genome[R,L)` deleted |
| FAR_DELETION | < −30 (up to 50 kb) | L1-mediated deletion |
| FAR_DUPLICATION | > 40 | L1-mediated duplication; for `gap >= dup_retained_min_span` the alt haplotype still carries both reference junctions |
| ONE_SIDED | `contig:L-oneside_L` (right end missing, real junction at L) / `contig:oneside_R-R` (left end missing, real junction at R) | only one junction is known |

### `combined.txt.gz` record format (what `--combined` provides)

```
@<locus>:R        # right junction: UPPERCASE genome[R-n, R) then lowercase INS 5' part
<seq>
+
<qual phred+33>
@<locus>:L        # left junction: lowercase INS 3' part (…poly-A) then UPPERCASE genome[L, L+m)
```
Case encodes aligned (upper) vs clipped/inserted (lower). Qualities are consensus confidences.
`insR` = lowercase of `:R`; `insL` = lowercase of `:L`. Both are in read orientation on the
reference-forward strand, i.e. exactly as they appear in a read spanning that junction.
One-sided loci have only the real side's record (built by `src/genotyping_contract_oneside.py`
as `<locus>:L`/`:R` in `genotyping.tprt.txt.gz` — the combined file has the two-sided name
`contig:L-R`... **A owner: check the one-sided case in `$E2E/combine/P1.combined.txt.gz` and
`$E2E/genotype/P1.genotyping.tprt.txt.gz` and document what you found in contract.rs.**)

## Per-locus haplotype model (`haplotype.rs`, owner A)

`F = cfg.flank` (default 300 ≥ read length + 100). All segments carry an `anchor`
(reference coordinate ↔ segment index) for the genome part so the aligner can seed a band.

* `REF`: if `|gap| <= F` one segment `genome[min(L,R)-F, max(L,R)+F)`; else two, `genome[L-F,L+F)` and `genome[R-F,R+F)`.
* `ALT_R` = `genome[R-F, R) ++ insR` (anchor: index 0 = R−F). Absent if the left end is open (`oneside_R-R` has the RIGHT junction real → ALT_R present, ALT_L absent; `L-oneside_L` → ALT_L present only).
* `ALT_L` = `insL ++ genome[L, L+F)` (anchor: index |insL| = L).
* `ALT_FULL` (replaces ALT_R/ALT_L) when `insR`'s suffix overlaps `insL`'s prefix by ≥ `merge_overlap_min` bases with ≤ `merge_overlap_max_mismatch_frac` mismatches and the overlap is not a homopolymer: `genome[R-F,R) ++ merge(insR, insL) ++ genome[L, L+F)`. This is the SHORT-insertion case (solo poly-A, short 5'-truncated L1) that the legacy `double_alt_is_artefact` rule mis-calls as chimeric.
* `alt_ref_junction_fraction` = 0.5 when kind is FAR_DUPLICATION and `gap >= dup_retained_min_span`, else 0.0 (the alt haplotype of a long duplication yields reference-junction reads too).
* Fetch windows: `|gap| <= F` → one window `[min-1, max+1)` (legacy semantics); else two `[L-1, L+1)`, `[R-1, R+1)`. One-sided: one window at the real breakpoint.
* Without `--combined`, insR/insL fall back to the 12-bp contract consensus (warn once).
* Missing a needed consensus → `LocusModel::error(reason)`; the driver emits an `error` row.

## Read likelihood (`align.rs`, `readlik.rs`, owner B)

Quality-aware affine-gap local alignment of the read against each segment (read may be
soft-clipped at either end at `ln(clip_prob)` per base, default ln(1/4); segment ends free).
Scores are natural-log likelihoods:

* match `ln(1-e)`, mismatch `ln(e/3)`, `e = min(0.75, e_read + e_cons)`, `e_x = 10^(-q/10)`,
  q clamped to `[base_q_min, base_q_max]`; `N` in read or segment: `ln(1/4)`.
* gap open `ln(10^(-gap_open_phred/10))`, extend `ln(10^(-gap_ext_phred/10))`; inside a
  haplotype homopolymer run ≥ `homopolymer_min_len` the open cost is `homopolymer_gap_open_phred`
  (poly-A length jitter is cheap, as in STR-aware pair-HMMs).
* Banded: diagonal from the read's BAM placement (`ref_start - leading_softclip` → segment
  index via the anchor), half-width `band_halfwidth`; if the banded score is below
  `expected − fallback_slack_nats` and `realign_fallback_full`, redo unbanded. Reads with no
  usable anchor (segment has no genome part covering the read) → unbanded.
* `ll_ref = max over REF segments`, `ll_alt = max over ALT segments`; with
  `alt_ref_junction_fraction = r > 0`: `ll_alt = logaddexp(ln(1-r) + ll_alt, ln r + ll_ref)`.
* `ReadObs.class`: `llr = ll_alt - ll_ref`; `> llr_informative` → Alt; `< -llr_informative` →
  Ref; else Uninformative. If the best hypothesis explains < `min_explained_frac` of the read's
  bases as aligned matches/mismatches (rest clipped), class = Unexplained (chimera / mismap)
  regardless of llr. `crosses_junction` = the best alignment covers a breakpoint column.

## Genotype model (`model.rs`, `output.rs`, owner C)

Dosage `g ∈ {0,1,2}`; purity grid `purity_grid` (weights uniform); alt-haplotype fraction
`φ(g,p)`: g=0 → `bg_alt_rate`; g=1 → `p/2 + (1-p/2)·bg`... precisely `φ1 = p/2` floored at
`bg_alt_rate`; g=2 → `max(p, 1-bg_alt_rate)`.
Per read `ll_r(φ) = logaddexp(ln(1-φ) + ll_ref, ln φ + ll_alt)`; Uninformative reads STAY in
the likelihood (not counted as votes; at a far duplication every reference-junction read has
llr = ln 0.5 and together they are the only evidence of absence); Unexplained reads are
excluded from GL but counted (`n_art`). With no voting read the locus is `no-coverage` unless
the uninformative reads move the posterior to ≥ `p_confident_absent`.
`GL_g = logmeanexp over p in grid of Σ_r ll_r(φ(g,p))`; `PL_g = round(-10·log10 e^(GL_g - max GL))`;
posterior with `prior` (default flat `1/3,1/3,1/3`); `GQ = min(99, -10·log10(1 - post_best))`.
`vaf` = MLE of φ on a 0..1 grid (step 0.01) of Σ_r ll_r(φ).

**No call strings (Jeremy, 2026-10-06: "don't keep legacy").** The per-colony output is numeric;
downstream (the joint step, annotate) decides from posteriors / likelihoods. A `status` word
marks rows without a model result: `no_reads` (nothing passed the gates), `high_coverage`
(depth gate), `error` (could not be modelled / fetched); their posterior columns are EMPTY.

Scores: `score_alt = Σ_{Alt reads} 10·log10(e)·llr`, `score_ref = Σ_{Ref reads} 10·log10(e)·(-llr)`
(Phred-scaled read evidence, integers), always both. `n_alt_l` / `n_alt_r` = Alt reads whose
best alt segment is ALT_L / ALT_R (an ALT_FULL read counts on both): the two-junction hallmark.

Output (gzip TSV, header exactly):
```
locus	kind	status	depth	n_alt	n_ref	n_uninf	n_art	n_disc	n_alt_l	n_alt_r	vaf	p_absent	p_het	p_hom	pl_absent	pl_het	pl_hom	gq	score_alt	score_ref
```
`kind` = TSD / BLUNT / TSD_DELETION / L1_MED_DELETION / L1_MED_DUPLICATION / ONE_SIDED.
`depth` = primary, mapped, non-duplicate reads overlapping the window(s) (capped at
`reads_for_high_coverage + 1` when the gate trips). The legacy consumers
(`src/combine_genotypes.py`, `annotate.py`'s `read_genotyping`, `tools/phylo`) read genotype
STRINGS and are NOT fed by this file; the joint step below is the per-patient consumer.

## Driver / I/O (`source.rs`, `read.rs`, `driver.rs`, owner D)

* `RegionSource` as legacy, but the `File` is wrapped in a `BufReader` with
  `io_buffer_bytes` (4 MiB) — Lustre latency, not CPU, bounds the legacy binary. Fills are
  adaptive: `io_fill_bytes` (256 KiB) right after a buffer miss, doubling per sequential fill up
  to the capacity (a fixed 4 MiB fill read 31 GB per PD37590 colony: loci are MiBs apart).
* Loci sorted by (contig order in header, min(L,R)); processed in that order; **one** indexed
  query per window: the same record stream feeds the depth count (early exit above
  `reads_for_high_coverage`, counting primary mapped non-dup records) and the evidence. Rows
  are streamed in processing order (downstream keys rows by name) with sync flush + stderr
  heartbeat every `heartbeat_every` loci. `--threads N` partitions the sorted loci into N
  contiguous chunks with one reader each; chunk outputs concatenated.
* Read gates (per record): primary (not secondary/supplementary), mapped, not qcfail, **not
  0x400** (Jeremy: dups are dropped), same contig, non-empty seq/qual; MAPQ ≥ `min_mapq` OR
  (MAPQ ≥ `min_mapq_clipped` AND a soft clip ≥ `min_clip_for_lowmapq` at a junction-facing
  end); a read must overlap a breakpoint (`[bp-1, bp+1)`) to be scored — others are not
  informative. Dedup by qname across the locus's windows (first wins); optional
  `lenient_dedup` (same start, end, strand of the fragment).
* Discordant anchors (counted, `n_disc`; weight `disc_weight_nats` into GL, default 0 =
  count only): primary, MAPQ ≥ `min_mapq`, mate unmapped / other contig / `|tlen| > disc_max_tlen`
  / same strand, FORWARD anchor ending in `[R - disc_span, R + 5]` or REVERSE anchor starting in
  `[L - 5, L + disc_span]`, qname not already an evidence read.
* `genotype_batch`: manifest `input<TAB>output`, consecutive, as legacy.

## Joint phylogenetic step (`newick.rs`, `joint.rs`, owner E)

`--step joint --tree P.tree --genotype-dir DIR (or --genotypes f1 f2 …) --out P.joint.tsv
--matrix P.joint_matrix.csv.gz [--root-prior 0.1] [--branch-prior length|uniform]`

Port of `tools/phylo/tree_fit.py` hypotheses, on the numeric per-colony files ONLY (no legacy
reader): per colony c, `P(d_c | absent) = 10^(-pl_absent/10)`, `P(d_c | present) =
½(10^(-pl_het/10) + 10^(-pl_hom/10))`; `status != ok` rows contribute log 1 to every hypothesis. Hypotheses: ROOT (every colony present), each branch
(clade below it), NOISE (none present, one shared alt-fraction fitted per locus — the
constant-allele-fraction artefact), INDEP (independent presence, π ~ U(0,1), DP over the
carrier count). Output per locus: best hypothesis id, log10 BF (tree vs max(NOISE, INDEP)),
carriers (tips), n_carriers, per-colony posterior P(carrier). Matrix: `;`-separated
(rows loci, cols colonies, empty first header cell) of NUMERIC per-colony P(carrier) (4
decimals); a colony whose row is `status != ok` or absent gets an empty cell. Tree tips are matched to genotype
file stems (`S1.txt.gz` → `S1`; farm: `<colony>.txt.gz`); unmatched tips are an error.

## Benchmark (owner F, after merge)

Local truth sets in `/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad/work/`:
`e2e_phylo` (10 colonies S1–S10 on `donor/tree.nwk`, WT control S11, 283 loci, 161 TP,
`phylo/truth.tsv`, `results_by_event.tsv`, reference `ref/reduced.fa`, contract
`genotype/P1.genotyping.tprt.txt.gz`, combined `combine/P1.combined.txt.gz`, legacy outputs
`genotype/tprt/S*.txt.gz`), `e2e_phylo2`, `e2e` (3 colonies + S4 WT). Scorers:
`test/e2e/score_genotypes.py`, `tools/phylo/tree_fit.py` + `calibration.py`. Report:
present/absent concordance per locus kind (legacy vs new), WT-control FPs, no-call rate,
score distributions (for `min_best_score`), joint-step accuracy vs `truth.tsv`
(clade/private/germline/nonclade/none), and ms/locus timing.
