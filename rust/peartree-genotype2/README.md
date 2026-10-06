# peartree-genotype2 — realignment-based, phylogeny-aware genotyper

Replaces `peartree-genotype` (the 12-bp junction vote model with VAF bands). Design and the
contract between its modules: `plans/genotype_v2/SPEC.md`.

```
peartree-genotype2 --step genotype --bam <bam|cram> --insertions <P.genotyping[.tprt].txt.gz> \
    --combined <P.combined.txt.gz> --reference <ref.fa|ref.2bit> --out <colony.txt.gz> \
    [--threads N] [--config cluster/config.genotype2.grch38]
peartree-genotype2 --step genotype_batch --manifest <input<TAB>output per line> --insertions .. --combined .. --reference ..
peartree-genotype2 --step joint --tree <patient.tree> --genotypes genotypes/*.txt.gz \
    --out <P.joint.tsv> --matrix <P.genotypes.csv.gz> [--root-prior 0.1] [--branch-prior length|uniform] \
    [--dropout 0.02] [--false-present 0] [--noise-max-frac 1.0] [--ref-bias off|auto|<b>] [--zygosity colony|locus]
```

## What it does

Per locus, two haplotypes are built from the reference genome (FASTA + `.fai`, or UCSC 2bit)
and the FULL junction consensus combine wrote to `combined.txt.gz` (lowercase = inserted part,
uppercase = flank): `REF = genome[min-F, max+F)`, `ALT_R = genome[R-F, R) ++ insR`,
`ALT_L = insL ++ genome[L, L+F)`, merged into one `ALT_FULL` when the two consensuses overlap
(short insertions). Every geometry falls out of `alt = genome[.., R) ++ INS ++ genome[L, ..)`:
TSD, blunt, target-site deletion, L1-mediated deletion / duplication (far pairs: one window
per breakpoint; a long duplication's alt haplotype still carries the reference junctions, so
reference-junction reads are weighted 0.5/0.5), one-sided loci (one alt segment).

Every primary, non-duplicate read overlapping a breakpoint is realigned against every segment
with a quality-aware affine-gap aligner (banded on the BAM diagonal, unbanded fallback;
homopolymer-aware gap opening so poly-A length jitter is cheap; soft-clipped read bases cost
`ln 1/4`). `ll_ref` / `ll_alt` are the best segment likelihoods; a read is Alt / Ref when the
log-likelihood ratio passes `llr_informative`, Unexplained when no haplotype aligns ≥ 80% of
it (chimera / mismap), otherwise Uninformative (still in the likelihood, not a vote). Reads
whose MAPQ reflects only a short flank are kept when they carry a junction-facing soft clip
(`min_mapq_clipped`), which removes the legacy reference bias.

Genotype likelihoods: dosage 0/1/2 with the alt-haplotype fraction marginalised over a colony
purity grid; PL (Phred), posterior (flat prior), GQ, and the VAF MLE. Optionally corrected for
reference bias (see "Reference bias" below; off by default). Discordant anchors are
counted (`n_disc`) and can be weighted (`disc_weight_nats`, default 0).

### Per-colony output (gzip TSV, numeric, no call strings)

```
locus kind status depth n_alt n_ref n_uninf n_art n_disc n_alt_l n_alt_r vaf p_absent p_het p_hom pl_absent pl_het pl_hom gq score_alt score_ref
```

`status` = `ok` / `no_reads` / `high_coverage` / `error` (non-ok rows have empty posteriors).
`n_alt_l` / `n_alt_r` = alt reads per junction (the two-junction hallmark). `score_*` = Phred
sums of per-read evidence. Then one `pl_f<‰>` column per value of `noise_frac_grid` (default
`pl_f010 pl_f020 pl_f050 pl_f100 pl_f200`): −10·log10 P(reads | alt fraction φ) on the PL scale
(negative when φ fits better than every dosage), the per-locus profile the joint step's NOISE
hypothesis is built from. With `ref_bias_grid` set, one `pl_het_b<‰>` column per grid value follows
(`pl_het_b400` .. `pl_het_b1000`): the het likelihood at reference bias b, same scale (see
"Reference bias"). Rows are in processing (coordinate) order; consumers key by `locus`
and by column name (the readers are header-driven).

### Joint step

Port of `tools/phylo/tree_fit.py`: per locus, hypotheses ROOT (every colony), each branch
(its clade), NOISE, INDEP; `P(d_c|present) = ½(10^-pl_het/10 + 10^-pl_hom/10)`,
`P(d_c|absent) = 10^-pl_absent/10`. Writes a per-locus table (best hypothesis, carriers,
log10 Bayes factor, per-colony P(carrier)) and the numeric P(carrier) matrix that
`tools/annotate_v2.py` reads (carrier at P ≥ 0.9). Tree tips without a genotype file are
missing data. Genotype-error terms: `P(d_c|present) = (1-ε₁) P1 + ε₁ P0` (`--dropout`, default
0.02) and `P(d_c|absent) = (1-ε₀) P0 + ε₀ P1` (`--false-present`, default 0). PD37590 showed why:
at a germline locus 1-3 of 44 colonies have 0-1 alt reads (their alt reads realign as
uninformative / unexplained) and a hard PL 20-50 "absent" per colony handed 1,364 all-carrier loci
to INDEP; with ε₁ = 0.02 ROOT (or a clade) tolerates ~3 such colonies. ε₀ is off by default so a
single strongly present colony stays a private event.

**Private events: judge them on P(carrier), not on the BF.** INDEP explains a single carrier up to
a combinatorial factor, so a private event's `log10_bf_tree` saturates however strong the reads
are: ~2.2 on PD37590's 44-colony tree, ~1.2 on the 10-colony bench tree (1.60 on the 44-tip test
tree, identical at 20, 40 or 80 reads; `joint::tests::private_bf_saturates`), and `post_best`
saturates too (~0.85 on the bench). The carrier's `p_<colony>` / matrix cell is not saturated
(55 of 57 bench privates ≥ 0.9): use it (the 0.9 carrier threshold of `tools/genotype2_io.py`).
`test/genotype2/score_joint.py` scores Rust-joint privates that way; no other consumer
(`genotype2_io.py`, `genotype2_compare.py`, `check_known.py`, `compare_arms.py`, `annotate_v2.py`)
thresholds on the BF, they only print it.
NOISE is "absent everywhere" or one alt fraction φ shared by every colony (mean over φ of the
`pl_f` profile columns): the signature of mismapped paralogous reads or slippage, which tree_fit
calls `noise` and which three-genotype PLs cannot express. PD37590 first run: 65 of 88 joint
"clade" calls were such diffuse loci. Files without profile columns fall back to absent-everywhere.

`--zygosity locus` (default `colony`, the original model): `P(d_c|present) = ½(het + hom)` per
colony costs ln 2 in every carrier whose reads clearly favour one dosage -- 44 ln 2 = 30 nats for a
ROOT locus on PD37590, which NOISE (one shared fraction) never pays, so a moderately low but
consistent alt fraction in every colony went to NOISE. In `locus` mode a tree hypothesis takes one
dosage for all its carriers (½ het everywhere + ½ hom everywhere, ln 2 once); INDEP keeps the
per-colony mixture. Bench (10 colonies): germline ROOT 17 → 18 of 21, clades 39 → 44 of 55,
nonclade 23/23 unchanged.

### Reference bias

True somatic hets on PD37590 read at an alt fraction of ~0.18-0.35 at purity ~0.835: purity alone
predicts ~0.42. The rest is reference bias -- alt reads lost or unassigned (mapping bias of
insertion-carrying reads, junction reads the realignment cannot place), varying by locus kind
(tree_fit's germline-het vote odds K: TSD 0.95, BLUNT 0.92, TSD_DELETION 0.79, L1_MED_DELETION 0.59,
L1_MED_DUPLICATION 0.64). Model: `b` = capture / assignment efficiency of an ALT-haplotype read
relative to a REF-haplotype read (= tree_fit's K), and the het fraction at purity p is
`h·b / (h·b + 1 − h)`, h = p/2 (hom: h = p, still floored at `1 − bg`). `b = 1` reproduces the
uncorrected model exactly.

* Fixed, per colony: `ref_bias = <b>` and/or `ref_bias_kind = TSD:0.95,L1_MED_DELETION:0.6`
  in the genotype config correct the per-colony calls themselves.
* Estimated, joint step: genotype with `ref_bias_grid = 0.4,0.5,0.6,0.7,0.8,0.9,1.0`
  (`cluster/config.genotype2.grch38.refbias` = the GRCh38 config + this key, via `include =`): the
  per-colony calls are unchanged, the rows gain the `pl_het_b<‰>` profile. `--step joint
  --ref-bias auto` then estimates b from germline-het loci (tree_fit's candidate rule on
  `n_alt`/`n_ref`: autosomal, pooled alt fraction 0.2-0.8, ≥ 10 votes, ≥ max(2, C/2) colonies with
  ≥ 3 votes and every one of them with an alt vote, no colony with ≥ 8 votes and 0 alt):
  global `b = (Σalt + ½)/(Σref + ½)`; per kind the same, shrunk to the global value with 200
  pseudo-votes and then on the log scale with weight `n_loci/(n_loci + 100)` (the bias is a locus
  property: PD37590's 17 L1-mediated deletions carried 25k votes at raw b 0.45, which unshrunk turned
  3/37-alt background cells into ~75 private calls), and the global value for a kind with < 5 candidate loci; NO correction (b = 1) for a kind whose
  candidate reads are > 50 % uninformative, and for any locus whose reads pooled over the colonies are
  > 50 % uninformative (far duplications: the reference-junction reads score ln ½ by construction, so
  their 2-3 votes do not describe the likelihood -- they never get a spurious ~0.6 of their own);
  per colony a factor `(Σalt + 200)/(Σ raw_kind·ref + 200)` (against the kind's own vote ratio, so it does not re-absorb the shrinkage) over the kind-estimated candidates. Each
  cell's het likelihood is the profile linearly interpolated at `b_kind × factor` (clamped to the
  grid); absent, hom and the NOISE profile are untouched. `--ref-bias <b>` plugs one fixed value.
  The estimates go to stderr and to `<out stem>.refbias.tsv` (`P.joint.tsv` → `P.joint.refbias.tsv`:
  scope, name, loci, votes, raw, b, source). Without `--ref-bias` the joint output is byte-identical.

Effect (unit tests, 44 colonies, every colony at k/20 alt reads, b = 0.6, `joint::tests::
sweep_bias_vs_noise`): the correction alone does not rescue a consistent low fraction from NOISE
under `--zygosity colony` (still NOISE up to 6/20 = 0.30: the per-colony ln 2 dominates the ~2 Phred
per colony the bias gains); with `--zygosity locus` ROOT from 6/20, and with `--noise-max-frac 0.1`
as well, ROOT from 4/20. 2/20 and 3/20 everywhere stay NOISE in every mode, and a 5-colony clade at
5/20 with 39 colonies at 0/20 is the clade in every mode (a shared φ cannot fit the 39 clean
colonies). Per colony: 5/20 at purity 0.85 and b = 0.6 is het (P > 0.99, higher GQ than
uncorrected), 0/20 stays absent (P(absent) ≈ 0.995 -- a lower expected het fraction makes each ref
read weaker absence evidence), 10/20 is het not hom. PD37590 joint-only ablation on the real data (2026-10-06, same
per-colony files, rows = tree_fit class): versus no correction, `--ref-bias auto` moves germline
ROOT 7048 → 7191 and tree_fit-private privates 76 → 92 (cost: tree_fit-noise → ROOT 938 → 1045, and
uninformative-depth privates 10 → 84 to be checked); `--zygosity locus` sends 685 germline loci to
INDEP (12 before: real colonies differ in dosage at germline loci -- LOH, copy number, sampling),
and `--noise-max-frac 0.1` sends tree_fit-noise loci to ROOT (938 → 1913 colony / 2826 locus). The
six known loci are identical in every mode. Recommended: `--ref-bias auto`, zygosity `colony`,
no NOISE cap. Bench (simulated, no real bias): the estimate
is b ≈ 1.05 (TSD 1.03), so the joint calls do not move; even a deliberately wrong fixed
`ref_bias = 0.6` on every colony gives present-ok 542 → 552 of 616, 0 wild-type-control false
carriers (absent calls become no-calls instead: truth-absent no-call 72 → 134).

Python consumers read both formats through `tools/genotype2_io.py` (auto-detection, the
0.9 / 0.8 / 0.1 thresholds in one place): `tools/phylo/tree_fit.py`, `tools/annotate_v2.py`,
`cluster/tprt/compare_arms.py` / `check_known.py`, `test/e2e/score_genotypes.py`. Where
`<P>.joint.tsv` sits beside the matrix they add `joint_*` columns. `src/combine_genotypes.py`
refuses numeric files (the joint step replaces it).

## Cluster

`cluster/pipeline.sh` uses it by default (`GENOTYPE_IMPL=v2`): phase 3 runs it per colony with
`insertions/<P>.combined.txt.gz` and `GENOME_2BIT`, phase 4 runs the joint step with the
patient's tree from `patients/*/<P>/*.tree` (`PATIENT_TREE` to override). `GENOTYPE_IMPL=legacy`
restores the old binaries. Standalone: `cluster/genotype_one.sh` / `submit_genotype.sh` honour
the same variables (`COMBINED`, `GENOME_2BIT`, `GENOTYPE2_BIN`, `GENO2_CFG`).

**Real-data comparison with the legacy genotyper** (`cluster/genotype2_farm_compare.sh <P>`): for a
patient whose legacy run exists under the tprt_ab kit (`$TPRT_ROOT/<P>/C_rust/<P>` by default;
`LEGACY_RUNDIR=` to point elsewhere) it genotypes the same staged BAMs against the same contract
(picked from `insertions/` by the legacy files' locus count) and combined consensus as an LSF array
`gt2_<P>`, then a chained `--evaluate` job runs the joint step (length and uniform branch prior)
and `cluster/genotype2_compare.py` (standard library only): legacy call × v2 bucket cross-tab per
(locus, colony), hard discordances with read counts, per-locus carrier counts against the joint
step's hypothesis, legacy `tree_fit` labels against the joint step, and a trace of
`patients/*/<P>/known_insertions.tsv`. A second report is written for every sibling legacy set
`genotypes.*` (e.g. `genotypes.het30`). Output `$TPRT_ROOT/<P>/V2/report/report.md`, copied to
`~/results/tprt_ab/<P>/genotype2_vs_legacy.md`. When the python has pandas/scipy it also runs
`tools/phylo/tree_fit.py` on the v2 files (`$TPRT_ROOT/<P>/V2/fit/`), whose `summary.md` ends with a
tree_fit-class × joint-class table: the independent Python check of the Rust joint step. `--dry-run` resolves every input and prints the
bsubs. `REFBIAS=1` makes it a reference-bias run: config `cluster/config.genotype2.grch38.refbias`,
output `$TPRT_ROOT/<P>/V2_refbias` (the existing V2 files lack the profile columns), joint step with
`--ref-bias auto` (+ `JOINT_EXTRA`), results copied with a `.refbias` tag. `cluster/build.sh` builds this crate together with the other two.

## Validation

`test/genotype2/bench.sh` runs legacy and v2 on the simulated truth sets (10-colony phylogeny
+ wild-type control, 3-colony set) and scores both (`test/e2e/score_genotypes.py` reads both
formats; `test/genotype2/score_joint.py` scores the joint step). 2026-10-06, 10 colonies at
15x, 161 true loci × 10 colonies: truth-present colonies called present 398 → 542 of 616,
called absent 38 → 19, wild-type-control false carriers 0 → 0. Timing ~5 ms/locus (page-cached;
the legacy binary on the farm was Lustre-latency bound at ~35 ms/locus with two index queries
per locus; v2 does one early-stopping query per window through a 4 MiB seek-aware buffer).
First farm run (PD37590, 17,264 loci, 44 colonies): a fixed 4 MiB fill after every buffer miss
read 31 GB per colony (loci are MiBs apart in a 100 GB BAM) and took 35-45 ms/locus, the legacy
speed; fills are now adaptive (`io_fill_bytes` 256 KiB after a miss, doubling while a query keeps
streaming), with byte-identical calls.

`cargo test` (unit tests in every module) must pass; `cargo build --release`.
