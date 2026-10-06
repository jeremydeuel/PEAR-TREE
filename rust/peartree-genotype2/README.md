# peartree-genotype2 — realignment-based, phylogeny-aware genotyper

Replaces `peartree-genotype` (the 12-bp junction vote model with VAF bands). Design and the
contract between its modules: `plans/genotype_v2/SPEC.md`.

```
peartree-genotype2 --step genotype --bam <bam|cram> --insertions <P.genotyping[.tprt].txt.gz> \
    --combined <P.combined.txt.gz> --reference <ref.fa|ref.2bit> --out <colony.txt.gz> \
    [--threads N] [--config cluster/config.genotype2.grch38]
peartree-genotype2 --step genotype_batch --manifest <input<TAB>output per line> --insertions .. --combined .. --reference ..
peartree-genotype2 --step joint --tree <patient.tree> --genotypes genotypes/*.txt.gz \
    --out <P.joint.tsv> --matrix <P.genotypes.csv.gz> [--root-prior 0.1] [--branch-prior length|uniform]
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
purity grid; PL (Phred), posterior (flat prior), GQ, and the VAF MLE. Discordant anchors are
counted (`n_disc`) and can be weighted (`disc_weight_nats`, default 0).

### Per-colony output (gzip TSV, numeric, no call strings)

```
locus kind status depth n_alt n_ref n_uninf n_art n_disc n_alt_l n_alt_r vaf p_absent p_het p_hom pl_absent pl_het pl_hom gq score_alt score_ref
```

`status` = `ok` / `no_reads` / `high_coverage` / `error` (non-ok rows have empty posteriors).
`n_alt_l` / `n_alt_r` = alt reads per junction (the two-junction hallmark). `score_*` = Phred
sums of per-read evidence. Rows are in processing (coordinate) order; consumers key by `locus`.

### Joint step

Port of `tools/phylo/tree_fit.py`: per locus, hypotheses ROOT (every colony), each branch
(its clade), NOISE, INDEP; `P(d_c|present) = ½(10^-pl_het/10 + 10^-pl_hom/10)`,
`P(d_c|absent) = 10^-pl_absent/10`. Writes a per-locus table (best hypothesis, carriers,
log10 Bayes factor, per-colony P(carrier)) and the numeric P(carrier) matrix that
`tools/annotate_v2.py` reads (carrier at P ≥ 0.9). Tree tips without a genotype file are
missing data.

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
`~/results/tprt_ab/<P>/genotype2_vs_legacy.md`. `--dry-run` resolves every input and prints the
bsubs. `cluster/build.sh` builds this crate together with the other two.

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
