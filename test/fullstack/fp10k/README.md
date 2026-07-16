# FP-stress full-stack harness (`test/fullstack/fp10k/`)

A companion to `test/fullstack/scale10k/`. Where scale10k is a **specificity gate** (10k real
young implants + a satellite compartment, tuned to yield ~1 genuine FP), this harness builds a
BAM that is deliberately **FP-rich**: **10 000 true insertions of all classes** and a large
population of **non-MEI decoys** tuned so that raw *discovery* throws **≥10 000 emergent false
positives**, for stress-testing the downstream MEI filtering (combine / annotate / coverage_mask).

Run: `OUT=~/Downloads/fullstack_fp10k bash test/fullstack/fp10k/run_fp10k.sh`
(bwa on the ~90 Mb donor is the long pole, ~20–30 min; whole run well under 2 h on a laptop).

## What is in the BAM

hs1-derived 150 bp reads mapped to **GRCh38** (`fp_stack.bam`). Two implanted populations,
spliced into large **contiguous euchromatin windows** so coverage is uniform ~25× and
`coverage_mask`'s median-of-populated-bins estimate behaves exactly as it would on real WGS.

### 1. TRUE insertions — 10 000, all classes (`role=TP`)
Canonical TPRT MEIs (element body + poly-A tail + target-site duplication) implanted at TTAAAA
(L1 endonuclease) sites in the 5 scale10k euchromatin TP windows (35 Mb). Classes:

| class | families | notes |
|---|---|---|
| young | L1HS, AluYa5, HERVK, SVA_E, SVA_F | full / 5′-trunc / 5′-inversion / 3′-transduction / provirus |
| older | L1PA2, AluSx, AluJb | older L1 & Alu subfamilies implanted as novel insertions |
| pseudogene | HNRNPA1/MALAT1/DUX4/CASP12/RPL21 | spliced mRNA + TPRT scar (from `test/genotyping/pseudogenes.fa`) |

These are genuine hs1→hg38 non-reference insertions — the true positives an MEI caller should recover.

### 2. DECOY false positives — ~14 000 (`role=FP`)
Non-MEI structural differences spliced into a **separate** set of ~18 euchromatin FP windows
(`gen_fp_windows.py`, ~36 Mb, spread across the autosomes, disjoint from TP windows and the
satellite compartment). Every decoy **omits the poly-A + TSD hallmarks**, so it is a
copy-number / non-canonical difference a good MEI caller must reject. The inserted material is
real, **genome-wide** hs1 repeat sequence, so the FP population *is* "all low-complexity + young
& old retrotransposon sites" as requested — only the genomic *context* is mappable euchromatin
(so each junction is a clean, unmasked clip cluster rather than a satellite MAPQ-0 pileup):

| decoy class | inserted sequence | source pool |
|---|---|---|
| `lcr` | extra tandem copies of a Low_complexity / Simple_repeat / Satellite tract | `lcr_sites.bed` |
| `young_rte` | a divergent fragment (150–600 bp) of a young element (AluY / L1P / SVA / ERVK) | `young_sites.bed` |
| `old_rte` | a divergent fragment of an old element (AluS/AluJ / MIR / L2 / CR1 / ERVL) | `old_sites.bed` |

The three `*_sites.bed` are one RepeatMasker pass (`hs1.repeatMasker.out.gz`), thinned to one
representative element per (contig, 50 kb bin) per category so the pools span the whole genome.

## Why cassettes were rejected (design note)

The first cut emitted each site as its own short `flank|insert|flank` cassette. That makes the
donor **sparse** (a few % of hg38, scattered), so many thinly-covered multi-map bins drag the
`coverage_mask` median (taken over *populated* bins) down to ~1–2 reads/bin — and the mask then
nukes the real 25× loci (1 call in the whole BAM). It also anchored decoys in their native
*repetitive* hs1 context, so only ~17 % fired. The window-splice geometry (this harness) fixes
both: uniform coverage → realistic `coverage_mask`, unique euchromatin flanks → ~97 % decoy
firing and ~85 % TP recall.

## Files

| file | role |
|---|---|
| `build_fp10k.py` | window-splice donor builder (TP MEIs + non-MEI decoys) |
| `gen_fp_windows.py` | euchromatin FP-decoy windows, disjoint from TP + satellite |
| `run_fp10k.sh` | build the BAM: windows → donor → wgsim → artefacts → bwa hg38 → markdup → discovery (recommended + plain configs) → chain-lift TP truth → score |
| `run_pipeline.sh` | run the FULL pipeline on that BAM: combine_insertions → genotype → combine_genotypes, then score the post-combine kept set (needs hg38 `.2bit`, bowtie2 hs1 index, chain). Annotate: follow the scale10k step-9 recipe on `test/fullstack/annotate/` (paired combined only — filter one-sided discordant calls first). |
| `score_fp10k.py` | TP recall (per class) + emergent-FP count, with decoy-explained vs unlabelled breakdown |

## End-to-end result (2026-07-16, recommended config)

| stage | TRUE recall | emergent FP | note |
|---|---|---|---|
| discovery | 97.1 % (9701/9994) | **21 061** | of 30 825 calls |
| + combine_insertions | 95.9 % (9582) | **2 277** | keeps 11 859; clean-remap-to-hs1 drops 18 330 (the one-sided decoys) → 89 % FP cut |
| + genotype | — | — | per-locus het/hom VAF calls (`fp.genotypes.txt.gz`) |
| + combine_genotypes | — | — | 0 (multi-sample panel step; N/A for one sample) |
| + annotate | — | — | recovers L1HS/AluY/SVA/HERVK/L1PA/pseudogene families; decoys → artefact / unknown / old-repeat family |

The FP count is highly config-dependent (plain `min_mapq=40` discovery = 2 536 emergent FP; the
recommended config's discordant one-sided anchoring is the discovery-FP driver, and those
one-sided calls are exactly what combine_insertions then removes).

Reused from siblings: `scale10k/gen_windows.py`, `scale10k/profile_satellite.sh`,
`scale10k/discovery_hs.config`, `fullstack/inject_artefacts.py`, `fullstack/lift_truth.py`,
`genotyping/pseudogenes.fa`.

## Outputs (`$OUT`)

- `fp_stack.bam` (+`.bai`) — the deliverable BAM (hs1 reads on hg38).
- `truth_hs1.tsv` / `truth_hg38.tsv` — every implant with `class` and `role` (TP/FP); the hg38
  file is chain-lifted and is what `score_fp10k.py` scores against.
- `discovery.txt.gz` / `discovery.plain.txt.gz` — discovery calls (recommended vs plain mq40).
- `score.txt` — recall + emergent-FP headline.

Tunables (env): `N_IMPLANTS` (10000), `N_DECOYS` (14000), `FP_MB` (35), `DEPTH` (25), `THREADS`.
The emergent-FP count scales with `N_DECOYS`; raise it to push further past the 10k floor.
