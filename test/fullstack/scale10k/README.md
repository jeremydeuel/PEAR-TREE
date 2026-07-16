# scale10k — 10 000-insertion full-stack specificity/recall gate

A scaled version of the `../` full-stack test: **10 000** real young elements implanted
into hs1 euchromatin (true positives) **plus a ~170 Mb centromere/telomere/alpha-satellite
FP compartment**, reads simulated and mapped hs1 → **GRCh38** with `bwa-mem`, then run
through **rust discovery + Python combine (step 2)** and scored in hg38 coordinates.

Where `../run_fullstack.sh` implants ~42 elements and produces a handful of organic
assembly-discordance FPs, this stresses specificity against the hardest real regions —
every chromosome's centromere and telomere, where hs1 and hg38 disagree most.

## Result (seed 1, 25×, 10 000 implants, 171.5 Mb FP compartment)

| stage | recall | genuine FP |
|---|---|---|
| rust discovery (default) | 96.0 % (9597/9994) | 383 |
| rust discovery (**`coverage_mask`**) | 96.0 % (9597/9994) | 235 |
| combine step 2 (default) | 95.4 % (9533/9994) | **1** |
| combine step 2 (**`coverage_mask`**) | 95.4 % (9533/9994) | **0** |

(genuine FP = calls matching no truth; the 1-to-1 scorer prints slightly higher — see note below.)

- 9994/10000 implants chain-lift into hg38 (6 have no hg38 homolog → excluded from scoring).
- Recall is uniform across every family/variant (94–97 %); the ~4.6 % misses are
  coverage / flank-mappability driven, not a class blind-spot.

### The one residual FP and the fix

The single genuine combine FP was a **pericentromeric classical-satellite mismap**
(`chr10:41.88 Mb`): reads from two different hs1 centromere regions both mis-align onto one
hg38 pericentromere at **MAPQ 60**, stacking to ~46× the genome-median depth; the consensus
is a `GGAAT` satellite tandem repeat, not element sequence. The MAPQ floor can't see it
(MAPQ 60) and the combine remap-to-hs1 can't (satellite maps everywhere).

The **SPEC-3 coverage mask** (`coverage_mask = true`, drop breakpoints above 5× median
depth) removes it — and 148 other discovery-level mismaps — with **byte-identical recall**,
because true insertions sit at ~1× median. `adaptive_evidence` (SPEC-4) also removes it but
costs ~45 low-VAF true insertions, so `coverage_mask` is the clean lever. The recommended
human discovery config is in `discovery_hs.config`.

> The 1-to-1 scorer prints a larger FP count (e.g. "11") because it also charges a 2nd call
> landing on an already-matched truth — a *duplicate call at a real insertion*, not an FP
> locus. `score_10k.py` now reports **genuine FP** (nearest-truth-beyond-window) separately.

## Files

| file | does |
|---|---|
| `profile_satellite.sh` | one pass over `hs1.repeatMasker.out.gz` → per-1 Mb satellite bp (`sat_bins.tsv`) |
| `gen_windows.py` | pick FP compartment (satellite bins ≥400 kb + telomere ends) + TP euchromatin windows |
| `build_donor_10k.py` | implant N elements at TTAAAA/EN sites + emit the FP compartment verbatim → `donor.fa`, `truth_hs1.tsv`, `flanks.fa` |
| `make_2bit.py` | minimal `.2bit` writer for the primary hg38 chromosomes (combine's genotyping loader needs a 2bit; no `faToTwoBit` required) |
| `discovery_hs.config` | rust discovery `--config` for human runs (min_mapq 60 + `coverage_mask`) |
| `config.local.py.template` | step-2 config shim (local samtools/bowtie2/hs1 index/chain/2bit paths) |
| `run_combine.py` | invoke `combine_insertions` with the shim shadowing `src/config.py` |
| `score_10k.py` | recall-by-family, precision, genuine-FP + FP-by-contig, vs the chain-lifted hg38 truth |
| `run_scale10k.sh` | orchestrates all of the above end-to-end |

## Run

```bash
# needs hs1.fa, bwa-indexed hg38.fa, bowtie2-indexed hs1, hs1ToHg38 chain, hs1 RepeatMasker
OUT=~/Downloads/fullstack_scale10k bash test/fullstack/scale10k/run_scale10k.sh
# env overrides: HS1= HG38= HS1_BT2= CHAIN= RMSK= THREADS= DEPTH= N_IMPLANTS=
```

`bwa-mem` over ~19 M read pairs is the long pole (a few hours on a workstation); rust
discovery is seconds, combine ~3 min.

## Caveat (inherited from `../README.md`)

Synthetic reads are clean and junctions exact — this proves the *code path* each region
hits, not real enzymatic-prep messiness. Green here is necessary, not sufficient; orthogonal
(long-read / PCR / IGV) confirmation of new calls on real data still stands.
