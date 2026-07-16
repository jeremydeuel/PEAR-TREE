# scale10k — 10 000-insertion full-stack specificity/recall gate

A scaled version of the `../` full-stack test: **10 000** real young elements implanted
into hs1 euchromatin (true positives) **plus a ~170 Mb centromere/telomere/alpha-satellite
FP compartment**, reads simulated and mapped hs1 → **GRCh38** with `bwa-mem`, then run
through **rust discovery + Python combine (step 2)** and scored in hg38 coordinates.

Where `../run_fullstack.sh` implants ~42 elements and produces a handful of organic
assembly-discordance FPs, this stresses specificity against the hardest real regions —
every chromosome's centromere and telomere, where hs1 and hg38 disagree most.

## Result (seed 1, 25×, 10 000 implants, 171.5 Mb FP compartment)

Recommended config = `discovery_hs.config` (`min_mapq 40` + `coverage_mask` + `mate_anchor_rescue`):

| stage | recall | genuine FP |
|---|---|---|
| rust discovery (min_mapq 60, no gates) | 94.1 % (9400/9994) | 199 |
| rust discovery (**recommended**) | 96.7 % (9662/9994) | 250 |
| combine step 2 (min_mapq 40 + coverage_mask) | 95.4 % (9533/9994) | 0 |
| combine step 2 (**recommended**) | **96.0 % (9595/9994)** | **0** |

(genuine FP = calls matching no truth; the 1-to-1 scorer prints slightly higher — see note below.)

Recall is uniform across classes at the recommended config (combine): L1 96.2 %, Alu 95.8 %,
SVA 96.0 %, HERV-K 95.7 %. Processed pseudogenes are not implanted (Feature B / `val1`).

### Sensitivity sweep (discovery level, all with `coverage_mask` on)

| config | recall | L1 | genuine FP |
|---|---|---|---|
| min_mapq 60 | 94.1 % | 92.8 % | 199 |
| + `mate_anchor_rescue` | 96.1 % | 96.1 % | 227 |
| min_mapq 40 | 96.0 % | 96.2 % | 235 |
| min_mapq 40 + `mate_anchor_rescue` | **96.7 %** | **96.9 %** | 250 |
| + `evidence_window`/`short_polya_clip`/`consensus_tolerant` | 96.7 % | 97.0 % | 258 |
| min_evidence = 1 | 96.8 % | 97.0 % | 290 |
| + `max_lowq_clip_ratio` guard | 95.9 % | 95.9 % | 240 |

**`mate_anchor_rescue` is the sensitivity lever** — the residual misses are insertions into
low-mapability-but-mate-unique flanks (biggest effect on L1). The SENS-1/7/8 toggles and a
lower evidence floor add ~nothing here (our misses aren't wobble/short-polyA/consensus/depth
limited), and the SENS-2 `max_lowq_clip_ratio` guard *removes* the rescued loci — do not pair
it with rescue. Combine cleans every extra discovery FP the relaxation admits (250 → 0).

- 9994/10000 implants chain-lift into hg38 (6 have no hg38 homolog → excluded from scoring).
- Recall is uniform across every family/variant (94–97 %); the ~4.6 % misses are
  coverage / flank-mappability driven, not a class blind-spot.

### Tried and rejected (no quick win in the residual ~4 %)

The ~400 remaining misses split into: (a) ~67 found by discovery but dropped in combine,
spread across *several* combine filters (clean-remap 246, clip-maps-near-breakpoint 36,
high-rate-region); and (b) ~330 never clustered by discovery — minus-strand insertions
whose 3′ poly-A tail reads out as a poly-T homopolymer clip while the 5′ clip is paralog-
ambiguous, so the pair never forms. Verified non-fixes:

| attempt | result |
|---|---|
| drop `coverage_mask` | 0 recall gain, +168 discovery FP — strictly worse |
| `reject_fully_mapping_reads=off` | +5 TP, +22 FP — bad trade |
| `short_polya_clip` / `evidence_window` / `consensus_tolerant` | +0.0–0.1 % recall |
| **hallmark-aware clean-remap rescue** (keep a clean-remapping call that carries a ≥12 bp poly-A/T tail) | **net +1 TP, +5 FP** — poly-A tracts are too common in genomic sequence; the addressable TP set (clean-remap-dropped AND poly-A-bearing) is ~1 |

Conclusion: the current recommended config (combine 96.0 % recall, 0 FP) is at the sensible
frontier. Further recall needs core discovery work (pairing homopolymer poly-A junctions with
paralog-ambiguous 5′ clips), not a filter tweak. A prototype poly-A Filter-A exception lived in
`combine_insertions.py` behind `clean_remap_keep_polya` (default off) and was reverted after
measuring the net-negative trade above.

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
