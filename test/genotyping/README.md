# Genotyping panel — end-to-end gate for the VAF genotyper

Validates the G1–G19 genotyping rework on **real BAMs**: the per-locus VAF caller
(`genotyping_insertion.py`), the multiprocessing driver + ordered writer (`genotype.py`),
and the panel filters (`combine_genotypes.py`). `../val1` and `../fullstack` only exercise
*discovery*; this is the only harness that drives genotyping on real reads.

Design rationale and decision log: `plans/improve_discovery/genotyping_simulation.md`.

## The idea in one paragraph

Genotyping needs **het** loci (VAF ≈ 0.5), but a single implanted donor is homozygous by
construction. So we build **two** haplotypes from the same hs1 windows — `donor_ins` (with the
implants) and `donor_ref` (the pre-insertion allele) — and mix `wgsim` reads per sample:
a **het** sample draws depth/2 from each donor (→ VAF 0.5), a **wild-type** sample draws only
from `donor_ref` (→ VAF 0). Reads map hs1→**hg38**, so assembly discordance yields
discovery-derived false positives, present **identically in all 10 samples**. A true het reads
wild-type in the 3 reference-only samples; an FP never does — that is how `combine_genotypes`
separates them, and it is the property under test.

## Panel: 10 samples

| files | depth | donor mix | expected at true loci |
|---|---|---|---|
| 1 | 100× | ins+ref | het — **discovery source** (first `HET_DEPTHS`) |
| 6 | 80/60/40/20/10/5× | ins+ref | het, degrading to `insertion?`/`no-coverage` at 5× |
| 3 | 40× | ref only | wild-type (earns `min_wild-types`) |

## Element classes

`build_haplotypes.py` implants L1HS / AluYa5 / HERVK / SVA_E/F **plus processed pseudogenes**
whose inserted sequence is a real mature mRNA (`pseudogenes.fa`: DUX4, MALAT1, HNRNPA1, CASP12,
RPL21 — NCBI RefSeq Select). Because that mRNA equals the parent transcript, exonic spanning
reads can bwa-mismap to the parent locus in hg38 → alt reads lost → depressed VAF. The scorer
reports pseudogene het-recall separately so this shows up.

## Files

| file | role |
|---|---|
| `build_haplotypes.py` | two-donor builder (+ pseudogene family from `pseudogenes.fa`) |
| `pseudogenes.fa` | 5 mature parent transcripts (RefSeq Select; provenance in the plan) |
| `config.panel.py.template` | panel config — 3 shipping defaults overridden (see below) |
| `run_genotyping.sh` | orchestrator: haplotypes → 10-sample mix+bwa → discover(100×) → combine_insertions → genotype×10 → combine_genotypes → lift → score |
| `run_step.py` | runs the genotype / combine_genotypes step with the panel config first on `sys.path` (a plain `python src/main.py` puts src/ first and shadows it); has the `__main__` guard multiprocessing-spawn needs |
| `score_genotyping.py` | het-recall vs depth per class, negative discriminator, FP funnel |

Reuses `../fullstack/{lift_truth.py,scale10k/{make_2bit,run_combine,discovery_hs.config}}`.

## Panel config overrides (do **not** copy into `src/config*.py`)

- `min_wild-types 20 → 2` — 20 is impossible with 10 samples.
- `reads_for_high_coverage 60 → 250` — else every 100× locus is flagged high-coverage NA.
  (`LOW_COV_GATE=1` flips it back to 60 to exercise the gate.)
- `max_artefact/max_na` sized for a 10-sample panel.

## Run

```bash
# needs ~/Downloads/{hs1.fa, hg38.fa(+bwa index), hs1 bowtie2 index, hs1ToHg38 chain}
# and the built rust discovery binary.
bash test/genotyping/run_genotyping.sh          # env: OUT= THREADS= N_IMPLANTS= HET_DEPTHS= ...
```

~30–45 min on a laptop at the 1000-implant default (bottleneck: 10× bwa). For a fast smoke
run: `N_IMPLANTS=150 HET_DEPTHS="60 20 5" WT_DEPTHS="40 40" bash …`. Odd depths halve with
integer rounding (5× → 2+2). Outputs land in `$OUT` (default `~/Downloads/genotyping_panel`);
`score.txt` holds the verdict.

## Discovery-recall ceiling (~88%, by design)

Reads map hs1→**hg38**, so wherever the hg38 homolog of a donor window is segdup/repeat-dense,
bwa places the junction reads at **MAPQ 0** and discovery (`min_mapq`) drops them — the locus
never reaches the contract. This is the *same* cross-assembly discordance `../fullstack` uses to
manufacture FPs; it caps contract recall at ~88% on the validated windows (the residual loss is
the chr21:22.1–22.2M LINE-dense patch no window avoids). It is a mapping property, **not** a
genotyping or depth effect — 60× and 100× lose the same MAPQ-0 loci. The genotyping metrics are
computed on the loci that *do* reach the contract and are unaffected. For ~full recall instead,
map to hs1 (`bwa index hs1.fa`), at the cost of zero assembly-discordance FPs.

## ⚠ Synthetic-only limitation

Same caveat as `../val1`/`../fullstack`: clean synthetic reads exercise the code path, not
real enzymatic-prep messiness. Green here is a gate, not a clearance.
