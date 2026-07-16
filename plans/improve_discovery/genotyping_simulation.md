# Genotyping simulation harness (`test/genotyping/`)

End-to-end validation of the G1–G19 genotyping rework on **real BAMs**: the per-locus
VAF caller (`genotyping_insertion.py`), the multiprocessing driver + ordered writer
(`genotype.py`), and the panel filters (`combine_genotypes.py`). The `val1`/`fullstack`
harnesses only exercise *discovery*; nothing today drives genotyping on real reads.

## Decisions (locked with the user)

- **hs1 two-haplotype, mapped to hg38**, reusing `fullstack`/`scale10k` machinery. Small
  scale — one hour on a laptop.
- **False positives are discovery-derived** (hs1→hg38 assembly discordance), not synthetic.
- **Discovery seeds the contract**: run real `discover → combine_insertions` on the 100× sample.

## The core idea: het via two haplotypes, and why the 3 negatives are the whole test

`build_donor.py` implants are homozygous by construction (every read carries the insertion
→ VAF≈1). To get **het (VAF 0.5)** we build two donors and mix reads:

- `donor_ins.fa` — TP windows **with** the 1000 implants + the FP compartment.
- `donor_ref.fa` — the **same** TP windows with **no** implants + the **same** FP compartment.

Per sample we run `wgsim` on each donor separately and concatenate the FASTQs:

| sample role | mixture | at a real implant | at an FP locus |
|---|---|---|---|
| **het** (depth D) | `ins`@D/2 + `ref`@D/2 | half the spanning reads carry the insertion → **VAF 0.5 → het** | identical to ref → assembly-discordance signal |
| **wild-type** (depth D) | `ref`@D only | no insertion reads → **wild-type** | identical to ref → same assembly-discordance signal |

**The discriminator:** a *real* het insertion is absent from the 3 reference-only samples
→ it reads **wild-type in exactly 3 samples** → clears `min_wild-types`. A *discovery-derived
FP* comes from assembly discordance present in **all 10** samples identically → it never
reads wild-type anywhere (or reads artefact everywhere) → **fails `min_wild-types` /
`min_insertions` / `max_artefact`** and is rejected. The 3 negative samples are precisely
what separates true insertions from FPs in `combine_genotypes` — that is the property under test.

## Element classes implanted (incl. processed pseudogenes)

`build_donor_10k.py` implants L1HS / AluYa5 / HERVK / SVA_E / SVA_F. The harness **adds a
processed-pseudogene class** so the panel covers a de-novo retrocopy, whose inserted sequence
is a **spliced mRNA** (a parent gene's exons concatenated, introns skipped) carrying the L1
TPRT scar (poly-A tail + TSD + EN-motif flank) — the "same poly-A/TSD hallmarks but spliced
exonic sequence" case from the RTE review §2.3.

- **Parent genes (user-chosen):** `DUX4`, `MALAT1`, `HNRNPA1`, `CASP12`, `RPL21` — a
  deliberately hard cross-mapping panel:
  - `HNRNPA1`, `RPL21` — among the most processed-pseudogene-rich genes in the genome
    (hundreds of existing retrocopies each) → an exonic clip maps to many paralogous loci →
    low MAPQ → dropped by `min_mapq` → alt reads lost → depressed VAF. The maximal-stress case.
  - `DUX4` — sits in the D4Z4 macrosatellite (chr4q/chr10q subtelomere), highly repetitive
    context and its own retrocopy family (DUXA-like); short (~1.7 kb) mRNA.
  - `CASP12` — itself a segregating pseudogene in most humans; the "parent" is already pseudogenic.
  - `MALAT1` — single-exon lncRNA (~8.7 kb): **no intron-skip signature** (splice-hallmark
    negative control), non-coding, and a *long* insert that stresses the junction geometry.
- **Construction:** `build_element` gains a `pseudogene` family that reads each mature
  transcript from a `pseudogenes.fa` (one record per gene), appends a poly-A tail, and implants
  with a TSD + EN-motif flank like a TPRT element — mirroring `val1/simulate.py`'s
  `--element-fasta` mechanism (there the exon set is synthetic; here the records are the real
  transcripts so clips are genuine exonic sequence bwa places against the parent locus in hg38).
  Transcript sequences are the real mature mRNAs, not hardcoded coordinates.
- **Provenance (`test/genotyping/pseudogenes.fa`, already fetched):** NCBI RefSeq Select
  transcripts via E-utilities (`esearch` `"refseq select"[Filter]` → `efetch rettype=fasta`):
  `HNRNPA1` NM_031157.4 (3744 bp), `MALAT1` NR_002819.5 (7472 bp, non-coding),
  `DUX4` NM_001306068.3 (1710 bp), `CASP12` NM_001191016.3 (3910 bp),
  `RPL21` NM_000982.4 (566 bp). All ACGT-only. `build_haplotypes.py` reads this file directly.
- **Why pseudogenes matter *for genotyping specifically* (not just class coverage):** the
  inserted mRNA is identical to the parent gene's transcript, so a spanning read whose
  soft-clip is exonic can be **mismapped by bwa to the parent-gene locus in hg38** instead of
  staying at the insertion site. That removes alt-supporting reads → depresses the apparent
  VAF → can push a true het toward `insertion?`/`wild-type`. No other implanted class has this
  failure mode. The scorer reports the pseudogene het-recall separately so we can see whether
  parent-gene cross-mapping erodes it.
- Count: a modest slice of the ~1000 implants (e.g. ~50–100), enough for a per-class recall
  number without dominating.

## The 10-sample matrix

| files | depth | donor mix | expected call at true loci |
|---|---|---|---|
| 1 | **100×** | ins+ref | het — **also the discovery source** |
| 6 | 80/60/40/20/10/**5×** | ins+ref | het, degrading toward `insertion?`/`no-coverage` at 5× |
| 3 | 40× | ref only | wild-type |

## Pipeline (exact `main.py` invocations)

```
# per sample K: build reads, map to hg38
wgsim donor_ins.fa ...  (depth D/2, het) ;  wgsim donor_ref.fa ... (depth D/2)
cat -> bwa mem hg38 -> fixmate|sort|markdup -> sampleK.bam
# (wt samples: wgsim donor_ref.fa only at depth D)

# 1. discover on the 100x sample  (Rust)
peartree-discovery --step discover --bam sample_100x.bam --out discovery.txt.gz --config discovery_hs.config

# 2. build the shared genotyping contract  (needs bowtie2/samtools/hg38 2bit — from scale10k)
python src/main.py --step combine_insertions --discovery_files discovery.txt.gz --out step2 --threads N
#   -> step2.genotyping.txt.gz   (>locus / @RIGHT_INSERTION / @RIGHT_REFERENCE / @LEFT_* )

# 3. genotype ALL 10 samples against that one contract
for K in 1..10: python src/main.py --step genotype --bam sampleK.bam \
    --out sampleK.genotypes.txt.gz --insertions step2.genotyping.txt.gz --threads N

# 4. combine into the panel matrix
python src/main.py --step combine_genotypes \
    --genotypes sample1.genotypes.txt.gz ... sample10.genotypes.txt.gz --out matrix.csv.gz
```

## New / reused files (`test/genotyping/`)

| file | role |
|---|---|
| `build_haplotypes.py` | wraps `build_donor_10k.py`: emit `donor_ins.fa` (implants, **incl. processed-pseudogene class**) **and** `donor_ref.fa` (same windows, implants suppressed) + `truth_hs1.tsv` with a `family` column that distinguishes `pseudogene` |
| `run_genotyping.sh` | orchestrator: haplotypes → per-sample wgsim mix + bwa hg38 (10 BAMs) → discover(100×) → combine_insertions → genotype×10 → combine_genotypes → lift truth → score |
| `config.panel.py.template` | **dedicated** genotyping config (see below) |
| `score_genotyping.py` | new scorer (see below) |
| reuse | `build_donor_10k.py`, `inject_artefacts.py`, `lift_truth.py`, `make_2bit.py`, `gen_windows.py`, `discovery_hs.config`, `run_combine.py` |

## Dedicated panel config (do **not** mutate the shipping human config)

Three shipping defaults break a 10-sample panel and must be overridden for the test:

1. `min_wild-types: 20 → 2` — 20 wild-types is impossible with 10 samples; every locus
   would fail. (2–3 keeps the 3-negative discriminator meaningful.)
2. `reads_for_high_coverage: 60 → 250` — at 100× a locus span exceeds 60 reads → the
   high-coverage gate fires → NA instead of a genotype. Raise it above 100× so the 100×
   sample yields real het calls (and add a lower-gate variant to test the gate itself).
3. `max_na`, `max_artefact`, `min_best_score` — set for a 10-sample panel, not a big tree.

## Scorer (`score_genotyping.py`) — what "pass" means

Match `matrix.csv.gz` rows to `truth_hg38.tsv` (±50 bp), then assert, per column:

- **True loci, 7 insertion samples:** fraction called het rises with depth. Report a
  het-recall curve vs depth; expect ~100% at ≥40×, degrading at 5–10× to
  `insertion?`/`wild-type?`/`no-coverage`. **This is the sensitivity floor the depth
  gradient is designed to measure — not a fixed 100%.**
- **Per element class** (L1/Alu/HERVK/SVA/**pseudogene**): het-recall broken out. Watch the
  pseudogene row — if it trails the others, parent-gene cross-mapping is stealing alt reads
  and the VAF bands / `min_supporting_reads` need attention.
- **True loci, 3 negative samples:** wild-type (this is what earns `min_wild-types`).
- **FP loci:** rejected by `combine_genotypes` (never emitted as a passing insertion).
- **Panel outcome:** every true locus passes the filters; every FP locus is dropped.
  Report the confusion matrix + which filter caught each FP.

## Scale & wall-clock (target: <1 h laptop, 8 threads)

- **TP window:** one ~8 Mb hs1 chr21 euchromatin window; `--min-site-gap ~7 kb` → ~1000 sites.
- **FP compartment:** a pericentromere/satellite window (scale10k `fp_windows`) sized to the
  FP count wanted. **FP count is emergent, not exactly 1000** — a bigger satellite
  compartment ⇒ more FPs ⇒ more bp ⇒ more reads. Dozens–hundreds of FPs is the cheap regime;
  ~1000 FPs needs the full satellite compartment (~30–40 min instead of ~15).
- Donor ≈ 11 Mb; total reads across all 10 bwa runs ≈ 16 M pairs; the 100× sample ≈ 3.7 M
  pairs. Bottleneck is the 10× `bwa` + `wgsim`; Rust discovery and per-locus genotyping are cheap.
- Everything is parameterised (`N_IMPLANTS`, depths, FP size) so the user can trade fidelity
  for wall-clock.

## Open risks

- **hg38 mapping loses some implants** (fullstack saw 30/34): a true locus that discovery
  misses on the 100× sample never enters the contract, so it can't be genotyped. Scored as a
  discovery miss, not a genotyping miss — report the two separately.
- **`combine_insertions` needs the bowtie2/hg38/2bit setup** from scale10k (`config.local.py.template`).
- **Pseudogene parent-gene cross-mapping** can depress alt-read counts (see element-class
  section) — a genotyping sensitivity concern the scorer's per-class row is meant to catch.
  Needs a real multi-exon gene with clean hs1 exon coordinates; if none is readily available,
  fall back to `val1`'s synthetic-exon construction (less faithful clips, same TPRT geometry).
- **VAF band edges** (0.10/0.30/0.85) and `min_supporting_reads=2` are the knobs the depth
  curve will stress; the scorer surfaces where they bite so they can be tuned on this synthetic
  before touching real data.
