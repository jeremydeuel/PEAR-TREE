# 9. References and sources

## 9.1 Primary literature reviewed (files in `~/Downloads`)

1. **Tubio JMC, Li Y, Ju YS, et al.** "Extensive transduction of nonrepetitive DNA mediated
   by L1 retrotransposition in cancer genomes." *Science* 2014;345(6196):1251343.
   → `science.1251343.pdf` (dup: `science.1251343-2.pdf`). Introduces **TraFiC**;
   positive/negative discordant-read clusters; 3′ transductions barcode source L1s
   (95% → 72 germline sources); PCR validation 98% TP.

2. **Zumalave S, et al.** "Concurrent L1 retrotransposition events promote reciprocal
   translocations in human tumorigenesis." *Science* 2026 (aee4513).
   → `science.aee4513.pdf`, `science.aee4513_sm.pdf`, `science.aee4513_tables_s1_to_s36.zip`.
   Introduces **MEIGA / MEIGA-SR**; seven MEI classes + RT-mediated rearrangements;
   long-read internal-architecture resolution; diagnostic-SNV + transduction source assignment.
   **Headline finding:** somatic L1 activity drives **152 RT-RGs** including **13 reciprocal
   translocations** from two concurrent L1 insertions, detected via BND meta-clusters —
   largely invisible to short reads. Whether PEAR-TREE can capture these is assessed in
   [§3.2.3](03_detecting_true_events.md) (verdict: not as rearrangements without a break-end
   arm) and logged as a gap in [§7.7 #10](07_peartree_code_review.md).

3. **Nam CH, Youk J, Kim JY, et al.** "Widespread somatic L1 retrotransposition in normal
   colorectal epithelium." *Nature* 2023;617:540–547 (s41586-023-06046-z).
   → `s41586-023-06046-z.pdf`, `41586_2023_6046_MOESM1_ESM.pdf`. Clonal-WGS design;
   four callers (MELT, TraFiC-mem, DELLY, xTea); poly‑A + TSD IGV confirmation; poly‑A
   dropout in LCM prep; layered filtering cascade; clonal VAF ≈ 0.5 as somatic evidence.
   Code: github.com/ju-lab/colon_LINE1.

4. **Chen L, et al.** "Characterization and mitigation of artifacts derived from NGS library
   preparation due to structure-specific sequences in the human genome."
   *BMC Genomics* 2024;25:227. → `12864_2024_Article_10157.pdf`. **PDSM** model;
   inverted-repeat/palindrome chimeras; enzymatic vs sonication; **ArtifactsFinder**
   blacklist (github.com/lilicai/ArtifactsFinder).

5. **Tanaka N, et al.** "Sequencing artifacts derived from a library preparation method
   using enzymatic fragmentation." *PLOS ONE* 2020;15(1):e0227427. → `plosone1.pdf`.
   SNV-centred palindromes; soft-clip ratio and positional-bias signatures;
   logistic-regression filter (specificity 0.914, sensitivity 0.979).

6. **Rausch T, Zichner T, Schlattl A, et al.** "DELLY: structural variant discovery by
   integrated paired-end and split-read analysis." *Bioinformatics* 2012;28(18):i333–i339.
   → `bioinformatics_28_18_i333.pdf`. PE clustering (graph/maximal clique, 3 SD) + split-read
   DP refinement (k=7 k-mer filter, ≥2 split reads, 10% size agreement).

## 9.2 Software

- **Delly** — github.com/dellytools/delly (PE+SR SV caller; `delly`, `delly filter`,
  `delly merge`, `delly cnv`).
- **PEAR-TREE** — `~/Documents/PEAR-TREE` (this review's subject code; author J. Deuel).
- **MELT** — MEI caller named in the literature (not supplied as code).
- Support tools referenced: BWA-MEM, samtools/pysam, Picard/samblaster (dedup), bowtie2,
  HMMER + DFAM, RepeatMasker, wtdbg2/racon/minimap2 (long-read consensus), NCBI dustmasker,
  velvet, BLAT/BLAST, MUSCLE + EMBOSS `cons`.

### 9.2.1 Repositories reviewed at code level ([§12](12_tool_implementations_compared.md))

The source (not just the papers) of these seven tools was read for [§12](12_tool_implementations_compared.md):

- **TraFiC / TraFiC-mem** — gitlab.com/mobilegenomesgroup/TraFiC (discordant-pair L1 caller;
  Perl/bash + Snakemake; TEIBA breakpoint assembler; active-source-L1 transduction tracing).
- **v‑TraFiC** — gitlab.com/mobilegenomesgroup/v-trafic (viral-integration fork; adds dustmasker
  low-complexity masking, both-ends-clip reject, and a local MAPQ/SMS pileup mask *in discovery*).
- **MEIGA-LR** — gitlab.com/mobilegenomesgroup/MEIGA_LR (long-read caller; shared `GAPI/` library).
- **MEIGA-SR** — gitlab.com/mobilegenomesgroup/MEIGA-SR (short-read; targeted VAF genotyper +
  discovery caller ending in a **logistic-regression classifier**, `modules/meiga.joblib`).
- **xTEA** — github.com/parklab/xTEA (clip **+** discordant caller; coverage-adaptive thresholds;
  clip↔disc consistency; centromere blacklist; **random-forest genotyper** over 15 features).
- **ArtifactsFinder** — github.com/lilicai/ArtifactsFinder (palindrome/inverted-repeat finder;
  emits per-base artefact **positions**, not an interval BED; params **D_LEN 8 / STEM_LEN 5 /
  S_LEN 2 / ±50 bp**; palindrome length gate ships disabled).
- **MEIsimulator** — gitlab.com/mobilegenomesgroup/MEIsimulator (spike-in generator: labelled MEIs
  + reads + seeded truth report; VAF via two-genome merge — the §8.4 recall benchmark engine).

## 9.3 Reference databases for germline/polymorphic MEI subtraction

- **1000 Genomes MEI** call set (polymorphic mobile-element insertions).
- **dbRIP** — database of retrotransposon insertion polymorphisms.
- **euL1db** — European database of L1HS retrotransposon insertions.
- **DFAM** — profile-HMM repeat library (used by PEAR-TREE `annotate`).
- Population-allele-frequency panels of normal genomes (e.g. the 2,860-genome panel in
  Nam 2023) for rc-L1 / germline frequency estimation.

## 9.4 Key numeric thresholds cited across sources (quick reference)

| Parameter | Value | Source |
|---|---|---|
| Discordant insert-size cutoff | 3 SD from library median | Delly |
| Min split reads | ≥ 2 (default) | Delly / MEIGA |
| Split k-mer filter | k = 7, k_min = 3 hits | Delly |
| PE/SR size agreement | within 10% | Delly |
| MAPQ floor | ≥ 20 (Delly/MEIGA); ≥ 40 (PEAR-TREE) | all |
| Poly‑A tail | ≥ 15 bp, ≥ 90% purity, ≤ 30 bp from *either* end (MEIGA); ≥ 12 bp (PEAR-TREE) | MEIGA / PEAR-TREE |
| TSD length | ~2–20 bp (commonly 10–20) | Tubio / MEIGA |
| Min supporting reads (somatic) | ≥ 10% of total; ≥ 3 tumour, 0 in normal | Nam / MEIGA |
| Germline SV flag (blood) | *many* discordant reads in matched blood (qualitative); ≥ 3 disc. pairs + SA tag is a *min-support* rule, not the germline test | Nam 2023 |
| Panel VAF → artefact/germline | ≥ 1% (0.01) VAF in PoN | Nam 2023 |
| Indel/clip proportion reject | > 70% | Nam 2023 |
| Excludable regions / read cap | telomere/centromere `-x`; ≤ 1000 split reads/interval | Delly |
| Artefact soft-clip ratio | ~50% (vs ~5% genuine) | Tanaka |
| Palindrome-to-read-edge | ≤ 30 bp in 90.4% of SCPs | Tanaka |
| IVR/PS blacklist | IVR ≥ 8 bp (≥ 5 bp spacer); PS ≥ 17 bp; ± 50 bp | Chen |
| MEIGA event-size floor | exclude < 50 bp | MEIGA |
| MEIGA final support | 3–500 reads; ≥ 40% insert assigned identity | MEIGA |

**From the code review ([§12](12_tool_implementations_compared.md)):**

| Parameter | Value | Source |
|---|---|---|
| Both-ends-clip (SMS) reject | CIGAR op 4/5 at first & last → drop | MEIGA / v‑TraFiC / xTEA |
| Local-pileup mask | drop if MAPQ<10 >30% or SMS >15% | v‑TraFiC |
| High-coverage island mask | local depth > 3× sample median | xTEA |
| Coverage-adaptive clip/disc counts | cov 30→(3,4,1), 100→(8,12,3), 300→(25,30,8) | xTEA |
| RepeatMasker self-mask | same-family reference copy, divergence <15–20% | xTEA / TraFiC |
| Clip MAPQ floor | ≥ 12 (+ low-q-clip ratio ≤ 0.65) | xTEA |
| Low-complexity (k-mer entropy) | drop if Σ_K\|freq−0.25^K\| > 0.83, K∈{1..4}; poly‑A>0.4 exempt | MEIGA-LR |
| TraFiC cluster / +/− gap | merge ≤ 200 bp; +/− pair ≤ 200 bp; ≥ 4 reads | TraFiC |
| TraFiC transduction | both mates MAPQ ≥ 37; source within 20 kb; 500 kb self-guard | TraFiC |
| ArtifactsFinder IVS | arm-pair ≥ 8 bp; spacer ≥ 5 bp; sub-arm ≥ 2 bp; ± 50 bp | ArtifactsFinder |
| MEIsimulator reads / VAF | art HS25, 150 bp PE, insert 350±10; VAF = two-genome merge ratio | MEIsimulator |

*(PEAR-TREE's own configured thresholds are listed in [§7](07_peartree_code_review.md) and
its `config.py`.)*
