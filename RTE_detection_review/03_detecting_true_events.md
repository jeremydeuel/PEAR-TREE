# 3. Detecting true retrotransposition events in short-read WGS (focus question 1)

*Sources: Tubio et al. 2014 (TraFiC); Zumalave et al. 2026 (MEIGA/MEIGA-SR); Nam et al.
Nature 2023 (colorectal soL1R); Rausch et al. 2012 (Delly). This section describes the
evidence and algorithms used to call a genuine RTE insertion, then places PEAR-TREE's
approach among them.*

---

## 3.1 The two orthogonal read-level signals

Every short-read MEI caller is built on one or both of these:

**(A) Discordant read pairs (the "cluster" signal).** A read pair in which **one mate maps
to the mobile-element consensus (L1/Alu/SVA) and the other maps uniquely to the flanking
genome**. Collect these, and reads pointing toward the insertion from the left form a
**positive cluster** and those from the right a **negative cluster**; a **positive/negative
cluster pair pointing at the same locus defines the insertion and its approximate
breakpoints** (TraFiC, Tubio Fig. 1E). This signal works even when **no single read spans
the junction**, so it is robust at low-to-moderate coverage — its great strength.

**(B) Split / soft-clipped reads (the "junction" signal).** A read that crosses the
breakpoint is clipped; the clipped portion is inserted (non-reference) sequence. This gives
**base-resolution breakpoints**, recovers the **poly‑A tail** and reconstructs the **TSD**.
Its weakness: it needs a read to actually span the junction, so it is coverage- and
insert-size-sensitive.

The mature tools **combine both**: cluster to localise, split to resolve. Delly is the
canonical PE→SR pipeline ([§6](06_delly_review.md)); TraFiC leads with clusters and confirms
with clips; **PEAR-TREE is almost purely signal (B)** ([§7](07_peartree_code_review.md)).

## 3.2 The reference pipelines

### 3.2.1 TraFiC (Transposon Finder in Cancer) — Tubio 2014
Detects **solo-L1, partnered transductions, orphan transductions** from PE WGS:
- discordant pairs (one mate on repeat consensus, one unique) → **positive/negative
  clusters** → insertion locus + breakpoints;
- **coverage increment** downstream of a source L1 flags transduction/amplification;
- **split/clipped reads** confirm poly‑A and TSD (partial 5′/3′ reconstruction for **89%**
  of insertions);
- **subfamily assignment (L1Hs) via diagnostic nucleotides** (86% unambiguous);
- **transductions traced to a source L1** by matching the unique downstream sequence back
  to a repertoire of full-length L1 loci — **95% of transductions → 72 germline sources**.
- Validation: PCR of 308 insertions → **98% true positive**; capillary breakpoints for 84%;
  in-silico specificity > 99%, sensitivity 73–83%.

> **The TraFiC-mem and v‑TraFiC source was reviewed for this report — see
> [§12.2](12_tool_implementations_compared.md).** From the code: clustering merges same-family
> reads within **200 bp**, a call needs a +/− cluster pair within **≤200 bp** with **≥4
> reads/cluster**, and the element read is explicitly the **MAPQ0 multimapper**; the filter
> cascade subtracts satellites, a **515k-entry family-matched germline MEI DB**, a **panel of
> normals**, and a **young-reference-copy RepeatMasker self-mask (≤20% divergence)**; source
> tracing links donor coordinates to a **124-locus active-L1 catalog**. The viral fork **v‑TraFiC**
> adds exactly the pieces PEAR-TREE lacks: **dustmasker low-complexity masking in discovery**, a
> **both-ends-clipped read reject**, and a **local MAPQ/SMS pileup mask**.

### 3.2.2 MELT, xTea — the community MEI callers
Not in the supplied PDFs but named as the Nature 2023 co-callers. MELT (Mobile Element
Locator Tool) is the 1000 Genomes standard: discordant-pair clustering against element
references + split-read breakpoint + poly‑A/TSD annotation, with a build-in reference-MEI
mask. xTea adds machine-learning-assisted filtering and long-read support. Both encode the
same hallmark logic as TraFiC. **xTea's actual code was reviewed for this report — see
[§12.4](12_tool_implementations_compared.md):** it is a clip **+** discordant caller with
**coverage-adaptive count thresholds**, a **clip↔disc geometric-consistency** filter, a
high-coverage-island mask (>3× median), a centromere blacklist, and a **random-forest genotyper
over 15 features**. It keeps clips at MAPQ ≥ 12, not ≥40.

### 3.2.3 MEIGA / MEIGA-SR — Zumalave 2026 (long-read + short-read screen)
The long-read successor, detecting **seven classes** (solo L1/Alu/SVA, partnered & orphan
transductions, processed pseudogenes, isolated poly(A/T)) plus retrotransposition-mediated
genomic rearrangements. Pipeline and defaults worth borrowing:
- recruit **spanning** and **split** reads from CIGAR; **MAPQ > 20**; **exclude events
  < 50 bp** (ONT indel noise);
- cluster split reads (**≥ 2 reads, same orientation, within 50 bp**); spanning insertions
  within 250 bp; **meta-cluster** requiring reciprocal distance ≤ 500 bp and **≥ 3 reads
  total**;
- **consensus assembly** (wtdbg2 → racon → minimap2) over insertion ± 2.5 kb flanks; then
  annotate repeat family, poly(A/T), transduction, pseudogene;
- final call filters: **3–500 supporting reads**, not entirely low-complexity (except
  poly(A/T)), unambiguous family, **≥ 40% of insertion assigned to an identity**, spanning/
  clip consistency; **somatic = zero supporting reads in matched normal**;
- **source assignment** two ways — transduction alignment to a **patient-specific DB of
  10-kb downstream intervals** of all full-length L1s, and **diagnostic-SNV inference**
  (100% of diagnostic SNVs shared, ≥ 75% total) — 91.7% concordant;
- benchmark precision **> 95%** (F1 up to 99.55), recall ~75–82% at 15% VAF. *(A separate
  MEIGA-SR re-validation vs TraFiC/xTea reports 99.9% precision / 95.7% recall; the ">99%"
  figure elsewhere in the paper is the **source-inference** specificity — sensitivity only
  47.6% — not the detection-benchmark precision.)*

> **The MEIGA-SR/LR source was reviewed for this report — see
> [§12.3–§12.4](12_tool_implementations_compared.md).** Corrections/additions from the code:
> the shipped MEIGA-SR is primarily a **targeted VAF genotyper**; discovery ends in a
> **logistic-regression classifier over 29 features** (poly‑A and clip support are the strongest
> positive weights); it rejects **both-ends-clipped (SMS)** reads and masks artefact regions from
> the **local MAPQ/SMS pileup**; and — notably — its **TSD detection is a commented-out stub with
> no length window and it scores no EN motif**, so PEAR-TREE's planned TSD-length + EN-motif
> scoring (R8) would lead it.

### 3.2.4 Delly — generic SV arm
No MEI mode; contributes single-breakpoint insertion signatures and (crucially)
**translocation/SV calls** that turn out to be RT-mediated. See [§6](06_delly_review.md).

## 3.3 Positive evidence that raises precision (use as many as possible)

Ranked by discriminating power (see [§2](02_biology_of_retrotransposition.md)):

1. **Poly‑A tail** at the 3′ junction (≥ 12–15 bp, high purity).
2. **TSD**: the same short direct repeat on both flanks; equivalently the `L < R`
   breakpoint ordering. A clean TSD is near-decisive.
3. **Both breakpoints recovered** and mutually consistent (positive + negative cluster, or
   left-clip + right-clip within the TSD window).
4. **Inserted sequence matches an element consensus** (L1/Alu/SVA/IAP) — establishes it is
   a MEI. In PEAR-TREE this is done post-hoc by the `annotate` tool via **DFAM HMMs +
   RepeatMasker**.
5. **3′ transduction traceable to a source L1** — simultaneously proves the event is real
   and names its origin.
6. **EN motif (TT/AAAA)** at the 5′ nick — per-cohort enrichment (190× in Nature 2023).
7. **Clonal VAF ≈ 0.5** across a clonally expanded sample — proves an in-vivo somatic event
   rather than a stochastic artefact.

## 3.4 The other half: classifying reference / germline / somatic

Recognising the TPRT scar is necessary but not sufficient; the event must be **classified**:

- **Reference/known:** mask against the assembly's annotated repeats and MEI polymorphism
  databases — **1000 Genomes MEI, dbRIP, euL1db**. (MELT/xTea/MEIGA all do this.)
- **Germline:** present in **matched normal/blood** or shared across unrelated donors.
  MEIGA-SR **excludes events shared by ≥ 2 donors**; Nam 2023 flagged germline SVs by
  **many discordant reads in matched blood** (the "≥ 3 discordant read pairs with an SA tag"
  figure is a *minimum-support* rule for retaining a call, **not** the germline test) and used
  a **2,860-genome population-allele-frequency panel**.
- **Somatic:** absent from matched normal; in clonal material, **VAF ≈ 0.5**. Requiring
  clonality (present in all cells of a crypt/organoid founder lineage) is what defines the
  event as a bona fide somatic mutation and enables developmental timing.

## 3.5 Validation is part of "detection"
Every reference study treats orthogonal validation as integral: **PCR + capillary**
(Tubio 98% TP), **long-read PacBio/ONT cross-platform** (MEIGA 92.5%, Nature 2023 PacBio),
**targeted capture**, **FISH/Micro-C** for RT-mediated rearrangements, and **mandatory IGV
inspection** of poly‑A + TSD for every accepted call. A short-read pipeline should output
enough evidence (both junctions, poly‑A, TSD, element identity) to make that inspection
possible.

## 3.6 Where PEAR-TREE sits

PEAR-TREE deliberately inverts the usual emphasis: it is a **TSD-first, clip-first,
poly‑A-aware** caller that uses **cross-sample population genotyping over a phylogenetic
tree** as its primary germline/somatic/artefact discriminator, in place of a single matched
normal.

- **Signal used:** split/clipped reads (B) with poly‑A as a first-class, rescuing signal;
  paired-end mates only *extend* the clipped consensus. It does **not** currently use the
  discordant-pair cluster signal (A). → high precision on spanning junctions, blind to
  events with no junction-spanning read.
- **Element identity:** deferred to the `annotate` tool (DFAM HMM + RepeatMasker +
  bowtie2 + liftover), not used to *gate* discovery — so discovery finds insertion junctions
  agnostically and classifies them later.
- **Classification:** replaces "matched normal" with "must vary across the cohort" — an
  insertion must have both **wild-type** and **confident insertion** samples and few
  **artefact/NA** samples ([§7](07_peartree_code_review.md) §7.5). For a phylogenetic study
  of many clones this is arguably *stronger* than one matched normal, because a recurrent
  artefact rarely produces a clean het/hom-vs-wt partition.

The consequence, developed in [§7](07_peartree_code_review.md) and
[§8](08_synthesis_and_recommendations.md): PEAR-TREE is well designed for **precision on
clonal cohorts** and its main improvement axes are (i) adding the discordant-pair signal for
sensitivity, (ii) turning on the coverage/blacklist artefact masks it already half-implements,
and (iii) exploiting the poly‑A/TSD/EN evidence more explicitly as a positive score.
