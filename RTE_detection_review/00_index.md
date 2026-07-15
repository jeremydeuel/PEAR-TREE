# Detecting retrotransposition in short-read WGS — a review

*A synthesis of the retrotransposition-detection and sequencing-artefact literature in
`~/Downloads`, the Delly repository, and the PEAR-TREE codebase — organised around three
questions: how to detect true retrotransposition events, how to reject artefacts and
non-RTE events at discovery, and what the common sequencing artefacts are.*

Prepared 2026-07-15 · scope and materials in [§1](01_overview.md).

---

## Table of contents

| # | Document | What's in it |
|---|----------|--------------|
| 0 | **[Index](00_index.md)** | This page |
| 1 | **[Overview and scope](01_overview.md)** | The problem, materials reviewed, how the pieces relate, one-line conclusions |
| 2 | **[Biology of retrotransposition](02_biology_of_retrotransposition.md)** | TPRT hallmarks — poly‑A, TSD, EN motif, 5′ truncation/inversion, transductions — the signatures detection exploits |
| 3 | **[Detecting true events](03_detecting_true_events.md)** *(focus Q1)* | Discordant-pair vs split-read signals; TraFiC / MELT / xTea / MEIGA / Delly; positive evidence; ref/germline/somatic classification; where PEAR-TREE sits |
| 4 | **[Rejecting artefacts at discovery](04_artefacts_in_discovery.md)** *(focus Q2)* | Per-artefact discovery-time tests, in cheapening order; mapped onto PEAR-TREE's actual filters |
| 5 | **[Common sequencing artefacts](05_sequencing_artefacts.md)** *(focus Q3)* | Structure-specific chimeras (PDSM), enzymatic fragmentation, poly‑A dropout, poly‑G, adapters, homopolymers, PCR/index-hopping/mapping |
| 6 | **[Delly review](06_delly_review.md)** | PE-clustering + split-read refinement; SV classes; no MEI mode; artefact controls; lessons for PEAR-TREE |
| 7 | **[PEAR-TREE code review](07_peartree_code_review.md)** | Architecture, the four steps, discovery algorithm, filter-to-artefact mapping, findings/gaps, strengths |
| 8 | **[Synthesis & recommendations](08_synthesis_and_recommendations.md)** | Best-practice pipeline; **18 prioritised recommendations (R1–R18)** for PEAR-TREE2; validation plan |
| 9 | **[References & sources](09_references.md)** | Full citations, software, databases, threshold quick-reference table |
| 10 | **[Empirical discovery patterns](10_empirical_discovery_patterns.md)** | Artefact patterns measured in **real** PEAR-TREE discovery output (560 JAK2/HNRNPA1 MPN files): homopolymer clips, recurrent Alu, palindrome/poly‑G — with prevalence and real examples |
| 11 | **[Cross-dataset & mouse ERV](11_cross_dataset_and_mouse_erv.md)** | Two more datasets (15-individual human, 32-individual mouse). Per-individual recurrence = germline-vs-artefact; **mouse ERVs make no poly‑A** — poly‑A must never be required; mouse has ~10× more palindrome + ~3.7× more poly‑G artefacts |
| 12 | **[Tool implementations compared](12_tool_implementations_compared.md)** | **Code-level review of 7 repos** (TraFiC, v‑TraFiC, MEIGA‑SR/LR, xTEA, ArtifactsFinder, MEIsimulator). Master comparison matrix; the universal **SMS both-ends-clip filter** (R15); local-pileup masking (R1); MEIGA/xTEA **ML classifiers** (R18); coverage-adaptive thresholds (R16); RepeatMasker self-mask (R17); ArtifactsFinder blacklist recipe; MEIsimulator recall benchmark |
| — | **[scans/](scans/README.md)** | The analysis scripts (`scans/scripts/`) and captured full-run outputs (`scans/outputs/`) behind §10–§11, with a README documenting each |

---

## The three focus questions, answered in one paragraph each

**1. How can we detect true retrotransposition events in WGS short-read data?**
By recognising the stereotyped scar of target-primed reverse transcription: a **poly‑A
tail** and a **target-site duplication (TSD)** flanking the insertion, with the inserted
sequence matching an **element consensus** (L1/Alu/SVA), and — where present — a **3′
transduction** whose unique downstream sequence traces the event to a specific source L1.
Two read-level signals carry this: **discordant read pairs** (one mate on the element
consensus, one on unique genome) clustered into a positive/negative pair that localises the
insertion, and **split/soft-clipped reads** that resolve the exact junction and read out the
poly‑A and TSD. Mature tools (TraFiC, MELT, xTea, MEIGA, and Delly as a generic-SV arm)
combine both; PEAR-TREE leads with the split-read + poly‑A signal and adds cross-sample
clonal genotyping. Somatic status is then established by **absence in matched normal** and a
**clonal VAF ≈ 0.5**. Full detail in [§2](02_biology_of_retrotransposition.md)–[§3](03_detecting_true_events.md).

**2. How can we detect common artefacts and non-RTE events already in discovery?**
Filter in cheapening order and demand positive hallmarks. Genome-free first: strip adapters
and poly‑G, reject low-complexity/homopolymer clips, reject **self-palindromic clips that
re-map locally in inverted orientation** (the top artefact for clipped-read callers),
require ≥ 2 concordant reads with an agreeing consensus, and mask coverage spikes. Then
genome-aware: reject clips that **map back within ~1 kb** of their breakpoint (local
rearrangement) or ends that **map end-to-end** (no novel junction), and drop candidates
inside **inverted-repeat/palindrome and telomere/centromere/segdup blacklists**. Finally,
separate non-somatic events by matched-normal/cohort recurrence and MEI polymorphism
databases. [§4](04_artefacts_in_discovery.md) maps every one of these onto PEAR-TREE's code.

**3. What are the common sequencing artefacts?**
Chiefly **structure-specific chimeras** — inverted-repeat/palindrome fold-backs whose
fill-in errors create clipped reads that mimic inversions, insertions and fusions,
massively amplified by **enzymatic fragmentation**; **poly‑A dropout** in low-input/enzymatic
prep (which specifically blinds L1 callers); **poly‑G** dark-cycle tails from two-colour
chemistry; **adapter read-through**; **homopolymer/microsatellite slippage**; **PCR
duplicates and chimeras**; **index hopping** on patterned flow cells; and **mapping/assembly**
artefacts near repeats and segmental duplications. Mechanisms, read-level signatures and
prevalence in [§5](05_sequencing_artefacts.md).

---

## Headline recommendations for PEAR-TREE (full list in [§8](08_synthesis_and_recommendations.md))

- **R1** Re-enable discovery-time high-coverage masking (currently `if False`).
- **R2** Wire in the already-written `is_low_complexity` / `has_well_defined_breakpoint` filters.
- **R6/R7** Add a self-palindrome filter and reference IVR/palindrome + telomere/centromere blacklists.
- **R8** Score poly‑A length, canonical TSD length and the EN motif as explicit *positive* evidence.
- **R11** Add discordant-read-pair discovery — the main sensitivity gap of a purely clip-based caller.
- **R15** *(cheapest new item)* Reject **both-ends-soft-clipped (SMS)** reads — a one-line filter MEIGA, v‑TraFiC and xTEA all use ([§12](12_tool_implementations_compared.md)).
- **R16/R17/R18** Coverage-adaptive thresholds, RepeatMasker young-copy self-mask, and an optional cohort-adapted ML classifier — all from the code-level review.

> **Empirical validation ([§10](10_empirical_discovery_patterns.md)):** a **full scan of all
> 560** JAK2/HNRNPA1 MPN discovery files (5.8M clipped sequences, 2.9M loci) confirms the
> artefact taxonomy — low-complexity clips **11.95%**, homopolymer-dominated **13.71%**,
> poly‑A/T **12.61%**, inverted-repeat/palindrome self-fold **0.65%**, poly‑G run **0.66%**;
> and **~45% of all clips are non-unique** across the cohort (Alu-consensus ~184k vs L1 ~7k
> instances). This directly justifies **R2** (low-complexity clips reaching output), **R6**
> (real fold-back chimeras present), and **R14** (recurrence-across-loci as an explicit filter).
>
> **Cross-dataset & mouse ([§11](11_cross_dataset_and_mouse_erv.md)):** across a 15-individual
> human set (50.7M clips) and a 32-individual mouse set (8.56M clips), the folder=individual
> layout makes recurrence a **germline-vs-artefact** test — clips in all individuals are
> reference-repeat (drop), clips private to one individual are germline (keep). **Mouse ERVs
> (IAP/LTR) make no poly‑A**, so poly‑A must never be required (refines **R8**) and recurrence
> must be scored **per individual, never "is-a-repeat"** (refines **R14**). Mouse also shows
> ~10× more palindrome and ~3.7× more poly‑G artefacts than human (raises **R6** priority).
>
> **Code-level tool review ([§12](12_tool_implementations_compared.md)):** the actual source of
> seven tools was read. Three findings recur across independent codebases — (1) **reject
> both-ends-clipped (SMS) reads** (MEIGA, v‑TraFiC, xTEA all do; one line; new **R15**); (2)
> **mask artefact regions from the local pileup**, no blacklist needed (concrete form of **R1**);
> (3) both mature short-read callers **end in a trained classifier** (MEIGA logistic reg / xTEA
> random forest; new **R18**). The review also **corrects ArtifactsFinder** (per-base positions,
> not intervals; real params D_LEN 8 / STEM_LEN 5 / ±50 bp) and adds a **MEIsimulator recall
> benchmark** for §8.4. PEAR-TREE's planned poly‑A/TSD/EN-motif scoring (**R8**) would *lead* the
> field: neither MEIGA nor xTEA scores the EN motif and MEIGA's TSD detection is a stub.
