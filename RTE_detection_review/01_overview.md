# 1. Overview and scope

## 1.1 Purpose

This review synthesises the literature and code supplied in `~/Downloads` and
`~/Documents/PEAR-TREE`, plus the Delly repository, around three questions:

1. **How can we detect true retrotransposition events in short-read WGS data?**
   → [§3](03_detecting_true_events.md), grounded in the biology of [§2](02_biology_of_retrotransposition.md).
2. **How can we detect common artefacts and non-retrotransposition events already at
   discovery time?** → [§4](04_artefacts_in_discovery.md).
3. **What are the common sequencing artefacts?** → [§5](05_sequencing_artefacts.md).

Delly is reviewed in [§6](06_delly_review.md); PEAR-TREE's own code is reviewed against all
of the above in [§7](07_peartree_code_review.md); and [§8](08_synthesis_and_recommendations.md)
turns the whole thing into concrete recommendations for PEAR-TREE2.

## 1.2 The problem in one paragraph

A retrotransposon insertion is a few hundred to a few thousand bases of mobile-element DNA
dropped into a new genomic site by target-primed reverse transcription (TPRT). In short-read
WGS it is visible only indirectly: as reads that **don't fit the reference** — soft-clipped
reads whose clipped part is the inserted sequence, and read pairs where one mate lands on a
repeat consensus and the other on unique genome. The difficulty is that **library-prep and
mapping artefacts produce the exact same signals** — clipped reads, split reads, discordant
pairs — often far more abundantly than real insertions. Reliable detection is therefore a
two-sided problem: **maximise the true TPRT signal (poly‑A tail, target-site duplication,
element identity, clonal allele fraction) while aggressively subtracting a well-characterised
set of artefacts.** Every tool reviewed here is a particular balance of those two forces.

## 1.3 Materials reviewed

**Retrotransposition detection literature**
- Tubio et al., *Science* 2014 — L1 3′ transductions in cancer; the **TraFiC** pipeline.
  (`science.1251343.pdf`)
- Zumalave et al., *Science* 2026 — concurrent L1 events driving reciprocal translocations;
  the **MEIGA / MEIGA-SR** long-read + short-read pipeline. (`science.aee4513.pdf`,
  `science.aee4513_sm.pdf`, tables zip)
- Nam et al., *Nature* 2023 — widespread somatic L1 retrotransposition in normal colorectal
  epithelium; a four-caller (MELT/TraFiC-mem/DELLY/xTea) clonal-WGS design with a detailed
  filtering cascade. (`s41586-023-06046-z.pdf`, `41586_2023_6046_MOESM1_ESM.pdf`)

**Sequencing-artefact literature**
- Chen et al., *BMC Genomics* 2024 — structure-specific (inverted-repeat/palindrome)
  artefacts; the **PDSM** mechanism and **ArtifactsFinder** blacklist tool.
  (`12864_2024_Article_10157.pdf`)
- Tanaka et al., *PLOS ONE* 2020 — enzymatic-fragmentation artefacts; SNV-centred
  palindromes and a logistic-regression filter. (`plosone1.pdf`)

**Structural-variant method**
- Rausch et al., *Bioinformatics* 2012 + `github.com/dellytools/delly` — the **Delly**
  integrated paired-end + split-read SV caller. (`bioinformatics_28_18_i333.pdf`)

**Code**
- **PEAR-TREE** (`~/Documents/PEAR-TREE`, branch `PEAR-TREE2`) — the user's own
  whole-genome retrotransposon-insertion pipeline, with its `PEAR-TREE2_PLAN.md` roadmap.
- **Seven comparator repositories, reviewed at source level** ([§12](12_tool_implementations_compared.md)):
  **TraFiC** and **v‑TraFiC** (gitlab.com/mobilegenomesgroup), **MEIGA-LR** and **MEIGA-SR**
  (shared `GAPI/` library; MEIGA-SR ships a `meiga.joblib` classifier), **xTEA**
  (github.com/parklab/xTEA), **ArtifactsFinder** (github.com/lilicai/ArtifactsFinder), and
  **MEIsimulator** — read as code, not papers, to extract concrete adoptable mechanisms and
  correct several literature-level assumptions.

> Not reviewed in depth: the `science.aee4513_tables_s1_to_s36.zip` supplementary tables
> (data tables, not methods) and several transient `~$*.docx` lock files present in the
> Downloads folder but unrelated to the topic.

## 1.4 How the pieces relate

```
        BIOLOGY (§2)  ──►  what a true event looks like
             │
   ┌─────────┴──────────┐
   ▼                    ▼
DETECTION (§3)     ARTEFACTS (§5)
 true-event signal   what mimics it
   │                    │
   └────────┬───────────┘
            ▼
  DISCOVERY-TIME FILTERING (§4)
  keep signal, drop mimics early
            │
   ┌────────┴─────────┐
   ▼                  ▼
 DELLY (§6)      PEAR-TREE (§7)
 PE+SR reference  clip-first, poly-A-aware,
                  cohort-genotyped
            │
            ▼
   RECOMMENDATIONS (§8)
```

## 1.5 One-line conclusions (details in each section)

- **True-event detection** rests on two hallmarks above all — the **poly‑A tail** and the
  **target-site duplication** — supported by element-consensus identity, 3′-transduction
  source tracing, and (for somatic calls) **clonal VAF ≈ 0.5** with matched-normal/cohort
  subtraction.
- **The most dangerous artefact** for a clipped-read caller is the **structure-specific
  chimera** from inverted repeats/palindromes, especially under **enzymatic fragmentation**;
  the decisive discovery-time test is "does the clipped part map back locally in inverted
  orientation?"
- **PEAR-TREE** already implements most recommended filters and adds a genuinely strong
  cohort-genotyping discriminator; its main opportunities are to **add the discordant-pair
  signal (sensitivity)**, **turn on coverage/blacklist masking (specificity)**, and
  **score the poly‑A/TSD/EN hallmarks explicitly (precision)**.
