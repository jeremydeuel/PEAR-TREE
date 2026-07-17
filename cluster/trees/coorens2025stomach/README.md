# Coorens et al. 2025 — normal gastric epithelium (STOMACH)

**"The somatic mutation landscape of normal gastric epithelium", Nature 2025**, PMC11981919,
doi 10.1038/s41586-025-08708-6. Martincorena / Stratton / Campbell.

**STOMACH IS AN ORGAN WE HAD ZERO COVERAGE OF.** 30 donors, 238 LCM WGS gastric glands.

Source: `github.com/TimCoorens/Stomach`, `Sequoia_output/<PD>_snv_tree_with_branch_length.tree`.
Fetched 2026-07-17. The repo also has `<PD>_indel_on_snv_tree_with_branch_length.tree` (same
topology, indel branch lengths) — NOT fetched; we want SNV branch lengths.

**Assembly: GRCh38** (paper states it explicitly; CNV/SV tables use `chr1`-style names).

## Data

- **EGAD00001015351 — WGS** (this is the one we want for MEI calling)
- EGAD00001015352 — targeted panel (829 samples; NOT usable for insertion discovery)

## Validation (2026-07-17)

30 trees, **240 tips total → 239 real samples**, vs **238 WGS rows** in Supplementary Table 2.
Per-donor tip counts match the supplement exactly on 28/30 donors. The two that don't:

- `PD41762` — tree has a **non-sample tip literally named `Ancestral`** (a root placeholder).
  5 real samples, matching the table. *** GOTCHA: 1 of the 30 trees carries `Ancestral`.
  Filter tips on `^PD` or you will ingest a fake sample. Same class of bug as the `Cl.N`
  liver tips and the PD48367b..h filename bug (commit 3f21262). ***
- `PD41756` — tree has 2 tips (`b_lo0001`, `b_lo0009`), supplement WGS sheet lists 1.
  Unexplained, +1. Not chased.

## Tree sizes are SMALL — this is not a colony cohort

Median ~6 tips/donor (range 2-23). These are LCM microbiopsies, not clonal colonies, so a
"phylogeny" here is a handful of glands. **Donors with 2-4 tips cannot support the monophyly
/ carrier-count logic our MEI filters assume** (see [[carrier-count-is-binomial-tail]]).
Usable-ish: PD41759(23), PD42787(20), PD41767(17), PD42789(16), PD40294(14), PD40293(12),
PD42790(11), PD41753/41763/41766(10).

## Two donors carry CANCER-scale mutation burdens

Median root->tip SNVs are ~450-3600 for most donors, but **PD41762 = 64,203** and
**PD41759 = 25,600**. These are gastric *cancer* patients and those clones look malignant
(hypermutator — likely MSI/POLE). Do not read those depths as "normal gastric epithelium".
Depth is tissue x age x disease. See [[cohort-paper-donor-map]].

## Cohort

21 of 30 are gastric-cancer patients (normal glands sampled adjacent to tumour); the rest are
non-cancer: obesity surgery (PD40297, PD42794), autopsy (PD41751-41758), transplant donors
(PD42787, PD42788, PD45518). Ages 23-85, both sexes. H pylori status in Supp Table 1.

**`PD45518` is NOT `PD45517`.** PD45517 = CB002, our cord-blood donor. One digit apart,
different study, different organ. See [[never-infer-tissue-from-pd-number]].
