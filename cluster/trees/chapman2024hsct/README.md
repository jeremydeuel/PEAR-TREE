# Spencer Chapman et al. 2024 — allogeneic HSCT, 10 sibling donor–recipient pairs

Nature 635:926-934 (2024), **doi 10.1038/s41586-024-08128-y**, PMID 39478227, PMC11602715.
20 trees (10 `trees_main/` + 10 `trees_no_dups/`), fetched 2026-07-17 from
`github.com/mspencerchapman/Clonal_dynamics_of_HSCT`.

2,824 colonies, sampled **9–31 years after HLA-matched sibling HCT**. EGA `EGAD00001010872`.

## The manifest is in the REPO, not the supplements

`tables/Demographics_table.xlsx` (the repo's copy of Extended Data Table 1). All three MOESM
supplements are PDFs (Supplementary Info, Reporting Summary, Peer Review) and contain **no PD
ids at all**. Per-donor clinical detail is in `../donor_clinical.tsv`.

## EVERY TREE IS TWO PEOPLE

Donor and recipient share an engrafted clone and are drawn on ONE tree. They are different
individuals with different ages and diagnoses. **Never treat a Pair as one donor.**

Rule, taken from the table's own `Donor or Recipient` column (NOT inferred from adjacency):
**lower PD id = DONOR, higher = RECIPIENT.** Holds for all 10 pairs.

Donors have a **blank Diagnosis cell** — that is the table marking them healthy, not a gap.

## THE PUBLISHED PAIR NUMBERS ARE NOT THE INTERNAL ONES

Our tree FILENAMES use the internal numbers. `data/metadata_files/Pair_metadata.csv`
(cols `Pair,Pair_new`) gives the renumbering:

    internal ->  published        internal ->  published
      11     ->  Pair_9             28     ->  Pair_10
      13     ->  Pair_7             31     ->  Pair_6
      21     ->  Pair_5             38     ->  Pair_8
      24     ->  Pair_1             40     ->  Pair_3
      25     ->  Pair_4             41     ->  Pair_2

**A LIVE COLLISION:** internal **Pair 3** (PD45790/91, unpublished) is NOT published
**Pair_3** (= internal Pair40 = PD45810/11). Reading a paper figure against our filenames will
silently pair the wrong people.

## Pair -> PD id

    file (internal)   published   DONOR      RECIPIENT   recipient dx
    Pair11            Pair_9      PD45792    PD45793     AML
    Pair13            Pair_7      PD45794    PD45795     NHL
    Pair21            Pair_5      PD45798    PD45799     AML   + CISPLATIN post-HSCT
    Pair24            Pair_1      PD45800    PD45801     secondary AML
    Pair25            Pair_4      PD45802    PD45803     AML
    Pair28            Pair_10     PD45804    PD45805     CML
    Pair31            Pair_6      PD45806    PD45807     AML   + CISPLATIN post-HSCT
    Pair38            Pair_8      PD45808    PD45809     AML
    Pair40            Pair_3      PD45810    PD45811     AML
    Pair41            Pair_2      PD45812    PD45813     AML   (SEX-MISMATCHED, F->M)

## Things that will bite

- **Tissue is peripheral blood CD34+ HSPC colonies for all 20.** The `Stem cell source`
  column (BM/PBSC) describes the **graft at transplant**, not what was sampled.
- **Siblings share ~half their germline.** The paper: "As the donors are siblings, recipients
  will share around half the same germline variants of the donor." Anything we do involving
  germline subtraction, contamination detection, or cross-donor "this variant is private"
  logic must account for this. In a pooled cohort contract these 20 are NOT 20 independent
  germlines — they are 10 sibships.
- **Two recipients had platinum AFTER the transplant**: PD45799 (cisplatin, lung carcinoma, 3y
  post-HSCT) and PD45807 (cisplatin, oesophageal cancer, post-HSCT). If they are used as
  "chemo-naive" controls, they are not.
- Each tree carries one extra **`Ancestral`** tip beyond the colony count.
- **No donor has multiple organs or timepoints** — one blood draw each. The only extra data
  are FACS-sorted mature subsets (granulocytes, monocytes, B and T cells) from the SAME draw,
  targeted-sequenced.

## These trees differ from the chapman2025 versions — do not interchange them

Colony counts here (`trees_main/`, excluding the `Ancestral` tip) vs `../chapman2025/MSC_BMT/`:

    pair      here   chapman2025
    Pair25     343       344
    Pair21     337       353
    Pair24     337       341
    Pair11     325       326      (agree closely)
    Pair13     298       299      (agree closely)
    Pair40     279       281
    Pair41     279       281
    Pair38     252       254
    Pair28     240       245
    Pair31     192       198

The 2025 repo carries re-derived trees. Our catalogue's tip counts follow **chapman2025**.

## Four recruited-but-dropped individuals

`data/metadata_files/sample_level_metadata.tsv` holds **12 pairs**; only 10 were published.
**PD45790/PD45791** (internal Pair 3) and **PD45796/PD45797** (internal Pair 18) have colonies
picked and plated (288/192 and 288/288 wells; BFU-E/CFU-GM/CFU-GEMM) but **no tree and no
demographics row**. The paper does not say why they were dropped. The repo's complete PD
inventory is exactly PD45790–PD45813, contiguous, 24 ids.

## An unresolved arithmetic discrepancy IN THE PAPER

The stated exclusions do not close: *"We excluded 46 colonies with low coverage, 58 technical
duplicates, 10 derived from a different germline (likely contamination) and 468 that were
non-clonal, leaving a final dataset of 2,824 genomes"* — but 3,399 − 582 = **2,817, not
2,824**. A 7-genome gap in the published text.

Separately, `trees_main/` totals 2,882 colonies; 2,882 − 58 technical duplicates = 2,824
exactly, which matches the headline — but that is OUR arithmetic, not a stated claim, and it
does not survive contact with `trees_no_dups/`, which removes **134** tips (total 2,748), with
Pair24 alone losing 85. Do not rely on either reconciliation.

## The six unplaceable donors are NOT here

PD49226, PD49227, PD49228, PD49229, PD53373, PD53374: **zero hits** in the paper body, all
three supplement PDFs, and the repo. (One apparent "53374" hit in Supplementary Information is
the genomic coordinate `22-41533746-G-A` — a false positive.) That is now **seven** papers
checked and empty.
