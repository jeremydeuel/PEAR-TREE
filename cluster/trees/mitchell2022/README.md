# Mitchell et al. 2022 — "Clonal dynamics of haematopoiesis across the human lifespan"

Nature 606:343-350 (2022), doi 10.1038/s41586-022-04786-y, PMC9177428.
10 Newick trees, fetched 2026-07-17 from **Mendeley `data.mendeley.com/datasets/np54zjkvxr/1`**.
Raw data: EGA `EGAD00001007851`.

## Mendeley has all 10; GitHub has only 8

`github.com/emily-mitchell/normal_haematopoiesis` (the Code availability link) is missing
**both cord bloods** (CB001, CB002). The Data availability link (Mendeley) is complete. This is
the second time in one day that a Nature paper's code repo held a subset and its data repo held
the whole thing — the first was Mitchell 2025 (Zenodo = 1 tree, Mendeley = 44). **Always check
the Data and Code DOIs separately.**

## THESE ARE NOT THE TREES WE USE

Our catalogue's trees come from `../chapman2025/EM/`, which re-analysed this cohort and built
DIFFERENT trees. Kept here as the primary-source reference:

    donor    PD id      this paper   chapman2025 (= ours)
    KX004    PD45534       451           922
    CB002    PD45517       390           379
    CB001    PD40315       216           195
    KX001    PD40521       407           407   (same)
    KX002    PD40667       380           380   (same)
    KX003    PD43974       328           328   (same)
    KX007    PD47738       315           315   (same)
    KX008    PD48402       367           367   (same)
    SX001    PD41048       362           362   (same)
    AX001    (none)        361           361   (same)

The published tip count equals the TOTAL colony count in `Summary_cut.csv` for every donor
except KX003 (343 colonies -> 328 tips). CB002 (95) and KX004 (99) are the two donors with
large **Progenitor** fractions, which is the likely axis of the 451-vs-922 difference.

## The complete donor list — 10, and there is no KX005 or KX006

Supp Table 1 "Demographic data and clinical details for donors used in the study", verbatim
(`Donor_ID / Age / Sex / Cause of death / Other diagnoses`):

    CB001   0    Female  NA                        NA
    CB002   0    Female  NA                        NA
    KX001   29   Male    Trauma - accident         Twisted bowel and appendectomy as infant.
    KX002   38   Male    Intracranial haemorrhage  Crohns disease diagnosed aged 28.
    SX001   48   Male    NA                        Selenoprotein deficiency.
    AX001   63   Male    NA                        Nothing of significance recorded.
    KX007   75   Male    Intracranial haemorrhage  Nothing of significance recorded.
    KX008   76   Female  Intracranial haemorrhage  Lyme disease aged 22.
    KX004   78   Female  Trauma - accident         Nothing of significance recorded.
    KX003   81   Male    Trauma - accident         Appendicectomy aged 6, ... Paget's aged 68, MI aged 78.

**KX004 is 77 OR 78** — Supp Table 1 says 78, the colony table `Summary_cut.csv` says 77. Both
are in the published record; they genuinely disagree. Not a transcription error on our side.

**KX002's `Previous chemotherapy` column reads `No`**, as it does for every donor. The main
text separately says KX002 "had inflammatory bowel disease treated with azathioprine" — Crohn's
IS an IBD and azathioprine is not chemotherapy, so these are consistent, not contradictory.

## AX001 has NO PD id in this paper

Its samples are named `BMH1_TG001_*` (Supp Table 3: `AX001  BMH1_TG001_P32_B03`). **`PD43976`
appears nowhere in this paper** — not in the text, any supplementary table, the Mendeley
archive, or any tree. This paper neither confirms nor refutes AX001=PD43976.

The AX001=PD43976 identification stands on the **Mitchell 2025** evidence instead: that paper's
`PD43976.tree` has tips named `BMH1_TG001_P31_A11`. See `../mitchell2025/README.md`.

## Sites: two donors have two, both blood. NO SPLEEN.

**The string "spleen" appears 0 times** in this paper's full text and 0 times in any
supplementary table. `sample_type` in `Summary_cut.csv` takes exactly three values: `BM`
(2,177), `PB` (809), `CB` (606).

    KX003 / PD43974   BM 282 + PB 61   "we sequenced both marrow and peripheral blood HSCs
                                        from the 81-year-old subject"
    KX001 / PD40521   BM 382 + PB 25

The spleen colonies for KX001/KX002/KX003 come from the SISTER paper, **Machado 2022**
(lymphocytes, doi 10.1038/s41586-022-05072-7). See donor_clinical.tsv for the unresolved
PD40667 donor-assignment conflict, which this paper CANNOT adjudicate: none of the disputed
spleen colonies appear in it.

What this paper DOES establish unambiguously: **PD40667 = KX002 = 38M, and KX001 = PD40521
only.** Supp Table 3 never pairs KX001 with PD40667; KX001's tree is 407/407 `PD40521*` tips
with zero PD40667 tips; KX002's tree is 380/380 `PD40667*`.

Multiple timepoints (same site): SX001 PB x3 (54/125/183 colonies), AX001 PB x2 (178/183).
