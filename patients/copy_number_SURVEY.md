# Clade-level copy-number losses / LOH — per-patient survey

Built 2026-10-18 for the `cn_masked` rule in `cluster/somatic_table.py` (commit 8310a02). A germline
het insertion on a homolog that a clade lost (loss, LOH, cnLOH; gains never mask) is absent from
that clade and present everywhere else, so it looks somatic. `patients/<organ>/<P>/copy_number.tsv`
lists those events; `cluster/pipeline.sh` passes it automatically.

How the files were made:
- Events were read from each paper's supplements or the authors' own data repos (per-colony
  CNA/LOH tables, PDD RDS `meta$CNA`/`meta$LOH`, SV tables), never from body text. Every row's
  `source` column names the table/figure, and a `#` comment above it gives the paper's own label/colony ids.
- Clades are tip names from OUR tree. Each file loads with `read_copy_number` against the patient
  tree in its colonies.tsv assembly (hs37d5 -> GRCh37, else GRCh38). The MRCA clade size equals the
  number of affected colonies on our tree (colonies the paper lists but our tree lacks are noted in comments).
- Coordinates are in the patient's assembly. Stomach (paper GRCh38 -> our hs37d5 BAMs) and liver
  (paper GRCh37 -> our GRCh38 BAMs) were lifted with UCSC chains; the PD37116 10q11 lift is approximate (69% mapped).
- Lower-confidence placements: Mitchell 2022 single-colony cnLOH (PD40521 14q, PD48402 chr4) and the
  PD45534 2p clade were placed from figure branch positions; Fabre 2022 (PD34493/PD41276/PD41305) from
  Fig. 3 branch labels by topology match. PD41276 LOY covers every tip (masks nothing; kept for record).
- Small SV deletions (down to ~1 kb, and 1-148 bp BCR::ABL1 junction deletions on CML clades) are
  listed as losses: real, but they will rarely overlap an insertion.
- Liver tips are anonymised `Cl.N` clones; deletions were kept only when exactly one clone's
  microdissection set equals the SV's (Ng 2021 Supplementary Code `x.structure.RData`).

Working notes (per cohort, incl. events found but unplaceable) are kept outside the repo in the
side-session notes folder (`cn_mpn.md`, `cn_sds.md`, `cn_mitchell.md`, `cn_chapman_cml_et.md`, `cn_solid.md`).

| organ | patient | events in copy_number.tsv | source |
|---|---|---|---|
| blood_bm_spleen | PD40521 | 1 rows (1 cnLOH); chr 14q | Mitchell 2022 Nature 606:343 |
| blood_bm_spleen | PD40667 | 1 rows (1 loss); chr Y | Mitchell 2022 Nature 606:343 |
| blood_bm_spleen | PD43974 | 63 rows (63 loss); chr Y | Mitchell 2022 Nature 606:343 |
| bone_marrow | PD42191 | 7 rows (7 loss); chr 5, X, 18, 20, 6p | Machado 2023 Nat Commun; Machado 2023 PDD.RDS SDS8 nodes |
| bone_marrow | PD43293 | 1 rows (1 loss); chr 7p | Machado 2023 Nat Commun |
| bone_marrow | PD43294 | none reported (Machado 2023 PDD.RDS meta LOH/CNA, cnv_mapped_to_trees, SV table, Fig 2) |  |
| bone_marrow | PD43295 | none reported (as SDS1) |  |
| bone_marrow | PD43975 | gain only (15q24-26 tetrasomy, 9 col) - not listed |  |
| bone_marrow | PD44706 | 4 rows (1 LOH, 3 loss); chr 7, 20 | Machado 2023 Nat Commun |
| bone_marrow | PD45534 | 3 rows (3 cnLOH); chr 2p, 12q, 19q | Mitchell 2022 Nature 606:343; Mitchell 2025 Nat Genet |
| bone_marrow | PD45886 | none reported (as SDS1) |  |
| bone_marrow | PD45887 | 2 rows (2 loss); chr 20 | Machado 2023 Nat Commun |
| bone_marrow | PD45889 | none reported (as SDS1) |  |
| bone_marrow | PD47738 | 2 rows (2 loss); chr Y | Mitchell 2022 Nature 606:343 |
| bone_marrow | PD48402 | 1 rows (1 cnLOH); chr 4 | Mitchell 2022 Nature 606:343 |
| bone_marrow | PD49236 | 1 rows (1 cnLOH); chr 20 | Mitchell 2025 Nat Genet |
| bone_marrow | PD49237 | none reported (Mitchell 2022/2025 CN tables) |  |
| bone_marrow | PD5182 | 3 rows (3 cnLOH); chr 2, 9 | Williams 2022 Nature 602:162 |
| bone_marrow | PD5847 | 1 rows (1 cnLOH); chr 9 | Williams 2022 Nature 602:162 |
| bone_marrow | SDS5_PD42190_PD45888 | 2 rows (2 cnLOH); chr 7 | Machado 2023 Nat Commun |
| bronchial_epithelium | PD26988 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD30160 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD34204 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD34205 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD34206 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD34207 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD34209 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD34210 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD34211 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD34215 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD37451 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD37452 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD37453 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD37454 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD37455 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| bronchial_epithelium | PD37456 | none placeable: Yoshida 2020 publishes per-colony CN counts only, no regions |  |
| colon | PD44887 | none placeable: Robinson 2022 gives per-crypt CNV counts only (all 6 'e' adenoma glands have >=1 CNV - likely a shared clade event, region unpublished) |  |
| colon | PD44890 | none placeable: Robinson 2022 gives per-crypt CNV counts only |  |
| colorectum | PD34200 | none placeable: Lee-Six 2019 gives SV counts only |  |
| colorectum | PD37449 | 17q cnLOH reported in PD37449e_lo0014 (Lee-Six 2019 ED Fig 7c) but that crypt is not on our tree |  |
| colorectum | PD37590 | none placeable: Lee-Six 2019 gives SV counts only |  |
| cord_blood | PD40315 | none reported (Mitchell 2022 CN_summ.csv) |  |
| cord_blood | PD45517 | none reported (Mitchell 2022 CN_summ.csv) |  |
| foetal_liver | PD41768 | none reported (Spencer Chapman 2021: diploid karyotype in every colony) |  |
| foetal_liver | PD43947 | none reported (Spencer Chapman 2021) |  |
| gastric_epithelium | PD40293 | 10 rows (1 cnLOH, 9 loss); chr 1, 2, 3, X, 17, 20 | Coorens 2025 Nature |
| gastric_epithelium | PD40294 | 2 rows (2 loss); chr 3, 9 | Coorens 2025 Nature |
| gastric_epithelium | PD40295 | 2 rows (2 loss); chr 3 | Coorens 2025 Nature |
| gastric_epithelium | PD40296 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD40297 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41751 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41752 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41753 | 1 rows (1 loss); chr 9 | Coorens 2025 Nature |
| gastric_epithelium | PD41754 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41756 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41757 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41758 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41759 | 1 rows (1 loss); chr X | Coorens 2025 Nature |
| gastric_epithelium | PD41760 | 4 rows (4 loss); chr 3, 7, 16, 19 | Coorens 2025 Nature |
| gastric_epithelium | PD41761 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41762 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41763 | 2 rows (2 loss); chr 9, 16 | Coorens 2025 Nature |
| gastric_epithelium | PD41764 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD41765 | 2 rows (2 loss); chr 3 | Coorens 2025 Nature |
| gastric_epithelium | PD41766 | 3 rows (3 loss); chr 3, 9, 14 | Coorens 2025 Nature |
| gastric_epithelium | PD41767 | 10 rows (10 loss); chr 3, 4, 9, X, 17 | Coorens 2025 Nature |
| gastric_epithelium | PD42787 | 3 rows (1 cnLOH, 2 loss); chr 15, 17, 20 | Coorens 2025 Nature |
| gastric_epithelium | PD42788 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD42789 | 1 rows (1 loss); chr 3 | Coorens 2025 Nature |
| gastric_epithelium | PD42790 | 2 rows (2 loss); chr 3, 4 | Coorens 2025 Nature |
| gastric_epithelium | PD42791 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD42792 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| gastric_epithelium | PD42793 | 1 rows (1 cnLOH); chr 11 | Coorens 2025 Nature |
| gastric_epithelium | PD42794 | 2 rows (2 loss); chr 3, X | Coorens 2025 Nature |
| gastric_epithelium | PD45518 | none reported (Coorens 2025 Supp Tables 4/5) or gains only |  |
| liver | PD36713 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD36714 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD36715 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD36717 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD36718 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37105 | 3 rows (3 loss); chr 3, 8, X | Ng 2021 Nature |
| liver | PD37107 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37110 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37111 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37113 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37114 | 1 rows (1 loss); chr 4 | Ng 2021 Nature |
| liver | PD37115 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37116 | 2 rows (2 loss); chr 4, 10 | Brunner 2019 Nature |
| liver | PD37118 | 2 rows (2 loss); chr 2, 18 | Ng 2021 Nature |
| liver | PD37230 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37231 | 1 rows (1 loss); chr 5 | Ng 2021 Nature |
| liver | PD37232 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37234 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37237 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37238 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37239 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37240 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37243 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37245 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37904 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37907 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37909 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37910 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37914 | 3 rows (3 loss); chr 3, 13, 18 | Ng 2021 Nature |
| liver | PD37915 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD37917 | 1 rows (1 loss); chr 5 | Ng 2021 Nature |
| liver | PD37918 | not placeable: Brunner 2019 Supp Table 3 CN/LOH is per biopsy; SV deletions not mapping to a unique clone |  |
| liver | PD48367 | not placeable: Ng 2021 Supp Table 5 lists 52 deletions per microdissection, no clone key for these trees |  |
| liver | PD48372 | not placeable (as PD48367) |  |
| peripheral_blood | AX001 | 21 rows (21 loss); chr Y | Mitchell 2022 Nature 606:343 |
| peripheral_blood | PD34493 | 1 rows (1 loss); chr Y | Fabre 2022 Nature |
| peripheral_blood | PD37580 | none reported: 0 CNA / 0 SV in 84 colonies (Mitchell 2025 CN_summ.csv, ED Fig 2d, Fig 5d); TP53 clade is biallelic by 2 SNVs, not 17p loss |  |
| peripheral_blood | PD41048 | 7 rows (7 loss); chr Y | Mitchell 2022 Nature 606:343 |
| peripheral_blood | PD41276 | 2 rows (2 loss); chr Y, 20q | Fabre 2022 Nature |
| peripheral_blood | PD41305 | 2 rows (2 loss); chr Y | Fabre 2022 Nature |
| peripheral_blood | PD43976 | 1 rows (1 loss); chr Y | Mitchell 2022 Nature 606:343 |
| peripheral_blood | PD44579 | 1 rows (1 cnLOH); chr 10q | Mitchell 2025 Nat Genet |
| peripheral_blood | PD47536 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47537 | 1 rows (1 cnLOH); chr 22 | Mitchell 2025 Nat Genet |
| peripheral_blood | PD47538 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47539 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47540 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47541 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47695 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47696 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47697 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47698 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47699 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47700 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47701 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47702 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD47703 | 3 rows (3 cnLOH); chr 9q, 14q | Mitchell 2025 Nat Genet |
| peripheral_blood | PD4781 | 5 rows (2 cnLOH, 3 loss); chr 7, 9 | Williams 2022 Nature 602:162 |
| peripheral_blood | PD50306 | none reported (Mitchell 2025 Chemo_CN_SV.csv) |  |
| peripheral_blood | PD50307 | none reported (Mitchell 2025 Chemo_CN_SV.csv) |  |
| peripheral_blood | PD50308 | none reported (Mitchell 2025 Chemo_CN_SV.csv) |  |
| peripheral_blood | PD5117 | 3 rows (3 cnLOH); chr 9 | Williams 2022 Nature 602:162 |
| peripheral_blood | PD5147 | none reported (Williams 2022 PD_MPN_CUT.RDS LOH/CNA lists empty; Supp Note 5) |  |
| peripheral_blood | PD5163 | 1 rows (1 loss); chr X | Williams 2022 Nature 602:162 |
| peripheral_blood | PD51632 | 2 rows (2 loss); chr 9, Y | Kamizela 2025 Nature |
| peripheral_blood | PD51633 | 5 rows (1 cnLOH, 4 loss); chr 3, 9, Y, 22 | Kamizela 2025 Nature |
| peripheral_blood | PD51634 | 1 rows (1 loss); chr 9 | Kamizela 2025 Nature |
| peripheral_blood | PD51635 | 2 rows (2 loss); chr 9, 22 | Kamizela 2025 Nature |
| peripheral_blood | PD5179 | 1 rows (1 LOH); chr 9 | Williams 2022 Nature 602:162 |
| peripheral_blood | PD54832 | none reported (Leongamornlert 2026 PDD_SCC_TN.RDS CNA/LOH empty; Fig 5) |  |
| peripheral_blood | PD56961 | 1 rows (1 loss); chr 22 | Kamizela 2025 Nature |
| peripheral_blood | PD57332 | none (BCR::ABL1 junction duplications only; Kamizela 2025 PDD RDS + Supp Table 2) |  |
| peripheral_blood | PD57333 | 3 rows (3 loss); chr Y | Kamizela 2025 Nature |
| peripheral_blood | PD57334 | 1 rows (1 loss); chr 9 | Kamizela 2025 Nature |
| peripheral_blood | PD57335 | 3 rows (3 loss); chr Y, 22 | Kamizela 2025 Nature |
| peripheral_blood | PD60009 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD60010 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD60011 | not assessed: Mitchell 2025 ran CN calling only on its 11 large-tree donors |  |
| peripheral_blood | PD6629 | 1 rows (1 cnLOH); chr 9 | Williams 2022 Nature 602:162 |
| peripheral_blood | PD6634 | none reported (as PD5147) |  |
| peripheral_blood | PD6646 | 2 rows (2 cnLOH); chr 9 | Williams 2022 Nature 602:162 |
| peripheral_blood | PD7271 | none reported (as PD5147) |  |
| peripheral_blood | PD9478 | 4 rows (4 cnLOH); chr 9 | Williams 2022 Nature 602:162 |
| peripheral_blood | Pair11_PD45792_PD45793 | 16 rows (5 cnLOH, 11 loss); chr 4, Y, 19, 20, 21, 16p | Spencer Chapman 2024 Nature |
| peripheral_blood | Pair13_PD45794_PD45795 | 5 rows (2 cnLOH, 3 loss); chr 3, Y, 7q, 9q | Spencer Chapman 2024 Nature |
| peripheral_blood | Pair21_PD45798_PD45799 | 3 rows (3 loss); chr X | Spencer Chapman 2024 Nature |
| peripheral_blood | Pair24_PD45800_PD45801 | none reported (Spencer Chapman 2024 github data/SV_and_CNA_data) |  |
| peripheral_blood | Pair25_PD45802_PD45803 | 9 rows (9 loss); chr Y | Spencer Chapman 2024 Nature |
| peripheral_blood | Pair28_PD45804_PD45805 | 2 rows (2 loss); chr 6, 11 | Spencer Chapman 2024 Nature |
| peripheral_blood | Pair31_PD45806_PD45807 | events reported (5q UPD, a deletion) but in colonies not on our tree |  |
| peripheral_blood | Pair38_PD45808_PD45809 | 2 rows (1 cnLOH, 1 loss); chr X, 22 | Spencer Chapman 2024 Nature |
| peripheral_blood | Pair40_PD45810_PD45811 | gain / translocation / inversion only - not listed |  |
| peripheral_blood | Pair41_PD45812_PD45813 | 4 rows (1 cnLOH, 3 loss); chr 1, 4, 11, 11p | Spencer Chapman 2024 Nature |
| tonsil | PD42775 | skipped - no tree |  |
| unknown_tissue_bcl | PD49226 | 5 rows (1 cnLOH, 4 loss); chr 14, 15, 16 | Spencer Chapman 2023 Nat Med |
| unknown_tissue_bcl | PD49227 | 2 rows (2 loss); chr 10, 11 | Spencer Chapman 2023 Nat Med |
| unknown_tissue_bcl | PD49228 | 8 rows (2 cnLOH, 6 loss); chr 3, 6, 8, X, 10, 12, 18, 19 | Spencer Chapman 2023 Nat Med |
| unknown_tissue_bcl | PD49229 | 11 rows (3 cnLOH, 8 loss); chr 1, 3, 5, 7, Y, 13, 14, 16, 20 | Spencer Chapman 2023 Nat Med |
| unknown_tissue_bcl | PD53373 | 5 rows (5 loss); chr 1, 4, 6, 10, 19 | Spencer Chapman 2023 Nat Med |
| unknown_tissue_bcl | PD53374 | 6 rows (1 cnLOH, 5 loss); chr 4, 7, 8, 15, 16, 14q | Spencer Chapman 2023 Nat Med |
