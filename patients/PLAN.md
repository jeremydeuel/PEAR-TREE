# patients/ — layout and the plan to complete it

Built 2026-07-17 by `cluster/trees/build_patients.py` (rerunnable; regenerates this tree).

## Layout

```
patients/<organ_group>/<patient>/
    <tree>.tree        one (or, for liver explants, several) Newick tree file(s)
    colonies.tsv       per-colony metadata — PLACEHOLDER for now (see below)
    NO_TREE_ON_DISK.txt only present when the tree isn't downloaded yet
patients/<organ_group>/<organ_group>.md   per-organ summary: patient, PD id(s), tips
```

**155 patients across 12 organ groups.** One patient = one individual, **except Chapman HSCT
pairs**, where donor+recipient twins share one folder (`PairNN_PDxxxxx_PDyyyyy`) because they
share one engrafted tree. Liver explants **PD48367 / PD48372** each carry their 8 per-region
trees (`_b`..`_i`) in a single patient folder.

## The colonies.tsv is a PLACEHOLDER

Every `colonies.tsv` currently holds only the header:

    donor  proj  ds  readlen  mapped  assembly  sample

**These will be filled from the per-tip header probe running on the farm now.** Status by cohort:

| cohort / organ | probe status | source to populate from |
|---|---|---|
| **Peripheral blood — Chapman HSCT** | **DONE — 2,882/2,882 tips WGS, 0 gaps.** TSVs populated (10 pair folders). | `~/Downloads/chapman_hsct.tips.bam_coverage.tsv` (`check_tree_bams.sh`) |
| Gastric (stomach) | spot-check done (2/donor×proj) | `~/coorens2025stomach.wgs.*` + full per-tip run |
| Liver (Brunner) | spot-check done | `~/probe_liver.txt` + full per-tip run |
| all other blood / bronchial / cord / foetal / BM / colon | **not probed** | needs `probe_donors.sh` per cohort |

**Chapman is filled** (`donor·proj·ds·readlen·assembly·sample`), with two caveats recorded in
each TSV's comment line: `ds=WGS` and `readlen=151` are cohort-confirmed (probe was unanimous)
not per-tip measured; `assembly` is per-project (hs37d5 vs hg19 varies per-BAM — see
`chapman2024-hsct-grch37-only`); **`mapped` is still blank** — it needs a per-tip `idxstats`
run (`check_tree_bams.sh` used the fast project-classification path and didn't read depth).

For the remaining cohorts, run `probe_donors.sh` / a per-tip pass and fill the same 7 columns.
**Populate the WGS BAM by DS, not by basename** — the same colony can have a WGS and a targeted
release under different projects (see `chapman2024-hsct-grch37-only` memory).

**Populate is a separate step from build.** `build_patients.py` writes placeholder TSVs but now
**preserves any populated `colonies.tsv` across rebuilds** (snapshots non-placeholder ones before
the wipe and restores them). So the order is: build → populate; reruns won't clobber filled data.

## Trees still missing on disk (1 patient)

Hunted 2026-07-17/18. **12 of 13 obtained**, all tip-count-verified against donor_trees.tsv:
- **2** from the web — PD44887, PD44890 (Robinson 2022 MUTYH colon) → `cluster/trees/robinson2022mutyh/`.
- **9** from the local working set `.../STEM_Green_Lab/HNRNPA1/spar_2ndrev/human/all_trees/`
  (117-tree collection) → the **6 BCL** donors (`cluster/trees/bcl_unpublished/`, under `BCL0NN`
  aliases) and the **3 Lee-Six colorectal** (`cluster/trees/leesix2019colorectal/`).
- **1** from the JAK2/HNRNPA1 paper repo — **PD6634** (Williams 2022 MPN, 36 tips) extracted
  2026-07-18 from `Fig1_trees/mpn_trees/PD_MPN_CUT.RDS` (`obj[["PD6634"]]$tree`, an `ape::phylo`)
  → `cluster/trees/williams2022mpn/`. It was absent from the Williams reanalysis GitHub repo and
  the local all_trees set, but present in that RDS.

**Only 1 remains missing — left blank for now (Jeremy, 2026-07-18):**

| organ | patient | tips | status / blocker |
|---|---|---|---|
| tonsil | PD42775 | 9 | Machado 2022 (10.1038/s41586-022-05072-7). No phylogeny objects in `machadoheather/lymphocyte_somatic_mutation`; not in the local all_trees set, not in the MPN RDS. |

It still has a folder + placeholder `colonies.tsv` + `NO_TREE_ON_DISK.txt`.

**Routes for the last 1:** (a) rebuild from the paper's per-donor mutation matrix with
`ape::write.tree` (R is available); (b) email the authors; (c) check other local trees folders.

> Note: `.../HNRNPA1/spar_2ndrev/human/all_trees/` is a near-complete 117-tree collection that
> overlaps almost entirely with what we assembled under `cluster/trees/`. Worth treating as a
> cross-check / canonical source (it uses `..._pval_post_mix` Chapman-pair variants, vs our
> `..._vaf_post_mix`).

## Not represented here

**NF1 multi-organ autopsy** (Oliver 2025; PD50297/51122/51123 — brain/viscera/spinal cord;
838 LCM WGS, 151bp GRCh38) has **no published trees** — they would have to be built from the
substitution data before this cohort can join. See `new-organ-candidates-2026-07` memory.

## Regenerate

    python3 cluster/trees/build_patients.py     # rebuilds patients/ from the tree dirs; preserves this PLAN.md
