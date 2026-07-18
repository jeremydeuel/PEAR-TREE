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

Status by cohort (as of 2026-07-18):

| cohort / organ | status | how |
|---|---|---|
| **118 / 155 patients — DONE** | tree-tip → BAM, per colony (13,508/13,510 rows WGS, 0 MISSING) | `cluster/populate_colonies_tsv.sh` (commit 2f7375d) + the 10 Chapman pairs |
| **Liver (34) + AX001/PD43976 (2)** | **donor-level inventory** — tree tips anonymised, see below | `cluster/populate_donor_level.sh` + `cluster/donor_level.manifest.tsv` |
| tonsil/PD42775 | left blank (no tree) | — |

**The two populate scripts (both header/index-only nst_links reads; WGS picked by @RG DS;
assembly per-BAM from @SQ; `mapped` = idxstats mapped reads; `readlen`=NA in default header
mode, real via `READLEN_MODE=record`):**

1. `populate_colonies_tsv.sh` — maps each **tree tip → its BAM**. Works for the 108 cohorts
   whose tips are real sample ids (`PDxxxxx[_lo]`, `_hum` stripped). Run: `mkdir -p logs &&
   bsub < cluster/populate_colonies_tsv.sh`. Skips patients with no PD-style tips.
2. `populate_donor_level.sh` — **donor-level colony inventory** for the 36 patients whose tree
   tips are anonymised and have NO public/local clone→sample map: liver `Cl.NN` (Chapman
   "Prolonged persistence" ships only donor-level `SN_samples.txt`) and AX001/PD43976 `BMH…`
   plate-well codenames (Mitchell/your own JAK2 data work in BMH space; no `BMH→PD` map exists).
   It lists **all of the donor's WGS colonies** from nst_links — real sample metadata, but NOT
   1:1 with the `Cl.NN`/`BMH` tree tips (the tree file is unchanged; a comment records this).
   Run: `mkdir -p logs && bsub < cluster/populate_donor_level.sh`.

**Populate the WGS BAM by DS, not by basename** — the same colony can have a WGS and a targeted
release under different projects (see `chapman2024-hsct-grch37-only` memory). True per-tip
population of liver/AX001 would need the raw-data holders' private manifests (Chapman for liver
`Cl→sample`; Mitchell/EGA for `BMH→PD`).

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
