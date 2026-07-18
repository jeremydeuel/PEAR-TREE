# williams2022mpn — PD6634 (recovered 2026-07-18)

`PD6634.tree` — Williams 2022 MPN colony phylogeny, **36 tips**
(`PD6634d_lo*`), matches `donor_trees.tsv` (36). Newick with branch lengths
(SNV-scaled edge lengths).

## Provenance
Extracted from the JAK2/HNRNPA1 paper repo's MPN tree collection:

    <OneDrive>/STEM_Green_Lab/HNRNPA1/JAK2_HNRNPA1_paper/Repo_080602025/
        Fig1_trees/mpn_trees/PD_MPN_CUT.RDS

`PD_MPN_CUT.RDS` is a named `list` of 12 MPN donors
(PD7271, PD5182, PD5163, PD5847, PD5179, **PD6634**, PD9478, PD6629, PD6646,
PD5117, PD4781, PD5147). Each element is `list(meta, cfg, tree)` where `tree`
is an `ape::phylo`. Extracted with:

    obj <- readRDS("PD_MPN_CUT.RDS")
    ape::write.tree(obj[["PD6634"]]$tree, "PD6634.tree")

This closes the last recoverable gap from PLAN.md's missing-trees hunt — PD6634
was absent from the Williams reanalysis GitHub repo and the local `all_trees`
set, but is present here. (The RDS also carries PD4781, whose tree we already
hold from chapman2025/NW; the two are independent copies.)

Only **PD42775** (Machado tonsil) now remains without a tree — left blank per
Jeremy 2026-07-18.
