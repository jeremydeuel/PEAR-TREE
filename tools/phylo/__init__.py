"""Phylogeny-based evaluation of the insertion caller (plans/tprt_hallmarks/PHYLO_EVAL.md).

tree                 Newick parsing, pruning, branch enumeration, random coalescent trees
genotype_likelihood  per locus x colony read-vote likelihoods, purity / background estimators
tree_fit             CLI: fit every locus onto the SNV tree (phylo_fit.tsv, violations.tsv, ...)
discrimination       CLI: ROC / PR of the caller's scores against phylo_consistent / violating
simulate_counts      count-level simulator (per-colony genotype files from a tree) for calibration
"""
