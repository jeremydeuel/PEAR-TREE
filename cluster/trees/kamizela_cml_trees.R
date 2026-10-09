suppressMessages({library(ape); library(phangorn)})
# Kamizela 2025 CML trees -> patients/peripheral_blood/<P>/<P>.tree (+ <P>.bcr_abl1.tsv).
# Usage: Rscript kamizela_cml_trees.R <outdir>. Inputs: nangalialab/CML @0991125 cache/PDD_A.RDS, cache/PDD_B.RDS,
# data/EGAD00001015353.sample_manifest.n1023.csv (local copies under ~/Documents/TP53_HSPC_trees/kamizela_cml).
# Tips missing from the manifest get the donor's single lo prefix; all names were checked against QC-passed rows.
src <- "~/Documents/TP53_HSPC_trees/kamizela_cml"
man <- read.csv(file.path(src, "_source/EGAD00001015353.sample_manifest.n1023.csv"))
man <- man[man$Sample_type == "Colony", ]
pdd <- c(readRDS(file.path(src, "PD51635/PDD_A.RDS")), readRDS(file.path(src, "_source/PDD_B.RDS")))
out <- commandArgs(TRUE)[1]
parse <- function(s) { m <- regmatches(s, regexec("^([a-z]*?)_?(lo)?([0-9]+)$", s))[[1]]
  if (length(m)) list(l = m[2], n = as.integer(m[4])) else list(l = s, n = NA) }
for (P in names(pdd)) {
  t <- pdd[[P]]$pdx$tree_ml
  t <- drop.tip(t, "zeros")
  samp <- man$sample[man$patient == P]; suf <- sub(P, "", samp, fixed = TRUE)
  ps <- lapply(suf, parse)
  new <- vapply(t$tip.label, function(tip) {
    hit <- samp[suf == tip]
    if (!length(hit)) { pt <- parse(tip)
      if (!is.na(pt$n)) hit <- samp[vapply(ps, function(q) !is.na(q$n) && q$n == pt$n &&
                                              (pt$l == "" || q$l == pt$l), TRUE)] }
    if (!length(hit) && grepl("^lo[0-9]+$", tip)) {   # tip absent from the manifest: use the donor's only lo prefix
      pre <- unique(sub("lo[0-9]+$", "", suf[grepl("_lo[0-9]+$", suf)]))
      if (length(pre) == 1) { hit <- paste0(P, pre, tip); message(P, ": ", tip, " not in manifest -> ", hit) } }
    if (length(hit) != 1) stop(P, ": tip ", tip, " -> ", length(hit), " matches: ", paste(hit, collapse = ","))
    hit }, "")
  stopifnot(!anyDuplicated(new))
  # driver clades on the ORIGINAL tree (node numbers refer to it)
  t0 <- pdd[[P]]$pdx$tree_ml; nd <- pdd[[P]]$nodes
  bcr <- nd$node[nd$driver == "BCR::ABL1"]
  bcr_tips <- if (length(bcr)) t0$tip.label[Descendants(t0, bcr[1], "tips")[[1]]] else character()
  t$tip.label <- unname(new)
  write.tree(t, file.path(out, paste0(P, ".tree")))
  write.table(data.frame(tip = unname(new), bcr_abl1 = names(new) %in% bcr_tips),
              file.path(out, paste0(P, ".bcr_abl1.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%s tips=%d BCR::ABL1=%d manifest_colonies=%d (QC passed %d)\n", P, length(new),
              sum(names(new) %in% bcr_tips), length(samp), sum(man$patient == P & man$Colony_QC == "Passed")))
}
