suppressMessages(library(ape))
base <- "/Users/jeremy/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents - STEM_Green_Lab/HNRNPA1"
rows <- list(); add <- function(src,donor,tr){
  if(!inherits(tr,"phylo")) return(invisible())
  dep <- if(is.null(tr$edge.length)) NA else median(node.depth.edgelength(tr)[seq_along(tr$tip.label)])
  rows[[length(rows)+1]] <<- data.frame(source=src, donor=donor, ntips=length(tr$tip.label),
      median_depth=dep, has_brlen=!is.null(tr$edge.length), tip=tr$tip.label, stringsAsFactors=FALSE)
}
# 1. the master rds
x <- readRDS(file.path(base,"spar_2ndrev/human/trees.rds"))
for (p in names(x)) add("trees.rds", p, x[[p]])
# 2. the MPN rds (nested: $tree)
y <- readRDS(file.path(base,"JAK2_HNRNPA1_paper/Repo_080602025/Fig1_trees/mpn_trees/PD_MPN_CUT.RDS"))
for (p in names(y)) add("PD_MPN_CUT.RDS", p, y[[p]]$tree)
# 3. loose .tree files (skip ignore/)
d <- file.path(base,"spar_2ndrev/human/all_trees")
for (f in list.files(d, pattern="\\.tree$", full.names=TRUE)) {
  tr <- tryCatch(read.tree(f), error=function(e) NULL)
  if (is.null(tr)) { message("UNPARSEABLE: ", basename(f)); next }
  if (inherits(tr,"multiPhylo")) tr <- tr[[1]]
  nm <- sub("\\.tree$","",basename(f)); nm <- sub("^tree_","",nm)
  nm <- sub("_(rmix|snp|adjusted|FINAL|v[0-9]|m40|[0-9]_01|[0-9]+_[0-9]+).*$","",nm)
  add(paste0("file:",basename(f)), nm, tr)
}
r <- do.call(rbind, rows)
write.table(r, "all_tips.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("sources=%d donors=%d tips=%d\n", length(unique(r$source)), length(unique(r$donor)), nrow(r)))
