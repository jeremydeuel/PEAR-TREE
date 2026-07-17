#!/usr/bin/env Rscript
# Pick the 10x10 benchmark: 10 blood/HSPC donors x 10 colonies, GRCh38.
#
#   Rscript cluster/trees/pick_clades.R > cluster/trees/picked.tsv
#
# DESIGN (hybrid, chosen by Jeremy 2026-07-17): per donor, 5 CLADE + 5 SPREAD.
#
#   CLADE  5 colonies under the deepest internal node with >=5 tips. They share that node's
#          entire root-to-node path, so every SOMATIC MEI on it must appear in all 5. The
#          path length in mutations IS the somatic true-positive signal available.
#   SPREAD 5 colonies chosen by greedy farthest-point sampling (max-min cophenetic distance),
#          seeded from the tip furthest from the clade. These coalesce with the clade only at
#          the root, so they share GERMLINE MEIs with it but essentially no somatic ones.
#
# WHY BOTH. Germline MEIs are shared by all 10 colonies of a donor by construction -> the
# primary TP set, and it needs no clade structure at all. Somatic MEIs live only on shared
# internal branches -> they need the clade. The spread half also keeps FP sampling honest:
# 10 tightly-related colonies would see the same artefacts together and could look like a
# consistent (i.e. "true") call.
#
# WHAT THE TREE SHAPE FORCES ON US (measured, do not re-litigate):
# Normal HSPC trees are STAR-SHAPED -- HSPCs coalesce in early embryogenesis, then each
# colony accrues a long private branch. So a clade of >=5 tips only exists near the root
# unless the donor had a CLONAL EXPANSION. Measured 5-tip clade depths:
#   PX001 3037 | KX008 1060 | KX004 1030 | KX003 933 | BMH1_TG 613 | KX002 174
#   SX001   43 | PX003   31 | KX001   28 | CB002   25
# The first five are the expansion/chemo donors and carry real somatic signal. The last four
# carry ~none: they are GERMLINE + FP donors, and that is biology, not bad picking.
# CB002 is cord blood -- its ENTIRE tree is ~47 mutations deep. It can never contribute a
# somatic TP. Do not "fix" this by loosening the clade criterion.
#
# Depth units are the tree's own (SNV counts here): comparable WITHIN a donor, NOT across
# donors of different age/exposure.
suppressMessages(library(ape))

BASE <- "/Users/jeremy/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents - STEM_Green_Lab/HNRNPA1"
N_CLADE  <- as.integer(Sys.getenv("N_CLADE",  "5"))
N_SPREAD <- as.integer(Sys.getenv("N_SPREAD", "5"))

# donor -> alias. All blood/HSPC, all GRCh38, all 151bp PE.
# Tissue for KX004/KX008/SX001 (blood) and CB002 (cord blood) CONFIRMED BY JEREMY 2026-07-17.
# Not inferred from the codename -- that inference is unsafe (it is what made me call `_lo`
# "clonal", which was wrong).
#
# EXCLUDED ON PURPOSE -- both exclusions are load-bearing, do not "restore" them:
#
#  PD44579/PX002. Blood, GRCh38, 174 tips: qualifies mechanically. But every filter we have
#    (min_dispersion, min_wild-types, the reference-context gate, the carrier-count band) was
#    developed on it. Scoring it measures fit to PD44579, not detection. POSITIVE CONTROL ONLY.
#
#  AX001/BMH1_TG. Its 361 tree tips are PLATE/WELL ids (BMH1_TG001_P31_F08), and NOTHING we
#    hold maps those to sample ids: zero catalogued samples begin with "BMH1_TG", and the
#    sibling donor PD43976 -- whose 10-tip tree in the chemo paper is made ENTIRELY of
#    BMH1_TG001_* tips, so it is likely the same donor -- has samples named PD43976aaa/xq.
#    An earlier claim that BMH1_TG001_* samples live in projects 2445/2318 is UNVERIFIED and
#    contradicted by the catalogue; do not act on it without a plate/well->sample key.
#    Re-add AX001 as a 10th donor once that key exists.
#
# So this is a 9x10 (90 colonies). Cohort gates in combine_genotypes assume >=2 donors, not 10;
# nothing downstream depends on the donor count being exactly 10.
WANT <- c(PD40521="KX001", PD40667="KX002", PD43974="KX003",
          PD47703="PX001", PD50307="PX003", PD45534="KX004", PD48402="KX008",
          PD45517="CB002", PD41048="SX001")

trees <- list()
x <- readRDS(file.path(BASE, "spar_2ndrev/human/trees.rds"))
for (p in names(x)) if (inherits(x[[p]], "phylo")) trees[[p]] <- x[[p]]
y <- readRDS(file.path(BASE, "JAK2_HNRNPA1_paper/Repo_080602025/Fig1_trees/mpn_trees/PD_MPN_CUT.RDS"))
for (p in names(y)) if (inherits(y[[p]]$tree, "phylo")) trees[[p]] <- y[[p]]$tree

# deepest internal node with >= n descendant tips; return its n tips nearest the node
deepest_clade <- function(tr, n, d) {
  ntip <- length(tr$tip.label)
  desc <- prop.part(tr)
  best <- NULL
  for (i in seq_along(desc)) {
    tips <- desc[[i]]
    if (length(tips) < n) next
    node <- ntip + i
    cand <- list(node = node, depth = d[node], ntips = length(tips), tips = tips)
    if (is.null(best) || cand$depth > best$depth ||
        (cand$depth == best$depth && cand$ntips < best$ntips)) best <- cand
  }
  if (is.null(best)) return(NULL)
  sel <- best$tips[order(d[best$tips] - best$depth)][seq_len(n)]
  list(node = best$node, depth = best$depth, size = best$ntips, sel = sel)
}

# greedy farthest-point: repeatedly take the tip maximising its MINIMUM distance to all
# already-chosen tips. Seeded with the clade, so pick 1 is the tip furthest from the clade.
spread_pick <- function(cop, avoid, pool, n) {
  chosen <- integer(0)
  ref <- avoid
  for (k in seq_len(n)) {
    if (!length(pool)) break
    mind <- apply(cop[pool, ref, drop = FALSE], 1, min)
    take <- pool[which.max(mind)]
    chosen <- c(chosen, take)
    ref <- c(ref, take)
    pool <- setdiff(pool, take)
  }
  chosen
}

cat("donor\talias\trole\tclade_node\tclade_depth\tclade_size\ttree_tips\ttip\ttip_depth\n")
for (donor in names(WANT)) {
  alias <- WANT[[donor]]
  key <- if (!is.null(trees[[donor]])) donor else if (!is.null(trees[[alias]])) alias else NA
  if (is.na(key)) { message("NO TREE: ", donor, " (", alias, ")"); next }
  tr <- trees[[key]]
  if (is.null(tr$edge.length)) { message("NO BRANCH LENGTHS: ", donor); next }

  d  <- node.depth.edgelength(tr)
  cl <- deepest_clade(tr, N_CLADE, d)
  if (is.null(cl)) { message("NO CLADE >= ", N_CLADE, ": ", donor); next }

  cop  <- cophenetic.phylo(tr)
  pool <- setdiff(seq_along(tr$tip.label), cl$tips)   # exclude the WHOLE clade, not just the 5
  sp   <- spread_pick(cop, cl$sel, pool, N_SPREAD)

  minsep <- min(cop[sp, cl$sel])
  message(sprintf("%-8s %-8s tips=%4d | clade d=%7.1f size=%3d n=%d | spread n=%d min_sep_from_clade=%7.1f",
                  donor, alias, length(tr$tip.label), cl$depth, cl$size, length(cl$sel),
                  length(sp), minsep))

  emit <- function(idx, role) for (i in idx)
    cat(sprintf("%s\t%s\t%s\t%d\t%.1f\t%d\t%d\t%s\t%.1f\n", donor, alias, role, cl$node,
                cl$depth, cl$size, length(tr$tip.label), tr$tip.label[i], d[i]))
  emit(cl$sel, "clade")
  emit(sp,     "spread")
}
