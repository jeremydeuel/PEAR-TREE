#!/usr/bin/env bash
# ============================================================================
# GRCh38 discovery ablation sweep — isolate WHICH of the 5 gating features that
# changed since the noslip baseline causes the ~30% common-MEI recall loss.
#
# Features (all `= true` in the shipped config.discovery.grch38; all absent in the
# noslip baseline the frozen 9x10 truth set was built on):
#   SLIP = slippage_filter        (SPEC-8,  mapped-side poly-A gate) [removal]
#   CLIP = clip_slippage_filter   (SPEC-8b, bilateral clip gate)     [removal]
#   COV  = coverage_mask          (SPEC-3,  5x local-median gate)    [removal, ratio]
#   DISC = discordant_anchor      (Feature A)                        [additive]
#   MATE = mate_anchor_rescue                                        [additive/relocate]
#
# 12-arm screening design (alloff + allon + 5 leave-one-out + 5 add-one-in).
# allon and loo_slip are already computed on THIS binary (insertions_grch38_speca8ba
# and insertions_grch38_noSPEC8), so this script runs the remaining 10 arms.
# Each arm = a 90-BAM discovery array + a dependent combine into its own contract.
#
# Self-contained: absolute paths, PATIENTS derived inline. Run on a farm22 head node.
# ============================================================================
set -euo pipefail

REPO=/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE-pd44579
BASE="$REPO/cluster/config.discovery.grch38"
FOFN=/nfs/users/nfs_j/jd43/catalogue/picked/all.bams.fofn
FOFNDIR=/nfs/users/nfs_j/jd43/catalogue/picked
STEM=mei9x10
THROTTLE=25          # per-array concurrency; 10 arrays => <=250 slots, LSF pends the rest
COMBINE_MEM=32000

cd "$REPO"
[ -s "$BASE" ] || { echo "no base config: $BASE" >&2; exit 1; }
[ -s "$FOFN" ] || { echo "no fofn: $FOFN" >&2; exit 1; }
mkdir -p logs/sweep

# anchored flag setter: matches "<key><sp>=<sp><word>" only; leaves *_multiplier,
# *_min_ref_run, *_min_reads, *_min_run etc. untouched.
set_flag() {
  local f=$1 key=$2 v=$3 val
  [ "$v" = 1 ] && val=true || val=false
  sed -E -i.bak "s/^(${key}[[:space:]]*=[[:space:]]*)[a-z]+/\1${val}/" "$f"
  rm -f "$f.bak"
}

# arm  SLIP CLIP COV DISC MATE   (allon 1 1 1 1 1 and loo_slip 0 1 1 1 1 already done)
ARMS="alloff:0:0:0:0:0
loo_clip:1:0:1:1:1
loo_cov:1:1:0:1:1
loo_disc:1:1:1:0:1
loo_mate:1:1:1:1:0
aoi_slip:1:0:0:0:0
aoi_clip:0:1:0:0:0
aoi_cov:0:0:1:0:0
aoi_disc:0:0:0:1:0
aoi_mate:0:0:0:0:1"

echo "=== generating + verifying 10 arm configs ==="
printf "%-9s %-5s %-5s %-5s %-5s %-5s\n" arm SLIP CLIP COV DISC MATE
for spec in $ARMS; do
  IFS=: read -r arm slip clip cov disc mate <<<"$spec"
  f="$REPO/cluster/config.discovery.grch38.$arm"
  cp "$BASE" "$f"
  set_flag "$f" slippage_filter      "$slip"
  set_flag "$f" clip_slippage_filter "$clip"
  set_flag "$f" coverage_mask        "$cov"
  set_flag "$f" discordant_anchor    "$disc"
  set_flag "$f" mate_anchor_rescue   "$mate"
  gv() { grep -E "^$1[[:space:]]*=" "$f" | head -1 | sed -E "s/.*=[[:space:]]*([a-z]+).*/\1/"; }
  printf "%-9s %-5s %-5s %-5s %-5s %-5s\n" "$arm" \
    "$(gv slippage_filter)" "$(gv clip_slippage_filter)" "$(gv coverage_mask)" \
    "$(gv discordant_anchor)" "$(gv mate_anchor_rescue)"
  # guardrail: the exon track splice_hallmark needs must still be present+readable
  grep -qE "^splice_hallmark[[:space:]]*=[[:space:]]*true" "$f" || { echo "FAIL $arm: splice_hallmark not true" >&2; exit 1; }
  ea=$(grep -E "^exon_annotation[[:space:]]*=" "$f" | sed -E "s/.*=[[:space:]]*//"); [ -r "$ea" ] || { echo "FAIL $arm: exon_annotation unreadable: $ea" >&2; exit 1; }
done

echo
echo "=== submitting discovery arrays (one per arm) ==="
for spec in $ARMS; do
  arm=${spec%%:*}
  OUTDIR="$REPO/discovery_grch38_$arm"
  CFG="$REPO/cluster/config.discovery.grch38.$arm"
  echo "--- $arm -> $OUTDIR"
  FOFN="$FOFN" OUTDIR="$OUTDIR" DISCOVER_CFG="$CFG" THROTTLE="$THROTTLE" \
    bash cluster/submit_discovery.sh
done

echo
echo "=== submitting dependent combine jobs ==="
# PATIENTS = every <stem>.bams.fofn in FOFNDIR except all.bams.fofn (what the assembly gate iterates)
for spec in $ARMS; do
  arm=${spec%%:*}
  DISCDIR="$REPO/discovery_grch38_$arm"
  OUTDIR="$REPO/insertions_grch38_$arm"
  DEP="ptdisc_discovery_grch38_$arm"     # job name submit_discovery.sh assigns to the array
  mkdir -p logs/sweep
  echo "--- combine $arm  (after ended($DEP)) -> $OUTDIR/$STEM.genotyping.txt.gz"
  bsub -J "combine_$arm" -w "ended($DEP)" \
       -n 16 -M "$COMBINE_MEM" \
       -R "select[mem>$COMBINE_MEM] rusage[mem=$COMBINE_MEM] span[hosts=1]" -q normal \
       -o "logs/sweep/combine_$arm.%J.log" -e "logs/sweep/combine_$arm.%J.err" \
       bash -c 'PATIENTS=$(for p in '"$FOFNDIR"'/*.bams.fofn; do b=$(basename "$p" .bams.fofn); [ "$b" = all ] || printf "%s " "$b"; done) \
                FOFNDIR='"$FOFNDIR"' DISCDIR='"$DISCDIR"' OUTDIR='"$OUTDIR"' STEM='"$STEM"' THREADS=16 \
                bash '"$REPO"'/cluster/combine_mei.sh'
done

echo
echo "=== submitted. watch with: bjobs -A ; bjobs -w | grep -E 'ptdisc_discovery_grch38_(alloff|loo_|aoi_)|combine_' ==="
echo "10 contracts will land at: $REPO/insertions_grch38_{alloff,loo_clip,loo_cov,loo_disc,loo_mate,aoi_slip,aoi_clip,aoi_cov,aoi_disc,aoi_mate}/$STEM.genotyping.txt.gz"
echo "(plus the 2 already done: insertions_grch38_speca8ba [allon] and insertions_grch38_noSPEC8 [loo_slip])"