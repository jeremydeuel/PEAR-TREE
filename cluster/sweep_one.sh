#!/usr/bin/env bash
# ============================================================================
# Run ONE discovery config against the 9x10 region slices and score it, fully on
# the farm. Chains: slice-discovery (median pinned) -> combine -> adjudicate, each
# LSF-dependent on the last, and appends a one-line verdict to verdicts.tsv.
#
# Because the slices are ~1% of each BAM, the whole chain finishes in minutes, so
# you can fire many configs back-to-back to fine-tune thresholds.
#
#   bash cluster/sweep_one.sh <arm-label> <config-file>
# e.g.
#   bash cluster/sweep_one.sh cov3 cluster/config.discovery.grch38.cov3
#
# Prereq: cluster/build_slices.sh has finished (slices + medians present).
# ============================================================================
set -euo pipefail

ARM="${1:?usage: sweep_one.sh <arm-label> <config-file>}"
CFG="${2:?usage: sweep_one.sh <arm-label> <config-file>}"

REPO=/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE-pd44579
SLICEFOFN="$REPO/slices/slices.fofn"
MEDDIR="$REPO/slices/medians"
FOFNDIR=/nfs/users/nfs_j/jd43/catalogue/picked
TRUTH="$REPO/analysis/mei9x10/frozen_truth.tsv"
STEM=mei9x10
THREADS=16
VERDICTS="$REPO/logs/sweep_slice/verdicts.tsv"

cd "$REPO"
[ -s "$SLICEFOFN" ] || { echo "no slice fofn: $SLICEFOFN  (run cluster/build_slices.sh)" >&2; exit 1; }
[ -d "$MEDDIR" ]    || { echo "no median dir: $MEDDIR" >&2; exit 1; }
[ -s "$CFG" ]       || { echo "no config: $CFG" >&2; exit 1; }
[ -s "$TRUTH" ]     || { echo "no truth set: $TRUTH" >&2; exit 1; }
NS=$(wc -l < "$SLICEFOFN"); NS=${NS// /}
NM=$(ls "$MEDDIR"/*.median 2>/dev/null | wc -l); NM=${NM// /}
[ "$NM" -ge "$NS" ] || { echo "only $NM/$NS medians present — build_slices not finished?" >&2; exit 1; }
mkdir -p logs/sweep_slice

DISCDIR="$REPO/discovery_slice_$ARM"
INSDIR="$REPO/insertions_slice_$ARM"
DISCJOB="ptdisc_discovery_slice_$ARM"       # submit_discovery derives this from basename(OUTDIR)
CMBJOB="combine_slice_$ARM"
ADJJOB="adj_slice_$ARM"

# PATIENTS = every <patient>.bams.fofn in FOFNDIR except all.bams.fofn
PATIENTS=$(for p in "$FOFNDIR"/*.bams.fofn; do b=$(basename "$p" .bams.fofn); [ "$b" = all ] || printf "%s " "$b"; done)

echo "=== [$ARM] slice-discovery ($NS colonies, median-pinned) -> $DISCDIR ==="
MEDIAN_DIR="$MEDDIR" FOFN="$SLICEFOFN" OUTDIR="$DISCDIR" DISCOVER_CFG="$CFG" THROTTLE="$NS" \
  bash cluster/submit_discovery.sh

echo "=== [$ARM] combine (after ended($DISCJOB)) -> $INSDIR/$STEM.genotyping.txt.gz ==="
bsub -J "$CMBJOB" -w "ended($DISCJOB)" -n "$THREADS" -M 32000 \
     -R "select[mem>32000] rusage[mem=32000] span[hosts=1]" -q normal \
     -o "logs/sweep_slice/$CMBJOB.%J.log" -e "logs/sweep_slice/$CMBJOB.%J.err" \
     "PATIENTS='$PATIENTS' FOFNDIR='$FOFNDIR' DISCDIR='$DISCDIR' OUTDIR='$INSDIR' STEM='$STEM' THREADS='$THREADS' bash '$REPO/cluster/combine_mei.sh'"

echo "=== [$ARM] adjudicate (after ended($CMBJOB)) -> $VERDICTS ==="
bsub -J "$ADJJOB" -w "ended($CMBJOB)" -n 1 -M 4000 \
     -R "select[mem>4000] rusage[mem=4000]" -q normal \
     -o "logs/sweep_slice/$ADJJOB.%J.log" -e "logs/sweep_slice/$ADJJOB.%J.err" \
     "python3 '$REPO/cluster/adjudicate_contract.py' '$INSDIR/$STEM.genotyping.txt.gz' '$TRUTH' '$ARM' | tee -a '$VERDICTS'"

echo
echo "verdict will land in: $VERDICTS  (watch: tail -f $VERDICTS)"
