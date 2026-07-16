#!/usr/bin/env bash
# Submit the discovery step as an LSF job array on farm22: one task per colony,
# 1 core each, throttled so we never storm shared Lustre with 185 full-BAM scans.
#
# Usage:  cluster/submit_discovery.sh
# Env (all optional):
#   FOFN      file list         (default: bams.fofn)
#   THROTTLE  max concurrent    (default: 50)
#   MEM       MB per task       (default: 16000 -- MEASURED: discovery peaked at 12.2 GB
#                                on PD44579 30x WGS. Do NOT lower to ~8 GB: every task
#                                dies with TERM_MEMLIMIT. Re-check with pilot.sh on new data.)
#   QUEUE     LSF queue         (default: normal; ~16-24 min/colony single-threaded)
#   GROUP     LSF fairshare group (-G), if your setup requires one
#
# SLURM equivalent (if farm22 ever moves to SLURM):
#   sbatch --array=1-$N%$THROTTLE --cpus-per-task=1 --mem=${MEM} \
#          --wrap 'cluster/discover_one.sh $SLURM_ARRAY_TASK_ID'
set -euo pipefail
cd "$(dirname "$0")/.."

FOFN="${FOFN:-bams.fofn}"
THROTTLE="${THROTTLE:-50}"
MEM="${MEM:-16000}"   # measured: 12.2 GB peak on PD44579 (8000 => TERM_MEMLIMIT)
QUEUE="${QUEUE:-normal}"
N="$(wc -l < "$FOFN" | tr -d ' ')"
[ "$N" -gt 0 ] || { echo "empty $FOFN" >&2; exit 1; }
mkdir -p logs discovery

GROUP_ARG=()
[ -n "${GROUP:-}" ] && GROUP_ARG=(-G "$GROUP")

echo "submitting discovery array: $N colonies, <=$THROTTLE concurrent, ${MEM}MB, queue=$QUEUE"
bsub \
    -J "ptdisc[1-${N}]%${THROTTLE}" \
    -o "logs/disc.%I.out" -e "logs/disc.%I.err" \
    -n 1 -q "$QUEUE" "${GROUP_ARG[@]}" \
    -R "select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" -M "${MEM}" \
    "FOFN='$FOFN' bash cluster/discover_one.sh \$LSB_JOBINDEX"

echo "watch with: bjobs -A ; tail -f logs/disc.1.out"
echo "when done: ls discovery/*.txt.gz | wc -l   # expect $N"
