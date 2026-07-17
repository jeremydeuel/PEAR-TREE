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
OUTDIR="${OUTDIR:-discovery}"
THROTTLE="${THROTTLE:-50}"
# MEASURED, twice — do not trim this back:
#   PD44579 30x WGS, 174 colonies : peak 12.2 GB  (8000 => every task TERM_MEMLIMIT)
#   PD41048b_lo0015 (proj 2073)   : peak 16.28 GB => TERM_MEMLIMIT at -M 16000, killed
#                                   at 811 s with Delta Memory -284 MB.
# Peak varies by cohort, not just by depth, and 16000 is NOT enough headroom: it lost
# 1/50 tasks on the 10x10 benchmark. 24000 covers the observed spread.
MEM="${MEM:-24000}"
QUEUE="${QUEUE:-normal}"
DISCOVER_CFG="${DISCOVER_CFG:-cluster/config.discovery.grch37}"
[ -s "$FOFN" ] || { echo "no such fofn: $FOFN (run cluster/build_fofn.sh?)" >&2; exit 1; }
N="$(wc -l < "$FOFN" | tr -d ' ')"
[ "$N" -gt 0 ] || { echo "empty $FOFN" >&2; exit 1; }
mkdir -p logs "$OUTDIR"

# Say out loud what we are about to do. Defaults here are PD44579's (bams.fofn ->
# discovery/), and if the caller's env does not reach this script those defaults submit a
# no-op array over an already-complete run: every task prints "output exists, skipping",
# the array vanishes in seconds and nothing new appears. That failure is silent, so print
# the resolved values and a sample of the input before submitting.
echo "  FOFN     = $FOFN  ($N BAMs)"
echo "  OUTDIR   = $OUTDIR  ($(ls "$OUTDIR" 2>/dev/null | wc -l | tr -d ' ') files already present)"
echo "  CONFIG   = $DISCOVER_CFG"
echo "  first    = $(head -1 "$FOFN")"

# PREFLIGHT: does the config's contig_allowlist actually name contigs in these BAMs?
#
# This is the single most expensive silent failure available to us. The allowlist is a set of
# NAMES; GRCh37/hs37d5 calls them 1..Y and GRCh38 calls them chr1..chrY. Point the grch37
# config at GRCh38 BAMs (the default DISCOVER_CFG below is grch37, so merely forgetting the
# env var does it) and NOTHING matches: every task runs to completion, exits 0, and writes an
# EMPTY output. "No reads on the allowlisted contigs" is indistinguishable from "no evidence
# found". 90 empty files that look like a clean run, ~24 min/colony of farm time, and the
# error only surfaces later as an empty contract.
#
# The Rust binary reports "contig allowlist: N contigs" but never intersects it with the BAM
# header, so it cannot catch this either. Do it here, once, before spending the array.
#
# Header read only (samtools view -H) -- allowed on any path under Sanger policy, and these
# are staged Lustre BAMs anyway.
FIRST_BAM="$(head -1 "$FOFN")"
if [ -s "$FIRST_BAM" ]; then
    module load samtools-1.19/python-3.12.0 2>/dev/null || true
    if command -v samtools >/dev/null 2>&1; then
        allow="$(awk -F'[[:space:]]*=[[:space:]]*' '/^[[:space:]]*contig_allowlist[[:space:]]*=/{print $2}' "$DISCOVER_CFG" | tr -d ' ')"
        if [ -n "$allow" ]; then
            hdr_sq="$(samtools view -H "$FIRST_BAM" 2>/dev/null | awk '$1=="@SQ"{for(i=2;i<=NF;i++) if($i ~ /^SN:/){sub(/^SN:/,"",$i); print $i}}')"
            n_allow=$(printf '%s' "$allow" | tr ',' '\n' | grep -c .)
            n_hit=$(printf '%s\n' "$allow" | tr ',' '\n' | grep -Fxc -f <(printf '%s\n' "$hdr_sq") - 2>/dev/null || true)
            n_hit=${n_hit:-0}
            echo "  allowlist: $n_hit/$n_allow allowlisted contigs present in $(basename "$FIRST_BAM")"
            if [ "$n_hit" -eq 0 ]; then
                echo >&2
                echo >&2 "REFUSING TO SUBMIT: none of the $n_allow allowlisted contigs exist in the BAM."
                echo >&2 "  config    : $DISCOVER_CFG"
                echo >&2 "  allowlist : $(printf '%s' "$allow" | cut -c1-60)..."
                echo >&2 "  BAM @SQ   : $(printf '%s\n' "$hdr_sq" | head -4 | paste -sd, -)..."
                echo >&2 "Every task would produce an EMPTY output and exit 0."
                echo >&2 "This is almost certainly an assembly/naming mismatch:"
                echo >&2 "  GRCh38 BAMs (chr-prefixed) -> DISCOVER_CFG=cluster/config.discovery.grch38"
                echo >&2 "  hs37d5 BAMs (numeric)      -> DISCOVER_CFG=cluster/config.discovery.grch37"
                exit 1
            elif [ "$n_hit" -lt "$n_allow" ]; then
                echo >&2 "WARNING: only $n_hit of $n_allow allowlisted contigs are in the BAM —"
                echo >&2 "         partial match. Check the config matches this cohort's assembly."
            fi
        fi
    else
        echo "  allowlist: NOT CHECKED (samtools not on PATH) — verify the config matches the BAM assembly" >&2
    fi
fi

GROUP_ARG=()
[ -n "${GROUP:-}" ] && GROUP_ARG=(-G "$GROUP")

echo "submitting discovery array: $N colonies, <=$THROTTLE concurrent, ${MEM}MB, queue=$QUEUE"
bsub \
    -J "ptdisc[1-${N}]%${THROTTLE}" \
    -o "logs/disc.%I.out" -e "logs/disc.%I.err" \
    -n 1 -q "$QUEUE" "${GROUP_ARG[@]}" \
    -R "select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" -M "${MEM}" \
    "FOFN='$FOFN' OUTDIR='$OUTDIR' DISCOVER_CFG='$DISCOVER_CFG' bash cluster/discover_one.sh \$LSB_JOBINDEX"

echo "watch with: bjobs -A ; tail -f logs/disc.1.out"
echo "when done: ls discovery/*.txt.gz | wc -l   # expect $N"
