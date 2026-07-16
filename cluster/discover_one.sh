#!/usr/bin/env bash
# Run Rust discovery on ONE colony BAM. Scheduler-agnostic: the submit wrapper
# passes the 1-based line number of bams.fofn to process. Safe to re-run: skips
# colonies whose output already exists, and writes atomically (tmp -> final) so a
# killed job never leaves a half-written .txt.gz that looks complete.
#
# Usage:            cluster/discover_one.sh <index>
# Env (all optional):
#   FOFN         file list                (default: bams.fofn)
#   OUTDIR       output dir               (default: discovery)
#   DISCOVER_BIN discovery binary         (default: rust/peartree-discovery/target/release/peartree-discovery)
#   DISCOVER_CFG rust discovery config    (default: cluster/config.discovery.grch37)
set -euo pipefail

IDX="${1:?usage: discover_one.sh <1-based-index-into-fofn>}"
FOFN="${FOFN:-bams.fofn}"
OUTDIR="${OUTDIR:-discovery}"
BIN="${DISCOVER_BIN:-rust/peartree-discovery/target/release/peartree-discovery}"
CFG="${DISCOVER_CFG:-cluster/config.discovery.grch37}"

BAM="$(sed -n "${IDX}p" "$FOFN")"
[ -n "$BAM" ] || { echo "no BAM at line $IDX of $FOFN" >&2; exit 1; }
# <id>/mapped_sample/<id>.sample.dupmarked.bam  ->  <id>
ID="$(basename "$(dirname "$(dirname "$BAM")")")"
OUT="$OUTDIR/${ID}.txt.gz"

mkdir -p "$OUTDIR"
if [ -s "$OUT" ]; then
    echo "[$IDX] $ID: output exists, skipping ($OUT)"
    exit 0
fi

TMP="$OUT.tmp.$$"
echo "[$IDX] $ID: discovering $BAM"
"$BIN" --step discover --bam "$BAM" --out "$TMP" --threads 1 --config "$CFG"
mv -f "$TMP" "$OUT"
echo "[$IDX] $ID: done -> $OUT"
