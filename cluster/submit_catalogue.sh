#!/usr/bin/env bash
# Submit PHASE 2 of the sample catalogue: header reads as a throttled LSF array.
#
#   bash cluster/submit_catalogue.sh          # after cluster/catalogue_scan.sh
#
# THROTTLE matters more than usual here. Every task does ~4 small reads per BAM over shared
# NFS (nst_links), and the whole point of the catalogue is that there are tens of thousands
# of them. Unthrottled this is a denial-of-service against a filesystem the whole institute
# uses. 20 concurrent tasks is deliberately modest — this is a background job, not a
# deadline. Raise only if you know the NFS server is idle.
#
# Env: CHUNK (rows/task, default 500), THROTTLE (default 20), MEM (default 2000), QUEUE
set -euo pipefail
cd "$(dirname "$0")/.."

OUT=${OUT:-$HOME/catalogue}
CHUNK=${CHUNK:-500}
THROTTLE=${THROTTLE:-20}
MEM=${MEM:-2000}          # header reads are tiny; this is samtools + bash overhead only
QUEUE=${QUEUE:-normal}
MAN="$OUT/manifest.tsv"

[ -s "$MAN" ] || { echo "no manifest: $MAN — run cluster/catalogue_scan.sh first" >&2; exit 1; }
ROWS=$(( $(wc -l < "$MAN") - 1 ))
[ "$ROWS" -gt 0 ] || { echo "manifest has no rows" >&2; exit 1; }
N=$(( (ROWS + CHUNK - 1) / CHUNK ))
mkdir -p logs "$OUT/parts"

echo "  MANIFEST = $MAN  ($ROWS samples)"
echo "  CHUNK    = $CHUNK rows/task  ->  $N tasks"
echo "  THROTTLE = $THROTTLE concurrent"
echo "  OUT      = $OUT/parts  ($(ls "$OUT/parts" 2>/dev/null | wc -l | tr -d ' ') parts already present)"

GROUP_ARG=(); [ -n "${GROUP:-}" ] && GROUP_ARG=(-G "$GROUP")

bsub \
    -J "ptcat[1-${N}]%${THROTTLE}" \
    -o "logs/cat.%I.out" -e "logs/cat.%I.err" \
    -n 1 -q "$QUEUE" "${GROUP_ARG[@]}" \
    -R "select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" -M "${MEM}" \
    "OUT='$OUT' CHUNK='$CHUNK' bash cluster/catalogue_headers.sh \$LSB_JOBINDEX"

echo
echo "watch:     bjobs -A"
echo "when done: bash cluster/catalogue_merge.sh"
