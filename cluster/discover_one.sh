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
# Derive the colony id from the FILENAME, not the directory. Both layouts we read from
# end in <id>.sample.dupmarked.bam, but their parent dirs differ:
#   staging  : <id>/mapped_sample/<id>.sample.dupmarked.bam       -> dir(dir()) = <id>
#   nst_links: <project>/<id>/<id>.sample.dupmarked.bam           -> dir(dir()) = <project>  !!
# The old dir-based rule silently collapsed every nst_links colony of a patient onto the
# PROJECT number, so all but the first "skipped, output exists" and you got 1 file per
# project instead of 1 per colony. Filename-based is correct for both.
ID="$(basename "$BAM")"; ID="${ID%.bam}"; ID="${ID%.cram}"; ID="${ID%.sample.dupmarked}"
[ -n "$ID" ] || { echo "could not derive id from $BAM" >&2; exit 1; }
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
# The binary writes sidecars next to --out ({out}.stats.json, and — when enabled —
# {out}.splice.tsv / {out}.hallmarks.tsv). Carry them through the atomic rename so
# they end up as $OUT.<ext>, not orphaned under the .tmp name.
for ext in stats.json splice.tsv hallmarks.tsv; do
    [ -e "$TMP.$ext" ] && mv -f "$TMP.$ext" "$OUT.$ext"
done
echo "[$IDX] $ID: done -> $OUT"
