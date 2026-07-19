#!/usr/bin/env bash
# Build ONE region-slice BAM + its exact full-BAM genome median. Called by build_slices.sh
# as an array task. Idempotent (skips work already done). Scheduler-agnostic.
#
# Usage: slice_one.sh <index> <full_fofn> <bed> <slicedir> <meddir> <disc_bin> <samtools_module>
set -euo pipefail

IDX="${1:?index}"; FOFN="${2:?fofn}"; BED="${3:?bed}"; SLICEDIR="${4:?slicedir}"
MEDDIR="${5:?meddir}"; BIN="${6:?bin}"; SAMTOOLS_MODULE="${7:?samtools module}"

module load "$SAMTOOLS_MODULE"

BAM="$(sed -n "${IDX}p" "$FOFN")"
[ -n "$BAM" ] || { echo "no BAM at line $IDX of $FOFN" >&2; exit 1; }
[ -s "$BAM" ] || { echo "missing BAM: $BAM" >&2; exit 1; }
ID="$(basename "$BAM")"; ID="${ID%.bam}"; ID="${ID%.cram}"; ID="${ID%.sample.dupmarked}"

SLICE="$SLICEDIR/${ID}.slice.bam"
MED="$MEDDIR/${ID}.median"

# 1) region slice (index-based multi-region iterator). Needs a coordinate index on the
#    full BAM; build one if absent (one-time). Reads only the ~1.1% of the BAM in the BED.
if [ ! -s "$SLICE" ] || [ ! -s "$SLICE.bai" ]; then
    if [ ! -e "$BAM.bai" ] && [ ! -e "$BAM.csi" ] && [ ! -e "${BAM%.bam}.bai" ]; then
        echo "[$IDX] $ID: no index on full BAM, building one"
        samtools index -@ 4 "$BAM"
    fi
    TMP="$SLICE.tmp.$$"
    echo "[$IDX] $ID: slicing -> $SLICE"
    samtools view -b -M -L "$BED" -@ 4 "$BAM" -o "$TMP"
    samtools index -@ 4 "$TMP"
    mv -f "$TMP" "$SLICE"; mv -f "$TMP.bai" "$SLICE.bai"
else
    echo "[$IDX] $ID: slice exists, skipping"
fi

# 2) exact genome median from the FULL BAM (one sequential pass; no index needed).
#    Store the bare number so discover_one can `cat` it straight into the env var.
if [ ! -s "$MED" ]; then
    echo "[$IDX] $ID: computing full-BAM coverage median"
    "$BIN" --step coverage-median --bam "$BAM" -@ 4 | awk -F'\t' '{print $2}' > "$MED.tmp.$$"
    [ -s "$MED.tmp.$$" ] || { echo "[$IDX] $ID: empty median!" >&2; rm -f "$MED.tmp.$$"; exit 1; }
    mv -f "$MED.tmp.$$" "$MED"
else
    echo "[$IDX] $ID: median exists, skipping"
fi

echo "[$IDX] $ID: done  slice=$SLICE  median=$(cat "$MED")"
