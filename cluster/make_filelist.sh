#!/usr/bin/env bash
# Build the discovery file list (fofn) for one donor's colonies.
#
# Picks the REAL per-colony BAMs (<id>/mapped_sample/<id>.sample.dupmarked.bam) and
# excludes the tmpExportData/progress/ staging duplicates. Writes one absolute path
# per line to bams.fofn.
#
# Usage:  cluster/make_filelist.sh <data-dir> [out.fofn]
#   e.g.  cluster/make_filelist.sh /lustre/scratch126/casm/staging/team273/jd43/2178
set -euo pipefail

DATADIR="${1:?usage: make_filelist.sh <data-dir> [out.fofn]}"
OUT="${2:-bams.fofn}"

find "$DATADIR" -type f \
     -path '*/mapped_sample/*.sample.dupmarked.bam' \
     -not -path '*/tmpExportData/*' \
  | sort > "$OUT"

n=$(wc -l < "$OUT" | tr -d ' ')
echo "wrote $n colony BAMs -> $OUT"
echo "sample colony IDs:"
head -3 "$OUT" | while read -r b; do basename "$(dirname "$(dirname "$b")")"; done
echo "  ..."
