#!/bin/bash
# Profile satellite/rRNA density per 1 Mb bin across hs1 -> sat_bins.tsv, the input
# gen_windows.py uses to pick the centromere/telomere FP compartment. One pass over the
# hs1 RepeatMasker .out.gz (columns: $5=contig $6=begin $7=end $11=class/family).
set -euo pipefail
RMSK="${1:?usage: profile_satellite.sh <hs1.repeatMasker.out.gz> <out/sat_bins.tsv>}"
OUT="${2:?out path}"
gzip -dc "$RMSK" | awk 'NR>3 && $11 ~ /Satellite|rRNA|acro|centr/ {
    c=$5; b=$6+0; e=$7+0; bin=int(b/1000000); sat[c"\t"bin]+=(e-b)
} END { for(k in sat) print k, sat[k] }' | sort -k1,1 -k2,2n > "$OUT"
echo "wrote $(wc -l < "$OUT") bins to $OUT"
