#!/usr/bin/env bash
# Pilot: run discovery on the FIRST colony as a single LSF job with generous memory
# and a wall-clock timer, so you can (a) confirm the whole path works on real GRCh37
# WGS before launching 185 jobs, and (b) read the actual peak RSS + CPU time from the
# LSF report to right-size MEM/QUEUE for the array.
#
# Usage:  cluster/pilot.sh
# Then:   check logs/pilot.out for "Max Memory" and "CPU time"; set MEM for submit_discovery.sh
set -euo pipefail
cd "$(dirname "$0")/.."

FOFN="${FOFN:-bams.fofn}"
BAM="$(sed -n '1p' "$FOFN")"
ID="$(basename "$(dirname "$(dirname "$BAM")")")"
mkdir -p logs discovery

echo "pilot colony: $ID"
echo "  bam: $BAM"
bsub -J "ptpilot" -o logs/pilot.out -e logs/pilot.err \
     -n 1 -q "${QUEUE:-normal}" \
     -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" -M 16000 \
     "FOFN='$FOFN' bash cluster/discover_one.sh 1"

echo "submitted. when it finishes:"
echo "  grep -E 'Max Memory|CPU time|Successfully' logs/pilot.out"
# discovery output is FASTQ (@<contig>:<l>-<r>:LEFT|RIGHT:...), NOT the '>'-headed
# genotyping contract; count distinct breakpoint loci from the record headers.
echo "  zcat discovery/${ID}.txt.gz | grep -oE '^@[^:]+:[0-9]+-[0-9]+' | sort -u | wc -l   # candidate breakpoints"
