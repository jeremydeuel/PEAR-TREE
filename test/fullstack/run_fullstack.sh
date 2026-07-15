#!/bin/bash
# Full-stack test-case builder: hs1 donor (+ implanted MEIs) -> wgsim reads -> library
# artefacts -> bwa-mem to GRCh38 -> markdup'd BAM, plus truth lifted into hg38 coords.
# Deliberately maps hs1-derived reads to hg38 so assembly discordance yields the false
# positives PEAR-TREE must filter, alongside the implanted true positives.
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
HS1="${HS1:-$HOME/Downloads/hs1.fa}"
HG38="${HG38:-$HOME/Downloads/hg38.fa}"
OUT="${OUT:-$HOME/Downloads/fullstack}"
PY="${PY:-/Users/jeremy/Documents/PEAR-TREE/venv/bin/python}"
CHAIN="${CHAIN:-$HOME/Downloads/hs1ToHg38.over.chain.gz}"  # hs1->hg38 liftOver chain for chain-based truth lift
THREADS="${THREADS:-8}"
DEPTH="${DEPTH:-30}"
mkdir -p "$OUT"
ts(){ date "+%H:%M:%S"; }

echo "[$(ts)] 1. build donor + implant MEIs"
"$PY" "$DIR/build_donor.py" --hs1 "$HS1" \
    --out-donor "$OUT/donor.fa" --out-truth "$OUT/truth_hs1.tsv" --out-flanks "$OUT/flanks.fa"
samtools faidx "$OUT/donor.fa"

BP=$(awk '{s+=$2} END{print s}' "$OUT/donor.fa.fai")
N=$(( BP * DEPTH / 300 ))
echo "[$(ts)] 2. wgsim ~${DEPTH}x over ${BP} bp -> $N pairs"
wgsim -N "$N" -1 150 -2 150 -d 320 -s 40 -e 0.005 -r 0 -R 0 -X 0 \
    "$OUT/donor.fa" "$OUT/r1.fq" "$OUT/r2.fq" > "$OUT/wgsim.log" 2>&1

echo "[$(ts)] 3. inject library artefacts"
"$PY" "$DIR/inject_artefacts.py" --in1 "$OUT/r1.fq" --in2 "$OUT/r2.fq" \
    --out1 "$OUT/r1a.fq" --out2 "$OUT/r2a.fq"

echo "[$(ts)] 4. bwa-mem -> hg38, fixmate + sort + markdup"
bwa mem -t "$THREADS" -R '@RG\tID:hs1sim\tSM:hs1sim\tPL:ILLUMINA' \
    "$HG38" "$OUT/r1a.fq" "$OUT/r2a.fq" 2>"$OUT/bwa.log" \
  | samtools fixmate -m -u -@ "$THREADS" - - \
  | samtools sort -u -@ "$THREADS" - \
  | samtools markdup -@ "$THREADS" - "$OUT/full_stack.bam"
samtools index "$OUT/full_stack.bam"

echo "[$(ts)] 5. lift truth into hg38 coords (chain-lift insertion point; map flanks as fallback)"
bwa mem -t "$THREADS" "$HG38" "$OUT/flanks.fa" 2>>"$OUT/bwa.log" \
  | samtools sort -o "$OUT/flanks.bam" -
samtools index "$OUT/flanks.bam"
CHAIN_ARG=(); [ -f "$CHAIN" ] && CHAIN_ARG=(--chain "$CHAIN") || echo "  (chain $CHAIN not found; flank-only lift)"
"$PY" "$DIR/lift_truth.py" --flanks-bam "$OUT/flanks.bam" \
    --truth-hs1 "$OUT/truth_hs1.tsv" --out "$OUT/truth_hg38.tsv" "${CHAIN_ARG[@]}"

echo "[$(ts)] cleanup fastq intermediates"
rm -f "$OUT/r1.fq" "$OUT/r2.fq" "$OUT/r1a.fq" "$OUT/r2a.fq"
echo "[$(ts)] DONE"
echo "  BAM   : $OUT/full_stack.bam"
echo "  truth : $OUT/truth_hg38.tsv  (hg38 coords; status=scoreable rows)"
samtools flagstat "$OUT/full_stack.bam" | head -5
