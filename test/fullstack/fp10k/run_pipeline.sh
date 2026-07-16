#!/bin/bash
# Run the FULL PEAR-TREE pipeline on the FP-stress BAM:
#   combine_insertions (step 2) -> genotype -> combine_genotypes, then score the
#   post-combine calls against truth to see how much of the ~21k discovery FP the
#   specificity chain removes (and how many true insertions survive).
# Assumes run_fp10k.sh already produced $OUT/{fp_stack.bam,discovery*.txt.gz,truth_hg38.tsv}.
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../../.." && pwd)"
SCALE="$REPO/test/fullstack/scale10k"
PY="${PY:-$REPO/venv/bin/python}"
HG38="${HG38:-$HOME/Downloads/hg38.fa}"
HS1_BT2="${HS1_BT2:-$HOME/Downloads/hs1}"
CHAIN="${CHAIN:-$HOME/Downloads/hs1ToHg38.over.chain.gz}"
SAMTOOLS="${SAMTOOLS:-$(command -v samtools)}"
BOWTIE2="${BOWTIE2:-$(command -v bowtie2)}"
OUT="${OUT:-$HOME/Downloads/fullstack_fp10k}"
THREADS="${THREADS:-8}"
DISC="${DISC:-discovery.txt.gz}"          # which discovery calls feed combine (recommended cfg)
ts(){ date "+%H:%M:%S"; }
cd "$OUT"

echo "[$(ts)] 0. hg38 primary .2bit"
[ -f hg38.primary.2bit ] || "$PY" "$SCALE/make_2bit.py" "$HG38" hg38.primary.2bit

echo "[$(ts)] 1. materialise combine config.py"
sed -e "s#@@HG38_2BIT@@#$OUT/hg38.primary.2bit#" -e "s#@@SAMTOOLS@@#$SAMTOOLS#" \
    -e "s#@@BOWTIE2@@#$BOWTIE2#" -e "s#@@HS1_BT2@@#$HS1_BT2#" -e "s#@@HS1_HG38_CHAIN@@#$CHAIN#" \
    "$SCALE/config.local.py.template" > config.py

echo "[$(ts)] 2. combine_insertions (step 2) on $DISC"
"$PY" "$SCALE/run_combine.py" "$DISC" step2 "$THREADS" "$OUT" "$REPO/src" >combine.log 2>&1
tail -3 combine.log
ls -la step2.combined.txt.gz step2.genotyping.txt.gz 2>&1 | awk '{print $5,$9}'

echo "[$(ts)] 3. genotype the sample against the step-2 contract"
"$PY" "$REPO/src/main.py" --step genotype --bam fp_stack.bam \
    --out fp.genotypes.txt.gz --insertions step2.genotyping.txt.gz --threads "$THREADS" >genotype.log 2>&1
tail -2 genotype.log

echo "[$(ts)] 4. combine_genotypes -> panel matrix (single sample)"
"$PY" "$REPO/src/main.py" --step combine_genotypes \
    --genotypes fp.genotypes.txt.gz --out matrix.csv.gz --threads "$THREADS" >combine_gt.log 2>&1
tail -2 combine_gt.log

echo "[$(ts)] 5. score discovery vs post-combine (KEPT insertions live in genotyping.txt.gz,"
echo "         not combined.txt.gz which is the clip FASTQ of every candidate)"
zcat < step2.genotyping.txt.gz | awk '/^>/{sub(/^>/,"@"); print $0":L"}' | gzip > kept_calls.txt.gz
"$PY" "$DIR/score_fp10k.py" truth_hg38.tsv \
    discovery=$DISC combine_kept=kept_calls.txt.gz --window 50 | tee score.pipeline.txt
echo "[$(ts)] PIPELINE DONE (matrix: $OUT/matrix.csv.gz)"
