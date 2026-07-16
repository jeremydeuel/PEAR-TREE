#!/bin/bash
# End-to-end 10k/10k full-stack gate:
#   satellite profile -> windows -> donor (10k implants + centromere/telomere FP compartment)
#   -> wgsim -> library artefacts -> bwa-mem hg38 -> markdup BAM
#   -> rust discovery (coverage_mask) -> chain-lift truth -> python combine -> score.
#
# Needs: hs1.fa, bwa-indexed hg38.fa, bowtie2-indexed hs1, hs1ToHg38 chain, hs1 RepeatMasker.
# All configurable by env var (defaults assume ~/Downloads and the repo venv).
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../../.." && pwd)"
PY="${PY:-$REPO/venv/bin/python}"
RUST="${RUST:-$REPO/rust/peartree-discovery/target/release/peartree-discovery}"

HS1="${HS1:-$HOME/Downloads/hs1.fa}"
HG38="${HG38:-$HOME/Downloads/hg38.fa}"
HS1_BT2="${HS1_BT2:-$HOME/Downloads/hs1}"                       # bowtie2 index prefix
CHAIN="${CHAIN:-$HOME/Downloads/hs1ToHg38.over.chain.gz}"
RMSK="${RMSK:-$REPO/hs1.repeatMasker.out.gz}"
SAMTOOLS="${SAMTOOLS:-$(command -v samtools)}"
BOWTIE2="${BOWTIE2:-$(command -v bowtie2)}"
OUT="${OUT:-$HOME/Downloads/fullstack_scale10k}"
THREADS="${THREADS:-8}"
DEPTH="${DEPTH:-25}"
N_IMPLANTS="${N_IMPLANTS:-10000}"
ts(){ date "+%H:%M:%S"; }
mkdir -p "$OUT"; cd "$OUT"

echo "[$(ts)] 0. hg38 primary .2bit (for combine genotyping loader)"
[ -f hg38.primary.2bit ] || "$PY" "$DIR/make_2bit.py" "$HG38" hg38.primary.2bit

echo "[$(ts)] 1. satellite profile + windows"
[ -f sat_bins.tsv ] || bash "$DIR/profile_satellite.sh" "$RMSK" sat_bins.tsv
"$PY" "$DIR/gen_windows.py" sat_bins.tsv "$HS1.fai" "$OUT"

echo "[$(ts)] 2. build donor ($N_IMPLANTS implants + FP compartment)"
"$PY" "$DIR/build_donor_10k.py" --hs1 "$HS1" --tp-bed tp_windows.bed --fp-bed fp_windows.bed \
    --n-implants "$N_IMPLANTS" --out-donor donor.fa --out-truth truth_hs1.tsv --out-flanks flanks.fa
samtools faidx donor.fa
BP=$(awk '{s+=$2} END{print s}' donor.fa.fai); N=$(( BP * DEPTH / 300 ))

echo "[$(ts)] 3. wgsim ~${DEPTH}x -> $N pairs, inject artefacts"
wgsim -N "$N" -1 150 -2 150 -d 320 -s 40 -e 0.005 -r 0 -R 0 -X 0 donor.fa r1.fq r2.fq >wgsim.log 2>&1
"$PY" "$REPO/test/fullstack/inject_artefacts.py" --in1 r1.fq --in2 r2.fq --out1 r1a.fq --out2 r2a.fq
rm -f r1.fq r2.fq

echo "[$(ts)] 4. bwa-mem -> hg38, fixmate/sort/markdup"
bwa mem -t "$THREADS" -R '@RG\tID:hs1sim\tSM:hs1sim\tPL:ILLUMINA' "$HG38" r1a.fq r2a.fq 2>bwa.log \
  | samtools fixmate -m -u -@ "$THREADS" - - | samtools sort -u -@ "$THREADS" - \
  | samtools markdup -@ "$THREADS" - full_stack.bam
samtools index full_stack.bam
bwa mem -t "$THREADS" "$HG38" flanks.fa 2>>bwa.log | samtools sort -o flanks.bam -; samtools index flanks.bam
rm -f r1a.fq r2a.fq

echo "[$(ts)] 5. rust discovery (coverage_mask gate)"
"$RUST" --step discover --bam full_stack.bam --out discovery.txt.gz --config "$DIR/discovery_hs.config" --threads "$THREADS"

echo "[$(ts)] 6. chain-lift truth to hg38"
"$PY" "$REPO/test/fullstack/lift_truth.py" --flanks-bam flanks.bam --truth-hs1 truth_hs1.tsv \
    --out truth_hg38.tsv --chain "$CHAIN"

echo "[$(ts)] 7. python combine (step 2)"
sed -e "s#@@HG38_2BIT@@#$OUT/hg38.primary.2bit#" -e "s#@@SAMTOOLS@@#$SAMTOOLS#" \
    -e "s#@@BOWTIE2@@#$BOWTIE2#" -e "s#@@HS1_BT2@@#$HS1_BT2#" -e "s#@@HS1_HG38_CHAIN@@#$CHAIN#" \
    "$DIR/config.local.py.template" > config.py
"$PY" "$DIR/run_combine.py" discovery.txt.gz step2 "$THREADS" "$OUT" "$REPO/src" >combine.log 2>&1
tail -2 combine.log

echo "[$(ts)] 8. score"
"$PY" "$DIR/score_10k.py" truth_hg38.tsv discovery=discovery.txt.gz combined=step2.combined.txt.gz --window 50 | tee score.txt

echo "[$(ts)] 9. annotate (family calls)"
# End-to-end annotation: nhmmscan family scan (build_hmm.sh library) + hs1 clip remap + RMSK +
# pseudogene exon track. scale10k has no genotyping step, so synthesise a 1-tip all-heterozygous
# genotypes file (every combined call gets annotated — a family-recovery check, not a genotype
# test). Needs the HMM library built once: test/fullstack/annotate/build_hmm.sh --from-hs1 $HS1
ANNOT="$DIR/../annotate"
if [ -f "$ANNOT/peartree_rte.hmm" ]; then
    "$PY" - "$OUT/step2.combined.txt.gz" "$OUT/step2.genotypes.csv.gz" <<'PYEOF'
import gzip, sys
comb, out = sys.argv[1], sys.argv[2]
keys = [l.strip()[1:-2] for l in gzip.open(comb, 'rt') if l.startswith('@') and l.rstrip().endswith('L')]
with gzip.open(out, 'wt') as o:
    o.write(';tip1\n')
    for k in keys:
        o.write(f'{k};heterozygous\n')
PYEOF
    HS1_BT2="$HS1_BT2" BOWTIE2="$BOWTIE2" PT_RMSK="$RMSK" \
        "$PY" "$ANNOT/run_annotate.py" --combined "$OUT/step2.combined.txt.gz" \
            --genotypes "$OUT/step2.genotypes.csv.gz" --out "$OUT/annotate.txt" --workdir "$OUT/annot"
    echo "  annotation report -> $OUT/annotate.txt"
else
    echo "  (skipped: build the HMM library first -> $ANNOT/build_hmm.sh --from-hs1 \$HS1)"
fi
echo "[$(ts)] DONE  (outputs in $OUT)"
