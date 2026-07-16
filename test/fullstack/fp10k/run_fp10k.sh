#!/bin/bash
# FP-stress full-stack harness: 10k true insertions (all classes) + a large non-MEI decoy
# population tuned so raw *discovery* throws >=10 000 EMERGENT false positives.
#
#   RepeatMasker -> LCR/young-RTE/old-RTE decoy site BEDs + euchromatin TP windows
#   -> donor (10k canonical MEIs + ~16k non-MEI decoy cassettes)
#   -> wgsim -> library artefacts -> bwa-mem hg38 -> markdup BAM
#   -> rust discovery (discovery_hs.config) -> chain-lift TP truth -> score.
#
# Needs: hs1.fa(.fai), bwa-indexed hg38.fa, hs1ToHg38 chain, hs1 RepeatMasker .out.gz.
# All configurable by env var (defaults assume ~/Downloads and the repo venv).
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../../.." && pwd)"
SCALE="$REPO/test/fullstack/scale10k"
PY="${PY:-$REPO/venv/bin/python}"
RUST="${RUST:-$REPO/rust/peartree-discovery/target/release/peartree-discovery}"

HS1="${HS1:-$HOME/Downloads/hs1.fa}"
HG38="${HG38:-$HOME/Downloads/hg38.fa}"
CHAIN="${CHAIN:-$HOME/Downloads/hs1ToHg38.over.chain.gz}"
RMSK="${RMSK:-$REPO/hs1.repeatMasker.out.gz}"
SAMTOOLS="${SAMTOOLS:-$(command -v samtools)}"
OUT="${OUT:-$HOME/Downloads/fullstack_fp10k}"
THREADS="${THREADS:-8}"
DEPTH="${DEPTH:-25}"
N_IMPLANTS="${N_IMPLANTS:-10000}"
N_DECOYS="${N_DECOYS:-14000}"
FP_MB="${FP_MB:-35}"
ts(){ date "+%H:%M:%S"; }
mkdir -p "$OUT"; cd "$OUT"

echo "[$(ts)] 1. euchromatin TP + FP windows + RepeatMasker decoy site BEDs"
[ -f sat_bins.tsv ] || bash "$SCALE/profile_satellite.sh" "$RMSK" sat_bins.tsv
"$PY" "$SCALE/gen_windows.py" sat_bins.tsv "$HS1.fai" "$OUT"   # writes tp_windows.bed + fp_windows.bed (satellite)
"$PY" "$DIR/gen_fp_windows.py" "$HS1.fai" tp_windows.bed fp_windows.bed fp_eu_windows.bed "$FP_MB"

# One RepeatMasker pass -> decoy site BEDs, thinned to one representative element per
# (contig, 50 kb bin) per category so sites are spread genome-wide but the BEDs stay small.
# RMSK cols: $5 contig $6 begin $7 end $9 strand $10 repeat $11 class/family.
if [ ! -f lcr_sites.bed ]; then
  gzip -dc "$RMSK" | awk 'NR>3 && $5 ~ /^chr[0-9XY]+$/ {
    c=$5; b=$6+0; e=$7+0; rep=$10; fam=$11; L=e-b; bin=int(b/50000);
    if (fam ~ /Low_complexity|Simple_repeat|Satellite/ && L>=50) {
      k=c"\t"bin"\tL"; if(!(k in sl)){sl[k]=1; print c,b,e,rep > "lcr_sites.bed"}
    } else if ((fam ~ /^SINE\/Alu/ && rep ~ /^AluY/ && L>=250) ||
               (fam ~ /^LINE\/L1/  && rep ~ /^L1(HS|P|PA1)/ && L>=1000) ||
               (fam ~ /Retroposon\/SVA/ && L>=800) ||
               (fam ~ /LTR\/ERVK/ && L>=800)) {
      k=c"\t"bin"\tY"; if(!(k in sy)){sy[k]=1; print c,b,e,rep > "young_sites.bed"}
    } else if ((fam ~ /^SINE\/Alu/ && rep ~ /^Alu(S|J)/ && L>=250) ||
               (fam ~ /SINE\/MIR/ && L>=180) ||
               (fam ~ /LINE\/L2/ && L>=400) ||
               (fam ~ /LINE\/CR1/ && L>=500) ||
               (fam ~ /LTR\/ERVL/ && L>=600)) {
      k=c"\t"bin"\tO"; if(!(k in so)){so[k]=1; print c,b,e,rep > "old_sites.bed"}
    }
  }'
fi
echo "    LCR sites:   $(wc -l < lcr_sites.bed)"
echo "    young sites: $(wc -l < young_sites.bed)"
echo "    old sites:   $(wc -l < old_sites.bed)"

echo "[$(ts)] 2. build donor ($N_IMPLANTS true implants + $N_DECOYS decoys)"
"$PY" "$DIR/build_fp10k.py" --hs1 "$HS1" --tp-bed tp_windows.bed --fp-bed fp_eu_windows.bed \
    --lcr-bed lcr_sites.bed --young-bed young_sites.bed --old-bed old_sites.bed \
    --pseudogenes "$REPO/test/genotyping/pseudogenes.fa" \
    --n-implants "$N_IMPLANTS" --n-decoys "$N_DECOYS" \
    --out-donor donor.fa --out-truth truth_hs1.tsv --out-flanks flanks.fa
samtools faidx donor.fa
BP=$(awk '{s+=$2} END{print s}' donor.fa.fai); N=$(( BP * DEPTH / 300 ))
echo "    donor ${BP} bp -> wgsim $N pairs at ${DEPTH}x"

echo "[$(ts)] 3. wgsim ~${DEPTH}x, inject library artefacts"
wgsim -N "$N" -1 150 -2 150 -d 320 -s 40 -e 0.005 -r 0 -R 0 -X 0 donor.fa r1.fq r2.fq >wgsim.log 2>&1
"$PY" "$REPO/test/fullstack/inject_artefacts.py" --in1 r1.fq --in2 r2.fq --out1 r1a.fq --out2 r2a.fq
rm -f r1.fq r2.fq

echo "[$(ts)] 4. bwa-mem -> hg38, fixmate/sort/markdup"
bwa mem -t "$THREADS" -R '@RG\tID:fp10k\tSM:fp10k\tPL:ILLUMINA' "$HG38" r1a.fq r2a.fq 2>bwa.log \
  | samtools fixmate -m -u -@ "$THREADS" - - | samtools sort -u -@ "$THREADS" - \
  | samtools markdup -@ "$THREADS" - fp_stack.bam
samtools index fp_stack.bam
bwa mem -t "$THREADS" "$HG38" flanks.fa 2>>bwa.log | samtools sort -o flanks.bam -; samtools index flanks.bam
rm -f r1a.fq r2a.fq

echo "[$(ts)] 5. rust discovery (recommended human config; + plain min_mapq40 for reference)"
"$RUST" --step discover --bam fp_stack.bam --out discovery.txt.gz --config "$SCALE/discovery_hs.config" --threads "$THREADS"
printf 'min_mapq = 40\n' > plain.config
"$RUST" --step discover --bam fp_stack.bam --out discovery.plain.txt.gz --config plain.config --threads "$THREADS"

echo "[$(ts)] 6. chain-lift TP truth to hg38"
"$PY" "$REPO/test/fullstack/lift_truth.py" --flanks-bam flanks.bam --truth-hs1 truth_hs1.tsv \
    --out truth_hg38.tsv --chain "$CHAIN"

echo "[$(ts)] 7. score"
"$PY" "$DIR/score_fp10k.py" truth_hg38.tsv \
    recommended=discovery.txt.gz plain_mq40=discovery.plain.txt.gz --window 50 | tee score.txt
echo "[$(ts)] DONE  (BAM: $OUT/fp_stack.bam, truth: $OUT/truth_hs1.tsv / truth_hg38.tsv)"
