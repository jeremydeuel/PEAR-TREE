#!/bin/bash
# GENOTYPING panel gate (test/genotyping).
#
#   two-haplotype hs1 donor (ins + ref, 1000 het implants incl. processed pseudogenes)
#   -> per-sample wgsim MIX (het = ins@D/2 + ref@D/2 ; wild-type = ref@D) -> bwa-mem hg38
#   -> 10 markdup'd BAMs -> rust discovery on the 100x sample -> combine_insertions (contract)
#   -> genotype all 10 against that contract -> combine_genotypes -> lift truth -> score.
#
# The 3 wild-type samples are the discriminator: a true het reads wild-type there, an
# assembly-discordance FP (present identically in all 10) does not -> combine_genotypes
# separates them. See plans/improve_discovery/genotyping_simulation.md.
#
# Needs: hs1.fa, bwa-indexed hg38.fa, bowtie2-indexed hs1, hs1ToHg38 chain, the built
# rust discovery binary, samtools/bwa/wgsim, and test/genotyping/pseudogenes.fa.
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../.." && pwd)"
SCALE="$REPO/test/fullstack/scale10k"
PY="${PY:-$REPO/venv/bin/python}"
RUST="${RUST:-$REPO/rust/peartree-discovery/target/release/peartree-discovery}"

HS1="${HS1:-$HOME/Downloads/hs1.fa}"
HG38="${HG38:-$HOME/Downloads/hg38.fa}"
HS1_BT2="${HS1_BT2:-$HOME/Downloads/hs1}"
CHAIN="${CHAIN:-$HOME/Downloads/hs1ToHg38.over.chain.gz}"
SAMTOOLS="${SAMTOOLS:-$(command -v samtools)}"
BOWTIE2="${BOWTIE2:-$(command -v bowtie2)}"
OUT="${OUT:-$HOME/Downloads/genotyping_panel}"
THREADS="${THREADS:-8}"          # bwa / rust discovery
MP="${MP:-1}"                    # python multiprocessing steps (genotype/combine); 1 = safe on macOS spawn
N_IMPLANTS="${N_IMPLANTS:-1000}"
N_PSEUDOGENE="${N_PSEUDOGENE:-100}"
SEED="${SEED:-1}"
# het sample depths (the FIRST is the discovery sample) + wild-type depths
HET_DEPTHS="${HET_DEPTHS:-100 80 60 40 20 10 5}"
WT_DEPTHS="${WT_DEPTHS:-40 40 40}"
LOW_COV_GATE="${LOW_COV_GATE:-0}"   # 1 -> reads_for_high_coverage=60 (exercise the gate)
ts(){ date "+%H:%M:%S"; }
mkdir -p "$OUT"; cd "$OUT"

echo "[$(ts)] 0. hg38 primary .2bit + panel config.py"
[ -f hg38.primary.2bit ] || "$PY" "$SCALE/make_2bit.py" "$HG38" hg38.primary.2bit
sed -e "s#@@HG38_2BIT@@#$OUT/hg38.primary.2bit#" -e "s#@@SAMTOOLS@@#$SAMTOOLS#" \
    -e "s#@@BOWTIE2@@#$BOWTIE2#" -e "s#@@HS1_BT2@@#$HS1_BT2#" -e "s#@@HS1_HG38_CHAIN@@#$CHAIN#" \
    "$DIR/config.panel.py.template" > config.py
if [ "$LOW_COV_GATE" = "1" ]; then
    sed -i.bak "s/'reads_for_high_coverage': 250/'reads_for_high_coverage': 60/" config.py
    echo "  (LOW_COV_GATE: reads_for_high_coverage=60)"
fi

echo "[$(ts)] 1. TP/FP windows (fullstack-validated high-mapability euchromatin + pericentromere)"
# hs1->hg38 mapping puts reads in low-mapability homologs at MAPQ 0 -> discovery drops them
# (min_mapq). These two windows are the fullstack-proven high-mapability euchromatin (~88%
# cross-assembly recall; the residual loss is the chr21:22.1-22.2M LINE-dense patch that no
# window choice avoids). Override with TP_WINDOWS / FP_WINDOWS (printf %b format) to retune.
printf "%b" "${TP_WINDOWS:-chr21\t20000000\t26000000\tTP\nchr22\t25000000\t31000000\tTP\n}" > tp_windows.bed
printf "%b" "${FP_WINDOWS:-chr21\t10000000\t13000000\tpericentromere\n}" > fp_windows.bed

echo "[$(ts)] 2. build two haplotypes ($N_IMPLANTS implants, $N_PSEUDOGENE pseudogenes)"
"$PY" "$DIR/build_haplotypes.py" --hs1 "$HS1" --tp-bed tp_windows.bed --fp-bed fp_windows.bed \
    --pseudogene-fasta "$DIR/pseudogenes.fa" --n-implants "$N_IMPLANTS" --n-pseudogene "$N_PSEUDOGENE" \
    --seed "$SEED" --out-donor-ins donor_ins.fa --out-donor-ref donor_ref.fa \
    --out-truth truth_hs1.tsv --out-flanks flanks.fa
samtools faidx donor_ins.fa; samtools faidx donor_ref.fa
INS_BP=$(awk '{s+=$2} END{print s}' donor_ins.fa.fai)
REF_BP=$(awk '{s+=$2} END{print s}' donor_ref.fa.fai)
echo "  donor_ins ${INS_BP} bp, donor_ref ${REF_BP} bp"

# map one FASTQ pair to hg38 -> $1.bam  (args: name r1 r2)
map_sample(){
  local name="$1" r1="$2" r2="$3"
  bwa mem -t "$THREADS" -R "@RG\tID:$name\tSM:$name\tPL:ILLUMINA" "$HG38" "$r1" "$r2" 2>>bwa.log \
    | samtools fixmate -m -u -@ "$THREADS" - - \
    | samtools sort -u -@ "$THREADS" - \
    | samtools markdup -@ "$THREADS" - "$name.bam"
  samtools index "$name.bam"
}
# wgsim `reads for coverage C over B bp` = C*B/(2*150)
npairs(){ python3 -c "print(int($1*$2/300))"; }
WG="wgsim -1 150 -2 150 -d 320 -s 40 -e 0.005 -r 0 -R 0 -X 0"

SAMPLES=(); DISC_BAM=""
echo "[$(ts)] 3. het samples (ins@D/2 + ref@D/2)"
first=1
for D in $HET_DEPTHS; do
  name="het_${D}x"
  $WG -N "$(npairs $((D/2)) "$INS_BP")" donor_ins.fa "${name}_i1.fq" "${name}_i2.fq" >>wgsim.log 2>&1
  $WG -N "$(npairs $((D/2)) "$REF_BP")" donor_ref.fa "${name}_r1.fq" "${name}_r2.fq" >>wgsim.log 2>&1
  cat "${name}_i1.fq" "${name}_r1.fq" > "${name}_1.fq"; cat "${name}_i2.fq" "${name}_r2.fq" > "${name}_2.fq"
  rm -f "${name}_i1.fq" "${name}_i2.fq" "${name}_r1.fq" "${name}_r2.fq"
  map_sample "$name" "${name}_1.fq" "${name}_2.fq"; rm -f "${name}_1.fq" "${name}_2.fq"
  SAMPLES+=("$name")
  if [ "$first" = 1 ]; then DISC_BAM="$name.bam"; first=0; fi   # first het depth = discovery sample
  echo "  [$(ts)] $name done"
done

echo "[$(ts)] 4. wild-type samples (ref@D only)"
WT_LETTERS=(a b c d e f g h)
wi=0
for D in $WT_DEPTHS; do
  name="wt_${D}x_${WT_LETTERS[$wi]}"; wi=$((wi+1))
  $WG -N "$(npairs "$D" "$REF_BP")" donor_ref.fa "${name}_1.fq" "${name}_2.fq" >>wgsim.log 2>&1
  map_sample "$name" "${name}_1.fq" "${name}_2.fq"; rm -f "${name}_1.fq" "${name}_2.fq"
  SAMPLES+=("$name")
  echo "  [$(ts)] $name done"
done

echo "[$(ts)] 5. discovery on the 100x sample ($DISC_BAM)"
# Feature B (processed-pseudogene splice): this harness implants pseudogenes, so enable the
# mate-splice hallmark. Build the hg38 exon track once (mates land in hg38), then append the
# splice block to the discovery config; combine emits <combined>.splice.tsv for step 11.
ANNOT="$REPO/test/fullstack/annotate"
EXON_HG38="$ANNOT/pseudogene_exons.hg38.bed"
[ -f "$EXON_HG38" ] || "$PY" "$ANNOT/build_exon_track.py" --pseudogene-fasta "$DIR/pseudogenes.fa" \
    --bwa "$(command -v bwa)" --ref "$HG38" --min-mapq 0 --out "$EXON_HG38"
cp "$SCALE/discovery_hs.config" "$OUT/discovery.config"
cat >> "$OUT/discovery.config" <<EOF

# Feature B — processed-pseudogene splice annotation (wired by run_genotyping.sh)
splice_hallmark = true
exon_annotation = $EXON_HG38
splice_min_exons = 2
EOF
"$RUST" --step discover --bam "$DISC_BAM" --out discovery.txt.gz \
    --config "$OUT/discovery.config" --threads "$THREADS"

echo "[$(ts)] 6. combine_insertions -> genotyping contract"
"$PY" "$SCALE/run_combine.py" discovery.txt.gz step2 "$MP" "$OUT" "$REPO/src" > combine.log 2>&1
tail -1 combine.log

echo "[$(ts)] 7. lift truth to hg38"
bwa mem -t "$THREADS" "$HG38" flanks.fa 2>>bwa.log | samtools sort -o flanks.bam -; samtools index flanks.bam
"$PY" "$REPO/test/fullstack/lift_truth.py" --flanks-bam flanks.bam --truth-hs1 truth_hs1.tsv \
    --out truth_hg38.tsv --chain "$CHAIN"

echo "[$(ts)] 8. genotype all ${#SAMPLES[@]} samples against the contract"
GT_FILES=()
for name in "${SAMPLES[@]}"; do
  "$PY" "$DIR/run_step.py" "$OUT" "$REPO/src" genotype "$name.bam" \
      "$name.genotypes.txt.gz" step2.genotyping.txt.gz "$MP" > "genotype_$name.log" 2>&1
  GT_FILES+=("$name.genotypes.txt.gz")
  echo "  [$(ts)] genotyped $name"
done

echo "[$(ts)] 9. combine_genotypes -> matrix"
"$PY" "$DIR/run_step.py" "$OUT" "$REPO/src" combine_genotypes matrix.csv.gz "$MP" \
    "${GT_FILES[@]}" > combine_genotypes.log 2>&1
tail -3 combine_genotypes.log

echo "[$(ts)] 10. score"
"$PY" "$DIR/score_genotyping.py" --truth truth_hg38.tsv --contract step2.genotyping.txt.gz \
    --matrix matrix.csv.gz --genotypes "${GT_FILES[@]}" \
    --het-depths "$HET_DEPTHS" --disc-bam "$DISC_BAM" --window 50 | tee score.txt

echo "[$(ts)] 11. annotate (family + pseudogene calls)"
# Family/pseudogene annotation of the combined calls. The splice sidecar from step 5/6 gives
# the mate-splice pseudogene calls; the HMM/RMSK/exon paths give RTE + single-locus pseudogene.
# No genotyping columns are needed for a family-recovery check, so synthesise a 1-tip all-het
# genotypes file. Needs the HMM library once: test/fullstack/annotate/build_hmm.sh
if [ -f "$ANNOT/peartree_rte.hmm" ]; then
    "$PY" - "$OUT/step2.combined.txt.gz" "$OUT/annot.genotypes.csv.gz" <<'PYEOF'
import gzip, sys
comb, out = sys.argv[1], sys.argv[2]
keys = [l.strip()[1:-2] for l in gzip.open(comb, 'rt') if l.startswith('@') and l.rstrip().endswith('L')]
with gzip.open(out, 'wt') as o:
    o.write(';tip1\n')
    for k in keys:
        o.write(f'{k};heterozygous\n')
PYEOF
    HS1_BT2="$HS1_BT2" BOWTIE2="$BOWTIE2" \
        "$PY" "$ANNOT/run_annotate.py" --combined "$OUT/step2.combined.txt.gz" \
            --genotypes "$OUT/annot.genotypes.csv.gz" --out "$OUT/annotate.txt" --workdir "$OUT/annot"
    echo "  annotation report -> $OUT/annotate.txt"
else
    echo "  (skipped: build the HMM library first -> $ANNOT/build_hmm.sh)"
fi
echo "[$(ts)] DONE (outputs in $OUT)"
