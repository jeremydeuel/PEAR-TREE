#!/bin/bash
# Genotyping stage of the TPRT-hallmark E2E (plans/tprt_hallmarks/E2E_REPORT.md, "Genotyping new
# locus kinds"). Needs a finished `test/e2e/run_e2e.sh` output in $OUT (colony BAMs, combine).
#
#   1. a wild-type control colony S<k> (reads from the event-free haplotype only, same simulator
#      and mapping) -> every TP locus is truth-absent there
#   2. contract + one-sided loci (src/genotyping_contract_oneside.py)
#   3. Rust genotyper per colony: legacy (config.genotype.grch38, legacy contract) and TPRT
#      (config.genotype.grch38.tprt, extended contract)
#   4. per-kind concordance vs simulator truth (test/e2e/score_genotypes.py)
#
# Usage: bash test/e2e/run_genotype_e2e.sh      (same env defaults as run_e2e.sh)
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../.." && pwd)"
SP="${SP:-/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad}"
PY="${PY:-$SP/venv/bin/python}"
OUT="${OUT:-$SP/work/e2e}"
SAMPLES="${SAMPLES:-3}"
DEPTH="${DEPTH:-15}"
PCR_DUP="${PCR_DUP:-0.10}"
PCR_DUP_JITTER="${PCR_DUP_JITTER:-3}"
SEED="${SEED:-7}"
THREADS="${THREADS:-8}"
GT="$OUT/genotype"
WT=$((SAMPLES + 1))
mkdir -p "$GT"
ts(){ date "+%H:%M:%S"; }

echo "[$(ts)] 0. build Rust genotyper"
(cd "$REPO/rust/peartree-genotype" && cargo build --release 2>&1 | tail -1)
BIN="$REPO/rust/peartree-genotype/target/release/peartree-genotype"

if [ ! -s "$GT/S$WT.bam.bai" ]; then
    echo "[$(ts)] 1. wild-type control colony S$WT (event-free haplotype only)"
    mkdir -p "$GT/wtdonor"
    for f in ref.hap.fa junctions.tsv slippage.tsv; do ln -sf "$OUT/donor/$f" "$GT/wtdonor/$f"; done
    printf 'sample\thap\tweight\n%s\tref.hap.fa\t1.0\n' "$WT" > "$GT/wtdonor/haps.tsv"
    : > "$GT/wtdonor/S$WT.molecules.fa"
    "$PY" "$REPO/test/fullstack/simulate_reads.py" --donor-dir "$GT/wtdonor" --sample "$WT" \
        --out-prefix "$GT/S$WT" --depth "$DEPTH" --seed "$SEED" \
        --pcr-dup-unflagged-frac "$PCR_DUP" --pcr-dup-jitter "$PCR_DUP_JITTER"
    bwa mem -t "$THREADS" -R "@RG\tID:S$WT\tSM:S$WT\tPL:ILLUMINA" "$OUT/ref/reduced.fa" \
        "$GT/S${WT}_R1.fq" "$GT/S${WT}_R2.fq" 2>"$GT/S$WT.bwa.log" \
      | samtools fixmate -m -u -@ "$THREADS" - - | samtools sort -@ "$THREADS" -o "$GT/S$WT.bam" -
    samtools index "$GT/S$WT.bam"
    rm -f "$GT/S${WT}_R1.fq" "$GT/S${WT}_R2.fq"
fi

echo "[$(ts)] 2. contract + one-sided loci"
"$PY" "$REPO/src/genotyping_contract_oneside.py" --contract "$OUT/combine/P1.genotyping.txt.gz" \
    --combined "$OUT/combine/P1.combined.txt.gz" --out "$GT/P1.genotyping.tprt.txt.gz" \
    --config-dir "$OUT/pyconf"

echo "[$(ts)] 3. genotype every colony (legacy + tprt)"
for mode in legacy tprt; do
    if [ "$mode" = legacy ]; then
        CFG="$REPO/cluster/config.genotype.grch38"; CON="$OUT/combine/P1.genotyping.txt.gz"
    else
        CFG="$REPO/cluster/config.genotype.grch38.tprt"; CON="$GT/P1.genotyping.tprt.txt.gz"
    fi
    mkdir -p "$GT/$mode"
    for i in $(seq 1 "$WT"); do
        BAM="$OUT/S$i.bam"; [ "$i" = "$WT" ] && BAM="$GT/S$WT.bam"
        "$BIN" --step genotype --bam "$BAM" --insertions "$CON" --out "$GT/$mode/S$i.txt.gz" \
            --threads "$THREADS" --config "$CFG" 2>"$GT/$mode/S$i.log"
    done
done

echo "[$(ts)] 4. concordance per locus kind"
for mode in legacy tprt; do
    python3 "$DIR/score_genotypes.py" --e2e-dir "$OUT" --geno-dir "$GT/$mode" --samples "$SAMPLES" \
        --wt-sample "S$WT" --label "($mode)"
    echo
done | tee "$GT/genotype_score.md"
echo "[$(ts)] DONE -> $GT/genotype_score.md"
