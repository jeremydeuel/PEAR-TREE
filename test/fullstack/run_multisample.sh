#!/bin/bash
# Multi-sample (one patient, n colonies) insertion-type catalogue test:
#   build_donor.py --types (hs1 windows + catalogue events, per-sample haplotypes)
#   -> simulate_reads.py per sample (poly-A jitter, unflagged PCR dups)
#   -> bwa-mem to the discovery reference -> fixmate/sort (NO markdup by default: the
#      independence rule must work from the reads, never from the 0x400 flag)
#   -> rust discovery per sample -> score_types.py (per-type recall, pooled + per sample)
#
# Reference: HG38=<bwa-indexed FASTA> (whole GRCh38), OR a reduced reference built on the fly
# from HG38_2BIT + HG38_CONTIGS (default chr22) + the hs1 transduction-source loci as decoy
# contigs (smoke-test mode; no whole-genome bwa index needed).
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../.." && pwd)"
HS1="${HS1:-$HOME/Downloads/hs1.fa}"                      # hs1 FASTA or .2bit
HS1_RMSK="${HS1_RMSK:-}"                                   # hs1 RepeatMasker .out.gz (fallback library)
RTE_LIBRARY="${RTE_LIBRARY:-$REPO/resources/rte_library}"  # used when populated
HG38="${HG38:-}"
HG38_2BIT="${HG38_2BIT:-}"
HG38_CONTIGS="${HG38_CONTIGS:-chr22}"
OUT="${OUT:-$HOME/Downloads/fullstack_multisample}"
PY="${PY:-python3}"
DISC_BIN="${DISC_BIN:-$REPO/rust/peartree-discovery/target/release/peartree-discovery}"
DISC_CONFIG="${DISC_CONFIG:-}"
THREADS="${THREADS:-8}"
SAMPLES="${SAMPLES:-2}"
DEPTH="${DEPTH:-15}"
TYPES="${TYPES:-all}"
N_PER_TYPE="${N_PER_TYPE:-3}"
REGION="${REGION:-chr22:26000000-30000000}"
SEED="${SEED:-1}"
PCR_DUP="${PCR_DUP:-0.05}"
PCR_DUP_JITTER="${PCR_DUP_JITTER:-0}"   # +-bp start/end jitter of unflagged PCR copies
POLYA_JITTER="${POLYA_JITTER:-1.0}"
MARKDUP="${MARKDUP:-0}"
WINDOW="${WINDOW:-30}"
mkdir -p "$OUT"
ts(){ date "+%H:%M:%S"; }

LIBARGS=()
[ -f "$RTE_LIBRARY/l1_intact.fa" ] && LIBARGS+=(--rte-library "$RTE_LIBRARY")
[ -n "$HS1_RMSK" ] && LIBARGS+=(--hs1-rmsk "$HS1_RMSK" --library-cache "$OUT/rte_cache")

echo "[$(ts)] 1. donor: $TYPES x $N_PER_TYPE, $SAMPLES samples, $REGION"
REGARGS=(); for r in $REGION; do REGARGS+=(--region "$r"); done
DARGS=(); for x in ${DONOR_ARGS:-}; do DARGS+=("$x"); done   # e.g. DONOR_ARGS="--tree random:10"
"$PY" "$DIR/build_donor.py" --types "$TYPES" --n-per-type "$N_PER_TYPE" --samples "$SAMPLES" \
    --hs1 "$HS1" "${LIBARGS[@]}" "${REGARGS[@]}" --seed "$SEED" --out-dir "$OUT/donor" ${DARGS[@]+"${DARGS[@]}"}
# tree mode (--tree): per-colony depth factors in donor/samples.tsv; its colony count must be SAMPLES
if [ -f "$OUT/donor/samples.tsv" ]; then
    NS=$(($(wc -l < "$OUT/donor/samples.tsv") - 1))
    [ "$NS" = "$SAMPLES" ] || { echo "donor has $NS colonies (tree tips) but SAMPLES=$SAMPLES"; exit 1; }
fi

if [ -z "$HG38" ]; then
    [ -n "$HG38_2BIT" ] || { echo "set HG38 (bwa-indexed) or HG38_2BIT for a reduced reference"; exit 1; }
    HG38="$OUT/ref/reduced.fa"
    if [ ! -f "$HG38.bwt" ]; then
        echo "[$(ts)] 1b. reduced reference: $HG38_CONTIGS from $HG38_2BIT + hs1 source decoys"
        mkdir -p "$OUT/ref"
        "$PY" - "$HG38_2BIT" "$HG38" "$OUT/donor/sources_hs1.fa" $HG38_CONTIGS <<'EOF'
import sys, py2bit
tb = py2bit.open(sys.argv[1]); out = open(sys.argv[2], "w")
for c in sys.argv[4:]:
    s = tb.sequence(c, 0, tb.chroms(c))
    out.write(f">{c}\n")
    for i in range(0, len(s), 80):
        out.write(s[i:i + 80] + "\n")
out.write(open(sys.argv[3]).read())
EOF
        samtools faidx "$HG38"
        bwa index "$HG38" > "$OUT/ref/bwa_index.log" 2>&1
    fi
fi

CALLS=(); SUPPORT=()
for i in $(seq 1 "$SAMPLES"); do
    S="S$i"
    SDEPTH="$DEPTH"
    if [ -f "$OUT/donor/samples.tsv" ]; then
        SDEPTH=$(awk -F'\t' -v s="$i" -v d="$DEPTH" '$1 == s { printf "%.3f", d * $4 }' "$OUT/donor/samples.tsv")
    fi
    echo "[$(ts)] 2. $S: reads (depth $SDEPTH) -> bwa-mem -> fixmate/sort"
    "$PY" "$DIR/simulate_reads.py" --donor-dir "$OUT/donor" --sample "$i" --out-prefix "$OUT/$S" \
        --depth "$SDEPTH" --seed "$SEED" --pcr-dup-unflagged-frac "$PCR_DUP" --pcr-dup-jitter "$PCR_DUP_JITTER" --polya-jitter "$POLYA_JITTER"
    bwa mem -t "$THREADS" -R "@RG\tID:$S\tSM:$S\tPL:ILLUMINA" "$HG38" "$OUT/${S}_R1.fq" "$OUT/${S}_R2.fq" \
        2>"$OUT/$S.bwa.log" \
      | samtools fixmate -m -u -@ "$THREADS" - - \
      | samtools sort -@ "$THREADS" -o "$OUT/$S.sorted.bam" -
    if [ "$MARKDUP" = "1" ]; then
        samtools markdup -@ "$THREADS" "$OUT/$S.sorted.bam" "$OUT/$S.bam"; rm "$OUT/$S.sorted.bam"
    else
        mv "$OUT/$S.sorted.bam" "$OUT/$S.bam"
    fi
    samtools index "$OUT/$S.bam"
    rm -f "$OUT/${S}_R1.fq" "$OUT/${S}_R2.fq"
    SUPPORT+=("$OUT/$S.support.tsv")
    if [ -x "$DISC_BIN" ]; then
        CFG=(); [ -n "$DISC_CONFIG" ] && CFG=(--config "$DISC_CONFIG")
        "$DISC_BIN" --step discover --bam "$OUT/$S.bam" --out "$OUT/$S.discovery.txt.gz" "${CFG[@]}" \
            > "$OUT/$S.discovery.log" 2>&1
        CALLS+=("$OUT/$S.discovery.txt.gz")
    fi
done

echo "[$(ts)] 3. lift truth (flank mapping) + score"
bwa mem -t "$THREADS" "$HG38" "$OUT/donor/flanks.fa" 2>/dev/null | samtools sort -o "$OUT/flanks.bam" -
samtools index "$OUT/flanks.bam"
if [ ${#CALLS[@]} -gt 0 ]; then
    "$PY" "$DIR/score_types.py" --truth "$OUT/donor/truth_types_hs1.tsv" --flanks-bam "$OUT/flanks.bam" \
        --calls "$(IFS=,; echo "${CALLS[*]}")" --support "$(IFS=,; echo "${SUPPORT[*]}")" \
        --window "$WINDOW" --out-tsv "$OUT/results_by_event.tsv" | tee "$OUT/score.txt"
else
    echo "  (no discovery binary at $DISC_BIN; skipping scoring)"
fi
echo "[$(ts)] DONE -> $OUT"
