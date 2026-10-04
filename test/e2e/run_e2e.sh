#!/bin/bash
# TPRT-hallmark end-to-end test (plans/tprt_hallmarks/E2E_REPORT.md):
#   fullstack multi-colony patient (all catalogue TP types + artefacts, unflagged PCR duplicates
#   with start/end jitter) -> bwa-mem to a reduced GRCh38 (chr22 + hs1 source decoys)
#   -> Rust discovery with the .tprt keys (evidence sidecar) per colony
#   -> combine_insertions with the .tprt python config (pooled >= 2 independent fragments,
#      lenient 2nd dedup, indel-aware consensus) -> annotate_v2 + tools/rte (resources/rte_library)
#   -> per-type scoring (test/e2e/score_e2e.py) + TPRT score calibration (tools/rte/calibrate.py)
#
# Everything big goes to $OUT. Re-running skips finished stages (delete $OUT to start over).
# Usage: bash test/e2e/run_e2e.sh            (env overrides below)
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../.." && pwd)"
SP="${SP:-/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad}"
PY="${PY:-$SP/venv/bin/python}"
GENOMES="${GENOMES:-$SP/genomes}"
OUT="${OUT:-$SP/work/e2e}"
SAMPLES="${SAMPLES:-3}"
N_PER_TYPE="${N_PER_TYPE:-8}"
DEPTH="${DEPTH:-15}"
PCR_DUP="${PCR_DUP:-0.10}"
PCR_DUP_JITTER="${PCR_DUP_JITTER:-3}"
REGION="${REGION:-chr22:19000000-33000000}"
THREADS="${THREADS:-8}"
SEED="${SEED:-7}"
mkdir -p "$OUT"
ts(){ date "+%H:%M:%S"; }

echo "[$(ts)] 0. build Rust discovery"
(cd "$REPO/rust/peartree-discovery" && cargo build --release 2>&1 | tail -1)
DISC_BIN="$REPO/rust/peartree-discovery/target/release/peartree-discovery"

# discovery config = cluster/config.discovery.grch38.tprt, minus what cannot run here: the
# GRCh38 exon track (splice_hallmark) lives on the farm. ignore_dup_flag=false (Jeremy
# 2026-10-04: 0x400 reads are dropped by discovery; the simulator never sets 0x400 anyway).
DCFG="$OUT/discovery.tprt.local"
sed -e 's/^splice_hallmark = true/splice_hallmark = false/' -e '/^exon_annotation = /d' \
    -e 's/^ignore_dup_flag = true/ignore_dup_flag = false/' \
    "$REPO/cluster/config.discovery.grch38.tprt" > "$DCFG"

if [ ! -f "$OUT/score.txt" ]; then
    echo "[$(ts)] 1-3. simulate $SAMPLES colonies, map, discover (test/fullstack/run_multisample.sh)"
    HS1="$GENOMES/hs1.2bit" HS1_RMSK="$GENOMES/hs1.repeatMasker.out.gz" HG38_2BIT="$GENOMES/hg38.2bit" \
    OUT="$OUT" N_PER_TYPE="$N_PER_TYPE" SAMPLES="$SAMPLES" DEPTH="$DEPTH" REGION="$REGION" \
    PCR_DUP="$PCR_DUP" PCR_DUP_JITTER="$PCR_DUP_JITTER" SEED="$SEED" PY="$PY" THREADS="$THREADS" \
    DISC_BIN="$DISC_BIN" DISC_CONFIG="$DCFG" \
        bash "$REPO/test/fullstack/run_multisample.sh"
fi

echo "[$(ts)] 4. combine/annotate references: bowtie2 index of the reduced reference, identity chain, rmsk"
REF="$OUT/ref/reduced.fa"
[ -f "$OUT/ref/reduced.1.bt2" ] || bowtie2-build --threads "$THREADS" "$REF" "$OUT/ref/reduced" > "$OUT/ref/bt2.log" 2>&1
"$PY" "$DIR/make_refs.py" --ref "$REF" --hg38-rmsk "$GENOMES/hg38.rmsk.txt.gz" --out-dir "$OUT/ref"
OVR=(); for o in ${CI_OVERRIDES:-}; do OVR+=(--override "$o"); done   # e.g. CI_OVERRIDES="slippage_reject=False"
"$PY" "$DIR/make_config.py" --repo "$REPO" --out "$OUT/pyconf/config.py" --ref-dir "$OUT/ref" \
    --hg38-2bit "$GENOMES/hg38.2bit" --hs1-2bit "$GENOMES/hs1.2bit" --workdir "$OUT/annot" ${OVR[@]+"${OVR[@]}"}
[ -f "$REPO/test/fullstack/annotate/peartree_rte.hmm.h3m" ] || bash "$REPO/test/fullstack/annotate/build_hmm.sh"

echo "[$(ts)] 5. combine_insertions (.tprt python config)"
CALLS=(); for i in $(seq 1 "$SAMPLES"); do CALLS+=("$OUT/S$i.discovery.txt.gz"); done
mkdir -p "$OUT/combine"
rm -f "$OUT/combine/P1".*
"$PY" "$DIR/run_combine.py" --config-dir "$OUT/pyconf" --out-stem "$OUT/combine/P1" --threads "$THREADS" \
    "${CALLS[@]}" > "$OUT/combine/combine.log" 2>&1
grep -E "independent-fragment gate|dropped:|SHORT overhang|evidence sidecars|indel-aware|wrote .* insertions" "$OUT/combine/combine.log" || true

echo "[$(ts)] 6. annotate_v2 + tools/rte"
mkdir -p "$OUT/annot"
rm -f "$OUT/annot/"*
"$PY" "$DIR/run_annotate_e2e.py" --config-dir "$OUT/pyconf" --combined "$OUT/combine/P1.combined.txt.gz" \
    --samples "$SAMPLES" --workdir "$OUT/annot" --table "$OUT/annot/P1.annotated.tsv" > "$OUT/annot/annotate.log" 2>&1
tail -2 "$OUT/annot/annotate.log"

echo "[$(ts)] 7. score per type + TPRT calibration (held-out half)"
"$PY" "$DIR/score_e2e.py" --out-dir "$OUT" --samples "$SAMPLES" --report "$OUT/e2e_tables.md" \
    --calib-truth "$OUT/calib_truth" | tee "$OUT/e2e_score.txt"
"$PY" "$DIR/classify_unexplained.py" --out-dir "$OUT" --samples "$SAMPLES" --hg38-2bit "$GENOMES/hg38.2bit" \
    --library "$REPO/resources/rte_library/consensus.fa" > "$OUT/unexplained_classes.md"
head -14 "$OUT/unexplained_classes.md"
for half in fit eval; do
    (cd "$REPO" && "$PY" -m tools.rte.calibrate --truth "$OUT/calib_truth.$half.tsv" --annot "$OUT/annot/P1.annotated.tsv" \
        --json "$OUT/calibrate.$half.json" > "$OUT/calibrate.$half.txt" 2>&1) || true
done
echo "[$(ts)] DONE -> $OUT (tables: $OUT/e2e_tables.md)"
