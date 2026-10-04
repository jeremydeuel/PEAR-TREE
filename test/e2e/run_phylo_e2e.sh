#!/bin/bash
# Phylogenetic evaluation E2E (plans/tprt_hallmarks/PHYLO_EVAL.md):
#   tree-mode fullstack patient (TP insertions on the branches of a random coalescent tree, carried
#   by exactly the clade at clonal het x colony purity 0.7-1.0, low-depth colonies, non-clade
#   decoys) -> run_e2e.sh (map, discovery, combine, annotate) -> run_genotype_e2e.sh (Rust
#   genotyper per colony, .tprt contract) -> tools/phylo/tree_fit.py -> calibration vs truth
#   (tools/phylo/calibration.py) -> tools/phylo/discrimination.py
#
# Usage: bash test/e2e/run_phylo_e2e.sh       (env overrides below; re-runs skip finished stages)
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../.." && pwd)"
export SP="${SP:-/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad}"
export PY="${PY:-$SP/venv/bin/python}"
export SAMPLES="${SAMPLES:-10}"
export TREE="${TREE:-random:$SAMPLES}"
export OUT="${OUT:-$SP/work/e2e_phylo}"
export N_PER_TYPE="${N_PER_TYPE:-10}"
export DEPTH="${DEPTH:-15}"
export SEED="${SEED:-21}"
export DONOR_ARGS="${DONOR_ARGS:---tree $TREE}"
ts(){ date "+%H:%M:%S"; }

echo "[$(ts)] A. simulate + discovery + combine + annotate ($SAMPLES colonies, $TREE)"
bash "$DIR/run_e2e.sh"
echo "[$(ts)] B. genotype every colony"
[ -s "$OUT/genotype/tprt/S$SAMPLES.txt.gz" ] || bash "$DIR/run_genotype_e2e.sh"

PH="$OUT/phylo"
mkdir -p "$PH"
echo "[$(ts)] C. truth table (genotyped loci -> simulated events)"
"$PY" "$DIR/phylo_truth.py" --e2e-dir "$OUT" --genotype-dir "$OUT/genotype/tprt" --out "$PH/truth.tsv"
echo "[$(ts)] D. tree_fit (the farm interface)"
"$PY" "$REPO/tools/phylo/tree_fit.py" --genotypes "$OUT/annot/P1.genotypes.csv.gz" \
    --genotype-dir "$OUT/genotype/tprt" --tree "$OUT/donor/tree.nwk" \
    --annotation "$OUT/annot/P1.annotated.tsv" --out "$PH/fit" --sex F ${TREE_FIT_ARGS:-}
echo "[$(ts)] E. calibration vs truth"
"$PY" "$REPO/tools/phylo/calibration.py" --fit "$PH/fit" --truth "$PH/truth.tsv" --out "$PH/calibration" \
    --label "E2E $SAMPLES colonies, $TREE, seed $SEED"
echo "[$(ts)] F. discrimination"
"$PY" "$REPO/tools/phylo/discrimination.py" --fit "$PH/fit/phylo_fit.tsv" --out "$PH/discrimination" \
    --truth "$PH/truth.tsv"
echo "[$(ts)] DONE -> $PH (fit/summary.md, calibration/report.md, discrimination/report.md)"
