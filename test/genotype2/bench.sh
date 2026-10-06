#!/bin/bash
# Benchmark peartree-genotype2 against the legacy genotyper on the local simulated truth sets
# (plans/genotype_v2/SPEC.md "Benchmark"). Per data set: genotype every colony with both binaries
# (timed), per-kind concordance vs truth (test/e2e/score_genotypes.py), score distributions, the
# joint phylogenetic step vs truth (test/genotype2/score_joint.py), and tree_fit on both.
#
# Usage: bash test/genotype2/bench.sh [e2e_phylo|e2e_phylo2|e2e ...]   (default: e2e_phylo e2e)
# Env: SP (old scratchpad with the simulations), V2_CFG (config for the new binary), THREADS
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../.." && pwd)"
SP="${SP:-/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad}"
PY="${PY:-$SP/venv/bin/python}"
THREADS="${THREADS:-1}"
V2_CFG="${V2_CFG:-}"
SETS=("$@"); [ ${#SETS[@]} -gt 0 ] || SETS=(e2e_phylo e2e)
ts(){ date "+%H:%M:%S"; }

echo "[$(ts)] build"
(cd "$REPO/rust/peartree-genotype" && cargo build --release 2>&1 | tail -1)
(cd "$REPO/rust/peartree-genotype2" && cargo build --release 2>&1 | tail -1)
OLD="$REPO/rust/peartree-genotype/target/release/peartree-genotype"
NEW="$REPO/rust/peartree-genotype2/target/release/peartree-genotype2"

for SET in "${SETS[@]}"; do
    OUT="$SP/work/$SET"
    GT="$OUT/genotype"
    case "$SET" in
        e2e) SAMPLES=3 ;;
        *)   SAMPLES=10 ;;
    esac
    WT=$((SAMPLES + 1))
    CON="$GT/P1.genotyping.tprt.txt.gz"
    COMB="$OUT/combine/P1.combined.txt.gz"
    REF="$OUT/ref/reduced.fa"
    REPORT="$GT/bench_v2.md"
    echo "[$(ts)] === $SET: $SAMPLES colonies + WT S$WT ==="
    mkdir -p "$GT/v2" "$GT/legacy_t"
    : > "$GT/timing.tsv"
    for i in $(seq 1 "$WT"); do
        BAM="$OUT/S$i.bam"; [ "$i" = "$WT" ] && BAM="$GT/S$WT.bam"
        t0=$(date +%s.%N)
        "$OLD" --step genotype --bam "$BAM" --insertions "$CON" --out "$GT/legacy_t/S$i.txt.gz" \
            --threads "$THREADS" --config "$REPO/cluster/config.genotype.grch38.tprt" 2>"$GT/legacy_t/S$i.log"
        t1=$(date +%s.%N)
        "$NEW" --step genotype --bam "$BAM" --insertions "$CON" --combined "$COMB" --reference "$REF" \
            --out "$GT/v2/S$i.txt.gz" --threads "$THREADS" ${V2_CFG:+--config "$V2_CFG"} 2>"$GT/v2/S$i.log"
        t2=$(date +%s.%N)
        printf 'S%s\t%.2f\t%.2f\n' "$i" "$(echo "$t1 - $t0" | bc)" "$(echo "$t2 - $t1" | bc)" >> "$GT/timing.tsv"
    done
    {
        echo "# genotype2 benchmark — $SET ($(date +%F))"
        echo
        echo "## Timing (s per colony, $THREADS thread): legacy vs v2"
        echo; echo '```'; cat "$GT/timing.tsv"; echo '```'; echo
        N=$(zcat < "$CON" | grep -c '^>')
        echo "contract loci: $N"
        echo
        echo "## Concordance vs truth (test/e2e/score_genotypes.py)"
        echo
        python3 "$REPO/test/e2e/score_genotypes.py" --e2e-dir "$OUT" --geno-dir "$GT/tprt" --samples "$SAMPLES" --wt-sample "S$WT" --label "(legacy tprt)"
        echo
        python3 "$REPO/test/e2e/score_genotypes.py" --e2e-dir "$OUT" --geno-dir "$GT/v2" --samples "$SAMPLES" --wt-sample "S$WT" --label "(v2)"
        echo
        echo "## score_genotype distribution of present calls (for combine_genotypes min_best_score)"
        echo
        for mode in tprt v2; do
            python3 - "$GT/$mode" <<'PYEOF'
import glob, gzip, os, statistics, sys
d = sys.argv[1]; vals = []
for f in glob.glob(os.path.join(d, 'S*.txt.gz')):
    with gzip.open(f, 'rt') as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if p[1] in ('heterozygous', 'homozygous', 'insertion'):
                vals.append(int(p[2]))
vals.sort()
q = lambda x: vals[min(len(vals) - 1, int(x * len(vals)))] if vals else 'NA'
print(f"- {os.path.basename(d)}: n={len(vals)} min={vals[0] if vals else 'NA'} p5={q(0.05)} p25={q(0.25)} median={q(0.5)} p95={q(0.95)}")
PYEOF
        done
        echo
        if [ -s "$OUT/donor/tree.nwk" ] && [ -s "$OUT/phylo/truth.tsv" ]; then
            echo "## Joint phylogenetic step vs truth"
            echo
            TIPS=$(seq -s, -f 'S%g' 1 "$SAMPLES")
            FILES=(); for i in $(seq 1 "$SAMPLES"); do FILES+=("$GT/v2/S$i.txt.gz"); done
            "$NEW" --step joint --tree "$OUT/donor/tree.nwk" --genotypes "${FILES[@]}" \
                --out "$GT/v2/P1.joint.tsv" --matrix "$GT/v2/P1.joint_matrix.csv.gz" 2>"$GT/v2/joint.log"
            FILESL=(); for i in $(seq 1 "$SAMPLES"); do FILESL+=("$GT/tprt/S$i.txt.gz"); done
            "$NEW" --step joint --tree "$OUT/donor/tree.nwk" --genotypes "${FILESL[@]}" \
                --out "$GT/tprt/P1.joint.tsv" --matrix "$GT/tprt/P1.joint_matrix.csv.gz" 2>"$GT/tprt/joint.log"
            python3 "$DIR/score_joint.py" --truth "$OUT/phylo/truth.tsv" --joint "$GT/v2/P1.joint.tsv" --tips "$TIPS" --label "(v2 genotypes)"
            python3 "$DIR/score_joint.py" --truth "$OUT/phylo/truth.tsv" --joint "$GT/tprt/P1.joint.tsv" --tips "$TIPS" --label "(legacy genotypes)"
            [ -s "$OUT/phylo/fit/phylo_fit.tsv" ] && python3 "$DIR/score_joint.py" --truth "$OUT/phylo/truth.tsv" --python-fit "$OUT/phylo/fit/phylo_fit.tsv" --tips "$TIPS" --label "(legacy genotypes, on disk)"
            echo
            echo "## tree_fit.py on v2 genotypes (legacy read-vote model on n_alt/n_ref)"
            echo
            if "$PY" "$REPO/tools/phylo/tree_fit.py" --genotypes "$OUT/annot/P1.genotypes.csv.gz" --genotype-dir "$GT/v2" \
                    --tree "$OUT/donor/tree.nwk" --annotation "$OUT/annot/P1.annotated.tsv" --out "$GT/v2/fit" --sex F >"$GT/v2/fit.log" 2>&1; then
                python3 "$DIR/score_joint.py" --truth "$OUT/phylo/truth.tsv" --python-fit "$GT/v2/fit/phylo_fit.tsv" --tips "$TIPS" --label "(v2 genotypes)"
            else
                echo "tree_fit.py failed on v2 genotypes (see $GT/v2/fit.log)"
            fi
        fi
    } | tee "$REPORT"
    echo "[$(ts)] -> $REPORT"
done
