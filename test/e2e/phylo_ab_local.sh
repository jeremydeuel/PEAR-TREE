#!/bin/bash
# Local proof of the farm evaluation chain (cluster/tprt/evaluate.sh -> tools/phylo/tree_fit.py ->
# tools/phylo/discrimination.py -> cluster/tprt/compare_arms.py) on a tree-mode E2E output
# (test/e2e/run_phylo_e2e.sh). Two stand-in arms from the same simulated patient:
#   A = legacy genotyper config + legacy contract   ($E2E/genotype/legacy)
#   B = .tprt genotyper config + extended contract  ($E2E/genotype/tprt)
# laid out as the kit's $TPRT_ROOT/<P>/<arm>/<P>/ rundirs (symlinks; genotypes.csv.gz = the raw
# per-colony calls, ungated), with patients/sim/<P>/{<P>.tree, colonies.tsv} from the simulator.
#
# Usage: bash test/e2e/phylo_ab_local.sh     (E2E=<run_phylo_e2e.sh OUT>, AB=<scratch layout>)
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$DIR/../.." && pwd)"
SP="${SP:-/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad}"
E2E="${E2E:-$SP/work/e2e_phylo}"
AB="${AB:-$SP/work/tprt_ab_local}"
P="${P:-SIMP1}"
PY="${PY:-$SP/venv/bin/python}"
rm -rf "$AB"
mkdir -p "$AB/patients/sim/$P"
cp "$E2E/donor/tree.nwk" "$AB/patients/sim/$P/$P.tree"
{ printf 'donor\tproj\tds\treadlen\tmapped\tassembly\tsample\n'
  awk -F'\t' -v p="$P" 'NR > 1 { printf "%s\t0\tWGS\t151\tNA\tGRCh38\t%s\n", p, $2 }' "$E2E/donor/samples.tsv"; } \
  > "$AB/patients/sim/$P/colonies.tsv"
for ARM in A B; do
    GM=legacy; CON="$E2E/combine/P1.genotyping.txt.gz"
    [ "$ARM" = B ] && { GM=tprt; CON="$E2E/genotype/P1.genotyping.tprt.txt.gz"; }
    RD="$AB/$P/$ARM/$P"
    mkdir -p "$RD/genotypes" "$RD/insertions" "$RD/discovery"
    for f in "$E2E/genotype/$GM"/S*.txt.gz; do ln -s "$f" "$RD/genotypes/$(basename "$f")"; done
    for f in "$E2E"/S*.discovery.txt.gz; do ln -s "$f" "$RD/discovery/$(basename "$f" .discovery.txt.gz).txt.gz"; done
    ln -s "$E2E/combine/P1.combined.txt.gz" "$RD/insertions/$P.combined.txt.gz"
    ln -s "$CON" "$RD/insertions/$P.genotyping$([ "$ARM" = B ] && echo .tprt).txt.gz"
    # annotate's table, renamed to the patient (tab-separated, as annotate_v2 writes it)
    gzip -c "$E2E/annot/P1.annotated.tsv" > "$RD/$P.annotated.csv.gz"
    # calls matrix in combine_genotypes' format (';', insertion;<colony>...), ungated
    "$PY" - "$RD/genotypes" "$RD/$P.genotypes.csv.gz" <<'EOF'
import sys, os, pandas as pd
d, out = sys.argv[1], sys.argv[2]
cols = {}
for f in sorted(os.listdir(d)):
    cols[f.split(".")[0]] = pd.read_csv(os.path.join(d, f), sep="\t", index_col=0)["genotype"]
df = pd.DataFrame(cols)
df.index.name = "insertion"
df.to_csv(out, sep=";")
EOF
done
echo "layout -> $AB"
TPRT_ROOT="$AB" PATIENTS_DIR="$AB/patients" VENV="$SP/venv" RESULTS_BASE="$AB/results" SEX=F \
    bash "$REPO/cluster/tprt/evaluate.sh" "$P"
