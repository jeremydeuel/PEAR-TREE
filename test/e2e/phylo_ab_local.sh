#!/bin/bash
# Local proof of the farm evaluation chain (cluster/tprt/evaluate.sh -> tools/phylo/tree_fit.py ->
# tools/phylo/discrimination.py -> cluster/tprt/compare_arms.py -> cluster/tprt/check_known.py) on a
# tree-mode E2E output (test/e2e/run_phylo_e2e.sh). Stand-in arms from the same simulated patient:
#   A = legacy genotyper config + legacy contract   ($E2E/genotype/legacy)
#   B = .tprt genotyper config + extended contract  ($E2E/genotype/tprt)
#   C = rust/peartree-genotype2 (numeric per-colony files $E2E/genotype/v2 + the joint step's
#       $E2E/genotype/v2_joint/P1.joint_matrix.csv.gz and P1.joint.tsv; test/genotype2/bench.sh),
#       only when those exist -- the GENOTYPE_IMPL=v2 run-dir layout
# laid out as the kit's $TPRT_ROOT/<P>/<arm>/<P>/ rundirs (symlinks; for A / B genotypes.csv.gz = the
# raw per-colony calls, ungated; for C the numeric P(carrier) matrix + <P>.joint.tsv), with
# patients/sim/<P>/{<P>.tree, colonies.tsv, known_insertions.tsv} from the simulator (known set =
# the truth clade / private loci of phylo/truth.tsv, when present).
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
if [ -s "$E2E/phylo/truth.tsv" ]; then
    { printf 'locus\ttier\tclass\tsubclass\tn_carriers\tcarriers\tnote\n'
      awk -F'\t' 'NR > 1 && ($2 == "clade" || $2 == "private") {
          n = split($3, c, ","); printf "%s\tknown\tsim\t%s\t%d\t%s\t%s\n", $1, $5, n, $3, $4 }' "$E2E/phylo/truth.tsv"; } \
      > "$AB/patients/sim/$P/known_insertions.tsv"
fi
ARMS="A B"
[ -s "$E2E/genotype/v2_joint/P1.joint.tsv" ] && [ -s "$E2E/genotype/v2/S1.txt.gz" ] && ARMS="A B C"
for ARM in $ARMS; do
    GM=legacy; CON="$E2E/combine/P1.genotyping.txt.gz"
    [ "$ARM" = B ] && { GM=tprt; CON="$E2E/genotype/P1.genotyping.tprt.txt.gz"; }
    [ "$ARM" = C ] && { GM=v2; CON="$E2E/genotype/P1.genotyping.tprt.txt.gz"; }
    RD="$AB/$P/$ARM/$P"
    mkdir -p "$RD/genotypes" "$RD/insertions" "$RD/discovery"
    for f in "$E2E/genotype/$GM"/S*.txt.gz; do ln -s "$f" "$RD/genotypes/$(basename "$f")"; done
    for f in "$E2E"/S*.discovery.txt.gz; do ln -s "$f" "$RD/discovery/$(basename "$f" .discovery.txt.gz).txt.gz"; done
    ln -s "$E2E/combine/P1.combined.txt.gz" "$RD/insertions/$P.combined.txt.gz"
    ln -s "$CON" "$RD/insertions/$P.genotyping$([ "$ARM" != A ] && echo .tprt).txt.gz"
    # annotate's table, renamed to the patient (tab-separated, as annotate_v2 writes it)
    gzip -c "$E2E/annot/P1.annotated.tsv" > "$RD/$P.annotated.csv.gz"
    if [ "$ARM" = C ]; then
        # GENOTYPE_IMPL=v2 phase 4 outputs: numeric P(carrier) matrix + the joint table
        cp "$E2E/genotype/v2_joint/P1.joint_matrix.csv.gz" "$RD/$P.genotypes.csv.gz"
        cp "$E2E/genotype/v2_joint/P1.joint.tsv" "$RD/$P.joint.tsv"
        continue
    fi
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
echo "layout -> $AB (arms $ARMS)"
TPRT_ROOT="$AB" PATIENTS_DIR="$AB/patients" VENV="$SP/venv" RESULTS_BASE="$AB/results" SEX=F \
    bash "$REPO/cluster/tprt/evaluate.sh" "$P"
