#!/usr/bin/env bash
# combine_genotypes over the POOLED 9x10 cohort -> the final call table.
# Runs INSIDE the LSF job; cluster/submit_combine_genotypes.sh dispatches it.
#
# THIS is where the cohort-level gates live (min_wild-types, min_dispersion, max_artefact,
# max_na -- src/config.py's combine_genotypes block), and pooling changes what they MEAN.
# In the single-donor PD44579 run a germline MEI was present in ~all colonies; here it is
# present in ~10 of 90 (one donor's worth) and WILD-TYPE in the other 80. So min_wild-types
# is now a CROSS-DONOR filter, not an allele-fraction one. None of these thresholds are
# validated for this cohort -- they are inherited from 174 colonies of one donor on hs37d5.
# Do not report the first call table's numbers as if the gates were tuned for it.
#
# Env: FOFN, GTDIR, OUT, THREADS, VENV, ALLOW_MISSING
set -euo pipefail
cd "$(dirname "$0")/.."
PT_ROOT="$PWD"

FOFN="${FOFN:-$HOME/catalogue/picked/all.bams.fofn}"
GTDIR="${GTDIR:-genotypes_grch38}"
OUT="${OUT:-mei9x10.genotypes.csv.gz}"
THREADS="${THREADS:-8}"
VENV="${VENV:-/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE/venv}"
ALLOW_MISSING="${ALLOW_MISSING:-0}"

[ -s "$FOFN" ] || { echo "no fofn: $FOFN" >&2; exit 1; }
[ -x "$VENV/bin/python" ] || { echo "no venv python at $VENV/bin/python" >&2; exit 1; }
if [ -s "$OUT" ]; then echo "call table exists, skipping: $OUT"; exit 0; fi

# Resolve genotype files from the FOFN, not a glob. A glob over GTDIR would silently pull in
# any stray colony left there by an earlier run, and in a POOLED table that contamination is
# invisible -- an extra genome just looks like more cohort, and it lands directly on the
# cross-donor gates the FP readout depends on.
files=(); missing=0
while read -r bam; do
    [ -n "$bam" ] || continue
    id="$(basename "$bam")"; id="${id%.bam}"; id="${id%.cram}"; id="${id%.sample.dupmarked}"
    f="$GTDIR/$id.txt.gz"
    if [ -s "$f" ]; then files+=("$f"); else echo "MISSING genotype: $f" >&2; missing=$((missing+1)); fi
done < "$FOFN"
echo "pooled cohort: ${#files[@]} genotype files, $missing missing (allowed: $ALLOW_MISSING)"
[ "${#files[@]}" -gt 0 ] || { echo "no genotype files — aborting" >&2; exit 1; }
[ "$missing" -le "$ALLOW_MISSING" ] || {
    echo "too many missing genotype files ($missing > $ALLOW_MISSING) — refusing to combine." >&2
    echo "A colony absent here is silently absent from every locus's wild-type count, which is" >&2
    echo "exactly what min_wild-types and min_dispersion gate on. Re-run the missing tasks." >&2
    exit 1; }

echo "combining ${#files[@]} colonies -> $OUT"
"$VENV/bin/python" "$PT_ROOT/src/main.py" --step combine_genotypes \
    --genotypes "${files[@]}" --out "$OUT" --threads "$THREADS"

[ -s "$OUT" ] || { echo "combine_genotypes produced no $OUT" >&2; exit 1; }
echo "calls: $OUT ($(zcat "$OUT" | tail -n +2 | wc -l) loci passing the cohort gates, over ${#files[@]} colonies)"
