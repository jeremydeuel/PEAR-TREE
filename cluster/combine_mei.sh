#!/usr/bin/env bash
# combine_insertions for the 10x10 MEI benchmark — ONE POOLED COHORT over all 50 colonies.
#
# WHY POOLED, not per-patient. The five patients are unrelated genomes, so it is tempting
# to build five contracts. That is wrong, and it breaks the benchmark:
#
#   Genotyping only ever visits loci that are IN the contract. With per-patient contracts,
#   a locus discovered in PD34200 but not PD37449 is never genotyped in PD37449's colonies
#   — so its absence there is uninterpretable: "artefact confined to one patient" and
#   "never tested outside that patient" produce identical output. Cross-patient recurrence
#   is THE readout of this benchmark, and it needs a common locus set across all 50.
#
# Pooling is what creates the comparison; it does not destroy patient independence. The
# independence is ours, held outside the algorithm in the patient->colony map. The caller
# is deliberately blind to it, and that blindness is the leverage:
#   - a somatic call shared across patients      -> artefact (unrelated genomes)
#   - a germline MEI shared across patients      -> expected (common polymorphism)
#   - a somatic call confined to one clade       -> candidate true positive
#
# Safe to pool statistically: the cohort-level gates (min_wild-types, min_dispersion,
# max_artefact, max_na) live in combine_genotypes, NOT here. The only cohort-ish knob in
# combine_insertions is exclude_files_with_many_insertions, which is per-file. Cohort size
# does not move this step's behaviour.
#
# Usage:  cluster/combine_mei.sh
# Env:
#   FOFNDIR   per-patient fofns + all.bams.fofn  (default: $HOME/mei10x10)
#   DISCDIR   discovery output                   (default: discovery_mei10x10)
#   OUTDIR    contract dir                       (default: insertions_mei10x10)
#   STEM      contract basename                  (default: mei10x10)
#   PATIENTS  space-separated list               (used for the per-patient assembly gate)
#   THREADS   combine cores                      (default: 16)
#   VENV      python venv (NOT in this run dir)  (default: the sibling PEAR-TREE checkout)
#   ALLOW_MISSING  tolerate N absent discovery files (default: 1 — PD41048b_lo0015)
#   COMBINE_IMPL   python (default) | rust (rust/peartree-combine, built by cluster/build.sh)
#   COMBINE_BIN    rust binary (default: rust/peartree-combine/target/release/peartree-combine)
set -euo pipefail
cd "$(dirname "$0")/.."
PT_ROOT="$PWD"

FOFNDIR="${FOFNDIR:-$HOME/mei10x10}"
DISCDIR="${DISCDIR:-discovery_mei10x10}"
OUTDIR="${OUTDIR:-insertions_mei10x10}"
STEM="${STEM:-mei10x10}"
THREADS="${THREADS:-16}"
VENV="${VENV:-/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE/venv}"
PATIENTS="${PATIENTS:-PD34200 PD37449 PD43947 PD41048 PD43974}"
ALLOW_MISSING="${ALLOW_MISSING:-1}"
COMBINE_IMPL="${COMBINE_IMPL:-python}"
COMBINE_BIN="${COMBINE_BIN:-$PT_ROOT/rust/peartree-combine/target/release/peartree-combine}"
case "$COMBINE_IMPL" in
    python) ;;
    rust) [ -x "$COMBINE_BIN" ] || { echo "COMBINE_IMPL=rust but missing: $COMBINE_BIN (run cluster/build.sh)" >&2; exit 1; } ;;
    *) echo "COMBINE_IMPL must be python or rust, got '$COMBINE_IMPL'" >&2; exit 1 ;;
esac

FOFN="$FOFNDIR/all.bams.fofn"
[ -s "$FOFN" ] || { echo "no pooled fofn: $FOFN (run cluster/build_fofn.sh)" >&2; exit 1; }
mkdir -p "$OUTDIR"

CONTRACT="$OUTDIR/$STEM.genotyping.txt.gz"
if [ -s "$CONTRACT" ]; then echo "contract exists, skipping: $CONTRACT"; exit 0; fi

# Resolve discovery files from the fofn, not a glob: a glob over DISCDIR would silently
# pull in any stray colony left there by an earlier run, and in a POOLED contract that
# contamination is invisible — an extra genome just looks like more cohort.
files=(); missing=0
while read -r bam; do
    [ -n "$bam" ] || continue
    id="$(basename "$bam")"; id="${id%.bam}"; id="${id%.cram}"; id="${id%.sample.dupmarked}"
    f="$DISCDIR/$id.txt.gz"
    if [ -s "$f" ]; then files+=("$f"); else echo "MISSING discovery: $f" >&2; missing=$((missing+1)); fi
done < "$FOFN"
echo "pooled cohort: ${#files[@]} discovery files, $missing missing (allowed: $ALLOW_MISSING)"
[ "${#files[@]}" -gt 0 ] || { echo "no discovery files — aborting" >&2; exit 1; }
[ "$missing" -le "$ALLOW_MISSING" ] || {
    echo "too many missing discovery files ($missing > $ALLOW_MISSING) — refusing to combine." >&2
    echo "A colony absent from the contract is absent from EVERY genotype, silently." >&2
    exit 1; }

# ASSEMBLY GATE — combine_insertions is assembly-specific (reads flanks from genome_2bit,
# lifts hs1 clip hits back through a chain). A mismatched config yields plausible, silently
# WRONG coordinates and no error.
#
# Check ONE BAM PER PATIENT, not just the first. nst_links assembly varies per BAM — even
# within a project (1903 holds hs37d5 and GRCh38 side by side). Pooled, a single
# wrong-assembly patient corrupts the shared contract for all 50 colonies, not just its own.
for P in $PATIENTS; do
    PF="$FOFNDIR/$P.bams.fofn"
    [ -s "$PF" ] || { echo "no fofn for $P: $PF" >&2; exit 1; }
    BAM1="$(head -1 "$PF")"
    [ -s "$BAM1" ] || { echo "[$P] first BAM unreadable: $BAM1" >&2; exit 1; }
    echo "[$P] verifying config against $BAM1"
    # install.sh defaults VENV to $PT_ROOT/venv, which is WRONG here: this run dir has no
    # venv — it lives in the sibling PEAR-TREE checkout. Pass ours through explicitly or
    # check-config dies with "no venv python at .../venv/bin/python" and we lose the
    # assembly gate to a path bug rather than to a real mismatch.
    VENV="$VENV" bash "$PT_ROOT/cluster/install.sh" check-config --bam "$BAM1" --pt-root "$PT_ROOT" \
        || { echo "[$P] assembly/config check FAILED — refusing to combine" >&2; exit 1; }
done

echo "combining ${#files[@]} colonies into one contract -> $OUTDIR/$STEM"
if [ "$COMBINE_IMPL" = rust ]; then
    echo "combine_insertions: rust ($COMBINE_BIN), config $PT_ROOT/src/config.py"
    PEARTREE_PYTHON="$VENV/bin/python" "$COMBINE_BIN" --step combine_insertions --config "$PT_ROOT/src/config.py" \
        --discovery_files "${files[@]}" --out "$OUTDIR/$STEM" --threads "$THREADS"
else
    "$VENV/bin/python" "$PT_ROOT/src/main.py" --step combine_insertions \
        --discovery_files "${files[@]}" --out "$OUTDIR/$STEM" --threads "$THREADS"
fi

[ -s "$CONTRACT" ] || { echo "combine_insertions produced no $CONTRACT" >&2; exit 1; }
echo "contract: $CONTRACT ($(zcat "$CONTRACT" | grep -c '^>') loci) over ${#files[@]} colonies"
