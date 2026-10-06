#!/usr/bin/env bash
# =============================================================================
# cluster/genotype2_farm_compare.sh — genotype a patient with peartree-genotype2 on farm22 and
# compare it with the legacy genotyping that already exists for the same patient (the tprt_ab
# kit's run dir: same contract, same combined consensus, same staged BAMs).
#
#   bash cluster/genotype2_farm_compare.sh PD37590              # submit the array + the chained evaluation
#   bash cluster/genotype2_farm_compare.sh PD37590 --dry-run    # resolve every input, print the bsubs, submit nothing
#   bash cluster/genotype2_farm_compare.sh PD37590 --evaluate   # joint step + comparison report (the chained job runs this)
#
# Resolved from files (every value an env override; defaults are farm22 absolute paths):
#   LEGACY_RUNDIR  $TPRT_ROOT/<P>/C_rust/<P>      the legacy pipeline run dir (insertions/, genotypes/, logs/)
#   LEGACY_GT      $LEGACY_RUNDIR/genotypes       legacy per-colony files <colony>.txt.gz
#   CONTRACT       the contract in $LEGACY_RUNDIR/insertions whose locus count matches the legacy files
#                  (<P>.genotyping.txt.gz = two-sided only, or <P>.genotyping.tprt.txt.gz = with one-sided loci)
#   COMBINED       $LEGACY_RUNDIR/insertions/<P>.combined.txt.gz
#   GENOME_2BIT    $JD/hg38.2bit
#   SAMPLES        $TPRT_ROOT/<P>/samples.tsv (sample<TAB>proj, written by cluster/tprt/run_ab.sh)
#   V2_ROOT        $TPRT_ROOT/<P>/V2             everything this script writes: bams.fofn genotypes/ logs/ joint/ report/
#   THROTTLE 12  MEM 4000 (MB; the reservation is what keeps LSF from packing the array onto one node)
#   QUEUE normal  EVAL_MEM 8000  RESULTS_BASE $HOME/results/tprt_ab (NFS copy of the report)
#   GENO2_THREADS 1  threads per genotype task (bsub -n + --threads)
#
# Jobs:  gt2_<P>[1-N]%THROTTLE   genotype_one.sh (GENOTYPE_IMPL=v2; skip-if-exists, atomic)
#        gt2_<P>_eval            ended(array) -> this script --evaluate: joint step (length + uniform branch
#                                prior), cluster/genotype2_compare.py -> $V2_ROOT/report/report.md,
#                                tools/phylo/tree_fit.py cross-check -> $V2_ROOT/fit/summary.md (needs pandas/scipy)
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=cluster/tprt/common.sh
source "$HERE/tprt/common.sh"          # JD, TPRT_ROOT, STAGING_ROOT, RESULTS_BASE, staged_bam, patient_dir, patient_tree, die, note
SELF="$HERE/$(basename "${BASH_SOURCE[0]}")"
PT_ROOT="$PT_ROOT_B"
count_files() { ls "$@" 2>/dev/null | wc -l | tr -d ' ' || true; }   # 0 when the glob matches nothing (pipefail-safe)

P=""; MODE="submit"
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) MODE="dry" ;;
        --evaluate) MODE="evaluate" ;;
        -h|--help) sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*) die "unknown option $1" ;;
        *) [ -z "$P" ] || die "one patient per call (got $P and $1)"; P="$1" ;;
    esac
    shift
done
[ -n "$P" ] || die "usage: genotype2_farm_compare.sh <PATIENT_ID> [--dry-run|--evaluate]"

THROTTLE="${THROTTLE:-12}"; MEM="${MEM:-4000}"; QUEUE="${QUEUE:-normal}"; EVAL_MEM="${EVAL_MEM:-8000}"
THREADS="${GENO2_THREADS:-1}"   # per-colony realignment is CPU-bound on the farm (~30 ms CPU/locus at 30x): -n THREADS + --threads
BIN="${GENOTYPE2_BIN:-$PT_ROOT/rust/peartree-genotype2/target/release/peartree-genotype2}"
CFG="${GENO2_CFG:-$PT_ROOT/cluster/config.genotype2.grch38}"
GENOME_2BIT="${GENOME_2BIT:-$JD/hg38.2bit}"
LEGACY_RUNDIR="${LEGACY_RUNDIR:-$TPRT_ROOT/$P/C_rust/$P}"
LEGACY_GT="${LEGACY_GT:-$LEGACY_RUNDIR/genotypes}"
COMBINED="${COMBINED:-$LEGACY_RUNDIR/insertions/$P.combined.txt.gz}"
SAMPLES="${SAMPLES:-$TPRT_ROOT/$P/samples.tsv}"
V2_ROOT="${V2_ROOT:-$TPRT_ROOT/$P/V2}"
FOFN="$V2_ROOT/bams.fofn"; OUTDIR="$V2_ROOT/genotypes"; LOGS="$V2_ROOT/logs"; JOINT="$V2_ROOT/joint"; REPORT="$V2_ROOT/report"

# --- resolve and check every input ----------------------------------------------------------
[ -x "$BIN" ] || die "no peartree-genotype2 binary at $BIN — build it on the head node: bash $PT_ROOT/cluster/build.sh"
[ -s "$CFG" ] || die "no config $CFG"
[ -s "$GENOME_2BIT" ] || die "reference $GENOME_2BIT missing"
if [ ! -d "$LEGACY_RUNDIR" ]; then
    echo "legacy run dir $LEGACY_RUNDIR does not exist. Candidates under $TPRT_ROOT/$P:" >&2
    ls -d "$TPRT_ROOT/$P"/*/"$P" 2>/dev/null >&2 || echo "  (none)" >&2
    die "set LEGACY_RUNDIR=<one of the above>"
fi
LEG_FILES=("$LEGACY_GT"/*.txt.gz)
[ -s "${LEG_FILES[0]}" ] || die "no legacy genotype files in $LEGACY_GT"
N_LEG_LOCI="$(( $(zcat "${LEG_FILES[0]}" | wc -l | tr -d ' ') - 1 ))"
if [ -z "${CONTRACT:-}" ]; then
    for c in "$LEGACY_RUNDIR/insertions/$P.genotyping.txt.gz" "$LEGACY_RUNDIR/insertions/$P.genotyping.tprt.txt.gz"; do
        [ -s "$c" ] || continue
        n="$(zcat "$c" | grep -c '^>' || true)"
        if [ "$n" -eq "$N_LEG_LOCI" ]; then CONTRACT="$c"; break; fi
        echo "contract $c has $n loci; the legacy files have $N_LEG_LOCI rows — not this one" >&2
    done
    [ -n "${CONTRACT:-}" ] || die "no contract in $LEGACY_RUNDIR/insertions matches the legacy files' $N_LEG_LOCI loci — set CONTRACT="
fi
[ -s "$CONTRACT" ] || die "no contract $CONTRACT"
[ -s "$COMBINED" ] || die "no combined consensus $COMBINED"
TREE="$(patient_tree "$P")"
PDIR="$(patient_dir "$P")"
[ -s "$SAMPLES" ] || die "no $SAMPLES (sample<TAB>proj) — cluster/tprt/run_ab.sh writes it; or: cluster/fleet.sh samples $P > $SAMPLES"
N="$(grep -c . "$SAMPLES" || true)"
[ "${N:-0}" -gt 0 ] || die "$SAMPLES is empty"
LEGACY_ARM="$(basename "$(dirname "$LEGACY_RUNDIR")")"
LEGACY_FIT="${LEGACY_FIT:-$TPRT_ROOT/$P/eval/$LEGACY_ARM/fit/phylo_fit.tsv}"
LEGACY_CALLS="${LEGACY_CALLS:-$LEGACY_RUNDIR/$P.genotypes.csv.gz}"
KNOWN="$PDIR/known_insertions.tsv"

note "patient       : $P   tree $TREE ($(tree_tips "$TREE" | wc -l | tr -d ' ') tips)"
note "legacy run    : $LEGACY_RUNDIR   genotypes $LEGACY_GT (${#LEG_FILES[@]} files, $N_LEG_LOCI loci each)"
note "contract      : $CONTRACT"
note "combined      : $COMBINED"
note "reference     : $GENOME_2BIT"
note "samples       : $SAMPLES ($N colonies)"
note "v2 binary     : $BIN   config $CFG"
note "v2 output     : $V2_ROOT"
[ -s "$LEGACY_CALLS" ] && note "legacy calls  : $LEGACY_CALLS" || note "legacy calls  : none at $LEGACY_CALLS (report without the combine_genotypes section)"
[ -s "$LEGACY_FIT" ] && note "legacy fit    : $LEGACY_FIT" || note "legacy fit    : none at $LEGACY_FIT (report without the tree_fit section)"
[ -s "$KNOWN" ] && note "known loci    : $KNOWN" || note "known loci    : none"

mkdir -p "$V2_ROOT" "$OUTDIR" "$LOGS" "$JOINT" "$REPORT"

# the BAM list: the staged copies the legacy arm genotyped (pipeline.sh bam_path convention)
build_fofn() {
    local missing=0 tmp="$FOFN.tmp.$$"
    : > "$tmp"
    while IFS=$'\t' read -r sample proj _; do
        [ -n "$sample" ] || continue
        local bam; bam="$(staged_bam "$sample" "$proj")"
        if [ ! -s "$bam" ]; then echo "  missing BAM: $bam" >&2; missing=$((missing + 1)); fi
        echo "$bam" >> "$tmp"
    done < "$SAMPLES"
    if [ "$missing" -gt 0 ]; then
        rm -f "$tmp"
        die "$missing staged BAM(s) missing — the legacy run's BAMs were cleaned up? Re-stage with cluster/tprt/run_ab.sh (stage only) or set STAGING_ROOT="
    fi
    mv -f "$tmp" "$FOFN"
    local n_leg_missing=0
    while IFS=$'\t' read -r sample _; do
        [ -n "$sample" ] || continue
        [ -s "$LEGACY_GT/$sample.txt.gz" ] || n_leg_missing=$((n_leg_missing + 1))
    done < "$SAMPLES"
    [ "$n_leg_missing" -eq 0 ] || note "WARNING: $n_leg_missing colonies have no legacy genotype file in $LEGACY_GT (compared on the rest)"
}

ENV_PASS="LEGACY_RUNDIR='$LEGACY_RUNDIR' LEGACY_GT='$LEGACY_GT' CONTRACT='$CONTRACT' COMBINED='$COMBINED' GENOME_2BIT='$GENOME_2BIT' SAMPLES='$SAMPLES' V2_ROOT='$V2_ROOT' GENOTYPE2_BIN='$BIN' GENO2_CFG='$CFG' LEGACY_FIT='$LEGACY_FIT' LEGACY_CALLS='$LEGACY_CALLS' TPRT_ROOT='$TPRT_ROOT' STAGING_ROOT='$STAGING_ROOT' PATIENTS_DIR='$PATIENTS_DIR' RESULTS_BASE='$RESULTS_BASE'"
GT_CMD="FOFN='$FOFN' OUTDIR='$OUTDIR' CONTRACT='$CONTRACT' COMBINED='$COMBINED' GENOME_2BIT='$GENOME_2BIT' GENOTYPE_IMPL=v2 GENOTYPE2_BIN='$BIN' GENO2_CFG='$CFG' GENO2_THREADS='$THREADS' bash '$PT_ROOT/cluster/genotype_one.sh' \$LSB_JOBINDEX"
EVAL_CMD="$ENV_PASS bash '$SELF' '$P' --evaluate"

# --- submit / dry-run -------------------------------------------------------------------------
if [ "$MODE" = submit ] || [ "$MODE" = dry ]; then
    build_fofn
    n_done="$(count_files "$OUTDIR"/*.txt.gz)"
    note "fofn          : $FOFN ($N BAMs; $n_done already genotyped -> skipped by genotype_one.sh)"
    ARRAY=(bsub -J "gt2_${P}[1-${N}]%${THROTTLE}" -o "$LOGS/gt.%I.log" -e "$LOGS/gt.%I.err" -n 1 -q "$QUEUE"
           -R "select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" -M "$MEM" "$GT_CMD")
    if [ "$MODE" = dry ]; then
        echo; echo "DRY-RUN:"; printf ' %q' "${ARRAY[@]}"; echo
        printf ' %q' bsub -J "gt2_${P}_eval" -w 'ended(<the array job id>)' -o "$LOGS/eval.%J.log" -e "$LOGS/eval.%J.err" -n 1 -q "$QUEUE" \
            -R "select[mem>${EVAL_MEM}] rusage[mem=${EVAL_MEM}] span[hosts=1]" -M "$EVAL_MEM" "$EVAL_CMD"; echo
        exit 0
    fi
    out="$("${ARRAY[@]}")"; echo "$out"
    JID="$(echo "$out" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+' | head -1)"
    [ -n "$JID" ] || die "could not parse the array job id from bsub's output"
    out2="$(bsub -J "gt2_${P}_eval" -w "ended($JID)" -o "$LOGS/eval.%J.log" -e "$LOGS/eval.%J.err" -n 1 -q "$QUEUE" \
            -R "select[mem>${EVAL_MEM}] rusage[mem=${EVAL_MEM}] span[hosts=1]" -M "$EVAL_MEM" "$EVAL_CMD")"; echo "$out2"
    JID2="$(echo "$out2" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+' | head -1)"
    printf 'genotype\t%s\nevaluate\t%s\n' "$JID" "${JID2:-?}" > "$V2_ROOT/jobids.tsv"
    echo
    echo "watch:   bjobs -A $JID ; tail -f $LOGS/gt.1.err"
    echo "done:    ls $OUTDIR/*.txt.gz | wc -l    # expect $N"
    echo "report:  $REPORT/report.md   (also $RESULTS_BASE/$P/genotype2_vs_legacy.md)"
    echo "re-run the evaluation by hand:  bash $SELF $P --evaluate"
    exit 0
fi

# --- evaluate -------------------------------------------------------------------------------
n_out="$(count_files "$OUTDIR"/*.txt.gz)"
note "genotype files: $n_out / $N"
if [ "$n_out" -lt "$N" ]; then
    note "missing colonies:"
    while IFS=$'\t' read -r sample _; do
        [ -n "$sample" ] && [ ! -s "$OUTDIR/$sample.txt.gz" ] && echo "  $sample" >&2
    done < "$SAMPLES" || true
    fails="$(grep -lE 'TERM_MEMLIMIT|TERM_RUNLIMIT|Exited with exit code' "$LOGS"/gt.*.log 2>/dev/null || true)"
    [ -z "$fails" ] || { note "failed tasks (LSF report):"; echo "$fails" | sed 's/^/  /' >&2; }
    note "re-submit the missing ones with:  bash $SELF $P   (genotype_one.sh skips finished colonies)"
fi
[ "$n_out" -gt 0 ] || die "nothing to evaluate"

note "joint step (length branch prior)"
"$BIN" --step joint --tree "$TREE" --genotype-dir "$OUTDIR" \
    --out "$JOINT/$P.joint.tsv" --matrix "$JOINT/$P.joint_matrix.csv.gz" 2> "$JOINT/joint.log"
tail -3 "$JOINT/joint.log" >&2 || true
note "joint step (uniform branch prior)"
"$BIN" --step joint --tree "$TREE" --genotype-dir "$OUTDIR" --branch-prior uniform \
    --out "$JOINT/$P.joint.uniform.tsv" --matrix "$JOINT/$P.joint_matrix.uniform.csv.gz" 2> "$JOINT/joint.uniform.log"

# python: the tprt kit's venv if it exists, else the system python3 (the report script is stdlib-only)
PY="${PY:-}"
[ -n "$PY" ] || { [ -x "$TPRT_ROOT/PEAR-TREE/venv/bin/python" ] && PY="$TPRT_ROOT/PEAR-TREE/venv/bin/python"; }
[ -n "$PY" ] || PY="$(command -v python3 || true)"
[ -n "$PY" ] || die "no python3"

run_report() {   # <legacy genotype dir> <out dir>
    local gt="$1" out="$2"
    mkdir -p "$out"
    local extra=()
    [ -s "$LEGACY_CALLS" ] && extra+=(--legacy-calls "$LEGACY_CALLS")
    [ -s "$LEGACY_FIT" ] && extra+=(--legacy-fit "$LEGACY_FIT")
    [ -s "$KNOWN" ] && extra+=(--known "$KNOWN")
    [ -d "$LEGACY_RUNDIR/logs" ] && extra+=(--legacy-logs "$LEGACY_RUNDIR/logs")
    "$PY" "$PT_ROOT/cluster/genotype2_compare.py" --patient "$P" --legacy-dir "$gt" --v2-dir "$OUTDIR" \
        --joint "$JOINT/$P.joint.tsv" --joint-uniform "$JOINT/$P.joint.uniform.tsv" --tree "$TREE" \
        --v2-logs "$LOGS" --out-dir "$out" "${extra[@]}" > "$out/report.stdout" 2> "$out/report.stderr" \
        || { cat "$out/report.stderr" >&2; die "genotype2_compare.py failed for $gt"; }
    [ -s "$out/report.stderr" ] && cat "$out/report.stderr" >&2
    note "report: $out/report.md"
}
run_report "$LEGACY_GT" "$REPORT"
# other legacy genotype sets kept next to the main one (e.g. genotypes.het30 = vaf_het_min 0.30): one report each
for d in "$LEGACY_RUNDIR"/genotypes.*; do
    [ -d "$d" ] || continue
    [ "$d" = "$LEGACY_GT" ] && continue
    ls "$d"/*.txt.gz >/dev/null 2>&1 || continue
    run_report "$d" "$REPORT/$(basename "$d")"
done

mkdir -p "$RESULTS_BASE/$P"
cp -f "$REPORT/report.md" "$RESULTS_BASE/$P/genotype2_vs_legacy.md"
[ -s "$REPORT/known.tsv" ] && cp -f "$REPORT/known.tsv" "$RESULTS_BASE/$P/genotype2_known.tsv"
cp -f "$JOINT/$P.joint.tsv" "$RESULTS_BASE/$P/" 2>/dev/null || true
note "copied to $RESULTS_BASE/$P/genotype2_vs_legacy.md"

# independent Python cross-check of the Rust joint step: tools/phylo/tree_fit.py (the read-vote model)
# on the v2 per-colony files + the numeric matrix; it joins <P>.joint.tsv beside the matrix and writes a
# tree_fit class x joint class table into summary.md. Needs pandas/scipy (the tprt kit venv has them).
FIT="$V2_ROOT/fit"
if "$PY" -c 'import pandas, scipy' 2>/dev/null; then
    note "tree_fit cross-check -> $FIT"
    rm -rf "$FIT"; mkdir -p "$FIT"
    if (cd "$PT_ROOT" && "$PY" tools/phylo/tree_fit.py --genotypes "$JOINT/$P.joint_matrix.csv.gz" \
            --genotype-dir "$OUTDIR" --tree "$TREE" --out "$FIT") > "$V2_ROOT/tree_fit.log" 2>&1; then
        cp -f "$FIT/summary.md" "$RESULTS_BASE/$P/genotype2_tree_fit_summary.md" 2>/dev/null || true
        sed -n '/^## Cross-check/,/^How to read/p' "$FIT/summary.md" >&2 || true
    else
        note "tree_fit FAILED (see $V2_ROOT/tree_fit.log) -- the report above stands without it"
        tail -5 "$V2_ROOT/tree_fit.log" >&2 || true
    fi
else
    note "$PY lacks pandas/scipy: tree_fit cross-check skipped (PY=<python with pandas+scipy> to enable)"
fi
sed -n '1,/^## Per locus/p' "$REPORT/report.md"
