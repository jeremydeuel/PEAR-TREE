#!/usr/bin/env bash
# =============================================================================
# cluster/tprt/evaluate.sh — phylogeny-based evaluation of both arms of one patient + A/B report.
#
#   bash cluster/tprt/evaluate.sh PD37449         # (run_ab.sh submits this after both annotates)
#
# Per arm with outputs ($TPRT_ROOT/<P>/<arm>/<P>/):
#   tools/phylo/tree_fit.py        --genotypes <P>.genotypes.csv.gz --genotype-dir genotypes/
#                                  --tree <patients/.../tree> --annotation <P>.annotated.csv.gz
#                                  --samples colonies.tsv --sex $SEX   -> eval/<arm>/fit/phylo_fit.tsv
#   tools/phylo/discrimination.py  --fit eval/<arm>/fit/phylo_fit.tsv  -> eval/<arm>/discrimination/
# then cluster/tprt/compare_arms.py -> eval/ab_report.md (+ a_only.tsv, b_only.tsv, matched*.tsv,
# stage_counts.tsv, lsf_resources.tsv), all copied to $RESULTS_BASE/<P>/.
#
# Idempotent: re-running recomputes the evaluation (cheap) from the arms' finished outputs. An arm
# without calls is skipped (one-arm report). SEX (default auto) is passed to tree_fit.
# =============================================================================
set -euo pipefail
# shellcheck source=cluster/tprt/common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

P="${1:?usage: evaluate.sh <PATIENT_ID>}"
SEX="${SEX:-auto}"
PY="$VENV/bin/python"
[ -x "$PY" ] || die "no venv python $PY"
TREE="$(patient_tree "$P")"
COLS="$(patient_colonies "$P")"
EV="$(eval_dir "$P")"
mkdir -p "$EV"

NHAVE=0; NFAILED=0
for ARM in "${ALL_ARMS[@]}"; do
    RD="$(arm_rundir "$P" "$ARM")"
    GT="$RD/$P.genotypes.csv.gz"; AN="$RD/$P.annotated.csv.gz"
    if [ ! -s "$GT" ] || [ ! -s "$AN" ]; then
        note "arm $ARM: no finished outputs ($GT / $AN) — skipped"
        continue
    fi
    NHAVE=$((NHAVE + 1))
    FIT="$EV/$ARM/fit"; DIS="$EV/$ARM/discrimination"
    rm -rf "$FIT" "$DIS"          # never report a previous run's labels if this tree_fit fails
    mkdir -p "$FIT" "$DIS"
    note "arm $ARM: tree_fit"
    if (cd "$PT_ROOT_B" && "$PY" tools/phylo/tree_fit.py \
            --genotypes "$GT" --genotype-dir "$RD/genotypes" --tree "$TREE" \
            --annotation "$AN" --out "$FIT" --sex "$SEX" --samples "$COLS") \
            > "$EV/$ARM/tree_fit.log" 2>&1; then
        note "arm $ARM: discrimination"
        (cd "$PT_ROOT_B" && "$PY" tools/phylo/discrimination.py --fit "$FIT/phylo_fit.tsv" --out "$DIS") \
            > "$EV/$ARM/discrimination.log" 2>&1 \
            || { note "arm $ARM: discrimination FAILED (see $EV/$ARM/discrimination.log) — report continues without it"; NFAILED=$((NFAILED + 1)); }
    else
        note "arm $ARM: tree_fit FAILED (see $EV/$ARM/tree_fit.log) — report continues without phylo labels"; NFAILED=$((NFAILED + 1))
        tail -5 "$EV/$ARM/tree_fit.log" >&2 || true
    fi
done
[ "$NHAVE" -gt 0 ] || die "neither arm of $P has finished outputs under $TPRT_ROOT/$P"

HAVE_A=0; [ -s "$(arm_rundir "$P" A)/$P.annotated.csv.gz" ] && HAVE_A=1
if [ "$HAVE_A" = 1 ] && [ -s "$(arm_rundir "$P" B)/$P.annotated.csv.gz" ]; then
    note "compare_arms A vs B"
    "$PY" "$TPRT_KIT_DIR/compare_arms.py" --patient "$P" \
        --rundir-a "$(arm_rundir "$P" A)" --rundir-b "$(arm_rundir "$P" B)" \
        --eval-a "$EV/A" --eval-b "$EV/B" --tree "$TREE" --out-dir "$EV"
else
    note "arm A or B has no annotated output — A vs B report skipped"
fi
if [ "$HAVE_A" = 1 ] && [ -s "$(arm_rundir "$P" C)/$P.annotated.csv.gz" ]; then
    note "compare_arms A vs C"
    mkdir -p "$EV/AC" "$RESULTS_BASE/$P"
    "$PY" "$TPRT_KIT_DIR/compare_arms.py" --patient "$P" --label-b C \
        --rundir-a "$(arm_rundir "$P" A)" --rundir-b "$(arm_rundir "$P" C)" \
        --eval-a "$EV/A" --eval-b "$EV/C" --tree "$TREE" --out-dir "$EV/AC"
    cp -f "$EV/AC/ac_report.md" "$RESULTS_BASE/$P/" 2>/dev/null || true
fi

mkdir -p "$RESULTS_BASE/$P"
# known insertions (patients/<organ>/<P>/known_insertions.tsv): trace each through every stage, per arm
KNOWN="$(patient_dir "$P")/known_insertions.tsv"
if [ -s "$KNOWN" ]; then
    for ARM in "${ALL_ARMS[@]}"; do
        RD="$(arm_rundir "$P" "$ARM")"
        [ -d "$RD/discovery" ] || continue
        note "arm $ARM: known insertions"
        if "$PY" "$TPRT_KIT_DIR/check_known.py" --known "$KNOWN" --rundir "$RD" --eval "$EV/$ARM" \
                --label "$ARM" --out-dir "$EV/$ARM/known" > "$EV/$ARM/known.log" 2>&1; then
            cp -f "$EV/$ARM/known/known_report.md" "$RESULTS_BASE/$P/known_report.$ARM.md" 2>/dev/null || true
            sed -n '1,/^## Potentially/p' "$EV/$ARM/known/known_report.md"
        else
            note "arm $ARM: check_known FAILED (see $EV/$ARM/known.log)"; NFAILED=$((NFAILED + 1))
        fi
    done
fi
cp -f "$EV"/ab_report.md "$EV"/*.tsv "$RESULTS_BASE/$P/" 2>/dev/null || true
for ARM in "${ALL_ARMS[@]}"; do
    if [ -s "$EV/$ARM/fit/phylo_fit.tsv" ]; then cp -f "$EV/$ARM/fit/phylo_fit.tsv" "$RESULTS_BASE/$P/phylo_fit.$ARM.tsv"; fi
done
note "reports: $EV/ab_report.md, $EV/AC/ac_report.md (whichever exist; copied to $RESULTS_BASE/$P/)"
if [ "$NFAILED" -gt 0 ]; then
    note "$NFAILED phylo step(s) FAILED — the report has no phylo labels for that arm; fix, then re-run this script"
    exit 2
fi
