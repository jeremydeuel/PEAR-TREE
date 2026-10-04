#!/usr/bin/env bash
# =============================================================================
# cluster/tprt/run_ab.sh — A/B the TPRT-hallmark pipeline against the current one, ONE patient.
#
#   bash cluster/tprt/run_ab.sh PD37449                 # both arms + evaluation (default)
#   bash cluster/tprt/run_ab.sh PD37449 --dry-run       # print every bsub, submit nothing
#   bash cluster/tprt/run_ab.sh PD37449 --arm B         # one arm only (the other from an earlier run)
#   bash cluster/tprt/run_ab.sh PD37449 --cleanup       # also delete the staged BAMs when both arms are done
#
# Options: --arm A|B|both (default both)   --dry-run   --cleanup   --no-eval
#          --skip-preflight (skip the quick preflight)   --force (resubmit although jobs of the
#          last submission of that arm are still pending/running)
#
# Resolves everything from files: patients/<organ>/<P>/ (the tree, colonies.tsv), the sample
# list through `cluster/fleet.sh samples <P>` (the fleet's GRCh38 & WGS rule), the per-sample
# iRODS project (staging = pipeline.sh's stageBam.pl convention). Each arm is one ordinary
# `cluster/pipeline.sh submit-list` DAG (stage+discover -> retry -> combine -> genotype ->
# retry -> combine_genotypes -> annotate), with its own WORKROOT, RESULTS_DIR, configs, job-name
# prefix and src/config.py (arm A runs from the PT_ROOT_A worktree, arm B from this checkout).
#
#   jobs                     depends on
#   <P>_A_*  (pipeline DAG)   -
#   <P>_B_*  (pipeline DAG)   B's stage+discover element i waits for A's element i (ended):
#                             A stages each BAM once, B reuses it (two stageBam.pl writing one
#                             file would corrupt it). Same-size LSF arrays => element-wise.
#   <P>_tprt_eval            done(<P>_A_annotate) && done(<P>_B_annotate) -> evaluate.sh
#   <P>_tprt_cleanup         done(<P>_A_cg) && done(<P>_B_cg)   (only with --cleanup)
#
# The pipelines' own phase-5 cleanup is DISABLED (PT_NO_CLEANUP=1): arm A finishing first would
# otherwise delete the BAMs arm B is still genotyping.
#
# Resources (MB): see the "per-arm resources" block; every value is an env override.
# =============================================================================
set -euo pipefail
# shellcheck source=cluster/tprt/common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

usage() { sed -n '2,34p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit "${1:-1}"; }

P=""; ARMS="both"; DRY=0; CLEANUP=0; EVAL=1; SKIP_PRE=0; FORCE=0
while [ $# -gt 0 ]; do
    case "$1" in
        --arm) ARMS="${2:?--arm needs A|B|both}"; shift 2 ;;
        --arm=*) ARMS="${1#--arm=}"; shift ;;
        --dry-run) DRY=1; shift ;;
        --cleanup) CLEANUP=1; shift ;;
        --no-eval) EVAL=0; shift ;;
        --skip-preflight) SKIP_PRE=1; shift ;;
        --force) FORCE=1; shift ;;
        -h|--help) usage 0 ;;
        -*) die "unknown option $1" ;;
        *) [ -z "$P" ] || die "one patient per call (got $P and $1)"; P="$1"; shift ;;
    esac
done
[ -n "$P" ] || usage
case "$ARMS" in A) RUN_ARMS=(A) ;; B) RUN_ARMS=(B) ;; both) RUN_ARMS=(A B) ;; *) die "--arm must be A, B or both" ;; esac

# --- per-arm resources (MB) -----------------------------------------------------
# discovery: post-jemalloc PD44579 cohort peak max 3.1 GB / p95 2.1 GB (job 954062) -> 8 GB for
#   arm A (the validated submit_discovery.sh default). Arm B adds the evidence sidecar,
#   fetch_all_mates, SHORT reads and a per-sample floor of ONE fragment (many more emitted
#   breakpoints held with full records) — unmeasured on real WGS, so 12 GB; the retry controller
#   escalates any TERM_MEMLIMIT once to 32 GB on `long`.
# combine_insertions: legacy 15.5 GB peak on 174 PD44579 files -> 24 GB (A). Arm B pools every
#   colony's sidecar reads per junction in Python (dedup, SHORT, mappy) — unmeasured, and combine
#   has NO retry controller, so 64 GB. Lower both after the pilot (ab_report.md prints the peaks).
# genotype: 140 MB measured; 4 GB reservation on purpose — it is what keeps LSF from packing the
#   whole array onto one 1.9 TB node (I/O starvation, genotyping-perf note); THROTTLE 12 per arm
#   (24 concurrent for both) for the same reason.
A_SD_MEM="${A_SD_MEM:-8000}";   B_SD_MEM="${B_SD_MEM:-12000}"; SD_MEM_T2="${SD_MEM_T2:-32000}"
A_CI_MEM="${A_CI_MEM:-24000}";  B_CI_MEM="${B_CI_MEM:-64000}"; CI_CORES="${CI_CORES:-8}"
GT_MEM_T1="${GT_MEM_T1:-4000}"; GT_MEM_T2="${GT_MEM_T2:-8000}"; GT_THROTTLE="${GT_THROTTLE:-12}"
STAGE_THROTTLE="${STAGE_THROTTLE:-20}"
AN_MEM="${AN_MEM:-32000}"; AN_CORES="${AN_CORES:-4}"
EV_MEM="${EV_MEM:-16000}"; EV_CORES="${EV_CORES:-2}"
QUEUE="${QUEUE:-normal}"

# --- resolve the patient from files ----------------------------------------------
PDIR="$(patient_dir "$P")"
TREE="$(patient_tree "$P")"
COLS="$(patient_colonies "$P")"
[ -s "$COLS" ] || die "no $COLS"
PD="$(pat_dir "$P")"
mkdir -p "$PD"
SAMPLES_TSV="$PD/samples.tsv"
patient_samples "$P" > "$SAMPLES_TSV.tmp.$$" || die "fleet.sh samples $P failed"
mv -f "$SAMPLES_TSV.tmp.$$" "$SAMPLES_TSV"
N="$(wc -l < "$SAMPLES_TSV" | tr -d ' ')"
[ "$N" -gt 0 ] || die "$P has no GRCh38 WGS colony in $COLS"
NTIPS="$(tree_tips "$TREE" | wc -l | tr -d ' ')"
NMISS="$(comm -23 <(tree_tips "$TREE") <(cut -f1 "$SAMPLES_TSV" | sort -u) | wc -l | tr -d ' ')"

note "patient  : $P  ($PDIR)"
note "tree     : $TREE  ($NTIPS tips; $NMISS without a GRCh38 WGS colony)"
note "samples  : $SAMPLES_TSV  ($N GRCh38 WGS colonies, projects: $(cut -f2 "$SAMPLES_TSV" | sort -u | paste -sd, -))"
note "arms     : ${RUN_ARMS[*]}   A=$PT_ROOT_A   B=$PT_ROOT_B"
[ "$NMISS" -eq 0 ] || note "WARNING  : $NMISS tree tip(s) have no colony row — run preflight.sh $P for the list"

# --- quick preflight (hard blockers only; the full one is preflight.sh <P>) --------
if [ "$SKIP_PRE" = 0 ]; then
    bash "$TPRT_KIT_DIR/preflight.sh" "$P" --quick || die "quick preflight failed (fix, or --skip-preflight)"
fi

# --- dry run: a bsub shim that prints instead of submitting ------------------------
if [ "$DRY" = 1 ]; then
    SHIM="$(mktemp -d "${TMPDIR:-/tmp}/tprt_dry.XXXXXX")"
    trap 'rm -rf "$SHIM"' EXIT
    echo 900000 > "$SHIM/counter"
    cat > "$SHIM/bsub" <<'SH'
#!/usr/bin/env bash
c="$(dirname "$0")/counter"; n=$(( $(cat "$c") + 1 )); echo "$n" > "$c"
printf 'DRY-RUN bsub'; for a in "$@"; do printf ' %q' "$a"; done; printf '\n'
echo "Job <$n> is submitted to queue <dry-run>."
SH
    chmod +x "$SHIM/bsub"
    export PATH="$SHIM:$PATH"
    note "DRY RUN: bsub is a printing shim; job ids below are fake (9000xx). Run dirs and samples.tsv ARE written (same content a real submit writes)."
fi

submit() {   # submit <bsub args...> -> echoes the job id
    local out
    out="$(bsub "$@" 2>&1)" || { echo "bsub FAILED: $out" >&2; exit 1; }
    echo "$out" >&2
    echo "$out" | sed -n 's/^Job <\([0-9]*\)>.*/\1/p'
}

jobid_of() { awk -F'\t' -v p="$2" '$1==p{id=$2} END{print id}' "$1"; }   # <jobids.tsv> <phase>

active_jobs() {   # ids of the last submission of an arm that are still pending/running
    local f="$1"
    [ -s "$f" ] && command -v bjobs >/dev/null 2>&1 || return 0
    cut -f2 "$f" | while read -r j; do
        bjobs -noheader -o stat "$j" 2>/dev/null | grep -qE 'PEND|RUN|PSUSP|USUSP|SSUSP' && echo "$j"
    done
    return 0
}

# shellcheck disable=SC2034  # JID_CG / JID_AN are read through the dep_all nameref
declare -A JID_SD JID_CG JID_AN
for ARM in "${RUN_ARMS[@]}"; do
    ROOT="$(arm_root "$ARM")"
    WR="$(arm_workroot "$P" "$ARM")"
    JF="$WR/jobids.tsv"
    mkdir -p "$WR" "$TPRT_ROOT/tmp/$ARM"
    if [ "$DRY" = 0 ] && [ "$FORCE" = 0 ]; then
        live="$(active_jobs "$JF" | paste -sd' ' -)"
        [ -z "$live" ] || die "arm $ARM of $P still has live jobs from its last submission ($live) — wait, bkill them, or --force"
    fi
    : > "$JF"
    if [ "$ARM" = A ]; then SDM="$A_SD_MEM"; CIM="$A_CI_MEM"; else SDM="$B_SD_MEM"; CIM="$B_CI_MEM"; fi
    WAIT=""
    if [ "$ARM" = B ] && [ -n "${JID_SD[A]:-}" ]; then WAIT="ended(${JID_SD[A]}[*])"; fi

    note "=== arm $ARM: $(arm_disc_cfg "$ARM" | xargs basename) + $(arm_geno_cfg "$ARM" | xargs basename) + src/config.py <- $(arm_py_base "$ARM")"
    env WORKROOT="$WR" \
        RESULTS_DIR="$(arm_results "$P" "$ARM")" \
        STAGING_ROOT="$STAGING_ROOT" \
        VENV="$VENV" \
        DISCOVER_BIN="$DISCOVER_BIN" GENOTYPE_BIN="$GENOTYPE_BIN" \
        DISC_CFG="$(arm_disc_cfg "$ARM")" GENO_CFG="$(arm_geno_cfg "$ARM")" \
        SAMTOOLS_MODULE="$SAMTOOLS_MODULE" \
        STAGE_THROTTLE="$STAGE_THROTTLE" GT_THROTTLE="$GT_THROTTLE" \
        SD_MEM_T1="$SDM" SD_MEM_T2="$SD_MEM_T2" GT_MEM_T1="$GT_MEM_T1" GT_MEM_T2="$GT_MEM_T2" \
        CI_MEM="$CIM" CI_CORES="$CI_CORES" AN_MEM="$AN_MEM" AN_CORES="$AN_CORES" QUEUE="$QUEUE" \
        PT_JOB_PREFIX="${P}_${ARM}" PT_NO_CLEANUP=1 PT_JOBIDS_FILE="$JF" PT_SD_WAIT="$WAIT" \
        bash "$ROOT/cluster/pipeline.sh" submit-list "$P" "$SAMPLES_TSV"
    JID_SD[$ARM]="$(jobid_of "$JF" sd)"
    # shellcheck disable=SC2034  # read through the dep_all nameref
    JID_CG[$ARM]="$(jobid_of "$JF" cg)"
    JID_AN[$ARM]="$(jobid_of "$JF" an)"
    [ -n "${JID_AN[$ARM]}" ] || die "arm $ARM: could not read its job ids from $JF"
done

dep_all() {   # dep_all <phase-array-name>: "done(a) && done(b)" over the arms submitted now
    local -n m="$1"; local d="" a
    for a in "${RUN_ARMS[@]}"; do d="${d:+$d && }done(${m[$a]})"; done
    echo "$d"
}

PASS="TPRT_ROOT='$TPRT_ROOT' PT_ROOT_A='$PT_ROOT_A' VENV='$VENV' RESULTS_BASE='$RESULTS_BASE' PATIENTS_DIR='$PATIENTS_DIR' STAGING_ROOT='$STAGING_ROOT'"
mkdir -p "$PD/logs"
if [ "$EVAL" = 1 ]; then
    jid_ev="$(submit -J "${P}_tprt_eval" -w "$(dep_all JID_AN)" \
        -o "$PD/logs/eval.%J.log" -e "$PD/logs/eval.%J.err" \
        -n "$EV_CORES" -q "$QUEUE" -M "$EV_MEM" -R "select[mem>$EV_MEM] rusage[mem=$EV_MEM] span[hosts=1]" \
        "$PASS bash '$TPRT_KIT_DIR/evaluate.sh' '$P'")"
    note "evaluation : $jid_ev  (after: $(dep_all JID_AN))"
fi
if [ "$CLEANUP" = 1 ]; then
    CROOT="$(arm_root "${RUN_ARMS[0]}")"
    jid_cu="$(submit -J "${P}_tprt_cleanup" -w "$(dep_all JID_CG)" \
        -o "$PD/logs/cleanup.%J.log" -e "$PD/logs/cleanup.%J.err" \
        -n 1 -q "$QUEUE" -M 1000 -R "select[mem>1000] rusage[mem=1000]" \
        "PT_RUNDIR='$(arm_rundir "$P" "${RUN_ARMS[0]}")' bash '$CROOT/cluster/pipeline.sh' cleanup")"
    note "cleanup    : $jid_cu  (staged BAMs of $P, after: $(dep_all JID_CG))"
else
    note "cleanup    : not submitted (staged BAMs are kept for re-runs; --cleanup, or later:"
    note "             PT_RUNDIR=$(arm_rundir "$P" "${RUN_ARMS[0]}") bash $(arm_root "${RUN_ARMS[0]}")/cluster/pipeline.sh cleanup )"
fi

note "watch      : bjobs -w | grep ${P}_ ;  for a in ${RUN_ARMS[*]}; do PT_JOB_PREFIX=${P}_\$a WORKROOT=$TPRT_ROOT/$P/\$a bash $PT_ROOT_B/cluster/pipeline.sh status $P; done"
note "report     : $(eval_dir "$P")/ab_report.md   (copied to $RESULTS_BASE/$P/)"
