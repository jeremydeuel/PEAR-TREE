#!/usr/bin/env bash
# =============================================================================
# PEAR-TREE end-to-end pipeline for ONE patient, on farm22 (LSF).
#
#   cluster/pipeline.sh submit      <PROJECT_ID> <PATIENT_ID>   # samples from irods.txt
#   cluster/pipeline.sh submit-list <PATIENT_ID> <SAMPLES_TSV>  # samples from a prebuilt list
#   cluster/pipeline.sh status      <PATIENT_ID>
#
# SAMPLES_TSV is two tab-separated columns: <sample><TAB><iRODS project id>. This
# is what cluster/fleet.sh builds from patients/<organ>/<patient>/colonies.tsv, and
# it lets a single patient span several projects (per-sample staging).
#
# Submits the whole DAG and returns immediately. One self-dispatching file: the
# LSF tasks re-invoke this same script with an internal subcommand.
#
#   phase                depends on                 what it does
#   -------------------  -------------------------  -----------------------------
#   0 resolve            (login node)               build samples.tsv (sample<TAB>proj)
#   1 sd[1-N]%K          -                          PER SAMPLE: stageBam.pl -> discover
#   1r sdr               ended(sd)                  RETRY oom/timeout discoveries at tier2
#   2 combine            done(sdr)                  combine_insertions -> contract
#   3 gt[1-N]%K          done(combine)              PER SAMPLE: index + genotype (+stats)
#   3r gtr               ended(gt)                  RETRY oom/timeout genotypes at tier2
#   4 combine_genotypes  done(gtr)                  -> <patient>.genotypes.csv.gz
#   5 cleanup            done(combine_genotypes)    delete this patient's staged BAMs
#   6 annotate           done(combine_genotypes)    annotate_v2 (runs beside cleanup)
#
# WHY stage+discover share a task: discovery for a sample starts the instant THAT
# sample lands — no waiting for the other 191, no LSF element-dependency tricks,
# and one sample's failure cannot touch another's.
#
# OOM / TIMEOUT ESCALATION (sdr, gtr controllers): an OOM (TERM_MEMLIMIT) or wall
# (TERM_RUNLIMIT) kill is a SIGKILL the task cannot trap, and -M is fixed at submit
# time. So a controller job runs AFTER each array, classifies each genuine failure
# from its own LSF log, and resubmits only the oom/timeout ones ONCE at tier2
# (2x-4x mem + the `long` queue). There is no tier3: a sample that still fails, or
# that failed for any non-oom/timeout reason, is FLAGGED and EXCLUDED
# (<phase>_excluded.tsv) and the phase proceeds without it. If EVERY sample fails,
# the controller writes FLEET_FATAL.<phase> and exits non-zero -> the patient aborts.
#
# ERROR TOLERANCE (by design — iRODS routinely lists samples with no data; for
# PD44579, 192 samples were listed and only 174 exist, and the phylogeny has
# exactly those 174 leaves):
#   * a sample with no BAM logs a marker and exits 0 — never fails the array
#   * fan-in barriers use ended(), not done(), so dead samples don't poison them
#   * downstream steps GLOB for what exists rather than assuming N files
#   * every task is idempotent -> re-running `submit` resumes where it stopped
#   * cleanup hangs off done() so a failed genotyping never deletes your BAMs
#   * no two tasks ever write the same file (Lustre multi-writer corrupts): the
#     "missing sample" record is one marker FILE per sample, not a shared list
#
# ASSEMBLY: defaults to the GRCh38 configs. combine_insertions is assembly-specific
# and is gated by install.sh check-config against a staged BAM header before it runs.
# =============================================================================
set -euo pipefail

SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
PT_ROOT="$(cd "$(dirname "$SELF")/.." && pwd)"

# --- site configuration (override via environment) ----------------------------
WORKROOT="${WORKROOT:-/lustre/scratch126/casm/teams/team273/users/jd43/pt_runs}"
STAGING_ROOT="${STAGING_ROOT:-/lustre/scratch126/casm/staging/team273/jd43}"
IRODS_TXT="${IRODS_TXT:-/lustre/scratch126/casm/teams/team273/users/jd43/pt_hu_trees/irods.txt}"
VENV="${VENV:-/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE/venv}"
RESULTS_DIR="${RESULTS_DIR:-$HOME/results}"

DISCOVER_BIN="${DISCOVER_BIN:-$PT_ROOT/rust/peartree-discovery/target/release/peartree-discovery}"
GENOTYPE_BIN="${GENOTYPE_BIN:-$PT_ROOT/rust/peartree-genotype/target/release/peartree-genotype}"
DISC_CFG="${DISC_CFG:-$PT_ROOT/cluster/config.discovery.grch38}"
GENO_CFG="${GENO_CFG:-$PT_ROOT/cluster/config.genotype.grch38}"
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-samtools-1.19}"

# --- resources (tuned from the PD44579 run) -----------------------------------
STAGE_THROTTLE="${STAGE_THROTTLE:-20}"   # concurrent stageBam.pl -> bounds iRODS + Lustre I/O
GT_THROTTLE="${GT_THROTTLE:-30}"
# tier1 budgets are MEASURED on the full PD44579 run (peak RSS from LSF, with headroom):
#   discovery          12.2 GB peak  -> 16 GB
#   combine_insertions 15.5 GB peak  -> 24 GB
#   genotype (Rust)    140 MB peak   ->  2 GB
#   combine_genotypes   1.3 GB peak  ->  8 GB
SD_MEM_T1="${SD_MEM_T1:-16000}"; SD_MEM_T2="${SD_MEM_T2:-32000}"   # discovery: tier1 / tier2
GT_MEM_T1="${GT_MEM_T1:-2000}";  GT_MEM_T2="${GT_MEM_T2:-8000}"    # genotype:  tier1 / tier2
CI_MEM="${CI_MEM:-24000}"; CI_CORES="${CI_CORES:-8}"
CG_MEM="${CG_MEM:-8000}"; CG_CORES="${CG_CORES:-8}"
AN_MEM="${AN_MEM:-32000}"; AN_CORES="${AN_CORES:-4}"
QUEUE="${QUEUE:-normal}"                  # tier1 wall: 12 h
RETRY_QUEUE="${RETRY_QUEUE:-long}"        # tier2 wall: 48 h
CTRL_QUEUE="${CTRL_QUEUE:-week}"          # controller waits out the tier2 array (up to 48 h)

# --- optional hooks for wrappers (cluster/tprt/run_ab.sh); all unset = default behaviour ---
#   PT_JOB_PREFIX   LSF job-name prefix (default: the patient id). Two runs of one patient
#                   (A/B arms) otherwise share job names, so bjobs/status cannot tell them apart.
#   PT_SD_WAIT      extra LSF dependency for the stage+discover array, e.g. "ended(123[*])"
#                   (a same-size array dependency is element-wise in LSF: element i waits for i).
#   PT_NO_CLEANUP=1 do not submit phase 5 (staged-BAM cleanup); the wrapper owns cleanup
#                   because another run still reads the same staged BAMs.
#   PT_JOBIDS_FILE  append "<phase><TAB><jobid>" for every submitted job (for dependencies).
PT_SD_WAIT="${PT_SD_WAIT:-}"
PT_NO_CLEANUP="${PT_NO_CLEANUP:-0}"
PT_JOBIDS_FILE="${PT_JOBIDS_FILE:-}"

log() { echo "[$(date +%H:%M:%S)] $*"; }

bam_path() { echo "$STAGING_ROOT/$1/$2/mapped_sample/$2.sample.dupmarked.bam"; }  # proj sample

# fields of a samples.tsv line: 1=sample 2=proj
samp_field() { sed -n "${1}p" "$SAMPLES" | cut -f"$2"; }

# Submit and echo the job id (stderr carries LSF's own message).
submit_job() {
    local out
    out="$(bsub "$@" 2>&1)" || { echo "bsub FAILED: $out" >&2; exit 1; }
    echo "$out" >&2
    echo "$out" | sed -n 's/^Job <\([0-9]*\)>.*/\1/p'
}

# Block until an LSF job (or array) has ended. Prefer bwait; fall back to polling.
wait_job() {
    local jid="$1"
    if command -v bwait >/dev/null 2>&1 && bwait -w "ended($jid)" 2>/dev/null; then return 0; fi
    while bjobs -a -noheader -o stat "$jid" 2>/dev/null | grep -qiE 'PEND|RUN|WAIT|PROV'; do sleep 30; done
}

load_env() {
    : "${PT_RUNDIR:?internal subcommands need PT_RUNDIR}"
    # shellcheck disable=SC1091
    source "$PT_RUNDIR/run.env"
    JOB_PREFIX="${JOB_PREFIX:-$PATIENT_ID}"   # run.env written before the hook existed has none
}

record_jobid() {   # record_jobid <phase> <jobid>  (PT_JOBIDS_FILE hook; no-op when unset)
    [ -n "$PT_JOBIDS_FILE" ] && printf '%s\t%s\n' "$1" "$2" >> "$PT_JOBIDS_FILE"
    return 0
}

# =============================================================================
# phase 0 — build samples.tsv + freeze run.env, then submit the DAG
# =============================================================================

# write run.env (called after RUNDIR + SAMPLES exist)
freeze_env() {
    local RUNDIR="$1"
    cat > "$RUNDIR/run.env" <<EOF
PATIENT_ID='$PATIENT_ID'
JOB_PREFIX='$JOB_PREFIX'
RUNDIR='$RUNDIR'
SAMPLES='$RUNDIR/samples.tsv'
PT_ROOT='$PT_ROOT'
STAGING_ROOT='$STAGING_ROOT'
VENV='$VENV'
RESULTS_DIR='$RESULTS_DIR'
DISCOVER_BIN='$DISCOVER_BIN'
GENOTYPE_BIN='$GENOTYPE_BIN'
DISC_CFG='$DISC_CFG'
GENO_CFG='$GENO_CFG'
SAMTOOLS_MODULE='$SAMTOOLS_MODULE'
STAGE_THROTTLE='$STAGE_THROTTLE'
SD_MEM_T2='$SD_MEM_T2'
GT_MEM_T2='$GT_MEM_T2'
RETRY_QUEUE='$RETRY_QUEUE'
CI_CORES='$CI_CORES'
CG_CORES='$CG_CORES'
EOF
}

cmd_submit() {
    local PROJECT_ID="${1:?usage: pipeline.sh submit <PROJECT_ID> <PATIENT_ID>}"
    PATIENT_ID="${2:?usage: pipeline.sh submit <PROJECT_ID> <PATIENT_ID>}"
    local RUNDIR="$WORKROOT/$PATIENT_ID"
    [ -r "$IRODS_TXT" ] || { echo "cannot read IRODS_TXT=$IRODS_TXT" >&2; exit 1; }
    preflight_binaries
    mkdir -p "$RUNDIR"/{discovery,genotypes,insertions,stats,logs,missing}

    # >PROJECT block of irods.txt, filtered to this patient -> sample<TAB>proj
    awk -v proj=">$PROJECT_ID" -v pat="$PATIENT_ID" -v p="$PROJECT_ID" '
        $0 ~ /^>/ { inblk = ($0 == proj); next }
        inblk && index($0, pat) { print $1 "\t" p }
    ' "$IRODS_TXT" | sort -u > "$RUNDIR/samples.tsv"
    finish_submit "$RUNDIR"
}

cmd_submit_list() {
    PATIENT_ID="${1:?usage: pipeline.sh submit-list <PATIENT_ID> <SAMPLES_TSV>}"
    local SRC="${2:?usage: pipeline.sh submit-list <PATIENT_ID> <SAMPLES_TSV>}"
    [ -s "$SRC" ] || { echo "empty/missing samples list: $SRC" >&2; exit 1; }
    preflight_binaries
    local RUNDIR="$WORKROOT/$PATIENT_ID"
    mkdir -p "$RUNDIR"/{discovery,genotypes,insertions,stats,logs,missing}
    # accept <sample> or <sample><TAB><proj>; a missing proj is an error (we can't stage it)
    awk -F'\t' 'NF>=2 && $1!="" && $2!="" {print $1"\t"$2}' "$SRC" | sort -u > "$RUNDIR/samples.tsv"
    local bad; bad="$(awk -F'\t' 'NF<2 || $1=="" || $2==""' "$SRC" | wc -l | tr -d ' ')"
    [ "$bad" -eq 0 ] || echo "WARNING: dropped $bad line(s) from $SRC lacking a project id" >&2
    finish_submit "$RUNDIR"
}

preflight_binaries() {
    for f in "$DISCOVER_BIN" "$GENOTYPE_BIN" "$DISC_CFG" "$GENO_CFG"; do
        [ -e "$f" ] || { echo "missing: $f (run cluster/build.sh?)" >&2; exit 1; }
    done
    [ -x "$VENV/bin/python" ] || { echo "missing venv python: $VENV/bin/python" >&2; exit 1; }
}

finish_submit() {
    local RUNDIR="$1"
    SAMPLES="$RUNDIR/samples.tsv"
    JOB_PREFIX="${PT_JOB_PREFIX:-$PATIENT_ID}"
    local N; N="$(wc -l < "$SAMPLES" | tr -d ' ')"
    [ "$N" -gt 0 ] || { echo "no samples for $PATIENT_ID" >&2; exit 1; }
    rm -f "$RUNDIR"/FLEET_FATAL.* "$RUNDIR"/*_excluded.tsv 2>/dev/null || true
    freeze_env "$RUNDIR"
    log "patient $PATIENT_ID: $N samples"
    log "run dir: $RUNDIR"
    log "NOTE: listed samples with no data are skipped, not errors."
    submit_dag "$RUNDIR" "$N"
    log "submitted. watch:  bjobs -A ;  $SELF status $PATIENT_ID"
}

submit_dag() {
    local RUNDIR="$1" N="$2"
    local R="span[hosts=1]"
    local W="PT_RUNDIR='$RUNDIR' bash '$SELF'"
    local jid_sd jid_sdr jid_ci jid_gt jid_gtr jid_cg jid_cu jid_an

    local SD_WAIT=(); [ -n "$PT_SD_WAIT" ] && SD_WAIT=(-w "$PT_SD_WAIT")
    jid_sd=$(submit_job -J "${JOB_PREFIX}_sd[1-$N]%$STAGE_THROTTLE" ${SD_WAIT[@]+"${SD_WAIT[@]}"} \
        -o "$RUNDIR/logs/sd.%I.log" -e "$RUNDIR/logs/sd.%I.err" \
        -n 1 -q "$QUEUE" -M "$SD_MEM_T1" -R "select[mem>$SD_MEM_T1] rusage[mem=$SD_MEM_T1] $R" \
        "$W stage-discover \$LSB_JOBINDEX")
    log "phase 1 stage+discover : $jid_sd  (tier1 ${SD_MEM_T1}MB/$QUEUE)${PT_SD_WAIT:+  waits: $PT_SD_WAIT}"; record_jobid sd "$jid_sd"

    jid_sdr=$(submit_job -J "${JOB_PREFIX}_sdr" -w "ended($jid_sd)" \
        -o "$RUNDIR/logs/sd_retry.%J.log" -e "$RUNDIR/logs/sd_retry.%J.err" \
        -n 1 -q "$CTRL_QUEUE" -M 1000 -R "select[mem>1000] rusage[mem=1000]" \
        "$W retry discover")
    log "phase 1r retry(disc)   : $jid_sdr  (tier2 ${SD_MEM_T2}MB/$RETRY_QUEUE, oom/timeout only)"; record_jobid sdr "$jid_sdr"

    jid_ci=$(submit_job -J "${JOB_PREFIX}_ci" -w "done($jid_sdr)" \
        -o "$RUNDIR/logs/combine.%J.log" -e "$RUNDIR/logs/combine.%J.err" \
        -n "$CI_CORES" -q "$QUEUE" -M "$CI_MEM" -R "select[mem>$CI_MEM] rusage[mem=$CI_MEM] $R" \
        "$W combine")
    log "phase 2 combine        : $jid_ci"; record_jobid ci "$jid_ci"

    jid_gt=$(submit_job -J "${JOB_PREFIX}_gt[1-$N]%$GT_THROTTLE" -w "done($jid_ci)" \
        -o "$RUNDIR/logs/gt.%I.log" -e "$RUNDIR/logs/gt.%I.err" \
        -n 1 -q "$QUEUE" -M "$GT_MEM_T1" -R "select[mem>$GT_MEM_T1] rusage[mem=$GT_MEM_T1] $R" \
        "$W genotype \$LSB_JOBINDEX")
    log "phase 3 genotype       : $jid_gt  (tier1 ${GT_MEM_T1}MB/$QUEUE)"; record_jobid gt "$jid_gt"

    jid_gtr=$(submit_job -J "${JOB_PREFIX}_gtr" -w "ended($jid_gt)" \
        -o "$RUNDIR/logs/gt_retry.%J.log" -e "$RUNDIR/logs/gt_retry.%J.err" \
        -n 1 -q "$CTRL_QUEUE" -M 1000 -R "select[mem>1000] rusage[mem=1000]" \
        "$W retry genotype")
    log "phase 3r retry(geno)   : $jid_gtr  (tier2 ${GT_MEM_T2}MB/$RETRY_QUEUE, oom/timeout only)"; record_jobid gtr "$jid_gtr"

    jid_cg=$(submit_job -J "${JOB_PREFIX}_cg" -w "done($jid_gtr)" \
        -o "$RUNDIR/logs/combine_gt.%J.log" -e "$RUNDIR/logs/combine_gt.%J.err" \
        -n "$CG_CORES" -q "$QUEUE" -M "$CG_MEM" -R "select[mem>$CG_MEM] rusage[mem=$CG_MEM] $R" \
        "$W combine-genotypes")
    log "phase 4 combine_gt     : $jid_cg"; record_jobid cg "$jid_cg"

    if [ "$PT_NO_CLEANUP" = 1 ]; then
        log "phase 5 cleanup        : NOT submitted (PT_NO_CLEANUP=1; the caller owns cleanup)"
    else
        jid_cu=$(submit_job -J "${JOB_PREFIX}_cleanup" -w "done($jid_cg)" \
            -o "$RUNDIR/logs/cleanup.%J.log" -e "$RUNDIR/logs/cleanup.%J.err" \
            -n 1 -q "$QUEUE" -M 1000 -R "select[mem>1000] rusage[mem=1000]" \
            "$W cleanup")
        log "phase 5 cleanup        : $jid_cu  (only on success)"; record_jobid cu "$jid_cu"
    fi

    jid_an=$(submit_job -J "${JOB_PREFIX}_annotate" -w "done($jid_cg)" \
        -o "$RUNDIR/logs/annotate.%J.log" -e "$RUNDIR/logs/annotate.%J.err" \
        -n "$AN_CORES" -q "$QUEUE" -M "$AN_MEM" -R "select[mem>$AN_MEM] rusage[mem=$AN_MEM] $R" \
        "$W annotate")
    log "phase 6 annotate       : $jid_an"; record_jobid an "$jid_an"
}

# =============================================================================
# phase 1 — stage THIS sample, then discover it (fault-isolated, idempotent)
# =============================================================================
cmd_stage_discover() {
    load_env
    local IDX="${1:?stage-discover <index>}"
    local SAMPLE PROJ
    SAMPLE="$(samp_field "$IDX" 1)"; PROJ="$(samp_field "$IDX" 2)"
    [ -n "$SAMPLE" ] || { echo "no sample at line $IDX"; exit 0; }
    [ -n "$PROJ" ]   || { echo "$SAMPLE: no project id in samples.tsv" >&2; exit 1; }

    local OUT="$RUNDIR/discovery/$SAMPLE.txt.gz"
    if [ -s "$OUT" ]; then log "$SAMPLE: discovery exists, skipping"; exit 0; fi

    local BAM; BAM="$(bam_path "$PROJ" "$SAMPLE")"
    if [ ! -s "$BAM" ]; then
        log "$SAMPLE: staging from iRODS project $PROJ"
        module load dataImportExport >/dev/null 2>&1 || true
        stageBam.pl --lustre 126 --types m --sample "$SAMPLE" \
            --project "$PROJ" -o "$STAGING_ROOT" -fo || \
            log "$SAMPLE: stageBam.pl returned non-zero (continuing)"
    fi

    if [ ! -s "$BAM" ]; then
        : > "$RUNDIR/missing/$SAMPLE"          # one marker file per sample: no shared-file writes
        log "$SAMPLE: no BAM after staging -> skipping (this is normal, not an error)"
        exit 0
    fi

    log "$SAMPLE: discovering"
    local TMP="$OUT.tmp.$$"
    if ! "$DISCOVER_BIN" --step discover --bam "$BAM" --out "$TMP" --threads 1 --config "$DISC_CFG"; then
        log "$SAMPLE: DISCOVERY FAILED"; rm -f "$TMP" "$TMP".*; exit 1
    fi
    mv -f "$TMP" "$OUT"
    for ext in stats.json splice.tsv hallmarks.tsv evidence.tsv.gz; do
        [ -e "$TMP.$ext" ] && mv -f "$TMP.$ext" "$OUT.$ext"
    done
    log "$SAMPLE: done -> $OUT"
}

# =============================================================================
# retry controller — escalate oom/timeout failures ONCE to tier2, then exclude
#   usage (internal): retry <discover|genotype>
# =============================================================================
cmd_retry() {
    load_env
    cd "$RUNDIR"
    local PHASE="${1:?retry <discover|genotype>}"
    local OUTDIR PREFIX WORKCMD MEM_T2 THROT need_disc
    case "$PHASE" in
        discover) OUTDIR=discovery; PREFIX=sd; WORKCMD=stage-discover; MEM_T2="$SD_MEM_T2"; THROT="$STAGE_THROTTLE"; need_disc=0 ;;
        genotype) OUTDIR=genotypes; PREFIX=gt; WORKCMD=genotype;       MEM_T2="$GT_MEM_T2"; THROT="$STAGE_THROTTLE"; need_disc=1 ;;
        *) echo "retry: bad phase '$PHASE'" >&2; exit 2 ;;
    esac
    local N; N="$(wc -l < "$SAMPLES" | tr -d ' ')"

    # classify: emit "idx<TAB>sample<TAB>reason" for every GENUINE failure.
    # reason from the sample's own LSF log: TERM_MEMLIMIT->oom, TERM_RUNLIMIT->timeout, else error.
    # prefers the tier2 log (<prefix>_r2.<i>.log) when present.
    classify() {
        local i S out logf errf reason
        for ((i=1; i<=N; i++)); do
            S="$(sed -n "${i}p" "$SAMPLES" | cut -f1)"; [ -n "$S" ] || continue
            out="$OUTDIR/$S.txt.gz"
            [ -s "$out" ] && continue                                  # succeeded
            [ -e "missing/$S" ] && continue                            # legit no-data
            if [ "$need_disc" = 1 ] && [ ! -s "discovery/$S.txt.gz" ]; then continue; fi  # never eligible to genotype
            if [ -s "logs/${PREFIX}_r2.${i}.log" ] || [ -s "logs/${PREFIX}_r2.${i}.err" ]; then
                logf="logs/${PREFIX}_r2.${i}.log"; errf="logs/${PREFIX}_r2.${i}.err"
            else
                logf="logs/${PREFIX}.${i}.log"; errf="logs/${PREFIX}.${i}.err"
            fi
            reason=error
            if   grep -qs TERM_MEMLIMIT "$logf" "$errf" 2>/dev/null; then reason=oom
            elif grep -qs TERM_RUNLIMIT "$logf" "$errf" 2>/dev/null; then reason=timeout
            fi
            printf '%s\t%s\t%s\n' "$i" "$S" "$reason"
        done
    }

    local fails; fails="$(classify)"
    if [ -n "$fails" ]; then
        local retry_idx; retry_idx="$(awk -F'\t' '$3=="oom"||$3=="timeout"{print $1}' <<<"$fails" | paste -sd, -)"
        if [ -n "$retry_idx" ]; then
            log "$PHASE: tier2 escalation (${MEM_T2}MB / $RETRY_QUEUE) for indices [$retry_idx]"
            local R="span[hosts=1]" jid
            jid=$(submit_job -J "${JOB_PREFIX}_${PREFIX}_r2[$retry_idx]%$THROT" \
                -o "$RUNDIR/logs/${PREFIX}_r2.%I.log" -e "$RUNDIR/logs/${PREFIX}_r2.%I.err" \
                -n 1 -q "$RETRY_QUEUE" -M "$MEM_T2" -R "select[mem>$MEM_T2] rusage[mem=$MEM_T2] $R" \
                "PT_RUNDIR='$RUNDIR' bash '$SELF' $WORKCMD \$LSB_JOBINDEX")
            log "$PHASE: waiting on tier2 array $jid"
            wait_job "$jid"
            fails="$(classify)"     # re-evaluate against tier2 logs/outputs
        else
            log "$PHASE: $(wc -l <<<"$fails") failure(s), none oom/timeout — not retrying"
        fi
    else
        log "$PHASE: no failures"
    fi

    # ---- after tier2: flag+exclude survivors; fatal only if NOTHING succeeded ----
    local exf="$RUNDIR/${PHASE}_excluded.tsv"
    if [ -n "$fails" ]; then
        { printf 'sample\treason\n'; awk -F'\t' '{print $2"\t"$3}' <<<"$fails"; } > "$exf"
        log "$PHASE: EXCLUDING $(($(wc -l <<<"$fails"))) sample(s):"
        sed 's/^/    /' "$exf"
    else
        rm -f "$exf"
    fi

    # find, not ls: under `set -o pipefail` an empty glob made `ls` fail and killed this
    # controller with exit 2 BEFORE the FATAL message below (PD37449 pilot, 2026-10-04).
    local succ; succ="$(find "$OUTDIR" -maxdepth 1 -name '*.txt.gz' 2>/dev/null | wc -l | tr -d ' ')"
    if [ "$succ" -eq 0 ]; then
        local nmiss; nmiss="$(find missing -maxdepth 1 -type f 2>/dev/null | wc -l | tr -d ' ')"
        { echo "FATAL: $PHASE produced ZERO outputs for $PATIENT_ID — all samples failed."
          [ "$nmiss" -gt 0 ] && echo "  $nmiss sample(s) have a missing/ marker = no BAM after staging; check logs/${PREFIX}.<i>.log for the stageBam.pl error"
          [ -n "$fails" ] && awk -F'\t' '{print $2"\t"$3}' <<<"$fails"; } | tee "$RUNDIR/FLEET_FATAL.$PHASE" >&2
        exit 1
    fi
    log "$PHASE: $succ output(s) present; proceeding"
}

# =============================================================================
# TPRT one-sided loci: the genotype contract
# =============================================================================
# combine_insertions leaves one-sided loci (`contig:L-oneside_L` / `contig:oneside_R-R`) out
# of <patient>.genotyping.txt.gz. When the genotype config enables `one_sided_loci`
# (cluster/config.genotype.grch38.tprt), the combine task appends them once, to
# <patient>.genotyping.tprt.txt.gz, and every genotype task reads that file instead. Built
# ONLY in the single combine task (no two tasks ever write the same file).
geno_one_sided() { grep -Eq '^[[:space:]]*one_sided_loci[[:space:]]*=[[:space:]]*(true|True|1)' "$GENO_CFG" 2>/dev/null; }

geno_contract() {   # the contract the genotype tasks read (path relative to nothing: absolute)
    if geno_one_sided; then
        echo "$RUNDIR/insertions/$PATIENT_ID.genotyping.tprt.txt.gz"
    else
        echo "$RUNDIR/insertions/$PATIENT_ID.genotyping.txt.gz"
    fi
}

ensure_geno_contract() {   # (combine task only) build the one-sided extension if enabled
    geno_one_sided || return 0
    local base="$RUNDIR/insertions/$PATIENT_ID.genotyping.txt.gz"
    local ext="$RUNDIR/insertions/$PATIENT_ID.genotyping.tprt.txt.gz"
    local comb="$RUNDIR/insertions/$PATIENT_ID.combined.txt.gz"
    if [ -s "$ext" ] && [ "$ext" -nt "$base" ]; then log "one-sided contract exists"; return 0; fi
    "$VENV/bin/python" "$PT_ROOT/src/genotyping_contract_oneside.py" \
        --contract "$base" --combined "$comb" --out "$ext" \
        || { echo "one-sided contract extension FAILED" >&2; exit 1; }
    log "contract (+one-sided): $ext ($(zcat "$ext" | grep -c '^>') loci)"
}

# =============================================================================
# phase 2 — combine_insertions over whatever discovery files exist
# =============================================================================
cmd_combine() {
    load_env
    cd "$RUNDIR"
    local CONTRACT="insertions/$PATIENT_ID.genotyping.txt.gz"
    if [ -s "$CONTRACT" ]; then log "contract exists, skipping"; ensure_geno_contract; exit 0; fi

    shopt -s nullglob
    local files=(discovery/*.txt.gz)
    shopt -u nullglob
    local n_missing; n_missing="$(ls "$RUNDIR/missing" 2>/dev/null | wc -l | tr -d ' ')"
    log "combining ${#files[@]} discovery files ($n_missing samples had no data)"
    [ "${#files[@]}" -gt 0 ] || { echo "no discovery files at all — aborting" >&2; exit 1; }

    # ASSEMBLY GATE — combine_insertions is assembly-specific; prove the config matches
    # a real staged BAM before spending cycles (a mismatch is silently wrong coordinates).
    local a_bam="" a_sample a_idx
    for a_sample in $(printf '%s\n' "${files[@]}" | sed -e 's/\.txt\.gz$//' -e 's#.*/##'); do
        a_idx="$(grep -nxF -m1 -- "$a_sample" <(cut -f1 "$SAMPLES") | cut -d: -f1)"
        [ -n "$a_idx" ] || continue
        a_bam="$(bam_path "$(samp_field "$a_idx" 2)" "$a_sample")"
        [ -s "$a_bam" ] && break || a_bam=""
    done
    if [ -n "$a_bam" ]; then
        log "verifying config against the data ($a_bam)"
        bash "$PT_ROOT/cluster/install.sh" check-config --bam "$a_bam" --pt-root "$PT_ROOT" \
            || { echo "assembly/config check FAILED — refusing to combine" >&2; exit 1; }
    else
        log "WARNING: no staged BAM left to verify the assembly against — proceeding unchecked"
    fi

    "$VENV/bin/python" "$PT_ROOT/src/main.py" --step combine_insertions \
        --discovery_files "${files[@]}" --out "insertions/$PATIENT_ID" --threads "$CI_CORES"

    [ -s "$CONTRACT" ] || { echo "combine_insertions produced no $CONTRACT" >&2; exit 1; }
    log "contract: $CONTRACT ($(zcat "$CONTRACT" | grep -c '^>') loci)"
    ensure_geno_contract
}

# =============================================================================
# phase 3 — genotype THIS sample against the contract (+ per-BAM stats sidecar)
# =============================================================================
cmd_genotype() {
    load_env
    local IDX="${1:?genotype <index>}"
    local SAMPLE PROJ
    SAMPLE="$(samp_field "$IDX" 1)"; PROJ="$(samp_field "$IDX" 2)"
    [ -n "$SAMPLE" ] || exit 0

    local OUT="$RUNDIR/genotypes/$SAMPLE.txt.gz"
    if [ -s "$OUT" ]; then log "$SAMPLE: genotype exists, skipping"; exit 0; fi

    local BAM; BAM="$(bam_path "$PROJ" "$SAMPLE")"
    if [ ! -s "$BAM" ]; then log "$SAMPLE: no BAM -> skipping (expected)"; exit 0; fi

    if [ ! -s "$BAM.bai" ] && [ ! -s "${BAM%.bam}.bai" ]; then
        log "$SAMPLE: indexing BAM"
        module load "$SAMTOOLS_MODULE" >/dev/null 2>&1 || true
        samtools index "$BAM" || { log "$SAMPLE: samtools index failed"; exit 1; }
    fi

    write_bam_stats "$SAMPLE" "$BAM"   # step 5: avg coverage / #reads / read length

    local CONTRACT; CONTRACT="$(geno_contract)"
    [ -s "$CONTRACT" ] || { log "$SAMPLE: contract $CONTRACT missing (re-run combine)"; exit 1; }
    local TMP="$OUT.tmp.$$"
    if ! "$GENOTYPE_BIN" --step genotype --bam "$BAM" \
            --insertions "$CONTRACT" \
            --out "$TMP" --threads 1 --config "$GENO_CFG"; then
        log "$SAMPLE: GENOTYPING FAILED"; rm -f "$TMP"; exit 1
    fi
    mv -f "$TMP" "$OUT"
    log "$SAMPLE: genotyped -> $OUT"
}

# per-BAM minimal stats: sample n_reads read_len mean_cov  (idxstats-based, no extra full pass)
write_bam_stats() {
    local S="$1" BAM="$2" sf="$RUNDIR/stats/$1.tsv"
    [ -s "$sf" ] && return 0
    module load "$SAMTOOLS_MODULE" >/dev/null 2>&1 || true
    local idx rl
    idx="$(samtools idxstats "$BAM" 2>/dev/null)" || { log "$S: idxstats failed (no stats)"; return 0; }
    rl="$(samtools view "$BAM" 2>/dev/null | head -1 | awk '{print length($10)}')"; [ -n "$rl" ] || rl=NA
    printf '%s\n' "$idx" | awk -v s="$S" -v rl="$rl" '
        $1 ~ /^(chr)?([0-9]+|X|Y)$/ { m+=$3; L+=$2 }
        END { cov=(L>0 && rl!="NA") ? sprintf("%.2f", m*rl/L) : "NA";
              printf "%s\t%d\t%s\t%s\n", s, m, rl, cov }' > "$sf.tmp.$$" && mv -f "$sf.tmp.$$" "$sf"
}

# =============================================================================
# phase 4 — combine_genotypes
# =============================================================================
cmd_combine_genotypes() {
    load_env
    cd "$RUNDIR"
    local CALLS="$PATIENT_ID.genotypes.csv.gz"
    if [ -s "$CALLS" ]; then log "calls exist, skipping"; exit 0; fi

    shopt -s nullglob
    local files=(genotypes/*.txt.gz)
    shopt -u nullglob
    log "combining ${#files[@]} genotype files"
    [ "${#files[@]}" -gt 0 ] || { echo "no genotype files — aborting" >&2; exit 1; }

    "$VENV/bin/python" "$PT_ROOT/src/main.py" --step combine_genotypes \
        --genotypes "${files[@]}" --out "$CALLS" --threads "$CG_CORES"

    [ -s "$CALLS" ] || { echo "combine_genotypes produced no $CALLS" >&2; exit 1; }

    # per-BAM stats summary (step 5/7): one table for the whole patient
    local STATS_SUM="$PATIENT_ID.bam_stats.tsv"
    { printf 'sample\tn_reads\tread_len\tmean_cov\n'; cat stats/*.tsv 2>/dev/null; } > "$STATS_SUM"

    mkdir -p "$RESULTS_DIR/$PATIENT_ID"
    cp -f "$CALLS" "$STATS_SUM" "insertions/$PATIENT_ID".* "$RESULTS_DIR/$PATIENT_ID/" 2>/dev/null || true
    for x in discover genotype; do
        [ -s "${x}_excluded.tsv" ] && cp -f "${x}_excluded.tsv" "$RESULTS_DIR/$PATIENT_ID/" || true
    done
    log "calls: $CALLS  (copied to $RESULTS_DIR/$PATIENT_ID)"
}

# =============================================================================
# phase 5 — delete THIS patient's staged BAMs (only reached on success)
# =============================================================================
cmd_cleanup() {
    load_env
    local freed=0 SAMPLE PROJ
    # per-sample deletion only; the staging project dir can hold OTHER patients.
    while IFS=$'\t' read -r SAMPLE PROJ; do
        [ -n "$SAMPLE" ] && [ -n "$PROJ" ] || continue
        local d="$STAGING_ROOT/$PROJ/$SAMPLE"
        if [ -d "$d" ]; then rm -rf "$d" && freed=$((freed+1)); fi
    done < "$SAMPLES"
    # rmdir each distinct project dir, only if now empty
    cut -f2 "$SAMPLES" | sort -u | while read -r PROJ; do
        [ -n "$PROJ" ] && rmdir "$STAGING_ROOT/$PROJ" 2>/dev/null || true
    done
    log "cleanup: removed $freed staged sample dirs under $STAGING_ROOT"
}

# =============================================================================
# phase 6 — annotate (element classification); needs no BAMs
# =============================================================================
cmd_annotate() {
    load_env
    cd "$RUNDIR"
    log "annotating $PATIENT_ID"
    "$VENV/bin/python" "$PT_ROOT/tools/annotate_v2.py" "$PATIENT_ID" "$PATIENT_ID.annotated.csv.gz"
    if [ -s "$PATIENT_ID.annotated.csv.gz" ]; then
        mkdir -p "$RESULTS_DIR/$PATIENT_ID"
        cp -f "$PATIENT_ID.annotated.csv.gz" "$RESULTS_DIR/$PATIENT_ID/"
        log "annotated -> $RESULTS_DIR/$PATIENT_ID/$PATIENT_ID.annotated.csv.gz"
    else
        log "annotate produced no output"; exit 1
    fi
}

# =============================================================================
# status
# =============================================================================
cmd_status() {
    local PATIENT_ID="${1:?usage: pipeline.sh status <PATIENT_ID>}"
    local RUNDIR="$WORKROOT/$PATIENT_ID"
    [ -d "$RUNDIR" ] || { echo "no run dir $RUNDIR" >&2; exit 1; }
    local n_s n_d n_m n_g
    n_s=$(wc -l < "$RUNDIR/samples.tsv" | tr -d ' ')
    n_d=$(ls "$RUNDIR"/discovery/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')
    n_m=$(ls "$RUNDIR"/missing 2>/dev/null | wc -l | tr -d ' ')
    n_g=$(ls "$RUNDIR"/genotypes/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')
    echo "patient   : $PATIENT_ID   ($RUNDIR)"
    echo "samples   : $n_s"
    echo "no data   : $n_m  (iRODS lists samples that were never sequenced)"
    echo "discovered: $n_d / $((n_s - n_m))"
    echo "genotyped : $n_g / $((n_s - n_m))"
    for x in discover genotype; do
        [ -s "$RUNDIR/${x}_excluded.tsv" ] && echo "excluded  : $x -> $(($(wc -l < "$RUNDIR/${x}_excluded.tsv")-1)) sample(s), see ${x}_excluded.tsv"
        [ -s "$RUNDIR/FLEET_FATAL.$x" ] && echo "FATAL     : $x — all samples failed (see FLEET_FATAL.$x)"
    done
    [ -s "$RUNDIR/insertions/$PATIENT_ID.genotyping.txt.gz" ] && echo "contract  : yes" || echo "contract  : no"
    [ -s "$RUNDIR/insertions/$PATIENT_ID.genotyping.tprt.txt.gz" ] && echo "contract+1: yes (one-sided loci appended)"
    [ -s "$RUNDIR/$PATIENT_ID.genotypes.csv.gz" ] && echo "calls     : yes" || echo "calls     : no"
    bjobs -J "${PT_JOB_PREFIX:-$PATIENT_ID}_*" -A 2>/dev/null || true
}

# --- dispatch -----------------------------------------------------------------
case "${1:-}" in
    submit)            shift; cmd_submit "$@" ;;
    submit-list)       shift; cmd_submit_list "$@" ;;
    status)            shift; cmd_status "$@" ;;
    stage-discover)    shift; cmd_stage_discover "$@" ;;
    retry)             shift; cmd_retry "$@" ;;
    combine)           shift; cmd_combine "$@" ;;
    genotype)          shift; cmd_genotype "$@" ;;
    combine-genotypes) shift; cmd_combine_genotypes "$@" ;;
    cleanup)           shift; cmd_cleanup "$@" ;;
    annotate)          shift; cmd_annotate "$@" ;;
    *) sed -n '2,35p' "$SELF"; exit 1 ;;
esac
