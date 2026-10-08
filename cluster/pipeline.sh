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
# Genotyping implementation (plans/genotype_v2/SPEC.md): v2 = peartree-genotype2, the
# realignment genotyper with numeric per-colony output + the phylogenetic joint step in phase 4
# (needs the patient's SNV tree, see patient_tree()); legacy = peartree-genotype +
# combine_genotypes.py (string calls). v2 is the default.
GENOTYPE_IMPL="${GENOTYPE_IMPL:-v2}"
GENOTYPE2_BIN="${GENOTYPE2_BIN:-$PT_ROOT/rust/peartree-genotype2/target/release/peartree-genotype2}"
GENO2_CFG="${GENO2_CFG:-$PT_ROOT/cluster/config.genotype2.grch38}"
# BAM assembly: GRCh38 (default) or GRCh37 (hs37d5). Exported into every job: the generated
# src/config.py (cluster/tprt/arm_config.py) reads it to pick genome_2bit + the hs1 chain.
PT_ASSEMBLY="${PT_ASSEMBLY:-GRCh38}"
case "$PT_ASSEMBLY" in GRCh38) _asm_2bit=hg38 ;; GRCh37) _asm_2bit=hg19 ;; *) echo "PT_ASSEMBLY must be GRCh38 or GRCh37, got '$PT_ASSEMBLY'" >&2; exit 1 ;; esac
# reference of the BAMs' assembly for the haplotype flanks (same file as config.py's genome_2bit;
# hg19.2bit serves hs37d5 loci: genotype2 toggles the chr prefix)
GENOME_2BIT="${GENOME_2BIT:-/lustre/scratch126/casm/teams/team273/users/jd43/$_asm_2bit.2bit}"
# the patient's SNV tree (Newick, tips = colony ids); default: the one file under patients/*/<id>/
PATIENT_TREE="${PATIENT_TREE:-}"
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-samtools-1.19}"
# combine_insertions implementation: python (default, src/main.py) or rust (the
# rust/peartree-combine port, built by cluster/build.sh; byte-identical outputs, reads the same
# src/config.py through $VENV/bin/python; refuses require_independent_fragments=True).
COMBINE_IMPL="${COMBINE_IMPL:-python}"
COMBINE_BIN="${COMBINE_BIN:-$PT_ROOT/rust/peartree-combine/target/release/peartree-combine}"
# Genotyping contract: one-sided loci on/off. Unset = GENO_CFG's `one_sided_loci` (legacy
# behaviour); 0 = two-sided loci only (PD37590 decision: only insertions with both ends are
# genotyped), 1 = with the one-sided extension.
GENO_ONE_SIDED="${GENO_ONE_SIDED:-}"
# extra arguments of the genotype2 joint step (phase 4), e.g. "--ref-bias auto" with a
# GENO2_CFG that writes the pl_het_b profile (cluster/config.genotype2.grch38.refbias)
JOINT_ARGS="${JOINT_ARGS:-}"
# phase 7 (report): tools/phylo/tree_fit.py + cluster/somatic_table.py after annotate; 0 = skip
PT_REPORT="${PT_REPORT:-1}"
# per-BAM header gate after staging (empty = off): GRCh38 = @SQ chr1 must be 248956422 bp,
# GRCh37 = @SQ 1 (hs37d5 naming) must be 249250621 bp; @RG DS must say WGS and the reads must be
# >= PT_MIN_READLEN (100) bp, else the sample is logged and treated as no-data (missing/ marker).
# For sample lists not built from BAM headers (hsc_run.sh populate's iRODS fallback, colonies.tsv
# rows with readlen NA).
PT_HEADER_GATE="${PT_HEADER_GATE:-}"
PT_MIN_READLEN="${PT_MIN_READLEN:-100}"

# --- resources (tuned from the PD44579 run) -----------------------------------
STAGE_THROTTLE="${STAGE_THROTTLE:-20}"   # concurrent stageBam.pl -> bounds iRODS + Lustre I/O
PT_RESTAGE="${PT_RESTAGE:-0}"            # 1 = drop each BAM after discovery, re-stage it for genotyping (drop_bam)
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
# the patient's SNV tree: $PATIENT_TREE if set, else the single *.tree under patients/*/<id>/
patient_tree() {
    if [ -n "$PATIENT_TREE" ]; then echo "$PATIENT_TREE"; return; fi
    local t; t=$(ls "$PT_ROOT"/patients/*/"$PATIENT_ID"/*.tree 2>/dev/null | head -1)
    echo "$t"
}

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
    export PT_ASSEMBLY="${PT_ASSEMBLY:-GRCh38}"   # read by src/config.py in the python steps
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
GENOTYPE_IMPL='$GENOTYPE_IMPL'
GENOTYPE2_BIN='$GENOTYPE2_BIN'
GENO2_CFG='$GENO2_CFG'
GENOME_2BIT='$GENOME_2BIT'
PATIENT_TREE='$PATIENT_TREE'
DISC_CFG='$DISC_CFG'
GENO_CFG='$GENO_CFG'
SAMTOOLS_MODULE='$SAMTOOLS_MODULE'
COMBINE_IMPL='$COMBINE_IMPL'
COMBINE_BIN='$COMBINE_BIN'
GENO_ONE_SIDED='$GENO_ONE_SIDED'
JOINT_ARGS='$JOINT_ARGS'
PT_HEADER_GATE='$PT_HEADER_GATE'
PT_MIN_READLEN='$PT_MIN_READLEN'
PT_ASSEMBLY='$PT_ASSEMBLY'
STAGE_THROTTLE='$STAGE_THROTTLE'
PT_RESTAGE='$PT_RESTAGE'
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
    local gt_bin="$GENOTYPE_BIN" gt_cfg="$GENO_CFG"
    if [ "$GENOTYPE_IMPL" = v2 ]; then
        gt_bin="$GENOTYPE2_BIN"; gt_cfg="$GENO2_CFG"
        [ -s "$GENOME_2BIT" ] || { echo "GENOME_2BIT $GENOME_2BIT missing (needed by peartree-genotype2)" >&2; exit 1; }
        [ -n "$(patient_tree)" ] || { echo "no SNV tree for $PATIENT_ID under $PT_ROOT/patients/*/$PATIENT_ID/ (set PATIENT_TREE)" >&2; exit 1; }
    fi
    for f in "$DISCOVER_BIN" "$gt_bin" "$DISC_CFG" "$gt_cfg"; do
        [ -e "$f" ] || { echo "missing: $f (run cluster/build.sh?)" >&2; exit 1; }
    done
    [ -x "$VENV/bin/python" ] || { echo "missing venv python: $VENV/bin/python" >&2; exit 1; }
    case "$COMBINE_IMPL" in
        python) ;;
        rust) [ -x "$COMBINE_BIN" ] || { echo "COMBINE_IMPL=rust but missing: $COMBINE_BIN (run cluster/build.sh)" >&2; exit 1; } ;;
        *) echo "COMBINE_IMPL must be python or rust, got '$COMBINE_IMPL'" >&2; exit 1 ;;
    esac
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
        -Q 99 "PT_GT_WATCHDOG=1 $W genotype \$LSB_JOBINDEX")
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

    if [ "$GENOTYPE_IMPL" = v2 ] && [ "$PT_REPORT" = 1 ]; then
        jid_rp=$(submit_job -J "${JOB_PREFIX}_report" -w "done($jid_an)" \
            -o "$RUNDIR/logs/report.%J.log" -e "$RUNDIR/logs/report.%J.err" \
            -n 1 -q "$QUEUE" -M 16000 -R "select[mem>16000] rusage[mem=16000]" \
            "$W report")
        log "phase 7 report         : $jid_rp"; record_jobid rp "$jid_rp"
    fi
}

# =============================================================================
# phase 1 — stage THIS sample, then discover it (fault-isolated, idempotent)
# =============================================================================
# Stage ONE sample's BAM from iRODS into bam_path (callers check the file with -s afterwards).
# Returns 0 (stageBam.pl ran), 2 = iRODS lists no files for this sample (legit no-data),
# 1 = stageBam.pl failure. Uses SAMPLE / PROJ / BAM from the caller.
stage_bam() {
    log "$SAMPLE: staging from iRODS project $PROJ -> $(dirname "$BAM")"
    # Three stageBam.pl traps, all hit on the PD37449 pilot (2026-10-04):
    #  1. ASYNC: it submits its own LSF transfer job ("Job <N> is submitted to queue <normal>.")
    #     and returns at once -> wait on that job, else every colony looks "missing".
    #  2. CONCURRENCY: 10 tasks calling it at the same second with the same -o crossed their
    #     requests (lo0006's task listed + transferred lo0016) -> 20 duplicate transfers,
    #     18 failed, 7 colonies never requested. So: a PRIVATE -o per sample, submissions
    #     serialised under flock, and the printed sample name is verified.
    #  3. PERL: the submitting shell's modules leak into the job (samtools-1.19 -> perl 5.38);
    #     dataImportExport's perl 5.36 then dies on "ListUtil.c: loadable library and perl
    #     binaries are mismatched". So: a clean subshell (module purge, PERL* unset).
    # Layout: stageBam.pl appends <proj>/<sample>/mapped_sample/ to -o; the file list it
    # PRINTS shows a flat <-o>/mapped_sample/ path, which is not where the transfer writes.
    local PRIV="$STAGING_ROOT/.stagebam/$SAMPLE"
    rm -rf "$PRIV"; mkdir -p "$PRIV"
    local so rc=0
    so="$( {
        command -v flock >/dev/null 2>&1 && exec 9>"$STAGING_ROOT/.stagebam.lock" && flock -w 1800 9
        unset PERL5LIB PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
        module purge >/dev/null 2>&1 || true
        module load dataImportExport >/dev/null 2>&1 || true
        stageBam.pl --lustre 126 --types m --sample "$SAMPLE" --project "$PROJ" -o "$PRIV" -fo
    } 2>&1 )" || rc=$?
    printf '%s\n' "$so"
    if grep -qE 'total files 0\b' <<<"$so"; then
        log "$SAMPLE: iRODS lists no files for project $PROJ -> skipping (not an error)"
        rm -rf "$PRIV"; return 2         # iRODS has nothing for this sample: legit no-data
    fi
    # the listed sample must be OURS (trap 2); stageBam prints it on its own line
    if ! grep -qx "$SAMPLE" <<<"$so"; then
        log "$SAMPLE: stageBam.pl did not list this sample (listed: $(grep -xE 'PD[0-9]+[a-z]+_?[a-z0-9]*' <<<"$so" | paste -sd, -)) -- refusing"
        return 1
    fi
    local tj; tj="$(sed -n 's/^Job <\([0-9]*\)> is submitted.*/\1/p' <<<"$so" | tail -1)"
    if [ -n "$tj" ]; then
        log "$SAMPLE: waiting for stageBam.pl transfer job $tj"
        wait_job "$tj"
    elif [ "$rc" -ne 0 ]; then
        log "$SAMPLE: stageBam.pl failed (rc=$rc) and submitted no transfer job"; return 1
    fi
    # find the published BAM in the private dir (never the tmpExportData/ progress copy),
    # allowing lustre a short grace period, then move the set into bam_path's directory
    local t=0 got=""
    while [ -z "$got" ] && [ "$t" -le "${STAGE_GRACE_S:-600}" ]; do
        for c in "$PRIV/$PROJ/$SAMPLE/mapped_sample/$SAMPLE.sample.dupmarked.bam" \
                 "$PRIV/mapped_sample/$SAMPLE.sample.dupmarked.bam"; do
            [ -s "$c" ] && { got="$c"; break; }
        done
        [ -n "$got" ] || { sleep 30; t=$((t+30)); }
    done
    if [ -n "$got" ]; then
        mkdir -p "$(dirname "$BAM")"
        for f in "$got" "$got.bai" "$got.bas" "$got.met.gz"; do
            [ -e "$f" ] && mv -f "$f" "$(dirname "$BAM")/"
        done
        rm -rf "$PRIV"
    fi
}

# PT_RESTAGE=1: a staged BAM lives only as long as ONE job needs it -- discovery deletes it when
# done and genotyping stages it again (and deletes it after), so the peak is ~(running sd +
# running gt) BAMs instead of the whole patient (PD49229: 722 colonies ~17 TB vs 8 TB free).
# Costs a second iRODS transfer per colony. Off = BAMs stay until phase 5 cleanup.
drop_bam() {   # drop_bam <proj> <sample>  (no-op unless PT_RESTAGE=1)
    [ "${PT_RESTAGE:-0}" = 1 ] || return 0
    rm -rf "${STAGING_ROOT:?}/$1/$2" && log "$2: PT_RESTAGE -> staged BAM removed"
    rmdir "$STAGING_ROOT/$1" 2>/dev/null || true
}

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
        local src=0; stage_bam || src=$?
        if [ "$src" -eq 2 ]; then : > "$RUNDIR/missing/$SAMPLE"; exit 0; fi
        [ "$src" -eq 0 ] || exit 1
    fi

    if [ ! -s "$BAM" ]; then
        log "$SAMPLE: STAGING FAILED -- no BAM at $BAM after stageBam.pl (see output above)"
        exit 1
    fi
    module load "$SAMTOOLS_MODULE" >/dev/null 2>&1 || true
    if command -v samtools >/dev/null 2>&1 && ! samtools quickcheck "$BAM"; then
        log "$SAMPLE: staged BAM fails samtools quickcheck (truncated?): $BAM -- delete it and rerun"
        exit 1
    fi
    if [ -n "${PT_HEADER_GATE:-}" ]; then
        # the assembly's chr1 under the name DISC_CFG's allowlist uses, WGS in @RG DS, and reads
        # long enough for clip discovery (max over the first 20k records; 75 bp releases exist)
        local hdr chr1 ds rl want_name want_len
        case "$PT_HEADER_GATE" in
            GRCh38) want_name=chr1; want_len=248956422 ;;
            GRCh37) want_name=1;    want_len=249250621 ;;     # hs37d5 (numeric contigs)
            *) log "PT_HEADER_GATE must be GRCh38 or GRCh37, got '$PT_HEADER_GATE'"; exit 1 ;;
        esac
        hdr="$(samtools view -H "$BAM")" || { log "$SAMPLE: cannot read the BAM header"; exit 1; }
        chr1="$(awk -F'\t' -v want="$want_name" '$1=="@SQ" { n=""; l=""; for (i=2;i<=NF;i++) { if ($i ~ /^SN:/) n=substr($i,4); if ($i ~ /^LN:/) l=substr($i,4) } if (n==want) { print l; exit } }' <<<"$hdr")"
        ds="$(grep '^@RG' <<<"$hdr" | tr '\t' '\n' | sed -n 's/^DS://p' | sort -u | paste -sd, -)"
        rl="$( { samtools view "$BAM" 2>/dev/null || true; } | head -n 20000 | awk '{ if (length($10) > m) m = length($10) } END { print m + 0 }')"
        if [ "$chr1" != "$want_len" ] || ! grep -q WGS <<<"$ds" || [ "$rl" -lt "${PT_MIN_READLEN:-100}" ]; then
            log "$SAMPLE: HEADER GATE -- $want_name length '${chr1:-none}' ($PT_HEADER_GATE = $want_len), @RG DS '${ds:-none}' (need WGS), read length $rl (need >= ${PT_MIN_READLEN:-100}) -> skipped as no-data"
            printf 'header_gate\t%s=%s\tDS=%s\treadlen=%s\n' "$want_name" "${chr1:-none}" "${ds:-none}" "$rl" > "$RUNDIR/missing/$SAMPLE"
            drop_bam "$PROJ" "$SAMPLE"
            exit 0
        fi
        log "$SAMPLE: header gate ok ($PT_HEADER_GATE, DS ${ds}, read length $rl)"
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
    drop_bam "$PROJ" "$SAMPLE"
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
geno_one_sided() {
    case "${GENO_ONE_SIDED:-}" in 0) return 1 ;; 1) return 0 ;; esac
    grep -Eq '^[[:space:]]*one_sided_loci[[:space:]]*=[[:space:]]*(true|True|1)' "$GENO_CFG" 2>/dev/null
}

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

# genotype2 extra pass (`gt_extra_reads`, rust/peartree-genotype2/src/extra.rs): every genotype
# task needs to know which loci its colony DISCOVERED (= has reads in combine's reads FASTA) so it
# collects extra reads only at the others. tools/genotype_extra_reads.py condenses the reads
# FASTA (PD51635: 269 MB) once into <P>.members.tsv.gz. Built by the combine task; a genotype
# task of an older run builds it under a lock (the others wait for it) when it is missing.
geno_members() { echo "$RUNDIR/insertions/$PATIENT_ID.members.tsv.gz"; }

ensure_geno_members() {
    [ "$GENOTYPE_IMPL" = v2 ] || return 0
    local fa="$RUNDIR/insertions/$PATIENT_ID.insertions.reads.fa.gz" mem lock t=0
    mem="$(geno_members)"; lock="$mem.lock"
    [ -s "$fa" ] || return 0                         # no combine reads: no extra pass
    # current = newer than the reads FASTA and carries the colony count (older tables lack it)
    if [ -s "$mem" ] && [ "$mem" -nt "$fa" ] \
            && { gzip -dc "$mem" 2>/dev/null || true; } | head -n 1 | grep -q '^#colonies'; then
        return 0
    fi
    # germline-skip denominator = ALL colonies of the run: the discovery files combine read (=
    # the colonies genotyped: a sample without a BAM has a missing/ marker and no discovery
    # file), not only those with reads at some locus
    local n_col
    n_col="$(find "$RUNDIR/discovery" -maxdepth 1 -name '*.txt.gz' 2>/dev/null | wc -l | tr -d ' ')"
    if mkdir "$lock" 2>/dev/null; then
        "$VENV/bin/python" "$PT_ROOT/tools/genotype_extra_reads.py" members --reads-fa "$fa" --out "$mem" \
            --n-colonies "$n_col" || log "WARNING: members table failed -> genotyping without the extra pass"
        rmdir "$lock"
    else
        while [ -d "$lock" ] && [ "$t" -lt 900 ]; do sleep 10; t=$((t + 10)); done
    fi
    return 0
}

# =============================================================================
# phase 2 — combine_insertions over whatever discovery files exist
# =============================================================================
cmd_combine() {
    load_env
    cd "$RUNDIR"
    local CONTRACT="insertions/$PATIENT_ID.genotyping.txt.gz"
    if [ -s "$CONTRACT" ]; then log "contract exists, skipping"; ensure_geno_contract; ensure_geno_members; exit 0; fi

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

    if [ "${COMBINE_IMPL:-python}" = rust ]; then
        # same config file python reads (src/config.py), dumped to JSON by $VENV's python
        log "combine_insertions: rust ($COMBINE_BIN)"
        PEARTREE_PYTHON="$VENV/bin/python" "$COMBINE_BIN" --step combine_insertions \
            --config "$PT_ROOT/src/config.py" \
            --discovery_files "${files[@]}" --out "insertions/$PATIENT_ID" --threads "$CI_CORES"
    else
        # -u: unbuffered, so the job log shows combine's progress while it runs (it can take hours)
        "$VENV/bin/python" -u "$PT_ROOT/src/main.py" --step combine_insertions \
            --discovery_files "${files[@]}" --out "insertions/$PATIENT_ID" --threads "$CI_CORES"
    fi

    [ -s "$CONTRACT" ] || { echo "combine_insertions produced no $CONTRACT" >&2; exit 1; }
    log "contract: $CONTRACT ($(zcat "$CONTRACT" | grep -c '^>') loci)"
    ensure_geno_contract
    ensure_geno_members
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
    if [ ! -s "$BAM" ] && [ "${PT_RESTAGE:-0}" = 1 ] && [ ! -e "$RUNDIR/missing/$SAMPLE" ] \
            && [ -s "$RUNDIR/discovery/$SAMPLE.txt.gz" ]; then
        log "$SAMPLE: PT_RESTAGE -> staging again for genotyping"
        local src=0; stage_bam || src=$?
        if [ "$src" -eq 2 ]; then : > "$RUNDIR/missing/$SAMPLE"; exit 0; fi
        { [ "$src" -eq 0 ] && [ -s "$BAM" ]; } || { log "$SAMPLE: RE-STAGING FAILED -- no BAM at $BAM"; exit 1; }
    fi
    if [ ! -s "$BAM" ]; then log "$SAMPLE: no BAM -> skipping (expected)"; exit 0; fi

    if [ ! -s "$BAM.bai" ] && [ ! -s "${BAM%.bam}.bai" ]; then
        log "$SAMPLE: indexing BAM"
        module load "$SAMTOOLS_MODULE" >/dev/null 2>&1 || true
        samtools index "$BAM" || { log "$SAMPLE: samtools index failed"; exit 1; }
    fi

    write_bam_stats "$SAMPLE" "$BAM"   # step 5: avg coverage / #reads / read length

    local CONTRACT; CONTRACT="$(geno_contract)"
    [ -s "$CONTRACT" ] || { log "$SAMPLE: contract $CONTRACT missing (re-run combine)"; exit 1; }
    local TMP="$OUT.tmp.$$" rc=0
    if [ "$GENOTYPE_IMPL" = v2 ]; then
        # realignment genotyper: full junction consensus from combine + reference flanks
        local COMBINED="$RUNDIR/insertions/$PATIENT_ID.combined.txt.gz"
        [ -s "$COMBINED" ] || { log "$SAMPLE: $COMBINED missing (re-run combine)"; exit 1; }
        # extra pass (a no-op unless GENO2_CFG sets gt_extra_reads): this colony's discovery
        # memberships; it writes $TMP.extra_reads.fa.gz next to the output
        local xargs=()
        ensure_geno_members
        [ -s "$(geno_members)" ] && xargs=(--members "$(geno_members)" --sample "$SAMPLE")
        run_genotyper "$GENOTYPE2_BIN" --step genotype --bam "$BAM" \
            --insertions "$CONTRACT" --combined "$COMBINED" --reference "$GENOME_2BIT" \
            --out "$TMP" --threads 1 --config "$GENO2_CFG" ${xargs[@]+"${xargs[@]}"} || rc=$?
    else
        run_genotyper "$GENOTYPE_BIN" --step genotype --bam "$BAM" \
            --insertions "$CONTRACT" \
            --out "$TMP" --threads 1 --config "$GENO_CFG" || rc=$?
    fi
    if [ "$rc" -eq 99 ]; then rm -f "$TMP" "$TMP.extra_reads.fa.gz"; exit 99; fi    # stalled: LSF requeues (bsub -Q 99)
    if [ "$rc" -ne 0 ]; then log "$SAMPLE: GENOTYPING FAILED"; rm -f "$TMP" "$TMP.extra_reads.fa.gz"; exit 1; fi
    # sidecar first: the output's presence is what marks the sample done
    [ -e "$TMP.extra_reads.fa.gz" ] && mv -f "$TMP.extra_reads.fa.gz" "$OUT.extra_reads.fa.gz"
    mv -f "$TMP" "$OUT"
    log "$SAMPLE: genotyped -> $OUT"
    drop_bam "$PROJ" "$SAMPLE"
}

# Run a genotyper. In the tier-1 array (PT_GT_WATCHDOG=1, submitted with bsub -Q 99) a run that
# prints nothing to stderr for GT_STALL_S seconds is killed and returns 99, so LSF requeues the
# element: same job id (downstream dependencies intact), a fresh slot. The genotypers print a
# progress line every 1000 loci, which on a normal node comes every 20-80 s (PD51635: 115
# colonies 311-1611 s for 21,729 loci); two elements on node-14-08 (233 jobs, load 107) slowed
# ~70x to 1000 loci per 1600 s and would have hit the 12 h wall. At most GT_MAX_REQUEUE stall
# requeues per sample, then the run goes to the end wherever it is. Elsewhere (tier-2 retry):
# a plain run.
run_genotyper() {
    if [ "${PT_GT_WATCHDOG:-0}" != 1 ]; then "$@"; return; fi
    local cnt="$RUNDIR/genotypes/.$SAMPLE.requeues" n=0 max="${GT_MAX_REQUEUE:-3}" stall="${GT_STALL_S:-900}"
    [ -s "$cnt" ] && n="$(cat "$cnt")"
    if [ "$n" -ge "$max" ]; then
        log "$SAMPLE: $n stall requeues already -> running without the watchdog"; "$@"; return
    fi
    local hb="$RUNDIR/logs/gt.$SAMPLE.progress"
    : > "$hb"
    "$@" 2> >(tee -a "$hb" >&2) &
    local pid=$! age
    while kill -0 "$pid" 2>/dev/null; do
        sleep "${GT_WATCH_EVERY_S:-30}"
        age=$(( $(date +%s) - $(stat -c %Y "$hb") ))
        if [ "$age" -gt "$stall" ] && kill -0 "$pid" 2>/dev/null; then
            kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null || true
            echo $((n + 1)) > "$cnt"
            log "$SAMPLE: no genotyper progress for ${age}s on $(hostname) -> exit 99, LSF requeue $((n + 1))/$max"
            return 99
        fi
    done
    wait "$pid"
}

# per-BAM minimal stats: sample n_reads read_len mean_cov  (idxstats-based, no extra full pass)
write_bam_stats() {
    local S="$1" BAM="$2" sf="$RUNDIR/stats/$1.tsv"
    [ -s "$sf" ] && return 0
    module load "$SAMTOOLS_MODULE" >/dev/null 2>&1 || true
    local idx rl
    idx="$(samtools idxstats "$BAM" 2>/dev/null)" || { log "$S: idxstats failed (no stats)"; return 0; }
    # read length from the first record. NOT `samtools view | head -1`: head closing the pipe
    # kills samtools with SIGPIPE, and under `set -euo pipefail` that killed the whole genotype
    # task (exit 141, all 10 PD37449 arm-A tasks, 2026-10-04). Mask samtools' status instead.
    rl="$( { samtools view "$BAM" 2>/dev/null || true; } | awk 'NR==1 {print length($10); exit}')"; [ -n "$rl" ] || rl=NA
    printf '%s\n' "$idx" | awk -v s="$S" -v rl="$rl" '
        $1 ~ /^(chr)?([0-9]+|X|Y)$/ { m+=$3; L+=$2 }
        END { cov=(L>0 && rl!="NA") ? sprintf("%.2f", m*rl/L) : "NA";
              printf "%s\t%d\t%s\t%s\n", s, m, rl, cov }' > "$sf.tmp.$$" && mv -f "$sf.tmp.$$" "$sf"
}

# =============================================================================
# phase 4 — combine_genotypes
# =============================================================================
# genotype2 extra pass: the per-colony sidecars, kept only where the joint step calls the colony a
# carrier (P >= genotype2_io.P_CARRIER), -> insertions/<P>.insertions.genotype_reads.fa.gz, which
# annotate_v2 (tools/rte) reads as classification evidence. Never junction evidence: the report's
# hard rules read insertions.reads.fa.gz only.
merge_genotype_reads() {
    local CALLS="$1" out="insertions/$PATIENT_ID.insertions.genotype_reads.fa.gz"
    [ "$GENOTYPE_IMPL" = v2 ] && [ -s "$CALLS" ] || return 0
    compgen -G "genotypes/*.extra_reads.fa.gz" >/dev/null || return 0
    if [ -s "$out" ] && [ "$out" -nt "$CALLS" ] && { [ ! -e "$(geno_members)" ] || [ "$out" -nt "$(geno_members)" ]; }; then return 0; fi
    "$VENV/bin/python" "$PT_ROOT/tools/genotype_extra_reads.py" merge --genotype-dir genotypes \
        --matrix "$CALLS" --out "$out" --members "$(geno_members)" \
        || log "WARNING: genotype reads merge failed (annotate runs without them)"
}

cmd_combine_genotypes() {
    load_env
    cd "$RUNDIR"
    local CALLS="$PATIENT_ID.genotypes.csv.gz"
    if [ -s "$CALLS" ]; then log "calls exist, skipping"; merge_genotype_reads "$CALLS"; exit 0; fi

    shopt -s nullglob
    local files=(genotypes/*.txt.gz)
    shopt -u nullglob
    log "combining ${#files[@]} genotype files"
    [ "${#files[@]}" -gt 0 ] || { echo "no genotype files — aborting" >&2; exit 1; }

    if [ "$GENOTYPE_IMPL" = v2 ]; then
        # phylogenetic joint genotyping: every locus is placed on the SNV tree (ROOT / branch /
        # NOISE / INDEP); $CALLS becomes the NUMERIC P(carrier) matrix (rows loci, cols
        # colonies; annotate_v2 reads it), $PATIENT_ID.joint.tsv the per-locus table.
        local TREE; TREE="$(patient_tree)"
        [ -s "$TREE" ] || { echo "no SNV tree for $PATIENT_ID (set PATIENT_TREE)" >&2; exit 1; }
        # shellcheck disable=SC2086  # JOINT_ARGS is a word list
        "$GENOTYPE2_BIN" --step joint --tree "$TREE" --genotypes "${files[@]}" \
            --out "$PATIENT_ID.joint.tsv" --matrix "$CALLS" ${JOINT_ARGS:-}
        [ -s "$PATIENT_ID.joint.tsv" ] || { echo "joint step produced no $PATIENT_ID.joint.tsv" >&2; exit 1; }
    else
        "$VENV/bin/python" "$PT_ROOT/src/main.py" --step combine_genotypes \
            --genotypes "${files[@]}" --out "$CALLS" --threads "$CG_CORES"
    fi

    [ -s "$CALLS" ] || { echo "phase 4 produced no $CALLS" >&2; exit 1; }
    merge_genotype_reads "$CALLS"     # before the copy below: it goes to RESULTS_DIR with insertions/

    # per-BAM stats summary (step 5/7): one table for the whole patient
    local STATS_SUM="$PATIENT_ID.bam_stats.tsv"
    { printf 'sample\tn_reads\tread_len\tmean_cov\n'; cat stats/*.tsv 2>/dev/null; } > "$STATS_SUM"

    mkdir -p "$RESULTS_DIR/$PATIENT_ID"
    cp -f "$CALLS" "$STATS_SUM" "insertions/$PATIENT_ID".* "$RESULTS_DIR/$PATIENT_ID/" 2>/dev/null || true
    [ -s "$PATIENT_ID.joint.tsv" ] && cp -f "$PATIENT_ID.joint.tsv" "$RESULTS_DIR/$PATIENT_ID/" || true
    [ -s "$PATIENT_ID.joint.refbias.tsv" ] && cp -f "$PATIENT_ID.joint.refbias.tsv" "$RESULTS_DIR/$PATIENT_ID/" || true
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
    # annotate_v2 imports `src.config` and `tools.rte` -> the repo root must be importable; run
    # from $RUNDIR, sys.path[0] is tools/ only (ModuleNotFoundError: src, PD37449 arm A 2026-10-04)
    # perl modules leaked from the submitting shell break dfamscan.pl (ListUtil.c mismatch)
    unset PERL5LIB PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
    PYTHONPATH="$PT_ROOT${PYTHONPATH:+:$PYTHONPATH}" \
        "$VENV/bin/python" -u "$PT_ROOT/tools/annotate_v2.py" "$PATIENT_ID" "$PATIENT_ID.annotated.csv.gz"
    if [ -s "$PATIENT_ID.annotated.csv.gz" ]; then
        mkdir -p "$RESULTS_DIR/$PATIENT_ID"
        cp -f "$PATIENT_ID.annotated.csv.gz" "$RESULTS_DIR/$PATIENT_ID/"
        log "annotated -> $RESULTS_DIR/$PATIENT_ID/$PATIENT_ID.annotated.csv.gz"
    else
        log "annotate produced no output"; exit 1
    fi
}

# =============================================================================
# phase 7 — report (genotype2 only): tree_fit cross-check + the annotated somatic table
# =============================================================================
cmd_report() {
    load_env
    cd "$RUNDIR"
    local TREE; TREE="$(patient_tree)"
    [ -s "$TREE" ] || { echo "no SNV tree for $PATIENT_ID (set PATIENT_TREE)" >&2; exit 1; }
    local CALLS="$PATIENT_ID.genotypes.csv.gz" FIT="$RUNDIR/fit"
    [ -s "$CALLS" ] || { echo "no $CALLS (phase 4 not done)" >&2; exit 1; }
    log "tree_fit -> $FIT"
    rm -rf "$FIT"; mkdir -p "$FIT"
    (cd "$PT_ROOT" && "$VENV/bin/python" tools/phylo/tree_fit.py --genotypes "$RUNDIR/$CALLS" \
        --genotype-dir "$RUNDIR/genotypes" --tree "$TREE" --out "$FIT") > "$RUNDIR/logs/tree_fit.log" 2>&1 \
        || { log "tree_fit FAILED (logs/tree_fit.log); the table is built without it"; tail -5 "$RUNDIR/logs/tree_fit.log" >&2; }
    local known; known="$(dirname "$TREE")/known_insertions.tsv"
    local args=(--patient "$PATIENT_ID" --joint "$PATIENT_ID.joint.tsv" --genotype-dir genotypes
                --insertions-dir insertions --genome "$GENOME_2BIT" --out "$PATIENT_ID.somatic.xlsx")
    [ -s "$FIT/phylo_fit.tsv" ] && args+=(--fit-v2 "$FIT/phylo_fit.tsv")
    [ -s "$PATIENT_ID.annotated.csv.gz" ] && args+=(--annotation "$PATIENT_ID.annotated.csv.gz")
    [ -s "$PATIENT_ID.joint.refbias.tsv" ] && args+=(--refbias "$PATIENT_ID.joint.refbias.tsv")
    [ -s "$known" ] && args+=(--known "$known")
    "$VENV/bin/python" "$PT_ROOT/cluster/somatic_table.py" "${args[@]}"
    mkdir -p "$RESULTS_DIR/$PATIENT_ID"
    cp -f "$PATIENT_ID.somatic.xlsx" "$RESULTS_DIR/$PATIENT_ID/"
    [ -s "$FIT/summary.md" ] && cp -f "$FIT/summary.md" "$RESULTS_DIR/$PATIENT_ID/$PATIENT_ID.tree_fit_summary.md" || true
    log "report -> $RESULTS_DIR/$PATIENT_ID/$PATIENT_ID.somatic.xlsx"
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
    report)            shift; cmd_report "$@" ;;
    *) sed -n '2,35p' "$SELF"; exit 1 ;;
esac
