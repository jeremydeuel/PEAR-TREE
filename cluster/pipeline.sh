#!/usr/bin/env bash
# =============================================================================
# PEAR-TREE end-to-end pipeline for ONE patient, on farm22 (LSF).
#
#   cluster/pipeline.sh submit <PROJECT_ID> <PATIENT_ID>
#   cluster/pipeline.sh status <PATIENT_ID>
#
# Submits the whole DAG and returns immediately. One self-dispatching file: the
# LSF tasks re-invoke this same script with an internal subcommand.
#
#   phase                depends on                 what it does
#   -------------------  -------------------------  -----------------------------
#   0 resolve            (login node)               irods.txt -> sample list
#   1 stagedisc[1-N]%K   -                          PER SAMPLE: stageBam.pl -> discover
#   2 combine            ended(stagedisc)           combine_insertions -> contract
#   3 genotype[1-N]%K    done(combine)              PER SAMPLE: index + genotype
#   4 combine_genotypes  ended(genotype)            -> <patient>.genotypes.csv.gz
#   5 cleanup            done(combine_genotypes)    delete this patient's staged BAMs
#   6 annotate           done(combine_genotypes)    annotate_v2 (runs beside cleanup)
#
# WHY stage+discover share a task: discovery for a sample starts the instant THAT
# sample lands — no waiting for the other 191, no LSF element-dependency tricks,
# and one sample's failure cannot touch another's.
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
# ASSEMBLY: runs discovery directly on the GRCh37/hs37d5 BAMs (no bwa remap to
# hs1). combine_insertions remaps clips to hs1 and reconciles them to GRCh37 via
# the hs1->hg19 chain, so final coordinates are GRCh37.
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
DISC_CFG="${DISC_CFG:-$PT_ROOT/cluster/config.discovery.grch37}"
GENO_CFG="${GENO_CFG:-$PT_ROOT/cluster/config.genotype.grch37}"
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-samtools-1.19}"

# --- resources (tuned from the PD44579 run) -----------------------------------
STAGE_THROTTLE="${STAGE_THROTTLE:-20}"   # concurrent stageBam.pl -> bounds iRODS + Lustre I/O
GT_THROTTLE="${GT_THROTTLE:-30}"
SD_MEM="${SD_MEM:-16000}"                # discovery peaked at 12.2 GB on PD44579
CI_MEM="${CI_MEM:-64000}"; CI_CORES="${CI_CORES:-8}"
GT_MEM="${GT_MEM:-8000}"
CG_MEM="${CG_MEM:-32000}"; CG_CORES="${CG_CORES:-8}"
AN_MEM="${AN_MEM:-32000}"; AN_CORES="${AN_CORES:-4}"
QUEUE="${QUEUE:-normal}"                 # 12 h wall

log() { echo "[$(date +%H:%M:%S)] $*"; }

bam_path() { echo "$STAGING_ROOT/$1/$2/mapped_sample/$2.sample.dupmarked.bam"; }  # project sample

# Submit and echo the job id (stderr carries LSF's own message).
submit_job() {
    local out
    out="$(bsub "$@" 2>&1)" || { echo "bsub FAILED: $out" >&2; exit 1; }
    echo "$out" >&2
    echo "$out" | sed -n 's/^Job <\([0-9]*\)>.*/\1/p'
}

load_env() {
    : "${PT_RUNDIR:?internal subcommands need PT_RUNDIR}"
    # shellcheck disable=SC1091
    source "$PT_RUNDIR/run.env"
}

# =============================================================================
# phase 0 — resolve samples + submit the DAG
# =============================================================================
cmd_submit() {
    local PROJECT_ID="${1:?usage: pipeline.sh submit <PROJECT_ID> <PATIENT_ID>}"
    local PATIENT_ID="${2:?usage: pipeline.sh submit <PROJECT_ID> <PATIENT_ID>}"
    local RUNDIR="$WORKROOT/$PATIENT_ID"

    [ -r "$IRODS_TXT" ] || { echo "cannot read IRODS_TXT=$IRODS_TXT" >&2; exit 1; }
    for f in "$DISCOVER_BIN" "$GENOTYPE_BIN" "$DISC_CFG" "$GENO_CFG"; do
        [ -e "$f" ] || { echo "missing: $f (run cluster/build.sh?)" >&2; exit 1; }
    done
    [ -x "$VENV/bin/python" ] || { echo "missing venv python: $VENV/bin/python" >&2; exit 1; }

    mkdir -p "$RUNDIR"/{discovery,genotypes,insertions,logs,missing}

    # sample list: the >PROJECT block of irods.txt, filtered to this patient.
    local SAMPLES="$RUNDIR/samples.txt"
    awk -v proj=">$PROJECT_ID" -v pat="$PATIENT_ID" '
        $0 ~ /^>/ { inblk = ($0 == proj); next }
        inblk && index($0, pat) { print $1 }
    ' "$IRODS_TXT" | sort -u > "$SAMPLES"
    local N; N="$(wc -l < "$SAMPLES" | tr -d ' ')"
    [ "$N" -gt 0 ] || { echo "no samples for $PATIENT_ID in project $PROJECT_ID of $IRODS_TXT" >&2; exit 1; }

    # freeze the run's parameters so every task sees identical settings
    cat > "$RUNDIR/run.env" <<EOF
PROJECT_ID='$PROJECT_ID'
PATIENT_ID='$PATIENT_ID'
RUNDIR='$RUNDIR'
SAMPLES='$SAMPLES'
PT_ROOT='$PT_ROOT'
STAGING_ROOT='$STAGING_ROOT'
VENV='$VENV'
RESULTS_DIR='$RESULTS_DIR'
DISCOVER_BIN='$DISCOVER_BIN'
GENOTYPE_BIN='$GENOTYPE_BIN'
DISC_CFG='$DISC_CFG'
GENO_CFG='$GENO_CFG'
SAMTOOLS_MODULE='$SAMTOOLS_MODULE'
CI_CORES='$CI_CORES'
CG_CORES='$CG_CORES'
EOF

    log "patient $PATIENT_ID / project $PROJECT_ID: $N samples listed in iRODS"
    log "run dir: $RUNDIR"
    log "NOTE: iRODS lists samples that have no data; those are skipped, not errors."

    local R="span[hosts=1]"
    local jid_sd jid_ci jid_gt jid_cg jid_cu jid_an

    jid_sd=$(submit_job -J "${PATIENT_ID}_sd[1-$N]%$STAGE_THROTTLE" \
        -o "$RUNDIR/logs/sd.%I.log" -e "$RUNDIR/logs/sd.%I.err" \
        -n 1 -q "$QUEUE" -M "$SD_MEM" -R "select[mem>$SD_MEM] rusage[mem=$SD_MEM] $R" \
        "PT_RUNDIR='$RUNDIR' bash '$SELF' stage-discover \$LSB_JOBINDEX")
    log "phase 1 stage+discover : $jid_sd"

    jid_ci=$(submit_job -J "${PATIENT_ID}_ci" -w "ended($jid_sd)" \
        -o "$RUNDIR/logs/combine.%J.log" -e "$RUNDIR/logs/combine.%J.err" \
        -n "$CI_CORES" -q "$QUEUE" -M "$CI_MEM" -R "select[mem>$CI_MEM] rusage[mem=$CI_MEM] $R" \
        "PT_RUNDIR='$RUNDIR' bash '$SELF' combine")
    log "phase 2 combine        : $jid_ci"

    jid_gt=$(submit_job -J "${PATIENT_ID}_gt[1-$N]%$GT_THROTTLE" -w "done($jid_ci)" \
        -o "$RUNDIR/logs/gt.%I.log" -e "$RUNDIR/logs/gt.%I.err" \
        -n 1 -q "$QUEUE" -M "$GT_MEM" -R "select[mem>$GT_MEM] rusage[mem=$GT_MEM] $R" \
        "PT_RUNDIR='$RUNDIR' bash '$SELF' genotype \$LSB_JOBINDEX")
    log "phase 3 genotype       : $jid_gt"

    jid_cg=$(submit_job -J "${PATIENT_ID}_cg" -w "ended($jid_gt)" \
        -o "$RUNDIR/logs/combine_gt.%J.log" -e "$RUNDIR/logs/combine_gt.%J.err" \
        -n "$CG_CORES" -q "$QUEUE" -M "$CG_MEM" -R "select[mem>$CG_MEM] rusage[mem=$CG_MEM] $R" \
        "PT_RUNDIR='$RUNDIR' bash '$SELF' combine-genotypes")
    log "phase 4 combine_gt     : $jid_cg"

    jid_cu=$(submit_job -J "${PATIENT_ID}_cleanup" -w "done($jid_cg)" \
        -o "$RUNDIR/logs/cleanup.%J.log" -e "$RUNDIR/logs/cleanup.%J.err" \
        -n 1 -q "$QUEUE" -M 1000 -R "select[mem>1000] rusage[mem=1000]" \
        "PT_RUNDIR='$RUNDIR' bash '$SELF' cleanup")
    log "phase 5 cleanup        : $jid_cu  (only on success)"

    jid_an=$(submit_job -J "${PATIENT_ID}_annotate" -w "done($jid_cg)" \
        -o "$RUNDIR/logs/annotate.%J.log" -e "$RUNDIR/logs/annotate.%J.err" \
        -n "$AN_CORES" -q "$QUEUE" -M "$AN_MEM" -R "select[mem>$AN_MEM] rusage[mem=$AN_MEM] $R" \
        "PT_RUNDIR='$RUNDIR' bash '$SELF' annotate")
    log "phase 6 annotate       : $jid_an"

    log "submitted. watch:  bjobs -A ;  $SELF status $PATIENT_ID"
}

# =============================================================================
# phase 1 — stage THIS sample, then discover it (fault-isolated, idempotent)
# =============================================================================
cmd_stage_discover() {
    load_env
    local IDX="${1:?stage-discover <index>}"
    local SAMPLE; SAMPLE="$(sed -n "${IDX}p" "$SAMPLES")"
    [ -n "$SAMPLE" ] || { echo "no sample at line $IDX"; exit 0; }

    local OUT="$RUNDIR/discovery/$SAMPLE.txt.gz"
    if [ -s "$OUT" ]; then log "$SAMPLE: discovery exists, skipping"; exit 0; fi

    local BAM; BAM="$(bam_path "$PROJECT_ID" "$SAMPLE")"
    if [ ! -s "$BAM" ]; then
        log "$SAMPLE: staging from iRODS project $PROJECT_ID"
        module load dataImportExport >/dev/null 2>&1 || true
        # never let a staging failure kill the task: many listed samples have no data
        stageBam.pl --lustre 126 --types m --sample "$SAMPLE" \
            --project "$PROJECT_ID" -o "$STAGING_ROOT" -fo || \
            log "$SAMPLE: stageBam.pl returned non-zero (continuing)"
    fi

    if [ ! -s "$BAM" ]; then
        # EXPECTED for a fair number of samples — record and succeed.
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
    # carry the sidecars ({out}.stats.json, {out}.splice.tsv for Feature B) to the final name
    for ext in stats.json splice.tsv hallmarks.tsv; do
        [ -e "$TMP.$ext" ] && mv -f "$TMP.$ext" "$OUT.$ext"
    done
    log "$SAMPLE: done -> $OUT"
}

# =============================================================================
# phase 2 — combine_insertions over whatever discovery files exist
# =============================================================================
cmd_combine() {
    load_env
    cd "$RUNDIR"
    local CONTRACT="insertions/$PATIENT_ID.genotyping.txt.gz"
    if [ -s "$CONTRACT" ]; then log "contract exists, skipping"; exit 0; fi

    shopt -s nullglob
    local files=(discovery/*.txt.gz)
    shopt -u nullglob
    local n_missing; n_missing="$(ls "$RUNDIR/missing" | wc -l | tr -d ' ')"
    log "combining ${#files[@]} discovery files ($n_missing samples had no data)"
    [ "${#files[@]}" -gt 0 ] || { echo "no discovery files at all — aborting" >&2; exit 1; }

    # ASSEMBLY GATE. Unlike discovery, combine_insertions is assembly-specific: it
    # reads flanks from genome_2bit and lifts hs1 clip hits back with a chain. A
    # mismatched config produces plausible, silently wrong coordinates — no error.
    # So prove the config matches the actual data before spending the cycles.
    local a_sample a_bam
    for a_sample in $(sed -e 's/\.txt\.gz$//' -e 's#.*/##' <<<"$(printf '%s\n' "${files[@]}")"); do
        a_bam="$(bam_path "$PROJECT_ID" "$a_sample")"
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
}

# =============================================================================
# phase 3 — genotype THIS sample against the contract
# =============================================================================
cmd_genotype() {
    load_env
    local IDX="${1:?genotype <index>}"
    local SAMPLE; SAMPLE="$(sed -n "${IDX}p" "$SAMPLES")"
    [ -n "$SAMPLE" ] || exit 0

    local OUT="$RUNDIR/genotypes/$SAMPLE.txt.gz"
    if [ -s "$OUT" ]; then log "$SAMPLE: genotype exists, skipping"; exit 0; fi

    local BAM; BAM="$(bam_path "$PROJECT_ID" "$SAMPLE")"
    if [ ! -s "$BAM" ]; then log "$SAMPLE: no BAM -> skipping (expected)"; exit 0; fi

    # genotyping needs a coordinate index (discovery did not)
    if [ ! -s "$BAM.bai" ] && [ ! -s "${BAM%.bam}.bai" ]; then
        log "$SAMPLE: indexing BAM"
        module load "$SAMTOOLS_MODULE" >/dev/null 2>&1 || true
        samtools index "$BAM" || { log "$SAMPLE: samtools index failed"; exit 1; }
    fi

    local TMP="$OUT.tmp.$$"
    if ! "$GENOTYPE_BIN" --step genotype --bam "$BAM" \
            --insertions "$RUNDIR/insertions/$PATIENT_ID.genotyping.txt.gz" \
            --out "$TMP" --threads 1 --config "$GENO_CFG"; then
        log "$SAMPLE: GENOTYPING FAILED"; rm -f "$TMP"; exit 1
    fi
    mv -f "$TMP" "$OUT"
    log "$SAMPLE: genotyped -> $OUT"
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
    mkdir -p "$RESULTS_DIR/$PATIENT_ID"
    cp -f "$CALLS" "insertions/$PATIENT_ID".* "$RESULTS_DIR/$PATIENT_ID/" 2>/dev/null || true
    log "calls: $CALLS  (copied to $RESULTS_DIR/$PATIENT_ID)"
}

# =============================================================================
# phase 5 — delete THIS patient's staged BAMs (only reached on success)
# =============================================================================
cmd_cleanup() {
    load_env
    local freed=0
    # Per-sample deletion only. The staging project dir can hold OTHER patients —
    # the legacy remapping.bsub.sh did `rm -rdf .../$PROJECT_ID/`, which also
    # deleted the inputs of still-running sibling tasks.
    while read -r SAMPLE; do
        local d="$STAGING_ROOT/$PROJECT_ID/$SAMPLE"
        if [ -d "$d" ]; then rm -rf "$d" && freed=$((freed+1)); fi
    done < "$SAMPLES"
    log "cleanup: removed $freed staged sample dirs under $STAGING_ROOT/$PROJECT_ID"
    rmdir "$STAGING_ROOT/$PROJECT_ID" 2>/dev/null || true   # only if now empty
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
    n_s=$(wc -l < "$RUNDIR/samples.txt" | tr -d ' ')
    n_d=$(ls "$RUNDIR"/discovery/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')
    n_m=$(ls "$RUNDIR"/missing 2>/dev/null | wc -l | tr -d ' ')
    n_g=$(ls "$RUNDIR"/genotypes/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')
    echo "patient   : $PATIENT_ID   ($RUNDIR)"
    echo "samples   : $n_s listed in iRODS"
    echo "no data   : $n_m  (expected: iRODS lists samples that were never sequenced)"
    echo "discovered: $n_d / $((n_s - n_m))"
    echo "genotyped : $n_g / $((n_s - n_m))"
    [ -s "$RUNDIR/insertions/$PATIENT_ID.genotyping.txt.gz" ] && echo "contract  : yes" || echo "contract  : no"
    [ -s "$RUNDIR/$PATIENT_ID.genotypes.csv.gz" ] && echo "calls     : yes" || echo "calls     : no"
    bjobs -J "${PATIENT_ID}_*" -A 2>/dev/null || true
}

# --- dispatch -----------------------------------------------------------------
case "${1:-}" in
    submit)            shift; cmd_submit "$@" ;;
    status)            shift; cmd_status "$@" ;;
    stage-discover)    shift; cmd_stage_discover "$@" ;;
    combine)           shift; cmd_combine "$@" ;;
    genotype)          shift; cmd_genotype "$@" ;;
    combine-genotypes) shift; cmd_combine_genotypes "$@" ;;
    cleanup)           shift; cmd_cleanup "$@" ;;
    annotate)          shift; cmd_annotate "$@" ;;
    *) sed -n '2,30p' "$SELF"; exit 1 ;;
esac
