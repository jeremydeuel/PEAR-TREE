#!/usr/bin/env bash
# =============================================================================
# PEAR-TREE fleet driver — run every GRCh38 patient in patients/ end-to-end.
#
#   bash cluster/fleet.sh plan            # dry run: worklist, skip-list, capacity — no jobs
#   bash cluster/fleet.sh run [PAT...]    # steps 1-9 for all GRCh38 patients (or the named ones)
#   bash cluster/fleet.sh status          # fleet-wide progress table
#
# Wraps cluster/pipeline.sh (the per-patient stage->discover->combine->genotype->
# combine_genotypes DAG). fleet.sh adds the pre-flight (clean scratch, capacity),
# builds each patient's sample list from patients/<organ>/<patient>/colonies.tsv,
# throttles how many patients stage at once, and rolls up the QC gate.
#
# SCOPE: only colonies.tsv rows with assembly==GRCh38 & ds~WGS. Patients with none
# are written to the skip-list for the future non-GRCh38 bwa-remap extension.
#
# STORAGE (locked with Jeremy 2026-07-18):
#   scratch (disposable) = STAGING_ROOT           -> wiped in step 1, per-patient in step 9
#   persistent lustre    = WORKROOT/<patient>/     -> discovery+genotype intermediates KEPT
#   backed-up finals     = RESULTS_DIR/<patient>/  -> small CSVs + stats + QC (NFS home)
# =============================================================================
set -euo pipefail

SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
PT_ROOT="$(cd "$(dirname "$SELF")/.." && pwd)"
PIPELINE="$PT_ROOT/cluster/pipeline.sh"

WORKROOT="${WORKROOT:-/lustre/scratch126/casm/teams/team273/users/jd43/pt_runs}"
STAGING_ROOT="${STAGING_ROOT:-/lustre/scratch126/casm/staging/team273/jd43}"
RESULTS_DIR="${RESULTS_DIR:-$HOME/results}"
PATIENTS_DIR="${PATIENTS_DIR:-$PT_ROOT/patients}"
FLEET_DIR="$WORKROOT/_fleet"

# fleet-level knobs
PATIENTS_INFLIGHT="${PATIENTS_INFLIGHT:-4}"   # max patients staging concurrently
STAGE_THROTTLE="${STAGE_THROTTLE:-20}"        # per-patient concurrent stageBam.pl (passed through)
AVG_BAM_GB="${AVG_BAM_GB:-45}"                # for the capacity estimate
LUSTRE_VOL="${LUSTRE_VOL:-/lustre/scratch126}"
LUSTRE_GROUP="${LUSTRE_GROUP:-team273}"       # binding quota is the team, not the user
export WORKROOT STAGING_ROOT RESULTS_DIR STAGE_THROTTLE

log() { echo "[$(date +%H:%M:%S)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

# ---- worklist: scan colonies.tsv -> per-patient GRCh38-WGS sample lists -------
# Classifies each patient by how much of its WGS tree is GRCh38:
#   pure  (all WGS colonies GRCh38)  -> patients.txt          run now
#   mixed (some GRCh38, some not)    -> mixed_partial.tsv     remap-queue (partial tree if forced)
#   none  (no GRCh38 WGS colony)     -> skip_non_grch38.tsv   remap-queue
# INCLUDE_MIXED=1 promotes the mixed patients' GRCh38 subset into patients.txt anyway.
build_worklist() {
    mkdir -p "$FLEET_DIR/worklist"
    : > "$FLEET_DIR/patients.txt"
    : > "$FLEET_DIR/skip_non_grch38.tsv"
    { printf 'patient\tgrch38_wgs\ttotal_wgs\tgrch38_pct\n'; } > "$FLEET_DIR/mixed_partial.tsv"
    local f pat g w
    while IFS= read -r f; do
        pat="$(basename "$(dirname "$f")")"
        # colonies.tsv cols: donor proj ds readlen mapped assembly sample
        awk -F'\t' 'NR>2 && $6=="GRCh38" && $3 ~ /WGS/ && $7!="" && $2!="" {print $7"\t"$2}' "$f" \
            | sort -u > "$FLEET_DIR/worklist/$pat.samples.tsv"
        g="$(wc -l < "$FLEET_DIR/worklist/$pat.samples.tsv" | tr -d ' ')"
        w="$(awk -F'\t' 'NR>2 && $3 ~ /WGS/ && $7!="" {print $7}' "$f" | sort -u | wc -l | tr -d ' ')"
        if [ "$g" -eq 0 ]; then
            rm -f "$FLEET_DIR/worklist/$pat.samples.tsv"
            local asm; asm="$(awk -F'\t' 'NR>2{print $6}' "$f" | sort -u | paste -sd, -)"
            printf '%s\t%s\n' "$pat" "${asm:-none}" >> "$FLEET_DIR/skip_non_grch38.tsv"
        elif [ "$g" -lt "$w" ]; then
            printf '%s\t%s\t%s\t%d%%\n' "$pat" "$g" "$w" "$(( g * 100 / w ))" >> "$FLEET_DIR/mixed_partial.tsv"
            if [ "${INCLUDE_MIXED:-0}" = 1 ]; then echo "$pat" >> "$FLEET_DIR/patients.txt"
            else rm -f "$FLEET_DIR/worklist/$pat.samples.tsv"; fi
        else
            echo "$pat" >> "$FLEET_DIR/patients.txt"     # pure GRCh38
        fi
    done < <(find "$PATIENTS_DIR" -mindepth 3 -maxdepth 3 -name colonies.tsv | sort)

    sort -u -o "$FLEET_DIR/patients.txt" "$FLEET_DIR/patients.txt"
}

# ---- step 1: clean the scratch (staging) area --------------------------------
clean_scratch() {
    [ -d "$STAGING_ROOT" ] || { log "step1: staging root $STAGING_ROOT absent — nothing to clean"; return 0; }
    # safety: refuse if a live LSF job references the staging root
    if command -v bjobs >/dev/null 2>&1 && bjobs -w -u "$USER" 2>/dev/null | grep -q "$STAGING_ROOT"; then
        die "step1: live LSF jobs reference $STAGING_ROOT — refusing to wipe. Drain them first (bjobs -w)."
    fi
    local used; used="$(du -sh "$STAGING_ROOT" 2>/dev/null | cut -f1)"
    log "step1: wiping scratch/staging $STAGING_ROOT (was ${used:-?}); re-stageable from iRODS"
    if [ "${DRYRUN:-0}" = 1 ]; then log "step1: [dry-run] would rm -rf $STAGING_ROOT/*"; return 0; fi
    rm -rf "${STAGING_ROOT:?}/"* 2>/dev/null || true
    mkdir -p "$STAGING_ROOT"
}

# ---- step 2: lustre capacity gate --------------------------------------------
# The binding limit on farm22 is the TEAM (group) quota — the user quota is unlimited
# (limit 0). We parse `lfs quota -g <group>`: on that filesystem the mountpoint prints
# on its own line and the numbers wrap to the next; fields are
# used quota limit grace files ... so used=$1, limit=$3 on the numeric line.
capacity_gate() {
    local n_pat="$1"
    local peak_gb=$(( PATIENTS_INFLIGHT * STAGE_THROTTLE * AVG_BAM_GB ))   # bounded footprint
    local need_gb=$(( peak_gb * 115 / 100 ))
    log "step2: checking lustre capacity on $LUSTRE_VOL (group $LUSTRE_GROUP)"

    local qline used_gb limit_gb free_gb
    qline="$(lfs quota -g "$LUSTRE_GROUP" "$LUSTRE_VOL" 2>/dev/null | awk '/^[[:space:]]+[0-9]/{u=$1; l=$3; print u" "l; exit}')"
    if [ -n "$qline" ]; then
        used_gb=$(( $(echo "$qline" | cut -d' ' -f1) / 1024 / 1024 ))
        local limit_kb; limit_kb="$(echo "$qline" | cut -d' ' -f2)"
        if [ "$limit_kb" -gt 0 ] 2>/dev/null; then
            limit_gb=$(( limit_kb / 1024 / 1024 ))
            free_gb=$(( limit_gb - used_gb ))
            log "step2: team273 quota ${used_gb}GB used / ${limit_gb}GB limit -> ${free_gb}GB free ($(( used_gb*100/limit_gb ))% full)"
        else
            log "step2: group quota is unlimited (limit 0) — skipping hard gate"
        fi
    else
        free_gb="$(df -BG "$LUSTRE_VOL" 2>/dev/null | awk 'NR==2{gsub(/G/,"",$4);print $4}')"
        log "step2: lfs quota unparsed; df free ~${free_gb:-unknown}GB"
    fi

    log "step2: est. peak staging ~${peak_gb}GB (+15% -> ${need_gb}GB) for $n_pat patients (bounded by cleanup + PATIENTS_INFLIGHT)"
    if [ -n "${free_gb:-}" ] && [ "$free_gb" -lt "$need_gb" ] 2>/dev/null; then
        die "step2: insufficient lustre headroom (${free_gb}GB free < ${need_gb}GB needed). Lower PATIENTS_INFLIGHT/STAGE_THROTTLE or free space."
    fi
    log "step2: capacity OK"
}

# ---- steps 3-9 for the fleet: submit each patient DAG, throttled --------------
# The heavy work is child LSF jobs; this loop only submits + throttles, so it is
# cheap and can itself run under LSF (basement) for a walk-away run.
run_fleet() {
    local -a PATS=("$@")
    if [ "${#PATS[@]}" -eq 0 ]; then mapfile -t PATS < "$FLEET_DIR/patients.txt"; fi
    [ "${#PATS[@]}" -gt 0 ] || die "no GRCh38 patients to run"
    log "fleet: ${#PATS[@]} patient(s); max $PATIENTS_INFLIGHT staging concurrently"

    local pat sl
    for pat in "${PATS[@]}"; do
        sl="$FLEET_DIR/worklist/$pat.samples.tsv"
        [ -s "$sl" ] || { log "$pat: no worklist (skip)"; continue; }
        # throttle: wait until fewer than N patients have a live stage (sd) array
        while [ "$(bjobs -w -u "$USER" 2>/dev/null | grep -cE '_sd(\[|r)? ' || true)" -ge "$PATIENTS_INFLIGHT" ]; do
            sleep 60
        done
        log "$pat: submitting DAG ($(wc -l < "$sl" | tr -d ' ') GRCh38 colonies)"
        DISC_CFG="${DISC_CFG:-$PT_ROOT/cluster/config.discovery.grch38}" \
        GENO_CFG="${GENO_CFG:-$PT_ROOT/cluster/config.genotype.grch38}" \
            bash "$PIPELINE" submit-list "$pat" "$sl" || log "$pat: submit FAILED (continuing)"
    done
    log "fleet: all patients submitted. QC rolls up as they finish: $SELF status"
}

# ---- step 8: QC roll-up (mark, don't abort) ----------------------------------
qc_rollup() {
    local qc="$FLEET_DIR/fleet_qc.tsv"
    printf 'patient\tn_listed\tn_missing\tn_disc\tn_geno\tcontract\tcalls\tverdict\n' > "$qc"
    local pat rd n_s n_m n_d n_g contract calls verdict
    while IFS= read -r pat; do
        rd="$WORKROOT/$pat"
        [ -d "$rd" ] || continue
        n_s=$(wc -l < "$rd/samples.tsv" 2>/dev/null | tr -d ' '); n_s=${n_s:-0}
        n_m=$(ls "$rd/missing" 2>/dev/null | wc -l | tr -d ' ')
        n_d=$(ls "$rd"/discovery/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')
        n_g=$(ls "$rd"/genotypes/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')
        [ -s "$rd/insertions/$pat.genotyping.txt.gz" ] && contract=yes || contract=no
        [ -s "$rd/$pat.genotypes.csv.gz" ] && calls=yes || calls=no
        verdict=OK
        if ls "$rd"/FLEET_FATAL.* >/dev/null 2>&1; then verdict="FAIL:all-failed"
        elif [ "$calls" = no ]; then verdict="WARN:no-calls"
        elif [ "$n_s" -gt 0 ] && [ $(( n_m * 100 )) -ge $(( n_s * 5 )) ]; then
            verdict="WARN:missing=$(( n_m * 100 / n_s ))%"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$pat" "$n_s" "$n_m" "$n_d" "$n_g" "$contract" "$calls" "$verdict" >> "$qc"
    done < <( [ -s "$FLEET_DIR/patients.txt" ] && cat "$FLEET_DIR/patients.txt" || ls "$WORKROOT" )
    column -t -s $'\t' "$qc"
    echo "---"
    awk -F'\t' 'NR>1{c[$8 ~ /^OK/ ? "OK" : ($8 ~ /^WARN/ ? "WARN" : "FAIL")]++}
        END{printf "OK=%d  WARN=%d  FAIL=%d\n", c["OK"], c["WARN"], c["FAIL"]}' "$qc"
    log "QC written: $qc"
}

# ---- step 9 (final sweep): report leftover staged dirs ------------------------
staging_final_sweep() {
    [ -d "$STAGING_ROOT" ] || return 0
    local left; left="$(find "$STAGING_ROOT" -mindepth 2 -maxdepth 2 -type d 2>/dev/null | wc -l | tr -d ' ')"
    if [ "$left" -gt 0 ]; then
        log "step9: $left staged sample dir(s) remain under $STAGING_ROOT (patients not yet finished, or failures)."
        log "step9: per-patient cleanup runs on each success; run 'bash cluster/staging_clean.sh --delete' for a careful sweep once done."
    else
        log "step9: staging area clean."
    fi
}

# --- subcommands --------------------------------------------------------------
case "${1:-}" in
    plan)
        DRYRUN=1
        mkdir -p "$FLEET_DIR"
        build_worklist
        n=$(wc -l < "$FLEET_DIR/patients.txt" | tr -d ' ')
        skip=$(wc -l < "$FLEET_DIR/skip_non_grch38.tsv" | tr -d ' ')
        mixed=$(( $(wc -l < "$FLEET_DIR/mixed_partial.tsv" | tr -d ' ') - 1 ))
        col=$(cat "$FLEET_DIR/worklist/"*.samples.tsv 2>/dev/null | wc -l | tr -d ' ')
        log "plan: $n patient(s) to run (INCLUDE_MIXED=${INCLUDE_MIXED:-0}), $col colonies"
        log "plan: $mixed mixed-assembly patient(s) -> remap queue (mixed_partial.tsv); $skip non-GRCh38 -> skip_non_grch38.tsv"
        capacity_gate "$n" || true
        echo "--- patients to run ---"; cat "$FLEET_DIR/patients.txt"
        echo "--- mixed (partial tree; deferred to remap) ---"; column -t -s $'\t' "$FLEET_DIR/mixed_partial.tsv"
        ;;
    run)
        shift
        mkdir -p "$FLEET_DIR"
        build_worklist
        n=$(wc -l < "$FLEET_DIR/patients.txt" | tr -d ' ')
        clean_scratch                     # step 1
        capacity_gate "$n"                # step 2 (hard stop)
        run_fleet "$@"                    # steps 3-7,9 (per-patient DAGs) submitted
        staging_final_sweep               # step 9 report
        ;;
    status)
        [ -d "$FLEET_DIR" ] || die "no fleet dir $FLEET_DIR — run 'fleet.sh plan' or 'run' first"
        qc_rollup                         # step 8
        ;;
    *)
        sed -n '2,20p' "$SELF"; exit 1 ;;
esac
