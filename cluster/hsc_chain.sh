#!/bin/bash
# =============================================================================
# cluster/hsc_chain.sh — submit several patients now, run them one after another (farm22).
#
#   bash <checkout>/cluster/hsc_chain.sh <UNIT> [<UNIT> ...]
#
# Each unit is an ordinary `hsc_run.sh submit`, but its stage+discover array waits (PT_SD_WAIT)
# for ended(<previous unit's gate job>, see below): the next patient stages BAMs only once the previous
# one's staged BAMs are deleted, so at most one patient's BAMs sit on Lustre. Cleanup hangs off
# done(combine_genotypes) and runs beside annotate, so annotate/report of patient i overlap
# discovery of patient i+1.
#
#   * a unit whose $HSC_ROOT/<U>/jobids.tsv exists is skipped (already submitted)
#   * a submit refused by the Lustre capacity check is retried once with PT_RESTAGE=1
#   * a submit that still fails is reported and skipped; the chain continues from the last good unit
#   * HSC_CHAIN_AFTER=<jobid>  make the first unit wait for ended(<jobid>) too
#   * HSC_CHAIN_DRY=1          list what would be submitted, submit nothing
#
# Each unit also gets a tiny GATE job (<U>_hsc_gate) and the next unit waits on ended(gate), not
# on cleanup. The gate fires when cleanup ends OR any controller step exits (sdr, ci, gtr, cg):
#   * cleanup DONE        -> unit succeeded, next unit starts
#   * a step EXITed       -> writes $HSC_ROOT/<U>/CHAIN_FAILED, bkills the unit's now-orphaned
#                            PEND jobs, deletes its staged BAMs (pipeline.sh cleanup; ci/cg/annotate
#                            reruns need no BAMs) and releases the next unit
#   * HSC_CHAIN_ON_FAIL=stop  keep the old behaviour: gate exits 1, next unit waits on done(gate),
#                            the chain halts with the failed unit's BAMs still staged
# Before this (2026-10-25) a dead sdr/ci/gtr/cg left done()-dependent jobs pending forever and the
# whole chain sat in PEND for days.
#
#   bash <checkout>/cluster/hsc_chain.sh --rearm <chain.tsv>
#       retrofit gates onto an already-submitted chain: every unit whose cleanup is still pending
#       gets a gate, and each next unit's stage array is bmod'ed onto ended(previous gate).
# A step you resubmit by hand must be appended to jobids.tsv under the SAME phase name (sdr, ci,
# gtr, cg; cg2 etc. also accepted): the gate and --rearm use the LAST job id per phase.
# Writes the table to $HSC_ROOT/chain.<timestamp>.tsv.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JD="${JD:-/lustre/scratch126/casm/teams/team273/users/jd43}"
HSC_ROOT="${HSC_ROOT:-$JD/hsc}"
export HSC_ROOT

die()  { echo "hsc_chain: $*" >&2; exit 1; }
note() { echo "[$(date +%H:%M:%S)] $*" >&2; }

# last job id recorded for a phase (phase, phase2, phase3... all count) in a unit's jobids.tsv
jid_of() { awk -F'\t' -v p="$2" '$1 ~ ("^" p "[0-9]*$") {j=$2} END {print j}' "$1"; }
stat_of() { [ -n "$1" ] || return 0; bjobs -a -noheader -o stat "$1" 2>/dev/null | head -1 | tr -d ' ' || true; }

submit_gate() {   # submit_gate <U> -> prints the gate job id, records it as "gate" in jobids.tsv
    local U="$1" J="$HSC_ROOT/$1/jobids.tsv" dep p j out gid
    dep="ended($(jid_of "$J" cu))"
    for p in sdr ci gtr cg; do
        j="$(jid_of "$J" "$p")"; [ -n "$j" ] && dep="$dep || exit($j)"
    done
    out="$(bsub -J "${U}_hsc_gate" -q "${HSC_GATE_QUEUE:-normal}" -n 1 -M 500 -R "select[mem>500] rusage[mem=500]" \
        -w "$dep" -o "$HSC_ROOT/$U/chain_gate.%J.log" -e "$HSC_ROOT/$U/chain_gate.%J.log" \
        "HSC_ROOT='$HSC_ROOT' HSC_CHAIN_ON_FAIL='${HSC_CHAIN_ON_FAIL:-continue}' bash '$HERE/hsc_chain.sh' --gate '$U'" 2>&1)" \
        || die "$U: gate bsub failed: $out"
    gid="$(sed -n 's/^Job <\([0-9]*\)>.*/\1/p' <<<"$out")"
    [ -n "$gid" ] || die "$U: gate bsub gave no job id: $out"
    printf 'gate\t%s\n' "$gid" >> "$J"
    note "$U: gate $gid waits on $dep"
    echo "$gid"
}

next_wait() {     # dependency the NEXT unit's stage array waits on
    if [ "${HSC_CHAIN_ON_FAIL:-continue}" = stop ]; then echo "done($1)"; else echo "ended($1)"; fi
}

# --- gate (runs as an LSF job) ---------------------------------------------------------------
if [ "${1:-}" = --gate ]; then
    U="${2:?--gate <UNIT>}"; J="$HSC_ROOT/$U/jobids.tsv"; R="$HSC_ROOT/$U/$U"
    cu="$(jid_of "$J" cu)"; cs="$(stat_of "$cu")"
    if [ "$cs" = DONE ] || { [ -z "$cs" ] && [ -s "$R/$U.genotypes.csv.gz" ]; }; then
        note "$U: cleanup $cu ${cs:-gone, genotypes present} -> unit OK, releasing the chain"; exit 0
    fi
    failed=""
    for p in sdr ci gtr cg cu; do
        j="$(jid_of "$J" "$p")"; [ "$(stat_of "$j")" = EXIT ] && failed="$failed $p($j)"
    done
    { echo "$(date '+%F %T') unit failed at:${failed:- unknown (cleanup $cu: ${cs:-not found})}"
      echo "logs: $R/logs"; } | tee "$HSC_ROOT/$U/CHAIN_FAILED" >&2
    if [ "${HSC_CHAIN_ON_FAIL:-continue}" = stop ]; then
        note "$U: HSC_CHAIN_ON_FAIL=stop -> chain halts here (bkill this gate's dependents to release it)"; exit 1
    fi
    # the rest of this unit can never run now: kill its PEND jobs so they don't sit in the queue
    # (not when only cleanup itself failed: annotate/report hang off cg and are still valid)
    only_cu=0; [ -n "$failed" ] && [ -z "$(sed 's/ cu([0-9]*)//' <<<"$failed")" ] && only_cu=1
    [ "$only_cu" = 1 ] || cut -f1 "$J" | sed 's/[0-9]*$//' | sort -u | while read -r p; do
        [ "$p" = gate ] && continue
        j="$(jid_of "$J" "$p")"; [ "$(stat_of "$j")" = PEND ] && { bkill "$j" >&2 || true; }
    done || true
    PT_RUNDIR="$R" bash "$HERE/pipeline.sh" cleanup || note "$U: staged-BAM cleanup failed (check by hand)"
    note "$U: failed unit released; the next unit starts"; exit 0
fi

# --- rearm an already-submitted chain ------------------------------------------------------
if [ "${1:-}" = --rearm ]; then
    T="${2:?--rearm <chain.tsv>}"; [ -s "$T" ] || die "no chain table $T"
    prev=""
    while IFS=$'\t' read -r U status _; do
        [ "$U" = unit ] && continue
        J="$HSC_ROOT/$U/jobids.tsv"; [ -s "$J" ] || { note "$U: no jobids.tsv -- skipped"; continue; }
        cs="$(stat_of "$(jid_of "$J" cu)")"
        if [ "$cs" != PEND ]; then note "$U: cleanup ${cs:-gone} -- not part of the live chain"; prev=""; continue; fi
        sd="$(jid_of "$J" sd)"
        if [ -n "$prev" ]; then
            if [ "$(stat_of "$sd")" = PEND ]; then
                bmod -w "$(next_wait "$prev")" "$sd" >&2 || die "$U: bmod -w on stage array $sd failed"
                note "$U: stage array $sd now waits on $(next_wait "$prev")"
            else
                note "$U: stage array $sd already $(stat_of "$sd") -- dependency left as is"
            fi
        fi
        g="$(jid_of "$J" gate)"
        if [ -n "$g" ] && [ "$(stat_of "$g")" = PEND ]; then note "$U: gate $g already armed"; else g="$(submit_gate "$U")"; fi
        prev="$g"
    done < "$T"
    exit 0
fi

[ "$#" -gt 0 ] || { awk 'NR>2 && /^# ====/ {exit} NR>2' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1; }
[ "${KEEP_BAMS:-0}" != 1 ] || die "KEEP_BAMS=1 submits no cleanup job, so there is nothing to chain on"

TABLE="$HSC_ROOT/chain.$(date +%Y%m%d%H%M%S).tsv"
[ "${HSC_CHAIN_DRY:-0}" = 1 ] && TABLE=/dev/null
printf 'unit\tstatus\tsd_jobid\tcleanup_jobid\twaits_on_gate\n' > "$TABLE"

prev="${HSC_CHAIN_AFTER:-}"
for U in "$@"; do
    W="$HSC_ROOT/$U"
    if [ -s "$W/jobids.tsv" ]; then
        note "$U: already submitted ($W/jobids.tsv) -- skipped"
        printf "%s\tskipped (already submitted)\t-\t-\t-\n" "$U" >> "$TABLE"; continue
    fi
    if [ "${HSC_CHAIN_DRY:-0}" = 1 ]; then
        note "$U: would submit, stage+discover waiting on ${prev:-nothing}"; prev="<$U gate>"; continue
    fi

    log="$(mktemp)"; rc=0
    note "$U: submitting${prev:+ (stage+discover waits on $(next_wait "$prev"))}"
    PT_SD_WAIT="${prev:+$(next_wait "$prev")}" bash "$HERE/hsc_run.sh" submit "$U" 2>&1 | tee "$log" || rc=$?
    if [ "$rc" -ne 0 ] && grep -q 'rerun with PT_RESTAGE=1' "$log" && [ "${PT_RESTAGE:-0}" != 1 ]; then
        note "$U: refused by the capacity check -- retrying with PT_RESTAGE=1"; rc=0
        PT_RESTAGE=1 PT_SD_WAIT="${prev:+$(next_wait "$prev")}" bash "$HERE/hsc_run.sh" submit "$U" 2>&1 | tee "$log" || rc=$?
    fi
    rm -f "$log"
    if [ "$rc" -ne 0 ]; then
        note "$U: SUBMIT FAILED (rc=$rc) -- skipped; the next unit waits on ${prev:-nothing}"
        printf "%s\tSUBMIT FAILED\t-\t-\t%s\n" "$U" "${prev:--}" >> "$TABLE"; continue
    fi

    sd="$(awk -F'\t' '$1=="sd" {j=$2} END {print j}' "$W/jobids.tsv")"
    cu="$(awk -F'\t' '$1=="cu" {j=$2} END {print j}' "$W/jobids.tsv")"
    [ -n "$cu" ] || die "$U: submitted but $W/jobids.tsv has no cleanup job -- cannot chain the rest"
    gate="$(submit_gate "$U")"
    printf '%s\tsubmitted\t%s\t%s\t%s\n' "$U" "${sd:--}" "$cu" "${prev:--}" >> "$TABLE"
    prev="$gate"
done

[ "${HSC_CHAIN_DRY:-0}" = 1 ] && exit 0
note "chain table: $TABLE"
column -t -s $'\t' "$TABLE" >&2
