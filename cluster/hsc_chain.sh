#!/bin/bash
# =============================================================================
# cluster/hsc_chain.sh — submit several patients now, run them one after another (farm22).
#
#   bash <checkout>/cluster/hsc_chain.sh <UNIT> [<UNIT> ...]
#
# Each unit is an ordinary `hsc_run.sh submit`, but its stage+discover array waits (PT_SD_WAIT)
# for ended(<previous unit's cleanup job>): the next patient stages BAMs only once the previous
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
# If a patient fails before combine_genotypes its cleanup job never becomes runnable and the
# chain stops there (its BAMs are still staged). `bkill <that cleanup jobid>` releases the next
# unit (ended() includes EXIT); the summary table lists every cleanup job id.
# Writes the table to $HSC_ROOT/chain.<timestamp>.tsv.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JD="${JD:-/lustre/scratch126/casm/teams/team273/users/jd43}"
HSC_ROOT="${HSC_ROOT:-$JD/hsc}"
export HSC_ROOT

die()  { echo "hsc_chain: $*" >&2; exit 1; }
note() { echo "[$(date +%H:%M:%S)] $*" >&2; }

[ "$#" -gt 0 ] || { sed -n '2,23p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1; }
[ "${KEEP_BAMS:-0}" != 1 ] || die "KEEP_BAMS=1 submits no cleanup job, so there is nothing to chain on"

TABLE="$HSC_ROOT/chain.$(date +%Y%m%d%H%M%S).tsv"
[ "${HSC_CHAIN_DRY:-0}" = 1 ] && TABLE=/dev/null
printf 'unit\tstatus\tsd_jobid\tcleanup_jobid\twaits_on\n' > "$TABLE"

prev="${HSC_CHAIN_AFTER:-}"
for U in "$@"; do
    W="$HSC_ROOT/$U"
    if [ -s "$W/jobids.tsv" ]; then
        note "$U: already submitted ($W/jobids.tsv) -- skipped"
        printf "%s\tskipped (already submitted)\t-\t-\t-\n" "$U" >> "$TABLE"; continue
    fi
    if [ "${HSC_CHAIN_DRY:-0}" = 1 ]; then
        note "$U: would submit, stage+discover waiting on ${prev:-nothing}"; prev="<$U cleanup>"; continue
    fi

    log="$(mktemp)"; rc=0
    note "$U: submitting${prev:+ (stage+discover waits on ended($prev))}"
    PT_SD_WAIT="${prev:+ended($prev)}" bash "$HERE/hsc_run.sh" submit "$U" 2>&1 | tee "$log" || rc=$?
    if [ "$rc" -ne 0 ] && grep -q 'rerun with PT_RESTAGE=1' "$log" && [ "${PT_RESTAGE:-0}" != 1 ]; then
        note "$U: refused by the capacity check -- retrying with PT_RESTAGE=1"; rc=0
        PT_RESTAGE=1 PT_SD_WAIT="${prev:+ended($prev)}" bash "$HERE/hsc_run.sh" submit "$U" 2>&1 | tee "$log" || rc=$?
    fi
    rm -f "$log"
    if [ "$rc" -ne 0 ]; then
        note "$U: SUBMIT FAILED (rc=$rc) -- skipped; the next unit waits on ${prev:-nothing}"
        printf "%s\tSUBMIT FAILED\t-\t-\t%s\n" "$U" "${prev:--}" >> "$TABLE"; continue
    fi

    sd="$(awk -F'\t' '$1=="sd" {j=$2} END {print j}' "$W/jobids.tsv")"
    cu="$(awk -F'\t' '$1=="cu" {j=$2} END {print j}' "$W/jobids.tsv")"
    [ -n "$cu" ] || die "$U: submitted but $W/jobids.tsv has no cleanup job -- cannot chain the rest"
    printf '%s\tsubmitted\t%s\t%s\t%s\n' "$U" "${sd:--}" "$cu" "${prev:--}" >> "$TABLE"
    prev="$cu"
done

[ "${HSC_CHAIN_DRY:-0}" = 1 ] && exit 0
note "chain table: $TABLE"
column -t -s $'\t' "$TABLE" >&2
