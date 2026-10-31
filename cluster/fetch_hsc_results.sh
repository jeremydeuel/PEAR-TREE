#!/bin/bash
# =============================================================================
# cluster/fetch_hsc_results.sh — run on the Mac: pull finished hsc runs from the farm's NFS results
# dir into ~/Documents/hsc_results, and print what is still missing.
#
#   bash cluster/fetch_hsc_results.sh            # copy every finished unit not yet local
#   FETCH_DRY=1 bash cluster/fetch_hsc_results.sh   # report only
#
# "Finished" = <U>.somatic.xlsx exists in $REMOTE_ROOT/<U> (the report step's last copy).
# A unit already present locally is re-copied only if the farm's somatic.xlsx is newer.
# Roster = run_unit column of patients/run_tracker.tsv.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PT_ROOT="$(cd "$HERE/.." && pwd)"
FARM_USER="${FARM_USER:-jd43}"
FARM_HOST="${FARM_HOST:-farm22-head1}"
REMOTE_ROOTS="${REMOTE_ROOTS:-/nfs/users/nfs_j/jd43/results/hsc_rediscover /nfs/users/nfs_j/jd43/results/hsc_foldback2 /nfs/users/nfs_j/jd43/results/hsc}"
LOCAL_ROOT="${LOCAL_ROOT:-$HOME/Documents/hsc_results}"
TSH="${TSH:-tsh}"

units=$(awk -F'\t' '!/^#/ && $1!="pd" && $6!="" {print $6}' "$PT_ROOT/patients/run_tracker.tsv" | sort -u)

# one remote round-trip: "<root>\t<unit>\t<mtime epoch>" for every somatic.xlsx
listing=$("$TSH" ssh "$FARM_USER@$FARM_HOST" "for r in $REMOTE_ROOTS; do for u in $(echo $units); do
    f=\$r/\$u/\$u.somatic.xlsx; [ -s \$f ] && printf '%s\t%s\t%s\n' \$r \$u \$(stat -c %Y \$f); done; done; true")

printf '%-28s %-10s %s\n' UNIT STATE SOURCE
for u in $units; do
    # newest remote copy across roots
    line=$(printf '%s\n' "$listing" | awk -F'\t' -v u="$u" '$2==u' | sort -t$'\t' -k3,3n | tail -1)
    lx="$LOCAL_ROOT/$u/$u.somatic.xlsx"
    if [ -z "$line" ]; then
        [ -s "$lx" ] && printf '%-28s %-10s %s\n' "$u" local "(not on farm NFS)" \
                     || printf '%-28s %-10s %s\n' "$u" MISSING "-"
        continue
    fi
    root=$(cut -f1 <<<"$line"); rt=$(cut -f3 <<<"$line")
    lt=$( [ -s "$lx" ] && stat -f %m "$lx" || echo 0)
    if [ "$lt" -ge "$rt" ]; then
        printf '%-28s %-10s %s\n' "$u" local "$root"
        continue
    fi
    state=$([ "$lt" -eq 0 ] && echo NEW || echo UPDATED)
    printf '%-28s %-10s %s\n' "$u" "$state" "$root"
    [ -n "${FETCH_DRY:-}" ] && continue
    "$TSH" scp -r "$FARM_USER@$FARM_HOST:$root/$u" "$LOCAL_ROOT/"
done
