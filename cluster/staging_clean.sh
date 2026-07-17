#!/bin/bash
# Classify staged project dirs by whether they are REDUNDANT with nst_links, and (only with
# --delete) remove the ones that provably are.
#
#   bash cluster/staging_clean.sh                 # dry run: report only, deletes nothing
#   bash cluster/staging_clean.sh --delete        # remove REDUNDANT + unreferenced dirs
#
# WHY THIS IS NOT JUST `rm -rf`. Staging holds ~8.3T of jd43's footprint against a team273
# quota at 90/100T, and every dir is 3+ months stale, so the cleanup is real. But three
# things make a blind rm wrong:
#
#   1. REDUNDANCY IS NOT UNIFORM. Staging came from iRODS, and some staged names have no
#      nst_links equivalent: PD44579b_lo0093_2 (a re-release), *_merged, *_10xil (10x
#      Genomics). Those may be the only local copy. Re-staging 3.6T is a slow undo.
#   2. SOMETHING MAY POINT AT IT. If a bams.fofn references /lustre/.../staging/..., then
#      deleting that dir silently breaks re-running that step. Checked below.
#   3. STAGING IS NOT ONLY THIS PROJECT'S. These dirs may serve other work. This script
#      only ever claims "redundant with nst_links", never "you don't need this".
#
# Recoverability: staged data originates in iRODS, so worst case is re-staging time, not
# data loss. Recoverable is not the same as free.
#
# Env: S (staging root), NST (nst_links), FOFNGLOB (files that might reference staging)
set -uo pipefail

S=${S:-/lustre/scratch126/casm/staging/team273/jd43}
NST=${NST:-/nfs/cancer_ref01/nst_links/live}
DELETE=0
[ "${1:-}" = "--delete" ] && DELETE=1

[ -d "$S" ] || { echo "no staging root: $S" >&2; exit 1; }
[ -d "$NST" ] || { echo "no nst_links: $NST" >&2; exit 1; }

# Anything that references staging by path makes its dir off-limits, whatever nst_links says.
echo "scanning for references to staging paths..."
REFS=$(grep -rl "$S" "$HOME" ~/scratch126_user 2>/dev/null | head -20)
if [ -n "$REFS" ]; then
    echo "  FILES REFERENCING STAGING (their projects will be kept):"
    printf '    %s\n' $REFS
fi

# nst_links index: sample -> exists. One pass, then O(1) lookups.
echo "indexing nst_links..."
declare -A HAVE
for projdir in "$NST"/*/; do
    for sdir in "$projdir"*/; do
        [ -d "$sdir" ] || continue
        s=$(basename "$sdir")
        [ -f "$sdir/$s.sample.dupmarked.bam" ] && HAVE["$s"]=1
    done
done
echo "  ${#HAVE[@]} samples with BAMs in nst_links"
echo

printf "%-6s %-7s %6s %6s  %s\n" PROJ SIZE SAMPLES INNST VERDICT
redundant=(); keep=()
for d in "$S"/*/; do
    [ -d "$d" ] || continue
    p=$(basename "$d")
    mapfile -t samples < <(ls "$d" 2>/dev/null)
    n=${#samples[@]}
    if [ "$n" -eq 0 ]; then
        printf "%-6s %-7s %6s %6s  %s\n" "$p" "empty" 0 0 "EMPTY -> rmdir"
        redundant+=("$d"); continue
    fi
    size=$(du -sh "$d" 2>/dev/null | cut -f1)
    inn=0; miss=()
    for s in "${samples[@]}"; do
        # try the staged name, and the name minus a _<n> re-release suffix
        base="${s%_[0-9]}"
        if [ -n "${HAVE[$s]:-}" ] || [ -n "${HAVE[$base]:-}" ]; then inn=$((inn+1)); else miss+=("$s"); fi
    done
    if printf '%s\n' $REFS | xargs -r grep -l "$p" 2>/dev/null | grep -q .; then
        v="KEEP: referenced by a fofn/script"
        keep+=("$d")
    elif [ "$inn" -eq "$n" ]; then
        v="REDUNDANT: all $n in nst_links"
        redundant+=("$d")
    else
        v="KEEP: ${#miss[@]} NOT in nst_links (e.g. ${miss[0]:-?})"
        keep+=("$d")
    fi
    printf "%-6s %-7s %6s %6s  %s\n" "$p" "$size" "$n" "$inn" "$v"
done

echo
echo "REDUNDANT (safe to remove): ${#redundant[@]} dirs"
echo "KEEP:                       ${#keep[@]} dirs"

if [ "$DELETE" -eq 0 ]; then
    echo
    echo "DRY RUN — nothing deleted. Re-run with --delete to remove the REDUNDANT dirs above."
    echo "Review the list first: this script only claims 'also present in nst_links'."
    echo "It does NOT know whether you need these for work other than the MEI benchmark."
    exit 0
fi

[ "${#redundant[@]}" -gt 0 ] || { echo "nothing to delete"; exit 0; }
echo
echo "DELETING ${#redundant[@]} redundant staging dirs in 10s — Ctrl-C to abort"
sleep 10
for d in "${redundant[@]}"; do
    echo "  rm -rf $d"
    rm -rf "$d"
done
echo "done. new usage:"
lfs quota -h -u "$(whoami)" /lustre/scratch126 2>/dev/null | head -4
