#!/bin/bash
# Which patient/project combinations are hs37d5 WGS with enough colonies?
#
# Two naming conventions, two modes:
#
#   bash survey_assay.sh PD48402 PD45534 PD40521 ...          # dir-per-patient (PD*)
#   bash survey_assay.sh --list ax001.samples                 # bare sample names
#
# --list takes "PROJECT SAMPLE" lines, the format nst_links surveys print for donors whose
# samples are NOT under a <patient>* directory (AX001: `2315 T1_A3`). Samples are grouped
# by project, and up to PROBE per project get a header read.
#
# Run this BEFORE build_fofn.sh. The 10x10 benchmark's first attempt burned a full
# discovery run on three patients that were TARGETED_ILLUMINA_short gene panels: they
# passed the assembly/.bai/quickcheck checks and contributed 11, 7 and 92 loci to the
# pooled contract against 8,833 for a real WGS patient. Assay type is declared in @RG DS:
# and costs one header read — cheap, next to a wasted run.
set -uo pipefail

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

NST=${NST:-/nfs/cancer_ref01/nst_links/live}
MINCOL=${MINCOL:-10}       # colonies needed in one project to be worth using
PROBE=${PROBE:-2}          # BAMs header-read per patient/project (spot check, not a guarantee)
MIN_MAPPED_M=${MIN_MAPPED_M:-50}
EXCLUDE_RE='_merged$|_ds[0-9]|_ManyCrypts|_hum$'

hdr_of() { samtools view -H "$1" 2>/dev/null; }

# Verdict for one BAM: prints "DS<TAB>MAPPEDM<TAB>ASM<TAB>VERDICT"
judge() {
    local b="$1" hdr ds asm m v
    hdr=$(hdr_of "$b")
    [ -z "$hdr" ] && { printf 'UNREADABLE\t-\t-\tREJECT: no header\n'; return; }
    ds=$(echo "$hdr" | grep '^@RG' | tr '\t' '\n' | grep '^DS:' | sed 's/^DS://' | sort -u | paste -sd, -)
    if echo "$hdr" | grep -qm1 '^@SQ.*SN:hs37d5'; then asm=hs37d5
    else asm=$(echo "$hdr" | grep -m1 '^@SQ' | tr '\t' '\n' | grep -m1 '^SN:' | sed 's/^SN://'); fi
    # idxstats reads the INDEX only, so it is fast — but that also means it reports what the
    # index CLAIMS, not what the BAM holds. OBSERVED: PD34200c_lo0001 returns 0 mapped here
    # while `samtools view` still returns records. Cause unknown — its index is NOT older
    # than the BAM, so "stale index" is ruled out. Rejecting on a 0 is still right (a sample
    # whose index disagrees with its data breaks discovery too), but do not read this number
    # as ground truth for what is in the file.
    m=$(samtools idxstats "$b" 2>/dev/null | awk '{r+=$3} END {printf "%.0f", (r+0)/1000000}')
    # First reason wins — assay type is the root cause when several checks fail together
    # (a panel is both not-WGS and shallow; "not WGS" is what you act on).
    v=USE
    case "$ds" in WGS*) ;; *) v="REJECT: not WGS";; esac
    [ "$v" = USE ] && [ "$asm" != hs37d5 ] && v="REJECT: $asm not hs37d5"
    [ "$v" = USE ] && ! [ "${m:-0}" -ge "$MIN_MAPPED_M" ] 2>/dev/null && v="REJECT: only ${m:-0}M mapped"
    printf '%s\t%sM\t%s\t%s\n' "${ds:-NONE}" "${m:-0}" "$asm" "$v"
}

report() {  # report <label> <proj> <ncol> <bam>
    local ds m asm v
    IFS=$'\t' read -r ds m asm v < <(judge "$4")
    printf "%-12s %-6s %-5s %-26s %-8s %-8s %s\n" "$1" "$2" "$3" "$ds" "$m" "$asm" "$v"
}

printf "%-12s %-6s %-5s %-26s %-8s %-8s %s\n" DONOR PROJ NCOL DS MAPPED ASM VERDICT

if [ "${1:-}" = "--list" ]; then
    LIST="${2:?usage: $0 --list <file with 'PROJECT SAMPLE' lines>}"
    [ -s "$LIST" ] || { echo "no such list: $LIST" >&2; exit 1; }
    # group by project, preserving first-seen order
    for proj in $(awk 'NF>=2{print $1}' "$LIST" | sort -u); do
        mapfile -t samples < <(awk -v p="$proj" 'NF>=2 && $1==p{print $2}' "$LIST" | grep -Ev "$EXCLUDE_RE")
        n=${#samples[@]}
        done_n=0
        for s in "${samples[@]}"; do
            [ "$done_n" -ge "$PROBE" ] && break
            b="$NST/$proj/$s/$s.sample.dupmarked.bam"
            [ -s "$b" ] || continue
            report "$s" "$proj" "$n" "$b"
            done_n=$((done_n+1))
        done
        [ "$done_n" -gt 0 ] || printf "%-12s %-6s %-5s %-26s %-8s %-8s %s\n" "-" "$proj" "$n" "-" "-" "-" "no readable BAM"
    done
else
    [ $# -gt 0 ] || { echo "usage: $0 <patient> [...]  |  $0 --list <file>" >&2; exit 1; }
    for p in "$@"; do
        # Distinguish the two ways a donor yields nothing: no project has enough colonies,
        # versus projects have colonies but no readable BAM. Collapsing them into one
        # message sends you hunting for the wrong problem.
        best_n=0; probed=0
        for projdir in "$NST"/*/; do
            proj=$(basename "$projdir")
            mapfile -t cols < <(ls -d "$projdir/$p"* 2>/dev/null | xargs -r -n1 basename | grep -Ev "$EXCLUDE_RE")
            n=${#cols[@]}
            [ "$n" -gt "$best_n" ] && best_n=$n
            [ "$n" -ge "$MINCOL" ] || continue
            done_n=0
            for s in "${cols[@]}"; do
                [ "$done_n" -ge "$PROBE" ] && break
                b="$projdir/$s/$s.sample.dupmarked.bam"
                [ -s "$b" ] || continue
                report "$p" "$proj" "$n" "$b"
                done_n=$((done_n+1)); probed=$((probed+1))
            done
            [ "$done_n" -gt 0 ] || printf "%-12s %-6s %-5s %-26s %-8s %-8s %s\n" \
                "$p" "$proj" "$n" - - - "no readable .sample.dupmarked.bam"
        done
        [ "$probed" -gt 0 ] || [ "$best_n" -ge "$MINCOL" ] || \
            printf "%-12s %-6s %-5s %-26s %-8s %-8s %s\n" "$p" - "$best_n" - - - \
                   "no project with >=$MINCOL colonies (best: $best_n)"
    done
fi

echo
echo "Feed only USE rows to build_fofn.sh as <patient>:<project>."
echo "PROBE=$PROBE BAMs checked per donor/project — build_fofn.sh re-checks every BAM it picks."
