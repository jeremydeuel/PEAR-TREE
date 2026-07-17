#!/bin/bash
# Build the 10-colonies-per-patient fofn for the MEI benchmark.
#
#   bash build_fofn.sh PD41048:2073 PD34200:1748 PD37449:1743 PD43947:2244 PD43974:2567
#
# *** SUPERSEDED for the GRCh38 9x10 benchmark -- use cluster/stage_picked.sh instead. ***
#
# *** THIS SCRIPT VIOLATES SANGER POLICY AS WRITTEN. *** It emits nst_links paths into a
# fofn that discovery then reads END-TO-END. Policy: BAM records are ALWAYS read from a
# STAGED Lustre copy, NEVER from nst_links -- those are symlinks into the iRODS resource
# servers, and streaming terabytes off them hammers shared archive infrastructure. Only
# HEADER reads (view -H / idxstats, as used by the gates below) are allowed there.
#
# It survives because its per-BAM gates (assay, assembly, depth, .bai, truncation) are sound
# and worth keeping. If you resurrect it for a new cohort, make it emit STAGED paths --
# resolve via nst_links/the catalogue, stage with stageBam.pl, then write the fofn from
# $STAGE/$PROJ/$SAMPLE/... See cluster/stage_picked.sh for the shape that does this right.
#
# Takes explicit <patient>:<project> pairs, NOT a bare patient. This is deliberate:
# nst_links re-releases the same colonies under several project ids and the assembly
# differs between (and even within) those releases -- project 1903 holds PD41048 on
# hs37d5 and PD40667 on GRCh38. Globbing live/*/<patient>* therefore mixes assemblies
# and double-counts genomes. Pick the project from cluster/survey_projects.sh output.
#
# Every BAM is still individually checked for hs37d5 / .bai / truncation: the project
# is a hint, never a guarantee.
#
# Writes ~/mei10x10/<patient>.bams.fofn and appends ~/mei10x10/all.bams.fofn.
set -uo pipefail

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

NST=${NST:-/nfs/cancer_ref01/nst_links/live}
N=${N:-10}                        # colonies per patient
OUT=${OUT:-$HOME/mei10x10}
mkdir -p "$OUT"
: > "$OUT/all.bams.fofn"

[ $# -gt 0 ] || { echo "usage: $0 <patient>:<project> [...]" >&2; exit 1; }

is_hs37d5() { samtools view -H "$1" 2>/dev/null | grep -qm1 '^@SQ.*SN:hs37d5'; }

# ASSAY GATE. MEASURED, and the reason the first 10x10 attempt wasted a full run:
# PD34200/PD37449/PD43947 passed every check below (hs37d5, .bai, quickcheck) and were
# TARGETED_ILLUMINA_short — gene panels, 3-13M mapped reads. Their pooled contract
# contribution was 11, 7 and 92 loci against 8,833 for a real WGS patient (PD41048, ~480M
# reads). That is not a weak signal, it is the correct answer to the wrong question: MEI
# discovery reads soft-clips at breakpoints, and a panel covers a sliver of the genome, so
# there is nowhere for a genome-wide breakpoint to appear. Every human carries >1000
# polymorphic MEIs, so <100 loci from 10 colonies means the INPUT is wrong, never biology.
# The assay is declared in @RG DS: — check it, cheaply, before burning a discovery run.
is_wgs() { samtools view -H "$1" 2>/dev/null | grep '^@RG' | tr '\t' '\n' | grep -qm1 '^DS:WGS'; }

# Depth backstop for the case DS: is absent or lies. WGS colonies here run 166-510M mapped
# reads; targeted ran 3-13M. 50M separates them by an order of magnitude on both sides.
# idxstats reads the INDEX only, so this is fast — but it also means it reports what the
# index CLAIMS, not what the BAM holds. OBSERVED: PD34200c_lo0001 returns 0 mapped here
# while `samtools view` still returns records; cause unknown (its index is NOT older than
# the BAM, so staleness is ruled out). Rejecting on 0 is still correct — a sample whose
# index disagrees with its data breaks discovery too — but if this gate ever rejects a
# colony you believe is good WGS, suspect the metric before suspecting the data.
MIN_MAPPED=${MIN_MAPPED:-50000000}
mapped_reads() { samtools idxstats "$1" 2>/dev/null | awk '{r+=$3} END {print r+0}'; }

# Not colonies: merged/bulk/downsampled/xenograft-filtered derivatives.
EXCLUDE_RE='_merged$|_ds[0-9]|_ManyCrypts|_hum$'

for SPEC in "$@"; do
    P=${SPEC%%:*}; PROJ=${SPEC##*:}
    if [ "$P" = "$PROJ" ]; then echo "SKIP $SPEC: need <patient>:<project>" >&2; continue; fi
    DIR="$NST/$PROJ"
    [ -d "$DIR" ] || { echo "SKIP $SPEC: no project dir $DIR" >&2; continue; }

    # candidate samples in THIS project only; prefer _lo colonies when the patient has them
    mapfile -t ALL < <(ls -d "$DIR/${P}"* 2>/dev/null | xargs -r -n1 basename | sort -u | grep -Ev "$EXCLUDE_RE")
    mapfile -t LO  < <(printf '%s\n' "${ALL[@]:-}" | grep -E '_lo[0-9]+' || true)
    if [ "${#LO[@]}" -gt 0 ]; then SAMPLES=("${LO[@]}"); NOTE="_lo colonies"
    else                          SAMPLES=("${ALL[@]:-}"); NOTE="NO _lo naming - verify these are colonies"; fi
    total=${#SAMPLES[@]}
    if [ "$total" -eq 0 ]; then echo "SKIP $SPEC: no samples" >&2; continue; fi

    # even spread across the sorted list so one plate cannot dominate
    picked=(); tried=0
    while [ "${#picked[@]}" -lt "$N" ] && [ "$tried" -lt "$total" ]; do
        idx=$(( tried * total / N )); [ "$idx" -ge "$total" ] && idx=$((total-1))
        s="${SAMPLES[$idx]}"; tried=$((tried+1))
        b="$DIR/$s/$s.sample.dupmarked.bam"
        printf '%s\n' "${picked[@]:-}" | grep -qxF "$b" 2>/dev/null && continue
        [ -s "$b" ]     || { echo "  $P/$PROJ: no bam   $s" >&2; continue; }
        [ -s "$b.bai" ] || { echo "  $P/$PROJ: no .bai  $s" >&2; continue; }
        samtools quickcheck "$b" 2>/dev/null || { echo "  $P/$PROJ: truncated $s" >&2; continue; }
        is_hs37d5 "$b" || { echo "  $P/$PROJ: NOT hs37d5 $s" >&2; continue; }
        is_wgs "$b"    || { echo "  $P/$PROJ: NOT WGS ($(samtools view -H "$b" 2>/dev/null | grep -m1 '^@RG' | tr '\t' '\n' | grep '^DS:')) $s" >&2; continue; }
        m=$(mapped_reads "$b")
        [ "$m" -ge "$MIN_MAPPED" ] || { echo "  $P/$PROJ: only $((m/1000000))M mapped reads (need $((MIN_MAPPED/1000000))M) $s" >&2; continue; }
        picked+=("$b")
    done
    # top up sequentially if the spread lost some to the checks
    if [ "${#picked[@]}" -lt "$N" ]; then
        for s in "${SAMPLES[@]}"; do
            [ "${#picked[@]}" -ge "$N" ] && break
            b="$DIR/$s/$s.sample.dupmarked.bam"
            printf '%s\n' "${picked[@]:-}" | grep -qxF "$b" 2>/dev/null && continue
            # same gates as above — a laxer top-up would re-admit exactly what the main
            # loop just rejected, and silently
            [ -s "$b" ] && [ -s "$b.bai" ] && samtools quickcheck "$b" 2>/dev/null \
                && is_hs37d5 "$b" && is_wgs "$b" \
                && [ "$(mapped_reads "$b")" -ge "$MIN_MAPPED" ] \
                && picked+=("$b")
        done
    fi

    if [ "${#picked[@]}" -lt "$N" ]; then
        echo "SKIP $P: only ${#picked[@]}/$N usable hs37d5 BAMs in project $PROJ" >&2
        continue
    fi
    printf '%s\n' "${picked[@]}" > "$OUT/$P.bams.fofn"
    printf '%s\n' "${picked[@]}" >> "$OUT/all.bams.fofn"
    printf "OK   %-9s proj %-5s %2d BAMs (of %3d %s) -> %s\n" \
        "$P" "$PROJ" "${#picked[@]}" "$total" "$NOTE" "$OUT/$P.bams.fofn"
done

echo
echo "total: $(wc -l < "$OUT/all.bams.fofn") BAMs -> $OUT/all.bams.fofn"
echo "sanity: every path below must be a DISTINCT genome (no colony twice, all hs37d5)"
