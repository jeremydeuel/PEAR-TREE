#!/bin/bash
# Probe the new-organ candidates: map MANIFEST samples -> nst_links project, then header-read
# a few per (donor x project) to get assay type, read length, assembly and depth.
#
#   bash cluster/probe_new_organs.sh cluster/manifests/oliver2025nf1.wgs.tsv
#   bash cluster/probe_new_organs.sh cluster/manifests/coorens2025stomach.wgs.tsv
#   PROBE=4 bash cluster/probe_new_organs.sh <manifest>     # more spot checks per group
#
# WHY NOT survey_assay.sh: it hard-REJECTs anything that is not hs37d5, which is correct for
# the old GRCh37 cohorts and WRONG here — the NF1 WGS release (project 2571) is GRCh38 and
# would be rejected wholesale. This script REPORTS assembly instead of judging on it; we
# remap to hs1 anyway, so source assembly is triage, not a gate.
#
# WHY MANIFEST-DRIVEN, not a donor glob. Measured 2026-07-17 on the stomach cohort: a donor
# glob returns ~87 distinct samples for PD40293 but the paper only ran **12** of them as WGS
# — the rest are the 829-sample TARGETED PANEL, and the panel samples use the SAME
# `PD40293c_lo0003` naming as the WGS ones. You CANNOT tell them apart by name. You can only
# tell them apart by which supplementary table they appear in. Good news: the WGS and TGS
# name sets are DISJOINT (238 vs 829, overlap 0), so selecting by manifest name is exact.
# survey_assay.sh's header says a previous 10x10 run burned a full discovery on three
# TARGETED_ILLUMINA panels that contributed 11/7/92 loci vs 8,833 for real WGS. Same trap.
#
# READ LENGTH IS THE THING THAT CAN KILL THIS. 75bp vs 151bp materially changes clip-based
# MEI discovery and remapping cannot fix short reads. Costs ~1000 records to measure.
set -uo pipefail

MAN="${1:?usage: $0 <manifest.tsv>}"
[ -s "$MAN" ] || { echo "no such manifest: $MAN" >&2; exit 1; }
NST=${NST:-/nfs/cancer_ref01/nst_links/live}
PROBE=${PROBE:-2}

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

echo "manifest: $MAN"
echo "resolving manifest samples -> nst_links projects (this is the slow bit)..."
# donor <TAB> sample  -> donor <TAB> sample <TAB> project(s)
tail -n +2 "$MAN" | awk -F'\t' 'NF>=2{print $1"\t"$2}' | sort -u > "$TMP/man.tsv"
: > "$TMP/found.tsv"
while IFS=$'\t' read -r d s; do
    for p in $(ls -d "$NST"/*/"$s" 2>/dev/null | awk -F/ '{print $(NF-1)}'); do
        printf '%s\t%s\t%s\n' "$d" "$s" "$p" >> "$TMP/found.tsv"
    done
done < "$TMP/man.tsv"

tot=$(wc -l < "$TMP/man.tsv" | tr -d ' ')
got=$(cut -f1,2 "$TMP/found.tsv" | sort -u | wc -l | tr -d ' ')
echo
echo "manifest samples: $tot   resolved in nst_links: $got   unresolved: $((tot-got))"
echo "(unresolved != absent — nst_links is an incomplete view of iRODS; check with iquest)"

# The per-sample map is the ACTUAL artefact. Counts are a summary of it and summaries lie:
# on NF1 the glob count for PD51122 in proj 2571 was 413, exactly equal to the manifest's
# 413, which looked like proof the release was complete. It was not — the manifest-restricted
# count is 410, and 3 samples live only in proj 2789. Two aggregates agreed while disagreeing
# member-by-member. Always select on this map, never on "the project whose total looks right".
MAP="${MAP:-$HOME/$(basename "$MAN" .tsv).sample_project.tsv}"
{ echo -e "donor\tsample\tproject"; sort -u "$TMP/found.tsv"; } > "$MAP"
echo "per-sample map written: $MAP"

echo
echo "### manifest samples per (donor x project) — a SUMMARY of the map above, not a selector"
awk -F'\t' '{print $1"\t"$3}' "$TMP/found.tsv" | sort | uniq -c \
    | awk '{printf "  %-9s proj %-6s %s manifest samples\n", $2, $3, $1}'

# Samples reachable from exactly ONE project: these are the ones a "pin to project X"
# selection silently drops.
echo
echo "### samples that exist in ONLY ONE project (a single-project selection would drop these)"
cut -f2,3 "$TMP/found.tsv" | sort -u | cut -f1 | sort | uniq -c | awk '$1==1{print $2}' > "$TMP/single.txt"
if [ -s "$TMP/single.txt" ]; then
    join -1 1 -2 2 <(sort "$TMP/single.txt") <(sort -k2,2 "$TMP/found.tsv" | awk '{print $1"\t"$2"\t"$3}' | sort -k2,2) 2>/dev/null \
        | awk '{print $3}' | sort | uniq -c | awk '{printf "  proj %-6s is the ONLY home of %s manifest sample(s)\n", $2, $1}'
else
    echo "  (none — every manifest sample is reachable from more than one project)"
fi

echo
echo "### header probe: $PROBE BAM(s) per (donor x project)"
printf "%-9s %-6s %-24s %-8s %-7s %-14s %s\n" DONOR PROJ DS READLEN MAPPED ASSEMBLY SAMPLE
echo "--------------------------------------------------------------------------------------"
awk -F'\t' '{print $1"\t"$3}' "$TMP/found.tsv" | sort -u | while IFS=$'\t' read -r d p; do
    n=0
    awk -F'\t' -v d="$d" -v p="$p" '$1==d && $3==p{print $2}' "$TMP/found.tsv" | while read -r s; do
        [ "$n" -ge "$PROBE" ] && break
        b="$NST/$p/$s/$s.sample.dupmarked.bam"
        [ -s "$b" ] || continue
        hdr=$(samtools view -H "$b" 2>/dev/null)
        [ -z "$hdr" ] && continue
        ds=$(printf '%s\n' "$hdr" | grep '^@RG' | tr '\t' '\n' | grep '^DS:' | sed 's/^DS://' | sort -u | paste -sd, -)
        if printf '%s\n' "$hdr" | grep -qm1 '^@SQ.*SN:hs37d5'; then asm=hs37d5_GRCh37
        else
            c1=$(printf '%s\n' "$hdr" | grep -m1 '^@SQ.*SN:chr1[[:space:]]' | tr '\t' '\n' | grep -m1 '^LN:' | sed 's/^LN://')
            case "$c1" in 248956422) asm=GRCh38 ;; 249250621) asm=hg19 ;; *) asm="other(${c1:-?})" ;; esac
        fi
        # read length: needs RECORDS, not just the header. ~1000 reads is cheap and is the
        # ONLY way to get this — 75bp vs 151bp decides whether this cohort is usable at all.
        rl=$(samtools view "$b" 2>/dev/null | head -1000 \
             | awk '{print length($10)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
        m=$(samtools idxstats "$b" 2>/dev/null | awk '{r+=$3} END {printf "%.0f", (r+0)/1000000}')
        printf "%-9s %-6s %-24s %-8s %-7s %-14s %s\n" "$d" "$p" "${ds:-NONE}" "${rl:-?}bp" "${m:-0}M" "$asm" "$s"
        n=$((n+1))
    done
done

cat <<'EOF'

### READING THIS — this is INVENTORY, it does not select a cohort
  DS must say WGS. A TARGETED_ILLUMINA / panel DS means the sample is NOT usable for
  insertion discovery, no matter how good the depth looks.
  READLEN 151bp = good. 75bp = clip-based MEI discovery is badly compromised and remapping
  cannot fix it. A donor can be on the farm, correctly assembled, and still worthless.
  ASSEMBLY is triage only (we remap to hs1); but note it, and re-check the @SQ of every BAM
  you actually use — assembly is a property of the BAM, not the project.

  Record the answer in cluster/trees/INVENTORY.md. Choosing what to progress is a SEPARATE
  decision, made after the inventory is complete — not implied by anything printed above.
EOF
