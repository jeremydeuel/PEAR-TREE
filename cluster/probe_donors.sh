#!/bin/bash
# Probe a list of DONORS by globbing nst_links (no manifest needed): for each donor, find
# every (project x sample), then header-read a few per (donor x project) to REPORT assay
# type, read length, assembly and mapped depth.
#
#   bash cluster/probe_donors.sh PD36713 PD36714 ...           # ids on the command line
#   bash cluster/probe_donors.sh -f cluster/liver.donors.txt   # one id per line
#   PROBE=3 bash cluster/probe_donors.sh PD36713               # more spot checks per group
#
# WHY NOT probe_new_organs.sh: that one is MANIFEST-driven (donor<TAB>sample) and is the
# right tool when the paper's WGS sample list is known, because a manifest is the only exact
# way to exclude a same-named targeted panel. Use THIS script when there is NO sample list —
# it globs the donor instead. The cost: if the cohort has a panel arm named like its WGS, the
# glob will sweep it in. That is fine here because we REPORT DS and read it off; it is NOT a
# selector. Never build a discovery fofn straight from this — confirm WGS per sample first.
#
# WHY NOT survey_assay.sh: it hard-REJECTs anything that is not hs37d5 (correct for the old
# GRCh37 colony cohorts, WRONG for GRCh38 releases) and it does not measure read length. This
# script reports assembly instead of judging on it, and it reads ~1000 records/BAM to get the
# read length — the one field that can silently kill clip-based MEI discovery (75bp vs 151bp)
# and that a header alone cannot give you.
#
# READ-ONLY, head-node safe: only `samtools view -H`, `idxstats`, and `head -1000` records.
# Never point a real (record-reading) job at nst_links — stage to Lustre first.
set -uo pipefail

DONORS=()
if [ "${1:-}" = "-f" ]; then
    [ -s "${2:-}" ] || { echo "no such donor file: ${2:-}" >&2; exit 1; }
    mapfile -t DONORS < <(grep -oE 'PD[0-9]+' "$2" | sort -u)
else
    [ "$#" -ge 1 ] || { echo "usage: $0 PD.. PD.. | -f donors.txt" >&2; exit 1; }
    DONORS=("$@")
fi

NST=${NST:-/nfs/cancer_ref01/nst_links/live}
PROBE=${PROBE:-2}

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

echo "donors: ${#DONORS[@]}    nst_links: $NST    PROBE=$PROBE per (donor x project)"
echo

# ---- resolve each donor -> project<TAB>sample (dedup on sample basename per project) ------
: > "$TMP/found.tsv"
for d in "${DONORS[@]}"; do
    ls -d "$NST"/*/"${d}"* 2>/dev/null \
        | awk -F/ -v d="$d" '{print d"\t"$(NF-1)"\t"$NF}' | sort -u >> "$TMP/found.tsv"
done

echo "### distinct samples per (donor x project) — a SUMMARY, not a selector"
if [ -s "$TMP/found.tsv" ]; then
    awk -F'\t' '{print $1"\t"$2}' "$TMP/found.tsv" | sort | uniq -c \
        | awk '{printf "  %-9s proj %-6s %s samples\n", $2, $3, $1}'
else
    echo "  (none of the ${#DONORS[@]} donors resolved in nst_links — check iRODS, see irods-cgp-zone-limits)"
fi

echo
echo "### header probe: $PROBE BAM(s) per (donor x project)"
printf "%-9s %-6s %-24s %-8s %-7s %-14s %s\n" DONOR PROJ DS READLEN MAPPED ASSEMBLY SAMPLE
echo "--------------------------------------------------------------------------------------"
awk -F'\t' '{print $1"\t"$2}' "$TMP/found.tsv" | sort -u | while IFS=$'\t' read -r d p; do
    n=0
    awk -F'\t' -v d="$d" -v p="$p" '$1==d && $2==p{print $3}' "$TMP/found.tsv" | while read -r s; do
        [ "$n" -ge "$PROBE" ] && break
        b="$NST/$p/$s/$s.sample.dupmarked.bam"
        [ -s "$b" ] || continue
        hdr=$(samtools view -H "$b" 2>/dev/null)
        [ -z "$hdr" ] && continue
        ds=$(printf '%s\n' "$hdr" | grep '^@RG' | tr '\t' '\n' | grep '^DS:' | sed 's/^DS://' | sort -u | paste -sd, -)
        if printf '%s\n' "$hdr" | grep -qm1 '^@SQ.*SN:hs37d5'; then asm=hs37d5_GRCh37
        else
            c1=$(printf '%s\n' "$hdr" | grep -m1 '^@SQ.*SN:chr1[[:space:]]' | tr '\t' '\n' | grep -m1 '^LN:' | sed 's/^LN://')
            o1=$(printf '%s\n' "$hdr" | grep -m1 '^@SQ.*SN:1[[:space:]]'    | tr '\t' '\n' | grep -m1 '^LN:' | sed 's/^LN://')
            case "${c1:-$o1}" in 248956422) asm=GRCh38 ;; 249250621) asm=hg19_GRCh37 ;; *) asm="other(${c1:-${o1:-?}})" ;; esac
        fi
        rl=$(samtools view "$b" 2>/dev/null | head -1000 \
             | awk '{print length($10)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
        m=$(samtools idxstats "$b" 2>/dev/null | awk '{r+=$3} END {printf "%.0f", (r+0)/1000000}')
        printf "%-9s %-6s %-24s %-8s %-7s %-14s %s\n" "$d" "$p" "${ds:-NONE}" "${rl:-?}bp" "${m:-0}M" "$asm" "$s"
        n=$((n+1))
    done
done

cat <<'EOF'

### READING THIS — INVENTORY, not a cohort selection
  DS must say WGS. TARGETED_ILLUMINA / panel => not usable for insertion discovery.
  READLEN 151bp good; 75bp badly compromises clip-based MEI discovery (remap cannot fix it).
  ASSEMBLY is triage only (we remap to hs1) but re-check the @SQ of every BAM you actually use.
  If DS varies within a donor, the cohort has a panel arm mixed in by name — resolve per sample.
  Record the answer in cluster/trees/INVENTORY.md.
EOF
