#!/bin/bash
# For a manifest of TREE TIPS (colony sample names), check nst_links for a usable WGS BAM per
# tip. Answers: "do we have every BAM for every tip, and is it the WGS one (not a targeted
# twin)?" — the same colony basename can have BOTH a WGS and a targeted release under different
# projects (Chapman: PD45808b_lo0003 is WGS in proj 2569 AND targeted in proj 2781), so a tip
# is only truly covered if it has a WGS BAM, not merely *a* BAM. See chapman2024-hsct-grch37-only.
#
#   bash cluster/check_tree_bams.sh cluster/manifests/chapman_hsct.tips.tsv
#   WGS_PROJ="2189 2256 ..." bash cluster/check_tree_bams.sh <manifest>   # override known WGS projects
#
# Manifest: any TSV; the colony sample is grepped as PD<num><letters>_lo<num> from each line.
# Cheap first: classify each tip's projects against the known WGS-project list (no header read).
# Definitive fallback: for any tip NOT found in a known-WGS project, header-read its BAM(s) to
# read @RG DS: — so a WGS release in a project we did not pre-list is still counted, and a
# targeted-only tip is correctly flagged. Header reads only (idxstats/-H) — head-node safe.
set -uo pipefail

MAN="${1:?usage: $0 <manifest.tsv>}"
[ -s "$MAN" ] || { echo "no such manifest: $MAN" >&2; exit 1; }
NST=${NST:-/nfs/cancer_ref01/nst_links/live}
# Chapman 2024 HSCT WGS projects, from ~/probe_chapman_hsct.txt. Override for other cohorts.
WGS_PROJ="${WGS_PROJ:-2189 2256 2448 2496 2569 2679 2693}"

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

is_known_wgs() { local p="$1" w; for w in $WGS_PROJ; do [ "$p" = "$w" ] && return 0; done; return 1; }

ds_of() { samtools view -H "$1" 2>/dev/null | grep '^@RG' | tr '\t' '\n' | grep '^DS:' | sed 's/^DS://' | sort -u | paste -sd, -; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mapfile -t SAMPLES < <(grep -oE 'PD[0-9]+[a-z]+_lo[0-9]+' "$MAN" | sort -u)
echo "manifest: $MAN    distinct tips: ${#SAMPLES[@]}    known WGS projects: $WGS_PROJ"
echo "checking nst_links ($NST) — this walks every tip, be patient..."
echo

: > "$TMP/rows.tsv"   # sample  status  detail
for s in "${SAMPLES[@]}"; do
    d=$(printf '%s' "$s" | grep -oE '^PD[0-9]+')
    mapfile -t projs < <(ls -d "$NST"/*/"$s" 2>/dev/null | awk -F/ '{print $(NF-1)}' | sort -u)
    if [ "${#projs[@]}" -eq 0 ]; then
        printf '%s\t%s\t%s\t%s\n' "$d" "$s" "MISSING" "no dir in nst_links" >> "$TMP/rows.tsv"; continue
    fi
    # cheap path: any known WGS project?
    hit=""; for p in "${projs[@]}"; do is_known_wgs "$p" && { hit="$p"; break; }; done
    if [ -n "$hit" ]; then
        printf '%s\t%s\t%s\t%s\n' "$d" "$s" "WGS" "proj $hit" >> "$TMP/rows.tsv"; continue
    fi
    # fallback: header-read to see if any project is WGS after all
    wgsp=""; anyds=""
    for p in "${projs[@]}"; do
        b="$NST/$p/$s/$s.sample.dupmarked.bam"; [ -s "$b" ] || continue
        ds=$(ds_of "$b"); anyds="${anyds:+$anyds; }$p:$ds"
        case "$ds" in WGS*) wgsp="$p"; break;; esac
    done
    if [ -n "$wgsp" ]; then
        printf '%s\t%s\t%s\t%s\n' "$d" "$s" "WGS" "proj $wgsp (unlisted; header-confirmed)" >> "$TMP/rows.tsv"
    elif [ -n "$anyds" ]; then
        printf '%s\t%s\t%s\t%s\n' "$d" "$s" "NON_WGS_ONLY" "$anyds" >> "$TMP/rows.tsv"
    else
        printf '%s\t%s\t%s\t%s\n' "$d" "$s" "NO_BAM" "dirs exist but no .sample.dupmarked.bam: ${projs[*]}" >> "$TMP/rows.tsv"
    fi
done

echo "### per-donor coverage (tips with a WGS BAM / total tips)"
awk -F'\t' '{tot[$1]++; if($3=="WGS")ok[$1]++} END{
  for(d in tot) printf "  %-9s %4d / %-4d WGS%s\n", d, ok[d]+0, tot[d], (ok[d]==tot[d]?"":"   <-- INCOMPLETE")
}' "$TMP/rows.tsv" | sort

echo
tot=$(wc -l < "$TMP/rows.tsv" | tr -d ' ')
okc=$(awk -F'\t' '$3=="WGS"' "$TMP/rows.tsv" | wc -l | tr -d ' ')
echo "### TOTAL: $okc / $tot tips have a WGS BAM in nst_links"
echo
echo "### tips WITHOUT a WGS BAM (empty = every tip is covered)"
awk -F'\t' '$3!="WGS"{printf "  %-9s %-20s %-14s %s\n",$1,$2,$3,$4}' "$TMP/rows.tsv" | sort | head -200
miss=$(awk -F'\t' '$3!="WGS"' "$TMP/rows.tsv" | wc -l | tr -d ' ')
[ "$miss" -gt 200 ] && echo "  ...($miss total; showing 200)"
echo
echo "full per-tip table: writing to \$HOME/$(basename "$MAN" .tsv).bam_coverage.tsv"
{ echo -e "donor\tsample\tstatus\tdetail"; sort "$TMP/rows.tsv"; } > "$HOME/$(basename "$MAN" .tsv).bam_coverage.tsv"

cat <<'EOF'

### READING THIS
  WGS          = a WGS BAM exists for this tip (usable).
  MISSING      = no nst_links dir for this colony at all — may still be in iRODS only
                 (nst_links is an incomplete view; see irods-cgp-zone-limits). Confirm with iquest.
  NON_WGS_ONLY = colony exists but ONLY as a targeted/RNA release (e.g. the proj 2781 twin) —
                 no WGS BAM. This is the case dedup-on-basename can silently mask.
  NO_BAM       = dir exists but holds no .sample.dupmarked.bam (analysis-only release).
EOF
