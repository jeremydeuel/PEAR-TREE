#!/bin/bash
#BSUB -n 12
#BSUB -M 4000
#BSUB -R "select[mem>4000] rusage[mem=4000] span[hosts=1]"
#BSUB -q normal
#BSUB -J populate_donor
#BSUB -o logs/populate_donor.%J.log
#BSUB -e logs/populate_donor.%J.err
# ^ LSF opens these logs at DISPATCH — `logs/` MUST exist. Submit with:
#       mkdir -p logs && bsub < cluster/populate_donor_level.sh
#
# populate_donor_level.sh — DONOR-LEVEL colony inventory for patients whose tree tips are
# anonymised and therefore cannot be mapped tip->BAM (liver Cl.NN clone-ids; AX001/PD43976
# BMH plate-well codenames). For each (patient, donor) in the manifest it lists EVERY WGS
# colony of the donor found in nst_links — NOT keyed to the tree tips — and writes them to
# patients/<patient>/colonies.tsv:  donor proj ds readlen mapped assembly sample.
#
# This is the companion to cluster/populate_colonies_tsv.sh (which maps tree tips -> BAM for
# cohorts with real sample-id tips). Same policy: HEADER/INDEX reads only (samtools view -H,
# idxstats) — no records off nst_links; the WGS BAM is selected by @RG DS (not project/name/
# read length); assembly is read per-BAM from @SQ. readlen_mode as in the sibling script.
#
# WHY donor-level: the published re-analyses anonymised these tips (Chapman "Prolonged
# persistence" ships liver trees as Cl.NN with only donor-level SN_samples.txt; Mitchell's
# AX001/BMH1_TG works entirely in BMH codename space). No public or local clone/codename ->
# nst_links-sample map exists, so a per-tip inventory is impossible; the donor's actual
# sequenced WGS colonies are the honest, obtainable substitute. See patients/PLAN.md.
#
# NOTE: rows are ALL of the donor's WGS releases in nst_links — usually exactly the LCM/colony
# set, but a matched bulk-normal WGS (if the donor has one) would also appear. Rows are NOT
# 1:1 with the Cl.NN/BMH tree tips and the tree file is left unchanged.
#
#   MANIFEST (default cluster/donor_level.manifest.tsv): lines  <patient_relpath>\t<donor>
#   env: NST, PAR, READLEN_MODE (header|record) — same meaning as populate_colonies_tsv.sh.
set -uo pipefail
export LC_ALL=C
READLEN_MODE="${READLEN_MODE:-header}"

command -v git >/dev/null || { echo "git not found" >&2; exit 1; }
REPO="${REPO:-$(git rev-parse --show-toplevel 2>/dev/null)}"
[ -n "$REPO" ] && [ -d "$REPO/patients" ] || { echo "run from the repo root; $REPO/patients missing" >&2; exit 1; }
PAT="$REPO/patients"
NST="${NST:-/nfs/cancer_ref01/nst_links/live}"
[ -d "$NST" ] || { echo "nst_links not visible: $NST" >&2; exit 1; }
PAR="${PAR:-12}"
MAN="${MAN:-$REPO/cluster/donor_level.manifest.tsv}"
[ -s "$MAN" ] || { echo "manifest not found: $MAN" >&2; exit 1; }

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH (module load failed)" >&2; exit 1; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/popdonor.XXXXXX")" || exit 1
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/idx"
echo "repo=$REPO  nst=$NST  par=$PAR  readlen_mode=$READLEN_MODE  manifest=$MAN"

assembly_from_header() {   # stdin: @SQ lines
    awk '
      /SN:hs37d5/        {d=1}
      /SN:chr1\t/ && /LN:248956422/ {g38=1}
      /SN:1\t/    && /LN:249250621/ {g37=1}
      /SN:chr1\t/ && /LN:249250621/ {hg19=1}
      END{ if(d)print"hs37d5_GRCh37"; else if(g38)print"GRCh38"; else if(hg19)print"hg19_GRCh37"; else if(g37)print"GRCh37"; else print"UNKNOWN" }'
}
index_donor() {                 # arg: donor -> $WORK/idx/<donor>.tsv of  sample\tproject
    local d="$1"
    ls -d "$NST"/*/"$d"[a-z]* 2>/dev/null | awk -F/ '{print $NF"\t"$(NF-1)}' | sort -u > "$WORK/idx/$d.tsv"
}
# probe one (donor,sample): keep ONLY if a WGS release exists; write a row, else nothing.
probe_ds() {                    # args: donor sample
    local donor="$1" sample="$2" idx="$WORK/idx/$1.tsv" outdir="$WORK/rows/$1"
    local proj bam ds hdr asm mapped raw rl="NA" chosen="" chosends=""
    local projs=()
    mapfile -t projs < <(awk -F'\t' -v s="$sample" '$1==s{print $2}' "$idx" 2>/dev/null)
    for proj in "${projs[@]}"; do
        bam="$NST/$proj/$sample/$sample.sample.dupmarked.bam"; [ -s "$bam" ] || continue
        ds="$(samtools view -H "$bam" 2>/dev/null | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^DS:/){sub(/^DS:/,"",$i); print $i}}' | sort -u | paste -sd, -)"
        case "$ds" in WGS*) chosen="$proj"; chosends="$ds"; break ;; esac
    done
    [ -z "$chosen" ] && return 0     # no WGS release -> not a WGS colony, drop it
    bam="$NST/$chosen/$sample/$sample.sample.dupmarked.bam"
    hdr="$(samtools view -H "$bam" 2>/dev/null)"
    asm="$(printf '%s\n' "$hdr" | grep '^@SQ' | assembly_from_header)"
    if raw="$(samtools idxstats "$bam" 2>/dev/null)"; then
        mapped="$(printf '%s\n' "$raw" | awk -F'\t' '{s+=$3} END{print s+0}')"
    else mapped="IDXERR"; fi
    if [ "$READLEN_MODE" = record ]; then
        rl="$(samtools view "$bam" 2>/dev/null | head -1 | awk -F'\t' '{print length($10)}')"
        [ -z "$rl" ] && rl="NA"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$donor" "$chosen" "$chosends" "$rl" "$mapped" "$asm" "$sample" > "$outdir/$sample.row"
}
export -f assembly_from_header index_donor probe_ds
export NST WORK READLEN_MODE

# ── read manifest -> donor set + patient->donor map ───────────────────────────
declare -A P2D                       # patient_relpath -> donor
donors=()
while IFS=$'\t' read -r prel donor; do
    [ -z "${prel:-}" ] && continue
    case "$prel" in \#*) continue ;; esac
    [ -d "$PAT/$prel" ] || { echo "  WARN: no patient dir $prel (skipping)"; continue; }
    P2D["$prel"]="$donor"; donors+=("$donor")
done < "$MAN"
mapfile -t UDONORS < <(printf '%s\n' "${donors[@]}" | sort -u)
echo "manifest: ${#P2D[@]} patients, ${#UDONORS[@]} distinct donors"

# ── index each donor ONCE, then enumerate its distinct samples ────────────────
printf '%s\n' "${UDONORS[@]}" > "$WORK/donors"
echo "indexing ${#UDONORS[@]} donors in nst_links..."
xargs -P "$PAR" -L1 bash -c 'index_donor "$1"' _ < "$WORK/donors"
for d in "${UDONORS[@]}"; do
    mkdir -p "$WORK/rows/$d"
    awk -F'\t' -v d="$d" '{print d"\t"$1}' "$WORK/idx/$d.tsv" | sort -u
done | sort -u > "$WORK/jobs"
njob=$(wc -l < "$WORK/jobs" | tr -d ' ')
echo "probing $njob distinct (donor,sample) pairs..."
[ "$njob" -eq 0 ] && { echo "no samples found for any donor in nst_links"; exit 0; }

# ── probe each (donor,sample); WGS-only rows land in $WORK/rows/<donor>/ ───────
xargs -P "$PAR" -L1 bash -c 'probe_ds "$1" "$2"' _ < "$WORK/jobs"

# ── assemble each patient's colonies.tsv from its donor's WGS rows ─────────────
HDR=$'donor\tproj\tds\treadlen\tmapped\tassembly\tsample'
CMT="# donor-level inventory, populated $(date +%F) by cluster/populate_donor_level.sh"
CMT="$CMT (nst_links, header/index reads; readlen_mode=$READLEN_MODE). ROWS ARE ALL OF THE DONOR'S"
CMT="$CMT WGS COLONIES, not keyed to the anonymised tree tips (Cl.NN / BMH codenames); tree unchanged."
written=0; empty=""
for prel in "${!P2D[@]}"; do
    d="${P2D[$prel]}"; rowdir="$WORK/rows/$d"
    if ls "$rowdir"/*.row >/dev/null 2>&1; then
        { echo "$HDR"; echo "$CMT"; cat "$rowdir"/*.row | sort; } > "$PAT/$prel/.colonies.tsv.tmp" \
            && mv -f "$PAT/$prel/.colonies.tsv.tmp" "$PAT/$prel/colonies.tsv"
        written=$((written+1))
    else
        empty="${empty}\n  $prel (donor $d): NO WGS colonies found in nst_links"
    fi
done
echo "wrote $written colonies.tsv files"

echo; echo "### WGS colonies per donor"
for d in "${UDONORS[@]}"; do
    n=$(ls "$WORK/rows/$d"/*.row 2>/dev/null | wc -l | tr -d ' ')
    printf "  %-9s %4d WGS colonies\n" "$d" "$n"
done | sort
[ -n "$empty" ] && { echo; echo "### donors with NO WGS colonies in nst_links (check the id):"; printf '%b\n' "$empty"; }
