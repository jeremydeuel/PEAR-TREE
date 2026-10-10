#!/bin/bash
#BSUB -n 12
#BSUB -M 4000
#BSUB -R "select[mem>4000] rusage[mem=4000] span[hosts=1]"
#BSUB -q normal
#BSUB -J populate_tsv
#BSUB -o logs/populate_tsv.%J.log
#BSUB -e logs/populate_tsv.%J.err
# ^ LSF opens these logs at DISPATCH, before the script runs — so `logs/` MUST already exist
#   in the submission CWD or the job never starts. Submit with the mkdir baked in (see below):
#       mkdir -p logs && bsub < cluster/populate_colonies_tsv.sh
#
# populate_colonies_tsv.sh — fill the placeholder patients/<org>/<pat>/colonies.tsv
# files IN PLACE, so you can `git add patients && git commit && git push`.
#
# What it writes per tip:  donor  proj  ds  readlen  mapped  assembly  sample
#
# DATA SOURCE = nst_links, HEADER/INDEX reads only (samtools view -H, idxstats) — the
# Sanger-policy-safe access (see hpc-sanger skill: "always stage BAMs; never read records
# from nst_links; the ONLY permitted direct access is reading the HEADER"). We never stage
# and never stream records for the mandatory columns. The one exception is `readlen`, gated
# OFF by default — see READLEN_MODE below.
#
# ── PREREQUISITES (run order) ─────────────────────────────────────────────────
#   1. LOCALLY: git add patients/ (tree files + placeholder colonies.tsv) ; commit ; push.
#      The farm checkout needs the .tree files to know each patient's tips.
#   2. ON FARM, from the repo root of the checkout on Lustre:
#         cd /lustre/.../PEAR-TREE-pd44579 && git pull
#         mkdir -p logs && bsub < cluster/populate_colonies_tsv.sh    # submit (12h normal queue)
#      or run directly on a test/interactive node:  bash cluster/populate_colonies_tsv.sh
#   3. LOCALLY (or on farm): git add patients/*/*/colonies.tsv ; commit ; push.
#
# ── WHAT IT DOES NOT DO — patients with no PD-style tips are SKIPPED (listed at the end) ──
#   • Liver (Brunner, 34 patients): tips are anonymised clone ids (Cl.27), NOT sequencing
#     sample names — no tree→BAM path without the paper's Cl.NN→PD map.
#   • Lab-codename cohorts (e.g. AX001 / PD43976: tips like BMH1_TG001_P31_A11): the tip is a
#     lab codename, not an nst_links sample id — needs a codename→sample map to populate.
#   These are reported under "patients SKIPPED" so nothing is silently left unfilled. Populate
#   them separately once the relevant id map is in hand.
#   • Already-populated TSVs (the 10 Chapman pairs, or a rerun) are left untouched.
#
# ── CONFIG (env overrides) ────────────────────────────────────────────────────
#   NST         nst_links root      (default /nfs/cancer_ref01/nst_links/live)
#   PAR         parallel probes     (default 12; keep modest — these are iRODS-backed NFS)
#   READLEN_MODE  header|record     (default header) — see below
#   ONLY        space-separated organ groups to restrict to (default: all mappable)
#   PREFER_ASSEMBLY  GRCh38|hs37d5_GRCh37|...: a colony with WGS releases on several assemblies takes
#               this one (default: the first WGS release found)
set -uo pipefail
export LC_ALL=C

# ── read length: the one column not derivable from header/index ───────────────
# There is no reliable read-length field in a BAM header. Two ways to fill `readlen`:
#   header  (DEFAULT): DO NOT touch records. Emit readlen=NA. Nothing is guessed, nothing is
#           streamed off nst_links — fully policy-clean. (Read length is uniform 151bp for
#           every WGS cohort probed so far, but this script refuses to hard-code a guess into
#           a per-tip column; fill it from the cohort probe if you want it.)
#   record  (opt-in, READLEN_MODE=record): read the FIRST alignment record of each selected
#           BAM (samtools view | head -1 → length($10)). This is a single BGZF block per BAM,
#           not a stage and not a scan — but it IS a record read off nst_links, so it is
#           OFF by default and you are opting into it knowingly. Use only if you accept that.
READLEN_MODE="${READLEN_MODE:-header}"
PREFER_ASSEMBLY="${PREFER_ASSEMBLY:-}"

command -v git >/dev/null || { echo "git not found" >&2; exit 1; }
REPO="${REPO:-$(git rev-parse --show-toplevel 2>/dev/null)}"
[ -n "$REPO" ] && [ -d "$REPO/patients" ] || { echo "run from the repo root; $REPO/patients missing" >&2; exit 1; }
PAT="$REPO/patients"
NST="${NST:-/nfs/cancer_ref01/nst_links/live}"
[ -d "$NST" ] || { echo "nst_links not visible: $NST" >&2; exit 1; }
PAR="${PAR:-12}"

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH (module load failed)" >&2; exit 1; }
command -v python3  >/dev/null || { echo "python3 not on PATH" >&2; exit 1; }

# Node-local scratch for the many tiny per-tip files — NEVER on Lustre (bad at small I/O,
# and parallel workers each write their OWN file so no two jobs touch one file).
WORK="$(mktemp -d "${TMPDIR:-/tmp}/poptsv.XXXXXX")" || exit 1
trap 'rm -rf "$WORK"' EXIT
echo "repo=$REPO  nst=$NST  par=$PAR  readlen_mode=$READLEN_MODE  work=$WORK"

mkdir -p "$WORK/idx"
# Per-DONOR nst_links index, built ONCE (see index_donors below, after the job list is known).
# One `ls` glob per donor (~100 total) instead of stat-ing every project for every tip — keeps
# NFS metadata load on the shared archive to a minimum. probe_sample then does zero filesystem
# lookups: it greps the donor's cached "sample<TAB>project" table.
index_donor() {                 # arg: donor (PD<digits>) -> $WORK/idx/<donor>.tsv of  sample\tproject
    local d="$1"
    # sample = donor digits followed by a lowercase letter, so ${d}[a-z]* avoids PD1234/PD12345
    # prefix collisions; -d only real dirs; awk pulls (sample, project) from .../project/sample.
    ls -d "$NST"/*/"$d"[a-z]* 2>/dev/null | awk -F/ '{print $NF"\t"$(NF-1)}' | sort -u > "$WORK/idx/$d.tsv"
}
export -f index_donor

# ── tip extraction: robust to trees WITH and WITHOUT branch lengths ────────────
# A Newick leaf label sits right after '(' or ',' and starts with a letter (this excludes
# bootstrap/support numbers, which start with a digit and follow ')'). Works whether or not a
# ':branchlength' follows. Emits one raw tip label per line.
extract_tips() { python3 - "$1" <<'PY'
import re,sys
s=open(sys.argv[1]).read()
for m in re.findall(r'[(,]\s*([A-Za-z][A-Za-z0-9._\-]*)', s):
    print(m)
PY
}

# candidate nst_links sample ids for a raw tip label (most-specific first):
#   strip a trailing '_hum' analysis tag; strip a trailing '_<n>' dup marker; keep base.
# Only PD-style ids are sequencing samples; anything else (Cl.27) yields nothing -> UNMAPPABLE.
candidates() {
    local t="$1" b
    t="${t%_hum}"
    case "$t" in
        PD[0-9]*) : ;;
        *) return 0 ;;                         # not a sample id (e.g. Cl.27)
    esac
    printf '%s\n' "$t"                          # exact tip (minus _hum)
    b="${t%_[0-9]}"; [ "$b" != "$t" ] && printf '%s\n' "$b"   # drop trailing _<n> dup marker
    # base PD<num><letters>[_lo<num>]  (in case tip carried extra suffixes)
    b="$(printf '%s' "$t" | grep -oE '^PD[0-9]+[a-z]+(_lo[0-9]+)?')"
    [ -n "$b" ] && [ "$b" != "$t" ] && printf '%s\n' "$b"
}

# read length is uniform 151 across probed WGS cohorts; assembly detection from @SQ.
assembly_from_header() {   # stdin: @SQ lines
    awk '
      /SN:hs37d5/        {d=1}
      /SN:chr1\t/ && /LN:248956422/ {g38=1}
      /SN:1\t/    && /LN:249250621/ {g37=1}
      /SN:chr1\t/ && /LN:249250621/ {hg19=1}
      END{
        if(d)        print "hs37d5_GRCh37";
        else if(g38) print "GRCh38";
        else if(hg19)print "hg19_GRCh37";
        else if(g37) print "GRCh37";
        else         print "UNKNOWN";
      }'
}

# probe one sample: pick the WGS BAM by @RG DS across its projects; emit a TSV row to $WORK/<key>
probe_sample() {                # args: patient_key  raw_tip
    local key="$1" tip="$2"
    local out="$WORK/$key.rows/$tip.row"   # NB: own line — a var can't reference an earlier one in the SAME `local`
    local donor sample proj bam ds hdr asm mapped raw rl="NA" chosen="" chosends=""
    donor="$(printf '%s' "$tip" | grep -oE '^PD[0-9]+')"
    # resolve the sample by grepping the donor's cached index (built once, zero fs lookups here);
    # try each candidate id (exact, _hum-stripped, dup-marker-stripped) until one has projects.
    local c projs=() idx="$WORK/idx/$donor.tsv"
    while IFS= read -r c; do
        [ -z "$c" ] && continue
        mapfile -t projs < <(awk -F'\t' -v s="$c" '$1==s{print $2}' "$idx" 2>/dev/null)
        [ "${#projs[@]}" -gt 0 ] && { sample="$c"; break; }
    done < <(candidates "$tip")
    if [ -z "${sample:-}" ]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$donor" "-" "MISSING" "NA" "NA" "NA" "$tip" > "$out"; return
    fi
    # choose the WGS project (DS discriminates assay — NOT project id, name, or read length).
    # PREFER_ASSEMBLY (e.g. GRCh38): among the colony's WGS releases take the one on that assembly,
    # else the first WGS one -- PD45534 has colonies released on GRCh38 AND hs37d5 in different projects
    for proj in "${projs[@]}"; do
        bam="$NST/$proj/$sample/$sample.sample.dupmarked.bam"
        [ -s "$bam" ] || continue
        hdr="$(samtools view -H "$bam" 2>/dev/null)" || continue
        ds="$(printf '%s\n' "$hdr" | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^DS:/){sub(/^DS:/,"",$i); print $i}}' | sort -u | paste -sd, -)"
        case "$ds" in
            WGS*)
                case "$chosends" in WGS*) : ;; *) chosen="$proj"; chosends="$ds" ;; esac
                [ -z "${PREFER_ASSEMBLY:-}" ] && break
                if [ "$(printf '%s\n' "$hdr" | grep '^@SQ' | assembly_from_header)" = "$PREFER_ASSEMBLY" ]; then
                    chosen="$proj"; chosends="$ds"; break
                fi ;;
            *)    [ -z "$chosends" ] && { chosen="$proj"; chosends="$ds"; } ;;   # remember a fallback
        esac
    done
    # No WGS release for this colony (only a targeted/RNA twin, or no BAM at all). Do NOT leak
    # the non-WGS BAM's assembly/mapped into the WGS columns — flag the row and keep them NA.
    case "$chosends" in
        WGS*) : ;;
        "")   printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$donor" "NO_BAM"  "NO_BAM"          "NA" "NA" "NA" "$sample" > "$out"; return ;;
        *)    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$donor" "NO_WGS"  "NO_WGS($chosends)" "NA" "NA" "NA" "$sample" > "$out"; return ;;
    esac
    bam="$NST/$chosen/$sample/$sample.sample.dupmarked.bam"
    hdr="$(samtools view -H "$bam" 2>/dev/null)"
    asm="$(printf '%s\n' "$hdr" | grep '^@SQ' | assembly_from_header)"
    if raw="$(samtools idxstats "$bam" 2>/dev/null)"; then      # index read (header/index-safe)
        mapped="$(printf '%s\n' "$raw" | awk -F'\t' '{s+=$3} END{print s+0}')"
    else
        mapped="IDXERR"
    fi
    if [ "$READLEN_MODE" = record ]; then
        rl="$(samtools view "$bam" 2>/dev/null | head -1 | awk -F'\t' '{print length($10)}')"
        [ -z "$rl" ] && rl="NA"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$donor" "$chosen" "$chosends" "$rl" "$mapped" "$asm" "$sample" > "$out"
}
export -f probe_sample candidates assembly_from_header
export NST WORK READLEN_MODE PREFER_ASSEMBLY

# ── build the work list: (patient_key, tree, tip) for every unpopulated, mappable patient ──
ONLY="${ONLY:-}"
: > "$WORK/skipped"; : > "$WORK/jobs"; declare -A PDIR
mapfile -t tsvs < <(find "$PAT" -mindepth 3 -maxdepth 3 -name colonies.tsv | sort)
for tsv in "${tsvs[@]}"; do
    pdir="$(dirname "$tsv")"; org="$(basename "$(dirname "$pdir")")"; pat="$(basename "$pdir")"
    [ -n "$ONLY" ] && ! printf '%s ' $ONLY | grep -qw "$org" && continue
    # skip already-populated (has a data row that is neither comment nor header)
    if awk '!/^#/ && !/^donor\t/ && NF>0 {found=1; exit} END{exit !found}' "$tsv"; then continue; fi
    key="$org//$pat"; PDIR["$key"]="$pdir"; mkdir -p "$WORK/$key.rows"
    mappable=0
    while IFS= read -r tip; do
        [ -z "$tip" ] && continue
        if [ -n "$(candidates "$tip")" ]; then
            printf '%s\t%s\n' "$key" "$tip" >> "$WORK/jobs"; mappable=$((mappable+1))
        fi
    done < <(for t in "$pdir"/*.tree; do [ -e "$t" ] && extract_tips "$t"; done | sort -u)
    if [ "$mappable" -eq 0 ]; then
        # why: no tree at all, vs tips present but none PD-style (liver Cl.NN / lab codenames)
        extip="$(for t in "$pdir"/*.tree; do [ -e "$t" ] && extract_tips "$t"; done | head -1)"
        [ -z "$extip" ] && reason="no .tree on disk" || reason="no PD-style tips (e.g. '$extip')"
        printf '%s\t%s\n' "$key" "$reason" >> "$WORK/skipped"
        unset 'PDIR[$key]'
    fi
done
njob=$(wc -l < "$WORK/jobs" | tr -d ' ')
echo "probing $njob tips across ${#PDIR[@]} patients (skipping $(wc -l < "$WORK/skipped" | tr -d ' ') patients with no mappable tips)"
[ "$njob" -eq 0 ] && { echo "nothing to do"; exit 0; }

# ── build the per-donor nst_links index ONCE (one `ls` glob per donor) ─────────────────────
awk -F'\t' '{print $2}' "$WORK/jobs" | grep -oE '^PD[0-9]+' | sort -u > "$WORK/donors"
echo "indexing $(wc -l < "$WORK/donors" | tr -d ' ') donors in nst_links..."
xargs -P "$PAR" -L1 bash -c 'index_donor "$1"' _ < "$WORK/donors"

# ── probe in parallel; each worker writes its own file (Lustre-safe: node-local, 1 file/tip) ──
awk -F'\t' '{print $1"\t"$2}' "$WORK/jobs" | xargs -P "$PAR" -L1 bash -c 'probe_sample "$1" "$2"' _

# ── assemble per-patient colonies.tsv in the repo (single-threaded) ───────────────────────
HDR=$'donor\tproj\tds\treadlen\tmapped\tassembly\tsample'
CMT="# populated $(date +%F) by cluster/populate_colonies_tsv.sh (nst_links, header/index reads;"
CMT="$CMT readlen_mode=$READLEN_MODE). ds/assembly are per-BAM (WGS selected by @RG DS); mapped=idxstats mapped reads."
written=0
for key in "${!PDIR[@]}"; do
    pdir="${PDIR[$key]}"; rowdir="$WORK/$key.rows"
    ls "$rowdir"/*.row >/dev/null 2>&1 || continue
    # write-then-rename so a kill mid-write can't truncate a good placeholder in the git tree
    { echo "$HDR"; echo "$CMT"; cat "$rowdir"/*.row | sort; } > "$pdir/.colonies.tsv.tmp" \
        && mv -f "$pdir/.colonies.tsv.tmp" "$pdir/colonies.tsv"
    written=$((written+1))
done
echo "wrote $written colonies.tsv files"

echo; echo "### coverage summary (WGS tips / total tips, per patient)"
for key in "${!PDIR[@]}"; do
    rowdir="$WORK/$key.rows"; ls "$rowdir"/*.row >/dev/null 2>&1 || continue
    awk -F'\t' -v k="$key" '{tot++; if($3 ~ /^WGS/)wgs++} END{printf "  %-40s %4d / %-4d WGS%s\n", k, wgs+0, tot, (wgs==tot?"":"  <-- gaps")}' "$rowdir"/*.row
done | sort
echo
echo "### patients SKIPPED (no mappable tips — populate separately once an id map exists):"
sort "$WORK/skipped" | sed 's/^/  /'
