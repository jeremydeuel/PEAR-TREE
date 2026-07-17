#!/bin/bash
# PHASE 1 of the sample catalogue: enumerate the in-scope samples in nst_links and record
# BAM presence. NO header reads — this is stat() only, so it is safe to run on a head node.
#
#   bash cluster/catalogue_scan.sh                                   # donors from cluster/donors.txt
#   bash cluster/catalogue_scan.sh --samples ~/ax001.samples         # + non-PD-named donors
#
# SCOPE: only the donors in cluster/donors.txt — OUR cohort, not all of iRODS. Scanning
# every project would mean tens of thousands of header reads against shared NFS for samples
# we will never touch.
#
# --samples takes "DONOR PROJECT SAMPLE" lines, for donors whose sample names carry no
# donor (AX001's T1_A3, BMH1_plate_*). The donor column is REQUIRED and deliberately not
# inferred: the AX001 sample list spans projects 2315/2314/2360/2097/3226 AND 2445/2318
# "BMH1_*" rows, and whether AX001 and BMH1 are one person or two is not knowable from the
# names. Guessing wrong merges unrelated genomes into one clade — which manufactures
# exactly the cross-donor sharing this cohort exists to measure. Decide it, then write it.
#
# Phase 2 (catalogue_headers.sh, an LSF array) reads the headers for the rows found here —
# but only for the rows in todo.tsv, which this script computes. See "RESUME" below.
#
# Env: NST (default /nfs/cancer_ref01/nst_links/live), OUT (default ~/catalogue),
#      DONORS (default cluster/donors.txt), RESCAN=1 to force a full re-read
set -uo pipefail

cd "$(dirname "$0")/.."
NST=${NST:-/nfs/cancer_ref01/nst_links/live}
OUT=${OUT:-$HOME/catalogue}
DONORS=${DONORS:-cluster/donors.txt}
SAMPLES=""
[ "${1:-}" = "--samples" ] && { SAMPLES="${2:?--samples needs a file}"; }
mkdir -p "$OUT"
MAN="$OUT/manifest.tsv"

[ -d "$NST" ] || { echo "no nst_links at $NST" >&2; exit 1; }
[ -s "$DONORS" ] || { echo "no donor list: $DONORS" >&2; exit 1; }

# REFUSE TO RUN UNDER A LIVE ARRAY.
#
# Phase 2 slices its input BY ROW NUMBER (sed -n START,ENDp). Rewriting that input while
# tasks are still reading it re-points every task that has not started yet: task 7 keeps
# writing part.7, but part.7 now holds different samples. Nothing errors. The parts still
# merge. The catalogue is simply wrong, in a way no row count detects, because the row count
# is right. Wait for the array, or bkill it — the work is resumable now, so killing it costs
# only the tasks in flight.
if [ -s "$OUT/jobid" ] && command -v bjobs >/dev/null 2>&1; then
    JID=$(cat "$OUT/jobid")
    if bjobs "$JID" 2>/dev/null | grep -qE '\b(RUN|PEND)\b'; then
        echo "ABORT: catalogue array job $JID is still RUN/PEND." >&2
        echo "Phase 2 reads its input by row number; rewriting it now silently mis-assigns" >&2
        echo "headers to the wrong samples. Wait for it, or: bkill $JID" >&2
        exit 1
    fi
fi
mapfile -t WANT < <(sed 's/#.*//' "$DONORS" | awk 'NF{print $1}')
[ "${#WANT[@]}" -gt 0 ] || { echo "$DONORS has no donors" >&2; exit 1; }
echo "scope: ${#WANT[@]} donors from $DONORS${SAMPLES:+ + explicit samples from $SAMPLES}"

# CROSS-CHECK --samples AGAINST THE TREE TIPS, BEFORE SCANNING ANYTHING.
#
# A --samples file is the one input nothing else can validate: every row is well-formed, the
# sample dirs exist, the scan reports success — and the rows can still be the wrong donor
# entirely. AX001 was catalogued TWICE from a stale list (210 samples, "T1_A1", project 2315)
# when its real colonies are the 361 "BMH1_TG001_*" tips in project 2133. Both runs looked
# clean. 58 array tasks of header reads, twice, on samples we did not want.
#
# The tree tips are the independent witness: they say what a donor's colonies are CALLED, and
# they were derived from the phylogeny, not from a filename guess. Zero overlap between the
# file and the donor tips means the file is about something else -> refuse, do not scan.
#
# Partial overlap is NORMAL and only warns: a donor legitimately has bulk samples and extra
# colonies that never entered the tree. Absence from the tree is not absence from the donor.
TIPS=${TIPS:-cluster/trees/tips.tsv}
if [ -n "$SAMPLES" ]; then
    [ -s "$SAMPLES" ] || { echo "no such samples file: $SAMPLES" >&2; exit 1; }
    if [ ! -s "$TIPS" ]; then
        echo "WARNING: no tips at $TIPS — cannot cross-check $SAMPLES; scanning it unverified" >&2
    else
        echo "cross-check: $SAMPLES vs tree tips ($TIPS)"
        awk -v tips="$TIPS" '
        BEGIN{
            while ((getline line < tips) > 0) {
                n = split(line, a, /[ \t]+/); if (n < 2 || a[1] == "canon") continue
                tip[a[1] SUBSEP a[2]] = 1; hastips[a[1]] = 1
            }
        }
        /^[ \t]*(#|$)/ { next }
        {
            n = split($0, a, /[ \t]+/); if (n < 3) next
            d = a[1]; tot[d]++
            if ((d in hastips) && ((d SUBSEP a[3]) in tip)) ok[d]++
        }
        END{
            bad = 0
            for (d in tot) {
                if (!(d in hastips)) {
                    printf("  %-8s %4d samples — donor has no tree tips, cannot cross-check\n", d, tot[d])
                    continue
                }
                printf("  %-8s %4d samples, %4d match this donor tree tips\n", d, tot[d], ok[d]+0)
                if (ok[d]+0 == 0) { printf("  ^^ ZERO overlap: these are not %s colonies\n", d); bad = 1 }
            }
            exit bad
        }' "$SAMPLES" || {
            echo "ABORT: a --samples donor shares NO sample names with its tree tips." >&2
            echo "The file is about a different set than the donor it claims. Fix it before scanning." >&2
            exit 1
        }
    fi
fi

# stat is not portable: -c%s is GNU (farm), -f%z is BSD (mac). Getting this wrong returns
# an empty string that becomes 0 — a silent lie in a file we will later trust. Pick once.
#
# -L IS NOT OPTIONAL. nst_links is what its name says: every entry is a SYMLINK to the real
# BAM elsewhere. GNU stat does not dereference by default, so without -L this records the
# LENGTH OF THE TARGET PATH STRING (~99) instead of the file size — measured: 99 for
# PD34200c_lo0001, whose BAM is 3,258,431 bytes. That is a plausible-looking number in a
# column we would later use to spot truncated files, which is the worst kind of wrong.
if stat -Lc%s / >/dev/null 2>&1; then
    fsize() { stat -Lc%s "$1" 2>/dev/null || echo 0; }
else
    fsize() { stat -Lf%z "$1" 2>/dev/null || echo 0; }
fi

printf 'project\tsample\tdonor\tbam\tbam_bytes\tbai\tbai_stale\n' > "$MAN"

nsamp=0; nbam=0
emit() {  # emit <project> <sample> <donor>
    local proj="$1" s="$2" donor="$3" sdir="$NST/$1/$2" bam="" bytes bai stale ix
    [ -d "$sdir" ] || return 1
    nsamp=$((nsamp+1))
    for cand in "$sdir/$s.sample.dupmarked.bam" "$sdir/$s.sample.dupmarked.cram"; do
        [ -f "$cand" ] && { bam="$cand"; break; }
    done
    if [ -z "$bam" ]; then
        printf '%s\t%s\t%s\t-\t0\t-\t-\n' "$proj" "$s" "$donor" >> "$MAN"
        return 0
    fi
    nbam=$((nbam+1))
    bytes=$(fsize "$bam")
    bai="-"; stale="-"
    for ix in "$bam.bai" "$bam.crai" "${bam%.bam}.bai"; do
        [ -f "$ix" ] && { bai="$ix"; break; }
    done
    # Index-older-than-BAM check. An index that predates its BAM cannot describe it, and
    # idxstats reads the INDEX only — so a wrong index yields a confident wrong answer
    # (typically 0 mapped) that discovery would read as an empty genome, silently.
    # MEASURED: 1 of 14,525 samples trips this (2243/PD43947d_lo0010).
    # NOT the explanation for PD34200c_lo0001, which reports 0 mapped via idxstats while
    # `samtools view` still returns records, yet whose index is NEWER than its BAM. That
    # one is still unexplained — do not conflate the two.
    if [ "$bai" != "-" ]; then
        [ "$bai" -ot "$bam" ] && stale=YES || stale=no
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$proj" "$s" "$donor" "$bam" "$bytes" "$bai" "$stale" >> "$MAN"
}

# PD-named donors: find them wherever they appear. A donor is re-released under several
# project ids and the assembly differs between (and within) releases, so we record EVERY
# project/sample pair and let the query decide -- we do not pick a project here.
#
# The glob MUST NOT be a bare "$d"* : donor ids are not prefix-free. PD5163 is a real donor
# (71 tips, PD5163d_lo0074) and also a strict prefix of PD51632/3/4/5, the CML donors. A bare
# glob gave PD5163 449 colonies -- its own 71 plus ~378 stolen from four other donors -- and
# since those donors are ALSO in donors.txt, emit() ran twice per sample and the manifest
# carried the same BAM under two donor labels. Cross-donor contamination in the very table we
# use to tell donors apart.
# Sanger sample names are <DONOR><lowercase letter>_lo#### (PD5163d_lo0074, PD51632b_lo0001):
# the character after the donor id is NEVER a digit. That is the whole rule.
for projdir in "$NST"/*/; do
    proj=$(basename "$projdir")
    for d in "${WANT[@]}"; do
        for sdir in "$projdir$d"*/; do
            [ -d "$sdir" ] || continue
            s=$(basename "$sdir")
            case "$s" in "$d"[0-9]*) continue;; esac   # PD5163 must not eat PD51632b_lo0001
            emit "$proj" "$s" "$d"
        done
    done
done

# explicitly-mapped samples (non-PD naming): donor comes from the file, never inferred
if [ -n "$SAMPLES" ]; then
    [ -s "$SAMPLES" ] || { echo "no such samples file: $SAMPLES" >&2; exit 1; }
    while read -r donor proj s _; do
        case "${donor:-}" in ''|'#'*) continue;; esac
        [ -n "${s:-}" ] || { echo "  bad row (need 'DONOR PROJECT SAMPLE'): $donor $proj ${s:-}" >&2; continue; }
        emit "$proj" "$s" "$donor" || echo "  no dir: $NST/$proj/$s" >&2
    done < "$SAMPLES"
fi

echo "samples: $nsamp   with BAM/CRAM: $nbam"
echo "manifest: $MAN"

# ---------------------------------------------------------------- RESUME
#
# Phase 1 (this script) is stat() only and always re-runs in full — it is cheap, and it is
# the thing that discovers NEW samples. Phase 2 is the expensive half: ~4 NFS reads per BAM,
# tens of thousands of them. That is what we skip.
#
# THE DONE KEY IS THE ROW'S CONTENT, NOT ITS POSITION: project+sample+donor+bam+bytes.
#
# Keying on position (part.N exists -> skip) is what the old `-s "$DEST" && exit 0` did, and
# it is only safe while the manifest never changes. It changed: donors.txt went 97 -> 150
# donors, so every row index shifted and every existing part described samples it no longer
# aligned with. Content-keying is immune to that — and it also re-reads a sample whose BAM
# was re-released at a different size, or whose donor label we corrected, both of which a
# position key would happily keep stale forever.
#
# UNREADABLE rows are NOT done — a header read that failed is a transient to retry.
# NO_BAM rows ARE done: "this sample has no BAM" is a finding, and if a BAM later appears
# the bytes change, so the key changes, so it re-reads anyway.
DONE="$OUT/done.tsv"
TODO="$OUT/todo.tsv"

if [ "${RESCAN:-0}" = "1" ]; then
    echo
    echo "RESCAN=1: discarding $(( $(wc -l < "$DONE" 2>/dev/null || echo 1) - 1 )) done rows; everything re-reads"
    rm -f "$DONE"; rm -f "$OUT"/parts/*.tsv 2>/dev/null
fi

# Fold any completed parts into done.tsv, then clear them: parts are indexed against the
# todo file we are about to replace, so they cannot survive as parts. Their CONTENT survives.
mkdir -p "$OUT/parts"
nharv=$(ls "$OUT"/parts/*.tsv 2>/dev/null | wc -l | tr -d ' ')
if [ "$nharv" -gt 0 ]; then
    cat "$OUT"/parts/*.tsv >> "$DONE"
    rm -f "$OUT"/parts/*.tsv
    echo
    echo "harvested $nharv completed part(s) from the previous run into $DONE"
fi

# dedup done.tsv on the key, keeping the LAST occurrence (a later read supersedes an earlier)
if [ -s "$DONE" ]; then
    awk -F'\t' '{ k=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4 SUBSEP $5; if (!(k in seen)) { seen[k]=++n } row[seen[k]]=$0 }
         END { for (i=1; i<=n; i++) print row[i] }' "$DONE" > "$DONE.tmp" && mv -f "$DONE.tmp" "$DONE"
fi

awk -F'\t' -v done_f="$DONE" '
BEGIN {
    ndone = 0
    while ((getline line < done_f) > 0) {
        n = split(line, a, "\t"); if (n < 7) continue
        if (a[7] == "UNREADABLE") continue          # transient failure: retry it
        done[a[1] SUBSEP a[2] SUBSEP a[3] SUBSEP a[4] SUBSEP a[5]] = 1; ndone++
    }
}
NR == 1 { print; next }
{
    k = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4 SUBSEP $5
    if (k in done) { skipped++; next }
    print; todo++
}
END {
    printf("resume: %d already read, %d to read", skipped+0, todo+0) > "/dev/stderr"
    printf(" (%d done rows loaded)\n", ndone) > "/dev/stderr"
}' "$MAN" > "$TODO"

TODOROWS=$(( $(wc -l < "$TODO") - 1 ))
echo "todo:     $TODO  ($TODOROWS samples need header reads)"
echo
if [ "$TODOROWS" -le 0 ]; then
    echo "nothing to read — next: bash cluster/catalogue_merge.sh"
else
    echo "next: bash cluster/submit_catalogue.sh   # header reads, as a throttled LSF array"
fi
