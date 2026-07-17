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
# Phase 2 (catalogue_headers.sh, an LSF array) reads the headers for the rows found here.
#
# Env: NST (default /nfs/cancer_ref01/nst_links/live), OUT (default ~/catalogue),
#      DONORS (default cluster/donors.txt)
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
mapfile -t WANT < <(sed 's/#.*//' "$DONORS" | awk 'NF{print $1}')
[ "${#WANT[@]}" -gt 0 ] || { echo "$DONORS has no donors" >&2; exit 1; }
echo "scope: ${#WANT[@]} donors from $DONORS${SAMPLES:+ + explicit samples from $SAMPLES}"

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
echo
echo "next: bash cluster/submit_catalogue.sh   # header reads, as a throttled LSF array"
