#!/bin/bash
# PHASE 2 of the sample catalogue: read BAM headers for one CHUNK of the manifest and emit
# one TSV row per sample. Run as an LSF array task (see submit_catalogue.sh).
#
#   bash cluster/catalogue_headers.sh <chunk-index>
#
# Env: OUT (default ~/catalogue), CHUNK (rows per task, default 500)
#
# EVERYTHING the header can tell us, because we keep getting burned by not knowing what the
# data IS. Each field below exists because a specific assumption turned out to be wrong:
#
#   assay     (@RG DS)   TARGETED vs WGS. Cost us a full 50-colony discovery run: three
#                        "patients" were gene panels (11, 7, 92 loci vs 8,833 for real WGS).
#                        Varies PER BAM: PD43947 has TARGETED and WGS releases side by side.
#   assembly  (@SQ SN)   hs37d5 / chr-style / bare GRCh37. Varies WITHIN a project (2149
#                        holds both). combine_insertions is assembly-specific and a mismatch
#                        yields plausible, silently WRONG coordinates.
#   asm_name  (@SQ AS)   the header's own name for the build. Reads "NCBI38", not "GRCh38".
#   ref_len   (sum LN)   total reference length = a build fingerprint that cannot be faked by
#                        a mislabelled AS: tag. hs37d5 and GRCh38 differ here.
#   chr1_md5  (@SQ M5)   md5 of the chr1/1 sequence: the DEFINITIVE build identity. Two files
#                        claiming "GRCh38" with different M5 are different references.
#   read_len             75bp vs 151bp materially changes clip-based MEI discovery.
#                        Remapping CANNOT fix this — check before planning a remap.
#   mapped/unmapped      WGS ~166-510M vs targeted 3-13M. Also catches dud libraries:
#                        PD34200c_lo0001 has 72k reads total (a failed sample, not an index bug).
#   median_insert        library fragment size; drives discordant-pair anchoring.
#   n_libraries/n_runs   >1 means the sample is a merge; batch/technology cohorts.
#   run_dates (@RG DT)   2017 HiSeq vs 2019 NovaSeq -> different error and clip profiles.
#   platform/model       ditto.
#   dup_tool  (@PG)      bamsormadup vs bammarkduplicates2 vs biobambam version.
#   bwa + ref (@PG)      the exact aligner version and reference FASTA path used.
#
# CRAM: idxstats and headers are reference-free, but RECORDS are not. Without REF_PATH,
# samtools fetches reference slices from EBI over the internet — 30 array tasks doing that is
# slow and antisocial. We emit CRAM_NOREF rather than hang. Set REF_PATH (module load
# sanger-samtools-refpath) to get read_len/insert for CRAMs.
set -uo pipefail

IDX="${1:?usage: catalogue_headers.sh <1-based chunk index>}"
OUT=${OUT:-$HOME/catalogue}
CHUNK=${CHUNK:-500}
MAN="$OUT/manifest.tsv"
mkdir -p "$OUT/parts"
DEST="$OUT/parts/part.$IDX.tsv"

[ -s "$MAN" ] || { echo "no manifest: $MAN (run catalogue_scan.sh)" >&2; exit 1; }
if [ -s "$DEST" ]; then echo "[$IDX] part exists, skipping"; exit 0; fi

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

NCOL=27
blank() { printf '%s\t' "$@"; for i in $(seq $(( NCOL - $# )) ); do printf -- '-\t'; done; printf '\n'; }

START=$(( (IDX-1)*CHUNK + 2 ))
END=$(( IDX*CHUNK + 1 ))
TMP="$DEST.tmp.$$"
: > "$TMP"
n=0
while IFS=$'\t' read -r proj sample donor bam bytes bai stale; do
    [ -n "${proj:-}" ] || continue
    n=$((n+1))

    if [ "$bam" = "-" ] || [ ! -f "$bam" ]; then
        blank "$proj" "$sample" "$donor" "${bam:--}" "${bytes:-0}" "${stale:--}" "NO_BAM" >> "$TMP"; continue
    fi
    hdr=$(samtools view -H "$bam" 2>/dev/null)
    if [ -z "$hdr" ]; then
        blank "$proj" "$sample" "$donor" "$bam" "$bytes" "${stale:--}" "UNREADABLE" >> "$TMP"; continue
    fi

    rg=$(printf '%s\n' "$hdr" | grep '^@RG')
    sq=$(printf '%s\n' "$hdr" | grep '^@SQ')
    pg=$(printf '%s\n' "$hdr" | grep '^@PG')
    # unique values of an @RG tag, comma-joined
    rgf() { printf '%s\n' "$rg" | tr '\t' '\n' | grep "^$1:" | sed "s/^$1://" | sort -u | paste -sd, - ; }
    rgn() { printf '%s\n' "$rg" | tr '\t' '\n' | grep -c "^$1:" ; }
    rgu() { printf '%s\n' "$rg" | tr '\t' '\n' | grep "^$1:" | sort -u | wc -l | tr -d ' ' ; }

    ds=$(rgf DS); pl=$(rgf PL); pm=$(rgf PM); cn=$(rgf CN); sm=$(rgf SM)
    dt=$(printf '%s\n' "$rg" | tr '\t' '\n' | grep '^DT:' | sed 's/^DT://; s/T.*//' | sort -u | paste -sd, -)
    nlb=$(rgu LB); npu=$(rgu PU); nrg=$(printf '%s\n' "$rg" | grep -c '^@RG')
    pu1=$(printf '%s\n' "$rg" | tr '\t' '\n' | grep -m1 '^PU:' | sed 's/^PU://')

    # Match against the variable, NOT through a pipe into `grep -q`. Under `set -o pipefail`,
    # grep -q exits on first match and SIGPIPEs the upstream printf, so the pipeline reports
    # FAILURE on a successful match. It bit exactly the branch that mattered: chr1 is the
    # FIRST @SQ line (grep quits early -> printf killed -> "no match" -> fell through to the
    # fallback and labelled 9,925 GRCh38 BAMs "chr1"), while hs37d5 is the LAST contig
    # (printf finishes first -> no SIGPIPE -> worked). A silent miss on one build only.
    if   [[ "$sq" == *"SN:hs37d5"* ]]; then asm=hs37d5
    elif [[ "$sq" == *$'SN:chr1\t'* ]]; then asm=chr-style   # \t so chr10/chr11 don't match
    else asm=$(printf '%s\n' "$sq" | head -1 | tr '\t' '\n' | grep -m1 '^SN:' | sed 's/^SN://'); fi
    nsq=$(printf '%s\n' "$sq" | grep -c '^@SQ')
    asname=$(printf '%s\n' "$sq" | tr '\t' '\n' | grep -m1 '^AS:' | sed 's/^AS://')
    reflen=$(printf '%s\n' "$sq" | tr '\t' '\n' | grep '^LN:' | sed 's/^LN://' | awk '{s+=$1} END{print s+0}')
    # M5 of the FIRST sequence (chr1 or 1): the definitive reference fingerprint
    md5=$(printf '%s\n' "$sq" | head -1 | tr '\t' '\n' | grep -m1 '^M5:' | sed 's/^M5://')
    so=$(printf '%s\n' "$hdr" | grep -m1 '^@HD' | tr '\t' '\n' | grep -m1 '^SO:' | sed 's/^SO://')

    bwacl=$(printf '%s\n' "$pg" | grep -m1 'PN:bwa.*bwa mem')
    bwav=$(printf '%s\n' "$bwacl" | tr '\t' '\n' | grep -m1 '^VN:' | sed 's/^VN://')
    ref=$(printf '%s\n' "$bwacl" | grep -oE '/[^[:space:]]*\.(fa|fasta)' | head -1)
    duptool=$(printf '%s\n' "$pg" | grep -oE 'PN:(bamsormadup|bammarkduplicates2|bamstreamingmarkduplicates)' | sed 's/^PN://' | sort -u | paste -sd, -)
    npg=$(printf '%s\n' "$pg" | grep -c '^@PG')

    mapped=$(samtools idxstats "$bam" 2>/dev/null | awk '{m+=$3; u+=$4} END {print (m+0)"\t"(u+0)}')
    nmap=$(printf '%s' "$mapped" | cut -f1); nunmap=$(printf '%s' "$mapped" | cut -f2)

    case "$bam" in
        *.cram) if [ -z "${REF_PATH:-}" ] && [ -z "${REF_CACHE:-}" ]; then rl=CRAM_NOREF; ins=CRAM_NOREF; else rl=""; fi;;
        *) rl="";;
    esac
    if [ -z "${rl:-}" ]; then
        rl=$(samtools view "$bam" 2>/dev/null | head -1000 | awk '{print length($10)}' \
             | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
        ins=$(samtools view -f 2 -F 0x900 "$bam" 2>/dev/null | head -2000 \
              | awk '$9>0 && $9<2000 {print $9}' | sort -n \
              | awk '{a[NR]=$1} END {if(NR)print a[int(NR/2)+1]; else print "-"}')
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$proj" "$sample" "$donor" "$bam" "$bytes" "${stale:--}" \
        "${ds:-NONE}" "${asm:--}" "${asname:--}" "${reflen:-0}" "${md5:--}" "${nsq:-0}" "${so:--}" \
        "${rl:--}" "${nmap:-0}" "${nunmap:-0}" "${ins:--}" \
        "${pl:--}" "${pm:--}" "${cn:--}" "${dt:--}" "${nrg:-0}" "${nlb:-0}" "${npu:-0}" "${pu1:--}" \
        "${duptool:--}" "${bwav:--}" >> "$TMP"
done < <(sed -n "${START},${END}p" "$MAN")

mv -f "$TMP" "$DEST"
echo "[$IDX] $n rows -> $DEST"
