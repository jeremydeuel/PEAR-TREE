#!/bin/bash
# Dump EVERYTHING knowable about one sample: full BAM header, symlink target, and whatever
# metadata the surrounding systems expose. This is a DISCOVERY tool — run it on a couple of
# known-good samples and read the output, rather than trusting my guesses about where the
# Sanger keeps tissue and library-prep information.
#
#   bash cluster/sample_provenance.sh PD44890f_lo0001
#   bash cluster/sample_provenance.sh /nfs/cancer_ref01/nst_links/live/2460/X/X.sample.dupmarked.bam
#
# WHAT I ACTUALLY KNOW IS IN THE HEADER (measured, 14,525 BAMs):
#   @RG DS:   assay class — WGS_ILLUMINA_short / TARGETED_ILLUMINA_short / RNA-Seq / WXS.
#             This is assay, NOT prep chemistry: it will not tell you PCR-free vs PCR-plus.
#   @RG LB:   a library IDENTIFIER (e.g. "188036"). An id, not a type — useful as a batch
#             key, useless as prep description without a LIMS lookup.
#   @RG PU:   run_lane (e.g. "28162_1"); @RG DT: run date; @RG PM: instrument model.
#   @PG CL:   contains the run folder, e.g. .../analysis/190131_HX5_28162_A_HW72LCCXY/
#             = date_instrument_run_flowcell. Also the canpipe lane dir: 1748_329853.
#
# WHAT IS **NOT** AVAILABLE, ANYWHERE WE CAN REACH (settled 2026-07, do not re-litigate):
#   TISSUE / cell type and LIBRARY PREP chemistry.
#   - Not in the BAM header.
#   - Not in iRODS: farm22 sees ONLY the `cgp` zone (`ils /` -> just `C- /cgp`), and
#     `ils /seq` -> "does not exist or user lacks access permission". NPG's `seq` zone holds
#     the LIMS fields (library_type, sample_common_name, tissue) and is unreachable.
#   - cgp's own AVUs are five, none descriptive: file_extension, id_analysis_proc, id_ifile,
#     md5, version. Collection-level AVUs: none.
#   => Tissue/prep must come from the paper/supplement or EGA. The tissue column in
#      cluster/trees/donor_table.tsv is PAPER-DERIVED, never measured.
set -uo pipefail

ARG="${1:?usage: sample_provenance.sh <sample-name|bam-path>}"
NST=${NST:-/nfs/cancer_ref01/nst_links/live}
module load samtools-1.19/python-3.12.0 2>/dev/null || true

if [ -f "$ARG" ]; then BAM="$ARG"; else
    BAM=$(ls "$NST"/*/"$ARG"/"$ARG".sample.dupmarked.bam 2>/dev/null | head -1)
    [ -n "$BAM" ] || { echo "no BAM found for $ARG under $NST" >&2; exit 1; }
fi
echo "BAM: $BAM"

echo
echo "############ 1. the symlink chain (nst_links entries are symlinks)"
ls -l "$BAM"
tgt=$(readlink -f "$BAM" 2>/dev/null); echo "resolved: ${tgt:-<none>}"
ls -lL "$BAM" 2>/dev/null

echo
echo "############ 2. @CO comment lines — free-text; the most likely home for tissue"
samtools view -H "$BAM" 2>/dev/null | grep '^@CO' || echo "  (no @CO lines)"

echo
echo "############ 3. EVERY distinct @RG tag present (do not assume the tag set)"
samtools view -H "$BAM" 2>/dev/null | grep '^@RG' | tr '\t' '\n' | grep -oE '^[A-Za-z][A-Za-z0-9]:' \
    | sort -u | tr -d ':' | tr '\n' ' '; echo
echo "  -- full @RG lines:"
samtools view -H "$BAM" 2>/dev/null | grep '^@RG' | sed 's/^/    /'

echo
echo "############ 4. @HD / @SQ summary + build fingerprint"
samtools view -H "$BAM" 2>/dev/null | grep '^@HD' | sed 's/^/    /'
samtools view -H "$BAM" 2>/dev/null | grep -c '^@SQ' | xargs echo "    @SQ lines:"
samtools view -H "$BAM" 2>/dev/null | grep '^@SQ' | head -1 | sed 's/^/    /'
samtools view -H "$BAM" 2>/dev/null | grep '^@SQ' | tr '\t' '\n' | grep '^LN:' | sed 's/^LN://' \
    | awk '{s+=$1} END{print "    total reference length:", s}'

echo
echo "############ 5. run provenance buried in @PG CL strings"
samtools view -H "$BAM" 2>/dev/null | grep '^@PG' \
    | grep -oE '[0-9]{6}_[A-Z0-9]+_[0-9]+_[A-Z]_[A-Z0-9]+' | sort -u | sed 's/^/    run folder: /'
samtools view -H "$BAM" 2>/dev/null | grep '^@PG' \
    | grep -oE 'canpipe/live/data/lane/[0-9]+_[0-9]+' | sort -u | sed 's/^/    canpipe lane: /'
samtools view -H "$BAM" 2>/dev/null | grep '^@PG' \
    | grep -oE '/[^[:space:]]*\.(fa|fasta)' | sort -u | sed 's/^/    reference: /'

echo
echo "############ 6. iRODS (cgp zone) — every object for this sample, and its AVUs"
# Query by DATA_NAME. `imeta qu -d sample = <X>` returns "No rows found" in the cgp zone —
# NPG attribute names do not apply here. iquest-by-name assumes nothing and works.
if command -v iquest >/dev/null 2>&1; then
    sm=$(basename "$BAM"); sm="${sm%%.*}"
    echo "  objects named ${sm}* in /cgp:"
    iquest --no-page "SELECT COLL_NAME, DATA_NAME WHERE DATA_NAME like '${sm}%'" 2>&1 \
        | grep -E '^(COLL_NAME|DATA_NAME)' | paste - - | sed 's/^/    /' | head -12
    obj=$(iquest --no-page "SELECT COLL_NAME, DATA_NAME WHERE DATA_NAME like '${sm}%.sample.dupmarked.bam'" 2>/dev/null)
    c=$(printf '%s\n' "$obj" | grep -m1 '^COLL_NAME' | sed 's/^COLL_NAME = //')
    o=$(printf '%s\n' "$obj" | grep -m1 '^DATA_NAME' | sed 's/^DATA_NAME = //')
    if [ -n "${c:-}" ] && [ -n "${o:-}" ]; then
        echo "  AVUs on $c/$o:"
        imeta ls -d "$c/$o" 2>&1 | sed 's/^/    /'
        echo "  (expect exactly: file_extension, id_analysis_proc, id_ifile, md5, version —"
        echo "   none descriptive. Anything MORE than that is new and worth telling Jeremy.)"
    fi
else
    echo "  iquest not on PATH — run: module load IRODS   (CAPITALISED; lowercase does not exist)"
fi
