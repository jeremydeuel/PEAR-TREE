#!/bin/bash
# Survey candidate patients for the 10x10 MEI benchmark: locate their colony BAMs in
# nst_links and — critically — determine which ASSEMBLY each is mapped to.
#
# The benchmark only works on GRCh37/hs37d5, because the 1000G Phase 3 MEI truth set
# (testdata/mei/1kg.sv.vcf.gz) is called on hs37d5 with numeric contigs. A GRCh38 BAM
# would need a different truth set and would silently produce wrong coordinates through
# combine_insertions (see cluster/install.sh's assembly gate).
#
# Read-only and fast — run it on a head node, no bsub needed:
#   bash cluster/survey_patients.sh | tee ~/mei10x10/survey.txt
set -uo pipefail

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH; module load samtools-1.19/python-3.12.0" >&2; exit 1; }

NST=${NST:-/nfs/cancer_ref01/nst_links/live}
PATIENTS=${PATIENTS:-"PD41048 PD48402 PD45534 PD43974 PD40667 PD40521 PD34200 PD37449 PD37590 PD44887 PD44890 PD43947 PD45517"}

# Classify a BAM's reference from its @SQ header.
#   hs37d5  : has an SN:hs37d5 decoy contig            -> GRCh37, what we need
#   GRCh37  : AS:NCBI37 / numeric contigs, no decoy
#   GRCh38  : AS:GRCh38 or SN:chr* with chr1 length 248956422
#   hg19    : SN:chr* with chr1 length 249250621
assembly_of() {
    samtools view -H "$1" 2>/dev/null | awk '
        /^@SQ/ {
            nsq++
            sn=""; ln=""; as=""
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SN:/) sn=substr($i,4)
                if ($i ~ /^LN:/) ln=substr($i,4)
                if ($i ~ /^AS:/) as=substr($i,4)
            }
            if (nsq==1) { first_sn=sn; first_ln=ln }
            if (as != "") seen_as=as
            if (sn=="hs37d5") decoy=1
            if (sn=="chr1")   chr1_ln=ln
            if (sn=="1")      one_ln=ln
        }
        END {
            if (nsq==0)                       { print "NO_HEADER"; exit }
            if (decoy)                        { verdict="hs37d5_GRCh37" }
            else if (chr1_ln==248956422)      { verdict="GRCh38_chr" }
            else if (chr1_ln==249250621)      { verdict="hg19_chr" }
            else if (one_ln==249250621)       { verdict="GRCh37_numeric" }
            else if (one_ln==248956422)       { verdict="GRCh38_numeric" }
            else                              { verdict="UNKNOWN" }
            printf "%s\tnSQ=%d\tfirstSN=%s\tAS=%s", verdict, nsq, first_sn, (seen_as==""?"-":seen_as)
        }'
}

printf "%-9s %-7s %-7s %-16s %s\n" PATIENT PROJECT COLONIES ASSEMBLY DETAIL
echo "--------------------------------------------------------------------------------"

for P in $PATIENTS; do
    # colony sample dirs for this patient, across every project
    mapfile -t DIRS < <(ls -d "$NST"/*/"${P}"[a-z]_lo* 2>/dev/null | sort)
    if [ "${#DIRS[@]}" -eq 0 ]; then
        # fall back to any sample dir with this stem (naming may differ)
        mapfile -t DIRS < <(ls -d "$NST"/*/"${P}"* 2>/dev/null | sort)
    fi
    n=${#DIRS[@]}
    if [ "$n" -eq 0 ]; then
        printf "%-9s %-7s %-7s %-16s %s\n" "$P" "-" 0 "NOT_FOUND" "no dir under $NST/*/${P}*"
        continue
    fi
    PROJ=$(basename "$(dirname "${DIRS[0]}")")

    verdict="NO_BAM"; detail=""
    for d in "${DIRS[@]}"; do
        s=$(basename "$d"); b="$d/$s.sample.dupmarked.bam"
        [ -s "$b" ] || continue
        out=$(assembly_of "$b")
        verdict=$(echo "$out" | cut -f1)
        detail=$(echo "$out" | cut -f2-)
        break
    done
    printf "%-9s %-7s %-7s %-16s %s\n" "$P" "$PROJ" "$n" "$verdict" "$detail"
done

echo
echo "Only 'hs37d5_GRCh37' patients are usable with testdata/mei/1kg.sv.vcf.gz."
echo "Next: bash cluster/build_fofn.sh <patient> [<patient> ...]"
