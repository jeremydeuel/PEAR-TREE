#!/bin/bash
# Survey candidate patients at (patient x project) granularity.
#
# WHY: assembly is a property of the PROJECT (release era), not the patient. A patient's
# colonies are re-released under several project ids -- PD41048 has dirs in 1903, 1957,
# 2073, 2314 and 2441 -- and those releases are NOT all the same assembly. Picking the
# first project alphabetically (what survey_patients.sh did) reports the assembly of a
# 10-colony backwater while the 300-colony release may be GRCh38.
#
# Also reports DISTINCT colony names per project, because the same colony is linked under
# multiple projects (PD43947: 1232 dirs, 639 distinct) and a naive fofn would list it twice.
#
# Read-only, head-node:
#   bash ~/survey_projects.sh | tee ~/mei10x10/survey_projects.txt
set -uo pipefail

module load samtools-1.19/python-3.12.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

NST=${NST:-/nfs/cancer_ref01/nst_links/live}
PATIENTS=${PATIENTS:-"PD41048 PD43974 PD34200 PD37449 PD43947"}
TMP=${TMP:-$HOME/mei10x10/.survey}
mkdir -p "$TMP"

assembly_of() {
    samtools view -H "$1" 2>/dev/null | awk '
        /^@SQ/ {
            nsq++; sn=""; ln=""
            for (i=1;i<=NF;i++){ if($i~/^SN:/)sn=substr($i,4); if($i~/^LN:/)ln=substr($i,4) }
            if (sn=="hs37d5") decoy=1
            if (sn=="chr1")   chr1=ln
            if (sn=="1")      one=ln
        }
        END {
            if (nsq==0)                  print "NO_HEADER"
            else if (decoy)              print "hs37d5_GRCh37"
            else if (chr1==248956422)    print "GRCh38_chr"
            else if (chr1==249250621)    print "hg19_chr"
            else if (one==249250621)     print "GRCh37_numeric"
            else if (one==248956422)     print "GRCh38_numeric"
            else                         print "UNKNOWN"
        }'
}

printf "%-9s %-8s %-9s %-16s %s\n" PATIENT PROJECT COLONIES ASSEMBLY EXAMPLE_SAMPLE
echo "---------------------------------------------------------------------------"

for P in $PATIENTS; do
    # every dir for this patient, as  project<TAB>sample
    ls -d "$NST"/*/"${P}"* 2>/dev/null \
        | awk -F/ '{print $(NF-1) "\t" $NF}' | sort -u > "$TMP/$P.tsv"
    if [ ! -s "$TMP/$P.tsv" ]; then
        printf "%-9s %-8s %-9s %-16s %s\n" "$P" "-" 0 "NOT_FOUND" "-"
        continue
    fi

    cut -f1 "$TMP/$P.tsv" | sort -u | while read -r proj; do
        n=$(awk -F'\t' -v p="$proj" '$1==p' "$TMP/$P.tsv" | wc -l | tr -d ' ')
        s=$(awk -F'\t' -v p="$proj" '$1==p{print $2; exit}' "$TMP/$P.tsv")
        b="$NST/$proj/$s/$s.sample.dupmarked.bam"
        if [ -s "$b" ]; then a=$(assembly_of "$b"); else a="NO_BAM"; fi
        printf "%-9s %-8s %-9s %-16s %s\n" "$P" "$proj" "$n" "$a" "$s"
    done
done

echo
echo "Pick, per patient, ONE project that is hs37d5_GRCh37 and has >=10 colonies."
echo "Sample names above also show the colony naming convention (PD43974 has no _lo dirs)."
