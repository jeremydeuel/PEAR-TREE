#!/bin/bash
# Hunt for the new-organ donors on the farm. READ-ONLY. Head node. ~10-20 min with iRODS.
#
#   tsh ssh -D 8080 jd43@farm22-head1
#   cd /lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE && git pull
#   bash cluster/hunt_new_organs.sh 2>&1 | tee ~/hunt_new_organs.txt
#
# WHAT IT ANSWERS, in order:
#   1. Is each donor on the farm at all (nst_links / staging / lustre / iRODS)?
#   2. For those that are: which project, how many samples, WHAT ASSEMBLY.
#   3. How many of the paper's manifest samples we can actually account for.
#
# WHY MANIFEST-DRIVEN. cluster/locate_donors.sh answers "is it here" by globbing the donor
# id. That is the right first question but the WRONG way to build a cohort: a donor is
# re-released under several project ids (PD51634: 1,178 BAMs vs a 151-sample manifest, ~8x
# duplication), so globbing double-counts. The paper manifests in cluster/manifests/ are the
# authority for WHICH samples are real. See cluster/locate_donors.sh header.
#
# ASSEMBLY IS PER-BAM, NOT PER-DONOR OR PER-PROJECT. Project 1903 holds both hs37d5 and
# GRCh38 BAMs. This script reports assembly per (donor x project) as a TRIAGE signal only —
# re-check the @SQ header of every BAM you actually use. A wrong-assembly BAM is silently
# wrong, never loud.
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
NST=${NST:-/nfs/cancer_ref01/nst_links/live}
MAN=$REPO/cluster/manifests
DONORS=$(grep -v '^#' "$REPO/cluster/donors.new_organs.txt" | awk 'NF{print $1}')

module load samtools-1.19/python-3.12.0 2>/dev/null || true
module load IRODS 2>/dev/null || true   # CAPITALISED — lowercase `irods` does not exist

echo "############ STEP 1 — is each donor on the farm at all?"
echo "(nst_links is an INCOMPLETE view of iRODS: PD57333 = 156 in iRODS vs 79 linked."
echo " 'absent from nst_links' NEVER means 'absent from the farm'.)"
echo
DO_IRODS=1 bash "$REPO/cluster/locate_donors.sh" $DONORS

echo
echo "############ STEP 2 — project x assembly triage for donors that ARE in nst_links"
assembly_of() {
    samtools view -H "$1" 2>/dev/null | awk '
        /^@SQ/ { nsq++; sn=""; ln=""
            for (i=1;i<=NF;i++){ if($i~/^SN:/)sn=substr($i,4); if($i~/^LN:/)ln=substr($i,4) }
            if (sn=="hs37d5") decoy=1; if (sn=="chr1") chr1=ln; if (sn=="1") one=ln }
        END { if (nsq==0) print "NO_HEADER"
              else if (decoy) print "hs37d5_GRCh37"
              else if (chr1==248956422) print "GRCh38_chr"
              else if (chr1==249250621) print "hg19_chr"
              else if (one==249250621)  print "GRCh37_numeric"
              else if (one==248956422)  print "GRCh38_numeric"
              else print "UNKNOWN" }'
}
printf "%-9s %-8s %-9s %-16s %s\n" DONOR PROJECT SAMPLES ASSEMBLY EXAMPLE
echo "---------------------------------------------------------------------------"
for d in $DONORS; do
    # PD ids are not prefix-free (PD5163 vs PD51632): drop hits where a DIGIT follows the id.
    ls -d "$NST"/*/"$d"* 2>/dev/null | awk -F/ -v d="$d" '{if ($NF !~ "^"d"[0-9]") print}' \
        | awk -F/ '{print $(NF-1) "\t" $NF}' | sort -u > /tmp/hunt.$d.tsv
    [ -s /tmp/hunt.$d.tsv ] || { printf "%-9s %-8s %-9s %-16s %s\n" "$d" - 0 NOT_IN_NST -; continue; }
    cut -f1 /tmp/hunt.$d.tsv | sort -u | while read -r proj; do
        n=$(awk -F'\t' -v p="$proj" '$1==p' /tmp/hunt.$d.tsv | wc -l | tr -d ' ')
        s=$(awk -F'\t' -v p="$proj" '$1==p{print $2; exit}' /tmp/hunt.$d.tsv)
        b="$NST/$proj/$s/$s.sample.dupmarked.bam"
        # header-only read: the ONLY thing nst_links may be read for. Never stream records.
        if [ -s "$b" ]; then a=$(assembly_of "$b"); else a="NO_BAM(analysis release?)"; fi
        printf "%-9s %-8s %-9s %-16s %s\n" "$d" "$proj" "$n" "$a" "$s"
    done
done

echo
echo "############ STEP 3 — manifest coverage: how many PAPER samples can we find?"
for m in "$MAN"/coorens2025stomach.wgs.tsv "$MAN"/oliver2025nf1.wgs.tsv; do
    [ -s "$m" ] || continue
    echo "--- $(basename "$m")"
    # column 2 = sample. NF1 manifest additionally has arm/site/organ_group/status.
    tail -n +2 "$m" | awk -F'\t' '{print $1"\t"$2}' | sort -u > /tmp/hunt.man.tsv
    tot=$(wc -l < /tmp/hunt.man.tsv | tr -d ' ')
    found=0
    while IFS=$'\t' read -r d s; do
        ls -d "$NST"/*/"$s" >/dev/null 2>&1 && found=$((found+1))
    done < /tmp/hunt.man.tsv
    echo "    manifest samples: $tot   present in nst_links: $found   missing: $((tot-found))"
    echo "    (missing != absent — check iRODS with iquest before concluding EGA-only)"
done

cat <<'EOF'

############ WHAT TO DO WITH THE ANSWER

If the donors ARE on the farm:
  - STAGE them, never compute off nst_links (Sanger policy; nst_links symlinks into the
    iRODS resource servers and streaming whole BAMs hammers shared archive infrastructure):
      module load dataImportExport
      stageBam.pl --lustre 126 --types m --sample <file-of-sample-names> \
                  --project <PROJECT> -o /lustre/scratch126/casm/staging/team273/jd43 -fo
  - Build the sample file from cluster/manifests/*.tsv, NOT from a glob.

PRIORITY: the 305 NORMAL BRAIN LCM whole genomes in oliver2025nf1.wgs.tsv:
    awk -F'\t' 'NR>1 && $3=="WGS_LCM" && $6=="BRAIN" && $7=="NORMAL" {print $2}' \
        cluster/manifests/oliver2025nf1.wgs.tsv > ~/nf1_brain.samples
Brain is the tissue where somatic L1 retrotransposition is actually reported. 193 of the 838
NF1 LCM samples are TUMOUR (glioma) and 34 have an unknown site — both are excluded by the
filter above. Do not let them into a "normal tissue" cohort.

If a donor is NOWHERE: it is EGA-only (EGAD00001015398 NF1 / EGAD00001015351 stomach /
EGAD00001009812 Wilms) and needs a data-access request. Tissue and library_type CANNOT be
resolved on the farm — the NPG `seq` zone is not visible from farm22, only `cgp`.
EOF
