#!/bin/bash
# PHASE 3: merge the array's parts into one catalogue TSV and report coverage.
#
#   bash cluster/catalogue_merge.sh      # -> ~/catalogue/catalogue.tsv
#
# Refuses to merge if parts are missing: a catalogue with silent holes is worse than no
# catalogue, because you would query it and believe the answer.
set -uo pipefail

OUT=${OUT:-$HOME/catalogue}
CHUNK=${CHUNK:-500}
MAN="$OUT/manifest.tsv"
DEST="$OUT/catalogue.tsv"

[ -s "$MAN" ] || { echo "no manifest: $MAN" >&2; exit 1; }
ROWS=$(( $(wc -l < "$MAN") - 1 ))
N=$(( (ROWS + CHUNK - 1) / CHUNK ))

missing=()
for i in $(seq 1 "$N"); do [ -s "$OUT/parts/part.$i.tsv" ] || missing+=("$i"); done
if [ "${#missing[@]}" -gt 0 ]; then
    echo "MISSING ${#missing[@]}/$N parts: ${missing[*]:0:20}${missing[20]:+ ...}" >&2
    echo "re-run just those: for i in ${missing[*]:0:20}; do bash cluster/catalogue_headers.sh \$i; done" >&2
    echo "(or resubmit the array — completed parts are skipped)" >&2
    exit 1
fi

printf 'project\tsample\tdonor\tbam\tbam_bytes\tbai_stale\tassay\tassembly\tasm_name\tref_len\tchr1_md5\tn_seqs\tsort_order\tread_len\tmapped_reads\tunmapped_reads\tmedian_insert\tplatform\tmodel\tcentre\trun_dates\tn_readgroups\tn_libraries\tn_runs\tfirst_run_lane\tdup_tool\tbwa_version\n' > "$DEST"
cat "$OUT"/parts/part.*.tsv >> "$DEST"

got=$(( $(wc -l < "$DEST") - 1 ))
echo "catalogue: $DEST  ($got rows of $ROWS manifest samples)"
[ "$got" -eq "$ROWS" ] || echo "WARNING: row count $got != manifest $ROWS — investigate before trusting queries" >&2
echo
echo "== assay =="
awk -F'\t' 'NR>1{c[$7]++} END{for(k in c) printf "  %-28s %6d\n", k, c[k]}' "$DEST" | sort -k2 -rn
echo "== assembly =="
awk -F'\t' 'NR>1{c[$8]++} END{for(k in c) printf "  %-28s %6d\n", k, c[k]}' "$DEST" | sort -k2 -rn | head
echo "== reference builds seen (asm_name / ref_len / chr1 M5) =="
# ref_len and M5 fingerprint the build; AS: is the header's CLAIM and can disagree.
awk -F'\t' 'NR>1 && $8!="-"{c[$8"\t"$9"\t"$10"\t"substr($11,1,8)]++}
     END{for(k in c) printf "  %-40s %6d\n", k, c[k]}' "$DEST" | sort -k5 -rn | head
# Select on ref_len ($10), the summed @SQ LN — a fingerprint of the actual reference, not a
# label we derived. An earlier version keyed on $8=="chr-style" and reported ZERO GRCh38
# donors while 9,925 GRCh38 BAMs sat in the file: a pipefail/SIGPIPE bug had mislabelled them
# "chr1". Fingerprint > label; the label is a convenience, this is the query that decides runs.
#   GRCh38_full_analysis_set_plus_decoy_hla = 3217346917   (AS: reads NCBI38 *or* GRCh38)
#   hs37d5                                  = 3137454505
echo "== WGS + GRCh38 donors with >=10 colonies (the benchmark target) =="
awk -F'\t' 'NR>1 && $7 ~ /^WGS/ && $10==3217346917 && $3!="-" && $15>=50000000 {c[$3"\t"$1]++; rl[$3"\t"$1]=$14}
     END{for(k in c) if(c[k]>=10) printf "  %-24s %4d  read_len=%s\n", k, c[k], rl[k]}' "$DEST" | sort -k2 -rn | head -40
echo "== WGS + hs37d5 donors with >=10 colonies (remap candidates; check read_len first) =="
awk -F'\t' 'NR>1 && $7 ~ /^WGS/ && $10==3137454505 && $3!="-" && $15>=50000000 {c[$3"\t"$1]++; rl[$3"\t"$1]=$14}
     END{for(k in c) if(c[k]>=10) printf "  %-24s %4d  read_len=%s\n", k, c[k], rl[k]}' "$DEST" | sort -k2 -rn | head -20
echo "== any WGS BAM whose ref_len matches NEITHER known build (investigate before use) =="
awk -F'\t' 'NR>1 && $7 ~ /^WGS/ && $10!=3217346917 && $10!=3137454505 && $10>0 {c[$8"\t"$9"\t"$10]++}
     END{for(k in c) printf "  %-44s %6d\n", k, c[k]}' "$DEST" | sort -k4 -rn | head
echo
echo "copy back:  tsh scp -l jd43 farm22-head1:$DEST ."
