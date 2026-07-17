#!/bin/bash
# PHASE 3: merge the array's parts into one catalogue TSV and report coverage.
#
#   bash cluster/catalogue_merge.sh      # -> ~/catalogue/catalogue.tsv
#
# Refuses to merge if parts are missing: a catalogue with silent holes is worse than no
# catalogue, because you would query it and believe the answer.
#
# COMPLETENESS IS CHECKED PER SAMPLE, NOT PER PART. The catalogue is assembled from two
# sources — done.tsv (earlier scans) and parts/ (this array) — so "all N parts present" no
# longer means "all manifest samples covered", and counting rows does not either: done.tsv
# can hold rows for samples the CURRENT manifest no longer asks for (a donor removed from
# donors.txt), which inflates the count while real samples are missing. So we ask the only
# question that matters: is every project+sample in the manifest present in the output?
set -uo pipefail

OUT=${OUT:-$HOME/catalogue}
MAN="$OUT/manifest.tsv"
DONE="$OUT/done.tsv"
DEST="$OUT/catalogue.tsv"

[ -s "$MAN" ] || { echo "no manifest: $MAN" >&2; exit 1; }
ROWS=$(( $(wc -l < "$MAN") - 1 ))

TMP="$DEST.tmp.$$"
cat "$DONE" "$OUT"/parts/part.*.tsv 2>/dev/null > "$TMP"
[ -s "$TMP" ] || { echo "no header rows at all: neither $DONE nor $OUT/parts/*" >&2; rm -f "$TMP"; exit 1; }

# Restrict to the current manifest and report what is missing. Dedup on project+sample,
# last row wins — a fresh part supersedes an older done.tsv row for the same sample.
printf 'project\tsample\tdonor\tbam\tbam_bytes\tbai_stale\tassay\tassembly\tasm_name\tref_len\tchr1_md5\tn_seqs\tsort_order\tread_len\tmapped_reads\tunmapped_reads\tmedian_insert\tplatform\tmodel\tcentre\trun_dates\tn_readgroups\tn_libraries\tn_runs\tfirst_run_lane\tdup_tool\tbwa_version\n' > "$DEST"

awk -F'\t' -v man="$MAN" -v dest="$DEST" '
BEGIN {
    while ((getline line < man) > 0) {
        n = split(line, a, "\t"); if (n < 3 || a[1] == "project") continue
        want[a[1] SUBSEP a[2]] = 1; nwant++
    }
}
{
    k = $1 SUBSEP $2
    if (!(k in want)) { extra++; next }        # stale: not in the current scope
    if (!(k in seen)) seen[k] = ++nrow
    row[seen[k]] = $0
}
END {
    for (i = 1; i <= nrow; i++) print row[i] >> dest
    nmiss = 0
    for (k in want) if (!(k in seen)) { nmiss++; if (nmiss <= 20) { split(k, b, SUBSEP); miss = miss "  " b[1] "/" b[2] "\n" } }
    printf("covered %d/%d manifest samples", nrow+0, nwant+0) > "/dev/stderr"
    if (extra) printf("; dropped %d row(s) not in the current manifest", extra) > "/dev/stderr"
    printf("\n") > "/dev/stderr"
    if (nmiss) {
        printf("MISSING %d sample(s) — header never read:\n%s", nmiss, miss) > "/dev/stderr"
        if (nmiss > 20) printf("  ... and %d more\n", nmiss - 20) > "/dev/stderr"
        printf("re-run: bash cluster/catalogue_scan.sh && bash cluster/submit_catalogue.sh\n") > "/dev/stderr"
        exit 1
    }
}' "$TMP"
INCOMPLETE=$?
rm -f "$TMP"

got=$(( $(wc -l < "$DEST") - 1 ))
echo "catalogue: $DEST  ($got rows of $ROWS manifest samples)"
if [ "$INCOMPLETE" -ne 0 ]; then
    echo "REFUSING to call this complete — a catalogue with silent holes is worse than none." >&2
    echo "The file above is written but INCOMPLETE. Do not query it as if it were the cohort." >&2
    exit 1
fi

# Fold this array's parts into done.tsv so the next scan skips them, then drop the parts:
# they are indexed against a todo.tsv that the next scan will replace.
if ls "$OUT"/parts/part.*.tsv >/dev/null 2>&1; then
    cat "$OUT"/parts/part.*.tsv >> "$DONE"
    rm -f "$OUT"/parts/part.*.tsv
    echo "done.tsv: $(wc -l < "$DONE" | tr -d ' ') samples banked for the next scan"
fi
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
