#!/bin/bash
# Stage the 9x10 benchmark colonies onto Lustre, then build the discovery fofn.
#
#   bash cluster/stage_picked.sh            # resolve + print the stageBam.pl commands
#   bash cluster/stage_picked.sh --fofn     # after staging: verify + emit the fofn
#
# SANGER POLICY (why this script exists at all): BAM records are ALWAYS read from Lustre,
# NEVER from nst_links. nst_links entries are symlinks into the iRODS resource servers
# (/nfs/irods-cgp-<server>-sdf/...); streaming whole BAMs off them hammers shared archive
# infrastructure. Header reads (view -H / idxstats) are the ONLY permitted direct access --
# which is exactly, and only, what catalogue_headers.sh does. "It's already in nst_links so
# we don't need to stage" is a WRONG conclusion: availability is not the question, which
# filesystem absorbs ~2-3 TB of sequential reads is.
#
# WHAT IS BEING STAGED. 9 donors x 10 colonies, all blood/HSPC, all GRCh38, all 151bp PE.
# Per donor: 5 CLADE colonies (deepest internal node with >=5 tips -> shared trunk carries
# the SOMATIC true positives) + 5 SPREAD colonies (greedy farthest-point -> coalesce with the
# clade only at the root, so they share GERMLINE MEIs but no somatic ones, and they keep FP
# sampling independent). See cluster/trees/pick_clades.R for the selection and its caveats.
#
# NOT INCLUDED, deliberately: PD44579/PX002 (every filter we have was tuned on it -- it is a
# positive control, not a subject) and AX001/BMH1_TG (its tree tips are plate/well ids that
# nothing maps to sample names). Both are documented in pick_clades.R.
set -uo pipefail

PT_ROOT=${PT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
PICKED=${PICKED:-$PT_ROOT/cluster/trees/picked.tsv}
CAT=${CAT:-$HOME/catalogue/catalogue.tsv}
STAGE=${STAGE:-/lustre/scratch126/casm/staging/team273/jd43}
WORK=${WORK:-$HOME/catalogue/picked}
MODE=${1:-plan}

[ -s "$PICKED" ] || { echo "no picks: $PICKED (run: Rscript cluster/trees/pick_clades.R > $PICKED)" >&2; exit 1; }
[ -s "$CAT" ]    || { echo "no catalogue: $CAT (run catalogue_scan.sh + submit_catalogue.sh + catalogue_merge.sh)" >&2; exit 1; }
mkdir -p "$WORK"

# SCHEMA GATE. An 18-col catalogue exists in the wild (pre stat -L fix) whose columns MEAN
# different things: col 8 held the raw first SN and col 15 was n_libraries, so the queries
# below would fail on every row and report "0 resolved" with no hint why. Check the header
# width rather than let a column-index mismatch masquerade as "these colonies don't exist".
# (NO_BAM rows legitimately carry a trailing tab -> NF=28; that is blank()'s padding, not a
# schema problem, so gate on the HEADER only.)
NCOL=$(head -1 "$CAT" | awk -F'\t' '{print NF}')
if [ "$NCOL" -ne 27 ]; then
    echo >&2 "WRONG CATALOGUE SCHEMA: $CAT header has $NCOL columns, expected 27."
    echo >&2 "  header: $(head -1 "$CAT" | tr '\t' ' ')"
    echo >&2 "Rebuild it:  rm -rf ~/catalogue/parts && bash cluster/submit_catalogue.sh"
    echo >&2 "             # wait for the array, then:"
    echo >&2 "             bash cluster/catalogue_merge.sh"
    exit 1
fi
[ "$(head -1 "$CAT" | cut -f10)" = "ref_len" ] && [ "$(head -1 "$CAT" | cut -f12)" = "n_seqs" ] || {
    echo >&2 "catalogue col 10/12 are not ref_len/n_seqs -- the GRCh38 fingerprint test below"
    echo >&2 "would read the wrong fields. Header: $(head -1 "$CAT" | tr '\t' ' ')"; exit 1; }

# Resolve each picked tip -> (project, sample) using the CATALOGUE, never a filesystem glob.
# A glob across projects is how you end up staging a TARGETED release or a stale re-release:
# PD51634 has 1,178 BAMs across seven projects for 151 real samples. The catalogue carries the
# measured assay/assembly per BAM; the filesystem carries only a name.
#
# Column map (27-col schema): 1 project, 2 sample, 3 donor, 4 bam, 7 assay, 8 assembly,
# 10 ref_len, 12 n_seqs, 14 read_len, 15 mapped_reads.
#
# WE KEY GRCh38 ON THE REFERENCE FINGERPRINT, NOT THE assembly LABEL. The label is written by
# catalogue_headers.sh and has meant different things across versions: the farm's catalogue
# (built 2026-07-17 from a pre-pull copy of the script) says "chr1", while the current repo
# says "chr-style". Matching the label made this script resolve 0 of 90 colonies and blame
# the cohort. ref_len + n_seqs come straight from the header's @SQ block and cannot be
# relabelled by a script bug:
#   GRCh38_full_analysis_set_plus_decoy_hla -> n_seqs 3366, ref_len 3217346917 (9,925 rows,
#     ALL identical -- one uniform build, no hg19 mixed in; hg19 would be chr1 LN 249250621)
#   hs37d5                                  -> n_seqs   86, ref_len 3137454505 (12,692 rows)
# Never match on the literal string "GRCh38" either: the header's AS: tag reads "NCBI38".
GRCH38_REFLEN=3217346917
GRCH38_NSEQS=3366
awk -F'\t' -v OFS='\t' -v RL="$GRCH38_REFLEN" -v NS="$GRCH38_NSEQS" '
    NR==FNR { if (FNR>1) { role[$8]=$3; donor[$8]=$1 } ; next }
    FNR==1  { next }
    ($2 in role) && $7 ~ /^WGS/ && $10==RL && $12==NS && $15>=50000000 {
        if ($2 in seen) { dup[$2]=dup[$2] "," $1; next }
        seen[$2]=1; print donor[$2], role[$2], $1, $2, $4
    }
    END { for (s in dup) print "DUPLICATE-PROJECT\t" s "\t" dup[s] > "/dev/stderr" }
' "$PICKED" "$CAT" | sort -k1,1 -k2,2 -k4,4 > "$WORK/resolved.tsv"

WANT=$(( $(wc -l < "$PICKED") - 1 ))
GOT=$(wc -l < "$WORK/resolved.tsv" | tr -d ' ')
echo "picked colonies: $WANT   resolved to a GRCh38 WGS BAM: $GOT"
if [ "$GOT" -ne "$WANT" ]; then
    echo >&2 "UNRESOLVED -- refusing to proceed. A silently short cohort is worse than none."
    echo >&2
    # NB: match on FILENAME, not the NR==FNR idiom. When resolved.tsv is EMPTY (the 0-resolved
    # case, i.e. exactly when you need this message most) NR==FNR stays true for every line of
    # PICKED, so all 90 land in the load branch and nothing is ever printed. That bug made a
    # real failure look like a silent one.
    awk -F'\t' -v RES="$WORK/resolved.tsv" '
        FILENAME==RES { r[$4]=1; next }
        FNR>1 && !($8 in r) { print "  unresolved: " $1 "\t" $8 }
    ' "$WORK/resolved.tsv" "$PICKED" >&2
    echo >&2
    # Say WHY. Each picked tip fails for a reason the catalogue can state; guessing between
    # "stale catalogue" and "tip isn't a sample" wasted a round-trip once already.
    echo >&2 "why -- what the catalogue says about each unresolved tip:"
    awk -F'\t' -v RES="$WORK/resolved.tsv" -v CAT="$CAT" '
        FILENAME==RES { r[$4]=1; next }
        FILENAME==CAT { if (FNR>1) { seen[$2]=1
                                     asm[$2]=asm[$2] $8 "(reflen=" $10 ",nseq=" $12 ") "; assay[$2]=assay[$2] $7 " "
                                     dep[$2]=dep[$2] $15 " " } ; next }
        FNR>1 && !($8 in r) {
            if (!($8 in seen)) print "  " $8 ": NO ROW in the catalogue at all (tip is not a sample name?)"
            else printf "  %s: has rows -- assay=[%s] assembly=[%s] mapped=[%s]\n", $8, assay[$8], asm[$8], dep[$8]
        }
    ' "$WORK/resolved.tsv" "$CAT" "$PICKED" >&2 | head -20
    echo >&2
    echo >&2 "Reading the above: we require assay ~ /^WGS/, the GRCh38 fingerprint"
    echo >&2 "(ref_len=3217346917, n_seqs=3366 -- NOT the assembly label), and mapped >= 50M."
    echo >&2 "  - \"NO ROW\" for every tip      -> the donor was not in the scan scope, or the tip"
    echo >&2 "                                   is a plate/well id (AX001 is the known case)."
    echo >&2 "  - rows exist but assembly is  -> wrong release picked up; check the projects."
    echo >&2 "    hs37d5 for all of them"
    exit 1
fi

echo
echo "== per donor =="
awk -F'\t' '{c[$1"\t"$2]++} END{for(k in c) print "  " k "\t" c[k]}' "$WORK/resolved.tsv" | sort

if [ "$MODE" = "--fofn" ]; then
    # Post-staging: verify every colony landed on Lustre, then emit the fofn from STAGED paths.
    # Build into a TEMP and only publish on success. Writing $FOFN directly would leave a
    # PARTIAL fofn on disk after a failed run -- 89 of 90 lines, looking complete to the next
    # person who runs discovery by hand. That is precisely the silent-short-cohort failure
    # this script exists to prevent, so the fofn must never exist unless it is whole.
    FOFN="$WORK/all.bams.fofn"; TMP="$FOFN.tmp.$$"; : > "$TMP"; miss=0
    while IFS=$'\t' read -r donor role proj sample src; do
        b="$STAGE/$proj/$sample/$sample.sample.dupmarked.bam"
        if [ -s "$b" ]; then printf '%s\n' "$b" >> "$TMP"
        else echo "  NOT STAGED: $proj/$sample" >&2; miss=$((miss+1)); fi
    done < "$WORK/resolved.tsv"
    if [ "$miss" -gt 0 ]; then
        rm -f "$TMP" "$FOFN"      # also drop any fofn from an earlier, now-superseded run
        echo >&2 "$miss/$GOT colonies missing from $STAGE -- no fofn written."
        echo >&2 "Re-run the stageBam.pl command for the affected project(s), then retry --fofn."
        exit 1
    fi
    mv -f "$TMP" "$FOFN"
    echo
    echo "fofn: $FOFN  ($(wc -l < "$FOFN") BAMs, all on Lustre)"
    echo "next: PT_ROOT=$PT_ROOT BAMS=$FOFN bash cluster/submit_combine_mei.sh"
    exit 0
fi

# PLAN mode: write one sample-list per project and print the staging commands.
# Printed, not run: staging is a separately-submitted step in this workflow (downstream jobs
# hang off it with -w ended(...)), and it moves terabytes -- it should be a deliberate act.
rm -f "$WORK"/proj.*.samples
awk -F'\t' -v w="$WORK" '{print $4 >> (w "/proj." $3 ".samples")}' "$WORK/resolved.tsv"

echo
echo "== staging commands (run these; each moves ~10 x 30GB) =="
echo "module load dataImportExport/1.60.8"
for f in "$WORK"/proj.*.samples; do
    p=$(basename "$f" .samples); p=${p#proj.}
    printf 'stageBam.pl --lustre 126 --types m --project %s --sample %s -o %s -fo   # %s colonies\n' \
        "$p" "$f" "$STAGE" "$(wc -l < "$f" | tr -d ' ')"
done
echo
echo "then: bash cluster/stage_picked.sh --fofn"
