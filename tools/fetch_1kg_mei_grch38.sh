#!/usr/bin/env bash
# Build the GRCh38 germline-MEI truth set: 1000 Genomes 30x (3,202 samples), sites only.
#
#   bash tools/fetch_1kg_mei_grch38.sh          # -> testdata/mei/1kg.mei.grch38.sites.vcf.gz
#
# Env: OUT (default testdata/mei/1kg.mei.grch38.sites.vcf.gz), URL, KEEP_TMP=1
#
# WHY THIS FILE EXISTS. Our recall canary (testdata/mei/1kg.sv.vcf.gz) is 1000G PHASE 3:
# called in 2015 at low coverage, on hs37d5. The benchmark cohort is GRCh38, so that set
# cannot score it. The obvious move is to liftOver phase 3 hg19->hg38; the better move is not
# to lift a 2015 callset at all.
#
# WHAT WE USE INSTEAD -- and why this is an upgrade, not just a re-coordinate:
#   1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz
#   - GRCh38-NATIVE (chr-prefixed). No liftover, so no lift failures and no silently
#     mis-mapped coordinates to reason about.
#   - 3,202 samples at ~34.5x, vs phase 3's 2,504 at ~7x. MEIs are exactly the class low
#     coverage misses.
#   - MEI calls come from MELT (ALGORITHMS=melt), the same caller lineage as phase 3's
#     *_umary_* MEI calls -- so this is the same measurement redone properly, not a
#     different thing.
#   - Carries AF, which the recall table needs (the phase-3 finding that no COMMON MEI is
#     ever touched by the slippage gate is the load-bearing half of it).
#
# WHY NOT THE 2025 LONG-READ SET (1KG_ONT_VIENNA, Schloissnig et al., Nature 644:442).
# It is more COMPLETE for MEIs -- and that is the problem. It is ONT long-read, so it
# contains insertions that short-read WGS structurally cannot see. Our BAMs are Illumina
# 151bp PE; scoring them against a long-read truth set measures the technology gap, not the
# pipeline. Keep it in mind as a second, clearly-labelled canary; do not mix it into this one.
#
# WHAT IS SELECTED, and what is deliberately NOT:
#   KEEP  ALT = <INS:ME:ALU> | <INS:ME:LINE1> | <INS:ME:SVA> | <INS:ME>
#   DROP  <DEL:ME:...>  -- a DELETION of an element already in the reference. It is not an
#                          insertion, and counting it would credit us for finding the
#                          reference. Phase 3 spelled these DEL_ALU / DEL_LINE1.
#   DROP  <INS> / <INS:UNK> -- novel insertions of unclassified origin. Not known to be an MEI.
#   DROP  everything non-PASS.
# Sites only (fields 1-8): 3,202 genotype columns are ~1.7 GB of nothing we need.
#
# Streams via tabix rather than downloading the whole VCF: same bytes over the wire, but the
# genotypes are discarded as they arrive instead of landing on disk.
set -uo pipefail
cd "$(dirname "$0")/.."

URL="${URL:-http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz}"
OUT="${OUT:-testdata/mei/1kg.mei.grch38.sites.vcf.gz}"

command -v tabix >/dev/null || { echo "tabix not on PATH (brew install htslib)" >&2; exit 1; }
command -v bgzip >/dev/null || { echo "bgzip not on PATH (brew install htslib)" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"

TMP="$(mktemp -d)"; [ "${KEEP_TMP:-0}" = 1 ] || trap 'rm -rf "$TMP"' EXIT
# tabix writes the remote .tbi into $PWD; keep that litter out of the repo root.
cd "$TMP"

echo "source: $URL"
echo "fetching header..."
tabix -H "$URL" > hdr.vcf 2>/dev/null || { echo "cannot read remote header" >&2; exit 1; }
[ -s hdr.vcf ] || { echo "empty header — remote unreadable?" >&2; exit 1; }

# Assert GRCh38, from chr1's LENGTH — never from a label. This is the same fingerprint
# install.sh and the catalogue use: GRCh38 chr1 is 248,956,422; hg19/GRCh37 chr1 is
# 249,250,621. A truth set silently on the wrong build would score every locus as a miss.
c1="$(grep -m1 '^##contig=<ID=chr1,' hdr.vcf | grep -oE 'length=[0-9]+' | cut -d= -f2)"
if [ "${c1:-0}" != "248956422" ]; then
    echo "REFUSING: chr1 length is '${c1:-<absent>}', expected 248956422 (GRCh38)." >&2
    echo "This source is not GRCh38 — do not build a truth set from it." >&2
    exit 1
fi
echo "  chr1 length 248956422 — GRCh38 confirmed"

chroms="$(for i in $(seq 1 22) X Y; do echo "chr$i"; done)"
{
    grep '^##fileformat' hdr.vcf
    echo "##source=1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF (1000G 30x, 3202 samples, GRCh38-native)"
    echo "##source_url=$URL"
    echo "##pt_selection=ALT in {INS:ME:ALU,INS:ME:LINE1,INS:ME:SVA,INS:ME} AND FILTER=PASS; DEL:ME and INS:UNK excluded; sites only"
    grep -E '^##(contig|ALT|INFO|FILTER)=' hdr.vcf
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > sites.vcf

echo "streaming $(echo "$chroms" | wc -w | tr -d ' ') contigs (~1.7 GB over the wire)..."
for c in $chroms; do
    n=$(tabix "$URL" "$c" 2>/dev/null \
        | awk -F'\t' -v OFS='\t' '$7=="PASS" && $5 ~ /^<INS:ME(:|>)/ {print $1,$2,$3,$4,$5,$6,$7,$8}' \
        | tee -a sites.vcf | wc -l | tr -d ' ')
    printf '  %-6s %6s MEIs\n' "$c" "$n"
done

bgzip -f sites.vcf && tabix -f -p vcf sites.vcf.gz
cd - >/dev/null
cp "$TMP/sites.vcf.gz" "$OUT"; cp "$TMP/sites.vcf.gz.tbi" "$OUT.tbi"

echo
echo "wrote $OUT"
echo "== by element class =="
gunzip -c "$OUT" | grep -v '^#' | cut -f5 | sort | uniq -c | sort -rn | sed 's/^/  /'
echo "== total =="
gunzip -c "$OUT" | grep -vc '^#' | sed 's/^/  /'
echo "== AF present? (needed for the common-vs-rare half of the recall table) =="
gunzip -c "$OUT" | grep -v '^#' | grep -c 'AF=' | sed 's/^/  records with AF=: /'
