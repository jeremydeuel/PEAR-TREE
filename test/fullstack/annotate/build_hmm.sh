#!/bin/bash
# Build the nucleotide HMM library annotate_v2.py scans inserted-sequence clips against
# (CONFIG['annotate']['hmm']). One model per implanted RTE family, built from the exact hs1
# source locus the test donor samples (test/genotyping/build_haplotypes.py ELEMENT_LOCI), so
# the library covers precisely the families the harness implants. nhmmscan --dfamtblout then
# classifies each clip; annotate's conclusion() turns a family hit + poly-A/reciprocal side
# into the L1/Alu/SVA/HERVK call.
#
# rte_elements.fa is committed, so the default path needs only HMMER (hmmbuild/hmmpress) — no
# hs1.fa. Pass --from-hs1 <hs1.fa> to re-extract the source sequences (e.g. to add a family).
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
FA="$DIR/rte_elements.fa"
HMM="$DIR/peartree_rte.hmm"

# hs1 source loci (1-based inclusive), same as build_haplotypes.py ELEMENT_LOCI.
LOCI=( "L1HS:chrX:11289796-11295826" "AluYa5:chrX:2759615-2759925"
       "HERVK:chr10:5039043-5047200" "SVA_E:chr3:148851419-148853373"
       "SVA_F:chr18:54545509-54547673" )

if [ "${1:-}" = "--from-hs1" ]; then
    HS1="${2:?usage: build_hmm.sh --from-hs1 <hs1.fa>}"
    : > "$FA"
    for spec in "${LOCI[@]}"; do
        name="${spec%%:*}"; region="${spec#*:}"
        seq="$(samtools faidx "$HS1" "$region" | tail -n +2 | tr -d '\n')"
        printf ">%s\n%s\n" "$name" "$seq" >> "$FA"
    done
    echo "re-extracted $(grep -c '>' "$FA") element sources into $FA"
fi

[ -f "$FA" ] || { echo "missing $FA (run with --from-hs1 <hs1.fa> once to create it)"; exit 1; }
samtools faidx "$FA"
tmp="$(mktemp -d)"
: > "$HMM"
for name in $(grep '>' "$FA" | tr -d '>'); do
    samtools faidx "$FA" "$name" > "$tmp/$name.fa"
    hmmbuild --dna -n "$name" "$tmp/$name.hmm" "$tmp/$name.fa" > /dev/null
    cat "$tmp/$name.hmm" >> "$HMM"
done
rm -rf "$tmp"
hmmpress -f "$HMM" > /dev/null
echo "built + pressed $(grep -c '^NAME' "$HMM") models -> $HMM"
grep '^NAME' "$HMM"
