#!/usr/bin/env bash
# Genotype ONE colony BAM against the POOLED contract. Scheduler-agnostic: the submit
# wrapper passes the 1-based line number of the fofn to process. Safe to re-run: skips
# colonies whose output already exists, and writes atomically (tmp -> final) so a killed
# job never leaves a half-written .txt.gz that looks complete.
#
# POOLED, NOT PER-PATIENT — and that is the point. Every colony is genotyped against the
# SAME contract built from all 9 donors, including loci discovered in other donors. That is
# what makes absence interpretable: a locus called in one donor and wild-type in the other
# 8 is a candidate somatic event; a locus called subclonally ACROSS unrelated genomes is an
# artefact. With per-donor contracts a locus is never even tested outside its own donor, so
# "artefact confined to one donor" and "never looked" produce identical output and the FP
# readout quietly disappears. See cluster/combine_mei.sh's header.
#
# Usage:            cluster/genotype_one.sh <index>
# Env (all optional):
#   FOFN         file list             (default: bams.fofn)
#   OUTDIR       output dir            (default: genotypes)
#   CONTRACT     pooled contract       (default: insertions_grch38/mei9x10.genotyping.txt.gz)
#   GENOTYPE_BIN genotype binary       (default: rust/peartree-genotype/target/release/peartree-genotype)
#   GENO_CFG     rust genotype config  (default: cluster/config.genotype.grch38)
#   GENOTYPE_IMPL v2|legacy            (default: v2 = peartree-genotype2, realignment genotyper,
#                                       numeric output; needs COMBINED + GENOME_2BIT)
#   COMBINED     <cohort>.combined.txt.gz (default: CONTRACT with .genotyping[.tprt].txt.gz -> .combined.txt.gz)
#   GENOME_2BIT  reference 2bit/fasta  (default: /lustre/scratch126/casm/teams/team273/users/jd43/hg38.2bit)
#   GENOTYPE2_BIN / GENO2_CFG           (defaults: rust/peartree-genotype2/target/release/peartree-genotype2,
#                                       cluster/config.genotype2.grch38)
set -euo pipefail

IDX="${1:?usage: genotype_one.sh <1-based-index-into-fofn>}"
FOFN="${FOFN:-bams.fofn}"
OUTDIR="${OUTDIR:-genotypes}"
CONTRACT="${CONTRACT:-insertions_grch38/mei9x10.genotyping.txt.gz}"
BIN="${GENOTYPE_BIN:-rust/peartree-genotype/target/release/peartree-genotype}"
CFG="${GENO_CFG:-cluster/config.genotype.grch38}"
IMPL="${GENOTYPE_IMPL:-v2}"
BIN2="${GENOTYPE2_BIN:-rust/peartree-genotype2/target/release/peartree-genotype2}"
CFG2="${GENO2_CFG:-cluster/config.genotype2.grch38}"
GENOME_2BIT="${GENOME_2BIT:-/lustre/scratch126/casm/teams/team273/users/jd43/hg38.2bit}"
COMBINED="${COMBINED:-$(echo "$CONTRACT" | sed -E 's/\.genotyping(\.tprt)?\.txt\.gz$/.combined.txt.gz/')}"

[ -s "$CONTRACT" ] || { echo "no contract: $CONTRACT" >&2; exit 1; }

BAM="$(sed -n "${IDX}p" "$FOFN")"
[ -n "$BAM" ] || { echo "no BAM at line $IDX of $FOFN" >&2; exit 1; }
[ -s "$BAM" ] || { echo "BAM missing/empty: $BAM" >&2; exit 1; }
# Derive the colony id from the FILENAME, not the directory -- same rule as discover_one.sh,
# and for the same reason: the staging layout puts <id> two dirs up but nst_links puts the
# PROJECT there, which silently collapses every colony of a project onto one output name.
ID="$(basename "$BAM")"; ID="${ID%.bam}"; ID="${ID%.cram}"; ID="${ID%.sample.dupmarked}"
[ -n "$ID" ] || { echo "could not derive id from $BAM" >&2; exit 1; }
OUT="$OUTDIR/${ID}.txt.gz"

mkdir -p "$OUTDIR"
if [ -s "$OUT" ]; then
    echo "[$IDX] $ID: output exists, skipping ($OUT)"
    exit 0
fi

# Genotyping needs a coordinate index; discovery did not, so a staged BAM may well not have
# one. Index in place: these are OUR staged Lustre copies, not nst_links (Sanger policy --
# and writing a .bai next to an nst_links symlink would target the iRODS resource server).
if [ ! -s "$BAM.bai" ] && [ ! -s "${BAM%.bam}.bai" ]; then
    echo "[$IDX] $ID: indexing $BAM"
    module load samtools-1.19/python-3.12.0 2>/dev/null || true
    command -v samtools >/dev/null 2>&1 || { echo "[$IDX] $ID: samtools not on PATH — cannot index" >&2; exit 1; }
    samtools index "$BAM" || { echo "[$IDX] $ID: samtools index FAILED" >&2; exit 1; }
fi

TMP="$OUT.tmp.$$"
echo "[$IDX] $ID: genotyping $BAM against $CONTRACT"
if [ "$IMPL" = v2 ]; then
    [ -s "$COMBINED" ] || { echo "[$IDX] $ID: combined consensus missing: $COMBINED" >&2; exit 1; }
    [ -s "$GENOME_2BIT" ] || { echo "[$IDX] $ID: reference missing: $GENOME_2BIT" >&2; exit 1; }
    "$BIN2" --step genotype --bam "$BAM" --insertions "$CONTRACT" --combined "$COMBINED" \
        --reference "$GENOME_2BIT" --out "$TMP" --threads 1 --config "$CFG2"
else
    "$BIN" --step genotype --bam "$BAM" --insertions "$CONTRACT" \
        --out "$TMP" --threads 1 --config "$CFG"
fi
mv -f "$TMP" "$OUT"
echo "[$IDX] $ID: done -> $OUT"
