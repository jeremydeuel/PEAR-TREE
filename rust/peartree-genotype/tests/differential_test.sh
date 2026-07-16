#!/usr/bin/env bash
# Differential test: run the Python and Rust genotyping steps on the same BAM +
# insertion contract and diff their (decompressed) outputs. Exit 0 iff identical.
#
# Usage:
#   tests/differential_test.sh <bam|cram> <ins.genotyping.txt.gz> [python] [rust-binary] [threads] [reference.fa]
#
# This is the equivalence oracle for the Rust genotyping port. The gzip *bytes*
# differ (compression level / header mtime); only the decompressed table is
# compared, which is the on-disk contract downstream reads. Pass a reference FASTA
# as the 6th argument for external-reference CRAM input.
set -euo pipefail

BAM="${1:?usage: differential_test.sh <bam|cram> <ins.genotyping.txt.gz> [python] [rust-binary] [threads] [ref.fa]}"
INS="${2:?usage: differential_test.sh <bam|cram> <ins.genotyping.txt.gz> [python] [rust-binary] [threads] [ref.fa]}"
PY="${3:-python3}"
RS="${4:-rust/peartree-genotype/target/release/peartree-genotype}"
THREADS="${5:-1}"
REF="${6:-}"
RS_REF=()
if [ -n "$REF" ]; then RS_REF=(--reference "$REF"); fi

# resolve repo root (three levels up from this script) so relative defaults work
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "[diff-test] python:  $PY"
echo "[diff-test] rust:    $RS"
echo "[diff-test] bam:     $BAM"
echo "[diff-test] contract:$INS"
echo "[diff-test] threads: $THREADS"

echo "[diff-test] running python genotype..."
time "$PY" src/main.py --step genotype --bam "$BAM" --insertions "$INS" \
    --out "$TMP/py.txt.gz" --threads "$THREADS" >"$TMP/py.log" 2>&1

echo "[diff-test] running rust genotype..."
time "$RS" --step genotype --bam "$BAM" --insertions "$INS" \
    "${RS_REF[@]}" --out "$TMP/rs.txt.gz" --threads "$THREADS" >"$TMP/rs.log" 2>&1

if diff <(gzip -dc "$TMP/py.txt.gz") <(gzip -dc "$TMP/rs.txt.gz") >"$TMP/diff.txt"; then
    n=$(gzip -dc "$TMP/py.txt.gz" | tail -n +2 | wc -l | tr -d ' ')
    echo "[diff-test] *** BYTE-IDENTICAL *** ($n loci)"
    exit 0
else
    echo "[diff-test] !!! MISMATCH !!! (first 40 diff lines below)"
    head -40 "$TMP/diff.txt"
    echo "[diff-test] full outputs copied to ./diff-test-fail.geno.{py,rs}.txt"
    gzip -dc "$TMP/py.txt.gz" > ./diff-test-fail.geno.py.txt
    gzip -dc "$TMP/rs.txt.gz" > ./diff-test-fail.geno.rs.txt
    exit 1
fi
