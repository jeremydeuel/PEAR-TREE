#!/usr/bin/env bash
# Differential test: run the Python and Rust discovery steps on the same BAM and
# diff their (decompressed) outputs. Exit 0 iff byte-identical.
#
# Usage:
#   tests/differential_test.sh <coordinate-sorted.bam> [python] [rust-binary]
#
# This is the equivalence oracle for the Rust port. Point it at a real WGS BAM
# on the cluster to satisfy the Stage 3 exit criterion.
set -euo pipefail

BAM="${1:?usage: differential_test.sh <bam> [python] [rust-binary]}"
PY="${2:-venv/bin/python}"
RS="${3:-rust/peartree-discovery/target/release/peartree-discovery}"

# resolve repo root (two levels up from this script) so relative defaults work
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "[diff-test] python: $PY"
echo "[diff-test] rust:   $RS"
echo "[diff-test] bam:    $BAM"

echo "[diff-test] running python discovery..."
time "$PY" src/main.py --step discover --bam "$BAM" --out "$TMP/py.txt.gz" >"$TMP/py.log" 2>&1

echo "[diff-test] running rust discovery..."
time "$RS" --step discover --bam "$BAM" --out "$TMP/rs.txt.gz" >"$TMP/rs.log" 2>&1

if diff <(gzip -dc "$TMP/py.txt.gz") <(gzip -dc "$TMP/rs.txt.gz") >"$TMP/diff.txt"; then
    n=$(gzip -dc "$TMP/py.txt.gz" | grep -c '^@' || true)
    echo "[diff-test] *** BYTE-IDENTICAL *** ($n records)"
    exit 0
else
    echo "[diff-test] !!! MISMATCH !!! (first 40 diff lines below)"
    head -40 "$TMP/diff.txt"
    echo "[diff-test] full outputs kept: copying to ./diff-test-fail.{py,rs}.txt"
    gzip -dc "$TMP/py.txt.gz" > ./diff-test-fail.py.txt
    gzip -dc "$TMP/rs.txt.gz" > ./diff-test-fail.rs.txt
    exit 1
fi
