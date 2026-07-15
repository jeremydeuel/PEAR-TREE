#!/usr/bin/env bash
# Build the synthetic BAMs and run the Python-vs-Rust differential on each, plus
# the checked-in test_data/test.bam. This is the local byte-identity guard for the
# byte-identical refactors (SPD-*/OBS-*); it does NOT prove speed or exercise the
# real-WGS branches — run differential_test.sh on a cluster WGS BAM for that.
#
# Usage: tests/differential_all.sh [python]
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
PY="${1:-venv/bin/python}"
cd "$ROOT"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# checked-in real test BAM + the three synthetic generators
BAMS=("test_data/test.bam")
for g in multicontig polya xa; do
    "$PY" "$HERE/gen_$g.py" "$TMP/$g.bam" >/dev/null
    BAMS+=("$TMP/$g.bam")
done

fail=0
for b in "${BAMS[@]}"; do
    if out="$(bash "$HERE/differential_test.sh" "$b" 2>&1)"; then
        printf '  %-28s %s\n' "$(basename "$b")" "$(echo "$out" | grep -oE 'BYTE-IDENTICAL \*\*\* \([0-9]+ records\)')"
    else
        printf '  %-28s MISMATCH\n' "$(basename "$b")"
        echo "$out" | grep -A20 MISMATCH || true
        fail=1
    fi
done

if [ "$fail" -eq 0 ]; then
    echo "[diff-all] all byte-identical"
else
    echo "[diff-all] FAILURES above"
    exit 1
fi
