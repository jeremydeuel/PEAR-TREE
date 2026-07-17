#!/usr/bin/env bash
# A/B the homopolymer-slippage gate on the 9x10 GRCh38 cohort: submit BOTH discovery arms.
#
#   bash cluster/submit_discovery_ab.sh
#
# Env: FOFN (default ~/catalogue/picked/all.bams.fofn), THROTTLE, MEM, QUEUE, GROUP
#
#   arm A  slippage_filter = true   -> $OUT_A  (the shipped setting, production candidate)
#   arm B  slippage_filter = false  -> $OUT_B  (baseline; the honest superset)
#
# The two configs differ in EXACTLY ONE LINE -- this script asserts that before submitting.
# If anything else drifts between them, the experiment stops measuring the gate and starts
# measuring the drift, and nothing downstream would tell you.
#
# WHY TWO FULL PASSES. The gate drops clip clusters inside discovery
# (discovery.rs: `retain(!is_slippage_clip(..))`); there is no flag-don't-drop mode, and
# hallmarks.tsv is an unrelated feature (TSD/poly-A/EN motif). So the arms cannot be derived
# from one pass -- 180 BAM scans, not 90. If this A/B becomes routine, the cheaper fix is a
# Rust change emitting a slippage flag column instead of dropping, so one pass yields both.
#
# WHAT NOT TO DO DOWNSTREAM: arm B is a SUPERSET of arm A (the gate only ever removes), so do
# NOT genotype both arms -- that is 90 colonies x ~30k loci twice for nothing. Genotype the
# UNION (arm B's contract) once, then evaluate arm A by subsetting those calls to arm A's
# loci. Identical genotype calls in both arms, differing only in the locus set, which is
# exactly the comparison we want. Verify the subset relation before relying on it: dropping
# evidence can in principle shift a consensus and move a breakpoint, so confirm A's loci are
# actually contained in B's rather than assuming it.
#
# WHAT THIS CAN AND CANNOT ANSWER. Without a truth set it compares contract SIZE and
# cross-donor consistency -- i.e. how much the gate cuts, and whether what it cuts looks like
# artefact. It CANNOT tell you what recall the gate costs: that needs the 1kG MEI truth set
# lifted from hs37d5 to GRCh38. The shipped 0.6% (7/1144) figure was measured on ONE donor
# (PD44579), on hs37d5, and does not transfer to a 9-donor GRCh38 cohort. Until the liftover
# exists, arm A vs arm B is "how much is cut", not "is cutting it right".
set -euo pipefail
cd "$(dirname "$0")/.."

FOFN="${FOFN:-$HOME/catalogue/picked/all.bams.fofn}"
CFG_A="${CFG_A:-cluster/config.discovery.grch38}"
CFG_B="${CFG_B:-cluster/config.discovery.grch38.noslip}"
OUT_A="${OUT_A:-discovery_grch38_slip}"
OUT_B="${OUT_B:-discovery_grch38_noslip}"

[ -s "$FOFN" ]  || { echo "no fofn: $FOFN (run cluster/stage_picked.sh --fofn)" >&2; exit 1; }
[ -s "$CFG_A" ] || { echo "no config: $CFG_A" >&2; exit 1; }
[ -s "$CFG_B" ] || { echo "no config: $CFG_B" >&2; exit 1; }

# The arms must differ ONLY by slippage_filter. Compare the SETTINGS (comments and blank
# lines stripped), not the files -- the headers legitimately differ.
d="$(diff <(grep -vE '^[[:space:]]*#|^[[:space:]]*$' "$CFG_A") \
          <(grep -vE '^[[:space:]]*#|^[[:space:]]*$' "$CFG_B") || true)"
nchg=$(printf '%s\n' "$d" | grep -c '^[<>]' || true)
nslip=$(printf '%s\n' "$d" | grep '^[<>]' | grep -c 'slippage_filter' || true)
if [ "$nchg" -eq 0 ]; then
    # The likeliest mistake, and it would silently produce two IDENTICAL arms and a null
    # result -- which reads as "the gate does nothing" rather than "the experiment did not run".
    echo >&2 "ARMS ARE IDENTICAL — refusing to submit."
    echo >&2 "  $CFG_A and $CFG_B have the same settings, so both would run the same gate and"
    echo >&2 "  the A/B would show no difference for the wrong reason."
    grep -H '^slippage_filter' "$CFG_A" "$CFG_B" >&2 | sed 's/^/  /'
    echo >&2 "Regenerate arm B:  sed 's/^slippage_filter = true/slippage_filter = false/' \\"
    echo >&2 "                       $CFG_A > $CFG_B"
    exit 1
fi
if [ "$nchg" -ne 2 ] || [ "$nslip" -ne 2 ]; then
    echo >&2 "ARMS DIFFER BY MORE THAN slippage_filter — refusing to submit."
    echo >&2 "A confounded A/B is worse than none: it produces a number that looks like the"
    echo >&2 "gate's effect and is not. Settings diff ($CFG_A vs $CFG_B):"
    printf '%s\n' "$d" | sed 's/^/  /' >&2
    echo >&2 "Regenerate arm B:  sed 's/^slippage_filter = true/slippage_filter = false/' \\"
    echo >&2 "                       $CFG_A > $CFG_B"
    exit 1
fi
echo "arms differ by exactly one setting (slippage_filter) — OK"
grep -H '^slippage_filter' "$CFG_A" "$CFG_B" | sed 's/^/  /'

# DOES THE BINARY ACTUALLY HAVE THE GATE? Verifying the configs differ is necessary and NOT
# sufficient, which this A/B learned the expensive way: the slippage gate was never committed,
# the farm built a binary without it, config.rs's unknown-key arm ignored `slippage_filter`
# with a warning nobody read, and both arms of a 180-BAM-scan experiment ran IDENTICAL code
# and reported "the gate cuts 0.0%" across all 90 colonies. The assertion above passed the
# whole time. A config asking for a feature is not evidence the binary implements it.
#
# main.rs prints "slippage filter: ON/OFF" in its banner, so the literal is in the binary iff
# the feature is compiled in. That makes this a free, exact check.
BIN="${DISCOVER_BIN:-rust/peartree-discovery/target/release/peartree-discovery}"
if [ ! -x "$BIN" ]; then
    echo >&2 "no discovery binary at $BIN — build it: bash cluster/build.sh"
    exit 1
fi
# grep the binary DIRECTLY (-a = treat as text). Do NOT write `strings "$BIN" | grep -q ...`:
# under `set -o pipefail` that returns FAILURE on a SUCCESSFUL match — grep -q exits at the
# first hit, strings is still writing, takes SIGPIPE, and pipefail reports the pipeline as
# failed. That is the identical bug that made catalogue_headers.sh label every GRCh38 BAM
# "chr1" (commit 24e0a9b), and it silently inverted this check when first written: the probe
# refused a binary that DID have the gate.
if ! LC_ALL=C grep -qaF 'slippage filter:' "$BIN" 2>/dev/null; then
    echo >&2
    echo >&2 "REFUSING TO SUBMIT: $BIN does not implement the slippage gate."
    echo >&2 "  It would IGNORE slippage_filter (config.rs tolerates unknown keys by design,"
    echo >&2 "  warning to stderr), both arms would run identical code, and the A/B would"
    echo >&2 "  report a 0.0% effect that means 'not tested', not 'no effect'."
    echo >&2 "  Rebuild from current source:  bash cluster/build.sh"
    echo >&2 "  Then confirm:  $BIN --step discover --bam <bam> --out /tmp/x.gz --config $CFG_A 2>&1 | grep 'slippage filter'"
    exit 1
fi
echo "binary implements the slippage gate — OK ($BIN)"

echo
echo "=============== arm A: slippage_filter = true  -> $OUT_A"
FOFN="$FOFN" OUTDIR="$OUT_A" DISCOVER_CFG="$CFG_A" bash cluster/submit_discovery.sh

echo
echo "=============== arm B: slippage_filter = false -> $OUT_B"
FOFN="$FOFN" OUTDIR="$OUT_B" DISCOVER_CFG="$CFG_B" bash cluster/submit_discovery.sh

N="$(wc -l < "$FOFN" | tr -d ' ')"
cat <<EOF

both arms submitted ($N colonies each, ${N}x2 = $((N*2)) BAM scans).

watch:      bjobs -A
when done:  ls $OUT_A/*.txt.gz | wc -l   # expect $N
            ls $OUT_B/*.txt.gz | wc -l   # expect $N
            bash cluster/compare_ab.sh   # arm sizes + the subset check

then, per arm:
  DISCDIR=$OUT_A ... bash cluster/combine_mei.sh     # contract A
  DISCDIR=$OUT_B ... bash cluster/combine_mei.sh     # contract B (the union)
and genotype ONLY contract B — see the header of this script.
EOF
