#!/usr/bin/env bash
# Compare the two discovery arms of the slippage-gate A/B.
#
#   bash cluster/compare_ab.sh
#
# Env: OUT_A (slippage_filter=true), OUT_B (=false), FOFN, NCHECK (subset spot-check size)
#
# Reports per colony and in total: LOCI in each arm, what the gate cut, and whether arm A is
# genuinely a SUBSET of arm B.
#
# THE OUTPUT FORMAT, because getting this wrong silently inflates every number ~10x.
# Discovery writes 4-line FASTQ-like records, NOT a TSV:
#     @<contig>:<left>-<right>:<SIDE>:<TYPE>
#     <sequence>
#     +
#     <quality>
# SIDE is LEFT|RIGHT; TYPE is CLIPPED | ALIGNED | CLIPPED_POLYA | MATE1 | MATE2. Coordinates
# may be `polyA_<n>` rather than a plain integer (poly-A-anchored end).
#
# A LOCUS IS `contig:left-right` AND IT SPANS MANY RECORDS -- one per side per evidence type.
# Measured on a real colony: 225,300 lines = 56,325 records = 6,532 loci. So:
#   - counting LINES     overstates by ~34x
#   - counting RECORDS   overstates by ~8.6x
#   - `grep -c '^@'`     counts records, not loci -- AND is unsafe in principle, because a
#                        FASTQ quality line may legitimately begin with '@' (it happens to be
#                        0 in the file checked, which is luck, not a guarantee).
# Take every 4th line and strip the trailing :SIDE:TYPE. Never grep for '^@'.
set -euo pipefail
cd "$(dirname "$0")/.."

OUT_A="${OUT_A:-discovery_grch38_slip}"
OUT_B="${OUT_B:-discovery_grch38_noslip}"
FOFN="${FOFN:-$HOME/catalogue/picked/all.bams.fofn}"
NCHECK="${NCHECK:-5}"

[ -d "$OUT_A" ] || { echo "no arm A dir: $OUT_A" >&2; exit 1; }
[ -d "$OUT_B" ] || { echo "no arm B dir: $OUT_B" >&2; exit 1; }

# gunzip -c, never zcat: on macOS zcat wants .Z and fails, and a `2>/dev/null` on it once
# turned six read errors into a tidy column of zeros that read as real data.
# Locus key: 4th line -> drop '@' -> join all but the last two ':' fields (SIDE, TYPE).
# Contig names have no ':', but rebuilding from the front rather than cutting a fixed field
# keeps this correct even if one ever does.
loci() { gunzip -c "$1" | awk 'NR%4==1{sub(/^@/,""); n=split($0,p,":"); if(n<3) next; k=p[1]; for(i=2;i<=n-2;i++) k=k":"p[i]; print k}' | sort -u; }

na=$(ls "$OUT_A"/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')
nb=$(ls "$OUT_B"/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')
nf=$([ -s "$FOFN" ] && wc -l < "$FOFN" | tr -d ' ' || echo '?')
echo "arm A ($OUT_A): $na colonies"
echo "arm B ($OUT_B): $nb colonies"
echo "fofn:           $nf colonies"
if [ "$na" != "$nb" ]; then
    echo >&2 "ARM SIZES DIFFER ($na vs $nb) — the arms are not comparable yet."
    echo >&2 "Wait for both arrays, or re-run the missing tasks. Comparing now would attribute"
    echo >&2 "a scheduling artefact to the gate."
    exit 1
fi
[ "$na" -gt 0 ] || { echo "no discovery output in either arm" >&2; exit 1; }

echo
printf '%-28s %10s %10s %10s %8s\n' colony loci_A loci_B cut cut_pct
tot_a=0; tot_b=0; nonsub=0; checked=0
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
for fa in "$OUT_A"/*.txt.gz; do
    id="$(basename "$fa" .txt.gz)"
    fb="$OUT_B/$id.txt.gz"
    [ -s "$fb" ] || { echo "  MISSING in arm B: $id" >&2; continue; }
    loci "$fa" > "$TMP/a"; loci "$fb" > "$TMP/b"
    a=$(wc -l < "$TMP/a" | tr -d ' '); b=$(wc -l < "$TMP/b" | tr -d ' ')
    cut=$(( b - a ))
    pct=$(awk -v c="$cut" -v b="$b" 'BEGIN{printf "%.1f", (b>0 ? 100*c/b : 0)}')
    printf '%-28s %10d %10d %10d %7s%%\n' "$id" "$a" "$b" "$cut" "$pct"
    tot_a=$((tot_a + a)); tot_b=$((tot_b + b))

    # Subset spot-check. The gate only removes clip clusters, so arm A SHOULD be contained in
    # arm B -- but dropping evidence can change a cluster's consensus and shift a breakpoint,
    # which would put an arm-A locus at coordinates absent from arm B. Measure; don't assume.
    if [ "$checked" -lt "$NCHECK" ]; then
        only_a=$(comm -23 "$TMP/a" "$TMP/b" | wc -l | tr -d ' ')
        if [ "$only_a" -gt 0 ]; then
            echo "    NOT A SUBSET: $only_a loci in arm A absent from arm B ($id)" >&2
            comm -23 "$TMP/a" "$TMP/b" | head -3 | sed 's/^/      e.g. /' >&2
            nonsub=$((nonsub + 1))
        fi
        checked=$((checked + 1))
    fi
done

echo
cut_tot=$(( tot_b - tot_a ))
echo "TOTAL   arm A: $tot_a loci   arm B: $tot_b loci   cut: $cut_tot ($(awk -v c="$cut_tot" -v b="$tot_b" 'BEGIN{printf "%.1f", (b>0 ? 100*c/b : 0)}')%)"
echo
echo "For reference the gate cut 30% of the PD44579 contract — ONE donor, hs37d5, and measured"
echo "on the CONTRACT (post-combine), not on per-colony discovery loci. It is an expectation to"
echo "test, not a target to hit, and not the same quantity as the number above."
echo
if [ "$nonsub" -gt 0 ]; then
    echo "SUBSET CHECK FAILED on $nonsub/$checked colonies sampled."
    echo "  => Do NOT genotype only arm B and subset: the arms have diverged in coordinates."
    echo "     Build a contract per arm and genotype each separately."
else
    echo "SUBSET CHECK PASSED on $checked of $na colonies (arm A ⊆ arm B)."
    echo "  => Genotype arm B's contract ONCE, then evaluate arm A by subsetting those calls to"
    echo "     arm A's loci: same calls, different locus set, and it halves the costliest step."
    echo "  NB SAMPLED, NOT PROVEN — $checked of $na colonies. Raise with NCHECK=$na to check all."
fi
