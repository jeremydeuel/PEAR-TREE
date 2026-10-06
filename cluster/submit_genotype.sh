#!/usr/bin/env bash
# Submit genotyping as an LSF job array on farm22: one task per colony, 1 core each,
# every colony against the SAME POOLED contract (all 9 donors in one locus set).
#
# Usage:  cluster/submit_genotype.sh
# Env (all optional):
#   FOFN      file list         (default: $HOME/catalogue/picked/all.bams.fofn)
#   OUTDIR    output dir        (default: genotypes_grch38)
#   CONTRACT  pooled contract   (default: insertions_grch38/mei9x10.genotyping.txt.gz)
#   GENO_CFG  rust config       (default: cluster/config.genotype.grch38)
#   THROTTLE  max concurrent    (default: 50)
#   MEM       MB per task       (default: 2000 -- MEASURED: the Rust genotyper peaked at
#                                140 MB on PD44579, so this is ~14x headroom already. Do not
#                                inflate it: rusage[mem] RESERVES, so a bigger number just
#                                throttles how many tasks the scheduler co-locates.)
#   QUEUE     LSF queue         (default: normal; ~17 min/colony single-threaded)
#   GROUP     LSF fairshare group (-G), if your setup requires one
#
# WHY ONLY ONE ARM IS GENOTYPED. Arm B (slippage gate OFF) is a SUPERSET of arm A: the gate
# only ever removes loci. So genotype arm B's contract ONCE and evaluate arm A by subsetting
# those calls to arm A's loci -- identical calls, different locus set, and it halves the
# costliest step (90 colonies x ~30k loci). This is only valid if arm A is genuinely
# contained in arm B; cluster/compare_ab.sh tests that rather than assuming it.
set -euo pipefail
cd "$(dirname "$0")/.."

FOFN="${FOFN:-$HOME/catalogue/picked/all.bams.fofn}"
OUTDIR="${OUTDIR:-genotypes_grch38}"
CONTRACT="${CONTRACT:-insertions_grch38/mei9x10.genotyping.txt.gz}"
GENO_CFG="${GENO_CFG:-cluster/config.genotype.grch38}"
THROTTLE="${THROTTLE:-50}"
MEM="${MEM:-2000}"
QUEUE="${QUEUE:-normal}"
BIN="${GENOTYPE_BIN:-rust/peartree-genotype/target/release/peartree-genotype}"
# GENOTYPE_IMPL=v2 (default) runs peartree-genotype2 (see genotype_one.sh for COMBINED / GENOME_2BIT)
if [ "${GENOTYPE_IMPL:-v2}" = v2 ]; then
    BIN="${GENOTYPE2_BIN:-rust/peartree-genotype2/target/release/peartree-genotype2}"
    GENO_CFG="${GENO2_CFG:-cluster/config.genotype2.grch38}"
fi

[ -s "$FOFN" ]     || { echo "no such fofn: $FOFN (run cluster/stage_picked.sh --fofn)" >&2; exit 1; }
[ -s "$GENO_CFG" ] || { echo "no config: $GENO_CFG" >&2; exit 1; }
[ -x "$BIN" ]      || { echo "no genotype binary at $BIN — build it: bash cluster/build.sh" >&2; exit 1; }
# The contract is produced by combine_insertions. In the chained case this genotype array is
# submitted with WAIT=ended(<combine>) BEFORE the contract exists, so a hard existence check
# here makes dependency-chained submission impossible (the same silent-blocker class as the
# combine submitters). Defer the contract-dependent checks (loci count + contig preflight) to
# run time when a WAIT is set; hard-fail only for a standalone submission, where the contract
# genuinely should already be on disk.
CONTRACT_DEFERRED=""
if [ ! -s "$CONTRACT" ]; then
    if [ -n "${WAIT:-}" ]; then
        CONTRACT_DEFERRED=1
        echo >&2 "note: $CONTRACT not present yet — deferring loci + contig preflight to run time"
        echo >&2 "      (queuing behind: $WAIT)"
    else
        echo >&2 "no contract: $CONTRACT"
        echo >&2 "combine_insertions has not produced one yet. Either wait for it, or submit with"
        echo >&2 "  WAIT=\"ended(<combine job>)\"  to queue this array behind combine."
        exit 1
    fi
fi

N="$(wc -l < "$FOFN" | tr -d ' ')"
[ "$N" -gt 0 ] || { echo "empty $FOFN" >&2; exit 1; }
mkdir -p logs "$OUTDIR"

if [ -n "$CONTRACT_DEFERRED" ]; then
    NLOCI="(pending combine)"
else
    NLOCI="$(zcat "$CONTRACT" | grep -c '^>' || true)"
    [ "$NLOCI" -gt 0 ] || {
        echo >&2 "REFUSING TO SUBMIT: the contract $CONTRACT has 0 loci."
        echo >&2 "Every task would genotype nothing and exit 0. Check combine_insertions' log."
        exit 1; }
fi

# Say out loud what we are about to do -- the defaults are a live hazard, not a convenience:
# if the caller's env does not reach this script it submits a plausible array against the
# wrong inputs, and that failure is silent.
echo "  FOFN     = $FOFN  ($N colonies)"
echo "  CONTRACT = $CONTRACT  ($NLOCI loci, pooled over all donors)"
echo "  OUTDIR   = $OUTDIR  ($(ls "$OUTDIR" 2>/dev/null | wc -l | tr -d ' ') files already present)"
echo "  CONFIG   = $GENO_CFG"

# PREFLIGHT: do the CONTRACT's contigs exist in these BAMs?
#
# This is the genotyping twin of submit_discovery.sh's allowlist check, and the same silent
# failure: contract locus names are "<contig>:<start>-<end>", so a contract built in hs37d5
# space (1..Y) against GRCh38 BAMs (chr1..chrY) matches NOTHING. Every task runs to
# completion, exits 0, and calls every locus NA -- indistinguishable from "no evidence".
# 90 colonies x ~17 min of farm time, and it surfaces later as an all-NA call table.
FIRST_BAM="$(head -1 "$FOFN")"
if [ -n "$CONTRACT_DEFERRED" ]; then
    echo "  contigs  : NOT CHECKED (contract not built yet — deferred with the WAIT dependency)" >&2
elif [ -s "$FIRST_BAM" ]; then
    module load samtools-1.19/python-3.12.0 2>/dev/null || true
    if command -v samtools >/dev/null 2>&1; then
        # Read the contract's contigs and the BAM's @SQ, then intersect. Header read only.
        ctg="$(zcat "$CONTRACT" | awk '/^>/{sub(/^>/,""); sub(/:[^:]*$/,""); print}' | sort -u)"
        hdr_sq="$(samtools view -H "$FIRST_BAM" 2>/dev/null | awk '$1=="@SQ"{for(i=2;i<=NF;i++) if($i ~ /^SN:/){sub(/^SN:/,"",$i); print $i}}')"
        n_ctg=$(printf '%s\n' "$ctg" | grep -c . || true)
        n_hit=$(printf '%s\n' "$ctg" | grep -Fxc -f <(printf '%s\n' "$hdr_sq") - 2>/dev/null || true)
        n_hit=${n_hit:-0}
        echo "  contigs  : $n_hit/$n_ctg contract contigs present in $(basename "$FIRST_BAM")"
        if [ "$n_hit" -eq 0 ]; then
            echo >&2
            echo >&2 "REFUSING TO SUBMIT: none of the contract's $n_ctg contigs exist in the BAM."
            echo >&2 "  contract : $CONTRACT"
            echo >&2 "  contigs  : $(printf '%s\n' "$ctg" | head -4 | paste -sd, -)..."
            echo >&2 "  BAM @SQ  : $(printf '%s\n' "$hdr_sq" | head -4 | paste -sd, -)..."
            echo >&2 "Every locus would be NA and every task would exit 0."
            echo >&2 "The contract was almost certainly built from the wrong cohort's discovery files."
            exit 1
        fi
    else
        echo "  contigs  : NOT CHECKED (samtools not on PATH)" >&2
    fi
fi

GROUP_ARG=()
[ -n "${GROUP:-}" ] && GROUP_ARG=(-G "$GROUP")

# Apply the WAIT as an actual LSF dependency. This was the bug: WAIT was used ONLY to defer the
# contract preflight above (so the array could be submitted before combine built the contract),
# but it was never passed to bsub -- so a chained genotype array had NO dependency, ran
# immediately against a not-yet-existent contract, and every task hit `[ -s "$CONTRACT" ]` in
# genotype_one.sh and EXITED. (submit_combine_genotypes.sh applies WAIT; this one forgot to.)
# Deferring the check without gating the job is the worst of both: it removes the guard AND
# does not wait. Both must move together.
WAIT_ARG=()
[ -n "${WAIT:-}" ] && WAIT_ARG=(-w "$WAIT")

# Namespace the logs and the job name by OUTDIR: two arrays running at once (as the A/B arms
# did) otherwise share logs/gt.<i>.out, and LSF's -o APPENDS -- so the files end up holding
# both runs interleaved and no failure can be attributed to a run.
TAG="$(basename "$OUTDIR")"
LOGDIR="logs/$TAG"
mkdir -p "$LOGDIR"
echo "  LOGS     = $LOGDIR/gt.<task>.{out,err}"

echo "submitting genotype array: $N colonies x $NLOCI loci, <=$THROTTLE concurrent, ${MEM}MB, queue=$QUEUE"
[ -n "${WAIT:-}" ] && echo "  waiting on: $WAIT"
bsub \
    -J "ptgt_${TAG}[1-${N}]%${THROTTLE}" \
    -o "$LOGDIR/gt.%I.out" -e "$LOGDIR/gt.%I.err" \
    -n 1 -q "$QUEUE" "${GROUP_ARG[@]}" "${WAIT_ARG[@]}" \
    -R "select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" -M "${MEM}" \
    "FOFN='$FOFN' OUTDIR='$OUTDIR' CONTRACT='$CONTRACT' GENO_CFG='$GENO_CFG' bash cluster/genotype_one.sh \$LSB_JOBINDEX"

echo "watch with: bjobs -A ; tail -f $LOGDIR/gt.1.out"
echo "when done: ls $OUTDIR/*.txt.gz | wc -l   # expect $N"
echo "failures:  grep -lE 'TERM_MEMLIMIT|Exited with exit code' $LOGDIR/gt.*.err"
echo "then:      bash cluster/submit_combine_genotypes.sh"
