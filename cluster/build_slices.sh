#!/usr/bin/env bash
# ============================================================================
# Phase 0 of the fast threshold-sweep harness: build minimal region-slice BAMs
# for the 9x10 benchmark, plus the exact per-BAM genome coverage median.
#
# WHY: discovery re-reads every full BAM twice (coverage pre-pass + extract). On
# the 9x10 cohort a full arm takes hours. But the frozen truth set only scores
# 6,896 loci (TP+FP) inside ~1.1% of the genome. A read-start region slice keeps
# FULL depth at those loci and drops everything else, so discovery on the slice
# runs in seconds -> dozens of threshold configs per hour.
#
# THE TRAP (handled): the genome coverage median feeds the local/median ratio
# gates (SPEC-3, discordant_coverage_max_mult). Estimated from a slice it is
# inflated and the gates silently stop firing. So we ALSO precompute the true
# full-BAM median here and pin it at sweep time via MEDIAN_DIR/PEARTREE_COVERAGE_MEDIAN.
#
# One-time cost: one full-BAM pass per colony for the median (index-based slice is
# cheap). Pay once, sweep forever. Run on a farm22 head node.
# ============================================================================
set -euo pipefail

REPO=/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE-pd44579
FOFN=/nfs/users/nfs_j/jd43/catalogue/picked/all.bams.fofn      # 90 full staged BAMs
TRUTH="$REPO/analysis/mei9x10/frozen_truth.tsv"                # shipped (untracked)
BED="$REPO/analysis/mei9x10/truth_regions.bed"
SLICEDIR="$REPO/slices"
MEDDIR="$SLICEDIR/medians"
SLICEFOFN="$SLICEDIR/slices.fofn"
BIN="$REPO/rust/peartree-discovery/target/release/peartree-discovery"
PAD=1000
THROTTLE=40
MEM=8000
SAMTOOLS_MODULE=samtools-1.19/python-3.12.0

cd "$REPO"
[ -s "$FOFN"  ] || { echo "no fofn: $FOFN" >&2; exit 1; }
[ -s "$TRUTH" ] || { echo "no truth set: $TRUTH  (tsh scp it to the farm first)" >&2; exit 1; }
[ -x "$BIN"   ] || { echo "no binary: $BIN  (module load rust/1.87.0 && bash cluster/build.sh)" >&2; exit 1; }
# GUARD (hpc-sanger lesson): a stale binary silently ignores coverage_median_override and
# the whole median-pin is a no-op. Refuse to run unless the rebuilt binary knows the step.
# Probe with an INVALID step: the new binary lists its supported steps (incl. coverage-median)
# in the error, whereas a valid step with no --bam only prints the generic usage line. The
# probe deliberately EXITS NON-ZERO (error path), so capture with `|| true` and match the
# string — a bare `$BIN ... | grep` trips `set -o pipefail` regardless of whether grep matches.
probe_out="$("$BIN" --step __probe__ 2>&1 || true)"
case "$probe_out" in
  *coverage-median*) : ;;
  *) echo "FATAL: $BIN does not implement --step coverage-median -> REBUILD (module load rust/1.87.0 && bash cluster/build.sh)" >&2; exit 1 ;;
esac

mkdir -p "$SLICEDIR" "$MEDDIR" logs/slices

# 1) BED: pad each scored/excluded locus by +/-PAD, clamp, sort, merge. Includes the
#    11,300 'exclude' loci too so combine/merge behaviour around each scored locus is
#    faithful. ~17k intervals, ~1.1% of the genome.
echo "=== building region BED from $TRUTH ==="
awk -F'\t' -v P="$PAD" 'NR>1{
    n=split($1,a,":"); split(a[2],b,"-");
    s=b[1]-P; if(s<0)s=0; print a[1]"\t"s"\t"b[2]+P
  }' "$TRUTH" \
 | LC_ALL=C sort -k1,1 -k2,2n \
 | awk 'BEGIN{OFS="\t"}{ if($1==c && $2<=e){ if($3>e)e=$3; next } if(NR>1)print c,s,e; c=$1;s=$2;e=$3 } END{if(c!="")print c,s,e}' \
 > "$BED"
echo "  intervals: $(wc -l < "$BED")   bp: $(awk '{t+=$3-$2}END{print t}' "$BED")"

# 2) deterministic slice fofn (derived from the full fofn, colony-id order preserved)
awk -v D="$SLICEDIR" '{
    n=split($0,p,"/"); id=p[n];
    sub(/\.sample\.dupmarked\.bam$/,"",id); sub(/\.bam$/,"",id); sub(/\.cram$/,"",id);
    print D"/"id".slice.bam"
  }' "$FOFN" > "$SLICEFOFN"
N=$(wc -l < "$FOFN"); N=${N// /}
echo "=== $N colonies -> slices in $SLICEDIR, medians in $MEDDIR ==="

# 3) array: per colony, build the region slice (index-based) + the exact genome median.
#    Per-BAM median FILES (never a shared file from parallel tasks — Lustre rule).
bsub -J "ptslice[1-$N]%$THROTTLE" \
     -n 4 -M "$MEM" -R "select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -q normal \
     -o "logs/slices/slice.%I.log" -e "logs/slices/slice.%I.err" \
     "bash '$REPO/cluster/slice_one.sh' \$LSB_JOBINDEX '$FOFN' '$BED' '$SLICEDIR' '$MEDDIR' '$BIN' '$SAMTOOLS_MODULE'"

echo
echo "watch:  bjobs -A ; grep -lE 'Exited|TERM_' logs/slices/slice.*.err"
echo "when DONE, sanity-check counts (expect $N each):"
echo "  ls $SLICEDIR/*.slice.bam | wc -l ; ls $MEDDIR/*.median | wc -l"
echo "then run a config with:  MEDIAN_DIR=$MEDDIR FOFN=$SLICEFOFN OUTDIR=... bash cluster/submit_discovery.sh"
