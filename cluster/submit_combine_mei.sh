#!/usr/bin/env bash
# Submit combine_insertions for the 10x10 MEI benchmark as ONE job over ONE POOLED cohort
# of all 50 colonies (see cluster/combine_mei.sh for why pooled, not per-patient).
#
# Usage:  cluster/submit_combine_mei.sh
# Env:
#   MEM       MB     (default: 32000)
#   CORES     count  (default: 16)
#   QUEUE     LSF queue (default: normal)
#   WAIT      job dependency, e.g. "ended(818591)"
set -euo pipefail
cd "$(dirname "$0")/.."

PATIENTS="${PATIENTS:-PD34200 PD37449 PD43947 PD41048 PD43974}"
# PD44579 peaked at 15.5 GB combining 174 discovery files; 50 files is well under that,
# but combine is the memory-hungry step and discovery already surprised us once (a task
# died at 16.28 GB against a 16 GB reservation). This is a single job — headroom is cheap.
MEM="${MEM:-32000}"
CORES="${CORES:-16}"
QUEUE="${QUEUE:-normal}"
FOFNDIR="${FOFNDIR:-$HOME/mei10x10}"
DISCDIR="${DISCDIR:-discovery_mei10x10}"
OUTDIR="${OUTDIR:-insertions_mei10x10}"
STEM="${STEM:-mei10x10}"
ALLOW_MISSING="${ALLOW_MISSING:-1}"
# The venv is NOT in this run dir — it lives in the sibling PEAR-TREE checkout.
VENV="${VENV:-/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE/venv}"

mkdir -p logs "$OUTDIR"
[ -x "$VENV/bin/python" ] || { echo "no venv python at $VENV/bin/python — set VENV=" >&2; exit 1; }
[ -s "$FOFNDIR/all.bams.fofn" ] || { echo "missing pooled fofn: $FOFNDIR/all.bams.fofn" >&2; exit 1; }
for p in $PATIENTS; do
    [ -s "$FOFNDIR/$p.bams.fofn" ] || { echo "missing fofn: $FOFNDIR/$p.bams.fofn" >&2; exit 1; }
done

NBAM="$(grep -c . "$FOFNDIR/all.bams.fofn")"
# NB: `|| NDISC=0` is load-bearing. Under `set -euo pipefail`, when DISCDIR has no
# *.txt.gz yet (the normal case: this job is submitted with WAIT=ended(discovery) while
# discovery is STILL RUNNING), the glob matches nothing, `ls` exits non-zero, pipefail
# propagates it, and set -e kills the script BEFORE bsub — silently, with no output. The
# `||` catches that so the count is just 0 and submission proceeds. Without it, dependency-
# chained submission is impossible.
NDISC="$(ls "$DISCDIR"/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')" || NDISC=0
echo "  PATIENTS = $PATIENTS"
echo "  FOFN     = $FOFNDIR/all.bams.fofn  ($NBAM colonies)"
echo "  DISCDIR  = $DISCDIR  ($NDISC discovery files present)"
echo "  OUTDIR   = $OUTDIR   (contract stem: $STEM)"
echo "  VENV     = $VENV"
echo "  ALLOW_MISSING = $ALLOW_MISSING"

WAIT_ARG=(); [ -n "${WAIT:-}" ] && WAIT_ARG=(-w "$WAIT")
GROUP_ARG=(); [ -n "${GROUP:-}" ] && GROUP_ARG=(-G "$GROUP")

echo "submitting pooled combine: ${CORES} cores, ${MEM}MB, queue=$QUEUE"
[ -n "${WAIT:-}" ] && echo "  waiting on: $WAIT"

# Pass every var explicitly. Env set in the submitting shell does NOT reliably reach the
# job, and the failure is silent: the script falls back to its defaults and does something
# plausible-looking against the wrong inputs. This bit us twice already.
# Namespace the job name (and logs) by the run's OUTDIR. A FIXED name like "ptcomb" is a trap
# for downstream `WAIT=ended(ptcomb)`: LSF matches the name against ALL the user's jobs, so an
# already-ended ptcomb from an earlier run satisfies the condition instantly and the dependent
# genotype array fires early against a contract that does not exist yet. A per-run name keeps
# `ended(<name>)` chaining safe as long as each run uses a fresh OUTDIR; for extra safety,
# prefer depending on the numeric job id this script prints.
CJTAG="$(basename "$OUTDIR")"
bsub \
    -J "ptcomb_$CJTAG" \
    -o "logs/comb.$CJTAG.out" -e "logs/comb.$CJTAG.err" \
    -n "$CORES" -q "$QUEUE" "${WAIT_ARG[@]}" "${GROUP_ARG[@]}" \
    -R "select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" -M "${MEM}" \
    "PATIENTS='$PATIENTS' FOFNDIR='$FOFNDIR' DISCDIR='$DISCDIR' OUTDIR='$OUTDIR' STEM='$STEM' THREADS='$CORES' VENV='$VENV' ALLOW_MISSING='$ALLOW_MISSING' bash cluster/combine_mei.sh"

echo "when done: zcat $OUTDIR/$STEM.genotyping.txt.gz | grep -c '^>'"
