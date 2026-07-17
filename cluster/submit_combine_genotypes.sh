#!/usr/bin/env bash
# Submit combine_genotypes for the POOLED 9x10 cohort as ONE job (see cluster/combine_gt.sh,
# which it dispatches, for what the cohort gates mean once the cohort is pooled).
#
# Usage:  cluster/submit_combine_genotypes.sh
# Env:
#   FOFN     colony list   (default: $HOME/catalogue/picked/all.bams.fofn)
#   GTDIR    genotype dir  (default: genotypes_grch38)
#   OUT      call table    (default: mei9x10.genotypes.csv.gz)
#   MEM      MB            (default: 16000 -- MEASURED: 1.3 GB peak / 456 s on 174 colonies)
#   CORES/QUEUE/GROUP/VENV/ALLOW_MISSING/WAIT
set -euo pipefail
cd "$(dirname "$0")/.."

FOFN="${FOFN:-$HOME/catalogue/picked/all.bams.fofn}"
GTDIR="${GTDIR:-genotypes_grch38}"
OUT="${OUT:-mei9x10.genotypes.csv.gz}"
MEM="${MEM:-16000}"
CORES="${CORES:-8}"
QUEUE="${QUEUE:-normal}"
VENV="${VENV:-/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE/venv}"
ALLOW_MISSING="${ALLOW_MISSING:-0}"

[ -s "$FOFN" ] || { echo "no fofn: $FOFN" >&2; exit 1; }
[ -x "$VENV/bin/python" ] || { echo "no venv python at $VENV/bin/python — set VENV=" >&2; exit 1; }
# src/config.py is gitignored and environment-specific, and it is where the cohort gates
# live. Without it the job dies inside python after queueing; catch it here instead.
[ -s "src/config.py" ] || {
    echo >&2 "no src/config.py — the cohort gates live there and it is gitignored."
    echo >&2 "  cp cluster/config.py.grch38 src/config.py"
    exit 1; }

NBAM="$(grep -c . "$FOFN")"
NGT="$(ls "$GTDIR"/*.txt.gz 2>/dev/null | wc -l | tr -d ' ')"
echo "  FOFN  = $FOFN  ($NBAM colonies)"
echo "  GTDIR = $GTDIR  ($NGT genotype files present)"
echo "  OUT   = $OUT"
echo "  VENV  = $VENV"
echo "  ALLOW_MISSING = $ALLOW_MISSING"
# combine_gt.sh re-resolves and refuses on its own; this is the early, cheap warning so a
# short cohort is visible BEFORE it sits in the queue.
[ "$NGT" -ge "$NBAM" ] || echo "  NOTE: fewer genotype files than colonies — combine_gt.sh will refuse unless ALLOW_MISSING covers it"

mkdir -p logs
WAIT_ARG=(); [ -n "${WAIT:-}" ] && WAIT_ARG=(-w "$WAIT")
GROUP_ARG=(); [ -n "${GROUP:-}" ] && GROUP_ARG=(-G "$GROUP")

echo "submitting combine_genotypes: $CORES cores, ${MEM}MB, queue=$QUEUE"
[ -n "${WAIT:-}" ] && echo "  waiting on: $WAIT"

# Pass every var explicitly. Env set in the submitting shell does NOT reliably reach the job,
# and the failure is silent: the script falls back to its defaults and does something
# plausible-looking against the wrong inputs. This bit us twice already.
bsub \
    -J "ptcombgt" \
    -o "logs/combgt.out" -e "logs/combgt.err" \
    -n "$CORES" -q "$QUEUE" "${WAIT_ARG[@]}" "${GROUP_ARG[@]}" \
    -R "select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" -M "${MEM}" \
    "FOFN='$FOFN' GTDIR='$GTDIR' OUT='$OUT' THREADS='$CORES' VENV='$VENV' ALLOW_MISSING='$ALLOW_MISSING' bash cluster/combine_gt.sh"

echo "when done: zcat $OUT | head -2 ; zcat $OUT | tail -n +2 | wc -l   # loci passing the gates"
