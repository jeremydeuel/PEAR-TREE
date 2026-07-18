#!/usr/bin/env bash
# submit_all.sh -- fire off the whole 4-stage PEAR-TREE pipeline for ONE run, correctly
# sequenced with LSF dependencies, in a single command. Resumable and idempotent: it checks
# what already exists and submits ONLY the stages that are missing, so re-running after a
# mid-flight failure never re-does completed work and never duplicates an in-flight array.
#
# THE FOUR STAGES (each feeds the next):
#   1. discovery          submit_discovery.sh        -> $DISCDIR/<id>.txt.gz   (one per colony)
#   2. combine_insertions combine_mei.sh (bsub'd)    -> $CONTRACT              (pooled contract)
#   3. genotype           submit_genotype.sh         -> $GTDIR/<id>.txt.gz     (one per colony)
#   4. combine_genotypes  submit_combine_genotypes.sh-> $FINAL_OUT             (final call table)
#
# ---------------------------------------------------------------------------------------------
# USAGE
#   TAG=<label> DISCOVER_CFG=<cfg> bash cluster/submit_all.sh
#
# The only two things it will NOT guess (you MUST supply them):
#   TAG           run label, e.g. noSPEC8  -> names every output dir (discovery_grch38_$TAG, ...)
#   DISCOVER_CFG  discovery config, e.g. cluster/config.discovery.grch38.noSPEC8
#                 (there is no safe default -- the wrong one silently writes empty output.)
#
# Everything else is INFERRED and printed; override any of it by exporting it:
#   FOFN          colony list        (default: /nfs/users/nfs_j/jd43/catalogue/picked/all.bams.fofn)
#   GENO_CFG      genotype config    (default: cluster/config.genotype.grch38)
#   STEM          contract basename  (default: mei9x10)
#   VENV          python venv        (default: .../PEAR-TREE/venv)
#   DISCDIR/INS_OUTDIR/GTDIR/FINAL_OUT   output dirs/files (default: *_grch38_$TAG)
#   THROTTLE QUEUE GROUP  ALLOW_MISSING
#   DISC_MEM GT_MEM COMBINE_MEM COMBGT_MEM COMBINE_CORES COMBGT_CORES
#
# ATTACH TO AN ALREADY-RUNNING STAGE instead of resubmitting it (e.g. discovery already in the
# queue): pass its numeric job id. The next stage is then chained onto that id.
#   DISC_JOBID=957812        (discovery already running)
#   COMBINE_JOBID=... GT_JOBID=...
# The script also auto-detects a running stage by job name, but an explicit *_JOBID is safer.
#
#   DRYRUN=1      print the resolved plan (what would be submitted, and the dependency chain)
#                 and exit WITHOUT submitting anything. Recommended for a first look.
#   SKIP_ASSEMBLY_CHECK=1   skip the install.sh check-config genome_2bit/chain validation.
#
# WHY chain on the numeric job id, not the job name: LSF's ended(<name>) matches ALL your jobs
# of that name, so a same-named job that already ended satisfies the dependency instantly and
# the dependent fires early against inputs that do not exist yet (this bug has bitten this repo
# -- see submit_combine_genotypes.sh's header). This script captures each bsub's numeric id and
# chains ended(<id>), which is unambiguous.
# ---------------------------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")/.."
PT_ROOT="$PWD"

# ---------- help:  -h | --help | HELP=1 ----------
usage() {
cat <<'HELPDOC'
submit_all.sh -- submit the full 4-stage PEAR-TREE pipeline for ONE run, correctly sequenced
with LSF dependencies, in a single command. Resumable and idempotent: only missing stages are
submitted, so re-running after a failure never re-does completed work or duplicates a job.

STAGES (each feeds the next):
  1 discovery          -> $DISCDIR/<id>.txt.gz        (one per colony)
  2 combine_insertions -> $CONTRACT                   (pooled contract)
  3 genotype           -> $GTDIR/<id>.txt.gz          (one per colony)
  4 combine_genotypes  -> $FINAL_OUT                  (final call table)

USAGE
  TAG=<label> DISCOVER_CFG=<cfg> bash cluster/submit_all.sh

REQUIRED (no defaults -- the script refuses to guess these):
  TAG           run label, e.g. noSPEC8 -> names every output dir (discovery_grch38_$TAG, ...)
  DISCOVER_CFG  discovery config, e.g. cluster/config.discovery.grch38.noSPEC8
                (the wrong assembly's config silently writes empty output -> no default.)

INFERRED (printed at startup; override by exporting):
  FOFN          colony list      (default /nfs/users/nfs_j/jd43/catalogue/picked/all.bams.fofn)
  GENO_CFG      genotype config  (default cluster/config.genotype.grch38)
  STEM          contract stem    (default mei9x10)
  VENV          python venv      (default .../PEAR-TREE/venv)
  DISCDIR INS_OUTDIR GTDIR FINAL_OUT   output dirs/file  (default *_grch38_$TAG)
  PATIENTS      donor list       (computed from the <donor>.bams.fofn files next to FOFN)

TUNABLES:
  THROTTLE(50) QUEUE(normal) GROUP  ALLOW_MISSING(1)
  DISC_MEM(8000) GT_MEM(2000) COMBINE_MEM(32000) COMBGT_MEM(16000)
  COMBINE_CORES(16) COMBGT_CORES(8)  DISCOVER_BIN GENOTYPE_BIN

ATTACH to an already-running stage (do not resubmit it); the next stage chains onto its id:
  DISC_JOBID=957812   COMBINE_JOBID=...   GT_JOBID=...

MODES:
  DRYRUN=1              print the resolved plan + dependency chain, submit nothing.
  SKIP_ASSEMBLY_CHECK=1 skip the install.sh check-config genome_2bit/chain validation.
  -h | --help | HELP=1  show this help.

EXAMPLES:
  # see the plan without submitting
  DRYRUN=1 TAG=noSPEC8 DISCOVER_CFG=cluster/config.discovery.grch38.noSPEC8 bash cluster/submit_all.sh
  # discovery already queued as 957812 -> attach and chain stages 2-4 behind it
  TAG=noSPEC8 DISCOVER_CFG=cluster/config.discovery.grch38.noSPEC8 DISC_JOBID=957812 bash cluster/submit_all.sh
  # resume after a partial failure: just run the same command again -- done stages are skipped.

Full manual: cluster/SUBMIT_ALL.md
HELPDOC
}
case "${1:-}" in -h|--help|help) usage; exit 0 ;; esac
[ -n "${HELP:-}" ] && { usage; exit 0; }

# ---------- required run identity (refuse to guess) ----------
TAG="${TAG:-}"
DISCOVER_CFG="${DISCOVER_CFG:-}"

# ---------- inferred / overridable ----------
FOFN="${FOFN:-/nfs/users/nfs_j/jd43/catalogue/picked/all.bams.fofn}"
GENO_CFG="${GENO_CFG:-cluster/config.genotype.grch38}"
STEM="${STEM:-mei9x10}"
VENV="${VENV:-/lustre/scratch126/casm/teams/team273/users/jd43/PEAR-TREE/venv}"
THROTTLE="${THROTTLE:-50}"
QUEUE="${QUEUE:-normal}"
GROUP="${GROUP:-}"
ALLOW_MISSING="${ALLOW_MISSING:-1}"
DISC_MEM="${DISC_MEM:-8000}"
GT_MEM="${GT_MEM:-2000}"
COMBINE_MEM="${COMBINE_MEM:-32000}"
COMBGT_MEM="${COMBGT_MEM:-16000}"
COMBINE_CORES="${COMBINE_CORES:-16}"
COMBGT_CORES="${COMBGT_CORES:-8}"
DISCOVER_BIN="${DISCOVER_BIN:-rust/peartree-discovery/target/release/peartree-discovery}"
GENOTYPE_BIN="${GENOTYPE_BIN:-rust/peartree-genotype/target/release/peartree-genotype}"
DISC_JOBID="${DISC_JOBID:-}"
COMBINE_JOBID="${COMBINE_JOBID:-}"
GT_JOBID="${GT_JOBID:-}"
DRYRUN="${DRYRUN:-}"
SKIP_ASSEMBLY_CHECK="${SKIP_ASSEMBLY_CHECK:-}"

# ---------- derived (deterministic; printed below) ----------
DISCDIR="${DISCDIR:-discovery_grch38_$TAG}"
INS_OUTDIR="${INS_OUTDIR:-insertions_grch38_$TAG}"
CONTRACT="$INS_OUTDIR/$STEM.genotyping.txt.gz"
GTDIR="${GTDIR:-genotypes_grch38_$TAG}"
FINAL_OUT="${FINAL_OUT:-$STEM.$TAG.genotypes.csv.gz}"
FOFNDIR="$(dirname "$FOFN")"

say()  { printf '%s\n' "$*"; }
die()  { printf '\n\033[1mCannot start -- something is missing:\033[0m\n%s\n' "$*" >&2; exit 1; }

# ============================================================================================
# 0. required-variable gate (kind, specific)
# ============================================================================================
problems=()
[ -n "$TAG" ]          || problems+=("TAG is not set. Give the run a label, e.g. TAG=noSPEC8 -- it names every output dir.")
[ -n "$DISCOVER_CFG" ] || problems+=("DISCOVER_CFG is not set. Point it at the discovery config for this run, e.g.
      DISCOVER_CFG=cluster/config.discovery.grch38.noSPEC8
      (there is no default: the wrong assembly's config silently writes empty output.)")
if [ ${#problems[@]} -gt 0 ]; then
    msg=""
    for p in "${problems[@]}"; do msg+="  - $p"$'\n'; done
    die "$msg"$'\n'"Then re-run:  TAG=<label> DISCOVER_CFG=<cfg> bash cluster/submit_all.sh"
fi

# ============================================================================================
# 1. resolve PATIENTS from the per-patient fofns (computed, never guessed)
#    combine_insertions needs one BAM per donor for its assembly gate; the donor set is exactly
#    the *.bams.fofn files sitting next to all.bams.fofn.
# ============================================================================================
PATIENTS=""
if [ -d "$FOFNDIR" ]; then
    for f in "$FOFNDIR"/*.bams.fofn; do
        [ -e "$f" ] || continue
        b="$(basename "$f" .bams.fofn)"
        [ "$b" = "all" ] && continue
        PATIENTS+="$b "
    done
fi
PATIENTS="${PATIENTS% }"

# ============================================================================================
# 2. resolved plan -- say everything out loud before touching the farm
# ============================================================================================
say "======================================================================"
say " PEAR-TREE full pipeline -- run TAG=$TAG"
say "======================================================================"
say "  FOFN         = $FOFN"
say "  FOFNDIR      = $FOFNDIR"
say "  PATIENTS     = ${PATIENTS:-<none found!>}"
say "  DISCOVER_CFG = $DISCOVER_CFG"
say "  GENO_CFG     = $GENO_CFG"
say "  discovery    -> $DISCDIR/"
say "  contract     -> $CONTRACT"
say "  genotypes    -> $GTDIR/"
say "  final calls  -> $FINAL_OUT"
say "  VENV         = $VENV"
say "  queue=$QUEUE throttle=$THROTTLE  mem(disc/gt/comb/combgt)=${DISC_MEM}/${GT_MEM}/${COMBINE_MEM}/${COMBGT_MEM}"
say "----------------------------------------------------------------------"

# ============================================================================================
# 3. PREFLIGHT -- check every file/reference we can, collect ALL problems, then report once.
#    Nothing is submitted if anything required is missing.
# ============================================================================================
miss=()
note() { miss+=("$1"); }

# fofn + colonies
if [ ! -s "$FOFN" ]; then
    note "colony list FOFN not found or empty: $FOFN
      (build it with cluster/stage_picked.sh --fofn, or set FOFN=...)"
    N=0
else
    N="$(grep -c . "$FOFN")"
    # every listed BAM must exist (discovery/genotype read records; a missing one silently
    # shortens the cohort and skews the cross-donor gates).
    nbad=0; firstbad=""
    while read -r bam; do
        [ -n "$bam" ] || continue
        if [ ! -s "$bam" ]; then nbad=$((nbad+1)); [ -z "$firstbad" ] && firstbad="$bam"; fi
    done < "$FOFN"
    [ "$nbad" -eq 0 ] || note "$nbad of $N BAMs in the FOFN are missing/empty (first: $firstbad).
      Stage them to Lustre before running (cluster/stage_picked.sh)."
fi

# per-patient fofns / donors
[ -d "$FOFNDIR" ] || note "FOFNDIR does not exist: $FOFNDIR (it must hold all.bams.fofn + one <donor>.bams.fofn per donor)."
[ -n "$PATIENTS" ] || note "found no per-patient fofns (<donor>.bams.fofn) in $FOFNDIR.
      combine_insertions' assembly gate needs one BAM per donor. Re-run cluster/stage_picked.sh."

# configs
[ -s "$DISCOVER_CFG" ] || note "discovery config not found: $DISCOVER_CFG"
[ -s "$GENO_CFG" ]     || note "genotype config not found: $GENO_CFG"

# binaries
[ -x "$DISCOVER_BIN" ] || note "discovery binary missing/not executable: $DISCOVER_BIN
      build it:  module load rust/1.87.0 && bash cluster/build.sh"
[ -x "$GENOTYPE_BIN" ] || note "genotype binary missing/not executable: $GENOTYPE_BIN
      build it:  module load rust/1.87.0 && bash cluster/build.sh"

# python venv + config.py (combine steps)
[ -x "$VENV/bin/python" ] || note "venv python missing: $VENV/bin/python  (set VENV= or build the venv)"
[ -s "src/config.py" ]    || note "src/config.py missing -- the combine steps read reference paths + cohort gates from it,
      and it is gitignored. Copy the assembly template:  cp cluster/config.py.grch38 src/config.py"

# discovery reference: exon_annotation track, only required when splice_hallmark is on
splice="$(awk -F'[[:space:]]*=[[:space:]]*' '/^[[:space:]]*splice_hallmark[[:space:]]*=/{print $2}' "$DISCOVER_CFG" 2>/dev/null | tr -d ' ' || true)"
if [ "$splice" = "true" ]; then
    exon="$(awk -F'[[:space:]]*=[[:space:]]*' '/^[[:space:]]*exon_annotation[[:space:]]*=/{print $2}' "$DISCOVER_CFG" 2>/dev/null | tr -d ' ' || true)"
    if [ -z "$exon" ]; then
        note "$DISCOVER_CFG sets splice_hallmark=true but names no exon_annotation track.
      The discovery binary exits(1) without it. Add exon_annotation=... or set splice_hallmark=false."
    elif [ ! -s "$exon" ]; then
        note "exon_annotation track not found: $exon
      (splice_hallmark=true requires it; build with tools/build_grch38_exon_track.py or set splice_hallmark=false)."
    fi
fi

# combine assembly gate: genome_2bit + hs1->hg38 chain (a wrong pair = silently wrong coords).
# This is the project's own authoritative check; run it once now so we fail fast rather than
# discovering it after the discovery array has spent farm time.
if [ -z "$SKIP_ASSEMBLY_CHECK" ] && [ -x "$VENV/bin/python" ] && [ -s "src/config.py" ] && [ -s "$FOFN" ] && [ "${nbad:-1}" -eq 0 ]; then
    BAM1="$(head -1 "$FOFN")"
    if ! VENV="$VENV" bash cluster/install.sh check-config --bam "$BAM1" --pt-root "$PT_ROOT" >/tmp/ptcheck.$$ 2>&1; then
        note "the assembly gate (install.sh check-config) FAILED on $BAM1 -- genome_2bit / hs1->hg38 chain
      do not match these BAMs (this would produce silently wrong coordinates). Details:
$(sed 's/^/        /' /tmp/ptcheck.$$ | tail -12)
      Fetch the references:  bash cluster/install.sh install   (set SKIP_ASSEMBLY_CHECK=1 to bypass)"
    fi
    rm -f /tmp/ptcheck.$$
fi

if [ ${#miss[@]} -gt 0 ]; then
    m=""
    for x in "${miss[@]}"; do m+="  - $x"$'\n'; done
    die "$m"
fi
say "preflight OK -- all inputs, configs, binaries and references are present."
say "----------------------------------------------------------------------"

# ============================================================================================
# 4. per-stage state: COMPLETE (skip) / RUNNING (attach) / TODO (submit)
# ============================================================================================
fofn_ids() {
    awk 'NF{n=$0; sub(/.*\//,"",n); sub(/\.bam$/,"",n); sub(/\.cram$/,"",n); sub(/\.sample\.dupmarked$/,"",n); print n}' "$FOFN"
}
n_outputs() { # count <id>.txt.gz present in $1 for the fofn's ids
    local dir="$1" id c=0
    while read -r id; do [ -s "$dir/$id.txt.gz" ] && c=$((c+1)); done < <(fofn_ids)
    echo "$c"
}
running_id() { # numeric id of a PEND/RUN job matching name $1, or empty
    command -v bjobs >/dev/null 2>&1 || { echo ""; return 0; }
    bjobs -noheader -o 'id' -J "$1" 2>/dev/null | grep -oE '^[0-9]+' | head -1 || true
}

DISC_DONE="$(n_outputs "$DISCDIR")"
GT_DONE="$(n_outputs "$GTDIR")"
disc_complete()   { [ "$DISC_DONE" -eq "$N" ] && [ "$N" -gt 0 ]; }
comb_complete()   { [ -s "$CONTRACT" ]; }
gt_complete()     { [ "$GT_DONE" -eq "$N" ] && [ "$N" -gt 0 ]; }
combgt_complete() { [ -s "$FINAL_OUT" ]; }

# ============================================================================================
# 5. submit the incomplete stages, chaining by captured numeric job id.
# ============================================================================================
REPLY_ID=""
run_capture() { # runs "$@", echoes its output, sets REPLY_ID to the parsed LSF job id
    local out
    if out="$("$@" 2>&1)"; then
        printf '%s\n' "$out" | sed 's/^/    /'
        REPLY_ID="$(printf '%s\n' "$out" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -1)"
        [ -n "$REPLY_ID" ] || die "could not parse an LSF job id from the submission above."
    else
        printf '%s\n' "$out" | sed 's/^/    /'
        die "submission command failed (see output above); nothing further was submitted."
    fi
}

LAST=""          # numeric job id the NEXT stage must wait on ("" = no dependency)
declare -a PLAN  # human summary lines

wait_note() { [ -n "$LAST" ] && echo "after ended($LAST)" || echo "immediately"; }

# ---- stage 1: discovery ----
if disc_complete; then
    PLAN+=("discovery          SKIP  ($DISC_DONE/$N outputs already in $DISCDIR)")
elif [ -n "$DISC_JOBID" ]; then
    LAST="$DISC_JOBID"; PLAN+=("discovery          ATTACH to running job $DISC_JOBID (DISC_JOBID)")
elif rid="$(running_id "ptdisc_$(basename "$DISCDIR")")"; [ -n "$rid" ]; then
    LAST="$rid";       PLAN+=("discovery          ATTACH to running job $rid (detected by name)")
else
    say ">> submitting stage 1/4 discovery ($N colonies)"
    if [ -n "$DRYRUN" ]; then PLAN+=("discovery          SUBMIT $N-task array, $(wait_note)"); LAST="<discovery>"
    else
        e=(env FOFN="$FOFN" OUTDIR="$DISCDIR" DISCOVER_CFG="$DISCOVER_CFG" THROTTLE="$THROTTLE" QUEUE="$QUEUE" MEM="$DISC_MEM")
        [ -n "$GROUP" ] && e+=(GROUP="$GROUP")
        run_capture "${e[@]}" bash cluster/submit_discovery.sh
        LAST="$REPLY_ID"; PLAN+=("discovery          SUBMITTED job $LAST ($N-task array)")
    fi
fi

# ---- stage 2: combine_insertions (bsub combine_mei.sh; not a submit script of its own) ----
if comb_complete; then
    PLAN+=("combine_insertions SKIP  ($CONTRACT already exists)")
elif [ -n "$COMBINE_JOBID" ]; then
    LAST="$COMBINE_JOBID"; PLAN+=("combine_insertions ATTACH to running job $COMBINE_JOBID (COMBINE_JOBID)")
elif rid="$(running_id "ptcomb_$(basename "$INS_OUTDIR")")"; [ -n "$rid" ]; then
    LAST="$rid";           PLAN+=("combine_insertions ATTACH to running job $rid (detected by name)")
else
    say ">> submitting stage 2/4 combine_insertions ($(wait_note))"
    CTAG="$(basename "$INS_OUTDIR")"; mkdir -p logs
    if [ -n "$DRYRUN" ]; then PLAN+=("combine_insertions SUBMIT bsub combine_mei.sh, $(wait_note)"); LAST="<combine>"
    else
        b=(bsub -J "ptcomb_$CTAG" -o "logs/comb.$CTAG.out" -e "logs/comb.$CTAG.err"
           -n "$COMBINE_CORES" -q "$QUEUE"
           -R "select[mem>${COMBINE_MEM}] rusage[mem=${COMBINE_MEM}] span[hosts=1]" -M "$COMBINE_MEM")
        [ -n "$LAST" ]  && b+=(-w "ended($LAST)")
        [ -n "$GROUP" ] && b+=(-G "$GROUP")
        b+=(bash -c "FOFNDIR='$FOFNDIR' DISCDIR='$DISCDIR' OUTDIR='$INS_OUTDIR' STEM='$STEM' PATIENTS='$PATIENTS' THREADS='$COMBINE_CORES' VENV='$VENV' ALLOW_MISSING='$ALLOW_MISSING' bash cluster/combine_mei.sh")
        run_capture "${b[@]}"
        LAST="$REPLY_ID"; PLAN+=("combine_insertions SUBMITTED job $LAST")
    fi
fi

# ---- stage 3: genotype (submit_genotype.sh supports WAIT + defers contract checks) ----
if gt_complete; then
    PLAN+=("genotype           SKIP  ($GT_DONE/$N outputs already in $GTDIR)")
elif [ -n "$GT_JOBID" ]; then
    LAST="$GT_JOBID"; PLAN+=("genotype           ATTACH to running job $GT_JOBID (GT_JOBID)")
elif rid="$(running_id "ptgt_$(basename "$GTDIR")")"; [ -n "$rid" ]; then
    LAST="$rid";      PLAN+=("genotype           ATTACH to running job $rid (detected by name)")
else
    say ">> submitting stage 3/4 genotype ($N colonies, $(wait_note))"
    if [ -n "$DRYRUN" ]; then PLAN+=("genotype           SUBMIT $N-task array, $(wait_note)"); LAST="<genotype>"
    else
        e=(env FOFN="$FOFN" OUTDIR="$GTDIR" CONTRACT="$CONTRACT" GENO_CFG="$GENO_CFG" THROTTLE="$THROTTLE" QUEUE="$QUEUE" MEM="$GT_MEM" GENOTYPE_BIN="$GENOTYPE_BIN")
        [ -n "$LAST" ]  && e+=(WAIT="ended($LAST)")
        [ -n "$GROUP" ] && e+=(GROUP="$GROUP")
        run_capture "${e[@]}" bash cluster/submit_genotype.sh
        LAST="$REPLY_ID"; PLAN+=("genotype           SUBMITTED job $LAST ($N-task array)")
    fi
fi

# ---- stage 4: combine_genotypes (submit_combine_genotypes.sh supports WAIT) ----
if combgt_complete; then
    PLAN+=("combine_genotypes  SKIP  ($FINAL_OUT already exists)")
else
    say ">> submitting stage 4/4 combine_genotypes ($(wait_note))"
    if [ -n "$DRYRUN" ]; then PLAN+=("combine_genotypes  SUBMIT bsub combine_gt.sh, $(wait_note)")
    else
        e=(env FOFN="$FOFN" GTDIR="$GTDIR" OUT="$FINAL_OUT" MEM="$COMBGT_MEM" CORES="$COMBGT_CORES" QUEUE="$QUEUE" VENV="$VENV" ALLOW_MISSING="$ALLOW_MISSING")
        [ -n "$LAST" ]  && e+=(WAIT="ended($LAST)")
        [ -n "$GROUP" ] && e+=(GROUP="$GROUP")
        run_capture "${e[@]}" bash cluster/submit_combine_genotypes.sh
        LAST="$REPLY_ID"; PLAN+=("combine_genotypes  SUBMITTED job $LAST")
    fi
fi

# ============================================================================================
# 6. summary
# ============================================================================================
say "----------------------------------------------------------------------"
if [ -n "$DRYRUN" ]; then say " DRY RUN -- nothing was submitted. Plan:"; else say " Submitted. Plan:"; fi
for line in "${PLAN[@]}"; do say "   $line"; done
say "----------------------------------------------------------------------"
if [ -z "$DRYRUN" ]; then
    say " watch:   bjobs -A"
    say " final:   zcat $FINAL_OUT | tail -n +2 | wc -l   # loci passing the cohort gates"
    say " re-run this same command any time -- completed stages are skipped, so it is a safe"
    say " way to resume after a partial failure."
fi
