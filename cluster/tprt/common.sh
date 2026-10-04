# shellcheck shell=bash
# shellcheck disable=SC2034  # variables are consumed by the scripts that source this file
# =============================================================================
# cluster/tprt/common.sh — shared paths + resolvers for the TPRT A/B kit.
# Sourced by setup.sh, preflight.sh, run_ab.sh, evaluate.sh. Never executed directly.
#
# Every value is either a farm22 convention (overridable by env, the same way pipeline.sh
# and fleet.sh do it) or DERIVED from files in the repo (patient dir, tree, sample list).
#
#   arm A = the current pipeline: cluster/config.discovery.grch38 + config.genotype.grch38 +
#           config.py.grch38, run from a second worktree (PT_ROOT_A) at the SAME commit, so the
#           only difference between the arms is the configuration (every TPRT key defaults off
#           in code; the legacy configs are byte-identical to the pre-TPRT pipeline).
#   arm B = the TPRT-hallmark pipeline: the .tprt configs, run from this checkout (PT_ROOT_B).
#
# Two checkouts, because src/config.py is a single gitignored module that main.py and
# annotate_v2.py import from their own src/ directory: one checkout can only hold one of them,
# and the arms run concurrently.
# =============================================================================

# --- checkouts ----------------------------------------------------------------
TPRT_KIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PT_ROOT_B="${PT_ROOT_B:-$(cd "$TPRT_KIT_DIR/../.." && pwd)}"

# --- farm22 conventions (env overrides, like pipeline.sh) -----------------------
JD="${JD:-/lustre/scratch126/casm/teams/team273/users/jd43}"
TPRT_ROOT="${TPRT_ROOT:-$JD/tprt_ab}"                       # everything this kit creates
PT_ROOT_A="${PT_ROOT_A:-$TPRT_ROOT/PEAR-TREE-A}"            # arm-A worktree (same commit)
# One venv for both arms (same code). It lives at <arm-B checkout>/venv (gitignored) and the arm-A
# worktree gets a venv -> symlink, because install.sh check-config (called by pipeline.sh's combine
# task WITHOUT passing VENV) defaults to $PT_ROOT/venv.
VENV="${VENV:-$PT_ROOT_B/venv}"
TPRT_RES="${TPRT_RES:-$TPRT_ROOT/resources}"                # hs1.2bit, hs1 gene model, hs1 .mmi
STAGING_ROOT="${STAGING_ROOT:-/lustre/scratch126/casm/staging/team273/jd43}"
RESULTS_BASE="${RESULTS_BASE:-$HOME/results/tprt_ab}"       # NFS (backed up) copies of finals
NST="${NST:-/nfs/cancer_ref01/nst_links/live}"              # header reads ONLY (Sanger policy)
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-samtools-1.19}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.12.3}"
RUST_MODULE="${RUST_MODULE:-rust/1.87.0}"
MINIMAP2_VERSION="${MINIMAP2_VERSION:-2.28}"
MINIMAP2="${MINIMAP2:-$TPRT_ROOT/bin/minimap2}"
UCSC="${UCSC:-https://hgdownload.soe.ucsc.edu/goldenPath}"
PATIENTS_DIR="${PATIENTS_DIR:-$PT_ROOT_B/patients}"

# Binaries: built once in PT_ROOT_B; arm A runs the SAME binaries (same commit; the legacy
# configs leave every new key off, so the code path is the pre-TPRT one).
DISCOVER_BIN="${DISCOVER_BIN:-$PT_ROOT_B/rust/peartree-discovery/target/release/peartree-discovery}"
GENOTYPE_BIN="${GENOTYPE_BIN:-$PT_ROOT_B/rust/peartree-genotype/target/release/peartree-genotype}"
BUILD_STAMP="$TPRT_ROOT/build.stamp"

# --- per-arm configs -------------------------------------------------------------
# arm C = arm B with >= 2 distinct fragments per junction WITHIN each colony at discovery
# (config.discovery.grch38.tprt2frag); same checkout, combine/genotype/annotate configs as B.
ALL_ARMS=(A B C)
arm_root()     { case "$1" in A) echo "$PT_ROOT_A" ;; B|C) echo "$PT_ROOT_B" ;; *) return 1 ;; esac; }
arm_disc_cfg() { case "$1" in A) echo "$PT_ROOT_A/cluster/config.discovery.grch38" ;; B) echo "$PT_ROOT_B/cluster/config.discovery.grch38.tprt" ;; C) echo "$PT_ROOT_B/cluster/config.discovery.grch38.tprt2frag" ;; esac; }
arm_geno_cfg() { case "$1" in A) echo "$PT_ROOT_A/cluster/config.genotype.grch38" ;; B|C) echo "$PT_ROOT_B/cluster/config.genotype.grch38.tprt" ;; esac; }
arm_py_base()  { case "$1" in A) echo "cluster/config.py.grch38" ;; B|C) echo "cluster/config.py.grch38.tprt" ;; esac; }

# --- per-patient layout ------------------------------------------------------------
#   $TPRT_ROOT/<P>/samples.tsv            sample<TAB>proj (fleet.sh samples <P>)
#   $TPRT_ROOT/<P>/<arm>/                 pipeline WORKROOT of the arm
#   $TPRT_ROOT/<P>/<arm>/<P>/             pipeline RUNDIR  (discovery/ insertions/ genotypes/ logs/ ...)
#   $TPRT_ROOT/<P>/<arm>/jobids.tsv       phase<TAB>jobid of the last submission
#   $TPRT_ROOT/<P>/eval/                  tree_fit + discrimination per arm, ab_report.md
pat_dir()      { echo "$TPRT_ROOT/$1"; }
arm_workroot() { echo "$TPRT_ROOT/$1/$2"; }            # <P> <arm>
arm_rundir()   { echo "$TPRT_ROOT/$1/$2/$1"; }         # <P> <arm>
arm_results()  { echo "$RESULTS_BASE/$1/$2"; }         # <P> <arm>  (pipeline appends /<P>)
eval_dir()     { echo "$TPRT_ROOT/$1/eval"; }

die()  { echo "ERROR: $*" >&2; exit 1; }
note() { echo "[$(date +%H:%M:%S)] $*"; }

# patients/<organ>/<P>/ — exactly one, else die
patient_dir() {
    local P="$1" hits
    hits="$(find "$PATIENTS_DIR" -mindepth 2 -maxdepth 2 -type d -name "$P" | sort)"
    [ -n "$hits" ] || die "no patients/*/$P directory under $PATIENTS_DIR"
    [ "$(printf '%s\n' "$hits" | wc -l | tr -d ' ')" -eq 1 ] || die "several patients/*/$P directories: $(printf '%s' "$hits" | tr '\n' ' ')"
    echo "$hits"
}

# the patient's newick tree — exactly one *.tree file in the patient dir, else die
patient_tree() {
    local d; d="$(patient_dir "$1")" || exit 1
    local trees=("$d"/*.tree)
    [ -e "${trees[0]}" ] || die "no *.tree in $d"
    [ "${#trees[@]}" -eq 1 ] || die "several *.tree in $d: ${trees[*]}"
    echo "${trees[0]}"
}

patient_colonies() { echo "$(patient_dir "$1")/colonies.tsv"; }

# tree tip labels, one per line (newick leaves: a label directly after '(' or ',')
tree_tips() {
    tr -d '\n\r ' < "$1" | grep -oE '[(,][^(),:;]+' | sed 's/^[(,]//' | sort -u
}

# sample<TAB>proj for the patient, through the fleet's one rule (fleet.sh samples)
patient_samples() {
    PATIENTS_DIR="$PATIENTS_DIR" bash "$PT_ROOT_B/cluster/fleet.sh" samples "$1"
}

# staged BAM path — the same convention pipeline.sh's bam_path() uses
staged_bam() { echo "$STAGING_ROOT/$2/$1/mapped_sample/$1.sample.dupmarked.bam"; }   # <sample> <proj>
nst_bam()    { echo "$NST/$2/$1/$1.sample.dupmarked.bam"; }                          # <sample> <proj>

# iRODS (cgp zone) is the authority for what exists; nst_links may be absent on a node
# (it was gone from farm22-head2 by 2026-10). DATA_PATH is the physical replica path,
# printed WITHOUT the /nfs prefix (e.g. /irods-cgp-sr13-sdf/intproj/...); one line per replica.
load_irods() {
    command -v iquest >/dev/null 2>&1 && return 0
    # shellcheck disable=SC1091
    [ -n "${MODULESHOME:-}" ] && [ -f "$MODULESHOME/init/bash" ] && . "$MODULESHOME/init/bash" 2>/dev/null
    module load IRODS >/dev/null 2>&1 || true   # CAPITALISED; lowercase `irods` does not exist
    command -v iquest >/dev/null 2>&1
}
irods_replicas() {   # <sample> <proj> -> physical replica paths of the dupmarked BAM (may be empty)
    iquest --no-page "%s" "select DATA_PATH where COLL_NAME = '/cgp/intproj/$2/sample/$1' and DATA_NAME like '$1.%sample.dupmarked.bam'" 2>/dev/null \
        | grep -v -e CAT_NO_ROWS -e '^$' || true
}
# first header-readable physical replica (tries the path as printed and under /nfs), else empty
irods_readable_bam() {   # <sample> <proj>
    local r c
    while IFS= read -r r; do
        for c in "$r" "/nfs$r"; do
            samtools view -H "$c" >/dev/null 2>&1 && { echo "$c"; return 0; }
        done
    done < <(irods_replicas "$1" "$2")
    return 1
}

load_samtools() {
    command -v samtools >/dev/null 2>&1 && return 0
    # shellcheck disable=SC1091
    [ -n "${MODULESHOME:-}" ] && [ -f "$MODULESHOME/init/bash" ] && . "$MODULESHOME/init/bash" 2>/dev/null
    module load "$SAMTOOLS_MODULE" >/dev/null 2>&1 || true
    command -v samtools >/dev/null 2>&1
}

md5_of() {   # portable md5 (Linux md5sum / macOS md5)
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi
}
