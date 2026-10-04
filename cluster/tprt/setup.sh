#!/usr/bin/env bash
# =============================================================================
# cluster/tprt/setup.sh — deploy the TPRT A/B kit on farm22 (HEAD NODE: needs internet).
#
#   bash cluster/tprt/setup.sh all          # every step below, in order (idempotent)
#   bash cluster/tprt/setup.sh worktree     # arm-A worktree at THIS checkout's commit
#   bash cluster/tprt/setup.sh build        # both Rust binaries (rust module, CARGO_HOME fix) + build.stamp
#   bash cluster/tprt/setup.sh venv         # $VENV: requirements.txt (+edlib, mappy) + scipy, matplotlib
#   bash cluster/tprt/setup.sh minimap2     # module, else the official static release binary
#   bash cluster/tprt/setup.sh resources    # hs1.2bit (link existing / download) + hs1 RefSeq gene model
#   bash cluster/tprt/setup.sh hmmer        # nhmmscan for annotate: old dir / module / PATH / build from source
#   bash cluster/tprt/setup.sh index        # bsub: hs1 minimap2 sr index (optional; novel-source locator)
#   bash cluster/tprt/setup.sh configs      # src/config.py stubs in both checkouts + tmp dirs
#
# Layout (all under $TPRT_ROOT = /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab):
#   PEAR-TREE/      this checkout = arm B (branch tprt-hallmarks), binaries built here
#   PEAR-TREE-A/    git worktree, detached at the same commit = arm A
#   PEAR-TREE/venv (shared; PEAR-TREE-A/venv -> it)  bin/minimap2  resources/{hs1.2bit,hs1.gene_model.tsv.gz,hs1.sr.mmi}  tmp/{A,B}/
#   build.stamp     commit + md5 of the binaries (preflight.sh checks it)
# =============================================================================
set -euo pipefail
# shellcheck source=cluster/tprt/common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

modinit() {
    if ! command -v module >/dev/null 2>&1 && [ -n "${MODULESHOME:-}" ] && [ -f "$MODULESHOME/init/bash" ]; then
        # shellcheck disable=SC1091
        . "$MODULESHOME/init/bash"
    fi
}

fetch() {   # fetch <url> <dest>   (atomic)
    mkdir -p "$(dirname "$2")"
    note "downloading $1"
    curl -fL --retry 3 -o "$2.part" "$1"
    mv -f "$2.part" "$2"
}

step_worktree() {
    local head; head="$(git -C "$PT_ROOT_B" rev-parse HEAD)"
    if [ -d "$PT_ROOT_A/.git" ] || [ -f "$PT_ROOT_A/.git" ]; then
        git -C "$PT_ROOT_A" checkout --detach "$head"
        note "arm A worktree $PT_ROOT_A moved to ${head:0:12}"
    else
        mkdir -p "$(dirname "$PT_ROOT_A")"
        git -C "$PT_ROOT_B" worktree add --detach "$PT_ROOT_A" "$head"
        note "arm A worktree created: $PT_ROOT_A @ ${head:0:12}"
    fi
    [ "$(git -C "$PT_ROOT_A" rev-parse HEAD)" = "$head" ] || die "worktree is not at $head"
    # install.sh check-config (run by pipeline.sh's combine task without VENV) defaults to
    # $PT_ROOT/venv — give the arm-A worktree the shared venv under that name
    [ "$VENV" = "$PT_ROOT_A/venv" ] || ln -sfn "$VENV" "$PT_ROOT_A/venv"
}

step_build() {
    modinit
    command -v cargo >/dev/null 2>&1 || module load "$RUST_MODULE"
    # rust/1.87.0 sets CARGO_HOME to its read-only install dir; cluster/build.sh redirects an
    # unwritable CARGO_HOME to $PT_ROOT_B/.cargo-home (lustre, visible to compute nodes).
    bash "$PT_ROOT_B/cluster/build.sh"
    mkdir -p "$TPRT_ROOT"
    { printf 'commit\t%s\n' "$(git -C "$PT_ROOT_B" rev-parse HEAD)"
      printf 'built\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
      for b in "$DISCOVER_BIN" "$GENOTYPE_BIN"; do printf 'md5\t%s\t%s\n' "$(md5_of "$b")" "$b"; done
    } > "$BUILD_STAMP"
    note "build.stamp: $BUILD_STAMP"; cat "$BUILD_STAMP"
}

step_venv() {
    modinit
    if [ ! -x "$VENV/bin/python" ]; then
        module load "$PYTHON_MODULE" >/dev/null 2>&1 || true
        note "creating $VENV with $(python3 --version 2>&1)"
        python3 -m venv "$VENV"
    fi
    "$VENV/bin/pip" install -q --upgrade pip
    "$VENV/bin/pip" install -q -r "$PT_ROOT_B/requirements.txt"     # incl. edlib, mappy
    "$VENV/bin/pip" install -q scipy matplotlib
    "$VENV/bin/python" -c "import pysam, pyliftover, py2bit, numpy, pandas, edlib, mappy, scipy, matplotlib; print('venv OK: pysam', pysam.__version__, 'mappy', mappy.__version__)"
}

step_minimap2() {
    modinit
    if [ -x "$MINIMAP2" ]; then note "minimap2 present: $MINIMAP2 ($("$MINIMAP2" --version))"; return 0; fi
    local mod; mod="$(module -t avail minimap2 2>&1 | grep -iE '^minimap2' | sort -V | tail -1 || true)"
    if [ -n "$mod" ] && module load "$mod" 2>/dev/null && command -v minimap2 >/dev/null 2>&1; then
        mkdir -p "$(dirname "$MINIMAP2")"
        ln -sf "$(command -v minimap2)" "$MINIMAP2"
        note "minimap2 from module $mod -> $MINIMAP2"
    else
        local t="minimap2-${MINIMAP2_VERSION}_x64-linux"
        fetch "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/$t.tar.bz2" "$TPRT_ROOT/bin/$t.tar.bz2"
        tar -xjf "$TPRT_ROOT/bin/$t.tar.bz2" -C "$TPRT_ROOT/bin"
        ln -sf "$TPRT_ROOT/bin/$t/minimap2" "$MINIMAP2"
        note "minimap2 static binary -> $MINIMAP2"
    fi
    "$MINIMAP2" --version
}

step_resources() {
    mkdir -p "$TPRT_RES"
    # hs1.2bit: reuse a copy already on lustre, else download it
    if [ ! -s "$TPRT_RES/hs1.2bit" ]; then
        local c found=""
        for c in "$JD/hs1.2bit" "$JD/pt_hu_trees/hs1.2bit" "$JD/pt_hu_trees/hs1/hs1.2bit" "$JD/hs1/hs1.2bit"; do
            [ -s "$c" ] && { found="$c"; break; }
        done
        if [ -n "$found" ]; then ln -sf "$found" "$TPRT_RES/hs1.2bit"; note "hs1.2bit -> $found"
        else fetch "$UCSC/hs1/bigZips/hs1.2bit" "$TPRT_RES/hs1.2bit"; fi
    fi
    # hs1 RefSeq gene model (pseudogene exons for annotate; E2E_REPORT "Annotate round 2")
    # valid = gzip-intact AND non-empty (an earlier kit wrote it uncompressed via a .gz.part name)
    local gm="$TPRT_RES/hs1.gene_model.tsv.gz"
    if ! { [ -s "$gm" ] && gzip -t "$gm" 2>/dev/null && [ "$(gzip -dc "$gm" | head -2 | wc -l)" -gt 1 ]; }; then
        [ -e "$gm" ] && { note "rebuilding invalid/empty gene model $gm"; rm -f "$gm"; }
        [ -s "$TPRT_RES/hs1.ncbiRefSeq.gtf.gz" ] && gzip -t "$TPRT_RES/hs1.ncbiRefSeq.gtf.gz" 2>/dev/null \
            || fetch "$UCSC/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz" "$TPRT_RES/hs1.ncbiRefSeq.gtf.gz"
        [ -x "$VENV/bin/python" ] || die "venv first: bash $TPRT_KIT_DIR/setup.sh venv"
        # the temp name must END in .gz: build_gene_model.py chooses gzip by the output extension
        "$VENV/bin/python" "$PT_ROOT_B/tools/build_gene_model.py" --curated \
            "$TPRT_RES/hs1.ncbiRefSeq.gtf.gz" "$TPRT_RES/hs1.gene_model.part.tsv.gz"
        gzip -t "$TPRT_RES/hs1.gene_model.part.tsv.gz" || die "gene model build produced an invalid gzip"
        mv -f "$TPRT_RES/hs1.gene_model.part.tsv.gz" "$gm"
    fi
    note "hs1 gene model: $TPRT_RES/hs1.gene_model.tsv.gz ($(gzip -dc "$TPRT_RES/hs1.gene_model.tsv.gz" | wc -l) rows; ~252,903 expected)"
}

# HMMER (nhmmscan, used by annotate's Dfam scan). The base configs point at
# $JD/hmmer-3.3.2/bin, which no longer exists. Resolve, in order: that dir if it has nhmmscan;
# ~/.local/bin (where jd43's HMMER lives); an `hmmer` module; nhmmscan already on PATH; else build HMMER from source on the head node.
# The chosen bin dir is written to $TPRT_RES/hmmer_bin, which arm_config.py uses for both arms.
HMMER_VERSION="${HMMER_VERSION:-3.4}"
step_hmmer() {
    mkdir -p "$TPRT_RES"
    local d="" m
    if [ -x "$JD/hmmer-3.3.2/bin/nhmmscan" ]; then d="$JD/hmmer-3.3.2/bin"
    elif [ -x "$HOME/.local/bin/nhmmscan" ]; then d="$HOME/.local/bin"        # jd43's existing install
    elif [ -x "$TPRT_ROOT/hmmer-$HMMER_VERSION/bin/nhmmscan" ]; then d="$TPRT_ROOT/hmmer-$HMMER_VERSION/bin"
    else
        modinit
        for m in hmmer "hmmer/$HMMER_VERSION" hmmer-3.4 hmmer-3.3.2; do
            module load "$m" >/dev/null 2>&1 && command -v nhmmscan >/dev/null 2>&1 && break
        done
        command -v nhmmscan >/dev/null 2>&1 && d="$(dirname "$(command -v nhmmscan)")"
    fi
    if [ -z "$d" ]; then
        local src="$TPRT_ROOT/src/hmmer-$HMMER_VERSION"
        [ -d "$src" ] || { fetch "http://eddylab.org/software/hmmer/hmmer-$HMMER_VERSION.tar.gz" "$TPRT_ROOT/src/hmmer-$HMMER_VERSION.tar.gz"
                           tar -xzf "$TPRT_ROOT/src/hmmer-$HMMER_VERSION.tar.gz" -C "$TPRT_ROOT/src"; }
        note "building HMMER $HMMER_VERSION from source (~5 min)"
        ( cd "$src" && ./configure --prefix="$TPRT_ROOT/hmmer-$HMMER_VERSION" >/dev/null && make -j4 >/dev/null && make install >/dev/null )
        d="$TPRT_ROOT/hmmer-$HMMER_VERSION/bin"
    fi
    [ -x "$d/nhmmscan" ] || die "no usable nhmmscan (tried $d)"
    printf '%s\n' "$d" > "$TPRT_RES/hmmer_bin"
    note "HMMER: $d ($("$d/nhmmscan" -h | sed -n 2p | sed 's/^# *//'))"
}

step_index() {
    [ -s "$TPRT_RES/hs1.sr.mmi" ] && { note "index present: $TPRT_RES/hs1.sr.mmi"; return 0; }
    [ -x "$MINIMAP2" ] || step_minimap2
    [ -s "$TPRT_RES/hs1.fa.gz" ] || fetch "$UCSC/hs1/bigZips/hs1.fa.gz" "$TPRT_RES/hs1.fa.gz"   # head node: internet
    mkdir -p "$TPRT_ROOT/logs"
    bsub -J tprt_hs1_mmi -n 4 -q normal -M 32000 -R "select[mem>32000] rusage[mem=32000] span[hosts=1]" \
        -o "$TPRT_ROOT/logs/hs1_mmi.%J.log" -e "$TPRT_ROOT/logs/hs1_mmi.%J.err" \
        "'$MINIMAP2' -x sr -t 4 -d '$TPRT_RES/hs1.sr.mmi.part' '$TPRT_RES/hs1.fa.gz' && mv -f '$TPRT_RES/hs1.sr.mmi.part' '$TPRT_RES/hs1.sr.mmi'"
    note "index job submitted (~30-60 min); preflight warns until $TPRT_RES/hs1.sr.mmi exists"
}

step_configs() {
    local ARM R
    for ARM in A B; do
        R="$(arm_root "$ARM")"
        [ -d "$R" ] || die "no checkout for arm $ARM at $R (setup.sh worktree first)"
        if [ -e "$R/src/config.py" ] && ! grep -q 'TPRT_AB_ARM' "$R/src/config.py"; then
            cp -f "$R/src/config.py" "$R/src/config.py.pre_tprt_ab.$(date +%Y%m%d%H%M%S)"
            note "kept the previous $R/src/config.py as a timestamped .pre_tprt_ab backup"
        fi
        cat > "$R/src/config.py" <<EOF
# GENERATED by cluster/tprt/setup.sh configs — arm $ARM of the TPRT A/B kit. Do not edit:
# the configuration lives in $(arm_py_base "$ARM") + cluster/tprt/arm_config.py.
import os
import sys
os.environ.setdefault('TPRT_ROOT', '$TPRT_ROOT')   # baked at generation: LSF jobs need no env
os.environ.setdefault('TPRT_RES', '$TPRT_RES')
sys.path.insert(0, '$R/cluster/tprt')
from arm_config import build  # noqa: E402
TPRT_AB_ARM = '$ARM'
CONFIG = build(TPRT_AB_ARM)
EOF
        mkdir -p "$TPRT_ROOT/tmp/$ARM"
        (cd "$R" && TPRT_ROOT="$TPRT_ROOT" TPRT_RES="$TPRT_RES" "$VENV/bin/python" -c \
            "import sys; sys.path.insert(0, 'src'); import config; sys.path.insert(0, 'cluster/tprt'); import arm_config; print(arm_config.describe(config.TPRT_AB_ARM))")
    done
    note "assembly gate: preflight.sh runs install.sh check-config for both arms against a BAM of the patient it is given"
}

case "${1:-}" in
    all)       step_worktree; step_build; step_venv; step_minimap2; step_resources; step_hmmer; step_configs; step_index ;;
    worktree)  step_worktree ;;
    build)     step_build ;;
    venv)      step_venv ;;
    minimap2)  step_minimap2 ;;
    resources) step_resources ;;
    hmmer)     step_hmmer ;;
    index)     step_index ;;
    configs)   step_configs ;;
    *) sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
