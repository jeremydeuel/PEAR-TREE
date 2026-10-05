#!/usr/bin/env bash
# Build the PEAR-TREE Rust binaries (release) on the farm.
#
# Needs a Rust toolchain (>= 1.87). On farm22, either:
#   module avail rust        # and `module load` a >=1.87 rust, OR
#   curl https://sh.rustup.rs -sSf | sh    # per-user rustup (head node has proxy internet)
# The build fetches the noodles crates from crates.io, so run it on a node with
# outbound internet (the head node), not a compute node.
set -euo pipefail
cd "$(dirname "$0")/.."

# Load the rust module ourselves. It was documented in the header above and left to the
# caller to remember, which is not a plan: there is no cargo on farm22's default PATH, so
# forgetting it fails the build. That failure is loud and harmless. The one that is NOT is
# skipping the rebuild entirely — the binary then silently predates the config it is handed,
# and config.rs tolerates unknown keys BY DESIGN (warns, continues). That is exactly how an
# A/B of the slippage gate ran 180 BAM scans with both arms executing identical code: the
# binary had no `slippage_filter`, said "ignoring unknown config key" 209 times into stderr,
# and reported a 0.0% effect that meant "not tested".
#
# Must come BEFORE the CARGO_HOME check below: it is the module that sets the read-only value.
RUST_MODULE="${RUST_MODULE:-rust/1.87.0}"
if ! command -v cargo >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1 || [ -n "${MODULESHOME:-}" ]; then
        echo "no cargo on PATH — module load $RUST_MODULE"
        # `module` is a shell function; in a non-interactive shell it may need sourcing first.
        [ -n "${MODULESHOME:-}" ] && [ -f "$MODULESHOME/init/bash" ] && . "$MODULESHOME/init/bash" || true
        module load "$RUST_MODULE" 2>/dev/null || echo "  module load $RUST_MODULE failed — trying anyway"
    fi
fi
command -v cargo >/dev/null 2>&1 || {
    echo "cargo: MISSING" >&2
    echo "No Rust toolchain. On farm22:  module load $RUST_MODULE" >&2
    echo "(module avail rust  — to see what is there; needs >= 1.87)" >&2
    echo "Or a per-user rustup:  curl https://sh.rustup.rs -sSf | sh" >&2
    exit 1; }

# Some rust modules (e.g. farm22 `module load rust/1.87.0`) set CARGO_HOME to their
# own read-only install dir, which makes crate downloads fail with EACCES. If the
# current CARGO_HOME is set but not writable, redirect it to a local cache next to
# the build (the repo lives on lustre, which has space). An unset CARGO_HOME keeps
# cargo's default (~/.cargo).
if [ -n "${CARGO_HOME:-}" ] && [ ! -w "${CARGO_HOME}" ]; then
    export CARGO_HOME="$PWD/.cargo-home"
    mkdir -p "$CARGO_HOME"
    echo "CARGO_HOME was read-only; redirected to $CARGO_HOME"
fi

echo "cargo: $(command -v cargo)"
cargo --version

cargo build --release --manifest-path rust/peartree-discovery/Cargo.toml
cargo build --release --manifest-path rust/peartree-genotype/Cargo.toml
# combine_insertions port (only used with COMBINE_IMPL=rust; the python combine does not need
# it, so a failure here must not block the discovery/genotype builds above)
cargo build --release --manifest-path rust/peartree-combine/Cargo.toml \
    || echo "WARNING: rust/peartree-combine failed to build — COMBINE_IMPL=rust unavailable (python combine unaffected)" >&2

echo
echo "built:"
ls -la rust/peartree-discovery/target/release/peartree-discovery
ls -la rust/peartree-genotype/target/release/peartree-genotype
ls -la rust/peartree-combine/target/release/peartree-combine 2>/dev/null || true

# Report which OPTIONAL, CONFIG-GATED features this binary actually implements. A stale binary
# does not announce itself: config.rs tolerates unknown keys by design, so an old build fed a
# new config warns to stderr and carries on doing something else. The slippage A/B lost 180
# BAM scans to exactly that. These markers are the banner literals main.rs prints, so their
# presence in the binary is an exact test of what is compiled in.
#
# grep the binary DIRECTLY (-a). NEVER `strings "$BIN" | grep -q ...`: under `set -o pipefail`
# that returns FAILURE on a SUCCESSFUL match (grep -q exits first, strings takes SIGPIPE) —
# the same trap that labelled every GRCh38 BAM "chr1" in catalogue_headers.sh.
echo
echo "discovery features compiled in:"
BIN=rust/peartree-discovery/target/release/peartree-discovery
for feat in 'slippage filter:' 'contig allowlist:' 'coverage mask:' 'adaptive evidence floor:' 'RM self-mask:'; do
    if LC_ALL=C grep -qaF "$feat" "$BIN" 2>/dev/null; then
        printf '  yes  %s\n' "${feat%:}"
    else
        printf '  NO   %s   <-- configs setting this would be SILENTLY IGNORED\n' "${feat%:}"
    fi
done
echo
echo "Before trusting any run, check the banner reports the feature you are measuring:"
echo "  $BIN --step discover --bam <bam> --out /tmp/x.gz --config <cfg> 2>&1 | head"
