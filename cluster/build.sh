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

echo "cargo: $(command -v cargo || echo MISSING)"
cargo --version

cargo build --release --manifest-path rust/peartree-discovery/Cargo.toml
cargo build --release --manifest-path rust/peartree-genotype/Cargo.toml

echo
echo "built:"
ls -la rust/peartree-discovery/target/release/peartree-discovery
ls -la rust/peartree-genotype/target/release/peartree-genotype
