# peartree-combine

Rust port of PEAR-TREE step 2, `python src/main.py --step combine_insertions`
(`src/combine_insertions*.py`, `src/indel_consensus.py`, ...). Same inputs, same config file, same
output files, byte-identical after decompression. The one exception is the pooled
≥2-independent-fragment gate, which is not ported (see "Differences from Python" below).
`SPEC.md` describes the behaviour, and `PLAN.md` covers the memory and parallelism design.

## Build

```bash
cargo build --release --manifest-path rust/peartree-combine/Cargo.toml
# -> rust/peartree-combine/target/release/peartree-combine
```

You need Rust >= 1.87 and a C/C++ compiler: the build compiles the vendored edlib 1.2.7 (the
library python-edlib wraps), minimap2 2.30 (through the `minimap2` crate) and jemalloc. On the
farm, `bash cluster/build.sh` builds it next to the discovery and genotype binaries. That script
loads `rust/1.87.0` and redirects a read-only `CARGO_HOME`. A failed combine build only prints a
warning, because the Python combine does not need this binary.

## Run

```bash
PEARTREE_PYTHON=/path/to/venv/bin/python \
  rust/peartree-combine/target/release/peartree-combine --step combine_insertions \
  --config src/config.py --discovery_files discovery/*.txt.gz --out insertions/PD12345 --threads 8
```

* `--config` takes the same `config.py` the Python reads. A `.py` config is executed by
  `$PEARTREE_PYTHON` (default `python3`), which dumps `CONFIG` to JSON. A `.json` file is read
  directly. A relative `rte_library` is resolved like Python: first against the cwd, then against
  `<config dir>/..` (the repo root when the config is in `src/`), `$PT_ROOT`, and the binary's
  repo.
* The CLI accepts the Python spellings (`--discovery_files`, `--out`, `--threads`). The output
  stem strips `.gz`, as main.py does.
* Outputs: `<stem>.combined.txt.gz`, `.genotyping.txt.gz`, `.fq.gz`, `.bam`,
  `.insertionsonly.bam`, `.combined.splice.tsv`, and with evidence sidecars
  `.insertions.evidence.tsv.gz` and `.insertions.reads.fa.gz`. Resume semantics match Python: an
  existing `.bam` or `.insertionsonly.bam` skips that bowtie2 run.
* Scratch: `<stem>.evidence_shards/` holds the routed sidecar rows (deflated), rendered reads and
  the absorb lookup index. It is removed at the end.
* `PEARTREE_MEMLOG=1` prints stage timestamps and peak RSS to stderr.

### Cluster kit

`cluster/pipeline.sh` (and through it `cluster/fleet.sh`) and `cluster/combine_mei.sh` take
`COMBINE_IMPL=python|rust`. The default is `python`. With `rust`, they run `$COMBINE_BIN` (default
`rust/peartree-combine/target/release/peartree-combine`) with `--config $PT_ROOT/src/config.py`
and `PEARTREE_PYTHON=$VENV/bin/python`. `pipeline.sh` freezes the choice into `run.env`, and its
preflight fails if the binary is missing.

## Differences from Python

* **No pooled ≥2-independent-fragment gate.** If `require_independent_fragments` is truthy, the
  binary exits with an error rather than silently producing different output. Discovery
  (`min_evidence_fragments_per_sample`) now enforces that rule. Every evidence column, including
  `n_independent` and `supported`, is still computed exactly as in Python. See SPEC.md §0.
* Log text on stdout is close to Python's but is not part of the contract.
* Performance only, with outputs unchanged:
  * Chunks are split once they exceed 2000 member loci, so the parsed rows in flight stay bounded
    as sidecars grow. Python uses `n//2000`, capped at 256 chunks.
  * `absorb_one_sided` indexes the rows of the members it can look up in one parallel pass. It
    then re-evaluates batches of independent one-sided loci in parallel (the batch members have
    disjoint candidate windows) and applies the results in Python order.

## Testing

```bash
cargo test --release --manifest-path rust/peartree-combine/Cargo.toml
bash rust/peartree-combine/tests/equiv.sh                 # python vs rust, 4 config variants
THREAD_CHECK=1 bash rust/peartree-combine/tests/equiv.sh  # + threads 4 vs 1
```

`tests/equiv.sh` runs the Python reference and the binary on the 3-colony TPRT E2E fixture with
the `tprt`, `tprt_more`, `tprt_basecfg` and `legacy` configs. It then diffs every output after
decompression (BAMs as sorted SAM). The binary reads the generated `config.py` through
`PEARTREE_PYTHON`, so the `.py` config path is exercised too.
