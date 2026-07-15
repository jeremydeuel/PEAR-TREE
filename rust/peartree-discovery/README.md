# peartree-discovery (Rust port of the discovery step)

A drop-in replacement for **step 1 (discovery)** of PEAR-TREE, reimplemented in
Rust. It reads a coordinate-sorted BAM and emits the **same** custom FASTQ-based
`.txt.gz` breakpoint format as `python src/main.py --step discover`, so it slots
straight into the existing 4-step pipeline (combine_insertions / genotype /
combine_genotypes stay in Python).

Pure-Rust BAM reading via [noodles] — no htslib build dependency, which makes it
easy to build and run on the cluster.

## Build

```bash
cd rust/peartree-discovery
cargo build --release
# binary: target/release/peartree-discovery
```

Requires rustc ≥ 1.87 (noodles versions are pinned accordingly in Cargo.toml).

## Run

```bash
peartree-discovery --step discover --bam <coordinate-sorted.bam> --out <out.txt.gz>
```

The output file is identical in content to the Python discovery output. The
input **must be coordinate-sorted** (the Python step assumes this too; the Rust
port reads the same order).

## Equivalence testing

The Python implementation is the correctness oracle. `tests/differential_test.sh`
runs both and diffs the decompressed output:

```bash
# from the repo root
rust/peartree-discovery/tests/differential_test.sh <bam>
```

Validated byte-identical on:

| input | what it exercises |
|-------|-------------------|
| `test_data/test.bam` | single real insertion: left/right consensus, mate extension, SA/exclude flag, adapter clipping, multi-read consensus |
| `tests/gen_multicontig.py` | multiple contigs + output ordering; `MT` and >5-char contig skips |
| `tests/gen_polya.py` | polyA-rescue path (LEFT polyA paired with a right breakpoint) |
| `tests/gen_xa.py` | the XA/SA "maps fully elsewhere" early filter (on and off) |

Regenerate the synthetic inputs with:

```bash
venv/bin/python rust/peartree-discovery/tests/gen_multicontig.py /tmp/multi.bam
venv/bin/python rust/peartree-discovery/tests/gen_polya.py     /tmp/polya.bam
```

> **Not yet validated on a real WGS BAM.** The Stage 3 exit criterion (see
> `PEAR-TREE2_PLAN.md`) is equivalence on real whole-genome data. Run
> `differential_test.sh` on a real BAM on the cluster before switching the
> pipeline over to the Rust binary. Until then, keep the Python discovery as the
> reference implementation.

## Early "maps fully elsewhere" filter

A clipped candidate read is dropped during discovery if it carries an `XA`
(bwa alternative hit) or `SA` (supplementary) alignment that spans essentially
the **whole** read (clip in that alt `< min_clip_len`). Such a read maps
contiguously elsewhere in the reference and is therefore not a genuine chimeric
junction — step 2 would remove it anyway; doing it here shrinks the candidate
set early, using evidence already in the BAM (no genome required). An alt that
covers only the *clipped* part is the real junction signal and is **kept**.

Implemented identically in the Python and Rust discovery so they stay
byte-identical. On by default; disable with `PEARTREE_KEEP_FULLMAP=1`
(Python honours the same env var, and the `discovery.reject_fully_mapping_reads`
config key).

## Scope & fidelity notes

- The port is a **faithful** reimplementation of the Python discovery, including
  its quirks (e.g. the CLIP_LEFT double-append in `join`, the polyA dict-key
  output ordering, the SA-start vs 0-based-start off-by-one in the exclude
  check). Byte-identical output is the goal; behavioural changes (like the XA/SA
  filter above) are made in **both** implementations at once so they can still
  be diffed against each other.
- `extend_mates()` is a no-op in the Python (it iterates an already-emptied
  list) and is intentionally omitted here.
- Config values (MAPQ, clip lengths, adapters, …) are compiled in from
  `src/config.rs`, matching `src/config.py`. Making them load from a file is a
  Stage 2 item.
- BAM decompression is currently single-threaded; `--threads` / the
  `bam_threads` field are reserved for the Stage 2 multithreaded-decode work.

[noodles]: https://github.com/zaeleus/noodles
