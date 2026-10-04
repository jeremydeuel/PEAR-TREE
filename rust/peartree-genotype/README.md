# peartree-genotype — Rust port of the PEAR-TREE genotyping step

A drop-in replacement for

```
python src/main.py --step genotype --bam <bam> --insertions <ins.genotyping.txt.gz> --out <out.txt.gz> [--threads N]
```

emitting the same 8-column gzip table. The **decompressed** output is byte-identical
to the Python genotyper under the generic `src/config.py` genotyping defaults (the
gzip *bytes* differ — compression level and header mtime only).

```
peartree-genotype --step genotype \
    --bam sample.bam \
    --insertions cohort.genotyping.txt.gz \
    --out sample.genotypes.txt.gz \
    [--threads N] [--config genotyping.cfg] [--reference ref.fa]
```

Input may be **BAM or CRAM** (indexed: `.bai`/`.csi`, or `.crai`). CRAM whose bases
are stored against an external reference needs `--reference ref.fa` (names matching
the CRAM `@SQ`); self-contained / embedded-reference CRAM does not consult it.

### Batch mode (many samples, one run)

A phylogenetic tree genotypes the same contract against hundreds of colonies.
`genotype_batch` parses the contract once and processes a manifest of files
**consecutively** — one file at a time, with `--threads` cores spread across that
file's loci. This is deliberate for the cluster: 4 cores over 20 files never opens
20 readers at once, so memory and file descriptors stay bounded and predictable.

```
peartree-genotype --step genotype_batch \
    --manifest samples.tsv \
    --insertions cohort.genotyping.txt.gz \
    [--threads N] [--config genotyping.cfg] [--reference ref.fa]
```

`samples.tsv` has one `input<TAB>output` line per sample (BAM or CRAM; blank / `#`
lines ignored). Every input is validated up front, so a bad path fails before any
work. Each output is byte-identical to the corresponding single-file run.

## What it does

For each locus in the contract: count spanning depth (the high-coverage gate),
fetch the spanning reads via the BAM index, apply the same read gates and
per-fragment dedup as `Insertion.genotype`, score each read on both junctions
(`query_index_at_ref` CIGAR walk → ±1 bp `qscore` register search with the
perfect-match short-circuit), and summarise to a call with the VAF band model
(`summarise_evidence`). Loci are independent, so the work is split across
`--threads` workers (each with its own indexed reader) and reassembled in contract
order — the output is thread-invariant.

The genotype vocabulary, the raw-evidence columns (`coverage n_alt n_ref n_art`),
`min_reads_for_zygosity`, and `recover_low_coverage_presence` all match the Python
implementation.

## Layout

| file | port of |
|------|---------|
| `evidence.rs` | `genotyping_evidence_read.py` + `genotype_qscore.py` (`query_index_at_ref`, `qscore`, `qleft`/`qright`) |
| `insertion.rs` | `genotyping_insertion.py` (`Insertion.import_file`, `_side_call`, `summarise_evidence`) |
| `read.rs` | per-record decode over `sam::alignment::Record` (shared BAM/CRAM view, like the discovery crate) |
| `source.rs` | `RegionSource` trait unifying indexed BAM + CRAM region queries |
| `genotype.rs` | `genotype.py` driver (region fetch, gates, dedup, threading, 8-column gzip writer) |
| `config.rs` | the `genotyping` block of `src/config.py` (key=value file + `PEARTREE_MIN_MAPQ`) |

## Configuration

Defaults mirror the **generic** `src/config.py` (`min_mapq = 40`,
`reads_for_high_coverage = 180`, …). Production configs use `min_mapq = 60` —
pin it with `--config` or `PEARTREE_MIN_MAPQ` on real runs. The `--config` file is
the same `key = value` format as the discovery crate, so one file can carry both
sections.

### TPRT locus kinds (all keys default off -> byte-identical legacy output)

Discovery's TPRT pairing modes (`plans/tprt_hallmarks/SPEC.md`) name loci whose geometry
differs from a TSD (`R - L` in 2..40). Per-read scoring is per junction (`qleft` scores the
read 5' of L against `LEFT_REFERENCE`/`LEFT_INSERTION`, `qright` the read 3' of R), so a
target-site deletion (`R < L`, `genome[R, L)` absent from the alt allele) and a blunt
junction (`R - L` in {0, 1}) genotype correctly unchanged. `cluster/config.genotype.grch38.tprt`
turns on what the other kinds need:

| key | locus kind | what it does |
|---|---|---|
| `one_sided_loci` | `contig:L-oneside_L`, `contig:oneside_R-R` | parse the `oneside_` token, score only the real junction (off: such a locus is an `error` row) |
| `one_sided_open_window`, `one_sided_open_min_clip` | one-sided | skip the missing end's junction reads (soft clip on the open side within the window of the real breakpoint) — otherwise they are false ref votes |
| `split_breakpoint_span` | L1-mediated deletion / duplication (`\|R - L\|` up to 50 kb) | depth-gate and fetch each breakpoint as its own window (else the whole span counts as depth -> `high-coverage`) |
| `halve_single_junction_ref` | one-sided, split far pairs | each ref vote counts half: the VAF bands assume alt from two junctions per reference span |
| `dup_ref_discount_min_span` | far duplication | `n_ref -= min(n_ref, n_alt)`: the alt haplotype keeps both reference junctions |

The contract for one-sided loci comes from `src/genotyping_contract_oneside.py` (combine
leaves them out of `<patient>.genotyping.txt.gz`); `cluster/pipeline.sh` builds
`<patient>.genotyping.tprt.txt.gz` when the genotype config sets `one_sided_loci = true`.
Validation: `test/e2e/run_genotype_e2e.sh`, `plans/tprt_hallmarks/E2E_REPORT.md`.

## Scope / limitations

- **BAM and CRAM**, both indexed (BAM: `.bai`/`.csi`; CRAM: `.crai`). CRAM decoding
  goes through `source.rs`, which unifies the two backends behind one
  `RegionSource` trait (both yield records implementing `sam::alignment::Record`).
- The `.txt.gz` contract is the equivalence boundary: any change to the emitted
  table is a bug, not an optimisation.

## Testing

```
cargo test                                   # unit tests (evidence + all summarise bands)
tests/differential_test.sh <bam> <ins.genotyping.txt.gz> [python] [rust-bin] [threads]
```

The differential harness is the equivalence oracle: it runs both implementations
on the same BAM + contract and diffs the decompressed tables. Validated
byte-identical (and thread-invariant) on `test_data/test.bam` across the golden
locus, a 282-locus scan exercising no-coverage / artefact / heterozygous / dedup,
and the high-coverage path (lowered `reads_for_high_coverage`). CRAM was validated
the same way against a self-contained (embedded-reference) copy of the test data:
Rust-CRAM == Python(pysam)-CRAM == Rust-BAM == Python-BAM. Batch mode was validated
by running 20 mixed BAM/CRAM samples on 4 threads and confirming every output
matches the single-file oracle.
