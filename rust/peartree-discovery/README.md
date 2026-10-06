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

### CRAM input

CRAM is read natively (via noodles-cram) — the production format. CRAM stores
read bases as differences against the reference, so a reference FASTA is required
and its `@SQ` names must match the CRAM header (use `--reference`):

```bash
peartree-discovery --step discover --bam <sorted.cram> --reference <ref.fa> --out <out.txt.gz>
```

`--bam` accepts either `.bam` or `.cram` (format is chosen by extension). CRAM
uses the single-threaded scan path; the `mate_fetch` (SPD-4) and `contig_threads`
(SPD-3) indexed-fetch optimisations are BAM-only and are ignored for CRAM. Reading
the same alignments from BAM and from CRAM yields byte-identical discovery output.

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
- Config values default to the constants in `src/config.rs` (matching the
  generic `src/config.py`) and can be overridden at runtime — see below.
- BAM decompression uses `--threads N` (or `PEARTREE_BAM_THREADS`) to decode
  BGZF blocks on a worker pool; `1` (the default) keeps the single-core path.

## Runtime configuration (`--config <file>`)

`--config` reads a `key = value` file (`#` starts a comment); env vars override
the file. With no config the run is **byte-identical** to the pre-config port.

| key | default | effect |
|---|---|---|
| `min_mapq` | 40 | anchor MAPQ floor. Generic config = 40, `config_hs/mm` = 60 — pin it on real runs. Env: `PEARTREE_MIN_MAPQ` |
| `min_evidence_reads_per_breakpoint` | 2 | consensus evidence floor |
| `min_good_bases` | 10 | min clipped bases for a QC-passing position |
| `exclude_same_contig_supplementary` | 1000 | supplementary-exclusion distance |
| `cluster_window` / `tsd_min` / `tsd_max` (`max_bp_window`) | 6 / 2 / 40 | clustering + TSD-pairing windows |
| `polya_near_dist` / `polya_far_dist` | 12 / 120 | polyA-rescue proximity band |
| `polya_rescue_min_fragments` | 0 | Poly-A rescue fragment floor (both insertion ends need >= N independent fragments, enforced in discovery). When > 0: (1) the `join()` poly-A rescue of a lone / sub-floor clip cluster needs >= N distinct fragments among the cluster's poly-A reads; (2) a Bp+polyA emission needs >= N distinct fragments among the poly-A-mate reads pooled at that end (`cluster_window`) — on failure the Bp is treated exactly as if no poly-A read were in range; (3) Feature A anchor and sub-floor partner need `n_frags` >= N. Fragment = qname hash (mates / supplementaries of one template count once). Rejections in `stats.json` → `polya_rescue_floor` (block written only when on). 0 = off, byte-identical. Set to 2 in the `.tprt` configs |
| `reject_fully_mapping_reads` | true | XA/SA full-map early reject. Env: `PEARTREE_KEEP_FULLMAP=1` |
| `contig_allowlist` / `contig_allowlist_file` | none | SPEC-5/SENS-4 primary-assembly allowlist (comma list, or one name per line). When set, replaces the `len(name) <= 5` + not-MT heuristic — recovers RefSeq/T2T names like `NC_000014.9` |
| `exclude_bed` | none | SPEC-5 BED of regions whose breakpoints are dropped |
| `coverage_mask` / `coverage_mask_multiplier` | false / 5.0 | SPEC-3 pileup mask: drop breakpoints whose local coverage exceeds N× the genome-wide median |
| `adaptive_evidence` | false | SPEC-4: scale the evidence floor by local/median coverage (never below `min_evidence_reads_per_breakpoint`) |
| `coverage_bin_size` / `coverage_sample_size` | 500 / 3000 | shared SPEC-3/4 coverage-estimator bin size and median subsample bound |
| `evidence_window` | 0 | SENS-1/OBS-3: count support within ±N bp of the modal breakpoint (0 = exact). Gate with SPEC-3/4 |
| `consensus_tolerant` | false | SENS-7: extend consensus while the best base strictly beats the second-best |
| `max_lowq_clip_ratio` / `lowq_mapq_threshold` | none / 40 | SENS-2 guard: drop a locus whose fraction of clipped reads below the MAPQ threshold exceeds the ratio (ship with a lowered `min_mapq`) |
| `hallmark_score` | false | SENS-5: write `<out>.hallmarks.tsv` (poly-A purity, TSD length, EN motif) per insertion. NON-GATING — main output unchanged |
| `short_polya_clip` / `short_polya_min_clip` | false / 7 | SENS-8: allow clips down to N bp when the clipped consensus is a pure poly-A/T terminus |
| `rm_self_mask` / `rm_track` / `rm_divergence_max` | false / none / 5.0 | SPEC-7: drop breakpoints inside a *young* RepeatMasker element (percent divergence ≤ max) from `rm_track` (`.out[.gz]`) — divergence-gated, not family membership. Track must match the BAM assembly; hs1 tracks are in-repo, supply GRCh38/GRCm39 tracks for those |
| `mate_fetch` | false | SPD-4: resolve mates by indexed coordinate fetch (+ SA-tag supplementary loci) instead of a second full pass. Needs a `.bai`. ⚠ byte-identical only on complete BAMs — **validate on the real-WGS differential before enabling** |
| `contig_threads` | 1 | SPD-3: process contigs in parallel across N worker threads (indexed fetch per contig, merged in contig order). Needs a `.bai`. Reproduces single-threaded output on local multi-contig BAMs; **validate on the real-WGS differential before relying on it** |
| `discordant_anchor` | false | **Feature A**: rescue a one-sided junction (e.g. a lone poly-A clip) by pairing a real breakpoint that has no reciprocal partner in its TSD window with a cluster of discordant read pairs on the missing side. Emits a `disc_<pos>` coordinate token (no reads for that end). ⚠ adds calls — validate on the real-WGS differential + orthogonal confirmation before defaulting on |
| `discordant_max_tlen` | 1000 | same-contig template length beyond which a mapped pair counts as discordant (different-contig pairs always do; mate strand is not in the record, so non-FR is not tested) |
| `discordant_min_reads` | 3 | minimum distinct discordant pairs to form an anchoring cluster |
| `discordant_window` | 6 | single-linkage cluster width for discordant observations |
| `discordant_rte_track` | none | **Feature A2**: RepeatMasker `.out[.gz]` scored against the *mate* landing site to establish RTE origin. May be the same track as `rm_track` |
| `discordant_rte_divergence_max` | 20.0 | keep RTE copies with percent divergence ≤ this for the mate-origin test |
| `discordant_rte_only` / `discordant_rte_min` | false / 0.5 | gate: only let a cluster act as a partner when its RTE-origin fraction ≥ min (drops random-SV discordant clusters). Off = label only |
| `splice_hallmark` | false | **Feature B**: write `<out>.splice.tsv` flagging candidates whose mate reads span ≥ `splice_min_exons` exons of one gene with the introns skipped (processed-pseudogene signature). NON-GATING — main output unchanged |
| `exon_annotation` | none | exon annotation (`contig begin end gene_id`, 0-based half-open; `.gz` ok) for `splice_hallmark`. Required when on; must match the BAM assembly (external input, not shipped) |
| `splice_min_exons` | 2 | distinct same-gene exons a candidate's mates must hit to be flagged |

SPEC-3/4 add one lightweight coverage pre-pass over the BAM; it only runs when one
of those gates is enabled, so the default path is unchanged.

**Feature A pipeline note.** `discordant_anchor` is the first discovery toggle that
changes the 4-step contract: a discordant end is written as a `disc_<pos>` token
(mirroring `polyA_<pos>`) with no reads. The Python `combine_insertions` step parses
it (types `TYPE_*_DISC`) and currently **parks** these calls exactly as it parks
poly-A insertions — full downstream genotyping of coordinate-only ends is deferred
pending a dedicated filter. `splice_hallmark` is sidecar-only and touches nothing
downstream.

A `<out>.stats.json` reject-counter sidecar (OBS-1) is written next to every
output, mirroring the per-side `Breakpoint.stats` field set.

[noodles]: https://github.com/zaeleus/noodles
