# peartree-rte — Rust port of `tools/rte`

Why: annotate (`tools/annotate_v2.py` → `run_rte` → `RteAnnotator.load_evidence` →
`inputs.read_reads_fa`) loads `<P>.insertions.reads.fa.gz` for every insertion into Python
objects before annotating anything. PD49229 (722 colonies) was killed at 32 GB. This crate
annotates the same records with bounded memory. Work split: [PORT_PLAN.md](PORT_PLAN.md).

Reference behaviour = python `tools/rte` at `eaa2718`. Target: byte-identical `row()` output.

---------------------------------------------------------------------------------------------

## 1. What python does per insertion (annotator.py)

`annotate_v2.run_rte` builds `RteAnnotator(CONFIG['annotate'], gene_model=Insertion.gene_model)`,
loads the sidecars, builds `InsertionInput.from_legacy(ins, element_class(ins.conclusion()))` for
every insertion, then `annotate_all(inputs)`:

1. `_annotate_key(key, inp)`: evidence of the locus (junction rows + combine reads); when the
   genotype_reads file has GT_* reads for it, pool them AFTER the combine reads, annotate, and
   (rte_gt_compare, default True) annotate again with the combine reads only →
   `gt_changed = gt_changes(base, rec)`.
2. `annotate(inp, ev)`: richer junction strings (`_richer_junction`) → `split_junction` →
   `polya_info` (strand hint) → `locate_site` / `target_site` (discovery genome) → SiteContext
   windows (±`local_window` 600, ±`rte_wide_window` 10 kb) → `_cap_reads(reads, rte_max_reads=400)`
   → `Assembler.assemble` → supported tails (`rte_polya_min_fragments`) → poly-A strand logic →
   `polya_reads` → `en_motif` / `slippage_context` → `ExonJunctionIndex.find` (pseudogene candidate
   genes only) → `structure.classify` (NovelSourceFinder, `_premrna_fn` (gene model), pseudogene
   structure fn, legacy class) → site tags (TSD_DELETION, L1_MED_*, EN_INDEPENDENT) → beyond-poly-A
   → support counts → `foldback` → `ScoreInput` → `score` → `RteRecord`.
   Any exception → `RteRecord(key)` with `detail["error"] = f"{type(e).__name__}: {e}"[:200]`.
3. Cross-locus, after all records exist:
   * `cohort_source_pass`: only if `genome_2bit == remap_2bit` and a locator (`remap_index`)
     exists: L1 TPRT/LIKELY_TPRT sites become `novel.cohort_l1`; every record tagged TD3P without a
     `TD3P_SOURCE=` tag is RE-ANNOTATED (reads re-loaded). Off in the current farm configs.
   * `recurrence_pass`: (element, structure, int(j5)//5) shared by > `rte_recurrence_max` (3) L1/
     ALU/SVA loci (not FULL_LENGTH / 5P_UNRESOLVED) → `score_input.recurrent = True`,
     `detail["recurrence"]`, re-score. Needs every record's ScoreInput (small; kept in memory).

What annotate_v2 consumes back (`rte_records`): `rec.row()` (the 17 `RteRecord.COLUMNS`),
`rec.gt_reads`, `rec.gt_changed` (columns only when the genotype_reads file existed),
`'EXON_JUNCTION' in rec.tags` → `ins.exon_junction_proven` (loci with candidate genes), and
`locus_class.locus_class(cls, rec, ann.lib, ...)` which reads `rec.element / tprt_call / tags /
consensus / structure / detail[inv, inv_junction, j5, novel_tier, td_end, exon_junction]` and the
library's `sources` / `flanks3`. `print()` reads element/structure/tags/tprt_call/score/points.
`locus_class.py` and `calibrate.py` stay python (they run on records / the output table).

## 2. Inputs and parity notes

| file | format | Rust |
|---|---|---|
| `<P>.insertions.reads.fa.gz` | `>locus\|SIDE\|ROLE\|sample\|frag\|r12`, locus may contain `\|` (last 5 fields are the tail; < 6 fields padded) | `stream::ReadStore` (spilled, indexed) |
| `<P>.insertions.genotype_reads.fa.gz` | same; only `GT_*` roles kept | `ReadStore` with `role_prefix="GT_"` |
| `<P>.insertions.evidence.tsv.gz` | one row per junction, by header name | `inputs::read_evidence_tsv` (in memory, wanted loci) |

Parity notes (kept on purpose, python is the reference):
* `JunctionEvidence.cross_sample_identical` reads a column of that name; combine writes
  `n_cross_sample_identical`, so it is always 0 → the `cross_sample_identical` score point never
  fires. (Python bug, preserved.)
* combine re-emits the reads of one-sided loci (`…-oneside_N`) several times (PD51635: 5,337 loci
  in 2–17 identical runs, millions of records apart); python appends every copy, so the read
  lists contain duplicates. Preserved (and worth a separate combine fix).
* Header = first whitespace token after `>`; side/role upper-cased; sequence lines stripped and
  concatenated; case of read sequences is NOT kept (every consumer upper-cases).
* `_int(v)` = `int(float(v))` else 0; NaN/inf → 0 (python raises on inf; never occurs).

## 3. Streaming design (implemented: `src/stream.rs`)

Grouping finding: combine writes reads per name in `names` order, but the file is **not**
grouped by locus in practice (duplicated one-sided loci, above), and the genotype2 merge writes
colony by colony (e2e_phylo: 73 of 90 GT loci in > 1 run). A grouped stream is therefore unsafe.

Design: ONE pass over each gzipped FASTA; every record of a wanted locus is appended (4-bit
packed bases + header fields) to an anonymous spill file (`--tmp-dir`, default the output's
directory, unlinked at creation so it never outlives the process); RAM holds only
`locus → [(offset, len)] runs + count`. Per locus, `ReadStore::reads` `pread`s its runs in file
order → exactly python's list. `cap_reads` = `_cap_reads`, applied where python applies it (in
`annotate`, after pooling). `map_loci` runs loci with rayon in bounded chunks; results come back
in input order. The cohort pass re-loads loci from the store (random access).

Measured (`peartree-rte scan`, PD51635 reads FASTA, 269 MB gz / 1.57 Gbases / 10,387,602 reads /
116,637 loci): **peak RSS 120 MB**, 12.5 s wall (11.7 s indexing + 0.6 s reading every locus
back), spill 1.39 GB; max 9,089 reads (1.37 Mb) at one locus. Python held all of it as objects.

Parity check (`tests/golden_foundation.rs::e2e_stream_and_cap_match_python`): for all 336
python `annotate` calls of the e2e_phylo run (90 pooled with GT reads), the store's pooled list
equals python's read for read and `cap_reads` yields python's capped reads in order.

annotate_v2's own whole-file loads (integration): with the Rust engine `read_evidence_clips`
streams the evidence TSV and keeps only the two clip strings of the called insertions (the
parsed sidecar is only built for the python engine), and `read_gt_core` streams the
genotype_reads FASTA (`_iter_gt_fasta`) keeping per locus only the `4 * gt_core_max_queries`
smallest distinct query candidates -- provably the same selection as `gt_query_seqs` over all
reads (one sequence has at most 4 keys: GT_MATE-or-not x side; test
`test_streamed_bounded_selection_equals_gt_query_seqs`).

## 4. Interface (annotate_v2 ↔ binary)

```
peartree-rte annotate --config CFG.json --inputs IN.jsonl[.gz] --out OUT.tsv[.gz]
    [--evidence P.insertions.evidence.tsv.gz] [--reads P.insertions.reads.fa.gz]
    [--gt-reads P.insertions.genotype_reads.fa.gz] [--tmp-dir DIR] [--threads N] [--chunk N]
peartree-rte scan --reads FA [--gt] [--tmp-dir DIR]      # grouping / memory diagnostics
```
Missing sidecar paths are skipped like python. Python side: `tools/rte/rust_bridge.py` —
`write_inputs`, `write_config`, `read_records` (a `RustRecord` per row: `row()`, element,
structure, tags, tprt_call/score/points, consensus, strand, site, detail, gt_reads, gt_changed).

**annotate_v2 switch** (`VariantAnnotationContainer.rte_engine()` / `_run_rte_rust`):
`CONFIG['annotate']['rte_engine']` = `auto` (default: rust when the binary exists, else python
with a loud warning) | `rust` (missing binary = error) | `python`; `rte_binary` = binary path
(default `rust/peartree-rte/target/release/peartree-rte` of the checkout); `rte_threads`
(default `AN_CORES` / `LSB_DJOB_NUMPROC` / the CPU count, bounded by the CPUs available).
Scratch: `mkdtemp` under `TPRT_ANNOT_TMP` (else the run directory) for the inputs JSONL, config
JSON, output TSV and the reads spill (unlinked at creation); removed after a successful run,
kept on failure. A non-zero exit raises (the annotate job fails) -- never a python fallback.

**Input JSONL** (one object per insertion, annotate_v2 order; duplicates are an error):
`{"locus", "left_seq", "right_seq", "pseudogene_genes": [...], "legacy_class": str|null,
"sv": [rank, desc, contig, pos]|null}` = `InsertionInput.from_legacy(ins, element_class(...))`
(only `sv[0]` is read).

**Config JSON** = `CONFIG['annotate']` with callables dropped and `rte_library` absolute
(`write_config`). Keys read (defaults as python; sub-dicts override module DEFAULTS key by key):
`rte_library` (required), `genome_2bit`, `remap_2bit`, `remap_index`, `remap_rmsk`,
`exon_annotation`, `rte_exon_annotation`, `gene_model` (+ `splice_donor_window`,
`splice_acceptor_window`, `splice_ppt_window`, `splice_branch_window`, `promoter_up`),
`young_consensus_regex`, `rte_gt_compare`, `rte_max_reads`, `rte_premrna_window`,
`rte_wide_window`, `rte_max_target_site_deletion`, `rte_polya_min_fragments`,
`rte_recurrence_max`, `rte_assembly`, `rte_structure`, `rte_transduction`, `rte_pseudogene`,
`rte_score{weights,thresholds}`. See `src/config.rs`. NOTE annotate_v2 passes its `GeneModel`
object; the binary loads `gene_model` itself (same file, same windows).

**Output TSV** (header + one row per input, input order), `io::output_columns()`:
`locus`, the 17 `RteRecord.COLUMNS` exactly as python `row()`, `gt_reads`, `gt_changed`
(`.` when empty), `consensus`, `strand`, `site_contig`, `site_L`, `site_R`, `rte_detail_json`
(the typed `detail` dict, for `locus_class` and `exon_junction_proven`). annotate_v2 keeps
writing `gt_reads/gt_changed` only when the genotype_reads file exists (it knows the path).

## 5. Dependencies

* **minimap2**: `minimap2` crate 0.1.31 (bundles minimap2 2.30) through its FFI, replicating
  mappy 2.31's `Aligner.__cinit__` / `map` call sequence (`src/mm.rs`; options per call site in
  its doc). Same pair peartree-combine validated hit-for-hit (616,541 identical hits, its SPEC
  §8). Rejected: calling the minimap2 binary (per-read process / temp files; mappy's in-memory
  `seq=` indices for every locus window) and a pure-Rust aligner (would not reproduce minimap2's
  chaining, the source of every `h.q_st/r_st/mlen/blen`).
* **edlib**: the vendored edlib 1.2.7 C++ (copied from peartree-combine, compiled by `build.rs`
  with `cc`) — the library python-edlib 1.3.9 wraps → identical distances/locations/paths.
  `align::{distance, locate, path}`, wrapped by `sequtil::{edlib_best, edlib_path}`.
* **2bit**: py2bit-semantics reader copied from `rust/peartree-combine/src/genome.rs` (Sync,
  `pread`), with tools/rte's clamping + chr-alias (`src/genome.rs`); FASTA genomes in memory.
  (genotype2's `refseq.rs` pads out-of-range with N — different semantics, not used.)
* `regex` only for a configured `young_consensus_regex` (the default needs a lookahead → coded
  by hand); `serde_json` (`preserve_order`: python dict order matters in golden data), `rayon`,
  `flate2`, `rustc-hash`, `tikv-jemallocator` (default feature).
* Build: `cluster/build.sh` builds this crate as a REQUIRED binary (`hsc_run.sh setup` /
  `submit` check it); crates from crates.io on the head node, or `PT_CARGO_OFFLINE=1` from the
  crates already in CARGO_HOME; `Cargo.lock` committed. Compiles with rustc 1.87.0 (the farm22
  `rust/1.87.0` module; highest dependency `rust-version` is 1.85: hashbrown / indexmap).
  Needs a C and C++ compiler (`cc`: minimap2 2.30 C sources of minimap2-sys, vendored edlib
  C++) and zlib (libz-sys links the system zlib via pkg-config, else builds its bundled copy) --
  the same toolchain peartree-combine already builds with on farm22.

## 6. Ordering and determinism

Output order = input order. Hash maps are lookup-only. Python ordering semantics to keep in
every port: dict insertion order (junctions, detail, Counter keys), `max(..., key=)` /
`Counter.most_common(1)` = FIRST maximal element, `sorted` is stable, `int()` truncates,
`round()` half-to-even on the exact value (`pyfmt::py_round`), float `sum()` is Neumaier on
python 3.12 (`pyfmt::py_sum`), formatting via `pyfmt::{py_repr, py_g}`.

## 7. Golden harness (`golden/`)

`golden/make_golden.py` installs `golden/recorder.py` (monkeypatches tools/rte in place; results
unchanged) and runs either the tools/rte pytest suite (`pytest` mode → `golden/data/pytest.jsonl.gz`,
committed, 1.6 MB, plus `golden/data/libs/` = a snapshot of a tmp library one test builds) or
annotate_v2 on the e2e_phylo simulation (`e2e` mode → 8.4 MB, kept in the session scratchpad;
its `e2e.manifest.json` lists inputs/config/sidecars; tests find it via `PEARTREE_RTE_E2E`, else
skip). Regenerate:

```
PY=<venv with mappy edlib py2bit pysam pytest>
$PY rust/peartree-rte/golden/make_golden.py pytest --out rust/peartree-rte/golden/data/pytest.jsonl.gz
$PY rust/peartree-rte/golden/make_golden.py e2e --out SCRATCH/golden/e2e.jsonl.gz \
    --config SCRATCH/e2e_annot_config.py --tmp SCRATCH/golden/tmp \
    --gt SCRATCH/e2e/ins/P1.insertions.genotype_reads.fa.gz [--cache SCRATCH/annot/new_gt]
```

Event = one JSON line `{"kind", "seq", "case", "parent"?, "locus"?, "lib"?, "genome"?, "remap"?,
"cfg"?, "in", "out"}`; `genome` events define genome ids (`path` (repo-relative when inside the
repo) or inline `regions`). Kinds (count in pytest golden): `hallmark.*` (split_junction 224,
edge_run 452, polya_info 111, locate_site 111, target_site 111, tsd_from_flanks 7, en_motif 109,
slippage_context 109, foldback 112, parse_locus 113, en_bin 7) · `score` 137 · `assemble` 109
(in: ctx, junction_seqs, reads, strand_hint, merged cfg; out: full AssemblyResult incl.
raw/sense layouts) · `classify` 109 (assembly by `assembly_ref`, every callback answer recorded) ·
`known_source` 11 · `novel_find` 15 (locator answers recorded) · `cons_identity` 4 ·
`exon_find` 4 / `exon_structure` 6 · `premrna` 1 (gene-model calls logged) · `annotate` 109
(InsertionInput, junctions, reads, capped read names → record) · `annotate_key` / `final`
(records with `row`). Converters: `src/golden.rs`; tests: `tests/golden_foundation.rs` (on),
`tests/golden_modules.rs` (one `#[ignore]`d test per work package).

## 8. Fidelity risks

1. minimap2 2.30 (crate) vs mappy 2.31 — validated by peartree-combine for the library index;
   tools/rte also uses `seq=` indices (local window, `sr` wide window / pre-mRNA) and `sr` on a
   path (locator): the golden `assemble` events are the check.
2. python object aliasing in `AssemblyResult._pileup` (strand ≥ 0 shares Segment objects
   between `layouts` and `raw_layouts`) — see assembly.rs header.
3. `structure._local_templates(exclude={id(x)})` (identity) vs `e in f5` (equality).
4. `max`/`most_common` tie order, dict order, Neumaier `sum` — §6.
5. Floats compared exactly in the golden tests; identities are `mlen/blen` or `1 - ed/len` with
   the same operands → expected bit-identical.

## 9. Integration results (WP-INT)

Parity:
* `tests/golden_annotate.rs`: all 445 python `annotate` events (109 pytest + 336 e2e_phylo)
  reproduced record-for-record (every field incl. score_input, detail order, row()); classify's
  novel-source / pre-mRNA answers replayed from the recording.
* `golden_e2e_binary`: the binary over the e2e_phylo sidecars (+ genotype_reads) == python's
  final rows + gt_reads / gt_changed (246 loci, 90 with GT reads, 16 calls changed by them).
* annotate_v2 end to end on e2e_phylo, `rte_engine` rust vs python: the annotated table AND the
  verbose stdout report are byte-identical (paths aside), with and without the genotype_reads
  file; the table also equals the golden run's (`e2e.annotated.tsv`, python at eaa2718).
* Real data, PD51635 (hg38, every 17th of the 8,729 annotated loci = 514, python
  `annotate_all` vs the binary on the same inputs/sidecars): 0 differing rows.

Memory / time, PD51635 (local M-series Mac; reads FASTA 269 MB gz, the 8,729 loci of the farm's
annotated table hold 6,629,225 reads; config = config.py.grch38.tprt annotate block, local
hg38.2bit + hs1 rmsk):

| engine | peak RSS | time |
|---|---|---|
| python tools/rte | **4.15 GB** (3.85 GB after `load_evidence`) | 0.90 s/locus wall (514 loci: 461 s) -> ~2.2 h single-threaded for 8,729 (extrapolated) |
| peartree-rte, 1 thread | 0.45 GB (514 loci) | 0.63 s/locus |
| peartree-rte, 8 threads | **0.48 GB** (all 8,729 loci) | 919 s wall, 5,775 s CPU (0.66 s/locus); 10 s indexing |

Python's memory grows with the patient's reads (PD49229: > 32 GB); the binary's is bounded by
the index (~ loci x runs) + the reads of the loci in flight (`--chunk` 64 x <= 400 capped reads).
Time is edlib-bound (`Assembler::layout` rescue, ~70% of samples) like python; the gain is
threads (`AN_CORES`).
