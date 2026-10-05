# peartree-combine — implementation plan

Phase 1 (done, architect): crate skeleton that compiles, `SPEC.md`, this plan, the
equivalence harness `tests/equiv.sh`, and the three foundation modules everyone uses:
`align.rs` (vendored edlib 1.2.7 FFI, parity-tested 3000/3000 against python-edlib 1.3.9),
`pyfmt.rs` (CPython-3.12 Neumaier `sum`, half-even `round`, dedup budget, `_median`),
`model.rs` / `context.rs` (core types, interner, locus tokens).

Phase 2: six work packages in parallel (below). Phase 3: integration + equivalence + scale run.

## Rules for every package

* Work from `SPEC.md` + the Python file(s) named in your module's doc comment. The Python is
  normative; if SPEC is wrong, fix SPEC in the same commit.
* Only edit the files you own. Shared signatures (anything `pub` in the skeleton) change only by
  agreement — note it in your commit message so the others can rebase.
* Port literally first (same control flow, same iteration order); optimise only where the
  plan says so, and never in a way that changes iteration order or float evaluation order.
* Unit tests live in your files (`#[cfg(test)]`); for numerical code, generate expected values
  from the Python with the reference venv
  (`/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad/venv/bin/python`,
  has edlib/mappy/pysam/pyliftover/py2bit) and check them in as literals.
* `cargo build --release` must stay green; commit with `git commit --no-gpg-sign`.

## Work packages

| WP | files owned | Python | model | depends on | can start |
|---|---|---|---|---|---|
| **P1 core I/O + intersect** | `config.rs`, `seq.rs`, `insertion.rs`, `intersect.rs`, `region_filter.rs`, `genome.rs`, `model.rs` (only the two `todo!()` consensus helpers) | INS, ISX, RF, GS, QS, SC:89, TF:316-349, config | sonnet | – | now |
| **P2 consensus** | `consensus.rs` | IC (all) | **opus** | align, pyfmt (done) | now |
| **P3 junction evidence** | `evidence/row.rs`, `evidence/dedup.rs`, `evidence/junction.rs` | EV:80-148, 178-755, 1265-1305 | **opus** | P2 for `evaluate_junction` tests (signature fixed) | now |
| **P4 TPRT filters** | `tprt.rs`, `library.rs`, `evidence/filters.rs` | TF:45-311, EV:358-363, 1294-1452 | **opus** | P1 `seq`/`genome` for integration tests only | now |
| **P5 evidence orchestration** | `evidence/mod.rs`, `evidence/store.rs`, `evidence/output.rs` | EV:150-176, 766-1243, 1255-1331 | **opus** | P3, P4 (signatures fixed; stub them in tests) | now (store/output first) |
| **P6 remap, outputs, driver** | `remap.rs`, `liftover.rs`, `genotyping_out.rs`, `splice.rs`, `pipeline.rs`, `main.rs` | CI (all), main.py:148-188, pyliftover | sonnet | P1 for an end-to-end run | now (liftover/remap/CLI first) |

Integration order: P1 → (legacy variant runs end-to-end once P6 is in) → P2+P3 → P4 → P5 →
TPRT variants. The `legacy` harness variant needs only P1 + P6; it is the first milestone.

### Hot spots per package

* **P1** — `intersect_insertions` output order is the order of every output file (SPEC §3.3).
  Use an insertion-ordered map (Vec of entries + `FxHashMap<key, index>`); string keys of disc reps
  come after all tuple keys. The parser must not allocate per-record Strings (intern contigs,
  `Tok` tokens, boxed byte slices) and must skip MATE records without storing them. `genome.rs`:
  positional reads (`FileExt::read_at`) so the handle is `Sync`; N blocks must be honoured.
* **P2** — the highest fidelity risk in the port. Small insertion-ordered maps for `_Col.w`
  (keys A,C,G,T,N,'-') and `_Col.ins` (RLE strings); `py_sum` for `others`; `py_round`;
  identical edlib calls (wildcard pad = `"N" * len(read)`). Test against `test/test_indel_consensus.py`
  cases + randomised reads diffed against the Python (dump a JSON of reads → consensus fields).
* **P3** — `independent_clusters` is O(n²) per sample per junction (capped at 200 fragments by
  discovery): keep it literal; avoid allocation in `_seq_close` (reuse buffers). The evaluation of
  one junction must be a pure function of its rows (it runs on rayon workers).
* **P4** — `LibraryMatcher` via the `minimap2` crate (feature `mappy`): set
  `idxopt.k=11, w=3`, `mapopt.min_cnt=1, min_chain_score=15, min_dp_max=20, best_n=3`, CIGAR on;
  iterate hits in minimap2 order; strictly-greater mlen wins. Verify hit-for-hit against mappy on a
  few thousand clips taken from the fixture's evidence TSV (SPEC §8 risk 1). The matcher must be
  usable from several threads (one thread buffer per thread; a cache is optional).
* **P5** — see "Memory / parallelism design" below; results must be identical for any chunking
  and thread count (`THREAD_CHECK=1 tests/equiv.sh`). `absorb_one_sided` must use a positional
  index (Python is O(one-sided × all)).
* **P6** — command strings byte-identical to Python (they are recorded in the BAM @PG header,
  which the harness ignores, but keep them identical anyway); skip a bowtie2 run when its BAM
  exists (resume semantics); read BAMs via `samtools view` text. `liftover.rs` must load the
  hs1→hg38 chain (~ tens of MB gz) quickly; an interval index per source contig.

## Memory / parallelism design

Target patient: 44–100 colonies, ~760k discovery records, sidecars of tens of GB (Python today:
15.5 GB legacy peak on 174 files; TPRT unmeasured but larger). Target peak: a few GB.

**Discovery records (held).** Parsed in parallel, one rayon task per file, into compact
`Insertion`s: interned contig (u32), two `Tok`s (16 B each), up to four `QualSeq`s as boxed byte
slices (seq + qual), `files`/`member_loci` small vecs; MATE records skipped (in Python they are the
bulk of the discovery FASTQ). ~0.7 KB/record → ~0.5 GB for 760k records. After intersect the
merged survivors replace the input (inputs consumed by `merge_into` / dropped), typically
200–400k.

**Sidecar rows (streamed, never all in memory).** `Store::build`:
1. Build `FxHashMap<(FileId, LocusKey), Vec<u32 chunk id>>` from the insertions' member loci
   (~1–2 M entries × ~48 B → < 100 MB).
2. Chunks = consecutive insertion index ranges of ~2000 insertions (Python's rule
   `max(4*threads, n//2000)` capped at 256 is fine; any choice is legal).
3. One rayon task per sidecar file streams it (gzip → lines), parses only the `locus` column
   first (skip unwanted rows cheaply), and appends wanted rows (raw line bytes, or a compact binary
   encoding) to per-chunk in-memory buffers that are flushed (when > ~1 MB) as blocks into ONE
   scratch file per sidecar, recording `(chunk → [(offset, len)])`. Open handles = one per file
   being read; no per-(file × chunk) files.
4. `chunk_rows(c)`: read every file's blocks of chunk c (positional reads), parse into
   `EvidenceRow`s, group by (member, side) keeping sidecar order. Memory per chunk in flight:
   ~2000 insertions × ~50 rows × ~0.4 KB ≈ 40 MB; with `threads` chunks in flight < 1 GB.
5. `lookup(m, side)` for `absorb_one_sided`: read the first chunk holding m (small LRU of parsed
   chunks, like Python's 4-entry cache).

**Records (held, compact).** After a chunk is judged, each `JunctionRecord` keeps its counters,
both `ConsensusResult`s (needed for the TSV and for `_replace_clips` in absorb) and `aligned`;
its rows are rendered to FASTA text, written to a per-chunk reads scratch file
(`put_reads` → `ReadsRef{chunk, offset, len}`), and dropped. ~1.5 KB per junction → ~1 GB for
600k junctions. (If that is too much: render the TSV row text eagerly and keep only
`n_independent` + the combined consensus.)

**Outputs (streamed).** evidence.tsv rows and reads are written in final name order, reading
reads text back by `ReadsRef` (sequential within a chunk, so mostly sequential I/O).

**Parallel stages.** discovery parsing (per file), sidecar routing (per file), chunk evaluation
(per chunk, `rayon::ThreadPool` sized to `--threads`), bowtie2 (`--threads`). Order-dependent
stages (import concatenation, intersect, region filter, remap filters, absorb, outputs) are
sequential and linear-time.

**Expected peak for 760k records / 44 colonies:** discovery records ~0.5 GB (freed progressively
by intersect) + member index 0.1 GB + chunks in flight ≤ 1 GB + records ~1 GB + liftover chain
~0.2 GB ≈ **2–3 GB**, vs tens of GB in Python. Scratch disk ≈ the wanted sidecar rows
(uncompressed or lz-light) + rendered reads, removed at the end.

## Phase 3 — integration (one agent, opus)

1. Remove `#![allow(unused…)]` from lib.rs, fix warnings.
2. `bash rust/peartree-combine/tests/equiv.sh` → all four variants IDENTICAL, then
   `THREAD_CHECK=1` (threads 4 vs 1).
3. Scale: run Python and Rust on one real patient on the farm (e.g. the PD37590 TPRT arm) and diff
   outputs; record wall time and peak RSS; then switch `cluster/combine_mei.sh` / `pipeline.sh`
   behind a `COMBINE_IMPL=rust` toggle.
