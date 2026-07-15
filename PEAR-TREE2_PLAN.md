# PEAR-TREE2 — Improvement Plan

Branch: `PEAR-TREE2`

Goal: make the **discovery** step dramatically faster and remove more artefacts
earlier, without changing the four-step pipeline contract (each step still reads
and writes the same file formats). Work proceeds in three stages of increasing
effort and risk. Do the cheap, high-ROI Python fixes first; only rewrite in Rust
once the algorithm and filters are settled.

The three stages are independent enough to ship as separate PRs. Each stage ends
with a benchmark on `test_data/test.bam` (correctness) and, where available, one
real WGS BAM (wall-time).

---

## Guiding principles

- **Preserve the on-disk contract.** Discovery must keep emitting the custom
  `.txt.gz` breakpoint format so `combine_insertions` is unaffected.
- **Every change is measured.** Record wall-time and the set of emitted
  breakpoints before/after. A change that alters the breakpoint set must be
  justified as a bug fix, not accepted silently.
- **Correctness gate:** the provided test BAM must keep producing the documented
  insertion call at every step.

---

## Stage 1 — Quick wins in Python (hours–days, low risk)

Pure-Python changes, no architectural shifts. Expected: a meaningful multiple on
discovery wall-time plus better artefact removal, for very little risk. Several
of these are one-liners or deletions of code that currently does nothing.

### 1.1 Speed
- [ ] **Multithread BAM decompression.** Pass `threads=N` to every
  `pysam.AlignmentFile(...)` in `discovery.py` (and genotype/combine). One-line
  change, helps the I/O-bound scan immediately.
- [ ] **Rewrite `revcomp`** (`src/revcomp.py`) using a module-level
  `str.translate` table + `[::-1]` instead of the per-base if/elif chain.
- [ ] **Cut `QualitySeq` allocation churn.** Store quality as `bytes`/`bytearray`
  (as pysam hands it over), slice lazily, and stop calling `.upper()` on data
  that is already uppercase.

### 1.2 Kill debug overhead in the hot path
- [ ] **Gate `DEBUG`** behind an env var / CLI flag; default off. Currently
  hardcoded `DEBUG = True` in `discovery.py`, `genotype.py`,
  `genotyping_evidence_read.py`.
- [ ] **Remove `Breakpoint.DEBUG_check_if_breakpoint_of_interest`** (13-branch
  string compare on hardcoded coordinates, runs on every `Breakpoint.join`).
- [ ] Replace per-region / per-breakpoint `print()` with a `logging` logger.

### 1.3 Fix disabled / dead logic (these change results — treat as bug fixes)
- [ ] **`extend_mates()` is a silent no-op.** It iterates
  `self.temporary_breakpoints`, which `extract_chimeric()`'s final `cleanup()`
  already emptied. Either point it at the final breakpoint lists or delete it.
  Decide with a benchmark on whether mate-extension actually helps calls.
- [ ] **Re-enable / implement high-coverage exclusion.** `max_read_count` is in
  config and the README but referenced nowhere in code. Genotyping's coverage
  check is behind `if False and ...` (`genotype.py`). Decide the intended
  behaviour and wire it in (see Stage 2 for the discovery-time version).
- [ ] **Fix `except TypeError: print(e)` fall-through** in `extract_chimeric`
  (cigartuples `None` case) — should `continue`, not reuse the previous read's
  `left_class/left_len`.
- [ ] **Add a sortedness guard.** Assert the BAM is coordinate-sorted
  (`header['HD']['SO']`) and fail loudly otherwise; the contig-change cleanup
  logic silently produces wrong results on unsorted/name-sorted input.

### 1.4 Dead-code cleanup
- [ ] Remove the unreachable discordant-mate branch after
  `continue #dont do any of this` (~25 lines in `extract_chimeric`).
- [ ] Remove the unreachable tail of `is_good_consensus` after `return True`.
- [ ] De-duplicate the two `clean_clipped_seq` implementations (adapter.py
  str-based vs sequence_checks.py QualitySeq-based).

**Exit criteria:** test BAM call unchanged; discovery wall-time recorded and
improved; no behavioural change except the documented bug fixes.

---

## Stage 2 — Structural Python improvements (days–weeks, medium risk)

Restructure so discovery scales with cores and rejects more artefacts before it
ever writes a candidate.

### 2.1 Parallelism
- [ ] **Contig-level parallelism.** One worker per contig via `fetch(contig)` in
  a `multiprocessing` pool; merge per-contig outputs. Discovery is embarrassingly
  parallel across chromosomes.

### 2.2 Single-pass / targeted mate resolution
- [ ] **Eliminate the second full BAM pass.** `find_mates()` currently rescans
  the entire file to grab mate sequences for a small set of qnames. Replace with
  either (a) a bounded in-memory buffer within the insert-size window during the
  first pass, or (b) targeted `fetch()` around known mate coordinates (available
  from the primary read's flags). This removes ~half the per-read work.

### 2.3 Artefact removal at discovery time (no genome required)
- [ ] **High-coverage region masking.** Running coverage estimate per window;
  drop breakpoints inside pileup spikes. Biggest single artefact-reduction win.
- [ ] **Wire in existing low-complexity / repeat filters.**
  `has_well_defined_breakpoint` and `is_low_complexity` exist in
  `sequence_checks.py` but are not called on the discovery path.
- [ ] **Cluster-level split-read rejection.** Aggregate `SA`-tagged split reads
  that all map locally (cruciform / microindel) at the cluster level, not just
  per read.
- [ ] **Orientation sanity.** Reject right-clip/left-clip pairs whose orientation
  within `max_bp_window` is inconsistent with a TSD-flanked insertion, before
  emitting them.

> Note: the artefact class that needs the reference genome (clipped part maps
> back near the breakpoint) stays in `combine_insertions` — that is correct.

### 2.4 Packaging
- [ ] Turn `src/` into an installable package (`pyproject.toml`, console entry
  point) so imports stop depending on the working directory and the broken
  `peartree` launcher shebang.
- [ ] Move config out of a gitignored `config.py` full of hardcoded absolute
  cluster paths into a `--config <file.yaml>` argument. Keep species templates.

**Exit criteria:** discovery scales roughly with core count; second BAM pass
removed; artefact count in combine measurably lower for equal true-positive
recall on the test/validation set.

---

## Stage 3 — Rust rewrite of discovery only (weeks, higher effort)

Once the algorithm and filters are settled in Stage 2, port **discovery only**
to Rust. Do not rewrite combine/genotype/annotate — they are subprocess- and
I/O-bound glue where Python is fine and where the fiddly, still-evolving
filtering logic lives.

- [ ] Implement discovery with [`rust-htslib`](https://github.com/rust-bio/rust-htslib)
  or [`noodles`](https://github.com/zaeleus/noodles) + [`rayon`](https://github.com/rayon-rs/rayon)
  for per-contig parallelism.
- [ ] **Emit the exact existing `.txt.gz` breakpoint format** so the Rust binary
  is a drop-in replacement for step 1 in the current pipeline.
- [ ] Port the Stage-2 discovery-time artefact filters.
- [ ] Ship as a standalone binary invoked by the existing orchestration.
- [ ] **Equivalence test:** Rust and Python discovery must produce the same
  breakpoint set (modulo documented improvements) on the test BAM and at least
  one real WGS BAM before switch-over.

Expected: ~20–50× over the current Python discovery.

**Exit criteria:** Rust discovery matches Python output on the validation set and
is adopted as the default step-1 implementation; Python discovery kept as a
reference/fallback until confidence is high.

---

## Sequencing & rollback

- Ship Stage 1 → Stage 2 → Stage 3 as separate PRs; each is independently
  useful and revertible.
- Keep the Python discovery path working through Stage 3 as the correctness
  oracle for the Rust port.
- Do not begin Stage 3 until Stage 2's filters are stable — porting a moving
  target to Rust wastes the rewrite.
