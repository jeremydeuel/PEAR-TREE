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
- [x] **Multithread BAM decompression.** `pysam.AlignmentFile(..., threads=BAM_THREADS)`
  in `discovery.py`, controlled by `PEARTREE_BAM_THREADS` (default 1 to preserve
  the documented single-core footprint; raise it together with `cpus-per-task`).
- [x] **Rewrite `revcomp`** (`src/revcomp.py`) with a `str.translate` table +
  `[::-1]`. Verified equivalent on 5k random cases; ~18× faster on 150 bp reads.
- [x] **Cut `QualitySeq` allocation churn.** Kept the **list** backing (consensus
  scores exceed the 0-255 phred range, and — measured — an `array('i')` backing
  was ~2× *slower* to build at read length, so `bytes`/`array` were rejected).
  Instead: store the quality list by reference instead of copying it on every
  slice/upper/revcomp result (only copy when handed a non-list such as pysam's
  `array('B')`); short-circuit `upper()`/`lower()` when the sequence is already
  in case; reverse with `[::-1]`. ~1.4× faster on the read-length hot path.
  Verified behaviourally identical to the original across 20k randomised cases
  (incl. out-of-phred-range scores and empty seqs); discovery output and the
  genotype call remain byte-identical.

### 1.2 Kill debug overhead in the hot path
- [x] **Gate `DEBUG`** behind `PEARTREE_DEBUG` (default off) in `discovery.py`,
  `genotype.py`, `genotyping_evidence_read.py`, `breakpoint.py`.
- [x] **Short-circuit `Breakpoint.DEBUG_check_if_breakpoint_of_interest`** when
  debug is off (early `return False`), so the 13-branch coordinate compare no
  longer runs on every `Breakpoint.join`. (Kept the function for now rather than
  deleting, since the coordinates document breakpoints of interest.)
- [ ] Replace remaining per-region `print()` with a `logging` logger.
  _Deferred:_ per-contig prints (~24/run) are cheap; full logging refactor is
  tidier done alongside Stage 2 packaging.

### 1.3 Fix disabled / dead logic (these change results — treat as bug fixes)
- [ ] **`extend_mates()` is a silent no-op.** _Deferred — needs your decision._
  It iterates `self.temporary_breakpoints`, already emptied by
  `extract_chimeric()`'s final `cleanup()`. Pointing it at the final breakpoint
  lists would change emitted clipped consensus sequences; the test BAM is too
  small to judge whether mate-extension helps or hurts real calls. Requires a
  real-WGS benchmark before enabling or deleting.
- [ ] **Re-enable / implement high-coverage exclusion.** _Deferred to Stage 2.1
  (discovery-time masking)._ `max_read_count` is in config/README but referenced
  nowhere in code; genotyping's check is behind `if False and ...`.
- [x] **Fix `except TypeError` fall-through** in `extract_chimeric` — now
  `continue`s on `cigartuples is None` instead of reusing the previous read's
  cigar values.
- [x] **Add a sortedness guard.** `_assert_coordinate_sorted` checks
  `header['HD']['SO'] == 'coordinate'` and fails loudly otherwise.

### 1.4 Dead-code cleanup
- [x] Removed the unreachable discordant-mate branch after
  `continue #dont do any of this`.
- [x] Removed the unreachable tail of `is_good_consensus` after `return True`
  (and its stray debug `print`).
- [x] Removed the dead str-based `clean_clipped_seq` from `adapter.py` (all
  callers use the QualitySeq version in `sequence_checks.py`) and the now-unused
  `is_good_consensus` / `has_well_defined_breakpoint` imports in `discovery.py`.

**Exit criteria:** test BAM call unchanged; discovery wall-time recorded and
improved; no behavioural change except the documented bug fixes.

**Status:** test BAM output byte-identical to baseline before/after; all
discovery-path modules import clean. Two items deferred with rationale above
(`QualitySeq` churn; the `extend_mates` decision — flagged for your call).

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
