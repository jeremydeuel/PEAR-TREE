# PEAR-TREE discovery — speed recommendations (Rust port)

*Scope: the **Rust** discovery implementation only —
`rust/peartree-discovery/` on branch `PEAR-TREE2`. Goal: make step 1 as fast as
possible **without changing the emitted `.txt.gz` breakpoint contract** (the port
is currently byte-identical to the Python oracle on the differential tests, and
must stay diffable against it). The Python discovery is out of scope here except
as the correctness oracle.*

*Every line/file reference below is to `rust/peartree-discovery/src/…`.*

---

## 0. TL;DR — where the Rust time goes and what to do

The port is a faithful single-threaded, single-pass-per-phase translation of the
Python control flow. Its cost, for a 30× WGS BAM (~1 billion reads), is dominated
by, in rough order:

1. **BGZF decompression of the whole file, done twice** — `extract_chimeric`
   (`discovery.rs:177`) and `find_mates` (`discovery.rs:344`) each stream the
   entire BAM via `build_from_path` (single-threaded decode).
2. **Per-read heap allocation in `BamRead::from_record`** (`read.rs:54–137`) —
   `seq: Vec<u8>`, `qual: Vec<i32>`, `query_name: String`, `reference_name:
   String`, and the SA/XA tag strings are all built for **every** record, then
   most are immediately discarded by cheap gates.
3. **No parallelism** — `bam_threads` is stored but never used (`discovery.rs:81`;
   README line 92–93); discovery runs on one core.

Everything downstream (clustering, `join`, consensus, output) is cold-path:
breakpoints are rare relative to reads, so `model.rs` / `filters.rs` are not worth
optimising. **All the wall-time is in the read loop and the decode.**

| # | Change | Est. gain | Effort | Risk |
|---|--------|----------|--------|------|
| **R1** | **Lazy / cheap per-read decode** — reorder `from_record` so heavy fields are built only for surviving clipped candidates | **large constant factor** (skips alloc on ~99% of reads) | Medium | Low |
| **R2** | **Multithreaded BGZF decode** — wire the unused `bam_threads` to noodles' multithreaded reader | **~2–4×** (I/O-bound) | Low | Low |
| **R3** | **Contig-level parallelism** with `rayon` + an indexed reader | **~N cores** | Medium | Medium |
| **R4** | **Eliminate the 2nd full BAM pass** (`find_mates`) via mate-coordinate fetch | **~1.5–2×** | Medium | Medium |
| **R5** | **Cheap internals** — integer contig compare, `u8` quals, faster hasher, `Record` reuse | **1.2–1.5×** stacked | Low–Med | Low |

R1, R2 and R5 are single-thread constant-factor wins that compose; R3 and R4
change the I/O structure. **Measure first** (§6) — the "20–50× over Python"
estimate in the Plan is on data too small to trust, and no real-WGS benchmark
exists yet.

---

## R1 — Lazy per-read decode (biggest single-thread win)

`BamRead::from_record` (`read.rs:54–137`) eagerly materialises, for **every**
record, fields most reads never need:

- `seq: Vec<u8>` (`read.rs:91`) — heap alloc, ~150 bytes.
- `qual: Vec<i32>` (`read.rs:92`) — heap alloc **at 4 bytes/base** (see R5a), built
  by mapping each `u8`→`i32`.
- `query_name: String` (`read.rs:94–97`) — heap alloc.
- `reference_name: String` (`read.rs:60–62`) via `from_utf8_lossy(...).into_owned()`
  — heap alloc **per read** just to compare against the current contig.
- `sa` / `xa` tag strings (`read.rs:99–107`) — `value_to_string` allocates.

But look at the gates that run first in `extract_chimeric` (`discovery.rs:185–233`):
the MAPQ gate, the secondary/qcfail/duplicate flag gates, the contig-name/`MT`
skip, the `has_cigar` check, and the clip-shape classification. **None of them need
seq, qual, name, or tags** — only flags, mapq, reference id, and the first/last
CIGAR op. On a typical BAM the overwhelming majority of reads are dropped by these
gates (not clipped, or duplicate, or low-MAPQ), so all that allocation is wasted.

**Fix: split decoding into two tiers.**

1. A cheap `from_record_header` that fills only `flags`, `mapq`, `reference_id`
   (integer), and the first/last CIGAR op (`kind`, `len`) — noodles exposes CIGAR
   ops without materialising seq/qual.
2. Run the gates. Only for a surviving **clipped candidate** (the `Some(clip)` arm,
   `discovery.rs:233`) call a second step that decodes `seq`, `qual`,
   `query_name`, and — only when needed — the SA/XA strings.

Note the low-MAPQ poly-A branch (`discovery.rs:185–190`) needs `seq` too, but only
for reads below `MIN_MAPQ` **and** not in a proper pair (`find_polya` bails
immediately on `is_proper_pair`, `polya.rs:132`) — so decode seq lazily there as
well, after the `is_proper_pair`/`mate_is_mapped` pre-checks, which only need
flags.

This removes per-read heap allocation from ~99% of records. It is the single
biggest constant-factor improvement available in the current Rust code and it is
pure mechanical refactoring — no algorithm change, output stays byte-identical.

---

## R2 — Multithreaded BGZF decompression

BGZF decode is the dominant CPU cost of a linear BAM scan and parallelises
trivially at the block level. The port already reserves the knob but does nothing
with it: `bam_threads` is stored in `Discovery` (`discovery.rs:81–82`), parsed from
`--threads` / `PEARTREE_BAM_THREADS` in `main.rs:31–42`, and then **never read**.
Both readers use `bam::io::reader::Builder::default().build_from_path(...)`
(`discovery.rs:177`, `:344`), which is single-threaded decode.

**Fix:** build the reader over noodles' multithreaded BGZF reader
(`noodles_bgzf::MultithreadedReader`, worker count = `bam_threads`) and wrap it in
`bam::io::Reader::new(...)`, or use whatever multithreaded-decode constructor the
pinned noodles version exposes. Add `noodles-bgzf` to `Cargo.toml` if not already
transitively pinned. This is a small, low-risk change with a large I/O-bound
payoff, and it finally makes `--threads` mean something.

Decode-thread scaling on a single file saturates (block reorder + downstream become
the bottleneck) around 4–8 workers — beyond that, prefer R3. See §5 for splitting a
core budget between R2 and R3.

---

## R3 — Contig-level parallelism (rayon)

Discovery is embarrassingly parallel across chromosomes: the scan, `cleanup`
clustering, and `join` consensus for one contig are independent of another. This is
the deferred item in the Plan (Stage 3: "`rayon` per-contig parallelism deferred").

**Fix:**
- Switch from the streaming `build_from_path` to an **indexed** reader
  (`bam::io::indexed_reader`, needs `.bai`/`.csi`), and `rayon`-map over the
  main-chromosome list. Each worker opens its own reader handle + file descriptor
  and `query()`s exactly its contig region.
- Collect per-contig `(final_left, final_right, polyA)` and concatenate. `output`
  already groups and sorts by reference name (`discovery.rs:387–438`), so the merge
  is a concatenation followed by the existing per-contig sort — no new logic.

**The one coupling is mate resolution (R4).** Mates can be on another contig, so
contig-parallel pass 1 must emit its breakpoints **plus** the list of needed mate
coordinates, and a second phase resolves them. This is exactly why R4's
mate-*coordinate* design should land with R3: the current qname-scan `find_mates`
does not parallelise per-contig, but a coordinate-driven fetch does (group the mate
targets by contig → parallel).

Expected scaling is near-linear in cores until shared-file decode or the merge
dominates. On a 16–32-core node this is the hours→minutes step.

---

## R4 — Eliminate the second full BAM pass

`find_mates` (`discovery.rs:342`) opens and streams the **entire BAM a second
time** (`discovery.rs:344`) only to pull mate sequences for the small qname sets
built by `get_mates`. That is roughly half of all discovery decode spent on a
lookup.

**Key fact:** the anchor read already carries its mate's location in `RNEXT`/`PNEXT`
(`next_reference_id` / `next_reference_start`), but `BamRead` never decodes it and
`find_mates` matches on `qname` instead. Decode the mate coordinate in pass 1, and
the second scan becomes a handful of short indexed queries.

**Fix (preferred): targeted fetch.**
1. Add `mate_reference_id: Option<usize>` and `mate_alignment_start: Option<i64>`
   to `BamRead` (cheap — no seq/qual needed).
2. When a read becomes a breakpoint (or a poly-A), record
   `(mate_ref, mate_pos, qname, target-slot)` — the slot is the `BpRef`
   (`Left/Right/PolyA(idx)`) that `get_mates` already computes.
3. After pass 1, sort targets by `(ref, pos)`, coalesce nearby ones, and `query()`
   only those regions with the indexed reader. Fetched reads still pass the same
   `is_secondary/qcfail/duplicate` and read1/read2 gates (`discovery.rs:349–358`)
   before feeding `bp.mate_seqs` / `set_mate`.

Mates that are unmapped (placed next to their mate in a coord-sorted BAM) or
interchromosomal both fall out naturally from fetching at `(mate_ref, mate_pos)`.
Output stays byte-identical because the same mate reads are visited — only the way
they're located changes. Composes directly with R3.

> If random-access seek overhead shows up in the profile, fall back to a bounded
> in-pass buffer (keep reads within one insert-size window during pass 1; most
> mates of proper pairs are already resident) — but start with fetch; it's simpler
> to keep correct.

---

## R5 — Cheap internal wins (stack these)

Small, safe, single-thread constant-factor improvements. Each is minor alone; they
compound and are mostly local.

- **R5a. `reference_name` as an integer compare.** Today every read allocates a
  `String` for `reference_name` (`read.rs:60–62`) purely to compare against the
  current contig and run the `MT` / `len > 5` skips (`discovery.rs:198–207`).
  Carry the integer `reference_sequence_id` instead; resolve the name to a `String`
  **once per contig change** (the `cleanup` boundary). Removes one heap allocation
  from every read. Pairs with R1.
- **R5b. Store raw qualities as `u8`, widen only at consensus.** `qual: Vec<i32>`
  (`read.rs:92`, `qseq.rs:22`) is 4 bytes/base and requires a `u8`→`i32` map on
  every kept read. Raw phred fits in `u8`; only `find_consensus`' *scores*
  (`filters.rs:117`) exceed 255. Keep read/`QualitySeq` quals as `u8` and use a
  separate `i32` score track only for consensus output. Cuts quality memory traffic
  ~4× on the hot path. Medium effort (touches `QualitySeq`), so gate it behind a
  differential re-test — but it's a real win given quals are allocated per kept
  read.
- **R5c. Faster hasher for the qname maps.** `get_mates`/`find_mates` use
  `HashSet<String>` / `HashMap<String, BpRef>` (`discovery.rs:298–301`) with the
  std default SipHash. qnames aren't adversarial input — swap to `FxHashMap` /
  `ahash`. One-line type aliases; measurable on large breakpoint sets.
- **R5d. Reuse the record buffer.** `reader.records()` (`discovery.rs:181`, `:346`)
  yields an owned `Record` per iteration. Use `reader.read_record(&mut record)`
  into a single reusable `bam::Record` to cut per-iteration allocation across the
  ~billion-read loop.
- **R5e. Skip non-main contigs at the reader, not after decode.** With the indexed
  reader from R3, simply **don't iterate** alt/decoy/unplaced/HLA contigs and `MT`
  at all, rather than decoding each read and discarding on
  `ref_name.len() > 5`/`== "MT"` (`discovery.rs:202–207`). On GRCh38-style BAMs the
  decoy+alt reads are a non-trivial slice you currently pay full decode for. (This
  also matches the "explicit main-chromosome allowlist" correctness recommendation
  in `RTE_detection_review/08…md` R5.)

Not worth touching (cold path, correctness-sensitive, byte-identical constraints):
`most_common_first`'s O(n²) scan (`model.rs:60` — n = reads in a <6 bp cluster,
tiny), the `join`/consensus clones (`model.rs`, `filters.rs`), the naive substring
search in `is_adapter`/`find_polya` (per-candidate only), and the output gzip level
(output is small).

---

## 6. Measure first — a Rust-specific profiling plan

Nothing here is benchmarked on real WGS. Before implementing, establish the
baseline so each change is attributable:

1. **Baseline wall-time + read count.** `samtools flagstat` for the read total;
   `/usr/bin/time -v target/release/peartree-discovery --step discover …` on one
   real WGS BAM. Record RSS (R3 raises it; R5b lowers it).
2. **Attribute the two passes.** `extract_chimeric` and `find_mates` are separate
   calls in `discovery()` (`discovery.rs:379–385`); time them individually. The
   ratio is R4's ceiling.
3. **Confirm the alloc hypothesis (R1/R5).** Build with a profiler
   (`cargo flamegraph`, or `perf record` + `perf report`) on a single-chromosome
   subset (`samtools view -b in.bam chr20`). Expect `from_record` allocation +
   BGZF decode to dominate. If the flamegraph shows most time in `Vec`/`String`
   allocation and `from_utf8_lossy`, R1/R5a/R5b are confirmed high-value; if it's
   almost entirely in `inflate`/BGZF, R2 is the lever.
4. **Decode-thread sweep (R2).** Once wired, run `--threads` 1/2/4/8 on the same
   BAM to find where decode scaling saturates — that sets the R2-vs-R3 core split.
5. **Real-WGS differential (correctness gate).** Run
   `tests/differential_test.sh <real.bam>` after **each** change and confirm
   byte-identical output vs the Python oracle. This is also the still-outstanding
   Stage 3 exit criterion (README: "Not yet validated on a real WGS BAM").

---

## 7. Recommended sequencing

Ship as independent, revertible changes, re-running `differential_test.sh` each
time:

1. **Measure** (§6) — baseline + flamegraph + pass split. *(hours)*
2. **R2** — wire `bam_threads` to the multithreaded reader. *(low effort, immediate
   I/O win, no output change)*
3. **R1 + R5a + R5d** — lazy decode, integer contig compare, record reuse. All
   local to `read.rs`/`discovery.rs`, all byte-identical. *(the big single-thread
   constant factor)*
4. **R5c** — faster hasher for the mate maps. *(one line)*
5. **R4** — mate-coordinate fetch; needs the indexed reader, so pull that in here.
   *(medium; removes the 2nd pass)*
6. **R3** — `rayon` per-contig parallelism over the indexed reader, plus **R5e**
   (skip non-main contigs by not iterating them). *(medium; scales with cores)*
7. **R5b** — `u8` quals with a separate consensus score track. *(do last; touches
   `QualitySeq`, so isolate it behind its own differential re-test)*

**Guardrails:** preserve the `.txt.gz` contract; every change stays byte-identical
to the Python oracle on `test_data/test.bam` **and** the real-WGS differential
(the documented call at `13:32992169-32992177` must survive); keep the Python
discovery as the reference/fallback until the Rust path is validated on real WGS.
Any change that alters the emitted breakpoint set is a bug, not an optimisation.

## 8. Out of scope

- Artefact-removal changes (SMS reject, low-complexity, coverage masking, etc.) —
  those are in `RTE_detection_review/08_synthesis_and_recommendations.md`. Some
  (e.g. skipping high-coverage pileups) would *also* speed discovery up, but they
  change the output and belong with the filter work, not here.
- `combine_insertions` / `genotype` / `combine_genotypes` — Python glue, explicitly
  outside the Rust discovery port (Plan Stage 3).
- `extend_mates()` — a no-op in the Python and intentionally omitted from the Rust
  port (`discovery.rs:382–383`); nothing to speed up until its fate is decided.
