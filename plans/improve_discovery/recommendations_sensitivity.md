# Improving retrotransposition-detection sensitivity in the Rust discovery step

Scope: this document reviews **only the Rust discovery crate**
(`rust/peartree-discovery/src/`) and proposes changes to raise the *recall*
(sensitivity) of the discovery step for genuine RTE insertion events. It does not
touch the Python code or the downstream combine/genotype stages. Rationale for
each item is grounded in the RTE_detection_review documents (referenced as
`review §NN`) and cross-checked against the actual Rust source.

The Rust port is a faithful, near-byte-identical reimplementation of
`src/discovery.py`, so every sensitivity limitation the review identified in the
Python is present here too. Line references below are to the Rust files.

---

## 0. How a read reaches an emitted insertion (the funnel)

Every true event must survive this sequential gauntlet. Each gate is a place a
real insertion can be silently lost. File:line refer to the Rust source.

| # | Gate | Location | Drops |
|---|------|----------|-------|
| 1 | `mapq < MIN_MAPQ` (40) → poly-A-only path | `discovery.rs:185` | junction reads anchored near/in repeats |
| 2 | secondary / qcfail / duplicate | `discovery.rs:191` | (correct) |
| 3 | `ref_name.len() > 5` | `discovery.rs:202` | **any contig whose name is >5 chars** (see R9) |
| 4 | must be single-sided soft-clip (or longer side wins; equal → dropped) | `discovery.rs:217-233` | both-ends-equal clips |
| 5 | `maps_fully_elsewhere` (XA/SA full-length) | `discovery.rs:237` | multi-mapping anchors |
| 6 | clustering into `<6 bp` groups | `discovery.rs:152` | — (grouping only) |
| 7 | per-read `clipped.len() >= MIN_GOOD_BASES` (10) to count | `model.rs:140` | short-clip evidence reads |
| 8 | **`n` reads at the *exact* modal breakpoint `>= 2`** | `model.rs:148-150` | jittered-junction events (see R2) |
| 9 | single-read rescue only if exact `AAAAAAAA`/`TTTTTTTT` 8-mer | `model.rs:85,95,154-155` | single-read events without a clean homopolymer terminus |
| 10 | `clipped_cons.len() <= MIN_CLIP_LEN` (12) | `model.rs:221` | error-truncated consensuses |
| 11 | `unclipped_cons.len() <= 40` | `model.rs:224` | error-truncated anchors |
| 12 | first-24-bp `n`-polymer reject | `model.rs:228-235` | MEIs anchored in microsatellite |
| 13 | output pairing: needs L+R (2–40 bp TSD) **or** clip+poly-A | `discovery.rs:441-469` | **all one-sided clipped events** (see R5) |

Two facts dominate everything below:

- **PEAR-TREE discovery is a pure split-read caller.** There is no
  discordant-read-pair leg at all (the Python's dead-code branch is simply absent
  in Rust, `discovery.rs:379-385`). Any event where no single read crosses the
  junction with a clip is invisible (review §07-finding5, §03).
- **The thresholds are compile-time constants** (`config.rs`), so none of this is
  tunable without a rebuild, and per-species tuning (the shipped `config_hs.py` /
  `config_mm.py` use `min_mapq: 60`, not the `40` hardcoded at `config.rs:4`) is
  impossible. See R0.

---

## Priority 1 — highest sensitivity return

### R0. Make discovery thresholds runtime parameters (unblocks everything else)

**Problem.** `config.rs` hardcodes `MIN_MAPQ=40`, `MIN_CLIP_LEN=12`,
`MIN_EVIDENCE_READS_PER_BREAKPOINT=2`, `MIN_GOOD_BASES=10`, `POLYA_CUTOFF=12`,
etc. as `const`. Every sensitivity experiment below requires a recompile, and the
Rust `MIN_MAPQ=40` does **not** match the `min_mapq: 60` in the production human
and mouse configs — so the port is only behaviourally identical to the *generic*
`config.py`, and any real-WGS validation must pin this explicitly.

**Change.** Load the `discovery` block from a config file / CLI flags (the code
already anticipates this: `config.rs:2` "Stage 2 will make these load from a
file"). Thread a `DiscoveryConfig` struct through `Discovery`, `join`, and
`PolyABreakpoint`. This is a prerequisite for tuning recall on real data rather
than guessing.

**Trade-off.** None; pure enabler. **Validation.** Set the values to the current
constants and confirm byte-identical output on the test BAM
(`13:32992169-32992177` IAPEz call must persist, review §08.4).

---

### R2. Count evidence within the cluster window, not at the exact modal position

**This is the single most impactful, lowest-risk recall fix in the Rust code.**

**Problem.** Breakpoints are *clustered* with a tolerance of `<6 bp`
(`discovery.rs:152`), but the 2-read evidence requirement counts reads only at the
**exact** modal coordinate:

```rust
let (best_bp, n) = most_common_first(&bps);   // model.rs:148  — exact-position mode
if n < MIN_EVIDENCE_READS_PER_BREAKPOINT { ... reject ... }  // model.rs:150
```

`most_common_first` (`model.rs:60-79`) counts identical `i64` positions. Genuine
TPRT junctions routinely place the soft-clip boundary 1–3 bp apart across reads
because of target-site microhomology, TSD jitter, and aligner-dependent clip
placement. Two reads that both support the same real junction but land at, say,
`pos` and `pos+2` are in the *same* `<6 bp` group yet each has count 1 → `n = 1`
→ the event is rejected (unless it happens to carry an exact poly-A terminus).
The clustering tolerance (6 bp) and the counting tolerance (0 bp) are mismatched.

**Change.** Count all precise, QC-passing reads in the group that fall within a
tolerance (e.g. `±max_bp_jitter`, default 3–5 bp) of the modal position, and use
that as `n`. Keep the modal coordinate as `best_bp` for consensus alignment (the
`delta_bp` machinery at `model.rs:174-216` already handles per-read offsets from
`best_bp`, so it is designed for exactly this). Concretely, replace the exact
`most_common_first` count with a windowed count:

```rust
let n = bps.iter().filter(|&&p| (p - best_bp).abs() <= max_bp_jitter).count();
```

**Impact.** Recovers low-VAF / low-coverage events whose two supporting reads
disagree by 1–3 bp — a large fraction of real 2-read events. Directly addresses
review R3 (the "6 bp vs 40 bp window" inconsistency) but is more precise: the bug
is 6 bp-cluster vs 0 bp-count.

**Trade-off.** Slightly relaxes the evidence definition; pairs naturally with R1
(a positive-hallmark or self-consistency check) to hold specificity. **Validation.**
Spike-in sweep across VAF 5/10/25/50 % at 15/30/60× (review §08.4,
MEIsimulator); expect the largest recall gain at low VAF/coverage.

---

### R3. Lower the MAPQ anchor floor, with compensating guards

**Problem.** `discovery.rs:185`: any read with `mapq < MIN_MAPQ` (40 in Rust, 60
in prod) is **never** usable as a junction anchor — it is only inspected for a
poly-A tail (`PolyABreakpoint::find_polya`) and otherwise discarded. MEIs insert
into and near repeats, exactly where anchoring reads carry MAPQ 20–40. The review
(§12.4) notes xTEA anchors clips at **MAPQ ≥ 12**, Delly ≥ 20, MEIGA ~15, and
states plainly that *"a hard MAPQ≥40 floor is stricter than necessary if the other
guards are present."*

**Change.** Drop the anchor floor to ~20–30 (tunable via R0). To keep
specificity, add xTEA-style compensators that the review recommends:
- a per-locus low-MAPQ-clip ratio cap (xTEA `MAX_LOWQ_CLIP_RATIO=0.65`): reject a
  breakpoint if too high a fraction of its supporting clips are low-MAPQ;
- keep the existing `maps_fully_elsewhere` XA/SA guard (`discovery.rs:237`), which
  already removes the worst multimappers.

**Trade-off.** More FP surface near repeats; must ship with the ratio guard and
ideally the R1 hallmark score. **Validation.** `check_germline_coverage.py`
against known germline heterozygous insertions (review §08.4) — these are the
events most often lost at a high MAPQ floor.

---

### R5. Emit one-sided clipped breakpoints as candidates

**Problem.** The output pairing loop only ever emits an event when it can pair a
left and a right breakpoint within a 2–40 bp TSD, or pair one clipped side with a
nearby poly-A (`discovery.rs:441-469`). Two structural losses follow:

1. **One-sided events are never emitted.** A 5′-truncated L1, or any insertion
   where only one junction produced clean clipped reads (the other end fell in
   unmappable element sequence or below threshold), yields a single L or R
   breakpoint with no partner and no poly-A → dropped entirely. This is the
   common case for the 5′ end of truncated L1s and for Alu/SVA (review §04, §08
   R9).
2. **The loop discards the tail of the longer list.** `while il < l.len() && ir <
   r.len()` (`discovery.rs:441`) terminates as soon as *either* side is exhausted;
   and in the matched branch it advances only `il` (`discovery.rs:467`). Remaining
   unpaired breakpoints in the longer list are never considered.

**Change.** After the L↔R / clip↔poly-A pairing pass, emit the *unpaired*
surviving `Breakpoint`s as single-anchor candidates (with a flag in the FASTQ
name, e.g. `:LEFT_ONLY:`), so downstream combine/annotate can decide. Gate them
on stronger per-side evidence (e.g. `n_reads >= 2` after R2, or a poly-A/TSD/EN
hallmark from R1) to avoid flooding output. At minimum, fix the loop so the tail
of the longer list is visited rather than silently truncated.

**Trade-off.** More candidates → heavier downstream filtering; the review's R9
explicitly anticipates a "heavy filter" for poly-A-only/one-sided loci rather
than dropping them. **Validation.** Spike-in sweep over 5′-truncation length
(review §08.4) — recall for truncated L1 should rise most.

---

## Priority 2 — hallmark capture and consensus robustness

### R1. Score positive TPRT hallmarks instead of hard single-read poly-A gating

**Problem.** The only way a sub-threshold (single-read) event survives is an
**exact** homopolymer terminus: `bp.clipped[-8:] == "AAAAAAAA"` (CLIP_LEFT) or
`bp.clipped[:8] == "TTTTTTTT"` (CLIP_RIGHT), at `model.rs:85,95,154-155`. A single
sequencing error, a 7-A run, or an interrupted poly-A (`AAAAAAGAAAAA`) fails the
literal match and the event is lost. There is no credit for a TSD or an
endonuclease (EN) motif.

**Change.** Replace the exact 8-mer test with a **fractional-purity poly-A/T
test** (e.g. ≥ 8 of the terminal 10 bases are A, MEIGA uses window 8 / purity
95 %; review §12.4), and add a small additive positive score for:
poly-A length × purity; canonical TSD length 2–20 bp; and the L1 EN nick motif
(`TT|AAAA`, 5′-`TTAAAA`) at the anchor (review §08 R8; the review notes PEAR-TREE
scoring the EN motif would *lead the field*). Let a strong hallmark score rescue a
1-read event and let a hallmark-free clip require the full R2 count.

**Critical caveat (recall-protective).** **Never require poly-A globally.** Mouse
ERV/LTR elements (IAP, MusD/ETn, MMERVK) have *no* poly-A tail and must be scored
by LTR-consensus identity + short (~6 bp) TSD instead (review §08 R8, §11.4). A
global poly-A filter would delete every genuine ERV insertion — including the
IAPEz correctness-test call.

**Trade-off.** Raises both recall and precision when done additively; the risk is
element-class bias, mitigated by the ERV caveat. **Validation.** Real-data ERV
recall (cannot use MEIsimulator, which emits no EN motif / no ERV — review §08.4).

---

### R6. Loosen poly-A *capture* constraints in `find_polya`

**Problem.** `PolyABreakpoint::find_polya` (`polya.rs:131-164`) drops a poly-A read
unless **all** of:
- `!is_proper_pair` (`polya.rs:132`) — a poly-A read flagged as a proper pair is
  never captured, yet real poly-A-containing junction reads are frequently in
  nominally proper pairs;
- `mate_is_mapped` (`polya.rs:135`);
- a **contiguous** run of ≥ `POLYA_CUTOFF` (12) A or T (`polya.rs:138-140`, exact
  `[b'A';12]` substring) — an interrupted or 11-bp tail fails;
- for poly-A, the run must start at `fi > 6` (`polya.rs:78`), so a tail beginning
  within the first 6 bp yields `clipped = None` → dropped at `polya.rs:160`.

**Change.** (a) Allow poly-A capture regardless of `is_proper_pair` (a true
insertion can leave both mates mapped in-window). (b) Replace the exact contiguous
12-mer with a purity-based run detector (`find_parts` already computes the longest
run, `polya.rs:24-65` — switch the accept test from "contiguous ≥12" to "≥N bases,
≥95 % purity in a ≥8 window"). (c) Make the `fi > 6` and `POLYA_CUTOFF` values
config-driven (R0) and consider lowering the minimum run toward xTEA's 7–9 bp for
*pure* poly-A clips (review §12.4, `MINIMUM_POLYA_CLIP=7`).

**Trade-off.** More poly-A candidates; guard downstream with R5's flagging and
R1's scoring. **Validation.** Alu/SVA and 5′-truncated-L1 spike-ins (poly-A end is
their strongest signal).

---

### R7. Make consensus building tolerant of a single disagreement

**Problem.** `find_consensus` (`filters.rs:83-124`) stops at the **first ambiguous
position**:

```rust
let delta_best = sorted[0].1 - sorted[1].1 - sorted[2].1 - sorted[3].1;  // filters.rs:115
if delta_best > 0 { ... } else { break; }                                // filters.rs:119-121
```

The winning base's quality sum must exceed the summed quality of *all three* other
bases. With only 2 evidence reads, a single disagreement (two reads, different
base, comparable Q) gives `delta_best = 0` → the consensus is **truncated at that
position**. A truncated consensus then trips `clipped_cons.len() <= 12`
(`model.rs:221`) or `unclipped_cons.len() <= 40` (`model.rs:224`) and the whole
breakpoint is rejected. So one sequencing error in one of two reads can delete a
real event.

**Change.** Use a majority/quality vote that (a) emits the best base with an
ambiguity/quality annotation and *continues* past isolated disagreements rather
than halting, or (b) at minimum, only halts after K consecutive ambiguous
positions, not the first. This keeps consensus length above the length gates in
the presence of normal error. (The `delta_best > 0` rule is unusually strict; a
simple "best base strictly beats second-best" would already help.)

**Trade-off.** Slightly noisier consensus tails; downstream bowtie2/DFAM
annotation tolerates a few mismatches, so net recall gain. **Validation.**
Differential run: count how many breakpoints move from `clipped_failed` /
`unclipped_failed` into `passed` (the `Breakpoint.stats` counters, Python side)
on a real BAM.

---

## Priority 3 — smaller / situational

### R9. Replace the `ref_name.len() > 5` contig filter (latent recall cliff)

`discovery.rs:202` and `discovery.rs:419,430` skip any contig whose **name** is
longer than 5 characters. On `chr1`…`chr22`/`chrX` this is fine, but on RefSeq
(`NC_000014.9`), Ensembl-style, or T2T/`hs1` scaffold naming it silently skips
**main chromosomes**, zeroing recall genome-wide with no error (review
§07-finding6). Replace with an explicit include-list / a
"primary assembly" check driven by the BAM header (drop `_alt`/`_random`/`chrUn`
by pattern, keep everything else), configurable via R0.

### R8. Allow short clips when the clip is pure poly-A

`MIN_CLIP_LEN=12` (`config.rs:5`) gates both the clipped consensus (`model.rs:221`)
and the poly-A rescue (`model.rs:87,97`). Review §10.1 notes ~3.86 % of clips sit
just under 15 bp — a real population at the boundary. Permit a shorter clip (down
to ~7 bp) *only* when it is a pure poly-A/T terminus (xTEA behaviour), keeping the
12 bp floor for arbitrary sequence.

### R4. Emit the `Breakpoint.stats` rejection counters from Rust

The Python tracks per-side `too_few / rescued_pA / excluded / too_few_after_filter
/ clipped_failed / unclipped_failed / polymer / passed`. The Rust `join`
(`model.rs`) drops reads at the same points but **keeps no counters** — recall
losses are invisible per run. Add an atomic counter struct incremented at each
`return None` in `join` and print it at end of `cleanup`. This is diagnostic, not
a recall change, but it is how you will *measure* every recommendation above
(review §07 R4).

### R10. (Architectural, largest ceiling) Add an SA-bridge / discordant leg

The only fix for junctionless, large-TSD, low-coverage, and unmappable-junction
events (review §08 R11, §03). Cheapest first step that fits the clip-first design:
adopt MEIGA's `SA_as_DISCORDANTS` trick — turn a read's SA-tag supplementary clip
into a *pseudo-discordant* so a single split read substitutes for a pair. The
`BamRead` already parses the `SA` tag (`read.rs:100-103`) and the `Breakpoint`
already carries a `bp_precise: bool` flag (`model.rs:19`) intended for imprecise
(discordant) breakpoints — the scaffolding is present. Require the pseudo-discordant
to converge with clipped/poly-A evidence before emitting (keeps specificity;
review insists the discordant leg be *complementary*, §08 R11), and add xTEA's
clip↔discordant geometric-consistency check as the confirmation gate.

---

## Suggested order of work

1. **R0** (config plumbing) — unblocks tuning; zero behaviour change.
2. **R2** (windowed evidence count) + **R9** (contig filter) — highest recall per
   line changed, low risk.
3. **R5** (emit one-sided) + **R7** (consensus tolerance) — recover truncated /
   error-affected events.
4. **R3** (MAPQ floor) + **R6/R8** (poly-A capture) with **R1** (hallmark score)
   as the specificity counterweight.
5. **R4** (stats) alongside all of the above to measure them.
6. **R10** (SA-bridge / discordant) — the biggest ceiling, largest effort; do last.

## Validation harness (applies to every item)

- **Correctness gate:** the test BAM must keep producing the IAPEz call at
  `13:32992169-32992177` (review §08.4). Any behaviour change necessarily breaks
  byte-identity with the Python port — that is expected once R0 is in; pin the old
  constants to reproduce the baseline.
- **Recall:** MEIsimulator spike-ins sweeping VAF (5/10/25/50/100 %), coverage
  (15/30/60×), 5′-truncation length, and element class.
- **Simulator blind spots** (validate on real data, not MEIsimulator): EN motif,
  5′ twin-priming inversions, and all mouse ERV/LTR insertions — the simulator
  emits none of these (review §08.4).
- Use the **R4 stats counters** to attribute every recovered event to the gate
  that previously dropped it.
