# Implementing the RTE-detection review in the Rust discovery step

*Scope: this note reads the Rust discovery port (`rust/peartree-discovery/src/*.rs`)
against the 18 recommendations (R1–R18) in
[`RTE_detection_review/08_synthesis_and_recommendations.md`](RTE_detection_review/08_synthesis_and_recommendations.md)
and the code review in [`07_peartree_code_review.md`](RTE_detection_review/07_peartree_code_review.md).
For each recommendation it answers: **does it belong in the discovery step, how would
we implement it in the Rust port concretely, and what is the effort/risk?***

*Speed is out of scope here — that is [`recommendations_speed.md`](recommendations_speed.md).
This is purely about the artefact-rejection / sensitivity changes, mapped onto the Rust
code as it stands on branch `PEAR-TREE2`.*

---

## 0. TL;DR

The review's recommendations split three ways relative to the Rust discovery step:

| Belongs in Rust discovery **now** | Belongs in discovery **later** (needs infra) | **Not** discovery (combine / genotype / separate tool) |
|---|---|---|
| **R15** SMS both-ends-clip reject *(one filter, highest ROI)* | R1 / R16 local-depth mask + coverage-adaptive thresholds | R9 poly-A-only fate (`combine_insertions`) |
| **R2** wire in low-complexity filters | R5 / R7 include-BED + blacklist BED | R12 Delly companion run |
| **R4** emit `Breakpoint.stats` sidecar | R8 positive hallmark scoring (poly-A/TSD/EN) | R13 MEI-polymorphism DB subtraction |
| **R3** reconcile the 6-bp vs 40-bp windows | R14 within-BAM clip-recurrence flag | R18 final ML classifier (combine) |
| **R10** decide `extend_mates` fate (already omitted in Rust) | R17 RepeatMasker young-copy self-mask | R11 (mostly) discordant-pair engine — but see cheap first step |

**Do first, in order: R15 → R2 → R4 → R3.** These are all local to the discovery
code, need no reference genome, no new file formats, and no architectural change. R15
is a ~10-line filter that three independent mature callers converge on; R2 and R4 are
"wire in / surface code that is already conceptually present." Everything else needs
either a config/BED plumbing layer, a coverage estimator, or crosses into `combine`.

A key structural fact: **the Rust port currently has no stats counters and no
low-complexity filter** (only the n-polymer check in `model.rs:227-235`). So R2 and R4
are *ports of Python logic that was never carried over*, not new algorithms.

---

## 0.1 Cross-checked against the corrected review (2026-07-15 edits)

The review docs were revised after an accuracy pass (`git diff RTE_detection_review/`).
Five assertions were corrected; I re-verified each against these recommendations:

- **Poly-A detection window is "within 30 bp of *either* end of the insertion," not the
  3′ end only** (§02, §09, MEIGA spec: ≥15 bp, ≥90% purity). **This one matters here** —
  it is folded into R8 (the concrete poly-A scoring target) and R2 (the escape hatch)
  below, and it *validates* PEAR-TREE's existing design: `polya.rs` already scans A-runs
  **and** T-runs on both `CLIP_LEFT`/`CLIP_RIGHT` sides, i.e. both ends, so no orientation
  assumption in the port needs changing.
- **MEIGA detection precision is >95% (F1 ≤ 99.55), not ">99%"** — the ">99%" was
  *source-inference* specificity (sensitivity only 47.6%) (§03). Not cited here; R18's
  claim is only that MEIGA *ends in a logistic-regression classifier*, which is unchanged.
- **Nam 2023's "≥3 discordant pairs in blood" is a minimum-support rule, not the germline
  test** (§03, §09). Not cited here (germline discrimination is a combine/genotype concern,
  R13, and PEAR-TREE uses cohort partition instead of a matched normal).
- **Chen's ">1,000 palindromes/sample" was a conflation** (real figure: median 115 artefact
  SNVs/sample, enzymatic prep) (§05). Not cited here; R6/R7/R15's palindrome rationale rests
  on the §10–§11 *PEAR-TREE-measured* prevalences, which are unchanged.
- **`discovery.py` cluster-window lines are 178,189 (not 136,147)** (§07 finding 2). My R3
  anchors are the **Rust** line numbers (`discovery.rs:152` cluster, `:443,454` TSD), so
  they are unaffected; noted so the Python cross-reference stays correct.

**Net effect: no discovery-step recommendation changes.** Only R8 and R2 gain a more precise
poly-A spec (below).

---

## 1. Tier-1 — implement now (local to discovery, no new infra)

### R15 — Reject both-ends-soft-clipped (SMS) reads  ★ do this first

**Why:** the read-level fingerprint of a reference palindrome / inverted repeat,
adapter dimer, or spurious multimap. MEIGA, v-TraFiC and xTEA all reject it; it keys on
CIGAR *shape*, not "is a repeat," so it is safe for genuine repeat-consensus clips
(incl. mouse ERV). Largely subsumes R6/R7 for the enzymatic-prep palindrome problem.

**Where it goes wrong today:** [`discovery.rs:217-231`](rust/peartree-discovery/src/discovery.rs)
— when `left_soft && right_soft`, the current code keeps the read and classifies it on
its *longer* clip. That is exactly the read SMS wants to drop.

**Implementation:** add a threshold to `config.rs` and one guard in the clip-classification
block.

```rust
// config.rs
pub const MAX_CLIP_CLIP_LEN: usize = 8; // xTEA MAX_CLIP_CLIP_LEN; both-ends clip reject
```

```rust
// discovery.rs, inside the `right_soft && left_soft` arm (currently lines 221-228)
if right_soft && left_soft {
    // SMS reject: a read genuinely soft-clipped on BOTH ends is a chimera/palindrome,
    // not a single junction. Tolerate a few stray bases (adapter/noise) on the short end.
    if left_len.min(right_len) >= MAX_CLIP_CLIP_LEN {
        continue; // R15
    }
    // otherwise keep the old "longer clip wins" behaviour
    if left_len > right_len { Some(CLIP_LEFT) }
    else if right_len > left_len { Some(CLIP_RIGHT) }
    else { None }
}
```

**Note:** keep the tolerance (`>= MAX_CLIP_CLIP_LEN` on the *shorter* clip) rather than
"both ends soft at all" — a real junction read often has a couple of clipped bases on
the anchored end. Gate on the short clip so only genuine dual-clips are dropped.

**Effort:** ~10 lines. **Risk:** low, but it *changes the emitted breakpoint set*, so it
must ship behind a toggle (like `PEARTREE_KEEP_FULLMAP`) and be measured against the test
BAM call at `13:32992169-32992177` and one real WGS BAM. Add `PEARTREE_KEEP_SMS=1` to
disable, mirroring the existing full-map toggle in `config.rs:16`.

### R2 — Port and wire in the low-complexity filters

**Why:** `is_low_complexity` / `has_well_defined_breakpoint` exist in
[`src/sequence_checks.py:22,43`](src/sequence_checks.py) but were **never ported to Rust**
and are not on the Python discovery path either. The empirical scan (§10) found ~12% of
clips are low-complexity reaching output.

**Implementation:** add to `filters.rs`, call from `model.rs::join` on the consensus,
next to the existing n-polymer check ([`model.rs:227-235`](rust/peartree-discovery/src/model.rs)):

```rust
// filters.rs — port of sequence_checks.is_low_complexity (≤2 distinct ATGC bases → low)
pub fn is_low_complexity(seq: &[u8]) -> bool {
    let (mut a, mut t, mut g, mut c) = (false, false, false, false);
    for &b in seq {
        match b.to_ascii_uppercase() { b'A'=>a=true, b'T'=>t=true, b'G'=>g=true, b'C'=>c=true, _=>{} }
    }
    (a as u8 + t as u8 + g as u8 + c as u8) <= 2
}
```

Call it in `join` after `unclipped_cons` is built and **before** returning `Some(b)`,
rejecting when *both* the clip tail and the unclipped head are low-complexity — but keep
a **poly-A escape hatch** (MEIGA/§10.2): do not reject if the clip is a clean poly-A/T
tail (a run of ≥15 bp at ≥90% purity within 30 bp of *either* end, per the corrected
MEIGA spec, §02/§09), or the `rescue_polya` path already handled it. `polya.rs` already
tests both ends, so the escape hatch reuses that logic rather than adding a 3′-only check.

> **Optional upgrade (medium):** the code-grounded form in §12.2–12.3 is v-TraFiC's
> `dustmasker` or MEIGA-LR's k-mer entropy (`Σ_K|freq − 0.25^K| > 0.83`, K∈{1..4}). The
> entropy version is a self-contained ~30-line function with no external dependency and
> is a strict improvement over the ≤2-distinct-bases heuristic. Do the cheap version
> first; upgrade to entropy only if the stats sidecar (R4) shows low-complexity clips
> still leaking through.

**Effort:** cheap version ~20 lines. **Risk:** low; still gate behind a toggle + measure.

### R4 — Emit a `Breakpoint.stats` sidecar JSON  ★ do this alongside R15/R2

**Why:** the Rust port tracks **no reject-reason counters at all**. Every filter above
(R15, R2, and the existing exclude/polymer/too-few paths in `model.rs`) should increment
a counter so that a specificity regression is visible run-to-run. This operationalises
the plan's "every change is measured" and is the cheapest way to *validate* R15/R2.

**Implementation:** a plain `struct DiscoveryStats { too_few, excluded, polymer,
low_complexity, sms_rejected, fullmap_rejected, adapter, rescued_polya, passed, ... }`
of `u64` counters, threaded through `Discovery` (it is single-threaded, so plain fields
suffice). Increment at each `continue`/`return None` site in `discovery.rs` and
`model.rs`. After `output()`, write `<out>.stats.json`. No external crate needed — hand-write
the JSON or add `serde_json`.

**Effort:** ~1–2 h. **Risk:** none (observability only; does not touch the breakpoint set).
**Do this first of all** — it is the measurement instrument for every other change.

### R3 — Reconcile the 6-bp clustering window with the 40-bp TSD window

**Why:** `cleanup` clusters clips when positions are `< 6` bp apart
([`discovery.rs:152`](rust/peartree-discovery/src/discovery.rs)) but the TSD pairing in
`output()` uses `tsd < 2` / `tsd > 40` ([`discovery.rs:443,454`](rust/peartree-discovery/src/discovery.rs)).
A 6–40 bp TSD can be clustered inconsistently.

**Implementation:** lift the hard-coded `6`, `2`, `40`, `12`, `120` into named `config.rs`
constants (`CLUSTER_WINDOW`, `TSD_MIN`, `MAX_BP_WINDOW`, `POLYA_PAIR_MIN`, `POLYA_PAIR_MAX`)
and document the relationship. This is a **refactor first, behaviour-change second**: pull
the constants out with *identical values* (byte-identical output, safe), then experiment
with reconciling them under the R4 stats + toggle.

**Effort:** low. **Risk:** the constant-extraction step is zero-risk; any value change is
measured via R4.

### R10 — `extend_mates`: already omitted in Rust; make the omission a decision

The Rust port deliberately omits `extend_mates` ([`discovery.rs:382-383`](rust/peartree-discovery/src/discovery.rs))
because the Python version is a silent no-op (iterates the already-emptied
`temporary_breakpoints`). **The Rust code is already in the "removed" state.** The only
action is to *decide*: either (a) implement real mate-extension against the final
breakpoint lists (longer clipped contigs → better DFAM/RepeatMasker annotation), or (b)
formally delete it and update the plan. Per R10 and Plan §1.3 this **needs a real-WGS
benchmark** before choosing — longer contigs can introduce new artefacts. No code change
until that benchmark exists; flagged here so it is not forgotten.

---

## 2. Tier-2 — implement after a small infra layer

These are genuinely valuable at discovery time but each needs one new capability the
port does not yet have (a coverage estimator, BED plumbing, or a bundled track loader).

### R1 + R16 — Local-depth mask and coverage-adaptive thresholds

**Why:** coverage spikes are one of the richest artefact sources; a fixed
`MIN_EVIDENCE_READS_PER_BREAKPOINT = 2` ([`config.rs:6`](rust/peartree-discovery/src/config.rs))
is simultaneously too lax in pile-ups and too strict at low coverage.

**Fit with the streaming design:** the BAM is coordinate-sorted and the port streams it
once. Maintain an **O(1)-per-read running depth**: increment a counter at each read's
`reference_start`, decrement (via a small min-heap / ring buffer keyed on `reference_end`)
as the scan passes each position. That yields exact local depth at every breakpoint
position with no per-chromosome coverage array. Histogram the depth to get the **sample
median** in the same pass.

Then:
- **R1 mask:** drop a candidate breakpoint whose local depth `> 3 × sample_median`
  (xTEA `MAX_COV_TIMES=3`). Wire the flag into `add_breakpoint` / `join`.
- **R16 adaptive:** replace the constant `MIN_EVIDENCE_READS_PER_BREAKPOINT` with a
  lookup indexed by local coverage (xTEA table: cov 5→1, 30→3, 100→8, 300→25). Interpolate.

The mature "local-pileup mask" form (§12.6) also uses **fraction of low-MAPQ reads** and
**fraction of SMS reads** in a `bkp±100` window (v-TraFiC `<30% MAPQ<10 & <15% SMS`;
MEIGA `percMAPQ`/`percSMS`). SMS fraction is free once R15 is in (we already classify
SMS reads). This is the more powerful version and needs the same windowed buffer.

**Effort:** medium (the depth estimator is the real work; ~half a day). **Risk:** medium,
measured via R4. **Sequencing note:** land R15 first — SMS fraction reuses it.

### R5 + R7 — Include-BED for contigs, and blacklist BEDs

**Why:** `ref_name.len() > 5` ([`discovery.rs:202,419`](rust/peartree-discovery/src/discovery.rs))
is a fragile proxy that breaks on `NC_000014.9`-style naming and cannot express a
blacklist. R7 wants telomere/centromere/segdup + inverted-repeat/palindrome BEDs.

**Implementation:** add a config-supplied **main-chromosome allowlist** (replacing the
`len>5` and `MT`/`chrM` checks) and an optional **exclude-interval index** (a sorted
`Vec<(start,end)>` per contig, binary-searched — no crate needed). Check breakpoint
position against the exclude set in `add_breakpoint`. R7's BEDs then just populate the
same exclude index. Build the palindrome BED **per species** (T2T human *and* GRCm39 —
mouse has ~10× the palindrome load, §11) via the ArtifactsFinder recipe in §12.8
(params corrected in §12.5: `D_LEN 8 / STEM_LEN 5 / S_LEN 2 / ±50 bp`, re-enable
`LEN≥17`).

**Effort:** BED plumbing low; *generating* the palindrome BED is a separate offline
pipeline (medium, one-off). **Risk:** low for the plumbing. **Note:** the port needs an
**indexed** reader to *skip* excluded contigs at iteration time rather than decode-then-drop
— that overlaps with the speed work (S5), so coordinate the two.

### R8 — Positive hallmark scoring (poly-A length × purity, TSD length, EN motif)

**Why:** today the port uses TSD *geometry* (`R − L`, `output()`) and poly-A *presence*
but rewards neither poly-A length×purity, canonical TSD length (2–20 bp), nor the EN
motif (`TT/AAAA`, 190× enriched, Nam 2023). An additive positive score lets us *keep*
well-supported events at low coverage while *distrusting* hallmark-free clips.

**Concrete scoring target (corrected MEIGA spec, §02/§09):** reward a poly-A/T run of
**≥15 bp at ≥90% purity found within 30 bp of *either* end of the insertion** (not 3′
only — the corrected review makes this explicit, and it matches PEAR-TREE's existing
both-ends scan in `polya.rs`); reward a TSD in the **~2–20 bp** canonical range (commonly
10–20); reward the EN nick motif `TT/AAAA` at the 5′ junction.

**Implementation:** compute the three features where the data already is — poly-A length
in `polya.rs::find_parts` (already returns `polya_len`), TSD length in `output()` (already
`= r.breakpoint − l.breakpoint`), EN motif by inspecting the reference-adjacent bases at
the nick. Emit an additive score.

**The real constraint is the on-disk contract:** discovery emits the fixed FASTQ-like
`.txt.gz`. Options: (a) encode the score in the locus name suffix (backward-compatible if
`combine` ignores unknown suffixes — verify), or (b) write it to the R4 stats sidecar keyed
by locus. Prefer (b) until `combine` is taught to read it.

**Critical caveat (§11.4):** make it **element-class-aware and never require poly-A
globally.** Mouse ERV/LTR (IAP, MusD/ETn, MMERVK) make **no poly-A** and are scored by
LTR-consensus identity + a short ~6 bp TSD. A poly-A *filter* would delete every genuine
ERV insertion. So R8 must be a *score*, never a gate, and poly-A weight must be gated on
element class.

**Effort:** medium. **Risk:** low if it stays a non-gating score written to the sidecar.

### R14 — Flag within-BAM clip-sequence recurrence

**Why:** §10.3 — ~45% of clips are non-unique across the cohort (poly-A/T + Alu-consensus
fragments). A clip identical at many independent loci is not a novel locus-specific
insertion. Discovery can catch the *within-BAM* share cheaply.

**Implementation:** a post-pass over the final breakpoint lists — `HashMap<Vec<u8>, u32>`
counting clipped-consensus sequences, then flag/drop those recurring `≥ K` times. Runs
after `cleanup`, before `output`.

**Two hard caveats that make this a Tier-2, not Tier-1, item:**
1. Filter on **sequence recurrence, not "is a repeat"** — a *locus-unique* Alu clip with a
   poly-A tail and TSD is a candidate somatic Alu insertion and must survive. Doubly true
   for mouse ERVs whose true-insertion clips *are* the IAP/LTR consensus.
2. **Cross-individual recurrence is a combine-step signal, not discovery** — a clip shared
   across many colonies of *one* individual is germline (keep); shared across many
   *individuals* is artefact (drop). Discovery is per-BAM and cannot see this; only the
   within-BAM count belongs here. The per-individual/folder recurrence stays in `combine`.

**Effort:** within-BAM version low. **Risk:** medium — must ship as a *flag in the stats
sidecar first* (observe the distribution) before it is allowed to *drop* anything.

### R17 — RepeatMasker young-copy self-mask

**Why:** a distinct artefact class — reads mis-donated by a *young* reference element of
the same family near the breakpoint, which a static IR/germline blacklist misses. TraFiC
drops calls in a same-family RM element with divergence ≤20%; xTEA `<15%`.

**Fit:** the tracks are **already in the repo** — `hs1.repeatMasker.out.gz` and
`hs1.rte.out.gz`. Load into the same interval index as R5/R7 but carry `(family,
divergence)`; drop a candidate whose breakpoint falls in a same-family element with
divergence below cutoff. **Gate on divergence, not family membership**, so a locus-unique
young insertion survives.

**Effort:** medium (parse RM `.out`, build the divergence-annotated index). **Risk:** low
if divergence-gated + measured.

---

## 3. Tier-3 — architectural / mostly not the discovery step

### R11 — Discordant-read-pair discovery (biggest sensitivity gap)

Full discordant clustering is a large addition and the plan already scopes the Rust work
to the clip-first design. **But there is a cheap first step that fits discovery:** MEIGA's
**`SA_as_DISCORDANTS`** — convert an SA-tag supplementary clip into a pseudo-discordant so
a *single split read* substitutes for a pair. The port already parses SA tags
([`read.rs:99-107`](rust/peartree-discovery/src/read.rs), used in `maps_fully_elsewhere`
and the exclude logic), so the plumbing exists. Then borrow xTEA's clip↔disc geometric
consistency (`_is_distance_consistency`: a left-clip confirmed only if right-side support
maps within one insert-size) — xTEA's single strongest FP filter. A full discordant engine
(Delly per-library 3-SD insert model) is a separate project. **Effort:** SA-pseudo-discordant
medium; full engine high. **Priority:** after Tier-1/2 land and are validated.

### Not the discovery step — note and route elsewhere

- **R9** (fate of poly-A-only loci) — `combine_insertions_intersect_insertions.py:84`,
  the unconditional `continue`. A combine-step decision; discovery already produces the
  poly-A rescue evidence it needs.
- **R12** (Delly/MELT/xTea companion run) — a separate parallel caller, orchestration not
  discovery.
- **R13** (1000G MEI / dbRIP / euL1db subtraction) — classification, combine/genotype.
- **R18** (cohort-adapted logistic/RF classifier) — the *final* combine classifier. If
  built, reuse MEIGA/xTEA feature schemas but **replace matched-normal features with the
  cohort/tree-derived R14 signals**, and gate poly-A features on element class (§11) or
  every true mouse ERV insertion is down-weighted. Discovery's role is only to *emit the
  features* (R8/R14 sidecar) the classifier will consume.

---

## 4. Proposed implementation order for the Rust port

1. **R4** stats sidecar — the measurement instrument. *(first, zero-risk)*
2. **R15** SMS reject — highest ROI, ~10 lines, behind `PEARTREE_KEEP_SMS`.
3. **R2** low-complexity filter (cheap ≤2-base version) with poly-A escape.
4. **R3** extract clustering/TSD-window constants (byte-identical refactor), then tune.
5. **R1 + R16** running-depth estimator → high-coverage mask + adaptive evidence threshold.
6. **R5/R7** contig allowlist + exclude-BED plumbing (coordinate with the indexed-reader
   speed work); generate the per-species palindrome BED offline.
7. **R17** RepeatMasker divergence self-mask (reuse the R5 interval index; tracks already
   bundled).
8. **R8** hallmark scoring → sidecar (non-gating, element-aware).
9. **R14** within-BAM clip-recurrence — flag in sidecar first, drop only after review.
10. **R11** SA-as-discordant first step, then evaluate a full discordant branch.
11. **R10** benchmark-driven decision on `extend_mates`.

Every step ships behind an env toggle (as `PEARTREE_KEEP_FULLMAP` already does,
`config.rs:16`) and is measured with the R4 sidecar.

## 5. Validation (per §8.4 and the plan)

1. **Correctness gate:** the test BAM must keep producing the IAPEz call at
   `13:32992169-32992177`. Any filter that is *on by default* must not remove it.
2. **Specificity:** diff the R4 stats sidecar + final insertion count on one real WGS BAM
   before/after each filter. A filter must be shown to remove *artefacts*, not true calls —
   check against known germline hets (`check_germline_coverage.py`).
3. **Recall (R8/R11):** MEIsimulator spike-in sweep over VAF (5/10/25/50/100%), coverage
   (15/30/60×), 5′-truncation, event class — but validate EN-motif, twin-priming inversion
   and **mouse ERV** recall on *real* data, since MEIsimulator models none of them (§12.9).
4. **Byte-identical baseline:** keep Python discovery as the oracle; each toggle **off**
   must reproduce the current byte-identical `.txt.gz` (the differential harness at
   `rust/peartree-discovery/tests/differential_test.sh`), so every behaviour change is
   isolated to its flag.
