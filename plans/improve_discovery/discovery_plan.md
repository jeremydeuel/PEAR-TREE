# PEAR-TREE2 discovery — consolidated recommendation register & decision sheet

*Merges every recommendation from `implement_review.md`, `recommendations_sensitivity.md`
and `recommendations_speed.md` into one namespace, annotates each with the adversarial
findings from `adverserial_scientific.md`, gives a recommended disposition, and leaves a
**Decision** field for you to set.*

*A second, code-level adversarial pass (line-reference audit + byte-identity pressure-test
against the actual crate and the Python originals) has been folded in; its additions are
tagged **`2nd-review`** so their provenance stays distinct from the `#N` findings above.
That pass confirmed the vast majority of `file:line` citations are exact, verified the
`min_mapq 40-vs-60` and `extend_mates` no-op claims against the Python source, and confirmed
both RepeatMasker tracks are on disk — but it also hardened two byte-identity claims (SPD-4,
SPD-5e) that were previously only flagged as "unproven."*

**How to use:** for each item set `Decision:` to one of **ADOPT** / **ADOPT WITH CHANGES** /
**DON'T ADOPT** / **DEFER**, and edit the `Notes:` line. My recommended disposition is
pre-filled as a starting point — override freely.

**Namespacing (finding #1).** The two behaviour plans reused `R2`, `R8`, `R1` for
*opposite-signed* changes. To stop that hazard this register prefixes every item:
`SPD-` (speed), `OBS-` (observability/refactor), `SPEC-` (specificity / artefact rejection),
`SENS-` (sensitivity / recall), `ARCH-` (architectural). Original IDs are shown in brackets.

---

## 0. Summary table (recommended disposition at a glance)

| ID | Item | Source | Original tier | Recommended disposition |
|---|---|---|---|---|
| SPD-1 | Lazy per-read decode | speed R1 | do early | **Adopt** |
| SPD-2 | Multithreaded BGZF decode | speed R2 | do early | **Adopt** |
| SPD-3 | Contig-level parallelism (rayon + indexed reader) | speed R3 | later | **Adopt with changes** |
| SPD-4 | Mate-coordinate fetch (kill 2nd BAM pass) | speed R4 | later | **Adopt with changes** |
| SPD-5 | Cheap internals (int contig, u8 quals, fxhash, record reuse, skip-at-reader) | speed R5a–e | anytime | **Adopt** 5a/c/d; **with changes** 5b; **5e → move to SPEC-5/SENS-4** (not byte-identical) |
| OBS-1 | `Breakpoint.stats` reject-counter sidecar | impl R4 / sens R4 | do first | **Adopt** |
| OBS-2 | Extract clustering/TSD-window constants (byte-identical) | impl R3 | Tier 1 | **Adopt** |
| OBS-3 | Reconcile 6 bp cluster vs 40 bp TSD window (value change) | impl R3 | Tier 1 | **Adopt with changes** |
| OBS-4 | Runtime config params (no recompile) | sens R0 | P1 | **Adopt** |
| OBS-5 | Decide `extend_mates` fate | impl R10 | Tier 1 decision | **Defer** |
| SPEC-1 | SMS both-ends-clip reject | impl R15 | do FIRST | **Adopt with changes** (off by default) |
| SPEC-2 | Low-complexity filter wired into discovery | impl R2 | Tier 1 | **Adopt with changes** (entropy version) |
| SPEC-3 | Local-depth / MAPQ / SMS pileup mask | impl R1 | Tier 2 | **Adopt with changes** |
| SPEC-4 | Coverage-adaptive evidence threshold | impl R16 | Tier 2 | **Adopt with changes** |
| SPEC-5 | Contig allowlist + exclude-BED plumbing | impl R5/R7 | Tier 2 | **Adopt** (plumbing) |
| SPEC-6 | Palindrome/IVR reference blacklist (per species) | impl R7 | Tier 2 | **Defer** |
| SPEC-7 | RepeatMasker young-copy self-mask (divergence-gated) | impl R17 | Tier 2 | **Adopt with changes** |
| SPEC-8 | Within-BAM clip-recurrence flag | impl R14 | Tier 2 | **Adopt with changes** (flag only, never drop) |
| SENS-1 | Windowed evidence count (not exact modal position) | sens R2 | P1 | **Adopt with changes** (needs a real gate) |
| SENS-2 | Lower MAPQ anchor floor 40→20–30 | sens R3 | P1 | **Adopt with changes** (guard ships same commit) |
| SENS-3 | Emit one-sided clipped candidates | sens R5 | P1 | **Don't adopt** (as scoped) / **Defer** |
| SENS-4 | Replace `len(name)>5` contig filter | sens R9 | P1 | **Adopt** |
| SENS-5 | Positive-hallmark scoring (poly-A/TSD/EN), non-gating | sens R1 / impl R8 | P2 | **Adopt with changes** |
| SENS-6 | Loosen `find_polya` capture | sens R6 | P2 | **Adopt with changes** |
| SENS-7 | Consensus tolerant of single disagreement | sens R7 | P2 | **Adopt with changes** |
| SENS-8 | Short clips allowed when pure poly-A | sens R8 | P3 | **Adopt with changes** |
| ARCH-1 | SA-bridge → pseudo-discordant leg | impl R11 / sens R10 | Tier 3 | **Defer** |
| ARCH-2 | Poly-A-only locus fate in combine | impl R9 | not discovery | **Defer** (blocks SENS-3) |
| ARCH-3 | Delly/MELT/xTEA companion run | impl R12 | not discovery | **Defer** |
| ARCH-4 | MEI-polymorphism DB subtraction | impl R13 | not discovery | **Defer** |
| ARCH-5 | Cohort-adapted final classifier | impl R18 | not discovery | **Defer** |
| VAL-1 | Validation harness (see finding #6) | §8.4 / all | — | **Adopt with changes** (blocker) |

---

## 1. Speed (`recommendations_speed.md`) — lowest controversy

### SPD-1 — Lazy per-read decode *(speed R1)*
Reorder `BamRead::from_record` so seq/qual/name/tags are built only for surviving clipped
candidates. Byte-identical; large single-thread win.
- **Finding:** none for correctness — the refactor is byte-identical on its own.
- **Finding (`2nd-review`):** not quite "pure mechanical / ~99% skip." The low-MAPQ poly-A
  branch (`discovery.rs:185-190`) runs *before* the secondary/dup gate and still needs `seq`
  for every read passing `!is_proper_pair` + `mate_is_mapped` — a non-trivial slice in
  repeat-rich low-MAPQ regions, so the "~99%" is data-dependent. More important, it **couples
  to SENS-6:** R1's savings on that branch depend on `find_polya` bailing on `is_proper_pair`
  (`polya.rs:132`), which is *exactly* the bail SENS-6(a) proposes to remove. If SENS-6 lands,
  every low-MAPQ read decodes `seq` and this saving evaporates — sequence SPD-1 before SENS-6.
- **Recommended:** **ADOPT** (still byte-identical standalone; note the SENS-6 ordering).
- `Decision:` ADOPT
- `Notes:`

### SPD-2 — Multithreaded BGZF decode *(speed R2)*
Wire the already-parsed-but-unused `bam_threads` to a multithreaded reader.
- **Finding:** none for correctness — decode block order is preserved, so byte-identical.
- **Finding (`2nd-review`):** "Low effort" is mildly optimistic. `Cargo.toml` pins
  `noodles-bam = "0.79"` / `noodles-sam = "0.75"` and **does not list `noodles-bgzf`**; wiring
  `MultithreadedReader` means adding it at the exact version noodles-bam 0.79 pins transitively
  (a mismatched explicit version type-conflicts against `bam::io::Reader::new`). Low *risk*
  stands; it is not the one-liner the source table implies.
- **Recommended:** **ADOPT.** Makes `--threads` mean something.
- `Decision:` ADOPT
- `Notes:`

### SPD-3 — Contig-level parallelism, rayon + indexed reader *(speed R3)*
- **Finding (#10):** correct in principle, but couples to mate resolution (SPD-4) and changes
  RSS. Merge order must reproduce the existing per-contig sort exactly.
- **Recommended:** **ADOPT WITH CHANGES** — land after SPD-4; prove byte-identical on the
  real-WGS differential (which does not yet exist — see VAL-1).
- `Decision:` ADPèT WITH CHANGES (see above)
- `Notes:`

### SPD-4 — Mate-coordinate fetch, eliminate 2nd BAM pass *(speed R4)*
- **Finding (#10):** "byte-identical because the same reads are visited" is **not proven** —
  coordinate-fetch + coalescing changes iteration order, and mates are pushed in encounter
  order (`mate_seqs`, the `MATE{i}` numbering at `discovery.rs:504,517`). Order-sensitivity
  must be demonstrated, not assumed.
- **Finding (`2nd-review`) — the claim is not just unproven, it is _false_:** the `find_mates`
  gate excludes secondary/qcfail/duplicate but **not supplementary** reads
  (`discovery.rs:349`), and it matches by `qname` across the *whole* BAM. A mate with
  supplementary alignments therefore contributes its primary record at PNEXT **plus every
  supplementary record at other loci**, each appended as its own `:MATE{i}:` line. A fetch at
  the single `(mate_ref, mate_pos)` from PNEXT captures only the primary → strictly fewer MATE
  lines → **provably not byte-identical.** Same mechanism corrupts poly-A mates: `set_mate`
  overwrites `reference_name`/`breakpoint` on every call (`polya.rs:112-127`), so "last matching
  record wins" — dropping a higher-coordinate supplementary changes the emitted breakpoint.
  This is material precisely because RTE mates land in repeats where supplementary alignments
  are common. So it is "the same reads visited" that is wrong, independent of iteration order.
- **Recommended:** **ADOPT WITH CHANGES** — pursue as a *speed* change, but **drop the
  byte-identical label**; either also fetch each mate's SA-tag supplementary loci, or accept a
  documented output delta and re-baseline the oracle. Gate behind the real-WGS differential;
  keep the bounded in-pass buffer fallback ready.
- `Decision:` ADOPT WITH CHANGES
- `Notes:` also fetch mates' SA-tag supplementary loci

### SPD-5 — Cheap internals *(speed R5a–e)*
Integer contig compare, `u8` quals, `FxHashMap`, record reuse, skip non-main contigs at the reader.
- **Finding (#10):** all safe except **R5b (`u8` quals)** — consensus scores (`filters.rs:117`,
  `delta_best`) sum qualities past 255, so the read-quality/score-track split must be airtight.
- **Finding (`2nd-review`):** 5a/5c/5d verified byte-identical (name resolved once per contig
  is identical to the per-read string; the `get_mates`/`find_mates` maps are used for membership
  only, and output order comes from BAM coordinate order + `union.sort()`, not hash order). But
  **R5e is self-contradictory.** The current filter is exactly `ref_name.len() > 5` plus
  `MT`/`chrM` (`discovery.rs:202-207`). "Don't iterate alt/decoy/unplaced/HLA/MT" is
  byte-identical **only if** it skips exactly that set — yet 5e cites SPEC-5/SENS-4's
  header-driven allowlist, whose whole purpose is to *stop* dropping `NC_000014.9`-style main
  chromosomes that `len>5` currently discards. So 5e is either byte-identical (replicate `len>5`
  precisely, but then it is *not* the allowlist) or it is the allowlist (adds calls on
  RefSeq/T2T BAMs, breaking the speed doc's own guardrail) — not both. The divergence is latent
  on the `chr13`-style test BAM, so the current differential harness would not catch it.
- **Recommended:** **ADOPT** 5a/5c/5d; **ADOPT WITH CHANGES** 5b (isolate behind its own
  differential re-test, do last) and **5e** (treat it as the *behaviour-changing* SPEC-5/SENS-4
  item, not a byte-identical speed win — do it once, there, behind a toggle; don't smuggle it in
  as "cheap internals").
- `Decision:` ADOPT
- `Notes:`

---

## 2. Observability & safe refactors — do first, low/zero risk

### OBS-1 — `Breakpoint.stats` reject-counter sidecar *(impl R4 / sens R4)*
Per-`continue`/`return None` counters → `<out>.stats.json`.
- **Finding:** none — this is the measurement instrument for everything else. **But note
  finding #6:** counters show *what changed*, not *whether truth was removed*; do not mistake a
  falling count for "artefacts removed."
- **Finding (`2nd-review`):** the two sources spec this incompatibly — reconcile before
  building. Python's real schema (`breakpoint.py:37-57`, printed at `discovery.py:202`) is
  per-side `{too_few, rescued_pA, excluded, too_few_after_filter, clipped_failed,
  unclipped_failed, polymer, passed}`; **sens R4 mirrors it, but impl R4 invents different names
  and omits three of these** — take the field set from Python. Also, impl R4's justification
  "single-threaded, so plain fields suffice" is **falsified by SPD-3** (rayon per-contig): if
  SPD-3 lands, use atomics or per-thread accumulation (sens R4's choice), so decide OBS-1's
  storage *after* SPD-3's disposition.
- **Recommended:** **ADOPT — do this first of all** (Python field names; atomics if SPD-3 is on).
- `Decision:` ADOPT
- `Notes:`

### OBS-2 — Extract cluster/TSD constants, byte-identical *(impl R3, step 1)*
Lift `6`, `2`, `40`, `12`, `120` into named `config.rs` constants with identical values.
- **Finding:** none (zero-risk refactor).
- **Recommended:** **ADOPT.**
- `Decision:` ADOPT
- `Notes:`

### OBS-3 — Reconcile the 6 bp cluster window vs the 40 bp TSD window *(impl R3, step 2)*
Actual value change so a 6–40 bp TSD is not clustered inconsistently.
- **Finding:** legitimate, but note single-linkage chaining (`discovery.rs:152`) already lets a
  group span >6 bp; widening the cluster window interacts with SENS-1 (windowed count) — decide
  them together, not independently.
- **Recommended:** **ADOPT WITH CHANGES** — tune under OBS-1 stats + a toggle, jointly with SENS-1.
- `Decision:` ADOPT WITH CHANGES (see above)
- `Notes:`

### OBS-4 — Runtime config parameters *(sens R0)*
Load the `discovery` block from file/CLI; thread a `DiscoveryConfig`. Pure enabler.
- **Finding:** none, but flags a real latent bug: Rust hardcodes `MIN_MAPQ=40` while prod
  `config_hs.py`/`config_mm.py` use `60`. Pin this explicitly in any real-WGS run or the port
  is only faithful to the *generic* config.
- **Recommended:** **ADOPT** — prerequisite for tuning SENS-* honestly.
- `Decision:` ADOPT
- `Notes:`

### OBS-5 — Decide `extend_mates` fate *(impl R10)*
Rust already omits it (Python is a silent no-op). Either implement real mate-extension or
formally delete.
- **Finding:** longer contigs can introduce new artefacts; needs a real-WGS benchmark first.
- **Recommended:** **DEFER** — no code change until the benchmark exists; keep flagged.
- `Decision:` DEFER
- `Notes:`

---

## 3. Specificity / artefact rejection (`implement_review.md`)

### SPEC-1 — SMS both-ends-clip reject *(impl R15)* ⚠ top real-detection risk
- **Finding (#5, and the single most dangerous change):** the "safe for repeat clips" proof
  covers the *insert* clip but not the **anchor** clip. A true junction read anchored in a
  divergent young repeat / LTR flank gets ≥8 bp incidental anchor-end soft-clip → looks SMS →
  dropped; under the 2-read floor that deletes the event. Worst for the mouse ERV/IAP class of
  the IAPEz anchor, in the genome with 10× palindrome load. Slated first + on-by-default. "Three
  callers converge" is direction-only; the `MAX_CLIP_CLIP_LEN=8` comparator is underspecified.
- **Finding (`2nd-review`) — the code is fine, the science is the risk:** the clip-classification
  description (`discovery.rs:217-231`: single-sided → that side; dual-soft → longer wins, equal →
  `None`) is exact, and the proposed snippet is semantically correct — `left_len.min(right_len)
  >= MAX_CLIP_CLIP_LEN → continue` drops only genuine dual-clips while preserving "longer clip
  wins" for a junction read with a small stray anchor-end clip. Two splice caveats: (i) it must
  go **inside** the existing `} else if right_soft && left_soft {` arm (`discovery.rs:221-228`),
  not as the bare `if` the doc prints, or it breaks the `if/else` chain; (ii) `continue` targets
  the `for` loop over records — valid inside the `let clip = if…` expression since `continue`
  is `!`. So the residual risk is entirely the anchor-clip science in finding #5, not the code.
- **Recommended:** **ADOPT WITH CHANGES** — ship **OFF by default**; gate on the *shorter* clip
  with the tolerance already proposed; **prove non-removal on the real IAPEz anchor reads**
  before it is ever defaulted on; add the anchor-side caveat.
- `Decision:` DON'T ADOPT
- `Notes:`

### SPEC-2 — Low-complexity filter wired into discovery *(impl R2)*
- **Finding (#9):** premise correct (n-polymer runs only on `unclipped_cons`, `model.rs:228`),
  but the cheap "≤2 distinct bases" test + the mandatory poly-A escape hatch **nearly cancel** —
  it mostly catches poly-A that the hatch then exempts. Near-no-op as specified.
- **Recommended:** **ADOPT WITH CHANGES** — skip the ≤2-base form; implement the k-mer-entropy
  version (`Σ_K|freq−0.25^K| > 0.83`) with the poly-A escape. ERV clips are GC-diverse so this is
  safe for the anchor.
- `Decision:` ADOPT WITH CHANGES use k-mer entropy version with poly-A escape.
- `Notes:`

### SPEC-3 — Local depth / MAPQ / SMS pileup mask *(impl R1)*
- **Finding (#2, #3):** genuinely one of the richest artefact sources; the *right* guard for the
  sensitivity relaxations. But it is Tier-2 while the relaxations it guards are Priority-1 — the
  ordering leaves the relaxations unguarded. SMS-fraction reuses SPEC-1.
- **Recommended:** **ADOPT WITH CHANGES** — promote to ship *alongside* SENS-1/2/3, not after;
  it is the missing real gate.
- `Decision:` ADOPT WITH CHANGES
- `Notes:` use > 5x median coverage, samples from random 3000 locations as cutoff.

### SPEC-4 — Coverage-adaptive evidence threshold *(impl R16)*
Replace fixed `MIN_EVIDENCE_READS_PER_BREAKPOINT` with an xTEA-style cov→count lookup.
- **Finding:** sound and directly counters SENS-1's relaxation in pile-ups. Depends on the
  SPEC-3 depth estimator.
- **Recommended:** **ADOPT WITH CHANGES** — bundle with SPEC-3; keep the low-coverage end from
  dropping the 2-read floor below 2 for hallmark-bearing events.
- `Decision:` ADOPT WITH CHANGES, see recommendation above.
- `Notes:`

### SPEC-5 — Contig allowlist + exclude-BED plumbing *(impl R5 / R7 plumbing)*
- **Finding:** same underlying fix as SENS-4; do once. Low risk.
- **Recommended:** **ADOPT** — build the allowlist + sorted interval index; coordinate with the
  indexed reader (SPD-3).
- `Decision:` ADOPT
- `Notes:`

### SPEC-6 — Palindrome/IVR reference blacklist, per species *(impl R7 blacklist)*
- **Finding:** the review itself downgrades this to "second line behind SPEC-1"; generating it is
  a separate offline O(N²) pipeline; ArtifactsFinder PS length gate ships disabled.
- **Recommended:** **DEFER** — only build if residual palindrome FPs survive SPEC-1 + SPEC-3.
- `Decision:` DEFER
- `Notes:`

### SPEC-7 — RepeatMasker young-copy self-mask *(impl R17)*
Drop a candidate whose breakpoint falls in a same-family element below a divergence cutoff.
Tracks (`hs1.repeatMasker.out.gz`, `hs1.rte.out.gz`) are already in-repo.
- **Finding:** must **gate on divergence, not family membership**, or it deletes locus-unique
  young insertions — and for mouse ERV the clip *is* the LTR consensus, so a family-membership
  cut would delete the target events. Human-only tracks are bundled; mouse needs its own.
- **Recommended:** **ADOPT WITH CHANGES** — divergence-gated only; measure; supply a GRCm39 track
  before enabling on mouse.
- `Decision:` ADOPT WITH CHANGES, see above. This has to work with GRCm39 and GRCh38 as well as hs1.
- `Notes:`

### SPEC-8 — Within-BAM clip-recurrence flag *(impl R14)*
- **Finding (#8):** within one mouse, K real IAP insertions share the identical LTR-consensus
  clip → within-BAM recurrence **cannot honour its own "ERVs must survive" caveat.** The
  discriminating signal (cross-individual recurrence) is unavailable within one BAM.
- **Recommended:** **ADOPT WITH CHANGES** — sidecar **flag only, never a drop**, for any cohort
  with LTR/ERV elements. Real recurrence filtering belongs in combine (cross-individual).
- `Decision:` DONT ADOPT.
- `Notes:`

---

## 4. Sensitivity / recall (`recommendations_sensitivity.md`)

### SENS-1 — Windowed evidence count *(sens R2)*
Count QC-passing reads within ±jitter of the modal position, not only at the exact mode.
- **Finding (#3):** premise verified (`most_common_first` counts exact positions; consensus
  machinery already delta-aligns). Real recall gain — but in the 2-read regime two unrelated
  artefacts within jitter now fabricate a breakpoint, and the named counterweight (SENS-5) is
  "never a gate." Interacts with OBS-3.
- **Recommended:** **ADOPT WITH CHANGES** — only alongside a real gate (SPEC-3 pileup mask
  and/or SPEC-4 adaptive floor); tie jitter to the reconciled cluster window (OBS-3).
- `Decision:` ADOPT
- `Notes:`

### SENS-2 — Lower MAPQ anchor floor 40→20–30 *(sens R3)*
- **Finding (#3):** MEIs live where anchors are MAPQ 20–40, so real recall lever. But the
  compensating `MAX_LOWQ_CLIP_RATIO=0.65` per-locus guard is more work than the relaxation and is
  effectively the SPEC-3 machinery; the "at minimum keep `maps_fully_elsewhere`" fallback is not
  enough on its own near repeats.
- **Recommended:** **ADOPT WITH CHANGES** — the low-MAPQ-clip-ratio guard ships in the *same
  commit*; validate against germline hets **and** a low-VAF truth set (VAL-1).
- `Decision:` ADOPT
- `Notes:`

### SENS-3 — Emit one-sided clipped candidates *(sens R5)*
- **Finding (#4):** scoped to discovery only, but combine already **unconditionally drops
  one-ended/poly-A-only loci** (the ARCH-2 `continue`), and combine's **≥4-per-100 bp density
  mask** can wipe the now-denser true loci. Recovered events dead-end or, worse, real two-sided
  events get wiped.
- **Recommended:** **DON'T ADOPT as scoped.** Revisit only after ARCH-2 lands and a density-mask
  exemption exists. (The "fix the loop tail truncation" sub-part of R5 is a separate, safe bug
  fix — see Notes.)
- `Decision:` DON't ADOPT, REVISIT LATER
- `Notes:` (consider adopting *only* the `while il<l && ir<r` tail-visit fix separately)

### SENS-4 — Replace `len(name) > 5` contig filter *(sens R9)*
- **Finding:** correct latent recall cliff (`NC_000014.9`-style naming zeroes recall silently).
  Same fix as SPEC-5.
- **Recommended:** **ADOPT** — explicit primary-assembly allowlist from the BAM header; do once
  with SPEC-5.
- `Decision:` ADOPT
- `Notes:`

### SENS-5 — Positive-hallmark scoring, non-gating *(sens R1 / impl R8)*
Additive poly-A×purity, TSD length 2–20, EN motif — element-class-aware, poly-A never required.
- **Finding (#7):** the ERV caveat (never require poly-A) is airtight and correct. **But** the EN
  motif's "190×" is *cohort-level enrichment*, rated "moderate, not per-event proof" by the
  review's own cheat-sheet — using it as per-event additive weight biases toward AT-rich loci and
  the motif is written 4 inconsistent ways. TSD-length reward is unvalidatable on MEIsimulator
  (fixed 15 bp) and on ERV (≈6 bp).
- **Recommended:** **ADOPT WITH CHANGES** — keep poly-A×purity and TSD-length; make EN a
  tie-breaker at most with one fixed motif definition; write to the sidecar, never gate. ADOPT WITH CHANGES, see recommendation
- `Decision:`
- `Notes:`

### SENS-6 — Loosen `find_polya` capture *(sens R6)*
Drop the `is_proper_pair` reject; purity instead of contiguous-12; config-drive `fi>6`/cutoff;
lower toward 7 bp.
- **Finding (#3, runner-up worst):** `find_polya` runs only on the **low-MAPQ path**
  (`discovery.rs:185`). Dropping `is_proper_pair` + purity + 7 bp floor **multiplies genomic
  poly-A false rescues**, worst in the A/T-rich mouse genome where §11.4 already flags the 18.48%
  poly-A-rescued rate as suspect. The plan cites §11.4 for the recall caveat but ignores its
  specificity warning about this exact path.
- **Recommended:** **ADOPT WITH CHANGES** — do **not** relax `is_proper_pair` and the run-length
  together; keep a real gate; **hold off entirely on the mouse config** until measured on a mouse
  truth set. Consider human-only initially.
- `Decision:` DONT ADOPT
- `Notes:`

### SENS-7 — Consensus tolerant of a single disagreement *(sens R7)*
- **Finding:** premise verified (`delta_best = s0−s1−s2−s3`, break on first tie, `filters.rs:115`).
  One error in one of two reads can truncate a real consensus below the length gate. Real recall
  fix; low FP risk (downstream DFAM/bowtie2 tolerate a few mismatches).
- **Finding (`2nd-review`) — narrower than stated:** with 2 reads the break fires only when the
  two quals are **exactly equal** (`s0 − s1 = 0`); a 1-point Q difference already passes and
  keeps the better base. So the "comparable Q" framing overstates the trigger — the loss is real
  but confined to exact ties, which slightly shrinks the expected recall gain. Doesn't change the
  disposition.
- **Recommended:** **ADOPT WITH CHANGES** — "best strictly beats second-best, continue past
  isolated ambiguity" rather than the K-consecutive variant; keep byte-identical toggle-off.
- `Decision:` ADOPT WITH CHANGES, see above.
- `Notes:`

### SENS-8 — Short clips allowed when pure poly-A *(sens R8)*
Permit clips down to ~7 bp only when a pure poly-A/T terminus; keep the 12 bp floor otherwise.
- **Finding (#3):** compounds with SENS-6 on the same poly-A false-rescue surface; safe only if
  "pure poly-A" stays strict and element-class-gated.
- **Recommended:** **ADOPT WITH CHANGES** — implement *after* SENS-6 is measured; strict purity;
  human first.
- `Decision:` ADOPT
- `Notes:`

---

## 5. Architectural / out-of-discovery-scope

### ARCH-1 — SA-bridge → pseudo-discordant leg *(impl R11 / sens R10)*
Convert an SA-tag supplementary clip into a pseudo-discordant; require convergence with
clip/poly-A evidence; add xTEA clip↔disc geometric-consistency gate.
- **Finding:** the biggest sensitivity ceiling, but real engineering; scaffolding (SA parse,
  `bp_precise`) exists.
- **Finding (`2nd-review`) — "scaffolding" is thinner than advertised:** SA-tag parsing is real
  (`read.rs:99-107`, consumed by `maps_fully_elsewhere` and the exclude logic), but `bp_precise`
  (`model.rs:19`) is **inert** — it defaults `true` (`model.rs:50`), is only ever *read*
  (`model.rs:140,171`), and is **never set `false` anywhere in the port.** It is a vestigial
  field, not working discordant machinery; a real SA-bridge must add the code that sets and
  honours it. Also note RNEXT/PNEXT (needed for the discordant geometry) are still undecoded in
  `BamRead` — same gap SPD-4 depends on. So this remains genuine engineering, not "wiring."
- **Recommended:** **DEFER** — after Tier-1/2 land and are validated.
- `Decision:` ADOPT
- `Notes:`

### ARCH-2 — Poly-A-only locus fate in combine *(impl R9)*
The unconditional `continue` that drops poly-A-only loci.
- **Finding (#4):** **blocks SENS-3** from ever producing output. Combine-step change.
- **Recommended:** **DEFER** (but sequence it *before* SENS-3 if SENS-3 is ever revived).
- `Decision:` DEFER
- `Notes:`

### ARCH-3 — Delly/MELT/xTEA companion run *(impl R12)*
- **Recommended:** **DEFER** — orchestration, not discovery.
- `Decision:` DEFER
- `Notes:`

### ARCH-4 — MEI-polymorphism DB subtraction *(impl R13)*
- **Recommended:** **DEFER** — combine/genotype; PEAR-TREE's cohort genotyping already removes
  non-varying germline, so additive not essential.
- `Decision:` DEFER
- `Notes:`

### ARCH-5 — Cohort-adapted final classifier *(impl R18)*
- **Finding (#7):** if built, replace matched-normal features with cohort/tree features and
  **gate poly-A features on element class** or every mouse ERV is down-weighted.
- **Recommended:** **DEFER** — depends on SENS-5/SPEC-8 feature emission first.
- `Decision:` DEFER
- `Notes:`

---

## 6. Validation — a cross-cutting blocker

### VAL-1 — The validation harness cannot currently prove specificity *(finding #6)*
The §8.4 plan (stats diff + insertion count on one real BAM + germline-het check) **cannot
distinguish artefact-removal from truth-removal:** no artefact/truth labels, one BAM can't
estimate a rate, germline hets are high-VAF (blind to the low-VAF regime the relaxations
affect), MEIsimulator is blind to EN/ERV/twin-priming and fixes TSD at 15 bp, and the named S1
LINE-call fallback (`chr3:146970723`, `chr5:40029918`) is **unresolved, not confirmed truth** —
circular. Byte-identical-when-off proves the OFF state only; the 11-toggle ON-state space is
never tested against truth.
- **Finding (`2nd-review`) — the two plans share one oracle and their postures collide:** the
  speed doc's guardrail ("every change stays byte-identical to the Python oracle") and the
  sensitivity doc's admission ("any behaviour change necessarily breaks byte-identity") both key
  off the *same* `tests/differential_test.sh` vs. Python `src/discovery.py`. The moment any
  SPEC-*/SENS-* item is defaulted **on**, that differential goes red and can no longer validate
  the speed refactors. The only coherent ordering is the one §7 already implies but should state
  explicitly: **land every byte-identical speed/refactor item first while the oracle is green,
  then fork a second baseline** (Python-with-the-same-toggles, or a frozen Rust reference) for
  the behaviour work. Do not run SPD-3/SPD-4 validation against a tree that has any behaviour
  toggle on by default.
- **Recommended:** **ADOPT WITH CHANGES — treat as a gate on every SPEC-*/SENS-* item.** Before
  any relaxation or filter is defaulted on: (a) a low-VAF labelled spike-in truth set (patch
  MEIsimulator TSD to vary; add an ERV/LTR model or use a real mouse ERV truth set), (b)
  orthogonal (long-read/PCR/IGV) confirmation of a sample of *new* calls, (c) at least one
  combined-toggle ON-state run, not just one-at-a-time.
- `Decision:` ADOPT with the recommended changes
- `Notes:`

---

## 7. Recommended sequencing (if you adopt the above dispositions)

1. **OBS-1** (stats sidecar) + **OBS-2** (constant extraction) + **OBS-4** (runtime config) —
   zero-risk instruments and enablers.
2. **SPD-1/2/5a,c,d** — byte-identical speed wins; establish the real-WGS differential (VAL prereq).
3. **SPEC-5 / SENS-4** (contig allowlist + exclude-BED, done once).
4. **VAL-1** — stand up the low-VAF truth set and orthogonal check *before* any behaviour change
   is defaulted on.
5. **SPEC-3 + SPEC-4** (pileup mask + adaptive floor) — the *gates* the relaxations need.
6. **SENS-1 + SENS-2 + SENS-7** — recall levers, each shipped *with* its gate and OFF by default.
7. **SPEC-1** (SMS) — OFF by default, proven on the IAPEz reads first.
8. **SENS-5** (hallmark score, sidecar, non-gating) + **SPEC-8** (recurrence flag, never drop).
9. **SENS-6 / SENS-8** (poly-A relaxations) — human-only until mouse-measured.
10. **SPEC-7** (RM self-mask), **SPEC-6** (blacklist) only if residual FPs survive.
11. **SPD-3 / SPD-4** (parallelism + mate fetch), **SPD-5b** — behind the differential.
12. Architectural (**ARCH-1..5**), **SENS-3**, **OBS-5** — deferred, benchmark-gated.

*Every behaviour change ships behind an env toggle (as `PEARTREE_KEEP_FULLMAP` already does),
defaults consistent with the current output until VAL-1 clears it, and is attributed via OBS-1.*

---

## 8. Feature A / B — discordant anchoring & processed-pseudogene annotation

Added after the original register (full plan: `discordant_and_pseudogene_plan.md`).
Both ship OFF by default; default output stays byte-identical (differential_all).

### DA-1 — discordant-mate anchoring of one-sided junctions
- **Problem:** a junction with a real breakpoint on only one flank (classically a lone
  poly-A clip) is dropped — output pairing requires a real reciprocal side.
- **Change:** collect discordant read pairs (`discordant_anchor`); a real breakpoint
  with no reciprocal partner in its TSD window is paired with a cluster of ≥
  `discordant_min_reads` discordant mates on the missing side, emitted with a
  `disc_<pos>` token (no reads for that end).
- **Contract:** first discovery change touching the 4-step pipeline; Python
  `combine_insertions` parses `disc_` (`TYPE_*_DISC`) and parks these calls like
  poly-A insertions (downstream genotyping of coordinate-only ends deferred).
- **Decision:** ADOPT, opt-in + flagged. VAL-1: recall 0.93→1.0 on the one-sided
  class. Gate before defaulting on: real-WGS differential + orthogonal confirmation.

### DA-2 — RTE-origin of the discordant mates
- **Change:** score the *mate* landing site against a RepeatMasker track (reuse the
  SPEC-7 loader, divergence-gated); label by default, opt-in `discordant_rte_only`
  gate requiring RTE-origin ≥ `discordant_rte_min`, which drops random-SV clusters.
- **Decision:** ADOPT, label by default / gate opt-in. VAL-1: the gate restores
  precision by rejecting the non-RTE artefact class while keeping RTE-origin calls.

### PG-1 — processed-pseudogene (splice) annotation
- **Change:** non-gating `<out>.splice.tsv` (`splice_hallmark` + `exon_annotation`)
  flagging candidates whose mate reads span ≥ `splice_min_exons` exons of one gene
  with the introns skipped. Discovery detects the *insertion*; this annotates the
  *processed-ness*. Main output unchanged.
- **Decision:** ADOPT, sidecar-only. VAL-1: flags all synthetic pseudogenes, none of
  the single-exon controls. Needs a real exon annotation matching the BAM assembly.
