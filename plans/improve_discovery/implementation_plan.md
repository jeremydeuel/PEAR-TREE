# PEAR-TREE2 discovery — implementation plan

Derived from the `Decision:` fields in [`discovery_plan.md`](discovery_plan.md), with the
four blocking decisions resolved (2026-07-15):

1. **Validation of behaviour changes is Rust-only, no equivalence reference.** The existing
   Python differential ([`tests/differential_test.sh`](rust/peartree-discovery/tests/differential_test.sh))
   remains the oracle **only** for byte-identical speed/refactor work (all behaviour toggles off).
   Behaviour toggles are validated by VAL-1 metrics + Rust unit tests + toggle-OFF regression
   snapshots. Python `src/discovery.py` is **not** kept in sync with the science.
2. **ARCH-1 deferred** (overrides its ADOPT mark) — dead-ends without ARCH-2, needs discordant
   infra + VAL-1.
3. **VAL-1 built in full by me**, synthetic. ⚠ Synthetic-only: blind to real EN/ERV/twin-priming
   and (unless modelled) fixes TSD at 15 bp — this is the known ceiling on what it can prove.
4. **SENS-5 = ADOPT WITH CHANGES.**

Global rule (from the register's closing note): every behaviour change ships behind a config
toggle (env-overridable, like `PEARTREE_KEEP_FULLMAP`), **defaults to current output**, and is
attributed via the OBS-1 stats sidecar. Nothing defaults ON until VAL-1 clears it.

---

## Phase 0 — Instruments & enablers (byte-identical, zero risk)

| Item | Task | Files | Validation |
|---|---|---|---|
| **OBS-2** | Lift `6, 2, 40, 12, 120` cluster/TSD windows into named constants, identical values | `config.rs`, `discovery.rs:152,442,444,447,454` | Python differential green |
| **OBS-4** | `DiscoveryConfig` struct loaded from `--config` file / CLI, threaded through `Discovery`. Surface the **`MIN_MAPQ` 40-vs-60** gap (Rust hardcodes 40; `config_hs.py`/`config_mm.py` use 60) and require it be pinned per run | `config.rs`, `main.rs`, `discovery.rs` | Differential green with defaults = current constants |
| **OBS-1** | `<out>.stats.json` reject-counter sidecar. **Python field set** ([`breakpoint.py:37-57`](src/breakpoint.py#L37)): per-side `{too_few, rescued_pA, excluded, too_few_after_filter, clipped_failed, unclipped_failed, polymer, passed}`. Plain fields now; atomics/per-thread deferred until SPD-3 | new `stats.rs`, `discovery.rs`, `filters.rs`, `polya.rs` | Counts reconcile against emitted record count on test BAM |

Do OBS-1 **first** — it is the measurement instrument for every behaviour change downstream.

---

## Phase 1 — Byte-identical speed wins (oracle stays green)

| Item | Task | Notes |
|---|---|---|
| **SPD-1** | Reorder `BamRead::from_record` so seq/qual/name/tags build only for surviving clipped candidates | Byte-identical standalone. (SENS-6 dropped, so the ordering hazard the register flagged is moot.) |
| **SPD-2** | Wire `bam_threads` to a multithreaded reader | Add `noodles-bgzf` at the version `noodles-bam 0.79` pins transitively (not an explicit mismatched version — type-conflicts against `bam::io::Reader::new`) |
| **SPD-5a** | Integer contig compare | Byte-identical |
| **SPD-5c** | `FxHashMap` for `get_mates`/`find_mates` maps | Membership-only maps; output order is BAM coord + `union.sort()`, so hash order is irrelevant |
| **SPD-5d** | Record reuse | Byte-identical |
| **VAL prereq** | Run `differential_test.sh` on a real WGS BAM on the cluster | Establishes the real-WGS byte-identity oracle the SPD-3/4 work later needs |

Deferred within SPD-5: **5b** (`u8` quals) → Phase 8, its own differential (consensus scores
exceed 255, [`filters.rs:115-117`](rust/peartree-discovery/src/filters.rs#L115)). **5e** is
**not** a byte-identical speed win → folded into SPEC-5/SENS-4 (Phase 2).

---

## Phase 2 — Contig handling, done once

| Item | Task | Files |
|---|---|---|
| **SPEC-5 + SENS-4 (= SPD-5e)** | Replace `len(name)>5` + `MT`/`chrM` skip with a **header-driven primary-assembly allowlist** + **exclude-BED** (sorted interval index). Recovers `NC_000014.9`-style RefSeq/T2T main chromosomes currently zeroed | `discovery.rs:202-207,419,421`, new `intervals.rs`, `config.rs` |

Behaviour-changing on RefSeq/T2T naming (not on `chr13`-style BAMs). Ships behind config;
coordinate the interval index with the indexed reader used by SPD-3.

---

## Phase 3 — VAL-1 truth infrastructure ⚠ BLOCKER

Gates every ON-default in Phases 4–7. Built by me, synthetic.

- **Spike-in truth set**: patch MEIsimulator to (a) vary TSD length (currently fixed 15 bp) and
  (b) add an ERV/LTR insertion model, so recall can be measured across element classes and low VAF.
- **Toggle-matrix runner**: run discovery over the labelled set across the toggle space
  (one-at-a-time **and** at least one combined-ON run), emitting recall/precision + OBS-1 deltas.
- **Orthogonal check hook**: structure for IGV/long-read/PCR confirmation of a sample of *new*
  calls (data supplied later; harness ready now).
- ⚠ **Documented limitation**: synthetic truth cannot prove specificity against real
  EN/ERV/twin-priming artefacts. Treat green VAL-1 as necessary, not sufficient.

---

## Phase 4 — The gates (land before / with the relaxations)

| Item | Task | Scope you set |
|---|---|---|
| **SPEC-3** | Local coverage mask: sample coverage at **3000 random locations**, drop breakpoints in windows **>5× median** | Coverage-only (SPEC-1 not adopted, so the SMS-fraction sub-part is dropped) |
| **SPEC-4** | Coverage-adaptive evidence floor (xTEA-style cov→count) built on SPEC-3's estimator | Never drops below the 2-read floor for hallmark-bearing events |

---

## Phase 5 — Recall levers (OFF by default, each with its gate)

| Item | Task | Coupling |
|---|---|---|
| **OBS-3 + SENS-1** | Reconcile 6 bp cluster vs 40 bp TSD window; windowed evidence count within ±jitter tied to the reconciled window | Decide jointly; gated by SPEC-3/4 |
| **SENS-2** | MAPQ floor 40→20–30 **with** `MAX_LOWQ_CLIP_RATIO≈0.65` per-locus guard in the **same commit** | Guard reuses SPEC-3 machinery |
| **SENS-7** | Consensus continues past *isolated* ties ("best strictly beats second-best") | Byte-identical toggle-off ([`filters.rs:115`](rust/peartree-discovery/src/filters.rs#L115)) |

---

## Phase 6 — Hallmark & poly-A (human-first)

| Item | Task |
|---|---|
| **SENS-5** | Additive poly-A×purity + TSD-length (2–20) to the **sidecar, non-gating**; EN motif a **tie-breaker only** with one fixed motif definition; **poly-A never required** (ERV-safe) |
| **SENS-8** | Short (~7 bp) clips allowed **only** on strict-pure poly-A terminus, element-class gated, human-first. (Decoupled from SENS-6, which is not adopted — so this is the conservative half only) |

---

## Phase 7 — RepeatMasker self-mask

| Item | Task |
|---|---|
| **SPEC-7** | Drop a candidate whose breakpoint falls in a same-family element **below a divergence cutoff** (divergence-gated, **not** family-membership — family-membership would delete mouse ERV targets). Off by default | 

Must support **GRCh38 + GRCm39 + hs1**. hs1 tracks in-repo (`hs1.repeatMasker.out.gz`,
`hs1.rte.out.gz`); **GRCh38 and GRCm39 `.out.gz` tracks need sourcing** before enabling on those.

---

## Phase 8 — Parallelism & mate fetch (behind the differential)

| Item | Task | Notes |
|---|---|---|
| **SPD-4** | Mate-coordinate fetch to kill the 2nd BAM pass, **+ also fetch each mate's SA-tag supplementary loci** (your change) | Restores MATE-line completeness the naive fetch would lose; **drop the byte-identical label**, re-baseline / document the delta. Adds RNEXT/PNEXT + SA decode to `BamRead` |
| **SPD-3** | Rayon per-contig + indexed reader, **after SPD-4** | Merge must reproduce the per-contig sort exactly; prove on the real-WGS differential. Triggers OBS-1 → atomics/per-thread |
| **SPD-5b** | `u8` quals, **last** | Isolate read-quality vs consensus-score tracks (scores >255); own differential |

---

## Independent — SENS-3 tail-visit fix only

Extract just the `while il<l && ir<r` tail-truncation fix ([`discovery.rs:441`](rust/peartree-discovery/src/discovery.rs#L441)),
**not** the full SENS-3 (one-sided emission — not adopted, dead-ends in combine). Verify against
Python first: if it changes output it's a gated behaviour change, not a free bug fix.

---

## Not in scope (your decisions)

- **Don't adopt**: SPEC-1 (SMS both-ends reject), SPEC-8 (clip-recurrence flag), SENS-6 (loosen `find_polya`).
- **Defer**: OBS-5 (`extend_mates`), SPEC-6 (palindrome blacklist), SENS-3-full, ARCH-1/2/3/4/5.

---

## Dependency ordering (critical path)

```
OBS-1/2/4  →  SPD-1/2/5acd + real-WGS differential  →  SPEC-5/SENS-4
                                                          │
                                    VAL-1 (blocker) ──────┤
                                                          ▼
                            SPEC-3 → SPEC-4  →  SENS-1/2/7 (+OBS-3)  →  SENS-5/8  →  SPEC-7
                                                          │
                                    SPD-4 → SPD-3 → SPD-5b (behind differential, any time after Phase 1)
```

Each phase is an independently revertible PR. Behaviour PRs (Phase 4+) do not begin defaulting
anything ON until VAL-1 is standing.
