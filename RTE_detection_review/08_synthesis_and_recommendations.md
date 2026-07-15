# 8. Synthesis and recommendations for PEAR-TREE2

*This section pulls [§3](03_detecting_true_events.md)–[§7](07_peartree_code_review.md)
together into (a) a reference best-practice pipeline and (b) a prioritised, concrete
recommendation list for PEAR-TREE2. It is written to complement, not replace,
`PEAR-TREE2_PLAN.md` — where they overlap, that is noted.*

---

## 8.1 The best-practice short-read RTE pipeline (distilled)

From the four detection studies and the two artefact studies, a robust pipeline has six
layers. PEAR-TREE already implements layers 1, 3 (partly), 4 and 6 strongly; its gaps are
in layers 2 and 5.

1. **Dual read-level signal.** Collect *both* discordant pairs (one mate on element
   consensus, one unique) *and* soft-clipped/split reads. Cluster to localise; split to
   resolve. → PEAR-TREE has the split half only.
2. **Genome-free artefact rejection first** (cheap): adapter, poly‑G, low-complexity/
   homopolymer, **both-ends-clipped (SMS) reject** (R15), self-palindrome/inverted-remap,
   ≥ 2 concordant reads, local-pileup MAPQ/SMS/depth mask (R1). → PEAR-TREE is missing the SMS
   reject and has the pileup mask disabled.
3. **Positive-hallmark scoring:** poly‑A (length × purity), TSD (present, canonical length),
   EN motif, element-consensus match. A candidate's score should *rise* with hallmarks, not
   only fall with artefacts.
4. **Genome-aware artefact rejection:** clipped-maps-locally, end-maps-entirely,
   reference IVR/palindrome blacklist, telomere/centromere/segdup blacklist.
5. **Classification:** matched-normal and/or cohort-recurrence + MEI polymorphism databases
   (1000G MEI, dbRIP, euL1db) + population-allele-frequency panel.
6. **Clonality / genotyping:** VAF ≈ 0.5 for somatic in clonal material; population
   genotype-matrix consistency (require wt + insertion partition); orthogonal validation
   (long-read/PCR/IGV).

## 8.2 Prioritised recommendations for PEAR-TREE

Ranked by (impact on precision/recall) ÷ (effort). "Plan §" cross-references
`PEAR-TREE2_PLAN.md`.

### Tier 1 — high impact, low effort (do first)

**R1. Turn discovery-time high-coverage masking back on.** `max_read_count` is unreferenced
and the genotype high-coverage check is `if False and ...` ([§7](07_peartree_code_review.md)
finding 1). Coverage spikes are one of the richest artefact sources ([§5.6/§4.2.6](05_sequencing_artefacts.md)).
Add a running per-window depth estimate in discovery and drop candidates in spikes.
*(= Plan §2.3, "biggest single artefact-reduction win".)* **▲ Code-grounded form
([§12.6](12_tool_implementations_compared.md)):** the mature version is a **local-pileup mask
that needs no blacklist** — recompute, per candidate locus, the fraction of low-MAPQ reads and
the fraction of SMS (both-ends-clipped) reads and/or local depth vs the sample median. Concrete
thresholds from the field: v-TraFiC drops a locus unless **<30 % reads MAPQ<10 and <15 % SMS**
(`filterSMSclusters.sh`); MEIGA thresholds `percMAPQ`/`percSMS` over `bkp±100`
(`bins_lowMAPQ_SMS.py` + `area` filter); xTEA masks **depth > 3× sample median**
(`MAX_COV_TIMES=3`, median from 3,000 random sites). Any of these is cheaper than a static
blacklist and directly targets the multimapping/palindrome pile-ups.

**R2. Wire in the low-complexity filters that already exist.** `is_low_complexity` and
`has_well_defined_breakpoint` are implemented in `sequence_checks.py` but never called on
the discovery path ([§7](07_peartree_code_review.md) finding 7). Call them before emitting a
breakpoint. Near-zero effort. *(= Plan §2.3.)* **▲ Code-grounded form
([§12.2–§12.3](12_tool_implementations_compared.md)):** two field-proven algorithms to wire in
*before* clustering — v-TraFiC runs NCBI **`dustmasker`** on candidate-read sequences and drops
wholly low-complexity reads; MEIGA-LR uses **k-mer entropy** (drop if `Σ_K|freq−0.25^K| > 0.83`
over K∈{1,2,3,4}) plus a microsatellite check (`|obs−0.25^len| > 0.2`). **Both keep a poly‑A
escape hatch** (MEIGA: `poly‑A>0.4` auto-passes) so genuine poly‑A tails are not deleted — the
same exemption [§10.2](10_empirical_discovery_patterns.md) insists on.

**R3. Reconcile the two clustering windows.** Discovery clusters clipped reads at a
hard-coded **6 bp** while the documented/TSD window is **`max_bp_window = 40`**
([§7](07_peartree_code_review.md) finding 2). A 6–40 bp TSD can be clustered inconsistently.
Make the 6 configurable and reconcile with the TSD logic in `output()`.

**R4. Emit the `Breakpoint.stats` audit per run.** The per-side counters (too_few / excluded /
polymer / clipped_failed / passed / rescued_pA …) already exist and are a ready-made
discovery-time artefact dashboard ([§7](07_peartree_code_review.md) §7.3.3). Write them to a
sidecar JSON so specificity regressions are visible run-to-run. Trivial effort, big
observability payoff — and it operationalises the Plan's "every change is measured" principle.

**R5. Replace the `len(name) > 5` contig heuristic with an explicit include/exclude BED.**
([§7](07_peartree_code_review.md) finding 6). This simultaneously fixes non-`chrN`
assemblies and gives you the blacklist mechanism R7 needs. *(= Plan §2.4 config work.)*

### Tier 2 — high impact, medium effort

**R6. Add a self-palindrome / inverted-repeat discovery filter.** The top artefact for a
clipped-read caller, especially on enzymatic-prep data ([§5.1](05_sequencing_artefacts.md),
[§4.2.1](04_artefacts_in_discovery.md)). Genome-free version: reject when the clipped segment
equals the reverse complement of the adjacent anchor (or has ~50% clip ratio with the
junction ≤ 30 bp from the read edge, at a palindrome centre). PEAR-TREE's `SA`-based
same-contig exclusion catches the *split* case but not the *single-read self-fold* case.
**▲ Largely superseded by R15 ([§12.6](12_tool_implementations_compared.md)):** the cheapest,
field-standard version of this filter is the read-level **SMS reject** (drop any read
soft-clipped on *both* ends), which MEIGA, v-TraFiC and xTEA all implement. Do R15 first; keep
the revcomp-of-anchor check as a complement for the single-clip self-fold case R15 does not
cover.

**R7. Add reference blacklists.** (a) A **telomere/centromere/segdup** BED (Delly's `-x`
idea). (b) An **inverted-repeat/palindrome** BED generated à la ArtifactsFinder — this raised
sonication/enzymatic concordance from 7.8% → 80.4% in Chen et al. Drop candidates inside either.
Depends on R5's BED plumbing. **▲ Corrected from the code
([§12.5/§12.8](12_tool_implementations_compared.md)):** ArtifactsFinder does **not** emit an
interval BED — it writes per-base artefact *positions*; the real IVS params are **arm‑pair ≥8 bp
(`D_LEN`), spacer ≥5 bp (`STEM_LEN`), sub‑arm ≥2 bp (`S_LEN`), ±50 bp**, and its palindrome
minimum-length gate ships **disabled** (re-enable `LEN≥17`). A reproducible genome-wide recipe
(tile → `bedReform` 200 bp windows → run both finders → collapse arms/sites → `bedtools slop 15`
+ merge) is in [§12.8](12_tool_implementations_compared.md). Treat this as a **second line behind
R15**, and build it **per species** (human T2T *and* mouse GRCm39 — mouse has ~10× the palindrome
load, [§11](11_cross_dataset_and_mouse_erv.md)).

**R8. Score the poly‑A/TSD/EN hallmarks explicitly — species- and element-aware.** Today
PEAR-TREE uses TSD *geometry* and poly‑A presence, but does not reward **poly‑A length ×
purity**, **canonical TSD length (2–20 bp)**, or the **EN motif (TT/AAAA)** at the 5′ nick
([§7](07_peartree_code_review.md) finding 8; [§3.3](03_detecting_true_events.md)). A small
additive positive score would let you *keep* well-supported events at lower coverage while
*distrusting* hallmark-free clips — raising both recall and precision. The EN-motif
enrichment (190× in Nam 2023) makes it a cheap, strong prior. **Critical caveat from
[§11.4](11_cross_dataset_and_mouse_erv.md): make this element-class-aware and never require
poly‑A globally.** poly‑A + `TT/AAAA` apply to human L1/Alu and mouse L1Md/SINE; **mouse
ERV/LTR elements (IAP, MusD/ETn, MMERVK) produce NO poly‑A** and are scored instead by
**LTR-consensus identity + a short (~6 bp) TSD**. A poly‑A-weighted score applied blindly
would penalise — and a poly‑A *filter* would delete — every genuine ERV insertion.

**R9. Decide the fate of poly‑A-only insertions in `combine`.** The unconditional `continue`
that skips all poly‑A-only loci ([§7](07_peartree_code_review.md) finding 4) discards the
sensitivity that the discovery-time poly‑A rescue works hard to produce — costly for
5′-truncated L1s and Alu/SVA. Implement the "heavy filter" the comment anticipates
(e.g. require poly‑A length + unique anchor + cohort recurrence) rather than dropping them
wholesale.

**R10. Fix or remove `extend_mates()`.** It currently iterates an emptied list (silent
no-op) ([§7](07_peartree_code_review.md) finding 3; Plan §1.3). Either point it at the final
breakpoint lists (longer clipped contigs → better DFAM/RepeatMasker annotation) or delete it
— but decide with a real-WGS benchmark, since longer contigs can introduce new artefacts.

### Tier 3 — high impact, higher effort (architectural)

**R11. Add discordant-read-pair discovery.** PEAR-TREE's single biggest *sensitivity* gap:
it sees only events with a junction-spanning clipped read
([§7](07_peartree_code_review.md) finding 5; [§3.1](03_detecting_true_events.md)). Adding
TraFiC/Delly-style positive/negative discordant clusters (one mate on element consensus, one
unique) would recover large-TSD, low-coverage and unmappable-junction events. Use Delly's
**per-library insert-size model (3 SD)** rather than a fixed window. This is the deferred
discordant branch removed from `discovery.py`; re-introduce it as a *complementary* signal
that must still converge with the clipped/poly‑A evidence to be reported. **▲ Cheapest first
step ([§12.3](12_tool_implementations_compared.md)):** before a full discordant engine, adopt
MEIGA's **`SA_as_DISCORDANTS`** trick — convert an SA-tag supplementary clip into a pseudo-
discordant so a *single split read* substitutes for a pair; this fits the clip-first design with
minimal new code. Then borrow xTEA's **clip↔disc geometric consistency** check
(`_is_distance_consistency`): a left-clip is confirmed only if right-side discordant support maps
within one insert-size — xTEA's single strongest FP filter.

**R12. Consider a Delly (or MELT/xTea) companion run for non-RTE SVs and hard cases.**
Nam 2023 ran four callers in parallel precisely because each has blind spots. A generic SV
caller alongside PEAR-TREE catches RT-mediated translocations and junctionless insertions;
PEAR-TREE then supplies the RTE-specific poly‑A/TSD classification Delly lacks
([§6.3](06_delly_review.md)).

**R13. Optional: MEI-polymorphism-database subtraction.** For studies where a clean somatic
call set matters, intersect against 1000G MEI / dbRIP / euL1db and a population panel
([§3.4](03_detecting_true_events.md)). PEAR-TREE's cohort-genotyping already removes
non-varying germline insertions (`min_wild-types`), so this is additive rather than essential.

**R14. Treat recurrence-across-loci as an explicit filter.** *(Surfaced by the empirical
scan — [§10.3](10_empirical_discovery_patterns.md).)* Across the **full 560-file**
JAK2/HNRNPA1 discovery output, **~45% of all clipped sequences are non-unique** (recur ≥5×
across the cohort), dominated by poly‑A/T homopolymers and **Alu-consensus** fragments
(reference/polymorphic Alu or mismapping between Alu copies; ~184k Alu vs ~7k L1 instances). A
clip that is identical at many independent loci is by definition **not a novel locus-specific
insertion**. PEAR-TREE currently removes most of these only in `combine` via the bowtie2
end-to-end remap; flagging high-recurrence clip sequences at discovery (a hash count, or a
small bundled repeat-consensus screen) would drop a large share of this recurrent load
(~45% of clips are non-unique; Alu-consensus alone ~184k instances) earlier and cheaper.
**Caveat 1:** filter on *sequence recurrence*, not on "is a repeat" — a **locus-unique** Alu
clip with a poly‑A tail and TSD is a candidate somatic Alu insertion and must survive. This
is **doubly essential for mouse ERVs** ([§11.4](11_cross_dataset_and_mouse_erv.md)), whose
true insertions have clipped ends that *are* the IAP/LTR consensus — "is a repeat" would
delete exactly the events sought.
**Caveat 2 — score recurrence per individual, using the folder/prefix structure**
([§11.3](11_cross_dataset_and_mouse_erv.md)): a clip shared across **many individuals** is
reference-repeat/artefact (drop); a clip shared across **many colonies of one individual** is
a **germline** call (keep). In the 15-individual human set, 3,692 clips appear in all
individuals (artefact) while 781,044 (~47%) are private to one individual (candidate
germline). A naïve *global*-recurrence filter would delete the germline events the study
exists to find.

### Tier 4 — surfaced by the code-level review of seven tools ([§12](12_tool_implementations_compared.md))

*These four recommendations come from reading the actual source of TraFiC, v‑TraFiC, MEIGA‑SR/LR
and xTEA (not their papers). R15 is the single cheapest high-value item in the whole review.*

**R15. Reject "both-ends-soft-clipped" (SMS) reads at discovery.** *(High impact, ~one line.)*
A read soft/hard-clipped on **both** flanks is the read-level fingerprint of a reference
palindrome/inverted repeat, an adapter-dimer, or a spurious multimap — a structure-specific
chimera ([§5.1](05_sequencing_artefacts.md)). **Three independent mature callers all reject it**:
MEIGA (`GAPI/bamtools.py:654`, CIGAR op 4/5 at first *and* last), v‑TraFiC (discovery-awk CIGAR
`^S*M*S$`), xTEA (`MAX_CLIP_CLIP_LEN=8`). PEAR-TREE currently keeps such a read by classifying it
on its *longer* clip (`discovery.py`, "longer clip wins"). Adding an SMS reject is the cheapest,
most convergent artefact filter in this review and largely subsumes R6/R7 for the enzymatic-prep
palindrome problem. **It keys on CIGAR shape, not "is a repeat," so it is safe for genuine
repeat-consensus clips** — including the mouse ERV events of [§11](11_cross_dataset_and_mouse_erv.md),
for which mouse's ~10× palindrome load makes R15 *higher* priority.
([§12.6](12_tool_implementations_compared.md).)

**R16. Make clip/evidence-count thresholds coverage-adaptive.** *(Medium effort.)* PEAR-TREE uses
fixed minimum read counts; xTEA scales them with local coverage via a lookup table
(`x_parameter.py`: Illumina germline cov 5→(1,3,0), 30→(3,4,1), 100→(8,12,3), 300→(25,30,8);
lower table for case-control). A fixed threshold is simultaneously too lax in high-coverage
pile-ups (false positives) and too strict at low coverage (missed real events). A
coverage-indexed `min_evidence_reads_per_breakpoint` fixes both, and pairs naturally with the R1
depth estimate. ([§12.4](12_tool_implementations_compared.md).)

**R17. Add a RepeatMasker "young reference copy" self-mask.** *(Medium effort.)* A distinct
artefact class from germline recurrence: reads **mis-donated by a *young* reference element of the
same family** near the breakpoint. TraFiC drops a call whose breakpoint falls in a same-family
RepeatMasker element with **divergence ≤20%** (`TEs_poly_cleaner.pl`); xTEA uses **<15%**
(`REP_DIVERGENT_CUTOFF`). This catches mismapping a static IR/germline blacklist misses. Gate on
*divergence*, not mere family membership, so a locus-unique young insertion still survives.
([§12.2/§12.4](12_tool_implementations_compared.md).)

**R18. Consider a cohort-adapted final classifier.** *(Higher effort; optional.)* Both mature
short-read callers end in a trained model — MEIGA a **logistic regression over 29 features**
(coefficients in [§12.4](12_tool_implementations_compared.md); poly‑A and clip support are the
strongest positives, `svsNormal`/`nbNormal`/`areaSMS` the strongest negatives), xTEA a
**random forest over 15 coverage-normalised features**. PEAR-TREE uses a hand-tuned
quality-weighted score. If a classifier is added, reuse these feature schemas but **replace the
matched-normal features (`svsNormal`, `nbNormal`, `germPerc`) with cohort/tree-derived features**
(cross-donor sharing fraction, tree-consistency, per-branch recurrence — the R14 signals). A
simple regularised logistic model with class weighting is enough; no deep model needed.
**Element-class caveat ([§11](11_cross_dataset_and_mouse_erv.md)):** MEIGA's model over-weights
poly‑A; gate poly‑A features on element class or every true mouse ERV insertion is down-weighted.

## 8.3 What PEAR-TREE should *not* change

- **Keep the TSD-first, clip-first, poly‑A-as-first-class design.** It matches the biology
  and yields base-resolution breakpoints natively.
- **Keep the layered genome-free → genome-aware → cohort ordering.** It is exactly the
  defence-in-depth the artefact literature argues for and is efficient.
- **Keep population/phylogenetic genotyping as the primary somatic discriminator.** For
  clonal cohorts it is at least as powerful as a single matched normal, and it is
  PEAR-TREE's distinctive strength.
- **Keep the on-disk step contract** while doing all of the above (as the Plan already
  insists) so each change is independently shippable and measurable.

## 8.4 Suggested validation for any change (per the Plan's "every change is measured")

1. **Correctness gate:** the provided test BAM must keep producing the documented IAPEz call
   at `13:32992169-32992177`.
2. **Specificity:** track `Breakpoint.stats` and final insertion count on one real WGS BAM
   before/after each filter (R1–R8). A filter that removes candidates must be shown to remove
   *artefacts*, not true calls — check against known germline heterozygous insertions
   (the `check_germline_coverage.py` tool is built for exactly this sensitivity check).
3. **Recall:** for R9/R11 (sensitivity changes), confirm recovery of a spiked or known
   insertion set that the current pipeline misses. **▲ Concrete spike-in engine
   ([§12.9](12_tool_implementations_compared.md)):** use **MEIsimulator** to generate a merged
   BAM + labelled truth report (seeded, reproducible), then sweep **VAF** (`--clonality`
   5/10/25/50/100 %, the key clip-first recall-vs-VAF curve), **coverage** (15/30/60×),
   **5′ truncation / element size**, and **event class** (solo vs partnered/orphan transduction
   vs SV). Caveats to document: its TSD is a fixed 15 bp, solo-insertion poly‑A is fixed, and it
   models **no EN motif, no 5′ twin-priming inversion, and no ERV/LTR** — so validate EN-motif,
   inverted-truncation and **mouse ERV** ([§11](11_cross_dataset_and_mouse_erv.md)) recall on real
   data instead.
4. **Orthogonal truth where available:** long-read/PCR/IGV confirmation of a sample of new
   calls, as every reference study does.
