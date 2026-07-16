# Adversarial scientific review of the PEAR-TREE2 discovery plans

*Hostile external referee report on `implement_review.md`, `recommendations_sensitivity.md`
and `recommendations_speed.md`, cross-checked against their cited source
(`RTE_detection_review/`) and against the actual Rust discovery source
(`rust/peartree-discovery/src/{discovery,model,filters,polya,config,read}.rs`).*

*Stance: refutational. For every substantive recommendation I tried to make the case that
it is a bad idea and report only what survived. The correctness anchor is the IAPEz (mouse
ERV) call at `13:32992169-32992177`. Findings are ranked by likelihood × impact on real
detection.*

---

## Bottom line up front

The three plans are internally literate and the code-grounding is mostly accurate — I
verified the funnel gates, the exact-8-mer poly-A rescue, the exact-position evidence count,
and the first-ambiguous-base consensus truncation all exist as described in the Rust source.
But the two behaviour-changing plans **pull the same gates in opposite directions, reuse
colliding `R#` labels for opposite-signed changes, and lean on "guards" that do not gate and
a validation harness that cannot tell truth-removal from artefact-removal.** Several
individually-plausible recommendations become dangerous specifically for the
mouse-ERV/IAPEz class that is the stated correctness anchor.

**Verification caveats I owe the reader:**
1. Every quantitative literature number (`MAX_CLIP_CLIP_LEN=8`, the cov→reads table,
   `MAX_LOWQ_CLIP_RATIO=0.65`, v-TraFiC `<30%/<15%`, RepeatMasker divergence ≤20%/<15%) is
   **internally consistent across the review docs but rests entirely on the author's private
   reading of seven repos**, presented with file:line precision I could not independently
   verify without the upstream source. The R15 comparator is flagged as underspecified below.
2. All code-behaviour claims were confirmed against the Rust source directly; those are solid.

---

## 1. The two plans collide on the same gates — and share `R#` labels for opposite changes
**Verdict: KILLS the "do X first" ordering as written. Confidence: certain (textual).**

`implement_review.md` uses the review's R1–R18; `recommendations_sensitivity.md` invents its
**own** R0–R10. They overlap lexically but mean opposite things:

| Label | implement_review (specificity ↑) | sensitivity (recall ↑) |
|---|---|---|
| **R1** | coverage / low-MAPQ **mask** (drop candidates) | hallmark **score** (rescue candidates) |
| **R2** | low-complexity **filter** (drop clips) | **windowed evidence count** (admit more) |
| **R8** | hallmark scoring | **short poly-A clips** (down to 7 bp) |

On the actual gates the directions are opposed: implement_review **R1/R16** add a
low-MAPQ-fraction and depth mask; sensitivity **R3** *lowers* the MAPQ floor 40→20.
implement_review **R16** makes the 2-read floor coverage-adaptive (raises it in pile-ups);
sensitivity **R2** makes the same floor windowed (lowers its effective value).
implement_review **R14/R2** drop recurrent/low-complexity candidates; sensitivity **R5**
emits *more* one-sided candidates.

Anyone told "do R2 first" by both docs implements a filter *and* a relaxation on the same
counting logic. **The two "do-first" lists are not independently shippable; they must be
merged into one signed plan with one namespace before any ordering is meaningful.**

## 2. Which direction wins is architectural — and implement_review has it backwards for a clonal-cohort caller
**Verdict: BOUNDS both plans. Confidence: high on the asymmetry, medium on the net call.**

PEAR-TREE's distinctive backstop is **cross-cohort population genotyping** (require a
wt+insertion partition; review §7.5). Therefore a discovery-stage **false positive is
recoverable** (combine/genotype removes it) but a discovery-stage **false negative is
permanent** (combine cannot resurrect a junction discovery never emitted). This asymmetry
argues the *sensitivity* direction (emit generously, filter late) is better matched to the
tool than implement_review's push to filter hard at discovery (R1 mask, R2, R14 drop).
implement_review's aggressive discovery-time rejection **fights the tool's own
defence-in-depth** by moving irreversible decisions upstream of the population filter that is
PEAR-TREE's actual strength.

**Bounding exception:** truth-preserving, cheap items (R4 stats sidecar, R3
constant-extraction) are safe either way, and the sensitivity direction has its own failure
mode (finding #4). **Reorder both plans around "cheap + reversible first, irreversible drops
last," not "filters first."**

## 3. The sensitivity relaxations compound, and their only "guard" (R1) is defined as never gating
**Verdict: BOUNDS R2/R3/R6/R1 (specificity collapse in the target regime). Confidence: high.**

Priority-1/2 of `recommendations_sensitivity.md` simultaneously: windows the evidence count
(**R2**), drops the MAPQ floor to 20–30 (**R3**), loosens `find_polya` (drop the
`is_proper_pair` reject, purity instead of contiguous-12, floor→7 bp: **R6**), relaxes the
single-read poly-A rescue from exact `AAAAAAAA` to fractional purity (**R1**), and admits
sub-12-bp poly-A clips (**R8**). The stated counterweight everywhere is **R1's hallmark
score — which the plan explicitly says is "never a gate"** (implement_review R8 repeats: "a
*score*, never a gate"). A guard that can reject nothing does not hold specificity.

Concrete failure, in the low-VAF/low-coverage regime the plan *targets*: an event survives on
exactly 2 reads. R2 now counts any 2 clips within ±3–5 bp of the mode as support (verified:
`most_common_first`, `model.rs:148`, counts exact `i64` positions today). Two *unrelated*
single-read artefacts 4 bp apart — abundant given review §10's ~45%-recurrent,
poly-A/Alu-dominated clip population — now manufacture a breakpoint, while R3 has
simultaneously admitted the low-MAPQ repeat-anchored reads that are the richest artefact
source. The relaxations multiply; nothing gates the product. **At least one of R2/R3/R6/R1
must be a real gate (e.g. the R3 low-MAPQ-clip-ratio cap `MAX_LOWQ_CLIP_RATIO=0.65` and the
R16 coverage-adaptive floor) and it must ship in the same commit as the relaxation, not be
deferred to a later tier.**

## 4. sensitivity R5 emits one-sided candidates into a combine step that cannot consume them and a density mask that will wipe them
**Verdict: KILLS R5 as scoped. Confidence: high on the dead-end, medium on the density-wipe.**

R5 ("emit unpaired one-sided breakpoints with a `:LEFT_ONLY:` flag … downstream
combine/annotate can decide") is scoped by its own document to **discovery only**, explicitly
excluding combine. But:

- `combine_insertions.intersect_insertions` has an **unconditional `continue` that already
  drops poly-A-only / one-ended loci** (review §7.7 finding 4 — the very thing R9 exists to
  fix, and R9 is out of scope in the sensitivity plan).
- combine wipes any **100 bp window with ≥4 candidate insertions** (review §7.4.2).

Flooding discovery output with one-sided reference-Alu-boundary clips (every read clipped at a
reference Alu edge is a one-sided breakpoint) means a *genuine* insertion locus — now emitting
L, R, poly-A, plus R2's jittered variants — can itself reach ≥4 candidates in 100 bp and be
**deleted by the density mask.** So R5 does not merely waste work; via the density mask it can
**remove real two-sided events it was meant to help.** **R5 is not a discovery-only change: it
requires R9 plus a density-mask exemption to land first.**

## 5. R15 (SMS reject) ships on-by-default and first, and its safety proof has a hole at the ERV/repeat-anchored class of the correctness anchor
**Verdict: BOUNDS — and this is the single highest real-detection risk (see final section).
Confidence: medium-high.**

The review repeats that SMS "keys on CIGAR shape, not 'is a repeat,' so it is safe for genuine
repeat-consensus clips … including mouse ERV." That defends the *clipped* (insert) side and
**ignores the anchor side.** MEIs insert into and beside repeats; a genuine junction read is
anchored in genomic sequence that is *itself* a divergent young repeat (young L1, IAP flank).
bwa-mem soft-clips ≥8 bp off the **anchor** end when the terminal bases diverge from
reference — producing a read soft-clipped on both ends (real insert clip one side, incidental
repeat-divergence clip the other). R15 as written
(`left_len.min(right_len) >= MAX_CLIP_CLIP_LEN → continue`, targeting the
`right_soft && left_soft` arm at `discovery.rs:221-231`) then **drops a true junction read.**
Under the hard 2-read floor (`MIN_EVIDENCE_READS_PER_BREAKPOINT = 2`, `config.rs:6`), losing
one of two junction reads deletes the whole event. Mouse carries ~10× the palindrome/structure
load (review §11.5), so incidental anchor-end clipping is *more* common exactly where poly-A
cannot corroborate the loss.

Aggravators:
1. R15 is slated **first and on-by-default** (`PEARTREE_KEEP_SMS` *disables* it), so it gates
   the IAPEz anchor by default; the plan asserts anchor survival but never inspects the
   anchor's actual reads to show they clip <8 bp on the far end.
2. The "three callers converge" claim oversells: MEIGA/v-TraFiC reject **any** both-ends clip;
   xTEA only rejects when clips exceed 8 bp (review §12.6). They converge on *direction*, not
   *threshold* — and the exact comparator behind `MAX_CLIP_CLIP_LEN=8` (min? max? either
   side?) is asserted, not shown, so the ported `min(l,r) >= 8` may not be xTEA's actual rule.

**Ship R15 OFF by default, measure it on the real IAPEz anchor reads, gate on the shorter clip
only with the stated tolerance, and add the anchor-side caveat to the "safe for repeats" claim.**

## 6. The validation harness cannot distinguish artefact-removal from truth-removal
**Verdict: KILLS the "every change is measured" specificity guarantee. Confidence: high.**

The specificity gate (review §8.4, echoed in both plans) is: "diff the stats sidecar + final
insertion count on one real WGS BAM before/after each filter; a filter must be shown to remove
*artefacts*, not true calls — check against known germline hets." This is not a specificity
test:

- **No artefact/truth labels** on the real BAM — "count went down" is equally consistent with
  removing truth.
- **`check_germline_coverage.py` probes germline hets (VAF≈0.5)** — the *easiest*, highest-VAF,
  best-corroborated events. It is silent on the low-VAF somatic regime where every relaxation
  (R2/R3/R5/R6) adds FPs and where recall actually matters.
- **One BAM cannot estimate a specificity rate**, and the ~11 toggles are only ever measured
  one-at-a-time; the multiplicative ON-state interaction (finding #3) is never measured
  against truth.
- **MEIsimulator is blind** to EN motif, twin-priming, ERV/LTR; uses **fixed 15 bp TSD** and
  **fixed solo poly-A** (review §12.9) — so it exercises neither the R8 TSD-length score
  (fixed 15) nor the R8 EN score (absent) nor any ERV.
- The proposed real-data fallback for those blind spots names "**the S1 LINE-call
  reanalysis**" (review §12.9). Per the project's own prior finding, the S1 calls
  (`chr3:146970723`, `chr5:40029918`) are **unresolved complex junctions, not confirmed
  insertions** — using them as a truth set for EN-motif / twin-priming validation is circular.

**The byte-identical-when-off baseline is a false comfort:** it proves the OFF state matches
Python, never that any ON state is correct, and the combinatorial ON-state space (11 toggles)
is never validated against truth. A real spike-in truth set covering low VAF, plus orthogonal
(long-read/PCR) confirmation of a sample of *new* calls, are prerequisites, not options.

## 7. EN-motif "190×, would lead the field" is ecological enrichment misused as a per-event score
**Verdict: CAVEAT bordering on kill for the R8/R1 EN component. Confidence: medium-high.**

The 190× figure (review §2.2, §3.3, Nam 2023) is **cohort-level enrichment of insertions at
`TT/AAAA`**, and the review's own cheat-sheet rates the motif "**Moderate — enrichment, not
per-event proof**" (§2.5). R8 (implement) and R1 (sensitivity) nonetheless "reward the EN nick
motif at the 5′ junction" **per event.** The motif is degenerate and ubiquitous in an AT-rich
genome; per-event presence carries a weak prior and **biases scoring toward AT-rich loci — the
same loci that generate the false poly-A pairings review §11.4 warns about.** Compounding: the
motif is written four different ways across the docs (`TT/AAAA`, `TT|AAAA`, `TTAAAA`,
`TTTT|R`), so the implemented pattern is ambiguous, and it applies to L1/Alu only — it must be
gated off for ERV (§11.4), which the "would lead the field" framing never foregrounds. **Keep
EN as a tie-breaker at most, never additive weight comparable to poly-A/TSD; fix a single
motif definition. "Leads the field" is novelty, not a demonstrated gain.**

## 8. R14 within-BAM clip-recurrence cannot honour its own mouse-ERV caveat
**Verdict: BOUNDS (saved only by "flag-first"). Confidence: medium-high.**

R14 flags clip sequences recurring ≥K within a BAM, with Caveat 1 "filter on sequence
recurrence, not 'is a repeat' … doubly essential for mouse ERVs whose clips *are* the IAP/LTR
consensus." But **within one mouse, multiple independent real IAP insertions produce the
identical LTR-consensus clip end** — so their clip sequence recurs ≥K *by the very mechanism
that makes them real*, and R14's within-BAM sequence-recurrence count **cannot distinguish K
real ERV insertions from K reference-Alu mismaps.** The caveat says these "must survive"; the
mechanism deletes them. The cross-*individual* recurrence that *would* discriminate (review
§11.3) is unavailable within one BAM by construction. **R14's within-BAM form is only ever
safe as a sidecar *flag*, never a drop, in any cohort containing LTR/ERV elements;
implement_review's "flag first" is load-bearing and must be permanent for mouse.**

## 9. R2 (cheap low-complexity, ≤2-distinct-bases) is nearly a no-op after its own poly-A escape hatch
**Verdict: CAVEAT. Confidence: medium.**

Verified: the n-polymer filter runs only on `unclipped_cons` (`model.rs:228-235`), so
low-complexity *clips* do pass — R2's premise is correct. But the cheap "≤2 distinct ATGC
bases" test over a 12–150 bp clip is true essentially only for pure homopolymer / pure
dinucleotide runs, and those are exactly what the mandatory **poly-A escape hatch then
exempts.** The two nearly cancel: R2-cheap removes little that is not a poly-A being spared.
The 11.95% low-complexity figure (review §10.1) is dominated by poly-A/T homopolymers (§10.2),
i.e. the exempted class. **If R2 is worth doing it must be the k-mer-entropy version
(`Σ_K|freq−0.25^K| > 0.83`); stop advertising the ≤2-base form as meaningful specificity.**

## 10. Speed plan is largely sound but ships one correctness-risky item and one false "byte-identical" promise
**Verdict: CAVEAT. Confidence: medium.**

`recommendations_speed.md` is the strongest of the three (R1 lazy decode, R2 multithread BGZF,
R5 internals are genuinely byte-identical). Two flags:

- **R4 (mate-coordinate fetch replacing the qname scan)** is asserted byte-identical because
  "the same mate reads are visited," but coordinate-fetch + coalescing changes *iteration
  order*, and `get_mates`/`find_mates` push into `mate_seqs`/`mates` in encounter order; if any
  downstream consumer is order-sensitive (mate consensus, the `MATE{i}` numbering at
  `discovery.rs:504,517`), output diverges. It needs the differential test to *prove*
  order-independence, not assume it.
- **R5b (`u8` quals)** is correctly flagged as touching `QualitySeq`, but consensus scores
  (`filters.rs:117`, `delta_best`) sum qualities and exceed 255, so the split must be airtight
  or consensus truncation shifts.

**R4/R5b are not free byte-identical wins; gate both behind the real-WGS differential, which
the plan itself admits does not yet exist.**

---

## The single change most likely to degrade real-world detection if shipped as written

**R15 (SMS both-ends-clip reject), shipped on-by-default and first, as `implement_review.md`
specifies.** It is the flagship "do this first," it defaults to *active*, and its safety
argument defends the wrong side of the read: the "keys on CIGAR shape, safe for repeat clips"
proof covers the inserted clip but not the **anchor end**, which is soft-clipped ≥8 bp
precisely when a true junction read is anchored in the divergent young-repeat / LTR flank that
MEIs favour. Under the hard 2-read floor, silently dropping one such read deletes the event —
and it does so most often for the mouse ERV/IAP class that includes the IAPEz correctness
anchor, in the genome (mouse) where poly-A cannot corroborate the loss and the palindrome load
is 10× higher. It is the rare filter that is simultaneously (a) prioritized first, (b) on by
default, (c) targeted at the anchor's own element class, and (d) validated only by an
assertion that the anchor "must persist," never by inspecting the anchor's actual reads. Ship
it OFF by default, prove non-removal on the real IAPEz reads, or it will quietly cost ERV and
repeat-anchored L1 recall the moment it meets real WGS.

**Runner-up:** sensitivity **R1**'s relaxation of the single-read poly-A rescue (exact
`AAAAAAAA` → fractional purity) in the mouse configuration — review §11.4's own data (18.48%
poly-A-rescued loci, already suspect as genomic-poly-A false pairings) says this path is
*noisier* in mouse, yet the plan cites §11.4 only for the recall caveat and never notices it is
arguing against R1's own specificity there.
