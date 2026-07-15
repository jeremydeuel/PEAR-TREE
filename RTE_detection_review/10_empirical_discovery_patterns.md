# 10. Empirical artefact patterns in real PEAR-TREE discovery output

*Data: `.../trees_mpn/discovery/` — **560** PEAR-TREE step‑1 discovery files
(`*.discovery.fq.gz`) from the JAK2/HNRNPA1 MPN colony trees (PD4781 series). This section
reports what artefact patterns actually occur in this output, confirms they are covered by
the artefact taxonomy in [§4](04_artefacts_in_discovery.md)–[§5](05_sequencing_artefacts.md),
and flags where the empirics validate a code-review finding from [§7](07_peartree_code_review.md).*

*Method: parsed the custom FASTQ-like discovery format (`@seqname:R-L:SIDE:FIELD`) across
**all 560 files** (**5,798,872** clipped sequences, **2,899,436** loci) and measured the
prevalence of each pattern in the **CLIPPED** fields (the putative inserted sequence, where
artefacts surface). Full scan ran in ~66 s. Scripts: `scan_full.py` (full set),
`scan_discovery.py` / `scan2.py` (seeded-sample versions used for the concrete example
records in §10.3–§10.5). Figures below are the **full-set** values.*

---

## 10.1 Headline numbers (all 560 files)

| Pattern | Prevalence among clipped seqs | Interpretation |
|---|---:|---|
| **Homopolymer-dominated** clip (longest run ≥50% of clip) | **13.71%** | Homopolymer/microsatellite |
| **Poly‑A/T-dominated** clip (≥12 bp homopolymer, >70% A/T) | **12.61%** | Mix of true poly‑A tails **and** mononucleotide-run artefacts |
| **Low-complexity** clip (≤2 distinct bases) | **11.95%** | Overlaps the above |
| **Recurrent identical** clip (same seq ≥5× across the cohort) | **45.27%** of all clips (52.65% of ≥20 bp clips) — 190,259 distinct seqs, 2,625,045 instances | Non-unique clips: shared homopolymers + repeat-element consensus. **Scale-dependent — see §10.3** |
| ↳ of which **Alu**-consensus | **183,824** instances | Reference/polymorphic Alu + mismap at Alu |
| ↳ of which **L1**-consensus | **7,362** instances | ~25× rarer than Alu; genuine L1 clips are locus-specific (good) |
| **Short** clip (<15 bp) | **3.86%** | Near the `min_clip_len` boundary |
| **Poly‑G run ≥9** anywhere | **0.66%** (3′ tail **0.28%**) | Two-colour dark-cycle artefact |
| **Inverted-repeat / palindrome self-fold** (clip = RC of adjacent flank) | **0.65%** (37,626 clips) | Cruciform / structure-specific chimera |
| **Adapter** substring present | **0.04%** | Already well handled by `is_adapter`/`clean_clipped_seq` |
| **polyA-rescued loci** | **2.96%** of loci | The poly‑A pairing path is active but a minority |
| clip GC content (mean) | **0.40** | — |

The dominant real-world patterns are, in order: **(1) homopolymer/low-complexity clips
(~12–14%)**, **(2) recurrent repeat-element (overwhelmingly Alu) clips**, and a long tail of
**(3) inverted-repeat/palindrome and poly‑G structure artefacts (~1%)**. All three are
already in the taxonomy of [§5](05_sequencing_artefacts.md); the empirics tell us their
**relative weight in this data**. The per-clip prevalences are essentially identical to a
25–40-file sample (within ~0.2%), confirming they are stable properties of the data, not
sampling noise. The **recurrence** figure is the exception — it rises with cohort size
(§10.3).

## 10.2 Pattern 1 — homopolymer / low-complexity clips (the dominant pattern)

The most common recurrent clipped sequences are pure poly‑A and poly‑T runs of every length
(the top ~20 recurrent sequences are all `AAAA…`/`TTTT…`). Example real record:

```
@chr10:178454-polyA_178521:LEFT:CLIPPED
TTTTTTTAAAAAAAAAAAAAAGGGGGGGGG      ← poly-T → poly-A → poly-G in one clip
@chr10:178454-polyA_178521:LEFT:ALIGNED
GGGGGGGGGGGCGGGGGGGGGGCACAGTGGCTTGTGCCTATAATCCCAGCACTTTGGG…  ← poly-G-contaminated anchor
```

**This class is dual-natured** and must not be filtered blindly:
- A poly‑A/T clip at the **3′ end of a real insertion** is the **poly‑A tail hallmark**
  ([§2](02_biology_of_retrotransposition.md)) — genuine signal, correctly routed to the
  `PolyABreakpoint` logic.
- A poly‑A/T clip elsewhere, or a mixed poly‑T/poly‑A/poly‑G clip like the example above, is
  **mononucleotide-slippage / dark-cycle artefact**.

**Code-review link:** that ~12% of *output* clips are low-complexity **empirically confirms
[§7](07_peartree_code_review.md) finding 7 / recommendation R2** — `is_low_complexity` and
`has_well_defined_breakpoint` are implemented but **not wired into the discovery path**, and
the n‑polymer filter only inspects the *unclipped* anchor, so low-complexity *clipped*
sequences pass through. Wiring those filters in (while exempting genuine poly‑A tails routed
through the poly‑A path) would remove a large fraction of these before they reach `combine`.

## 10.3 Pattern 2 — recurrent Alu-consensus clips (the "non-RTE / germline" pattern)

Across the full 560-file cohort, **45.3% of all clipped sequences are non-unique** — the
exact same sequence recurs ≥5× (190,259 distinct sequences, 2.6M instances). **This figure
is deliberately scale-dependent:** in a 25-file subset the same threshold captured only
~8.9%, because a sequence must reach 5 copies to count and there are far more copies to find
across 560 files. The right reading is not a fixed percentage but the *shape*: **a large and
growing fraction of discovery clips are shared across many unrelated loci and samples**, and
that fraction is dominated by (a) poly‑A/T homopolymers (§10.2) and (b) repeat-element
consensus. Among the **non-homopolymer** recurrent clips the sequence is **overwhelmingly
Alu** (183,824 instances) versus **L1** (7,362) — a ~25:1 ratio. The most common
non-homopolymer clips (each ~3,000–4,700×) are the canonical AluY arms:

```
GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGA   ← AluY left arm (5')
CCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGC   ← Alu right arm (3')
```

These are reads clipped at the boundary of **reference/polymorphic Alu elements**, or
**mismapping between Alu copies** ([§5.9](05_sequencing_artefacts.md)). They recur because
Alu is ~10% of the genome, so the same consensus fragment appears at the clip boundary of
many independent, unrelated loci. Two consequences:

1. **Recurrence itself is a usable artefact/germline signal.** A clipped sequence identical
   across many loci and many samples is either a mobile-element consensus (reference/
   polymorphic, i.e. **not a novel somatic event**) or a recurrent library artefact. This is
   the same logic the artefact papers use for enzymatic palindromes and the reference
   pipelines use for MEI-polymorphism-database subtraction ([§3.4](03_detecting_true_events.md)).
2. **PEAR-TREE already removes most of these downstream, not at discovery.** The
   `combine_insertions` bowtie2 **end-to-end** remap deletes candidates whose end maps
   entirely to the reference (a reference Alu will), and the cross-file consensus + density
   mask remove others. So these are expected discovery-stage candidates that are filtered in
   step 2 — but they make up a large share of discovery output and cost combine-time work.
   That the L1 consensus is **~25× less** recurrent than Alu (7,362 vs 183,824 instances) is
   reassuring: genuine L1-insertion clips are rare and locus-specific, exactly as they should
   be.

**Note for a somatic study:** the goal is *novel* insertions, so recurrent reference-Alu
clips are noise here; but Alu is also an active *trans*-mobilised element
([§2.1](02_biology_of_retrotransposition.md)), so a *non-recurrent, locus-unique* Alu clip
with a poly‑A tail and TSD is a candidate somatic Alu insertion and must **not** be filtered
by "is Alu" alone — only by "maps end-to-end to a *reference* Alu locus."

## 10.4 Pattern 3 — inverted-repeat / palindrome / microsatellite chimeras

**0.65%** of clips (37,626 sequences) are **self-folding**: the clip's leading bases are the
reverse complement of the adjacent aligned flank — the structure-specific chimera of
[§5.1](05_sequencing_artefacts.md). Confirmed real examples:

```
chr1:144642607 RIGHT  clip  TATATATACATACACA
                      flankRC …TATATATATATATACACACAG     ← (TA)n/(CA)n microsatellite inverted repeat
chr11:37591676 RIGHT  clip  CTTGTTCAATATAAAATATATTGAACAATATAAA   ← internal palindrome TTGTTCAA…TTGAACAA
chr1:52982163  LEFT   clip  TGTGAGTGTACCGAGGCAAGTCGTGCGAGTGTACAC ← folds back on flank
chr1:65727956  LEFT   clip  GAATTAAAACCACAAAGTGATATTGCATCAAAC    ← clip[:14]=RC(flank)
```

Two sub-flavours are visible: **(a)** dinucleotide-microsatellite inverted repeats
(`(TA)n`, `(CA)n`) — abundant, low-complexity, and doubly caught by a complexity filter; and
**(b)** true unique-sequence fold-backs where the clip reverse-complements nearby genomic
sequence — the classic cruciform artefact that a complexity filter would **miss** and that
needs the explicit inverted-remap / self-palindrome test of **recommendation R6**
([§8](08_synthesis_and_recommendations.md)). At **0.65%** of clips (≈37.6k across the cohort)
this is a real, non-negligible contributor and the empirical justification for adding R6.
(This is a lower bound — it only counts clips whose *first 14 bp* match the RC of the
same-side flank; longer-range and offset fold-backs are not counted.)

## 10.5 Pattern 4 — poly‑G / poly‑C two-colour artefacts

Confirmed real examples, frequently **combined** with poly‑A/T or poly‑C:

```
chr1:740098    LEFT  GGGCGGGGGGGGG
chr1:2529649   LEFT  GCCCCCCCCCCCCCGGGGGGGGG          ← poly-C + poly-G, also self-folds
chr11:33906865 RIGHT TTTTTTGGGGGGGGGAAGAGGGGGGGGGGG   ← poly-T + poly-G
chr10:133577812 LEFT GGCTTGGGGGCCGGGGGCCCCGTGCCCTTGGGGGGCCCCCCCCCCCCCCC
```

`clean_clipped_seq` trims a *3′* poly‑G tail, which handles the cleanest cases, but mixed
poly‑G/poly‑C/poly‑T low-complexity clips like these are better removed by the complexity
filter of R2. Prevalence is modest (poly‑G run ≥9 = **0.66%**; 3′ poly‑G tail = **0.28%**)
but they are unambiguous artefacts.

## 10.6 What the empirics say about the documentation and the code

- **Coverage check:** every empirical pattern maps to an artefact class already documented
  in [§4](04_artefacts_in_discovery.md)/[§5](05_sequencing_artefacts.md) — homopolymer/STR
  slippage (§5.6), mapping/repeat artefacts (§5.9), structure-specific chimeras (§5.1),
  poly‑G (§5.4), adapters (§5.5). **No new artefact class was found that the review was
  missing**; the empirics refine the *weighting* (homopolymer and recurrent-Alu clips
  dominate this data) and supply concrete real examples.
- **Two recommendations are now empirically justified, not just literature-derived:**
  - **R2 (wire in low-complexity filters):** **11.95%** of discovery-*output* clips are
    low-complexity (≤2 bases) — direct evidence the filter is not running on the clipped path.
  - **R6 (self-palindrome / inverted-remap filter):** **0.65%** of clips are genuine fold-back
    chimeras, including unique-sequence ones a complexity filter cannot catch.
- **One reassurance:** recurrent **L1** clips are **~25× rarer** than recurrent **Alu**
  clips (7,362 vs 183,824 instances), and genuine insertion clips are locus-unique — so the
  dominant recurrent noise is reference-Alu (removed downstream by the end-to-end remap), not
  a failure to distinguish real L1 events.
- **One refinement to add to filtering logic** (fed into [§8](08_synthesis_and_recommendations.md)):
  treat **recurrence across loci/samples** as an explicit discovery-or-combine signal (a
  clip identical at many unrelated loci is reference-repeat or artefact). With **~45% of all
  clips non-unique** across the cohort, doing this at discovery would remove a large fraction
  of the load earlier and cheaper than the current per-candidate bowtie2 pass.

## 10.7 Reproducing this scan

Scripts and captured outputs are bundled in [`scans/`](scans/README.md):

```bash
python3 scans/scripts/scan_full.py       # FULL 560-file scan (~66 s): prevalence, recurrence,
                                         # repeat-consensus fraction, top-30 recurrent clips
python3 scans/scripts/scan_discovery.py  # seeded 40-file sample: prevalence + top recurrent
python3 scans/scripts/scan2.py           # seeded 25-file sample: palindrome + poly-G examples
# captured result: scans/outputs/human_mpn_trees_full.out
```
Point `DIR` at the discovery folder. `scan_full.py` parses all 560 files (5.8M clips) in
about a minute; the sampled scripts are kept because they produced the concrete example
records quoted above. Per-clip prevalences agree to within ~0.2% between sample and full
set; only the recurrence figure is scale-dependent (§10.3).
