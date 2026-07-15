# 4. Rejecting artefacts and non-retrotransposition events at discovery (focus question 2)

*This section is the practical bridge between the artefact catalogue ([§5](05_sequencing_artefacts.md))
and the detection logic ([§3](03_detecting_true_events.md)): for each artefact/non-RTE
class, what read-level test rejects it **as early as possible**, ideally without needing the
reference genome. It then maps those tests onto PEAR-TREE's actual filters.*

---

## 4.1 Principle: filter in cheapening order, positive-evidence first

The reference pipelines and the artefact papers converge on a **defence-in-depth** ordering.
Apply the cheapest, genome-free tests first so the expensive genome-aware tests run on far
fewer candidates:

1. **Read-intrinsic sequence tests** (no genome): adapter, poly‑G, homopolymer/low-complexity,
   self-palindrome.
2. **Cluster-consistency tests** (no genome): ≥ N independent reads, breakpoint convergence,
   consensus agreement, TSD geometry, orientation sanity.
3. **Local re-mapping tests** (needs genome): does the clipped part map back nearby? does an
   end map end-to-end?
4. **Cross-sample / population tests**: recurrence pattern, matched-normal, blacklist,
   panel-of-normals, VAF/clonality.

Equally important is to **require positive RTE evidence** (poly‑A, TSD, element match) rather
than only subtracting negatives — a candidate with a clean poly‑A + TSD is worth keeping even
against weak coverage, while one with neither hallmark near a repeat is worth distrusting.

## 4.2 Artefact-by-artefact discovery-time tests

### 4.2.1 Structure-specific chimeras (inverted repeats / palindromes) — the top priority
*The most dangerous mimic for a clipped-read caller ([§5.1](05_sequencing_artefacts.md)).*
Discovery-time tests:
- **Both-ends-soft-clipped (SMS) read reject — the cheapest, most convergent test.** A read
  whose CIGAR is soft/hard-clipped at *both* ends (`^S…M…S$`) is the read-level fingerprint of a
  reference palindrome/inverted repeat or adapter-dimer. **MEIGA, v‑TraFiC and xTEA all reject it**
  ([§12.6](12_tool_implementations_compared.md)); it is genome-free, one line, and keys on CIGAR
  shape (so it does not touch genuine repeat-consensus clips). PEAR-TREE does not yet do this
  (recommendation R15).
- **Self-palindrome / local-inverted-remap check.** If the clipped segment is the
  **reverse complement of the adjacent mapped sequence**, or re-maps within a short window
  in **inverted** orientation, reject. This can be done *without* a genome for the
  self-fold case (compare clip vs revcomp of the anchor), and *with* a genome for the
  distal case. **This is the single most valuable discovery-time filter for enzymatic-prep
  data.**
- **Split-read same-locus rejection.** A read with an `SA` (supplementary) alignment to the
  **same contig within ~1 kb** is a local rearrangement/foldback, not an insertion.
- **Palindrome-centre / positional-bias check.** Clip ratio ~50% and junction ≤ 30 bp from
  the read edge are red flags; a variant sitting at the exact centre of an odd-length
  palindrome is almost certainly a PDSM chimera.
- **Reference blacklist.** Pre-compute an inverted-repeat/palindrome BED and drop candidates
  inside it. Filtering raised sonication/enzymatic concordance from **7.8% → 80.4%** with no
  true-variant loss. *(From the ArtifactsFinder code, [§12.5](12_tool_implementations_compared.md):
  the real inverted-repeat params are **arm-pair ≥8 bp, spacer ≥5 bp, sub-arm ≥2 bp, ±50 bp**; the
  tool emits per-base artefact **positions**, not intervals, and its palindrome length gate ships
  disabled — recipe to build a genome-wide BED in [§12.8](12_tool_implementations_compared.md).)*
  Given the SMS test above, treat the blacklist as a **second line**.

### 4.2.2 Adapter read-through
Match the clipped tail against the **kit's adapter sequences** (and 1-bp-error variants) and
clip/reject. Genome-free. Must be configured per library kit/platform.

### 4.2.3 Poly‑G / dark-cycle
Trim high-quality 3′ **poly‑G** before using a clip; a clip that is mostly G is a
two-colour-chemistry artefact, not sequence. Genome-free.

### 4.2.4 Homopolymer / microsatellite / low-complexity
- Reject clips/flanks that are **low complexity** (≤ 2 distinct bases) or a short-period
  (1–4 bp) repeat.
- Cap homopolymer length immediately adjacent to the breakpoint.
- Require a **well-defined breakpoint**: the sequence straddling the junction must not be so
  repetitive that the exact position is ambiguous.
Genome-free. Note the **exception**: a *pure poly‑A/T* clip is not "low-complexity junk" —
it is the poly‑A hallmark and must be routed to the poly‑A logic, not discarded.

### 4.2.5 Random / non-reproducible chimeras (PCR, template switch)
- Require **≥ 2 independent evidence reads** converging on the same breakpoint.
- Require **consensus agreement** among those reads (a real junction has one sequence; a
  random chimera does not).
- Require **reproducibility across independent samples** where applicable.
- Remove **PCR duplicates upstream** (Picard/samblaster) so duplicates don't fake support.

### 4.2.6 High-coverage pile-up regions
Regions with anomalously high depth (centromere/telomere/segdup/rDNA/artefact hotspots)
generate clipped reads in bulk. **Mask by coverage** at discovery time (running per-window
depth; drop candidates in spikes) and/or by a **telomere/centromere/segdup blacklist BED**.
Delly's `-x` excludable-regions + per-interval read cap (L = 1000) is the reference design.

### 4.2.7 Local rearrangements masquerading as insertions
The clipped part of a true insertion is **novel** sequence; the clipped part of a local
inversion/deletion **maps back near the breakpoint**. Re-map the clipped segment
(local mode) and reject if it lands within ~1 kb of its own locus. Needs genome; run after
the cheap filters.

### 4.2.8 Reference structure-specific / mis-assembly / segdup
Re-map the **full aligned+clipped consensus end-to-end**; if either end maps entirely to
the reference, there is **no novel junction** and the candidate cannot be chimeric — reject.
Needs genome.

### 4.2.9 Non-RTE true events (real SVs that are not retrotransposition)
Real inversions, deletions, tandem duplications and translocations also produce clipped
reads. Distinguish them from RTE insertions by the **absence of RTE hallmarks** (no poly‑A,
no TSD, no element-consensus match) and by **breakpoint geometry** (reciprocal partner on
the far side, orientation). A candidate lacking every RTE hallmark is either a generic SV or
an artefact; either way it should not be reported as a retrotransposition. (Caveat from
[§3](03_detecting_true_events.md): some *genuine* reciprocal translocations are
RT-mediated — disambiguate with poly‑A + EN motif + interchanged TSDs at the bridge.)

### 4.2.10 Germline / reference / polymorphic insertions
Not artefacts, but non-*somatic*. Subtract via matched normal, MEI polymorphism databases
(1000 Genomes MEI / dbRIP / euL1db), panel-of-normals, and cross-donor recurrence.

## 4.3 Mapping the tests to PEAR-TREE (what already runs, what doesn't)

| Discovery-time test (§4.2) | PEAR-TREE implementation | Status |
|---|---|---|
| **Both-ends-clipped (SMS) read reject** | keeps the read on its *longer* clip (`discovery.py`) | ❌ **not done — R15**, though MEIGA/v‑TraFiC/xTEA all do ([§12.6](12_tool_implementations_compared.md)) |
| Self-palindrome / local inverted remap | `SA` same-contig-within-1000 bp → cluster poisoning (discovery); clipped `--local` remap + liftover (combine) | ✅ (self-fold via SA); ✅ genome case in combine |
| Reference IVR/palindrome blacklist BED | — | ❌ not implemented (R7; recipe [§12.8](12_tool_implementations_compared.md)) |
| Local-pileup artefact mask (MAPQ/SMS fraction) | — | ❌ not done — R1; v‑TraFiC/MEIGA/xTEA do it reference-free ([§12.6](12_tool_implementations_compared.md)) |
| Adapter read-through | `is_adapter`, `clean_clipped_seq` (NebNext) | ✅ |
| Poly‑G trim | `clean_clipped_seq` 3′ poly‑G | ✅ |
| Low-complexity / short-period repeat | `is_low_complexity`, `has_well_defined_breakpoint` exist; n-polymer unclipped filter runs | ⚠️ partial — the first two are **not wired into discovery** |
| Homopolymer cap at breakpoint | `max_homopolymer_len` (config) | ⚠️ config present; enforcement limited |
| Poly‑A routed to hallmark logic (not discarded) | `PolyABreakpoint`, poly‑A rescue/pairing | ✅ (a strength) |
| ≥ 2 independent reads + consensus agreement | `min_evidence_reads_per_breakpoint`; `find_consensus`; cross-file `sequence_matching_score ≥ 0.6` | ✅ |
| PCR duplicate removal | skips `is_duplicate` (relies on upstream marking) | ✅ (prerequisite) |
| High-coverage pile-up mask | `max_read_count` unreferenced; genotype check `if False`; 100 bp density mask in combine | ❌ discovery; ✅ partial (combine density) |
| Blacklist / weird-contig exclusion | `len(name) > 5`, `MT`/`chrM` skip | ⚠️ crude proxy, not a real BED |
| Clipped-maps-locally (genome) | bowtie2 `--local` + liftover, ≤ 1000 bp → reject | ✅ |
| End-maps-entirely (genome) | bowtie2 `--end-to-end` → reject | ✅ |
| Ref/alt similarity (nothing to genotype) | `sequence_matching_score` ref vs alt | ✅ |
| Germline/somatic classification | population genotyping: require wt + insertion samples; `max_artefact`, `max_na` | ✅ (cohort-based, in place of matched normal) |
| MEI polymorphism DB subtraction | — (annotate step gives element identity, not a polymorphism filter) | ❌ optional add |
| Element-consensus identity (positive) | `annotate` tool: DFAM HMM + RepeatMasker | ✅ post-hoc, not a discovery gate |
| TSD geometry (positive) | `L < R` ordering, `max_bp_window = 40` | ✅ (core design) |
| EN motif / TSD-length prior (positive) | — | ❌ optional add (cheap precision) |

> **Empirically grounded.** [§10](10_empirical_discovery_patterns.md) measures these classes
> across **all 560** real PEAR-TREE discovery files (5.8M clipped sequences): low-complexity
> clips **11.95%**, homopolymer-dominated **13.71%**, inverted-repeat/palindrome self-fold
> **0.65%**, poly‑G **0.66%**, and **~45% of clips non-unique** (recurrent, Alu-dominated)
> across the cohort — confirming which of the tests below matter most for this data and that
> the low-complexity filter (row "Low-complexity") is not yet running on the clipped path.

**Reading of the table.** PEAR-TREE already covers the majority of the recommended
discovery-time tests, and its genome-aware artefact removal (clipped-maps-locally +
end-maps-entirely) is well designed. The concrete gaps are: **(a)** no reference
inverted-repeat/palindrome blacklist — the highest-value missing filter for
enzymatic-prep data; **(b)** high-coverage masking effectively disabled at discovery;
**(c)** two implemented low-complexity filters not wired into the discovery path; **(d)** a
crude contig filter standing in for a proper blacklist; **(e)** no explicit positive score
for poly‑A length / EN motif / TSD length. These become the recommendation list in
[§8](08_synthesis_and_recommendations.md).
