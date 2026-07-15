# 6. Delly: integrated paired-end + split-read SV discovery

*Sources: Rausch et al. (2012) Bioinformatics 28(18):i333–i339; the current
`github.com/dellytools/delly` repository. Delly is included here because it is the
reference implementation of the discordant-pair→split-read paradigm that PEAR-TREE
deliberately does **not** use, and because it is one of the four callers the Nature 2023
colorectal-L1 study ran in parallel ([§3](03_detecting_true_events.md)).*

---

## 6.1 The two-stage algorithm

**Stage 1 — paired-end clustering.** Delly models each sequencing **library** separately,
estimating read-pair orientation and the insert-size distribution (median, SD). It
collects **uniquely mapping discordant pairs** — abnormal orientation, or insert size
beyond the expected range (**default cutoff: 3 SD from the median**). Pairs are binned per
chromosome and organised into a weighted graph *G(V,E)*: nodes are discordant pairs, edges
join pairs consistent with the *same* SV (same orientation, endpoint offsets within the
insert-size range). Edge weight = the difference in implied SV size. Each connected
component ideally is one SV; because data are noisy, Delly extracts a **maximal clique**
per component (seeded on the smallest-weight edge and greedily extended). **Singletons are
discarded.** The clique's extremal coordinates give an approximate SV interval.

**Stage 2 — split-read refinement.** Each PE cluster becomes a candidate breakpoint
interval that Delly screens for **split reads** to reach single-base resolution. It
gathers **single-anchored pairs** (one mate mapped, one unmapped) and optionally
soft-clipped reads, assigning any read within **2 SD** of a breakpoint to that SV. A fast
**k-mer filter (default k = 7)** counts hits along alignment diagonals; diagonals with
< **k_min = 3** hits are dropped, and ≥ 2 supported diagonals are required per read. Delly
then demands **≥ 2 split reads** (default), builds a **majority-vote consensus**, and
re-aligns it to the SV reference with **Gotoh affine-gap dynamic programming** (forward +
reverse matrices, an AGE-like double-DP). This pinpoints both breakpoints, **tolerates
non-templated microinsertions at the junction**, and requires the split-read SV size to
agree with the PE prediction within **10%**.

## 6.2 SV classes and their paired-end signatures

Delly calls **deletions, tandem duplications, inversions, translocations** and (in current
versions) **insertions**. The PE signatures:

| SV class | Discordant-pair signature |
|---|---|
| Deletion | pairs with **larger-than-expected** insert size, normal orientation |
| Tandem duplication | first/second read order **swapped**, orientation retained |
| Inversion | one mate **flipped**; separate left- and right-spanning clusters |
| Translocation | mates on **different chromosomes** (4 sub-types by sort/inversion) |
| Insertion | **smaller-than-expected** insert size, or a single-breakpoint split/clip cluster with no reciprocal partner |

For inversions/duplications/translocations Delly rewrites the reference so a standard
"deletion-type" split-read search resolves the junction.

## 6.3 What Delly does **not** do: dedicated MEI calling

Neither the paper nor the current repository provides a mobile-element/retrotransposon
mode — there is **no `delly mei` subcommand**, no LINE/Alu/SVA classification, no
RepeatMasker/consensus integration. A mobile-element insertion reaches Delly only as a
**generic single-breakpoint insertion signature**: clustered discordant/single-anchored
pairs plus soft-clipped split reads at one locus, with no matching second breakpoint.
Delly can resolve that junction and recover the microinserted sequence, but **the user
must classify the inserted sequence separately** (e.g. against a repeat consensus) and
must supply poly‑A/TSD reasoning externally. This is precisely why the dedicated MEI tools
(TraFiC, MELT, xTea, MEIGA — [§3](03_detecting_true_events.md)) exist, and why the Nature
2023 study used Delly as the **generic-SV arm** alongside three MEI-specific callers rather
than on its own.

## 6.4 Delly's artefact controls (directly relevant to PEAR-TREE)

- **Unique mapping only**; repo FAQ recommends `-q 20` mapping-quality and `-s 15` filters.
- **Support minima:** singletons discarded; ≥ 2 split reads and ≥ 2 good diagonals; the
  quality-filtered PE set in the paper required **avg MAPQ ≥ 20 and ≥ 3 supporting pairs**.
- **Insert-size modelling per library** at 3 SD — adapts the discordant cutoff to each
  library rather than a fixed number.
- **Excludable-regions blacklist (`-x`)**: telomere/centromere and unplaced contigs are
  excluded, because these repeat regions produce massive read pile-ups; Delly also **caps
  split reads per SV interval at L = 1000** to bound pile-up cost.
- **Merging thresholds:** inversion clusters merge at **≥ 80% reciprocal overlap**;
  translocation segments merge within **z = 300 bp**.
- **Germline/somatic separation & genotyping:** `delly filter -f somatic` (matched
  tumour/control via `samples.tsv`), `delly filter -f germline`, multi-sample `delly merge`
  + re-genotyping (`-v sites.bcf`), and read-depth/CNV via `delly cnv`/`classify`.

On simulated data Delly(PE) recovered > 90% of SVs; Delly(SR) traded a little sensitivity
for very high PPV (PCR-validated FDR ≤ 9.1% on 1000 Genomes deletions).

## 6.5 Lessons Delly offers PEAR-TREE

1. **Per-library insert-size modelling** (a data-driven discordant cutoff) is more robust
   than the fixed windows PEAR-TREE uses, and is the natural way to add discordant-pair
   discovery (PEAR-TREE's biggest missing sensitivity source — [§7](07_peartree_code_review.md) finding 5).
2. **An explicit excludable-regions blacklist (`-x`)** plus a **per-interval read cap** is
   the mature version of PEAR-TREE's `len(name)>5` contig filter and disabled
   `max_read_count`. Adopt a real blacklist BED (telomere/centromere/segdup/known-artefact).
3. **k-mer-filtered, DP-based split consensus** with junction-microinsertion tolerance is a
   more principled breakpoint refiner than first-ambiguous-base consensus truncation;
   worth considering if PEAR-TREE wants sub-base-accurate, indel-tolerant junctions.
4. **Delly is complementary, not a competitor:** running a generic SV caller alongside
   PEAR-TREE (as Nature 2023 did) catches events with no junction-spanning clipped read,
   at the cost of needing external MEI classification.
