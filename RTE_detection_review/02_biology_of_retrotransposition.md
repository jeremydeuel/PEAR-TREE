# 2. The biology of retrotransposition — the signatures detection exploits

*Why this section exists: every reliable way to detect a true retrotransposon (RTE)
insertion, and every way to reject an artefact, is grounded in the stereotyped molecular
scar that target-primed reverse transcription (TPRT) leaves in the genome. This section
lays out that scar so that [§3](03_detecting_true_events.md) and
[§4](04_artefacts_in_discovery.md) can refer to it.*

---

## 2.1 The mobile elements that matter in human/mouse WGS

Only a handful of element families are still active and therefore generate *de novo*
(germline or somatic) insertions detectable as variants:

- **LINE-1 (L1):** the only autonomous human retrotransposon. It encodes ORF1p (RNA
  chaperone) and ORF2p (endonuclease + reverse transcriptase). A full-length L1Hs is
  ~6 kb. The overwhelming majority of somatic retrotransposition, and essentially all
  transduction activity, is L1-driven.
- **Alu (SINE):** ~300 bp, non-autonomous, mobilised in *trans* by the L1 machinery. Most
  frequent germline MEI polymorphism.
- **SVA:** ~1–2 kb composite element, also *trans*-mobilised.
- **Processed pseudogenes:** mRNAs reverse-transcribed and inserted by L1 machinery; share
  the same poly‑A/TSD hallmarks but carry spliced exonic sequence.
- In **mouse**, the active families split into **two mechanistic classes** (this matters for
  detection — see below and [§11.4](11_cross_dataset_and_mouse_erv.md)):
  - **TPRT / poly‑A elements:** **L1Md** (LINE) and **B1 / B2 SINEs** — same poly‑A + TSD
    scar as human L1/Alu.
  - **LTR / ERV elements:** **IAP** (note PEAR-TREE's test data is an IAPEz insertion),
    **MusD/ETn**, **MMERVK**, etc. — see §2.1b.

### 2.1b Mouse ERV/LTR elements — a different scar (no poly‑A)
Endogenous retroviruses and other LTR retrotransposons integrate by a **retroviral,
integrase-mediated** mechanism rather than TPRT. Their scar is therefore different:
- **A short target-site duplication** — typically **~4–6 bp** (classically 6 bp for IAP) —
  so TSD geometry still applies.
- **Long terminal repeats (LTRs)** at both ends; the 5′ and 3′ LTRs are identical, so the
  left- and right-clipped ends of the insertion are both **LTR consensus sequence**.
- **No poly‑A tail.** This is the key point: **poly‑A is not a universal retrotransposition
  signal, and requiring it would miss every ERV/IAP event.** Detection of mouse ERVs relies
  on **LTR-consensus identity + the short TSD**, not poly‑A. The mouse dataset in
  [§11](11_cross_dataset_and_mouse_erv.md) shows exactly this — lower poly‑A clip content and
  recurrent LTR/GC-rich consensus clips instead of poly‑A tails.

## 2.2 TPRT and the hallmarks it leaves

L1 ORF2p endonuclease nicks the genomic bottom strand at a degenerate consensus, then uses
the exposed 3′-OH to prime reverse transcription of the L1 (or Alu/SVA/mRNA) poly‑A RNA.
Second-strand synthesis and repair complete the insertion. This mechanism deposits five
signatures, in decreasing order of usefulness for short-read detection:

1. **Poly‑A tail (3′).** The single most specific sequence signature. A homopolymeric
   adenine tract (on the `+` strand; poly‑T on the reverse read) sits at the 3′ junction.
   Detection tools require it to be reasonably long and pure — MEIGA uses **≥15 bp,
   ≥90% purity, within 30 bp of *either* end of the insertion**; PEAR-TREE uses a ≥12 bp A/T run. Its loss to
   library prep ("poly‑A dropout") is a major failure mode (see [§5](05_sequencing_artefacts.md)).
2. **Target-site duplication (TSD).** The staggered EN nick means the insertion is flanked
   by a short **direct repeat**, typically **~2–20 bp** (commonly 10–20 bp). Detecting the
   *same* short sequence on both flanks, in direct orientation, is the geometric proof of a
   TPRT event. In breakpoint coordinates, the two junctions of a TSD‑flanked insertion are
   ordered so that the 3′-side mapped position precedes the 5′-side clipped position — this
   is exactly the `L < R` ordering PEAR-TREE keys on ([§7](07_peartree_code_review.md)).
3. **Endonuclease cleavage motif (5′‑TT/AAAA‑3′).** L1 EN prefers a degenerate
   `5'-TTAAAA` nick site. Enrichment of insertions at this motif is a strong orthogonal
   sanity check — the Nature 2023 colorectal study saw **190‑fold enrichment** at
   `TTTT|R` motifs among true somatic events.
4. **5′ truncation and 5′ inversion.** Reverse transcription is frequently incomplete, so
   most somatic L1s are **5′‑truncated** (only ~5% are full-length). "Twin priming"
   produces **5′ inversions** (the 5′ portion of the insert is inverted). In the Nature
   2023 data, **29.5%** of events carried an intra‑RT (twin‑priming) inversion. The
   practical consequence for short reads: the 5′ junction is often in repetitive L1
   sequence and hard to place, while the 3′ (poly‑A) junction is cleaner — so callers lean
   on the 3′ end.
5. **3′ transductions.** L1 has a weak polyadenylation signal, so transcription often reads
   **through into unique downstream genomic DNA**, which is co‑mobilised:
   - **Partnered transduction:** L1 body **+** downstream unique sequence **+** poly‑A.
   - **Orphan transduction:** the downstream unique sequence **alone** (the L1 body is
     lost to 5′ truncation), still ending in poly‑A.
   Because the transduced sequence is locus‑specific, it **barcodes the source ("master")
   L1 element**. Tubio et al. (2014) showed **95% of transductions trace to just 72
   germline source loci**, two "hot‑L1s" (22q12, 6p24.1) driving over a third of them.
   Transductions are the strongest positive evidence a caller can have — and a notorious
   mimic of translocations (see below).

## 2.3 Structural variety of real insertions (what a caller must tolerate)

Long-read studies (MEIGA/Zumalave 2026) resolved the internal architecture of true events
and show detection must not assume a clean, full-length insert:

- ~51% canonical TPRT with 5′ deletion (median breakpoint ~5.5 kb into L1Hs);
- ~45% carry a 5′ inversion;
- plus internal deletions, duplications, templated insertions;
- **microtransductions < 50 bp** (median 24 bp) make up ~28% of transductions and are
  systematically **misclassified as solo‑L1 by short reads** — a known short‑read blind
  spot.

Insert **lengths** span < 100 bp (5′‑truncated) to full ~6 kb L1; transduced tails are
usually < 1 kb but reach ~12 kb.

## 2.4 Somatic vs germline vs polymorphic — a definitional point

- **Reference/known insertions** are already in the assembly or in MEI polymorphism
  databases (1000 Genomes MEI, dbRIP, euL1db). A read pattern that matches these is not a
  novel event.
- **Germline (polymorphic) insertions** are present in the individual's inherited genome
  (VAF ~0.5 or ~1.0 in *every* clone, and present in matched normal/blood).
- **Somatic insertions** are acquired post‑zygotically: present in some clones/tissues and
  absent from matched normal. In clonally expanded material they sit at **VAF ≈ 0.5**
  (heterozygous in the founder), which is itself the key evidence of a genuine in‑vivo
  event ([§3](03_detecting_true_events.md), Nature 2023).

The detection problem is therefore two-layered: (1) recognise the TPRT scar at all, and
(2) classify the event as reference/germline/somatic using matched samples and population
panels. PEAR-TREE addresses layer 2 with **population genotyping across a phylogenetic
tree of clones** rather than a single matched normal — see [§7](07_peartree_code_review.md).

## 2.5 Signature → detectability cheat-sheet

| Hallmark | Short-read signal | Specificity for true RTE |
|---|---|---|
| Poly‑A tail | clipped read of A/T homopolymer + adjacent unique/L1 seq | **Very high** (if not dropped by prep) |
| TSD | same short direct repeat on both flanks; `L<R` breakpoint order | **Very high** |
| EN motif TT/AAAA | reference sequence at 5′ nick | Moderate (enrichment, not per-event proof) |
| L1/Alu/SVA body | clipped/mate sequence matching repeat consensus | High for "is it a MEI", low for locus |
| 3′ transduction | second poly‑A + unique downstream seq matching a source L1 | **Very high + gives source** |
| 5′ inversion | split read with inverted L1 segment | Confirmatory |
| Clonal VAF ≈ 0.5 | balanced alt/ref spanning reads across clones | High (somatic, in-vivo) |
