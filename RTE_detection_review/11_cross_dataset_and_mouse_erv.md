# 11. Cross-dataset artefact comparison — multi-individual human and mouse ERV data

*Extends [§10](10_empirical_discovery_patterns.md) to two further PEAR-TREE discovery
datasets from `spar_2ndrev/`: a **15-individual human** set and a **32-individual mouse**
set. The mouse set carries a biology-specific caveat — **mice have active ERVs (IAP and
other LTR retrotransposons), which are inserted by a mechanism that produces NO poly‑A
tail** — so poly‑A cannot be treated as a universal retrotransposition signal there. Both
new sets are **organised so that each folder / filename-prefix is one individual**, which
lets recurrence be measured per-individual and turns "shared across individuals" vs "private
to one individual" into a germline-vs-artefact discriminator.*

*Full scans (not sampled). Scripts: `scan_mouse.py`, `scan_human2.py` (scratchpad).*

---

## 11.1 The datasets now analysed

| Dataset | Ref | Individuals | Files | Loci | Clipped seqs |
|---|---|---:|---:|---:|---:|
| Human MPN trees (`trees_mpn`) | [§10](10_empirical_discovery_patterns.md) | 1 cohort (unsplit) | 560 | 2.90 M | 5.80 M |
| **Human `spar_2ndrev`** | this section | **15** (folders) | 3,612 | 25.36 M | 50.72 M |
| **Mouse `spar_2ndrev`** | this section | **32** (name prefix `MD####`) | 2,074 | 4.28 M | 8.56 M |

## 11.2 Artefact-pattern prevalence: side-by-side

Prevalence among **clipped** sequences (the putative inserted sequence).

| Pattern | Human MPN (§10) | Human spar_2ndrev | **Mouse** spar_2ndrev |
|---|---:|---:|---:|
| Poly‑A/T-dominated clip | 12.61% | 10.00% | **6.99%** |
| Low-complexity (≤2 bases) | 11.95% | 10.00% | 8.82% |
| Homopolymer-dominated (≥50%) | 13.71% | 11.12% | 6.49% |
| Clip contains poly‑A/T run ≥12 | — | — | **9.52%** |
| **Poly‑G run ≥9** | 0.66% | 0.65% | **2.41%** |
| **Inverted-repeat / palindrome self-fold** | 0.65% | 0.27% | **2.64%** |
| Short clip (<15 bp) | 3.86% | 4.15% | 4.12% |
| Adapter substring | 0.04% | 0.02% | 0.03% |
| polyA-rescued loci (of all loci) | 2.96% | 3.16% | **18.48%** |

Two columns of the mouse row stand out and are discussed below: **structure artefacts
(palindrome ~10×, poly‑G ~3.7× the human rate)** and the **poly‑A profile** (lower poly‑A
clip content but a much higher poly‑A-*rescued-locus* rate).

## 11.3 Human, 15 individuals — recurrence separates germline from artefact

With one folder per individual, a clipped sequence can be scored by **how many individuals
it appears in** (counting only clips seen ≥2× within an individual, so single-read noise is
excluded). The result operationalises exactly the rule you stated — *a variant seen across
the colonies of one individual but not in others is a true germline event*:

- **3,692 clips appear in ALL 15 individuals** → reference-repeat / recurrent artefact
  (present everywhere, cannot be a private variant).
- **229,631 clips appear in ≥50% of individuals** → reference-repeat / very common polymorphism.
- **781,044 clips (47.4% of tracked recurrent clips) are private to exactly ONE individual**
  → **candidate germline / somatic events** — recurrent across that individual's colonies,
  absent from every other individual.

The clips present in all 15 individuals are precisely the ones §10 flagged: **poly‑A/T
homopolymers and Alu-consensus arms** (`GGCCGGGCGCGGTGGCTCAC…`, `CCTCGTGATCCGCCCGCC…`), each
seen 17,000–41,000×. This is the key refinement to **recommendation R14**: recurrence must
be scored **per individual, not globally** — a clip shared across many *individuals* is
reference-repeat/artefact and can be dropped, whereas a clip shared across many *colonies of
one individual* is a germline call and must be **kept**. A naïve global-recurrence filter
would delete exactly the germline events the study is looking for.

*(Cohort-wide, 73% of ≥20 bp clips recur ≥5× — higher than the 560-file set purely because
there are more files; see the scale-dependence note in [§10.3](10_empirical_discovery_patterns.md).)*

## 11.4 Mouse — the ERV / poly‑A caveat and what it means

Mouse retrotransposition is not L1-only. The active families split into two mechanistic
classes with **different molecular scars**:

- **LINE-1 (L1Md) and SINEs (B1, B2):** TPRT, **poly‑A tail + TSD** — same signature as
  human L1/Alu ([§2](02_biology_of_retrotransposition.md)). Poly‑A applies.
- **ERV / LTR retrotransposons (IAP, MusD/ETn, MMERVK, etc.):** replicate through an
  integrase-mediated, **LTR-flanked** mechanism (like a retrovirus), **not** TPRT. They
  produce **a short target-site duplication (~4–6 bp; classically 6 bp for IAP) but NO
  poly‑A tail.** Both ends of the insertion are **LTR sequence**, identical to each other.

The empirics are consistent with this being a real, sizeable fraction of the mouse data:

- **Poly‑A clip content is markedly lower than human** (poly‑A/T-dominated 6.99% vs
  10–12.6%; clips with a poly‑A run ≥12 only 9.52%). A large share of true mouse
  retrotransposition (the ERV component) simply **carries no poly‑A to detect**.
- The **top recurrent clips are GC-rich repeat/LTR-like consensus and microsatellites, not
  poly‑A** — e.g. `GCGCCGGCCGAGGCGAGGCGCCGCGCGGAAAACCGCGGCCCGGGGG…`,
  `GTGAGCTCTCGCTGGCCCTTGAAAATCCGGGG…`, `(AGTG)n`, `(TGAG)n` — each present in **all 32
  individuals** (the mouse reference-repeat load, the analogue of human Alu; these are the
  clipped LTR/SINE/LINE consensus ends of reference and polymorphic elements).

**Consequences for detection and filtering (important):**

1. **Never use poly‑A presence as a required positive signal or as a filter in mouse.** A
   poly‑A-gated caller would miss **every ERV/IAP insertion**. PEAR-TREE's core discovery is
   TSD-first / clipped-read-based, so it *can* find ERVs (as two-sided, LTR-clipped,
   short-TSD loci) — but its **poly‑A rescue and poly‑A-pairing paths** (`PolyABreakpoint`,
   the poly‑A rescue in `Breakpoint.join`, the 12–120 bp poly‑A pairing in `output()`) help
   only the L1Md/SINE fraction. Do not let poly‑A become a *requirement* anywhere in the
   mouse configuration.
2. **Recommendation R8 (score poly‑A / TSD / EN-motif as positive evidence) must be made
   species- and element-aware.** For mouse ERVs the positive evidence is **LTR-consensus
   identity + a short (~6 bp) TSD**, and the L1 endonuclease `TT/AAAA` motif does **not**
   apply. A poly‑A-weighted score would wrongly penalise genuine ERV insertions.
3. **"Looks like a mobile-element consensus" must never be a rejection criterion** — even
   more so than in human. A novel IAP insertion's clipped ends *are* the IAP LTR consensus;
   rejecting clips because they match a repeat would delete the very events sought. Reject on
   **locus/individual recurrence and reference end-to-end mapping**, not on sequence identity
   (this is R14 again, and it is essential for ERV calling).
4. **The high poly‑A-*rescued-locus* rate (18.48% vs ~3% human) deserves scrutiny, not
   trust.** Combined with the *lower* poly‑A clip content, it suggests many mouse poly‑A
   pairings arise from the A/T-rich mouse genome's abundant genomic poly‑A / simple-repeat
   tracts rather than from bona fide L1/SINE 3′ ends — i.e. in mouse the poly‑A path is both
   **less useful** (blind to ERVs) **and noisier** (more genomic-poly‑A false pairings). The
   low-complexity/poly-A-tract filters (R2) matter correspondingly more here.

## 11.5 Mouse — elevated structure-specific and poly‑G artefacts

Two artefact classes are several-fold higher in the mouse set than in either human set:

- **Inverted-repeat / palindrome self-fold: 2.64%** — ~**10×** the human `spar_2ndrev`
  rate (0.27%). The top recurrent clips include many self-folding and microsatellite-inverted
  sequences. This is the structure-specific chimera of [§5.1](05_sequencing_artefacts.md).
- **Poly‑G run ≥9: 2.41%** — ~**3.7×** the human rate; frequently combined with poly‑A in
  the same clip (e.g. `GGGGGGGGGGAAAAAAAAAGGG`, `…GTTGGGGGGGAAAAAAAATGG…`).

Whether this reflects the mouse **genome** (different inverted-repeat/simple-repeat content)
or the **library prep** (these mouse samples may be lower-input or more enzymatically
fragmented — [§5.2](05_sequencing_artefacts.md)) cannot be settled from discovery output
alone, but the operational conclusion is the same: **the self-palindrome / inverted-remap
filter (R6) and poly‑G handling are more important for this mouse data than for the human
data**, and should be enabled before the mouse calls are trusted. Because ERV detection can
*not* fall back on poly‑A to corroborate a locus, suppressing these chimeric artefacts up
front matters more here than in human.

## 11.6 Net effect on the recommendations

- **R2 (wire in low-complexity filters):** reaffirmed by both new sets (10% / 8.8%
  low-complexity clips in output); extra weight in mouse, where genomic poly‑A/simple-repeat
  tracts are abundant and cannot be corroborated by an ERV poly‑A tail.
- **R6 (self-palindrome / inverted-remap filter):** **substantially more important for
  mouse** (2.64% vs 0.27%). Prioritise for the mouse pipeline.
- **R8 (positive-hallmark scoring):** **must be species/element-aware.** poly‑A + `TT/AAAA`
  for human L1/Alu and mouse L1Md/SINE; **LTR identity + ~6 bp TSD, no poly‑A** for mouse
  ERV/IAP. Never require poly‑A globally.
- **R14 (recurrence as a filter):** **refined and strengthened** — compute recurrence
  **per individual** using the folder/prefix structure. Shared-across-individuals ⇒
  reference-repeat/artefact (drop); private-to-one-individual-but-shared-across-its-colonies
  ⇒ germline (keep). Filter on **locus/individual recurrence, never on "is a repeat"** —
  mandatory for mouse ERV, where true insertions have repeat-consensus clipped ends.

## 11.7 Reproducing

Scripts and captured outputs are bundled in [`scans/`](scans/README.md):

```bash
python3 scans/scripts/scan_mouse.py    # 2,074 mouse files, 32 individuals (prefix MD####/MX####)
python3 scans/scripts/scan_human2.py   # 15 human individuals (one folder each), 3,612 files
# captured results: scans/outputs/mouse_spar2ndrev_full.out, human_spar2ndrev_full.out
```
Both do full scans, report the prevalence table, global and **per-individual** recurrence,
and the top recurrent clipped sequences with the number of individuals each appears in.
