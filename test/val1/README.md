# VAL-1 — synthetic truth harness

The specificity/recall gate for every behaviour toggle (SPEC-\*/SENS-\*). A
behaviour change is **not defaulted on** until it clears VAL-1: it must not lose
true insertions or add unexplained calls versus baseline.

Because the plan chose *Rust-only, no equivalence reference*, this synthetic
truth set — plus per-toggle Rust unit tests — is the only automated check on the
toggle-**on** path. The Python differential only validates toggle-**off**.

## Pieces

| script | does |
|---|---|
| `simulate.py` | writes a coordinate-sorted BAM of labelled spike-in insertions + a truth TSV. **Every element class and structural variant** (L1Hs / Alu / SVA / ERV-K, plus HERV-K113 & HERV-K117 positive controls), variable TSD, variable VAF (via alt-read count), reference coverage reads, **and the full artefact catalogue** of `RTE_detection_review/05_sequencing_artefacts.md`. See "What simulate.py emits" below. |
| `score.py` | runs discovery under one config, matches calls to truth (±window), reports recall/precision overall, per class, per VAF bin, + the OBS-1 stats sidecar. |
| `run_matrix.py` | runs the scorer over `baseline` + every `configs/*.txt`, tabulating recall/precision. This is the one-at-a-time **and** combined-ON matrix the plan requires. |
| `new_calls.py` | diffs two discovery outputs and writes the calls new in the candidate as a BED — the list for orthogonal (IGV / long-read / PCR) confirmation on real data. |

## Use

```bash
# 1. simulate a labelled set
venv/bin/python test/val1/simulate.py --out-bam sim.bam --out-truth truth.tsv --n-l1 30 --n-erv 20

# 2. score one config (omit --config for baseline)
venv/bin/python test/val1/score.py --bam sim.bam --truth truth.tsv --config my_toggle.txt

# 3. run the whole matrix (baseline + configs/*.txt)
venv/bin/python test/val1/run_matrix.py --bam sim.bam --truth truth.tsv

# 4. when a toggle adds calls on REAL data, list them for orthogonal review
venv/bin/python test/val1/new_calls.py --baseline base.txt.gz --candidate toggle.txt.gz \
    --out-bed new.bed [--truth truth.tsv]
```

Adding a new toggle (later phases): drop a `configs/<name>.txt` enabling it, and
extend `combined.txt`. No harness code changes — `run_matrix.py` picks it up.

## What `simulate.py` emits

The truth `class` column is `<ELEMENT>_<variant>`, so `score.py` reports recall per
element **and** per variant.

**Calibrated to real data.** The defaults are measured from a real GRCm38 (mm10) mouse
WGS CRAM (`45521#13`, PCR-free NovaSeq), decoded against `GCF_000001635.20_GRCm38` via
an MD5-keyed `REF_CACHE`. Two profiling passes:

- *Geometry* (1.8M reads): 151 bp reads; soft-clip lengths drawn from the empirical
  distribution (p05/p50/p95 = 15/59/116 bp) so junction reads span each breakpoint at
  realistic, varying offsets; ~320 bp insert size; read-level artefact prevalences the
  defaults reflect (SMS 0.36 %, cruciform/SA-same-contig<1kb 2.47 %, maps-elsewhere
  3.2 %, duplicates 3.7 %). `--junction-lowmapq-frac` reproduces the ~12 % of clipped
  reads that anchor in a repeat and drop to MAPQ 0 (off by default so the recall gate
  stays deterministic).
- *Clip sequence composition* (168k clips): adapter read-through 5.5 % (the exact
  `AGATCGGAAGAGC…` adapter is a top recurrent clip), poly-G 2.8 %, low-complexity 2.6 %,
  and proximal poly-A/T run length median 7 bp — so TPRT poly-A tails are now sampled
  (13..`--polya-len`, above the 12 bp `POLYA_CUTOFF`) rather than a fixed clean 18.

`--read-len`, `--polya-len`, and `--element-fasta` expose the rest.

**Real insertions (in truth).** Each TPRT element (L1Hs, Alu, SVA) is emitted in
every structural variant — full-length, 5′ truncation, 5′ inversion (twin priming),
partnered / orphan 3′ transduction, 5′ transduction, 3′ (poly-A) deletion, internal
inv/del — with the poly-A tail and an L1 endonuclease `TTAAAA` flank motif. The
**transduction** variants carry a non-RTE **tag** co-mobilised from a unique source
locus: the reads crossing the tag-bearing junction have their **mates at that source
locus on another contig** (the barcode that traces the event to its master element).
This makes the transductions the *true-positive* mirror of the `--n-translocation`
mimic — both join the insertion site to a distant locus, so a transduction/translocation
discriminator must key on the RTE hallmarks (poly-A + TSD + EN motif, present here,
absent in the translocation), not merely on the cross-locus pointer. Processed
pseudogenes (Feature B) are the fully non-RTE case: an mRNA retrotransposed by L1
machinery, flagged by mates spanning ≥2 exons. ERV/LTR elements (generic
ERV-K/IAP-like) are emitted as full provirus, solo-LTR and 5′-inversion, with a short
TSD and **no poly-A**. **HERV-K113 and HERV-K117** are inserted as full-length HML-2
proviruses — positive controls so that if one ever retrotransposed we would see it.
Junction *clip content* is drawn from per-element consensus termini; the built-ins are
structural stand-ins (correct length/GC and the `TG…CA` LTR boundary) — pass
`--element-fasta` (records `<ELEMENT>_end5/_end3/_body/_ltr`) with the real Dfam
consensus to make the clips annotate-detectable.

**Sensitivity artefacts (in truth, MISSED at baseline; opt-in, default 0):**
`--n-wobble` (SENS-1), `--n-lowmapq` (SENS-2), `--n-shortpolya` (SENS-8),
`--n-polya-dropout` (§5.3 — a genuine L1 whose poly-A junction reads are lost to the
library prep, so it cannot be paired), and `--n-discordant-only` (the breakpoints fall in
the **unsequenced insert gap** between read1 and read2, so no read soft-clips anywhere and
the event is supported only by discordant pairs — PEAR-TREE's clip-first discovery has no
discordant-pair discovery, §7.7#5, so it is missed at baseline *and* by the one-sided
Feature-A rescue). Contrast Feature A `--n-discordant`, which models the same gap on **one**
junction but keeps a real soft-clip on the other for the `discordant_anchor` rescue.

**False-positive artefacts (NOT in truth; a call here is a FP; on by default):** all of
`RTE_detection_review` §5 — adapter read-through, poly-G/dark-cycle, homopolymer /
low-complexity, PCR duplicates, PCR chimeras, segdup/mismap (XA full-length), and
cruciform/inverted-repeat (SA same-contig <1 kb). Existing filters remove every one of
these (0 FPs expected). Some remain a **known discovery-time gap** and therefore surface
as baseline false positives until a filter lands:

- `--n-sms` (both-ends-soft-clipped chimeras) and `--n-palindrome` (self-fold
  palindromes) — the structure-specific-chimera gap (review R15);
- `--n-translocation` (§4.2.9) — a **chromosomal translocation between two
  retrotransposon loci**, the most dangerous mimic. Breakpoint microhomology gives a
  LEFT+RIGHT clip pair with a fake TSD, and because the loci are RTEs the clips are LTR
  consensus with no poly-A, so it is identical to an ERV insertion at discovery. Each
  read carries the real translocation signatures — a split alignment (`SA`) to the
  partner **contig** (so the same-contig cruciform check never fires) and a discordant
  mate there — and each translocation emits both reciprocal breakpoints. No discovery
  filter rejects it; only the combine-step end-to-end remap (the clip maps fully to the
  partner locus) or a reciprocal-partner / cross-chromosome-`SA` filter can. This is the
  target for validating any such filter.

`--n-artefacts` adds high-coverage pile-ups for the SPEC-3/4 coverage gates to remove.

**`--n-mismap-ambiguous` — mapping-ambiguity FP (models the full-stack assembly-discordance
calls).** A clipped pair whose anchor reads are *mostly* ambiguously placed (low MAPQ,
`XS`≈`AS`, and a **partial** `XA` alt so `reject_fully_mapping_reads` does not fire) with a
unique minority (`--mismap-ambiguous-uniqfrac`, default **0.85** — calibrated to the real
hs1→hg38 FPs, whose shifted hg38 homolog is locally unique so most reads still look
confident). It clears the ≥2-per-side call floor → a baseline FP. This is the target class
for the **anchor-uniqueness filter (#2)**: a true insertion has an all-unique anchor set,
this does not. Because the realistic unique fraction (0.85) sits just below a true
insertion's (~0.97 on real data, 1.0 in this synthetic), the separating threshold is narrow
— the filter is a genuine precision/recall trade, not a free win. A depth-normalised
*fraction* of unique-or-mate-anchored reads separates it; an absolute count does not (it
conflates anchor quality with read depth and erodes low-VAF recall).

**`--n-materescue` — mate-anchored rescue target (validates discovery's `mate_anchor_rescue`).**
A real L1 insertion into a low-mapability-but-mate-unique flank: the junction clips carry a
below-floor own MAPQ (`--materescue-mapq`, default 20) but each is a proper pair whose mate
maps uniquely (`MQ` tag = 60). **Missed at baseline** (own MAPQ < `min_mapq`); recovered by
`mate_anchor_rescue`, which trusts the unique mate. This is the full-stack `id=32` scenario
(an L1 3′-transduction whose hg38 homolog is repeat-adjacent). Note the interaction with
filter #2: these TPs have *non-unique own anchors*, so the anchor-uniqueness filter must be
**mate-aware** (count a read as anchored if its own anchor OR its mate is unique) or it
re-drops exactly what mate-rescue recovered.

**Scaling to large sets.** `--contig-len` (default 20 Mb) and `--step` (default 100 kb, the
inter-insertion spacing) let the 3-contig reference hold many thousands of insertions — e.g.
`--contig-len 30000000 --step 5000` fits ~10 000 TPs. Measured on such a run (10 252 TPs,
40 translocations, 300 cruciform, 200 mapping-ambiguity, seed 1): baseline recall
10 052/10 252 with **200 misses all in the `materescue` class → 100 % recall with
`mate_anchor_rescue`, 0 added FPs**; the 380 baseline FPs are 200 mapping-ambiguity + 80
translocation (40×2 reciprocal) + 50 SMS + 50 palindrome (cruciform/adapter/poly-G/PCR/
low-complexity/full-`XA` mismap all → 0). The mate-aware anchored *fraction* removes the 200
mapping-ambiguity FPs at 0 TP loss but is structurally blind to the 180 translocation/SMS/
palindrome FPs (unique anchors — those need reciprocal-`SA` / SMS / self-fold filters).

**Feature A — discordant anchoring (`discordant_anchor`):** `--n-discordant` emits a
one-sided junction (a real clip on one flank only) whose missing side is supplied by a
cluster of discordant mates pointing into an RTE band — recalled only when the toggle is
on (in truth). `--n-disc-artefact` is the same geometry but the mates point to random
genome (NOT in truth) — rescued by `discordant_anchor` alone, dropped once
`discordant_rte_only` is set. `--out-rmsk` writes the matching RepeatMasker track for
`discordant_rte_track`.

**Feature B — processed-pseudogene (`splice_hallmark`):** `--n-pseudogene` emits
insertions whose mate reads map into ≥2 exons of a synthetic gene (introns skipped) —
flagged in `<out>.splice.tsv`; `--n-single-exon` is the negative control (mates in one
exon, must not flag). `--out-exons` writes the matching exon annotation for
`exon_annotation`. These are annotation checks (the sidecar), scored by inspecting
`<out>.splice.tsv`, not by the recall/precision table.

Every count is tuneable; set any `--n-*` to 0 to drop that class. `--n-l1` / `--n-erv`
still drive the canonical full-length counts (backward compatible).

## ⚠ Synthetic-only limitation (do not skip)

This proves **necessity, not sufficiency**. The simulator now emits a synthetic
member of every artefact class (§5), so it exercises the *code path* each real
artefact would hit — but the reads are still clean and the junctions exact, so it
cannot reproduce the messiness of real enzymatic-prep chimeras. It still cannot
fully prove specificity against:

- real **EN/ERV/twin-priming** artefacts (the mouse ERV/IAP class especially),
- the **low-VAF** specificity regime,
- realistically noisy palindrome / cruciform / microindel chimeras.

So VAL-1 green is a gate, not a clearance. Before any toggle is defaulted on,
requirement **(b)** still stands: orthogonal (long-read / PCR / IGV) confirmation
of a sample of the *new* calls it produces on real data (`new_calls.py`).
