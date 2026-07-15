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

**Real insertions (in truth).** Each TPRT element (L1Hs, Alu, SVA) is emitted in
every structural variant — full-length, 5′ truncation, 5′ inversion (twin priming),
partnered / orphan 3′ transduction, 3′ (poly-A) deletion, internal inv/del — with the
poly-A tail and an L1 endonuclease `TTAAAA` flank motif. ERV/LTR elements (generic
ERV-K/IAP-like) are emitted as full provirus, solo-LTR and 5′-inversion, with a short
TSD and **no poly-A**. **HERV-K113 and HERV-K117** are inserted as full-length HML-2
proviruses — positive controls so that if one ever retrotransposed we would see it.
Junction *clip content* is drawn from per-element consensus termini; the built-ins are
structural stand-ins (correct length/GC and the `TG…CA` LTR boundary) — pass
`--element-fasta` (records `<ELEMENT>_end5/_end3/_body/_ltr`) with the real Dfam
consensus to make the clips annotate-detectable.

**Sensitivity artefacts (in truth, MISSED at baseline; opt-in, default 0):**
`--n-wobble` (SENS-1), `--n-lowmapq` (SENS-2), `--n-shortpolya` (SENS-8), and
`--n-polya-dropout` (§5.3 — a genuine L1 whose poly-A junction reads are lost to the
library prep, so it cannot be paired).

**False-positive artefacts (NOT in truth; a call here is a FP; on by default):** all of
`RTE_detection_review` §5 — adapter read-through, poly-G/dark-cycle, homopolymer /
low-complexity, PCR duplicates, PCR chimeras, segdup/mismap (XA full-length), and
cruciform/inverted-repeat (SA same-contig <1 kb). Existing filters remove every one of
these (0 FPs expected). Two remain a **known gap (review R15)** and therefore surface as
baseline false positives until a filter lands: `--n-sms` (both-ends-soft-clipped
chimeras) and `--n-palindrome` (self-fold palindromes). `--n-artefacts` adds
high-coverage pile-ups for the SPEC-3/4 coverage gates to remove.

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
