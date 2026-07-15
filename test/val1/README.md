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
| `simulate.py` | writes a coordinate-sorted BAM of labelled spike-in insertions + a truth TSV. L1 and ERV classes, variable TSD (2–40), variable VAF (via alt-read count), reference coverage reads. |
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

## ⚠ Synthetic-only limitation (do not skip)

This proves **necessity, not sufficiency**. The simulator makes clean clips,
exact junctions, and no chimeric / PCR / mapping artefacts, so it cannot prove
specificity against the real artefact classes the relaxations risk:

- real **EN/ERV/twin-priming** artefacts (the mouse ERV/IAP class especially),
- the **low-VAF** specificity regime,
- palindrome / cruciform / microindel chimeras.

So VAL-1 green is a gate, not a clearance. Before any toggle is defaulted on,
requirement **(b)** still stands: orthogonal (long-read / PCR / IGV) confirmation
of a sample of the *new* calls it produces on real data (`new_calls.py`).
