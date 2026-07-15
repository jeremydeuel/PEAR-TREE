# Empirical scan scripts and outputs

Scripts and raw results behind [§10](../10_empirical_discovery_patterns.md) and
[§11](../11_cross_dataset_and_mouse_erv.md) of the review. Each script parses PEAR-TREE
step‑1 discovery output (custom FASTQ-like format, `@seqname:R-L:SIDE:FIELD`) and measures
artefact-pattern prevalence in the **CLIPPED** fields (the putative inserted sequence),
plus sequence recurrence.

## Layout

```
scans/
├── scripts/   the analysis scripts (plain Python 3, stdlib only)
└── outputs/   the captured stdout of the full-set runs
```

## Scripts

| Script | Dataset it targets (`DIR`/`ROOT` constant at top) | What it computes |
|---|---|---|
| `scan_discovery.py` | Human MPN `trees_mpn/discovery` (seeded 40-file sample) | Prevalence table + top-25 recurrent clipped sequences. First-pass sample. |
| `scan2.py` | Human MPN `trees_mpn/discovery` (seeded 25-file sample) | Recurrence concentration, repeat-consensus fraction, and the concrete **palindrome self-fold** and **poly‑G** example records quoted in §10. |
| `scan_full.py` | Human MPN `trees_mpn/discovery` (**all 560 files**) | Full prevalence table, recurrence, Alu/L1 consensus fraction, top-30 recurrent clips. Produced `outputs/human_mpn_trees_full.out`. |
| `scan_mouse.py` | Mouse `spar_2ndrev/.../discovery` (**all 2,074 files, 32 individuals** by `MD####`/`MX####` prefix) | Prevalence table incl. poly‑A-run content, per-individual recurrence (germline vs artefact), top-40 recurrent clips with #individuals. Produced `outputs/mouse_spar2ndrev_full.out`. |
| `scan_human2.py` | Human `spar_2ndrev/.../discovery` (**all 15 individuals**, one folder each; skips `genotyping`) | Prevalence table, per-individual and cross-individual recurrence, top-35 recurrent clips with #individuals. Produced `outputs/human_spar2ndrev_full.out`. |

## Outputs

| File | Dataset | Files / individuals | Clipped seqs |
|---|---|---:|---:|
| `human_mpn_trees_full.out` | Human MPN trees | 560 / 1 cohort | 5.80 M |
| `human_spar2ndrev_full.out` | Human spar_2ndrev | 3,612 / 15 | 50.72 M |
| `mouse_spar2ndrev_full.out` | Mouse spar_2ndrev | 2,074 / 32 | 8.56 M |

## Running

Each script hard-codes the dataset path in a `DIR` / `ROOT` constant at the top; edit that to
re-point it, then:

```bash
python3 scripts/scan_full.py    > outputs/human_mpn_trees_full.out
python3 scripts/scan_mouse.py   > outputs/mouse_spar2ndrev_full.out
python3 scripts/scan_human2.py  > outputs/human_spar2ndrev_full.out
```

Python 3, standard library only (`gzip`, `glob`, `collections`). No external dependencies.
The full runs take ~1–2 min each and hold a `Counter` of distinct clipped sequences in
memory (a few GB for the 50 M-clip human set).

## Method notes / caveats

- **Prevalence** is per clipped sequence; a clip may fall in several categories (e.g. a
  poly‑A run is both low-complexity and poly‑A/T-dominated), so columns are **not** mutually
  exclusive.
- **Recurrence** counts identical clipped sequences (≥20 bp) and is **scale-dependent** —
  the same "≥5×" threshold yields a higher fraction on a larger file set (see §10.3). Compare
  recurrence only within a dataset, not across datasets of different size.
- **Per-individual recurrence** (mouse/human `spar_2ndrev`) counts, for each distinct clip
  seen ≥2× within an individual, how many individuals it appears in. "All individuals" =
  reference-repeat/artefact; "one individual" = candidate germline. This encodes the rule
  that a variant private to one individual (but shared across its colonies) is germline.
- **Repeat classification** (`Alu`, `L1`) uses a few hard-coded consensus 15-mers and is a
  **lower bound** — divergent subfamilies won't match. Mouse recurrent clips are described by
  composition (GC-rich LTR/repeat consensus, microsatellite) rather than formally annotated;
  a Dfam/RepeatMasker pass would label them precisely.
- Scripts assume the discovery files are the standard PEAR-TREE `.txt.gz` / `.discovery.fq.gz`
  4-line records; malformed files are skipped with an `ERR` line to stderr.
