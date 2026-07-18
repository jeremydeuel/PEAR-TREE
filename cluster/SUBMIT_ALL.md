# `submit_all.sh` — one-command pipeline submission

Submit the entire 4-stage PEAR-TREE pipeline for **one run** to farm22 (LSF), correctly
sequenced with job dependencies, in a single command. The script is **resumable and
idempotent**: it inspects what already exists and submits only the stages that are missing, so
re-running after a mid-flight failure never re-does completed work and never duplicates an
in-flight array.

```
TAG=<label> DISCOVER_CFG=<cfg> bash cluster/submit_all.sh
```

Quick help any time: `bash cluster/submit_all.sh --help` (also `-h` or `HELP=1`).

---

## The four stages

| # | stage | wraps | output | LSF job name |
|---|-------|-------|--------|--------------|
| 1 | discovery | `submit_discovery.sh` (array, 1 task/colony) | `$DISCDIR/<id>.txt.gz` | `ptdisc_<DISCDIR>` |
| 2 | combine_insertions | `combine_mei.sh` (bsub'd by this script) | `$CONTRACT` (pooled contract) | `ptcomb_<INS_OUTDIR>` |
| 3 | genotype | `submit_genotype.sh` (array, 1 task/colony) | `$GTDIR/<id>.txt.gz` | `ptgt_<GTDIR>` |
| 4 | combine_genotypes | `submit_combine_genotypes.sh` | `$FINAL_OUT` (final call table) | `ptcombgt_<GTDIR>` |

Each stage feeds the next. Stage *N+1* is submitted **held** (`bsub -w "ended(<id>)"`) behind
stage *N*'s numeric job id, so the whole chain can be launched at once and LSF releases each
stage as its predecessor finishes.

> **Why the numeric job id, not the name.** LSF's `ended(<name>)` matches *all* your jobs of
> that name — so a same-named job that already ended satisfies the dependency instantly and the
> dependent fires early against inputs that don't exist yet (a bug that has bitten this repo;
> see `submit_combine_genotypes.sh`'s header). `submit_all.sh` captures each `bsub`'s numeric id
> and chains `ended(<id>)`, which is unambiguous.

---

## Required variables (the script refuses to guess these)

| var | meaning |
|-----|---------|
| `TAG` | run label, e.g. `noSPEC8`. Names every output dir: `discovery_grch38_$TAG`, `insertions_grch38_$TAG`, `genotypes_grch38_$TAG`, and `$STEM.$TAG.genotypes.csv.gz`. |
| `DISCOVER_CFG` | discovery config, e.g. `cluster/config.discovery.grch38.noSPEC8`. **No default** — the wrong assembly's config silently produces empty output. |

If either is unset the script stops immediately with a specific, friendly message and submits
nothing.

## Inferred variables (printed at startup; override by exporting)

| var | default | notes |
|-----|---------|-------|
| `FOFN` | `/nfs/users/nfs_j/jd43/catalogue/picked/all.bams.fofn` | the cohort colony list (confirmed for the 9×10 cohort from job 954062's accounting). |
| `FOFNDIR` | `dirname $FOFN` | holds `all.bams.fofn` + one `<donor>.bams.fofn` per donor. |
| `PATIENTS` | *computed* from the `<donor>.bams.fofn` files in `FOFNDIR` | never guessed; used for combine's per-donor assembly gate. |
| `GENO_CFG` | `cluster/config.genotype.grch38` | genotype config (assembly-specific, not run-specific — same across A/B arms). |
| `STEM` | `mei9x10` | contract basename. |
| `VENV` | `.../PEAR-TREE/venv` | python venv for the combine steps. |
| `DISCDIR` | `discovery_grch38_$TAG` | override to point at an existing discovery dir. |
| `INS_OUTDIR` | `insertions_grch38_$TAG` | contract dir; `CONTRACT = $INS_OUTDIR/$STEM.genotyping.txt.gz`. |
| `GTDIR` | `genotypes_grch38_$TAG` | genotype output dir. |
| `FINAL_OUT` | `$STEM.$TAG.genotypes.csv.gz` | final call table. |

## Tunables

`THROTTLE`(50) `QUEUE`(normal) `GROUP` `ALLOW_MISSING`(1)
`DISC_MEM`(8000) `GT_MEM`(2000) `COMBINE_MEM`(32000) `COMBGT_MEM`(16000)
`COMBINE_CORES`(16) `COMBGT_CORES`(8) `DISCOVER_BIN` `GENOTYPE_BIN`

Memory numbers are MB per task/job. `DISC_MEM=8000` reflects the post-jemalloc measurement
(job 954062: peak 0.4–3.1 GB across 90 BAMs); see `submit_discovery.sh`'s header.

## Modes

| var | effect |
|-----|--------|
| `DRYRUN=1` | run the preflight, print the resolved plan + dependency chain, and **submit nothing**. Recommended before a real run. |
| `SKIP_ASSEMBLY_CHECK=1` | skip the `install.sh check-config` genome_2bit/chain validation in preflight (combine re-validates per-donor regardless). |
| `-h` / `--help` / `HELP=1` | print usage and exit. |

## Attaching to an already-running stage

If a stage is already in the queue (e.g. you launched discovery by hand), don't resubmit it —
tell the script its numeric job id and it chains the next stage onto it:

| var | attach to |
|-----|-----------|
| `DISC_JOBID` | a running discovery array |
| `COMBINE_JOBID` | a running combine_insertions job |
| `GT_JOBID` | a running genotype array |

The script also **auto-detects** a running stage by job name, but an explicit `*_JOBID` is
safer and always wins.

---

## What preflight checks (before anything is submitted)

All checks run first; **every** problem found is reported together in one block, and nothing is
submitted unless all required inputs are present:

- `FOFN` exists and is non-empty, and **every** BAM it lists exists on disk.
- `FOFNDIR` holds at least one `<donor>.bams.fofn` (so `PATIENTS` is non-empty).
- `DISCOVER_CFG` and `GENO_CFG` files exist.
- discovery and genotype binaries exist and are executable (else: `module load rust/1.87.0 && bash cluster/build.sh`).
- `VENV/bin/python` exists; `src/config.py` exists (gitignored — `cp cluster/config.py.grch38 src/config.py`).
- **discovery reference:** if `DISCOVER_CFG` has `splice_hallmark=true`, its `exon_annotation` track exists (the binary `exit(1)`s without it).
- **combine reference / assembly gate:** `install.sh check-config` validates `genome_2bit` + the `hs1→hg38` chain against the first BAM — a wrong pair yields *silently wrong coordinates*, so this fails fast (bypass with `SKIP_ASSEMBLY_CHECK=1`).

---

## How resume / idempotency works

The script classifies each stage as one of:

- **COMPLETE** → skip. Detected by outputs on disk: discovery/genotype need `<id>.txt.gz` for
  every colony in the FOFN; combine needs the contract; combine_genotypes needs the final table.
- **RUNNING** → attach (via `*_JOBID` or name auto-detect); the next stage chains onto it.
- **TODO** → submit, held behind whatever the previous non-skipped stage resolved to.

Because the per-task scripts (`discover_one.sh`, `genotype_one.sh`) already skip colonies whose
output exists and write atomically (`tmp → final`), **resubmitting a partially-complete array
just tops up the missing tasks** — it never recomputes finished ones. So the safe way to recover
from any partial failure is simply to **run the same command again**.

---

## Worked examples

**See the plan without submitting:**
```bash
DRYRUN=1 TAG=noSPEC8 DISCOVER_CFG=cluster/config.discovery.grch38.noSPEC8 \
  bash cluster/submit_all.sh
```

**Fresh full run (all four stages, chained):**
```bash
TAG=noSPEC8 DISCOVER_CFG=cluster/config.discovery.grch38.noSPEC8 \
  bash cluster/submit_all.sh
```

**Discovery already queued as 957812 — attach and chain stages 2–4 behind it:**
```bash
TAG=noSPEC8 DISCOVER_CFG=cluster/config.discovery.grch38.noSPEC8 DISC_JOBID=957812 \
  bash cluster/submit_all.sh
```

**Resume after a partial failure** — just re-run the original command; completed stages are
skipped and only the missing work is submitted.

**A different cohort / non-default inputs:**
```bash
TAG=myrun DISCOVER_CFG=cluster/config.discovery.grch38 \
  FOFN=/lustre/.../mycohort/all.bams.fofn \
  bash cluster/submit_all.sh
```

---

## After it finishes

```bash
bjobs -A                                        # watch the arrays
zcat <FINAL_OUT> | tail -n +2 | wc -l           # loci passing the cohort gates
```

The contract (`$CONTRACT`) is the artefact to diff against the frozen truth set
(`analysis/mei9x10/frozen_truth.tsv`) for the survived-vs-dropped TP/FP comparison; the final
call table (`$FINAL_OUT`) is the output of the cohort gates in `combine_genotypes`.

---

## Relationship to the individual scripts

`submit_all.sh` **orchestrates** the existing per-stage scripts; it does not replace them. Each
stage script remains independently runnable and keeps its own guards:

- `submit_discovery.sh` — discovery array + contig-allowlist preflight.
- `combine_mei.sh` — pooled combine_insertions + per-donor assembly gate (not a submitter; this
  script `bsub`s it).
- `submit_genotype.sh` — genotype array; supports `WAIT=ended(<combine>)` and defers
  contract-dependent checks when chained.
- `submit_combine_genotypes.sh` → `combine_gt.sh` — cohort gates; supports `WAIT`.

Use `submit_all.sh` for the whole run; reach for an individual script when you want to re-run or
debug a single stage in isolation.
