# Fleet workflow — process every GRCh38 patient in `patients/`

Design for a one-command, walk-away run over the whole `patients/` tree on farm22 (LSF).
Reuses `cluster/pipeline.sh` (the per-patient stage→discover→combine→genotype→combine_genotypes DAG)
and wraps it with pre-flight, per-BAM stats, archival, and a QC gate.

Decisions locked with Jeremy (2026-07-18):
- **Scratch = staging root** `/lustre/scratch126/casm/staging/team273/jd43` — disposable, wiped at start (step 1) and per-patient after success (step 9). Nothing precious ever lives here.
- **Persistent lustre working tree** = `pt_runs/<patient>/` — discovery + genotype **intermediates are kept here** ("lustre, not scratch"). Regenerable; may be purged.
- **Backed-up finals** = `~/results/<patient>/` (NFS) — small CSVs only: calls, contract, annotated, per-BAM stats summary, QC line.
- **Scope = GRCh38 WGS only.** Of the 50 patients with any GRCh38 WGS colony: **44 pure-GRCh38 run now** (3996 colonies); **6 mixed-assembly** (PD45534, PD48367, PD48372, AX001, PD43976, PD44579 — only 1–48% of their WGS tree is GRCh38) are routed to the remap queue (`mixed_partial.tsv`), not run partially, unless `INCLUDE_MIXED=1`. **105 non-GRCh38 patients** → `skip_non_grch38.tsv`. Both lists feed the future bwa-remap extension.

## Command surface (minimal input)

```bash
bash cluster/fleet.sh plan          # dry run: build worklist + skip-list, capacity estimate, no jobs
bash cluster/fleet.sh run           # steps 1–9 for all 50 GRCh38 patients
bash cluster/fleet.sh run PD44887   # one (or a few) patients by id
bash cluster/fleet.sh status        # fleet-wide progress table (wraps pipeline.sh status)
```

`run` does steps 1–2 synchronously on the head node (seconds), then submits **one long-running
orchestrator LSF job** (queue `basement`, low-CPU) that loops patients, throttled so only a bounded
number stage concurrently. You run one command and walk away.

## The 9 steps

### 1. Clean the scratch (staging) area
Wipe stale staged BAMs from previous runs before starting.
- Safety gate: refuse if any live LSF job references the staging root (`bjobs -w | grep`), so a
  concurrent run is never clobbered. Otherwise `rm -rf "$STAGING_ROOT"/*`.
- Re-stageable from iRODS, so worst case is re-staging time, not data loss.

### 2. Check lustre capacity
- `lfs quota -u $USER /lustre/scratch126` (user) + team273 quota.
- Estimate peak footprint = `STAGE_THROTTLE × FLEET_PATIENTS_INFLIGHT × avg_BAM(~40 GB)` for staging
  (bounded, because pipeline cleans each patient's BAMs on success) **+** cumulative `pt_runs`
  intermediates (discovery/genotype `.txt.gz`, ~small). Discovery/genotype text outputs are tiny;
  staged BAMs dominate and are transient.
- **Gate:** if free < estimate + 15% headroom → abort with a clear message (this is the one hard stop
  before compute; everything downstream is tolerant).

### 3. Stage BAMs + discover (pipelined)
- Per patient, build `samples.tsv` = **`sample<TAB>proj`** from `colonies.tsv` rows where
  `assembly==GRCh38 && ds~WGS`. **Per-sample project is required** — PD40667/PD43974/PD37111/AX001/PD43976
  each span 2–3 projects.
- `pipeline.sh` phase 1 (`stagedisc[1-N]%STAGE_THROTTLE`): each array element stages ITS bam
  (`stageBam.pl --project <that row's proj>`) then runs discovery the instant it lands.
  Config = `config.discovery.grch38`, `config.py.grch38`.
- Fault-tolerant by design: a listed sample with no BAM drops a marker in `missing/` and exits 0.
- **OOM/timeout escalation** (see "Failure escalation" below): a `sd_retry` controller runs after the
  array and re-runs only the OOM/timeout-killed samples at a higher memory/wall tier before combine.

### 4. combine_insertions
- `pipeline.sh` phase 2 over `discovery/*.txt.gz`. Keeps the **assembly-config gate**
  (`install.sh check-config` against a staged BAM header) — refuses to combine if the config doesn't
  match GRCh38, since combine is assembly-specific and a mismatch is silently wrong.
- Output contract → `pt_runs/<patient>/insertions/<patient>.genotyping.txt.gz`.

### 5. Genotype + per-BAM stats
- `pipeline.sh` phase 3 (`genotype[1-N]%GT_THROTTLE`): index (if needed) + genotype vs the contract.
  Config = `config.genotype.grch38`.
- **NEW — per-BAM stats sidecar** `pt_runs/<patient>/stats/<sample>.tsv`, written in the genotype task
  (BAM is already staged + indexed, so it's nearly free):
  `sample  n_reads  read_len  mean_cov` via one `samtools coverage` + `samtools idxstats`
  (mean_cov = Σ(meandepth·len)/Σlen over primary chromosomes; n_reads/read_len from idxstats/`stats`).
- These are per-BAM and one-file-per-sample (no shared-file writes — Lustre-safe).
- **OOM/timeout escalation**: a `gt_retry` controller mirrors the discovery one, re-running only the
  OOM/timeout-killed samples at a higher tier before combine_genotypes.

### 6. combine_genotypes
- `pipeline.sh` phase 4 over `genotypes/*.txt.gz` → `pt_runs/<patient>/<patient>.genotypes.csv.gz`.

### 7. Save patient data
- **Intermediates stay** in `pt_runs/<patient>/{discovery,insertions,genotypes,stats}` (persistent
  lustre working tree — the "lustre, not scratch" location).
- **Finals copied to NFS** `~/results/<patient>/`: `<patient>.genotypes.csv.gz`, the contract,
  `<patient>.annotated.csv.gz`, plus a **stats summary** `<patient>.bam_stats.tsv` (concatenation of the
  per-BAM `stats/*.tsv`) and the QC line. Small enough for the backed-up NFS quota.

### 8. QC gate (mark, don't abort)
Per patient, after combine_genotypes, compute and append one row to `pt_runs/_fleet/fleet_qc.tsv`:
`patient  n_listed  n_missing  n_discovered  n_genotyped  contract?  calls?  verdict`
- `verdict = OK` if `n_missing/n_listed < 0.05` **and** contract+calls exist **and**
  `n_genotyped ≈ n_discovered`.
- `verdict = WARN:<reason>` otherwise (e.g. `WARN:missing=8%`) — **recorded, run continues.** Never aborts.
- Fleet-level summary at the end: counts of OK / WARN / failed patients.

### 9. Clean the staging area
- `pipeline.sh` phase 5 (`cleanup`, hangs off `done(combine_genotypes)` so a failure never deletes BAMs):
  deletes **only this patient's** staged sample dirs (per-sample `rm`, never `rm -rf $PROJECT/`, which
  would hit sibling patients sharing a project).
- Fleet end: a final `staging_clean.sh`-style sweep of anything left behind.

## Failure escalation (OOM / timeout → respawn bigger, then proceed)

An OOM or wall-clock kill is a `SIGKILL` LSF sends when the task exceeds `-M` (→ `TERM_MEMLIMIT`) or
the queue's wall limit (→ `TERM_RUNLIMIT`). The task **cannot trap it to resubmit itself**, and `-M`
is fixed at submit time. So escalation is driven by a **retry-controller job** inserted into the DAG
between an array phase and the phase that consumes it:

```
sd_array (tier1) ─▶ sd_retry (controller) ─▶ combine        # discovery
gt_array (tier1) ─▶ gt_retry (controller) ─▶ combine_genotypes   # genotype
```

The controller (single low-mem job, depends `ended(array)`) loops escalation tiers until clean or tiers
exhausted:

1. **Find genuine failures** = samples in the list whose output `<sample>.txt.gz` is absent **and** which
   have **no** `missing/<sample>` marker (a missing marker = legit no-data, not a failure).
2. **Classify each** by grepping that element's own LSF log (`logs/sd.<I>.log`/`.err`, one per index):
   - `TERM_MEMLIMIT` → OOM → escalate **memory**.
   - `TERM_RUNLIMIT` → timeout → escalate **wall** (queue).
   - neither (real crash / bad data) → **not retried** — left for the step-8 QC gate to flag. This is the
     safeguard's boundary: only OOM/timeout respawn; genuine bugs don't loop.
3. **Resubmit only those indices** once, as a sparse array (`-J <patient>_sd_r2[3,17,42-44]`, same
   `sed -n ${LSB_JOBINDEX}p` body → identical work, idempotent: existing outputs skip). Tier 2 raises
   **both** mem and queue (a memlimit victim that also would have timed out is covered). **There is no
   tier 3.**

   | phase | tier1 (base) | tier2 (only escalation) |
   |---|---|---|
   | discovery | 16 GB / `normal` (12h) | 32 GB / `long` (48h) |
   | genotype | 2 GB / `normal` | 8 GB / `long` |

4. **Block on the retry array** (`bwait -w "ended(<retry>)"`, or poll `bjobs` if `bwait` is absent), then
   re-classify from the tier-2 logs.
5. **After tier 2 — no further retries:**
   - A sample that *still* fails (or failed tier1 for a non-OOM/timeout reason) is **flagged and excluded**:
     its `(sample, reason)` is written to `pt_runs/<patient>/<phase>_excluded.tsv`, and the phase proceeds
     with the surviving samples. `combine`/`combine_genotypes` glob existing outputs, so excluded samples
     just drop out; the QC gate (step 8) counts them.
   - **If ALL samples fail** (zero outputs produced): the controller writes `FLEET_FATAL.<phase>` with the
     collected reasons and **exits non-zero → aborts the patient** (the `done()` dependency stops `combine`
     and everything downstream). The fleet loop reports it and `bkill`s the patient's orphaned jobs.

Notes: tiers are env-overridable; base values are the measured PD44579 peaks (discovery 12.2 GB, genotype
140 MB) so tier1 already has headroom and escalation should be rare.

## What has to be built / changed

| item | change |
|---|---|
| `cluster/fleet.sh` | **new** orchestrator: steps 1, 2, patient loop w/ throttle, step 8 rollup, step 9 final sweep, `plan`/`run`/`status`. |
| `cluster/pipeline.sh` | (a) new `submit-list <PATIENT_ID> <SAMPLES_TSV>` entry that takes a prebuilt `sample<TAB>proj` list instead of resolving irods.txt; (b) **per-sample project** in `bam_path`/stage/cleanup (read col 2), replacing the single `PROJECT_ID`; (c) default configs → grch38; (d) `stats/<sample>.tsv` write in `cmd_genotype`; (e) **retry-controller** subcommand + `sd_retry`/`gt_retry` jobs wired into the DAG, with the tier table and log-based `TERM_MEMLIMIT`/`TERM_RUNLIMIT` classifier. |
| worklist builder | in `fleet.sh`: scan `patients/*/*/colonies.tsv`, emit per-patient `samples.tsv` + global `skip_non_grch38.tsv`. |
| stats summary + QC | in `fleet.sh` step 7/8, plus a `cmd_qc`/`cmd_archive` in `pipeline.sh` if cleaner as a DAG leaf. |

## Deliberately deferred (per Jeremy)
- **Non-GRCh38 remap extension.** When enabled: for `assembly!=GRCh38` BAMs, insert a bwa-mem→GRCh38
  remap array (samtools `fixmate` + `markdup`) before discovery, matching the exact aligner + post-align
  steps the Sanger used on a true GRCh38 BAM (read `@PG` bwa/samtools versions & args, plus `@RG`, from a
  known-good GRCh38 header and mirror them). Skip-list from step 3 is the input queue for this.
