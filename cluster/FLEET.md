# FLEET — operator manual

One-command sweep over every GRCh38 patient in `patients/` on farm22. Wraps the per-patient DAG in
`pipeline.sh`. For the full design + rationale see [`FLEET_PLAN.md`](FLEET_PLAN.md).

## TL;DR

```bash
cd <the GRCh38 PEAR-TREE checkout>          # rust binaries + venv must be built
bash cluster/fleet.sh plan                   # dry run: worklist + capacity, submits NOTHING
bash cluster/fleet.sh run PD44887            # smoke test: one pure-GRCh38 patient
bash cluster/fleet.sh run                    # full sweep: all 44 pure-GRCh38 patients
bash cluster/fleet.sh status                 # QC table (OK / WARN / FAIL per patient)
```

`run` does steps 1–2 on the head node (seconds), then submits one per-patient LSF DAG each, throttled so
only `PATIENTS_INFLIGHT` (default 4) stage at once. Walk away; re-run `status` to watch it drain.

## What it does (9 steps)

1. **Clean scratch** — wipes `STAGING_ROOT` (re-stageable from iRODS). Refuses if a live LSF job
   references it.
2. **Capacity gate** — reads the **team273 group** quota (`lfs quota -g team273`; the user quota is
   unlimited). Hard-stops only if free < ~4.1 TB estimate + 15%. This is the one hard stop.
3. **Stage + discover** — per BAM, `stageBam.pl` then discovery the instant it lands (`config.discovery.grch38`).
4. **combine_insertions** — with the assembly-config gate (refuses on a GRCh38 mismatch).
5. **Genotype + per-BAM stats** — writes `sample / n_reads / read_len / mean_cov` per colony.
6. **combine_genotypes** → `<patient>.genotypes.csv.gz`.
7. **Save** — intermediates KEPT in `pt_runs/<patient>/`; small finals + `<patient>.bam_stats.tsv` +
   QC copied to `~/results/<patient>/` (NFS, backed up).
8. **QC gate** — per-patient row in `pt_runs/_fleet/fleet_qc.tsv`; `WARN` if >5% BAMs missing, **never aborts**.
9. **Clean staging** — per-patient on success; final sweep reports leftovers.

## Scope

Only `colonies.tsv` rows with `assembly==GRCh38 && ds~WGS`.

| bucket | count | where | action |
|---|---|---|---|
| pure-GRCh38 | 44 patients / 3996 colonies | `pt_runs/_fleet/patients.txt` | run now |
| mixed-assembly | 6 (PD45534, PD48367, PD48372, AX001, PD43976, PD44579) | `mixed_partial.tsv` | remap queue; `INCLUDE_MIXED=1` forces partial run |
| non-GRCh38 | 105 | `skip_non_grch38.tsv` | remap queue |

The two queues feed the future bwa-remap-to-GRCh38 extension (not implemented).

## OOM / timeout retries

Each array is followed by a controller (`sdr`, `gtr`) that classifies failures from each element's LSF
log and resubmits **only OOM (`TERM_MEMLIMIT`) / timeout (`TERM_RUNLIMIT`) ones, once, at tier 2**:

| phase | tier 1 | tier 2 (only escalation) |
|---|---|---|
| discovery | 16 GB / `normal` | 32 GB / `long` |
| genotype | 2 GB / `normal` | 8 GB / `long` |

**No tier 3.** A sample still failing after tier 2, or failing for any other reason, is flagged and
**excluded** (`<phase>_excluded.tsv`); the phase proceeds on the survivors. If **every** sample fails,
the controller writes `FLEET_FATAL.<phase>` and exits non-zero → the patient aborts and downstream stops.

## Knobs (env overrides)

| var | default | meaning |
|---|---|---|
| `PATIENTS_INFLIGHT` | 4 | max patients staging concurrently |
| `STAGE_THROTTLE` | 20 | per-patient concurrent `stageBam.pl` |
| `INCLUDE_MIXED` | 0 | 1 = also run the 6 mixed patients' GRCh38 subset |
| `SD_MEM_T1/T2` | 16000/32000 | discovery tier1/tier2 MB |
| `GT_MEM_T1/T2` | 2000/8000 | genotype tier1/tier2 MB |
| `RETRY_QUEUE` | `long` | tier2 wall queue |
| `AVG_BAM_GB` | 45 | capacity estimate per BAM |
| `WORKROOT` / `STAGING_ROOT` / `RESULTS_DIR` | jd43 paths | run tree / scratch / NFS finals |

## Troubleshooting

- **Patient stuck / partial:** `bash cluster/pipeline.sh status <PAT>` — shows discovered/genotyped counts,
  `*_excluded.tsv`, `FLEET_FATAL.*`, and `bjobs` for that patient. Re-run `fleet.sh run <PAT>` — every task
  is idempotent and resumes.
- **`FLEET_FATAL.discover`:** all discoveries failed. Check `pt_runs/<PAT>/logs/sd*.{log,err}`; usually a
  bad config, unbuilt binary (`cluster/build.sh` after any `rust/` pull — see the discovery banner), or a
  wrong-assembly BAM the config gate rejected.
- **Capacity abort:** team273 quota is shared and near-full — lower `PATIENTS_INFLIGHT`/`STAGE_THROTTLE`,
  or run `bash cluster/staging_clean.sh --delete` first.
- **`bwait` absent:** the controllers fall back to polling `bjobs`; no action needed.
