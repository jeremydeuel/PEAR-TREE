# TPRT A/B kit — test the TPRT-hallmark pipeline on colon colonies (farm22)

One command per patient runs **two complete pipelines side by side** on the same staged BAMs and
then scores both against the patient's phylogeny:

| arm | what | discovery / genotype config | `src/config.py` from | checkout |
|---|---|---|---|---|
| **A** | current pipeline | `cluster/config.discovery.grch38` / `config.genotype.grch38` | `cluster/config.py.grch38` | `tprt_ab/PEAR-TREE-A` (git worktree) |
| **B** | TPRT-hallmark pipeline | `cluster/config.discovery.grch38.tprt` / `config.genotype.grch38.tprt` | `cluster/config.py.grch38.tprt` | `tprt_ab/PEAR-TREE` (clone, branch `tprt-hallmarks`) |

Same commit, same binaries, same BAMs — only the configuration differs (every TPRT key defaults
off in code, so arm A is the pre-TPRT pipeline). Two checkouts because `src/config.py` is one
gitignored module that `main.py` / `annotate_v2.py` import from their own `src/`, and the arms
run concurrently. Both arms are annotated with the **same** annotate settings (rte_library + hs1
gene model; `cluster/tprt/arm_config.py`): annotation is the measuring instrument, not part of
the pipeline under test. Each arm gets its own annotate scratch dir (`tprt_ab/tmp/{A,B}`) — the
base configs share `jd43/tmp/<patient>.*`, which two concurrent annotates would overwrite.

Everything lives under `/lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab`:

```
tprt_ab/PEAR-TREE/            arm B checkout (binaries built here, used by both arms)
tprt_ab/PEAR-TREE-A/          arm A worktree, detached at the same commit
tprt_ab/PEAR-TREE/venv/       shared venv (PEAR-TREE-A/venv is a symlink to it)
tprt_ab/bin/minimap2  tprt_ab/build.stamp
tprt_ab/resources/            hs1.2bit  hs1.gene_model.tsv.gz  hs1.sr.mmi (optional)
tprt_ab/tmp/{A,B}/            annotate scratch
tprt_ab/<P>/samples.tsv       sample<TAB>project (from fleet.sh samples <P>)
tprt_ab/<P>/{A,B}/<P>/        each arm's pipeline run dir (discovery/ insertions/ genotypes/ logs/ ...)
tprt_ab/<P>/{A,B}/jobids.tsv  phase<TAB>LSF job id of the last submission
tprt_ab/<P>/eval/             tree_fit + discrimination per arm, ab_report.md
~/results/tprt_ab/<P>/        NFS copies: each arm's finals (A/, B/) + the report
```

## Which patient first

| patient | colonies | tree | staging (est.) | role |
|---|---|---|---|---|
| **PD37449** (colorectum) | 10 | 10 tips, resolved clades | ~0.5 TB | **pilot**: small, informative tree, measures arm B's real memory/runtime |
| PD34200 (colorectum) | 5 | 5 tips, poorly supported | ~0.25 TB | optional smoke test (fastest; tree too small to judge much) |
| **PD44890** (colon) | 49 | deep, long-branch clades | ~4.5 TB (one colony has 3.5 G reads) | main test #1 |
| **PD37590** (colorectum) | 44 | adjusted-subs tree | ~4.1 TB | main test #2 |
| PD44887 (colon) | 32 | yes | ~2.9 TB | third |

All five: every tree tip is a GRCh38 WGS colony in `colonies.tsv` (checked by preflight).
Run PD37449 first, read its `ab_report.md` (memory table!), lower/raise the budgets, then
PD44890 and PD37590 **one at a time with `--cleanup`** (team273 quota: ~15 TB free in July).

## 0. Deploy (head node — has internet)

This needs the kit commit on GitHub (push `tprt-hallmarks` from the Mac first).

```bash
mkdir -p /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab
git clone -b tprt-hallmarks git@github.com:limebutterfly/PEAR-TREE.git /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE
bash /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/cluster/tprt/setup.sh all
```

`setup.sh all` = these steps, each also runnable alone and idempotent:

| step | does | check by hand |
|---|---|---|
| `worktree` | `git worktree add --detach tprt_ab/PEAR-TREE-A <HEAD>` (re-run after every `git pull` to move it along) | `git -C /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE-A log -1 --oneline` |
| `build` | `module load rust/1.87.0` + `cluster/build.sh` (redirects the read-only `CARGO_HOME`), writes `build.stamp` (commit + md5) | `cat /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/build.stamp` |
| `venv` | `module load python/3.12.3`; venv with `requirements.txt` (incl. `edlib`, `mappy`) + `scipy`, `matplotlib` | `/lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/venv/bin/python -c "import edlib, mappy, scipy, matplotlib, pysam"` |
| `minimap2` | a `minimap2` module if one exists, else the official static binary v2.28 from github.com/lh3/minimap2 | `/lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/bin/minimap2 --version` |
| `resources` | links an existing `hs1.2bit` (or downloads it from UCSC); downloads `hs1.ncbiRefSeq.gtf.gz` and builds the gene model (below) | `zcat /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/resources/hs1.gene_model.tsv.gz \| wc -l` (≈252,903) |
| `configs` | writes the generated `src/config.py` into both checkouts (arm A / arm B), creates `tmp/{A,B}` | `head -12 /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/src/config.py` |
| `index` | downloads `hs1.fa.gz`, **bsubs** `minimap2 -x sr -d hs1.sr.mmi` (4 cores, 32 GB, ~30-60 min) | `ls -la /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/resources/hs1.sr.mmi` |

The gene-model build `resources` runs is exactly (E2E_REPORT "Annotate round 2"):

```bash
curl -fL -o /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/resources/hs1.ncbiRefSeq.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz
/lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/venv/bin/python /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/tools/build_gene_model.py --curated /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/resources/hs1.ncbiRefSeq.gtf.gz /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/resources/hs1.gene_model.tsv.gz
```

Unchanged farm resources the configs already name (checked by preflight): `jd43/hg38.2bit`,
`jd43/hs1.hg38.all.chain.gz`, the hs1 bowtie2 index `jd43/pt_hu_trees/hs1/hs1`, bowtie2 2.5.4,
the GRCh38 discovery exon track `jd43/grch38.exons.bed.gz`, hs1 RepeatMasker, Dfam HMM +
dfamscan + HMMER. `resources/rte_library` ships with the repo (md5-checked against its
`manifest.tsv`). The hs1 minimap2 index is optional: without it annotate's novel-source locator
is off (known transduction sources still work); preflight warns.

**src/config.py is generated, never hand-edited.** It is a stub that loads the canonical
`cluster/config.py.grch38` (A) / `cluster/config.py.grch38.tprt` (B) and applies the short override
list in `cluster/tprt/arm_config.py` (annotate tmp, rte_library, hs1 gene model, hs1.2bit,
optional .mmi). An existing hand-made `src/config.py` is kept as `src/config.py.pre_tprt_ab.<time>`.
Print what each arm resolves to:
`/lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/venv/bin/python /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/cluster/tprt/arm_config.py`

## 1. Preflight (head node, per patient)

```bash
bash /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/cluster/tprt/preflight.sh PD37449
```

Checks, each `OK` / `WARN` / `FAIL` (exit 1 on any FAIL): arm-A worktree at the same commit;
binaries built from this commit (`build.stamp`, mtime vs the last `rust/` commit) **and accepting
every key of all four Rust configs** (each binary is run on each config and must print no
`ignoring unknown config key` — config.rs only warns, so a stale binary silently runs arm B as
arm A; a pre-TPRT build ignores 20 keys of the .tprt discovery config); venv imports; each `src/config.py` is the right arm;
`tools/phylo/{tree_fit,discrimination}.py --help`; minimap2; rte_library md5s; hs1 gene model /
hs1.2bit / .mmi; the discovery exon tracks and every combine/annotate path in both configs;
tree tips vs GRCh38 WGS colonies (both directions); each BAM staged + `samtools quickcheck` +
indexed, else header-readable in nst_links (header reads only — Sanger policy); and the assembly
gate `install.sh check-config --bam <a real BAM>` for **both** checkouts. `run_ab.sh` runs the
`--quick` subset (no BAM loop, no check-config) itself and refuses on a FAIL.

## 2. Submit both arms

```bash
bash /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/cluster/tprt/run_ab.sh PD37449 --dry-run
bash /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/cluster/tprt/run_ab.sh PD37449
```

`--dry-run` prints every `bsub` with fake job ids (a printing `bsub` shim on `PATH`; the arm run
dirs and `samples.tsv` are still written — same content a real submit writes). Other options:
`--arm A|B|both` (default both), `--cleanup` (delete the patient's staged BAMs after **both** arms'
combine_genotypes succeeded), `--no-eval`, `--force` (resubmit although jobs of the arm's last
submission are still pending/running), `--skip-preflight`.

What it resolves from files: `patients/<organ>/<P>/` (exactly one `*.tree`, `colonies.tsv`), the
sample list through `cluster/fleet.sh samples <P>` (the fleet's GRCh38 & WGS rule — the same rows
`fleet.sh run` would use), per-sample iRODS project → staged BAM
`/lustre/scratch126/casm/staging/team273/jd43/<proj>/<sample>/mapped_sample/<sample>.sample.dupmarked.bam`
(pipeline.sh's stageBam.pl convention). Each arm is one ordinary `cluster/pipeline.sh submit-list`
DAG with its own `WORKROOT`, `RESULTS_DIR`, configs and job-name prefix (`PD37449_A_*`, `PD37449_B_*`):

```
PD37449_A_sd[1-N]%20 -> _sdr -> _ci -> _gt[1-N]%12 -> _gtr -> _cg -> _annotate ─┐
PD37449_B_sd[1-N]%20 -> _sdr -> _ci -> _gt[1-N]%12 -> _gtr -> _cg -> _annotate ─┴─> PD37449_tprt_eval
   (B_sd element i waits ended(A_sd element i): A stages each BAM once, B reuses it)
```

pipeline.sh hooks used (all unset = unchanged default behaviour, verified: identical `bsub` lines):
`PT_JOB_PREFIX`, `PT_SD_WAIT`, `PT_NO_CLEANUP=1` (the pipelines' own cleanup would let the first
arm to finish delete the BAMs the other is still genotyping), `PT_JOBIDS_FILE`.

### Resources (MB; env overrides in `run_ab.sh`)

| stage | arm A | arm B | basis |
|---|---|---|---|
| stage+discover | 8000 | 12000 | post-jemalloc peak max 3.1 GB / p95 2.1 GB (job 954062, 90 BAMs). B adds the evidence sidecar, all mates, SHORT reads and a per-colony floor of ONE fragment — unmeasured on real WGS. TERM_MEMLIMIT → tier-2 retry at 32000 on `long` |
| combine_insertions (8 cores) | 24000 | 64000 | legacy 15.5 GB on 174 files. B pools every colony's sidecar reads in Python — unmeasured, and **combine has no retry controller** |
| genotype | 4000, THROTTLE 12 | 4000, THROTTLE 12 | 140 MB measured; the 4 GB reservation + low throttle stop LSF from packing the array onto one big node (the I/O-starvation "runs forever" tail) |
| combine_genotypes / annotate / eval | 8000 / 32000 / 16000 | same | measured 1.3 GB / — / small |

Run time (planning): discovery ~24 min per 30x colony (scales with reads: PD44890 averages ~800 M
reads, one colony 3.5 G → hours), genotype ~18 min/colony, combine 15-60 min, annotate ~0.5-2 h.
PD37449 ≈ 4-6 h wall, ~20 CPU-h for both arms; PD44890 ≈ 12-16 h wall, ~150 CPU-h.
Override example: `B_CI_MEM=96000 bash /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/cluster/tprt/run_ab.sh PD37449 --arm B`.

### Watch

```bash
bjobs -w | grep PD37449_
for a in A B; do PT_JOB_PREFIX=PD37449_$a WORKROOT=/lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PD37449/$a bash /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/cluster/pipeline.sh status PD37449; done
```

Re-running `run_ab.sh` resumes (every pipeline task skips existing outputs). A failed
combine (e.g. TERM_MEMLIMIT in arm B) leaves the rest of that arm pending on `done()`: `bkill` the
arm's pending jobs, raise `B_CI_MEM`, re-run with `--arm B`.

## 3. Evaluation (submitted automatically; re-runnable by hand)

```bash
bash /lustre/scratch126/casm/teams/team273/users/jd43/tprt_ab/PEAR-TREE/cluster/tprt/evaluate.sh PD37449
```

Per arm: `tools/phylo/tree_fit.py --genotypes <P>.genotypes.csv.gz --genotype-dir genotypes/
--tree <patients tree> --annotation <P>.annotated.csv.gz --samples colonies.tsv --sex auto`
→ `eval/<arm>/fit/phylo_fit.tsv`, then `tools/phylo/discrimination.py --fit … --out
eval/<arm>/discrimination`, then `cluster/tprt/compare_arms.py` → `eval/ab_report.md` with:
loci per stage (discovery per colony / pooled, combined, contract, +one-sided, calls, annotated;
B's combine-gate fail reasons), locus kinds from the name geometry, class / element / tprt_call,
carrier classes (private / shared / germline-like) × phylo buckets (consistent / violating /
other) per arm, fuzzy ±10 bp one-to-one A↔B overlap at contract and call level (`a_only.tsv`,
`b_only.tsv` with phylo label, tprt_score, carriers; `matched_*.tsv`), and peak memory / run time
/ CPU / TERM_* per stage from the LSF reports in each arm's `logs/*.log`. Everything is copied to
`~/results/tprt_ab/<P>/`. `SEX=F` (or `M`) overrides tree_fit's `--sex auto`.

`tools/phylo/` is being written in parallel; until it is on the branch, preflight warns, the arms
run, and the eval job fails at tree_fit — `git pull`, `setup.sh worktree`, then run `evaluate.sh`
by hand (compare_arms still writes the report without phylo labels when tree_fit fails; the
script then exits 2 so the LSF job shows as failed rather than done).

## Files

| file | role |
|---|---|
| `setup.sh` | deploy: worktree, build + stamp, venv, minimap2, resources, generated configs, hs1 index |
| `preflight.sh` | head-node checks (above); `--quick` for run_ab |
| `run_ab.sh` | one patient: both arm DAGs + eval (+ cleanup) with dependencies; `--dry-run` |
| `evaluate.sh` | tree_fit + discrimination per arm, compare_arms, copy to results |
| `compare_arms.py` | the A/B report |
| `arm_config.py` | builds each arm's `CONFIG` from the canonical config + overrides |
| `common.sh` | paths (env-overridable) and the patient / tree / sample resolvers |
