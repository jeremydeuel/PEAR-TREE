# Running PEAR-TREE on farm22 — donor PD44579

Runbook for genotyping one donor's colonies end-to-end on the Sanger farm (LSF).
Worked example: **PD44579**, ~185 colony WGS BAMs under
`/lustre/scratch126/casm/staging/team273/jd43/2178`.

## Assembly: these BAMs are GRCh37 / hs37d5 (already fully supported)

The dupmarked BAMs are **bwa-mem mapped to `1000Genomes_hs37d5`** (`@SQ … AS:NCBI37`,
contigs `1..22,X,Y,MT` + `GL000*`/`NC_007605`/`hs37d5` decoys).

- **Discovery** is assembly-agnostic → `min_mapq=60` + a GRCh37 contig allowlist
  (`config.discovery.grch37`).
- **combine_insertions** always remaps clips to **hs1 (T2T)** and reconciles them to the
  BAM assembly with a pyliftover chain; the `genome_2bit` wrapper maps `1`→`chr1`
  automatically. So GRCh37 needs **no new index** — only the hs1 index + the
  **hs1→hg19** chain, both already on the farm (step 5). hg19 coordinates equal
  GRCh37/hs37d5 for `1..22,X,Y`.

---

## Modules on farm22

Tools come from environment modules with **versioned** names, e.g.:

```bash
module load samtools-1.19      # provides samtools (needed for quickcheck + combine_insertions)
module load bowtie2            # (or a versioned bowtie2-x.y.z — check `module avail bowtie2`)
```

Use `which samtools` / `which bowtie2` after loading to get the absolute paths for
`config.py['combine_insertions']`.

## 0. Deploy + build (head node, has internet)

```bash
git clone -b PEAR-TREE2 git@github.com:limebutterfly/PEAR-TREE.git
cd PEAR-TREE
module load rust/1.87.0        # farm22: provides cargo 1.87 (see `module avail rust`)
bash cluster/build.sh          # -> peartree-discovery + peartree-genotype release binaries
```

> **CARGO_HOME trap (farm22):** the `rust/1.87.0` module sets `CARGO_HOME` to its own
> **read-only** install dir, so `cargo build` fails with `Permission denied` /
> `failed to create directory …/registry/cache`. `build.sh` auto-detects this and
> redirects `CARGO_HOME` to a writable local cache. If you build cargo by hand, first:
> `export CARGO_HOME=/lustre/scratch126/casm/teams/team273/users/jd43/.cargo`
> (a writable lustre path — home quota is too small for the crate cache). The build
> fetches crates from crates.io, so run it on the head node (compute nodes have no
> outbound internet).

> `src/config.py` is intentionally **not** in the repo (gitignored). The Rust discovery
> step does not use it. The Python combine/genotype steps do — create it in step 5.

## 1. File list (exclude the staging duplicates)

```bash
bash cluster/make_filelist.sh /lustre/scratch126/casm/staging/team273/jd43/2178
# -> bams.fofn : one <id>/mapped_sample/<id>.sample.dupmarked.bam per line
wc -l bams.fofn
```

The `tmpExportData/progress/` copies are excluded automatically.

## 2. Pilot one colony (smoke test + right-size resources)

```bash
bash cluster/pilot.sh
# when it finishes:
grep -E 'Max Memory|CPU time|Successfully' logs/pilot.out
zcat discovery/<firstID>.txt.gz | grep -c '^>'    # candidate breakpoints found
```

This confirms the full path works on real GRCh37 WGS before you launch 185 jobs, and the
LSF report gives the peak RSS / CPU time to set `MEM` for the array. Discovery is
single-threaded and streaming (no BAM index needed).

## 3. Discovery — LSF array, 1 core/file, throttled

```bash
MEM=<from pilot, e.g. 6000> THROTTLE=50 bash cluster/submit_discovery.sh
bjobs -A                              # array progress
ls discovery/*.txt.gz | wc -l         # expect the bams.fofn count
```

One core per colony (file-level parallelism is linear; intra-file threading is not).
`THROTTLE` caps concurrency so 185 full-BAM scans don't saturate Lustre — the bottleneck
here is shared I/O, not CPU. Re-running the array is safe: finished colonies are skipped.

---

## 4–7. Combine → genotype → combine (Python; needs GRCh37 indices)

These steps use the Python driver and therefore `src/config.py`. First create it:

```bash
cp src/config_hs.py src/config.py     # then edit the combine_insertions block (below)
```

### 5. Point config.py['combine_insertions'] at the existing farm files

No index build needed — everything is already under
`/lustre/scratch126/casm/teams/team273/users/jd43/`. Set:

```python
'combine_insertions': {
    'genome_2bit':      '/lustre/scratch126/casm/teams/team273/users/jd43/hg19.2bit',
    'bowtie2_index':    '/lustre/scratch126/casm/teams/team273/users/jd43/pt_hu_trees/hs1/hs1',
    'bowtie2_index2':   '/lustre/scratch126/casm/teams/team273/users/jd43/pt_hu_trees/hs1/hs1',
    'bowtie2_index2_lo':'/lustre/scratch126/casm/teams/team273/users/jd43/hs1.hg19.all.chain.gz',
    'samtools_executable': '<`module load samtools`; which samtools>',
    'bowtie2_executable':  '<`module load bowtie2`;  which bowtie2>',
    'exclude_files_with_many_insertions': 1_000_000,
    'clean_remap_max_insertion': 12,
    'clean_remap_min_as': -15,
},
```

Why these are correct for GRCh37 BAMs:
- `bowtie2_index` / `bowtie2_index2` = **hs1 (T2T)** — clips always remap to the most
  complete reference; assembly-independent.
- `bowtie2_index2_lo` = **hs1→hg19 chain** — lifts hs1 clip-hits to hg19 coordinates,
  which equal the GRCh37/hs37d5 BAM coordinates for `1..22,X,Y`. This is the piece that
  makes it GRCh37-correct.
- `genome_2bit` = **hg19.2bit** — supplies the reference-flank consensus at each
  breakpoint (BAM coords); `combine_insertions_get_sequence.get_sequence` maps the
  numeric BAM names `1`→`chr1` for the fetch.

> Validate combine on the pilot colony's discovery output before the full run.

Optional (element classification, post-genotyping `tools/annotate_v2.py`): the Dfam
`hs` HMM (`…/jd43/Dfam-curated_only-hs.hmm`, already hmmpressed), `dfamscan.pl`, and
`hmmer-3.3.2` are all present under `…/jd43/` if you want it later — not needed for the
core call table.

### 6. combine_insertions (build the shared genotyping contract)

```bash
python src/main.py --step combine_insertions \
    --discovery_files discovery/*.txt.gz \
    --out PD44579 --threads 8
# -> PD44579.genotyping.txt.gz  (the contract genotyped against every colony)
```

### 7. Genotype every colony against the contract (Rust batch — CONSECUTIVE)

```bash
# manifest: one "bam<TAB>output" per colony
awk '{ id=$0; sub(/.*\//,"",id); sub(/\.sample.*/,"",id);
       print $0 "\t" "genotypes/" id ".txt.gz" }' bams.fofn > genotype.manifest.tsv
mkdir -p genotypes

# BAMs need an index for genotyping (unlike discovery). Index any that lack one:
#   for b in $(cut -f1 genotype.manifest.tsv); do [ -f "$b.bai" ] || samtools index "$b"; done

rust/peartree-genotype/target/release/peartree-genotype \
    --step genotype_batch --manifest genotype.manifest.tsv \
    --insertions PD44579.genotyping.txt.gz --threads 4 \
    --config cluster/config.genotype.grch37      # min_mapq=60 (create like the discovery cfg)
```

Runs one colony at a time, 4 cores across its loci — bounded memory/handles. Wrap the
whole command in a single `bsub -n 4` job. (Or, for max throughput, submit it as an LSF
array over the manifest with 1 core each, like discovery.)

### 8. combine_genotypes (final call table)

```bash
python src/main.py --step combine_genotypes \
    --genotypes genotypes/*.txt.gz --out PD44579.calls.csv.gz --threads 8
```

---

## Scheduler note

Scripts here use LSF (`bsub`) for farm22. The per-file work scripts
(`discover_one.sh`) are scheduler-agnostic; only the submit wrappers are LSF-specific.
The SLURM equivalent is in the header of `submit_discovery.sh`.
