# Running PEAR-TREE on farm22 — donor PD44579

Runbook for genotyping one donor's colonies end-to-end on the Sanger farm (LSF).
Worked example: **PD44579**, ~185 colony WGS BAMs under
`/lustre/scratch126/casm/staging/team273/jd43/2178`.

## ⚠️ Assembly: these BAMs are GRCh37 / hs37d5

The dupmarked BAMs are **bwa-mem mapped to `1000Genomes_hs37d5`** (`@SQ … AS:NCBI37`,
contigs `1..22,X,Y,MT` + `GL000*`/`NC_007605`/`hs37d5` decoys). This drives two things:

- **Discovery** is assembly-agnostic → runs as-is with `min_mapq=60` and a GRCh37
  contig allowlist (see `config.discovery.grch37`). Ready now.
- **combine_insertions is NOT** → it remaps clipped consensuses with bowtie2 and lifts
  coordinates via a chain file, and the repo `config.py` points all of that at **hs1
  (T2T)**. Against GRCh37 reads that is wrong. Steps 2–4 are blocked until you build a
  **GRCh37/hs37d5 bowtie2 index + 2bit** (see step 5). Do NOT run combine against the
  hs1 config on this data.

---

## 0. Deploy + build (head node, has internet)

```bash
git clone -b PEAR-TREE2 git@github.com:limebutterfly/PEAR-TREE.git
cd PEAR-TREE
# Rust >= 1.87 (module load a rust, or per-user rustup). Build fetches noodles from crates.io.
bash cluster/build.sh          # -> peartree-discovery + peartree-genotype release binaries
```

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

### 5. Build the GRCh37 references combine_insertions needs (one-off)

Locate `hs37d5.fa` on the farm (the BAM header points at
`…/1000Genomes_hs37d5/all/fasta/hs37d5.fa`; the canpipe copy is
`/lustre/scratch119/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/genome.fa`).
Then, once:

```bash
bowtie2-build --threads 8 hs37d5.fa /path/to/idx/hs37d5      # ~1-2 h, ~4 GB output
faToTwoBit hs37d5.fa /path/to/idx/hs37d5.2bit                # for reference-flank fetch
```

Point `src/config.py['combine_insertions']` at them:

| key | set to |
|-----|--------|
| `genome_2bit` | `/path/to/idx/hs37d5.2bit` |
| `bowtie2_index` | `/path/to/idx/hs37d5` (assembly the BAMs are on — used to reject reference-matching clips) |
| `bowtie2_index2` | `/path/to/idx/hs37d5` (same assembly for a first pass; a newer annotated genome only if you want element classification) |
| `bowtie2_index2_lo` | the chain from `bowtie2_index2` → the BAM assembly; **identity when both are GRCh37** — leave as-is only if `index2` == `bowtie2_index` |
| `samtools_executable`, `bowtie2_executable` | valid farm paths (module load, or absolute) |

> Validate combine on the pilot colony's discovery output before the full run — it's the
> step most sensitive to the assembly/config mismatch.

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
