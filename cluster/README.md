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

### install.sh — check what's present, fetch what isn't

```bash
bash cluster/install.sh check   --bam <a real BAM>   # verify only; exits 1 if anything is missing
bash cluster/install.sh install --bam <a real BAM>   # same, but downloads what's missing
bash cluster/install.sh install --bam <BAM> --annotate   # also the step-3 (annotate) resources
bash cluster/install.sh detect-assembly <BAM>        # -> hg19 | hg38 | mm10 | mm39
bash cluster/install.sh build-exons <asm> --bam <BAM>    # Feature B exon model (~2 min)
```

Passing `--bam` is what makes the check meaningful: the assembly is read from the
BAM rather than assumed. It checks the toolchain (samtools, bowtie2, the two Rust
binaries, the venv) and the three assembly-specific references, printing a ready-to-paste
`config['combine_insertions']` block when everything lines up. Downloads come from UCSC
goldenPath; run it on the **head node** (compute nodes have no outbound internet).

**Supported BAM assemblies**, and what each needs:

| BAM assembly | `genome_2bit` | clips remap to | `bowtie2_index2_lo` (chain) |
|---|---|---|---|
| hg19 (= GRCh37/hs37d5) | `hg19.2bit` | hs1 | `hs1ToHg19.over.chain.gz` |
| hg38 | `hg38.2bit` | hs1 | `hs1ToHg38.over.chain.gz` |
| mm10 (= GRCm38) | `mm10.2bit` | mm39 | `mm39ToMm10.over.chain.gz` |
| mm39 (= GRCm39) | `mm39.2bit` | mm39 | *identity chain, synthesised* |

Assemblies are identified by **chr1 length**, not contig name — so a GRCh37 BAM naming
it `1` and an hg19 BAM naming it `chr1` both resolve to `hg19`, which is right: they are
the same coordinates. hs37d5 only adds decoy contigs, so it takes the hg19 resources too.

For **mm39** the remap target *is* the BAM assembly, so there is nothing to lift and no
chain to download; `install.sh` synthesises a per-contig identity chain from
`mm39.chrom.sizes` (pyliftover always wants a chain). Verified: every coordinate maps to
itself.

The **bowtie2 index is built, not downloaded** (UCSC doesn't ship one). If it's missing:
`bash cluster/install.sh build-index hs1` submits a 16-core / 64 GB / `long` job (~3 h).
Once per species, not per run.

### Feature B (splice / processed-pseudogene) needs an exon model

`config.discovery.*` ships with `splice_hallmark = true`, which makes discovery emit a
per-colony `<out>.splice.tsv` of processed-pseudogene candidates
(`contig breakpoint side gene n_exons intron_bp span_bp`). It requires `exon_annotation`
— and **without it every discovery job dies at startup**, so `install.sh check` now
verifies it (only when the config actually enables the feature).

```bash
bash cluster/install.sh build-exons hg19 --bam <a real BAM>   # ~2 min, once per assembly
```

Source is the **Ensembl** GTF (release-75 for GRCh37, 102 for GRCm38), because Ensembl
contig names are numeric (`1,2,..,X,Y`) and so match hs37d5/GRCh37 BAMs directly;
`build-exons` reads the BAM header and adds a `chr` prefix if the BAM is UCSC-style.
The model is **per-gene merged exons** (every transcript's exons collapsed into
non-overlapping blocks) — raw transcript exons overlap each other and blur the
"mates span ≥ N exons, intron skipped" signature. PD44579/GRCh37 → 354,207 blocks.

> The naming must match the BAM because the model is looked up by the **mate read's
> contig**. This is also why the discovery exon file is GRCh37 while
> `config.py['annotate']['exon_annotation']` must be **hs1**-based — annotate works in
> remap space, discovery works in BAM space. They are not interchangeable.

> **Why this gate exists.** Discovery and genotyping only ever read the BAM, so they don't
> care about the assembly. `combine_insertions` does: it reads reference flanks from
> `genome_2bit` and lifts hs1 clip hits back onto the BAM's coordinates through the chain.
> Point either at the wrong assembly and **you get no error — you get plausible, silently
> wrong coordinates.** So `pipeline.sh` re-derives the assembly from a staged BAM and
> asserts `src/config.py` matches *before* combine runs (`install.sh check-config --bam …`).
> It compares basenames, so where you keep your references doesn't matter. Escape hatch for
> an unusual layout: `PT_SKIP_ASSEMBLY_CHECK=1`.

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
# discovery output is FASTQ (@contig:l-r:LEFT|RIGHT:...), not the '>'-headed contract:
zcat discovery/<firstID>.txt.gz | grep -oE '^@[^:]+:[0-9]+-[0-9]+' | sort -u | wc -l   # candidate breakpoints
```

This confirms the full path works on real GRCh37 WGS before you launch the array, and the
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

## Measured resources (PD44579: 174 colonies × 30,025 loci, GRCh37 30× WGS)

Peak RSS from the LSF reports of a real full run — `pipeline.sh` defaults are set from
these, not guessed. Re-measure on new data (`pilot.sh` reports discovery's peak).

| step | peak RSS | time | cores | budget |
|---|---|---|---|---|
| discovery (Rust) | **12.2 GB** | ~16–24 min/colony | 1 | `SD_MEM=16000` |
| combine_insertions | **15.5 GB** | 854 s | 8 | `CI_MEM=24000` |
| genotype (Rust) | **140 MB** | ~17 min/colony (1 core) | 1 | `GT_MEM=2000` |
| combine_genotypes | **1.3 GB** | 456 s | 8 | `CG_MEM=8000` |

Two things worth internalising:

- **Discovery really does need ~12 GB** (the mate-qname sets). An 8 GB budget kills every
  task with `TERM_MEMLIMIT` — this was the old `submit_discovery.sh` default and is fixed.
- **The Rust genotyper is memory-trivial (140 MB).** Don't over-provision it: `rusage[mem]`
  *reserves* memory, so an inflated `GT_MEM` throttles how many array elements the
  scheduler co-locates. Over-provisioning costs throughput, not just politeness.

Scale: the whole donor (discovery → calls) is a few hours wall at ~50-wide.

## Operational gotchas (learned the hard way on PD44579)

- **Discovery output is FASTQ** (`@contig:l-r:LEFT|RIGHT:...`); the genotyping contract is
  `>`-headed. `grep -c '^>'` on a discovery file returns **0** and looks exactly like a
  failed run. Count loci with
  `grep -oE '^@[^:]+:[0-9]+-[0-9]+' | sort -u | wc -l`.
- **LSF spools stdout** — `logs/*.out` is empty until the job *completes*. Use output-file
  counts (`ls discovery/*.txt.gz | wc -l`) as the live progress meter, not `tail`.
- **Discovery needs no BAM index** (it streams); **genotyping does** (per-locus region
  fetch). `pipeline.sh` indexes on demand.
- **Sidecars must ride the atomic rename.** The binary writes `{out}.splice.tsv` /
  `{out}.stats.json` next to `--out`, which is the `.tmp` path — they must be moved to the
  final name or Feature B output is silently orphaned.
- **`samtools quickcheck`** on every BAM before a run: it verifies the EOF block, which is
  exactly the truncation failure mode. All 174 PD44579 BAMs passed.
- Re-running any array is safe — every step skips work whose output already exists.

## Scheduler note

Scripts here use LSF (`bsub`) for farm22. The per-file work scripts
(`discover_one.sh`) are scheduler-agnostic; only the submit wrappers are LSF-specific.
The SLURM equivalent is in the header of `submit_discovery.sh`.
