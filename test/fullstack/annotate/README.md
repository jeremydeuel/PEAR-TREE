# Annotation resources + local runner (step 3 of the pipeline)

Runs `tools/annotate_v2.py` on a workstation / the test harness, and holds the reference
resources it needs. Annotation classifies each combined insertion's inserted sequence into an
RTE family (L1 / Alu / SVA / HERVK) or a processed pseudogene, from three signals:

1. **Family HMM scan** — `nhmmscan --dfamtblout` of the inserted clips against a nucleotide
   HMM library. On the cluster this is the full Dfam set via `dfamscan.pl`; here it is a small
   per-family library built from the exact hs1 source loci the donor samples.
2. **Clip remap + RepeatMasker** — bowtie2 of the clips to hs1, then the hs1 RepeatMasker track
   names the repeat at the landing site. Subfamily-resolving and DFAM-independent.
3. **Pseudogene exon track** — a processed pseudogene is a spliced mRNA, so its clips map into a
   parent gene's exons. Clips hitting ≥2 distinct exons of one gene (or 1 exon + a poly-A tail)
   are called a pseudogene of that gene.

## Files

| file | what | build |
|---|---|---|
| `rte_elements.fa` | the 5 implanted RTE source sequences (hs1 loci) | committed |
| `peartree_rte.hmm(.h3*)` | hmmpress'd family HMM library (config `annotate.hmm`) | `build_hmm.sh` |
| `pseudogene_exons.hs1.bed` | hs1 exon intervals of the pseudogene parent genes (config `annotate.exon_annotation`) | `build_exon_track.py` |
| `run_annotate.py` | local runner (patches `CONFIG['annotate']` to these paths + nhmmscan back-end) | — |

The `.hmm*` binaries are **not committed** (rebuild them once):

```bash
# family HMM library — needs only HMMER (rebuilds from the committed rte_elements.fa)
test/fullstack/annotate/build_hmm.sh
# to re-extract the source sequences from hs1 (e.g. to add a family):
test/fullstack/annotate/build_hmm.sh --from-hs1 ~/Downloads/hs1.fa

# pseudogene exon track — tiles the parent mRNAs through bowtie2/hs1 (same mapper the clips use)
./venv/bin/python test/fullstack/annotate/build_exon_track.py \
    --pseudogene-fasta test/genotyping/pseudogenes.fa \
    --bowtie2 "$(command -v bowtie2)" --index ~/Downloads/hs1 \
    --out test/fullstack/annotate/pseudogene_exons.hs1.bed
```

## Run

`run_scale10k.sh` calls this as step 9 automatically once `peartree_rte.hmm` exists. Standalone:

```bash
HS1_BT2=~/Downloads/hs1 ./venv/bin/python test/fullstack/annotate/run_annotate.py \
    --combined step2.combined.txt.gz --genotypes step2.genotypes.csv.gz \
    --out annotate.txt --workdir annot/
```

## Validated recovery (1k genotyping-smoke harness, 108 combined calls)

| class | correct | notes |
|---|---|---|
| L1HS | 32/33 | subfamily drifts (L1HS↔L1PAx) via RMSK-remap to paralogs — expected |
| AluYa5 | 38/39 | |
| HERVK | 6/6 | |
| SVA (E+F) | 10/11 | **needs the family HMM** — RMSK-remap alone lands SVA clips on L1/LTR |
| MALAT1 (pseudogene) | 4/4 | single-exon lncRNA → 1 exon + poly-A |
| HNRNPA1 / RPL21 / CASP12 / DUX4 | 5/15 | multi-copy parents; terminal clips multi-map or don't map |

**0 false pseudogene calls.** The remaining multi-copy pseudogene misses (clip unmapped or a
lone exon) are what **Feature B mate-splice** targets — mates spanning ≥2 exons of one gene,
independent of the terminal clip. It is implemented in Rust discovery (`splice_hallmark` +
`exon_annotation`, in the read-mapping genome = hg38) and surfaced by `annotate_v2.read_splice`;
enable it in discovery and rerun to recover those. Not exercised here because the smoke run's
discovery had it off and the source BAM is not retained.
