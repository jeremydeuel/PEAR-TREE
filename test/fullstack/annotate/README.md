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
| `pseudogene_exons.hs1.bed` | hs1 exon intervals of the pseudogene parents — **annotate** clip track (config `annotate.exon_annotation`) | `build_exon_track.py` |
| `pseudogene_exons.hg38.bed` | hg38 exon intervals — **discovery** mate-splice track (Feature B, `exon_annotation`) | `build_exon_track.py --bwa` |
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

## Feature B: mate-splice pseudogene detection

Multi-copy pseudogene parents (HNRNPA1/CASP12/DUX4/RPL21) have terminal clips that multi-map
or don't map, so the clip path alone misses them. **Feature B** catches them from the mates
instead: discovery (`splice_hallmark = true` + `exon_annotation = <hg38 track>`) flags a
breakpoint whose mates span ≥ `splice_min_exons` exons of one gene with the introns skipped,
writes `<out>.splice.tsv`; combine re-keys it to `<combined>.splice.tsv`; `annotate_v2`
consumes it. `run_genotyping.sh` wires this automatically (builds the hg38 track, enables
splice, then annotates as step 11).

## Validated recovery (genotyping harness, het_60x, ~106 combined calls)

| class | clip path only | + Feature B mate-splice |
|---|---|---|
| L1HS | 32/33 | 32/33 |
| AluYa5 | 38/39 | 37/38 |
| HERVK | 6/6 | 6/6 |
| SVA (E+F) | 10/11 | 10/11 |
| MALAT1 pseudogene | 4/4 | 4/4 |
| HNRNPA1 / CASP12 / DUX4 / RPL21 | 5/15 | **9/14** |
| **overall** | **88.0%** | **92.5%** |

**0 false pseudogene calls** in both. SVA needs the family HMM (RMSK-remap alone lands SVA
clips on L1/LTR); L1/Alu subfamily labels drift (L1HS↔L1PAx) via RMSK-remap to paralogs —
expected. Remaining pseudogene misses are 5′-truncated copies whose mates don't reach a second
exon at 60×.
