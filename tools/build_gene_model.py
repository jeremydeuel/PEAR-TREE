#!/usr/bin/env python3
"""Build the gene-model track that annotate_v2's insertion-SITE annotator (GeneModel) consumes.

Input : an Ensembl GTF (numeric contigs 1..22,X,Y; 1-based inclusive coords).
Output: `contig<TAB>start<TAB>end<TAB>gene<TAB>strand` rows, 0-based half-open, per-gene merged
        exons, sorted -- the exact format tools/annotate_v2.py:GeneModel reads. gzipped when the
        output path ends in .gz.

This differs from build_grch38_exon_track.py (the discovery / processed-pseudogene exon track)
in two ways that the site annotator needs: it carries the gene STRAND (donor vs acceptor and the
TSS are strand-dependent), and it is a track of the SAMPLE's DISCOVERY/BAM genome (where the
insertion locus lives), not of bowtie2_index2 (the clip-remap genome). GeneModel reconstructs the
gene span, the TSS and every intron/exon boundary from these merged-exon rows, so this one file
drives exon / splice-site / intron / promoter classification.

Contig naming must match the title-locus contigs. GeneModel's lookup is chr-prefix tolerant
(chr8 <-> 8), so either convention works, but pick the one your discovery reference uses:
    # GRCh38 discovery (chr-prefixed contigs):
    tools/build_gene_model.py Homo_sapiens.GRCh38.112.gtf.gz grch38.gene_model.tsv.gz
    # GRCh37/hs37d5 discovery (numeric contigs):
    tools/build_gene_model.py --no-chr Homo_sapiens.GRCh37.87.gtf.gz grch37.gene_model.tsv.gz

Only protein_coding + lncRNA genes are kept by default (the biotypes whose disruption is
interpretable); pass --biotypes to change. Contigs are restricted to 1..22,X,Y (MT, alts and
scaffolds dropped -- insertion loci are called on the primary assembly). The gene label is
gene_name when present, else gene_id; exons are merged within each (gene, contig).
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict

CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}
GID = re.compile(r'gene_id "([^"]+)"')
GNAME = re.compile(r'gene_name "([^"]+)"')
GBIO = re.compile(r'gene_biotype "([^"]+)"')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gtf", help="Ensembl GTF (.gz ok)")
    ap.add_argument("out", help="output gene-model track (.gz ok)")
    ap.add_argument("--biotypes", default="protein_coding,lncRNA",
                    help="comma list of gene_biotypes to keep (default: protein_coding,lncRNA; "
                         "'all' keeps every biotype)")
    ap.add_argument("--no-chr", action="store_true",
                    help="emit numeric contigs (1..Y) instead of chr-prefixed (chr1..chrY); "
                         "use for a GRCh37/hs37d5 discovery reference")
    a = ap.parse_args()
    keep = None if a.biotypes == "all" else set(a.biotypes.split(","))

    # (gene, contig) -> [strand, [(begin0, end), ...]]
    spans = defaultdict(lambda: [None, []])
    n_exon = 0
    opener = gzip.open if a.gtf.endswith(".gz") else open
    with opener(a.gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon" or f[0] not in CHROMS:
                continue
            attr = f[8]
            if keep is not None:
                mb = GBIO.search(attr)
                if not mb or mb.group(1) not in keep:
                    continue
            gn = GNAME.search(attr)
            gene = gn.group(1) if gn else GID.search(attr).group(1)
            contig = f[0] if a.no_chr else "chr" + f[0]
            strand = f[6] if f[6] in ("+", "-") else "+"
            rec = spans[(gene, contig)]
            rec[0] = strand
            rec[1].append((int(f[3]) - 1, int(f[4])))
            n_exon += 1

    rows = []
    for (gene, contig), (strand, ivs) in spans.items():
        ivs.sort()
        cb, ce = ivs[0]
        for b, e in ivs[1:]:
            if b <= ce:
                ce = max(ce, e)
            else:
                rows.append((contig, cb, ce, gene, strand)); cb, ce = b, e
        rows.append((contig, cb, ce, gene, strand))

    def chrkey(c):
        s = c[3:] if c.startswith("chr") else c
        return int(s) if s.isdigit() else (100 if s == "X" else 101)
    rows.sort(key=lambda r: (chrkey(r[0]), r[1], r[2]))

    out_opener = gzip.open if a.out.endswith(".gz") else open
    with out_opener(a.out, "wt") as o:
        for contig, b, e, gene, strand in rows:
            o.write(f"{contig}\t{b}\t{e}\t{gene}\t{strand}\n")

    n_genes = len({g for g, _ in spans})
    print(f"kept {n_exon} exons ({a.biotypes}); {n_genes} genes; {len(rows)} merged rows -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
