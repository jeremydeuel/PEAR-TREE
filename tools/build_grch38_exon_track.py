#!/usr/bin/env python3
"""Build the GRCh38 exon track that discovery Feature B (splice_hallmark) needs.

Input : an Ensembl GRCh38 GTF (numeric contigs 1..22,X,Y; 1-based inclusive coords).
Output: `contig begin end gene` rows, chr-prefixed (chr1..chrY), 0-based half-open,
        per-gene merged exons, sorted, gzipped -- the exact format the Rust reader
        (rust/peartree-discovery/src/exons.rs) and annotate_v2 consume.

Only protein_coding + lncRNA genes are kept by default (the biotypes a processed
pseudogene's parent transcript realistically belongs to); pass --biotypes to change.
Contigs are restricted to the discovery allowlist (chr1..chr22,chrX,chrY); MT, alts,
and scaffolds are dropped -- discovery only ever looks up mate landings on the primary
assembly.

    tools/build_grch38_exon_track.py Homo_sapiens.GRCh38.112.gtf.gz grch38.exons.bed.gz

The gene label is gene_name when present, else gene_id. Exons are merged within each
(gene, contig) so the track is non-overlapping per gene, which is what exons.rs's
lookup assumes.
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
    ap.add_argument("gtf", help="Ensembl GRCh38 GTF (.gz ok)")
    ap.add_argument("out", help="output exon track (.gz ok)")
    ap.add_argument("--biotypes", default="protein_coding,lncRNA",
                    help="comma list of gene_biotypes to keep (default: protein_coding,lncRNA; "
                         "'all' keeps every biotype)")
    a = ap.parse_args()
    keep = None if a.biotypes == "all" else set(a.biotypes.split(","))

    spans = defaultdict(list)   # (gene, chrom) -> [(begin0, end)]
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
            spans[(gene, "chr" + f[0])].append((int(f[3]) - 1, int(f[4])))
            n_exon += 1

    rows = []
    for (gene, chrom), ivs in spans.items():
        ivs.sort()
        cb, ce = ivs[0]
        for b, e in ivs[1:]:
            if b <= ce:
                ce = max(ce, e)
            else:
                rows.append((chrom, cb, ce, gene)); cb, ce = b, e
        rows.append((chrom, cb, ce, gene))

    def chrkey(c):
        s = c[3:]
        return int(s) if s.isdigit() else (100 if s == "X" else 101)
    rows.sort(key=lambda r: (chrkey(r[0]), r[1], r[2]))

    out_opener = gzip.open if a.out.endswith(".gz") else open
    with out_opener(a.out, "wt") as o:
        for chrom, b, e, gene in rows:
            o.write(f"{chrom}\t{b}\t{e}\t{gene}\n")

    n_genes = len({g for g, _ in spans})
    print(f"kept {n_exon} exons ({a.biotypes}); {n_genes} genes; {len(rows)} merged rows -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
