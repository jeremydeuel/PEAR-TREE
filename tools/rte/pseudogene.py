"""Processed-pseudogene proof: an exon-exon junction covered by a read.

annotate_v2's `_pseudogene()` flags a candidate when junction clips land in exons of one gene;
that is not proof (a clip inside an exon is also what a genomic fragment / mis-mapped read
gives). A retrotransposed spliced mRNA must contain exon i joined DIRECTLY to exon i+1 (intron
removed). We build, for each candidate gene, the junction "cores" -- the last `overhang` bp of
exon i + the first `overhang` bp of exon j (j = i+1, and i+2 for an exon-skipping isoform) from
the remap genome -- and search them (both strands, edlib, <= `max_edits`) in every clipped read,
mate and junction consensus of the insertion. A hit means the read spans the splice junction
with >= overhang bp on both sides => tag EXON_JUNCTION. Without it the call stays
`PSEUDOGENE_CANDIDATE` (tag only, element is not PSEUDOGENE).

Exon track: contig, start, end, gene[, strand] (0-based half-open; plain or gzipped) -- the same
file as CONFIG['annotate']['exon_annotation'] (remap-genome coordinates), or a GeneModel track.
"""
from __future__ import annotations

import gzip

from .sequtil import edlib_best

DEFAULTS = {"exon_junction_overhang": 20, "exon_junction_max_edits": 2,
            "exon_junction_min_intron": 30, "exon_junction_skip": 1}


def load_exons_by_gene(path):
    genes = {}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            genes.setdefault(f[3], []).append((f[0], int(f[1]), int(f[2])))
    for g in genes:
        genes[g] = _merge(sorted(genes[g]))
    return genes


def _merge(ivs):
    out = []
    for c, s, e in ivs:
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1] = (c, out[-1][1], max(out[-1][2], e))
        else:
            out.append((c, s, e))
    return out


class ExonJunctionIndex:
    def __init__(self, exons_by_gene: dict, genome, cfg=None):
        self.cfg = dict(DEFAULTS)
        self.cfg.update(cfg or {})
        self.exons = exons_by_gene
        self.genome = genome
        self._cores = {}

    def cores(self, gene):
        if gene in self._cores:
            return self._cores[gene]
        ov = self.cfg["exon_junction_overhang"]
        out = []
        ex = self.exons.get(gene, [])
        if self.genome is not None:
            for i in range(len(ex)):
                for j in range(i + 1, min(len(ex), i + 2 + self.cfg["exon_junction_skip"])):
                    (ci, si, ei), (cj, sj, ej) = ex[i], ex[j]
                    if ci != cj or sj - ei < self.cfg["exon_junction_min_intron"]:
                        continue
                    if ei - si < ov or ej - sj < ov:
                        continue
                    core = self.genome.fetch(ci, ei - ov, ei) + self.genome.fetch(cj, sj, sj + ov)
                    if len(core) == 2 * ov and "N" not in core:
                        out.append((f"{gene}:e{i + 1}-e{j + 1}", core))
        self._cores[gene] = out
        return out

    def find(self, genes, seqs):
        """genes: candidate gene ids; seqs: [(name, sequence)]. Returns [(junction_label, name)]."""
        k_frac = self.cfg["exon_junction_max_edits"] / (2 * self.cfg["exon_junction_overhang"])
        hits = []
        for g in genes:
            for label, core in self.cores(g):
                for name, s in seqs:
                    if len(s) < len(core):
                        continue
                    r = edlib_best(core, s, max_frac=k_frac + 1e-9)
                    if r is not None:
                        hits.append((label, name))
        return hits
