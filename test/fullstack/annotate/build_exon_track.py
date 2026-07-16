#!/usr/bin/env python3
"""Build an hs1 exon BED for the processed-pseudogene parent genes, for annotate_v2's
clip-level pseudogene detection (CONFIG['annotate']['exon_annotation']).

The processed-pseudogene test implants mature mRNAs (test/genotyping/pseudogenes.fa). When
an inserted mRNA's terminal clip is remapped to hs1 it lands on the parent gene's exon; the
annotator calls a pseudogene when clips hit >= 2 distinct exons of one gene (or 1 exon +
polyA). That needs the parent exon intervals in hs1 coordinates.

We reconstruct them with the SAME aligner the clips use (bowtie2 end-to-end vs the hs1
index): tile each mRNA into overlapping windows, map the tiles, keep confident (MAPQ >=
MIN_MAPQ) placements, and merge tiles whose genomic gap is small into one block. A gap larger
than MAX_MERGE_GAP (an intron) starts a new block, so each block is one genomic exon. Output
is the exon-track BED format annotate/discovery expect: contig<TAB>start<TAB>end<TAB>gene
(0-based half-open), one row per exon.

The clip-remap track (annotate) is built against hs1; the discovery mate-splice track
(Feature B) is built against the read-mapping genome (hg38) where the mates land — same
format, different reference. Pick the aligner you have an index for.

Usage:
  # hs1 (annotate clip track), bowtie2 index:
  build_exon_track.py --pseudogene-fasta test/genotyping/pseudogenes.fa \
      --bowtie2 $(command -v bowtie2) --index ~/Downloads/hs1 --out pseudogene_exons.hs1.bed
  # hg38 (discovery mate-splice track), bwa index (the .fa with its bwa index alongside):
  build_exon_track.py --pseudogene-fasta test/genotyping/pseudogenes.fa \
      --bwa $(command -v bwa) --ref ~/Downloads/hg38.fa --out pseudogene_exons.hg38.bed
"""
import argparse, os, subprocess, sys, tempfile

WIN, STEP = 100, 25          # tile window / step (bp)
MAX_MERGE_GAP = 40           # genomic gap up to this joins tiles into one exon; larger = intron
MIN_EXON = 30                # drop blocks shorter than this (spurious single-tile noise)
# MIN_MAPQ default 20 = confident, ~unique placement (right for the hs1 clip track). For the
# hg38 mate-splice track, lower it (e.g. 0): a processed pseudogene of a multi-copy gene has
# reads that multi-map, so a MAPQ floor drops exactly those parents (RPL21/DUX4). Discovery
# captures the mate landing regardless of the mate's MAPQ, so an inclusive track is what lets
# the splice check see those genes.


def read_fasta(path):
    name, seq, out = None, [], {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name:
        out[name] = "".join(seq)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pseudogene-fasta", required=True)
    ap.add_argument("--bowtie2", help="bowtie2 executable (use with --index)")
    ap.add_argument("--index", help="bowtie2 index prefix")
    ap.add_argument("--bwa", help="bwa executable (use with --ref)")
    ap.add_argument("--ref", help="reference fasta with a bwa index alongside")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-mapq", type=int, default=20,
                    help="tile placement MAPQ floor (20 for the hs1 clip track; 0 for the hg38 mate-splice track)")
    a = ap.parse_args()
    if not ((a.bowtie2 and a.index) or (a.bwa and a.ref)):
        ap.error("give either --bowtie2 + --index or --bwa + --ref")

    genes = read_fasta(a.pseudogene_fasta)
    # tile every mRNA; tile id encodes gene + within-mRNA offset
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        tiles_fa = tf.name
        for gene, seq in genes.items():
            for i in range(0, max(1, len(seq) - WIN + 1), STEP):
                tf.write(f">{gene}|{i}\n{seq[i:i + WIN]}\n")

    # map tiles, best hit each, same aligner family the clips/mates use
    if a.bowtie2:
        cmd = [a.bowtie2, "-x", a.index, "-f", tiles_fa, "--end-to-end", "-k", "1",
               "-p", str(a.threads), "--quiet"]
    else:
        cmd = [a.bwa, "mem", "-t", str(a.threads), a.ref, tiles_fa]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.exit(1)
    os.remove(tiles_fa)

    # collect confident placements per (gene, contig): list of (start, end)
    placements = {}
    for line in proc.stdout.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        qname, flag, contig, pos, mapq = f[0], int(f[1]), f[2], int(f[3]), int(f[4])
        if flag & 0x4 or mapq < a.min_mapq or contig == "*":
            continue
        gene = qname.split("|")[0]
        start = pos - 1                                  # SAM 1-based -> 0-based
        end = start + WIN
        placements.setdefault((gene, contig), []).append((start, end))

    # merge tiles into exon blocks (gap <= MAX_MERGE_GAP joins; larger splits at introns)
    rows = []
    for (gene, contig), ivs in placements.items():
        ivs.sort()
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce + MAX_MERGE_GAP:
                ce = max(ce, e)
            else:
                if ce - cs >= MIN_EXON:
                    rows.append((contig, cs, ce, gene))
                cs, ce = s, e
        if ce - cs >= MIN_EXON:
            rows.append((contig, cs, ce, gene))

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(a.out, "w") as out:
        for contig, s, e, gene in rows:
            out.write(f"{contig}\t{s}\t{e}\t{gene}\n")
    by_gene = {}
    for _, _, _, g in rows:
        by_gene[g] = by_gene.get(g, 0) + 1
    print(f"wrote {len(rows)} exon blocks to {a.out}")
    for g in sorted(by_gene):
        print(f"  {g}: {by_gene[g]} exons")


if __name__ == "__main__":
    main()
