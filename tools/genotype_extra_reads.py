#!/usr/bin/env python3
"""Genotyping-time extra evidence (peartree-genotype2 `gt_extra_reads`, rust/peartree-genotype2/src/extra.rs).

Two small steps around the per-colony genotyping:

  members  <P>.insertions.reads.fa.gz -> <P>.members.tsv.gz   (`locus<TAB>sample,sample,...`)
           The colonies combine kept reads from per locus = the discovery members. Built once per
           patient (cluster/pipeline.sh, combine phase) so every genotype job reads a small table
           instead of the full reads FASTA (`--members`).

  merge    genotypes/<colony>.txt.gz.extra_reads.fa.gz + <P>.genotypes.csv.gz
           -> insertions/<P>.insertions.genotype_reads.fa.gz
           Keeps a colony's GT_* reads of a locus ONLY when the joint step calls that colony a
           carrier of it (matrix P(carrier) >= genotype2_io.P_CARRIER). The colony is the sidecar's
           file stem (= the matrix column); records are copied unchanged, colonies in sorted order.

The GT_* reads are classification evidence for annotate_v2 (tools/rte), never junction evidence:
cluster/somatic_table.py's hard rules read only insertions.reads.fa.gz and its combine roles.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genotype2_io import P_CARRIER, read_matrix, colony_stem  # noqa: E402

SIDECAR_SUFFIX = ".extra_reads.fa.gz"


def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def _locus_sample(header):
    """(locus, sample) of a reads.fa header `locus|SIDE|ROLE|sample|frag|r12` (without '>'); a
    locus name containing '|' keeps the last five fields as the tail."""
    parts = header.rstrip("\n").split("|")
    if len(parts) < 6:
        return None, None
    return "|".join(parts[:-5]), parts[-3]


def build_members(reads_fa, out):
    """locus -> sorted samples with reads in combine's reads FASTA; written as a TSV."""
    members = {}
    with _open(reads_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                locus, sample = _locus_sample(line[1:])
                if locus is not None:
                    members.setdefault(locus, set()).add(sample)
    tmp = f"{out}.tmp.{os.getpid()}"
    with (gzip.open(tmp, "wt") if out.endswith(".gz") else open(tmp, "w")) as o:
        o.write("locus\tsamples\n")
        for locus, ss in members.items():
            o.write(f"{locus}\t{','.join(sorted(ss))}\n")
    os.replace(tmp, out)
    return members


def carriers_from_matrix(matrix):
    """{locus: set(colonies with P(carrier) >= P_CARRIER)} of the numeric joint matrix."""
    fmt, colonies, rows = read_matrix(matrix)
    if fmt is None:
        raise SystemExit(f"no matrix {matrix}")
    if fmt != "numeric":
        raise SystemExit(f"{matrix} is not a numeric genotype2 matrix (format {fmt})")
    out = {}
    for locus, cells in rows.items():
        cs = set()
        for c, v in zip(colonies, cells):
            try:
                if float(v) >= P_CARRIER:
                    cs.add(c)
            except ValueError:
                pass
        if cs:
            out[locus] = cs
    return out


def sidecars(genotype_dir):
    """{colony: sidecar path} of `<colony>.txt.gz.extra_reads.fa.gz` files (temporaries skipped)."""
    out = {}
    for f in sorted(os.listdir(genotype_dir)):
        if f.endswith(SIDECAR_SUFFIX) and ".tmp" not in f:
            out[colony_stem(f[: -len(SIDECAR_SUFFIX)])] = os.path.join(genotype_dir, f)
    return out


def merge(genotype_dir, matrix, out):
    """Write the carrier-filtered union of the sidecars; returns (kept, dropped, colonies)."""
    carriers = carriers_from_matrix(matrix)
    files = sidecars(genotype_dir)
    kept = dropped = 0
    tmp = f"{out}.tmp.{os.getpid()}"
    with gzip.open(tmp, "wt") as o:
        for colony, path in files.items():
            with _open(path) as fh:
                head = None
                for line in fh:
                    if line.startswith(">"):
                        head = line
                        continue
                    if head is None:
                        continue
                    locus, _sample = _locus_sample(head[1:])
                    if locus is not None and colony in carriers.get(locus, ()):
                        o.write(head)
                        o.write(line if line.endswith("\n") else line + "\n")
                        kept += 1
                    else:
                        dropped += 1
                    head = None
    os.replace(tmp, out)
    return kept, dropped, len(files)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("members", help="discovery members per locus from combine's reads FASTA")
    m.add_argument("--reads-fa", required=True)
    m.add_argument("--out", required=True)
    g = sub.add_parser("merge", help="merge the per-colony sidecars, joint carriers only")
    g.add_argument("--genotype-dir", required=True)
    g.add_argument("--matrix", required=True, help="numeric joint matrix <P>.genotypes.csv.gz")
    g.add_argument("--out", required=True, help="<P>.insertions.genotype_reads.fa.gz")
    a = ap.parse_args(argv)
    if a.cmd == "members":
        mem = build_members(a.reads_fa, a.out)
        print(f"members: {len(mem)} loci, {sum(map(len, mem.values()))} (locus, colony) pairs -> {a.out}")
    else:
        kept, dropped, n = merge(a.genotype_dir, a.matrix, a.out)
        print(f"genotype reads: {kept} kept (joint carriers, P >= {P_CARRIER}), {dropped} dropped, "
              f"from {n} colony sidecar(s) -> {a.out}")


if __name__ == "__main__":
    main()
