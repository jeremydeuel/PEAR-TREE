#!/usr/bin/env python3
"""Genotyping-time extra evidence (peartree-genotype2 `gt_extra_reads`, rust/peartree-genotype2/src/extra.rs).

Two small steps around the per-colony genotyping:

  members  <P>.insertions.reads.fa.gz -> <P>.members.tsv.gz   (`#colonies N`, then
           `locus<TAB>sample,sample,...`)
           The colonies combine kept reads from per locus = the discovery members. Built once per
           patient (cluster/pipeline.sh, combine phase) so every genotype job reads a small table
           instead of the full reads FASTA (`--members`). `--n-colonies N` = ALL colonies of the run
           (pipeline: the discovery files), the denominator of the germline skip: a locus
           discovered in > `max_member_frac` (genotype2 `gt_extra_max_member_frac`, 0.5) of them
           gets no extra pass.

  merge    genotypes/<colony>.txt.gz.extra_reads.fa.gz + <P>.genotypes.csv.gz
           -> insertions/<P>.insertions.genotype_reads.fa.gz
           Keeps a colony's GT_* reads of a locus ONLY when the joint step calls that colony a
           carrier of it (matrix P(carrier) >= genotype2_io.P_CARRIER). The colony is the sidecar's
           file stem (= the matrix column); records are copied unchanged, colonies in sorted order.
           With `--members`, germline loci (as above) are dropped too, so sidecars written before
           the germline skip obey it.

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
MAX_MEMBER_FRAC = 0.5      # = genotype2 gt_extra_max_member_frac


def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def _locus_sample(header):
    """(locus, sample) of a reads.fa header `locus|SIDE|ROLE|sample|frag|r12` (without '>'); a
    locus name containing '|' keeps the last five fields as the tail."""
    parts = header.rstrip("\n").split("|")
    if len(parts) < 6:
        return None, None
    return "|".join(parts[:-5]), parts[-3]


def build_members(reads_fa, out, n_colonies=None):
    """locus -> samples with reads in combine's reads FASTA; written as a TSV, after a
    `#colonies N` line when the run's colony count is given."""
    members = {}
    with _open(reads_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                locus, sample = _locus_sample(line[1:])
                if locus is not None:
                    members.setdefault(locus, set()).add(sample)
    tmp = f"{out}.tmp.{os.getpid()}"
    with (gzip.open(tmp, "wt") if out.endswith(".gz") else open(tmp, "w")) as o:
        if n_colonies is not None:
            o.write(f"#colonies {int(n_colonies)}\n")
        o.write("locus\tsamples\n")
        for locus, ss in members.items():
            o.write(f"{locus}\t{','.join(sorted(ss))}\n")
    os.replace(tmp, out)
    return members


def is_germline(n_members, n_colonies, max_frac=MAX_MEMBER_FRAC):
    """n_members / n_colonies > max_frac (strict; = extra.rs is_germline)."""
    return n_colonies > 0 and n_members > max_frac * n_colonies


def germline_loci(members_tsv, max_frac=MAX_MEMBER_FRAC):
    """Loci of a members table discovered in > max_frac of its `#colonies`; None without that line."""
    n, counts = None, {}
    with _open(members_tsv) as fh:
        for line in fh:
            if line.startswith("#colonies"):
                n = int(line.split()[1])
            elif line.startswith("#") or line.startswith("locus\t") or "\t" not in line:
                continue
            else:
                locus, ss = line.rstrip("\n").split("\t", 1)
                counts[locus] = len([x for x in ss.split(",") if x])
    if n is None:
        return None
    return {loc for loc, k in counts.items() if is_germline(k, n, max_frac)}


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


def merge(genotype_dir, matrix, out, members=None, max_frac=MAX_MEMBER_FRAC):
    """Write the carrier-filtered union of the sidecars (minus germline loci when `members` has
    a colony count); returns (kept, dropped non-carrier, dropped germline, colonies)."""
    carriers = carriers_from_matrix(matrix)
    germline = (germline_loci(members, max_frac) if members and os.path.exists(members) else None) or set()
    files = sidecars(genotype_dir)
    kept = dropped = dropped_germ = 0
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
                    if locus in germline:
                        dropped_germ += 1
                    elif locus is not None and colony in carriers.get(locus, ()):
                        o.write(head)
                        o.write(line if line.endswith("\n") else line + "\n")
                        kept += 1
                    else:
                        dropped += 1
                    head = None
    os.replace(tmp, out)
    return kept, dropped, dropped_germ, len(files)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("members", help="discovery members per locus from combine's reads FASTA")
    m.add_argument("--reads-fa", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--n-colonies", type=int, help="ALL colonies of the run (germline-skip denominator)")
    g = sub.add_parser("merge", help="merge the per-colony sidecars, joint carriers only")
    g.add_argument("--genotype-dir", required=True)
    g.add_argument("--matrix", required=True, help="numeric joint matrix <P>.genotypes.csv.gz")
    g.add_argument("--out", required=True, help="<P>.insertions.genotype_reads.fa.gz")
    g.add_argument("--members", help="<P>.members.tsv.gz: drop germline loci (needs its #colonies line)")
    g.add_argument("--max-member-frac", type=float, default=MAX_MEMBER_FRAC,
                   help="germline: discovered in more than this fraction of the colonies (default 0.5)")
    a = ap.parse_args(argv)
    if a.cmd == "members":
        mem = build_members(a.reads_fa, a.out, a.n_colonies)
        ng = sum(1 for ss in mem.values() if a.n_colonies and is_germline(len(ss), a.n_colonies))
        print(f"members: {len(mem)} loci, {sum(map(len, mem.values()))} (locus, colony) pairs, "
              f"{a.n_colonies if a.n_colonies is not None else '?'} colonies ({ng} loci germline at > "
              f"{MAX_MEMBER_FRAC}) -> {a.out}")
    else:
        kept, dropped, dg, n = merge(a.genotype_dir, a.matrix, a.out, a.members, a.max_member_frac)
        print(f"genotype reads: {kept} kept (joint carriers, P >= {P_CARRIER}), {dropped} dropped as non-carrier, "
              f"{dg} dropped at germline loci, from {n} colony sidecar(s) -> {a.out}")


if __name__ == "__main__":
    main()
