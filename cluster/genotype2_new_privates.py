#!/usr/bin/env python3
"""List loci the joint step calls private in one run but not in another (e.g. --ref-bias auto vs
off), with the carrier's per-colony reads and the alt reads in every other colony.

  python cluster/genotype2_new_privates.py --base base.joint.tsv --new bias.joint.tsv \
      --genotype-dir V2_refbias/genotypes [--fit V2_refbias/fit/phylo_fit.tsv]
"""
import argparse
import collections
import csv
import gzip
import os


def rows(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def joint_class(r):
    b = r["best"]
    if b in ("ROOT", "NOISE", "INDEP"):
        return b
    return "clade" if int(r["n_carriers"] or 0) >= 2 else "private"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="joint table of the reference run")
    ap.add_argument("--new", required=True, help="joint table of the run whose extra privates are listed")
    ap.add_argument("--genotype-dir", required=True, help="per-colony v2 files <colony>.txt.gz")
    ap.add_argument("--fit", help="tree_fit phylo_fit.tsv (adds its class column)")
    a = ap.parse_args()

    base = {r["locus"]: r for r in rows(a.base)}
    new_rows = rows(a.new)
    tf = {r["locus"]: r["class"] for r in rows(a.fit)} if a.fit else {}
    extra = [r for r in new_rows if joint_class(r) == "private" and joint_class(base[r["locus"]]) != "private"]
    loci = {r["locus"] for r in extra}
    colonies = [k[2:] for k in new_rows[0] if k.startswith("p_")] if new_rows else []
    gt = {}
    for c in colonies:
        p = os.path.join(a.genotype_dir, c + ".txt.gz")
        if os.path.exists(p):
            gt[c] = {g["locus"]: g for g in rows(p) if g["locus"] in loci}

    hdr = ["locus", "kind", "tree_fit", "base_class", "carrier", "post_best", "p_carrier", "depth", "n_alt",
           "n_ref", "n_uninf", "vaf", "other_colonies_with_alt"]
    print("\t".join(hdr))
    for r in sorted(extra, key=lambda r: r["locus"]):
        loc, c = r["locus"], r["carriers"]
        g = gt.get(c, {}).get(loc, {})
        others = [f"{o}:{gt[o][loc]['n_alt']}" for o in gt if o != c and loc in gt[o] and int(gt[o][loc]["n_alt"] or 0) > 0]
        print("\t".join([loc, g.get("kind", "NA"), tf.get(loc, "NA"), joint_class(base[loc]), c, r["post_best"],
                         r.get("p_" + c, ""), g.get("depth", ""), g.get("n_alt", ""), g.get("n_ref", ""),
                         g.get("n_uninf", ""), g.get("vaf", ""), ",".join(others) or "-"]))
    print(f"# n = {len(extra)}; by tree_fit class: {dict(collections.Counter(tf.get(r['locus'], 'NA') for r in extra))}; "
          f"base class: {dict(collections.Counter(joint_class(base[r['locus']]) for r in extra))}; "
          f"carrier alt reads: {dict(sorted(collections.Counter(gt.get(r['carriers'], {}).get(r['locus'], {}).get('n_alt', 'NA') for r in extra).items()))}")


if __name__ == "__main__":
    main()
