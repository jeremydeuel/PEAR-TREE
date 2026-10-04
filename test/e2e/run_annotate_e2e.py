#!/usr/bin/env python3
"""E2E helper: annotate_v2 (+ tools/rte) on the combined insertions.

Genotyping is not part of this test: a stub genotypes file marks every combined insertion
heterozygous in every colony, so annotate_v2 keeps all of them (it drops loci with no
het/hom tip). Everything else -- Dfam scan (nhmmscan, local family HMM), bowtie2 clip remap,
rmsk, the tools/rte pass with resources/rte_library and the combine sidecars -- is the real
annotate_v2 code path."""
import argparse
import gzip
import importlib.util
import os
import sys
import types

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config-dir", required=True)
    ap.add_argument("--combined", required=True)
    ap.add_argument("--samples", type=int, required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--table", required=True)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    sys.path.insert(0, os.path.abspath(a.config_dir))
    sys.path.insert(1, REPO)
    import config
    import src
    sys.modules["src.config"] = config          # annotate_v2 does `from src.config import CONFIG`
    src.config = config
    titles = []
    with gzip.open(a.combined, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0 and line.startswith("@"):
                t = line[1:].strip()[:-2]
                if t not in titles:
                    titles.append(t)
    sample = os.path.basename(a.combined).split(".")[0]
    gt = os.path.join(a.workdir, f"{sample}.genotypes.csv.gz")
    with gzip.open(gt, "wt") as fh:
        fh.write(";".join(["insertion"] + [f"S{i + 1}" for i in range(a.samples)]) + "\n")
        for t in titles:
            fh.write(";".join([t] + ["heterozygous"] * a.samples) + "\n")
    A = config.CONFIG["annotate"]
    A["insertions_file"] = lambda s: a.combined
    A["genotyping_file"] = lambda s: gt
    spec = importlib.util.spec_from_file_location("annotate_v2", os.path.join(REPO, "tools", "annotate_v2.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    f = m.VariantAnnotationContainer(sample, a.table)
    f.write_table(a.table)
    with open(os.path.join(a.workdir, f"{sample}.annotate_report.txt"), "w") as out:
        old = sys.stdout
        sys.stdout = out
        try:
            f.print()
        finally:
            sys.stdout = old


if __name__ == "__main__":
    main()
