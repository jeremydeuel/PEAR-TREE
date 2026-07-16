#!/usr/bin/env python3
"""Run the genotype / combine_genotypes step with a local config shim.

Running `python src/main.py ...` puts src/ first on sys.path (script-dir precedence),
which shadows the panel $OUT/config.py with the cluster src/config.py. So, like
run_combine.py, we insert the config-dir FIRST and call the pipeline function directly,
never importing main.py.

The `if __name__ == "__main__"` guard is REQUIRED: collect_genotype() uses a
multiprocessing Pool, and under macOS spawn each worker re-imports this module -- without
the guard that re-runs the step and spawns workers recursively (hang).

Usage:
  run_step.py <config_dir> <repo_src> genotype <bam> <out> <insertions> <threads>
  run_step.py <config_dir> <repo_src> combine_genotypes <out> <threads> <gt1> <gt2> ...
"""
import sys


def main():
    config_dir, repo_src, step = sys.argv[1], sys.argv[2], sys.argv[3]
    sys.path.insert(0, config_dir)   # panel config.py wins for `from config import CONFIG`
    sys.path.append(repo_src)        # src/ pipeline modules
    rest = sys.argv[4:]

    if step == "genotype":
        bam, out, insertions, threads = rest[0], rest[1], rest[2], int(rest[3])
        from genotype import genotype
        genotype(insertions, bam, out, threads)
    elif step == "combine_genotypes":
        out, threads = rest[0], int(rest[1])
        gt_files = rest[2:]
        from combine_genotypes import collect_genotype
        collect_genotype(gt_files, out, threads)
    else:
        sys.exit(f"unknown step {step!r}")


if __name__ == "__main__":
    main()
