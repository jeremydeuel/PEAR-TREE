#!/usr/bin/env python3
"""Run PEAR-TREE step 2 (combine_insertions) with a local config shim.

The combine chain does `from config import CONFIG`; we shadow src/config.py by putting
the config-dir (which holds a materialised config.py — see config.local.py.template)
first on sys.path, then src/ for the pipeline modules.

Usage: run_combine.py <disc.txt.gz> <out_stem> <threads> <config_dir> <repo_src>
"""
import sys

disc, stem, threads, config_dir, repo_src = (
    sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])
sys.path.insert(0, config_dir)   # holds config.py (local paths) -> wins for `import config`
sys.path.append(repo_src)        # src/ pipeline modules

from combine_insertions import combine_insertions
combine_insertions([disc],
                   f"{stem}.genotyping.txt.gz",
                   f"{stem}.combined.txt.gz",
                   f"{stem}.fq.gz",
                   f"{stem}.bam",
                   threads)
print("combine done")
