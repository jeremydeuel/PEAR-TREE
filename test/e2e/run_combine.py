#!/usr/bin/env python3
"""E2E helper: run combine_insertions with the generated config (test/e2e/make_config.py)."""
import argparse
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config-dir", required=True)
    ap.add_argument("--out-stem", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    sys.path.insert(0, os.path.abspath(a.config_dir))
    sys.path.insert(1, os.path.join(REPO, "src"))
    import config  # noqa: F401  (the generated one, first on sys.path)
    assert os.path.abspath(config.__file__).startswith(os.path.abspath(a.config_dir)), config.__file__
    from combine_insertions import combine_insertions
    s = a.out_stem
    combine_insertions(a.files, f"{s}.genotyping.txt.gz", f"{s}.combined.txt.gz", f"{s}.fq.gz",
                       f"{s}.bam", a.threads)


if __name__ == "__main__":
    main()
