#!/usr/bin/env python3
"""Generate golden data for the Rust port of tools/rte (rust/peartree-rte).

Runs the CURRENT python tools/rte with golden/recorder.py installed and writes every module-
boundary call as JSON lines (schema: SPEC.md "Golden harness", parser: src/golden.rs).

  # 1. the tools/rte pytest suite (synthetic insertions; committed as golden/data/pytest.jsonl.gz)
  make_golden.py pytest --out rust/peartree-rte/golden/data/pytest.jsonl.gz

  # 2. the e2e_phylo simulation through annotate_v2 (real combine sidecars + genotype reads; big,
  #    kept out of git). Needs the e2e config (copied to src/config.py for the run and removed
  #    again) and writes, next to --out, `e2e.inputs.jsonl` + `e2e.config.json` +
  #    `e2e.manifest.json` (the paths the Rust end-to-end test feeds the binary).
  make_golden.py e2e --out SCRATCH/golden/e2e.jsonl.gz --config SCRATCH/e2e_annot_config.py \\
      --tmp SCRATCH/golden/tmp --gt SCRATCH/e2e/ins/P1.insertions.genotype_reads.fa.gz [--cache DIR]

Python: a venv with mappy, edlib, py2bit, pysam, pytest (the e2e needs annotate_v2's deps and
bowtie2 / hmmscan unless --cache holds annotate_v2's tmp files of an earlier identical run).
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import recorder  # noqa: E402

TESTS = ["test_rte_hallmarks.py", "test_rte_inverted5p.py", "test_rte_polya_support.py",
         "test_rte_round2.py", "test_rte_score.py", "test_rte_structure.py", "test_rte_td_sources.py",
         "test_rte_transduction.py", "test_rte_locus_class.py"]


def run_pytest(a):
    import pytest
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    recorder.install(a.out)
    sys.path.insert(0, os.path.join(REPO, "test"))
    rc = pytest.main([os.path.join(REPO, "test", t) for t in (a.tests or TESTS)]
                     + ["-q", "-p", "no:cacheprovider", "-x" if a.stop else "-q"])
    recorder.close()
    print(f"pytest exit {rc}; golden -> {a.out} ({os.path.getsize(a.out)} bytes)")
    return int(rc)


def run_e2e(a):
    out = os.path.abspath(a.out)
    d = os.path.dirname(out)
    os.makedirs(d, exist_ok=True)
    tmp = os.path.abspath(a.tmp)
    os.makedirs(tmp, exist_ok=True)
    if a.cache:
        for f in glob.glob(os.path.join(a.cache, "P1.*")):
            if not f.endswith(".tsv"):
                shutil.copy(f, tmp)
    cfg_dst = os.path.join(REPO, "src", "config.py")
    if os.path.exists(cfg_dst):
        sys.exit(f"{cfg_dst} exists -- refusing to overwrite a deployment config")
    os.environ.update({"PT_REPO": REPO, "PT_TMP": tmp, "PT_GT": a.gt or ""})
    shutil.copy(a.config, cfg_dst)
    try:
        sys.path.insert(0, os.path.join(REPO, "tools"))
        from tools.rte import annotator as A, rust_bridge
        from src.config import CONFIG
        orig_all = A.RteAnnotator.annotate_all
        side = {}

        def annotate_all(self, inputs):
            # the Rust binary's inputs, written exactly where annotate_v2 would hand them over
            rust_bridge.write_inputs(inputs, os.path.join(d, "e2e.inputs.jsonl"))
            rust_bridge.write_config(self.cfg, os.path.join(d, "e2e.config.json"))
            side["n"] = len(inputs)
            return orig_all(self, inputs)
        A.RteAnnotator.annotate_all = annotate_all
        recorder.install(out, case_fn=lambda: "e2e", full_reads=not a.no_reads)
        import runpy
        sys.argv = ["annotate_v2.py", "P1", os.path.join(d, "e2e.annotated.tsv")]
        runpy.run_path(os.path.join(REPO, "tools", "annotate_v2.py"), run_name="__main__")
        recorder.close()
        ins = CONFIG["annotate"]["insertions_file"]("P1")
        ev, rd = A.default_sidecars(ins)
        import json
        with open(os.path.join(d, "e2e.manifest.json"), "w") as fh:
            json.dump({"inputs": os.path.join(d, "e2e.inputs.jsonl"), "config": os.path.join(d, "e2e.config.json"),
                       "evidence": ev, "reads": rd, "gt_reads": a.gt or None, "golden": out,
                       "annotated": os.path.join(d, "e2e.annotated.tsv"), "n_inputs": side.get("n")}, fh, indent=1)
    finally:
        os.remove(cfg_dst)
    print(f"e2e golden -> {out} ({os.path.getsize(out)} bytes)")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pytest")
    p.add_argument("--out", required=True)
    p.add_argument("--tests", nargs="*")
    p.add_argument("--stop", action="store_true")
    e = sub.add_parser("e2e")
    e.add_argument("--out", required=True)
    e.add_argument("--config", required=True, help="e2e annotate config (becomes src/config.py for the run)")
    e.add_argument("--tmp", required=True, help="annotate_v2 tmp dir (PT_TMP)")
    e.add_argument("--gt", help="genotype_reads FASTA (PT_GT)")
    e.add_argument("--cache", help="copy annotate_v2 tmp files of an earlier run (dfam / sam caches)")
    e.add_argument("--no-reads", action="store_true", help="annotate events without the read sequences")
    a = ap.parse_args(argv)
    return run_pytest(a) if a.cmd == "pytest" else run_e2e(a)


if __name__ == "__main__":
    sys.exit(main())
