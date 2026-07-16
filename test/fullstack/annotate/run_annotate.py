#!/usr/bin/env python3
"""Run annotate_v2 on a local (workstation / test-harness) host.

annotate_v2 reads its paths from src.config.CONFIG['annotate'], which ships with the cluster
layout (dfamscan.pl + full Dfam HMM library). This driver patches that dict for a local run:
the family HMM library built by build_hmm.sh, the nhmmscan-direct DFAM back-end (no dfamscan.pl),
the hs1 RepeatMasker track, the hs1 pseudogene exon track, and the local bowtie2 + hs1 index.
It then runs the annotation on an explicit combined / genotypes pair (no fixed genotyping_test/
layout required), writing the tmp dfam/sam/fasta next to the output.

Usage:
  run_annotate.py --combined step2.combined.txt.gz --genotypes step2.genotypes.csv.gz \
      --out annotate.txt [--workdir DIR]
Env (defaults in brackets):
  PT_HMM   [<scriptdir>/peartree_rte.hmm]        HMM library (build_hmm.sh)
  PT_EXON  [<scriptdir>/pseudogene_exons.hs1.bed] pseudogene exon track (build_exon_track.py)
  PT_RMSK  [<repo>/hs1.repeatMasker.out.gz]        hs1 RepeatMasker
  BOWTIE2  [$(command -v bowtie2)]                 bowtie2 executable
  HS1_BT2  [~/Downloads/hs1]                        bowtie2 hs1 index prefix (clip remap target)
"""
import argparse, importlib.util, os, sys

SCRIPT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(SCRIPT, "..", "..", ".."))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", required=True)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", default=None, help="dir for tmp dfam/sam/fa (default: alongside --out)")
    a = ap.parse_args()
    workdir = a.workdir or (os.path.dirname(os.path.abspath(a.out)) or ".")
    os.makedirs(workdir, exist_ok=True)

    sys.path.insert(0, REPO)
    import src.config as cfg
    A = cfg.CONFIG["annotate"]; CI = cfg.CONFIG["combine_insertions"]
    A["hmm"] = os.environ.get("PT_HMM", os.path.join(SCRIPT, "peartree_rte.hmm"))
    A["exon_annotation"] = os.environ.get("PT_EXON", os.path.join(SCRIPT, "pseudogene_exons.hs1.bed"))
    A["rmsk"] = os.environ.get("PT_RMSK", os.path.join(REPO, "hs1.repeatMasker.out.gz"))
    A["dfamscan"] = None          # force the nhmmscan --dfamtblout back-end
    A["hmmer"] = None             # nhmmscan already on PATH
    # explicit combined/genotypes; tmp files (dfam/sam/fa.gz) land in workdir keyed by sample
    A["insertions_file"] = lambda s: a.combined
    A["genotyping_file"] = lambda s: a.genotypes
    A["tmp"] = lambda ext: (lambda s: os.path.join(workdir, f"{s}.{ext}"))
    CI["bowtie2_executable"] = os.environ.get("BOWTIE2") or _which("bowtie2")
    CI["bowtie2_index2"] = os.environ.get("HS1_BT2", os.path.expanduser("~/Downloads/hs1"))

    spec = importlib.util.spec_from_file_location("annotate_v2", os.path.join(REPO, "tools", "annotate_v2.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    sample = os.path.basename(a.combined).split(".")[0]
    f = m.VariantAnnotationContainer(sample, a.out)   # setup logging -> console
    with open(a.out, "w") as fh:                       # per-insertion report -> --out
        old = sys.stdout
        sys.stdout = fh
        try:
            f.print()
        finally:
            sys.stdout = old
    print(f"wrote annotation report to {a.out}")


def _which(x):
    from shutil import which
    return which(x) or x


if __name__ == "__main__":
    main()
