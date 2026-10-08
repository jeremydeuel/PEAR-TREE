# PEAR-TREE - Dfam short-clip N-padding regression test
#
# Locks in Insertion._dfam_pad / get_dfam_fasta in tools/annotate_v2.py: a clip of 30..49 bp is
# N-padded to 50 bp for the Dfam scan only (the HMM null model under-scores short queries), with
# the Ns on the 3' end so Dfam ali_start/ali_end stay offsets into the unpadded clip. Clips under
# 30 bp, clips already >= 50 bp and (near-)pure poly-A/T clips go in unchanged; the bowtie2
# fasta (get_fasta) is never padded.
#
# Run with:  python test/test_annotate_dfam_pad.py     (from the repo root)
#        or:  pytest test/test_annotate_dfam_pad.py

import os
import sys
import types
import importlib.util

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)

# annotate_v2 imports pysam at module load, but nothing on this path uses it.
try:
    import pysam  # noqa: F401
except Exception:
    sys.modules["pysam"] = types.ModuleType("pysam")

_spec = importlib.util.spec_from_file_location(
    "annotate_v2", os.path.join(REPO, "tools", "annotate_v2.py"))
annotate_v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(annotate_v2)
Insertion = annotate_v2.Insertion

ALU = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAG"


def _pad(seq):
    return Insertion._dfam_pad(seq)


def test_short_clip_padded_on_3prime_end():
    clip = ALU[:35]
    assert _pad(clip) == clip + "N" * 15


def test_boundaries():
    assert _pad(ALU[:29]) == ALU[:29]              # < 30 bp: unchanged
    assert _pad(ALU[:30]) == ALU[:30] + "N" * 20   # 30 bp: padded
    assert _pad(ALU[:49]) == ALU[:49] + "N"        # 49 bp: padded
    assert _pad(ALU[:50]) == ALU[:50]              # >= 50 bp: unchanged
    assert _pad("") == ""


def test_polya_polyt_unchanged():
    assert _pad("A" * 40) == "A" * 40
    assert _pad("T" * 40) == "T" * 40
    near = "A" * 37 + "GCA"                        # 38/40 A = 95%: near-pure
    assert _pad(near) == near
    mixed = ALU[:20] + "A" * 20                    # element body + tail: padded
    assert _pad(mixed) == mixed + "N" * 10


def test_dfam_fasta_pads_but_bowtie2_fasta_does_not():
    # left junction = insert + REFERENCE flank; right = REFERENCE flank + insert
    left = ALU[:35].lower() + "ACGTACGTAC"
    right = "TTGACCATGA" + ("a" * 32)              # poly-A right clip: unchanged
    ins = Insertion("chr1:100-115", left, right)
    dfam = ins.get_dfam_fasta()
    assert dfam == (f">chr1:100-115:R\n{'A' * 32}\n"
                    f">chr1:100-115:L\n{ALU[:35]}{'N' * 15}\n")
    assert "N" not in ins.get_fasta()


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
