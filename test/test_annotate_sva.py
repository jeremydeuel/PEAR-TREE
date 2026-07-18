# PEAR-TREE - SVA composite-element annotation regression test
#
# Locks in the SVA recall fix in tools/annotate_v2.py (Insertion._sva_conclusion). SVA is a
# composite element (CCCTCT hexamer + Alu-like region + VNTR + SINE-R + poly-A); before the
# fix annotate mislabelled or dropped most real SVAs because (a) its poly-A tail can sit on
# either junction regardless of which strand the SVA HMM hits, while the Alu/L1 acceptance was
# strand-gated, (b) when the SVA hit was on the LEFT clip with a bare poly-A on the right, the
# whole Alu/L1 block (nested under `len(right_dfams)>0`) never ran, and (c) an SVA's Alu-like
# body scores an Alu model, so it was out-competed to ALU.
#
# These cases are REAL clips + Dfam bit scores from the 9x10 benchmark
# (analysis/mei9x10/mei9x10.combined.txt.gz + the cached nhmmscan table), so the poly-A and
# CCCTCT/AGAGGG hexamer detectors run on genuine sequence. Pure conclusion()-level logic, so
# no nhmmscan / bowtie2 / pysam runtime is needed.
#
# Run with:  python test/test_annotate_sva.py     (from the repo root)
#        or:  pytest test/test_annotate_sva.py

import os
import sys
import types
import importlib.util

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)

# annotate_v2 imports pysam at module load, but nothing on the conclusion() path uses it. Stub
# it so this logic-only test runs even in an environment without pysam installed.
try:
    import pysam  # noqa: F401
except Exception:
    sys.modules["pysam"] = types.ModuleType("pysam")

_spec = importlib.util.spec_from_file_location(
    "annotate_v2", os.path.join(REPO, "tools", "annotate_v2.py"))
annotate_v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(annotate_v2)
Insertion = annotate_v2.Insertion
Dfam_Annotation = annotate_v2.Dfam_Annotation


def _dfam(model, bits, strand="-"):
    """A Dfam_Annotation with a given model/bits/strand (other fields are filler; the SVA
    classifier only reads .model, .bits, .strand)."""
    line = [model, "ACC", "target", str(bits), "1e-10", "0.0",
            "1", "100", strand, "1", "100", "1", "100", "1400", "desc"]
    return Dfam_Annotation(line)


def _make(left_seq, right_seq, left_dfams=(), right_dfams=()):
    ins = Insertion("test", left_seq, right_seq)
    ins.left_dfams = list(left_dfams)
    ins.right_dfams = list(right_dfams)
    return ins


# --- real 9x10 SVA clips (kg_class == SVA, label_phylo == TP) ------------------------------

# chr17:78175213 - SVA_E on the LEFT clip (61.6 b), a bare poly-A on the RIGHT. Pre-fix this
# returned 'unknown': the left-dfam logic is nested under `len(right_dfams)>0`, which is empty
# here, so it was never consulted.
SVA_LEFT_POLYA_RIGHT = _make(
    "gaacaaaggtctctggttttcctaggcagaggaccctgcggccttggccttccgcagtgtttgtgtccctg"
    "GGTGGGTTTTTCTTTTAAACCCTCTGATTCACAAGTACTGAAGCCAGGAGAGTTTGCTCTGACCTTAGTCC"
    "TTCCCCCAACCCCCAGAGCACCATGCCCGGCTCTGTCCCTTTGCCAATTCAAAGCACCTGCCGCAACAGATGGCA",
    "GCTTCCCGGAAGCACTCCATGCACACCGCCGTCAGCTTCACCACGCTCTCGGCCAGCGGCACCAGGTTCAGG"
    "ATGGCCCCAAATGGCTGCGCTCCCATGGAAAGACAGAGGGCGCTCAGCCACCTTACCAGGTGGGTTTTTCTTTT"
    "tttttttttttttttttttttacggggggaagttttatttcaggaaat",
    left_dfams=[_dfam("SVA_E", 61.6, "-"), _dfam("LTR5_Hs", 32.7, "-")])

# chr10:69213563 - SVA_B on the LEFT (66.7 b), poly-A on the RIGHT. Same skipped-left-dfam bug.
SVA_B_LEFT_POLYA_RIGHT = _make(
    "ctcacctctgacgatgggcggccaggcagagacgctcctcacttcccagacggggtggcagccgggcagagg"
    "ctgcaatcttAGGCCACAGGTTTTTAGATCTCCTTAGTTTTTACCCAGTGGCCTTCTTCTCTTTTGGATCCC"
    "ATCCAAGACACCACATGACATTTAGTCGTCTCCTTAAGCTCCTCTTGGCTGTGACAGTTTCTCAGACTTGCTGTAGATGATCTC",
    "ATACAAAGAATTCCCATATACCCAACACATTTCCCCTATTAACATCTTACATTTCCCCCATTAACATCTTAC"
    "ATTAGTATGGTACGTTTGTTATAATTAATGAACCAAGACTGATACATTATTATTAACTAAAGGCCACAGGTTTTT"
    "ttttttttttttttttttttttattttgtttattttctttttttttttattttttttttttttttttttttttttttttt",
    left_dfams=[_dfam("SVA_B", 66.7, "-")])

# chr13:41828841 - the SVA-diagnostic (AGAGGG)n hexamer on the LEFT (SVA HMM only trails at
# 3.6 b, below the dominance floor), poly-A on the RIGHT. Exercises the hexamer-tandem path.
SVA_HEXAMER = _make(
    "ggagacggagacggagacggagagggagagggagagggagagggagagggagagggagagggagagggagc"
    "AACTTTTCTAATTTTTTTATAGCCGTTTTTATTTTGTTTCCTCCTCCAGGATTCAGTCAAGTCTCACCATT"
    "GCATTTGGTTGTTACATCTCCACAAGAAATGTTTAAATAAGCTACAATGGGTAATCCGATATCAGATGGCATAC",
    "CTTTATCCCTTAAGTATCTCAACATGCATCTCATAAGAATAGGAATATTCTCCTTCATAACCACAATATCAT"
    "TATCATACTTAAGGAAATTAGTAATTCAATAATATCATATAACATAGTCCATATTTAAACTTTTCTAATTTTTTT"
    "ttttttttttttttttttttttatcaagtttttttttttttttttttttttttttttttttttttttatatttttggttttttttgagttgttggt",
    left_dfams=[_dfam("SVA_F", 3.6, "-")])

# --- specificity guard: a real Alu must NOT be stolen to SVA -------------------------------
# chr8:127613044 - AluSc dominates the RIGHT clip (43.3 b) with only a trace SVA hit (8.4 b)
# and a single (non-tandem) CCCTCT; poly-A on the LEFT. The SVA classifier must decline
# (return None) so the Alu call stands. Protects Alu recall against SVA's Alu-like body.
ALU_NOT_SVA = _make(
    "ccaatcccccgcccaaacacccactgggacaaacacccaccaaataaaaaaaaaaaaaaaaattaaaaaaaa"
    "aaaaaaaaaaaaaaaaaaaaaaaaAAGATAGCTGTTTTTTTCCACTAAGCTTTGGGGTAGTTTGCATCACTG"
    "CTGTAGGCACCAAATCATGCCTTATCTAAAGGAGCAATTACTTCCTAGCTCGTGCCAACTGTTTCCCCATGGAAATAAGCACTCAGGTTAGCTAAATA",
    "CAGCCCCAAGCTGTTAGAGTACCTCCTAAGAGTCTTCCTGGCTCAGGCCCCAGACATTGCAGAGCAGAGACA"
    "AGCCATCCCCTCTGTGTCCTATTCAAGTTCCTGCCCCATGAAATCTGTGAGCATAATAAGATAGCTGTTTTT"
    "ggagccgagattgcaccactgcactccacccagacaacagagcaagactccgtctcaaaaaataaaataa",
    right_dfams=[_dfam("AluSc", 43.3, "+"), _dfam("SVA_B", 8.4, "-")])


def test_sva_left_dfam_polya_right():
    # bug (b): left SVA hit + bare poly-A right was 'unknown' (left-dfam block skipped).
    assert "SVA" in SVA_LEFT_POLYA_RIGHT.conclusion()


def test_sva_b_left_dfam_polya_right():
    assert "SVA" in SVA_B_LEFT_POLYA_RIGHT.conclusion()


def test_sva_hexamer_tandem():
    # the CCCTCT/AGAGGG hexamer is SVA-diagnostic; a weak SVA HMM hit + hexamer + poly-A -> SVA.
    c = SVA_HEXAMER.conclusion()
    assert "SVA" in c


def test_alu_not_relabelled_sva():
    # specificity: an Alu-dominant clip with only a trace SVA hit stays ALU.
    assert SVA_HEXAMER._sva_conclusion() is not None  # sanity: helper distinguishes cases
    assert ALU_NOT_SVA._sva_conclusion() is None
    assert "SVA" not in ALU_NOT_SVA.conclusion()


def test_element_class_maps_to_sva():
    VAC = annotate_v2.VariantAnnotationContainer
    assert VAC.element_class(SVA_LEFT_POLYA_RIGHT.conclusion()) == "SVA"


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL {name}: {e}")
    sys.exit(1 if fails else 0)
