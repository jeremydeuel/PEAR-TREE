# PEAR-TREE - proviral / lentiviral-vector LTR annotation regression test
#
# Locks in the provirus extension of tools/annotate_v2.py (Insertion._ltr_terminus,
# _ltr_u3_intact, _provirus_conclusion and the PROVIRUS_LTR class).
#
# An integrated retrovirus or lentiviral gene-therapy vector is NOT a retrotransposon, but its
# junctions look like an insertion's: a 5 bp target-site duplication (integrase, not TPRT) and a
# clip whose inserted part is the provirus. Unrecognised it classed as `unknown`, or as `ALU`
# when it landed in a reference Alu -- 204 unknown + 9 ALU out of PD49229's 214 vector
# integrations, every one confirmed against the published integration-site list of Spencer
# Chapman 2023 (the BCL11A sickle-cell gene-therapy cohort).
#
# The test is anchored: the LTR terminus must abut the insert/flank boundary. That is what makes
# it specific -- an unanchored 9-12-mer at 3-5 mismatches hits Alus all over the genome. On real
# data: 210 of PD49229's 214 with ZERO false positives, 0 of 232k loci in PD51635 + PD56961,
# and PD37580's single hit flagged as contamination rather than an integration.
#
# Sequences below are real clip_consensus strings from those runs (insert lower-case, aligned
# reference flank upper-case, as the combine evidence sidecar stores them).
#
# Run with:  python test/test_annotate_provirus.py     (from the repo root)
#        or:  pytest test/test_annotate_provirus.py

import os
import sys
import types
import importlib.util

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)

# pure sequence logic -- no pysam / nhmmscan / bowtie2 runtime needed
sys.modules.setdefault("pysam", types.ModuleType("pysam"))
_spec = importlib.util.spec_from_file_location("annotate_v2", os.path.join(REPO, "tools", "annotate_v2.py"))
av2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(av2)

Insertion = av2.Insertion
element_class = av2.VariantAnnotationContainer.element_class

# PD49229 1:8545322-8545328 -- SIN vector, provirus on the minus strand, inside a reference Alu
# (the upper-case flank is AluY: the locus classed as ALU before this check existed)
MINUS_L = ("ttttccacactgactaaaagggtctgagggatctctagttaccagagtcacacaacagacgggcacacactactttgagcactcaagg"
           "caagctttattgaggcttaagcagtgggttccctagttagccagagagctcccaggctcagatctggtctaaccagagagacccagta"
           "caggcaaaaagcagatcttgtcttctttgggagtgaattagcccttcc"
           "AGGCAGGAGAATCACTTGAACCTGGGAGGCGGAGGTTGTGGTGAGCCGAGATCGTGCCATCGCACTCCAGCCTGGGGAACAACAGCGA")
MINUS_R = ("CCAAGGTGGGCAGATCACCTGATGTAGGGAGTTCAAGACCAGCCTGACCAATATGGAGAAACCCCGTCTCTACTAAGAATACAAAATT"
           "CACTGGACGTGGTGCCGCATGTCTGTAATCCTAGCTACTTGGGAGGCTGAGGCAG"
           "tgctagagattttccacactgactaaaagggtctgagggatctctagttaccagagtcacacaacagacgggcacacactactttgag"
           "cactcaaggcaagctttattgaggcttaagcagtggg")

# PD49229 10:73602633-73602638 -- the other orientation (provirus on the plus strand). Both
# junctions resolve, and the right clip reaches R, so U3 is seen to be deleted.
PLUS_L = ("agggaacccactgcttaagcctcaataaagcttgccttgagtgctcaaagtagtgtgtgcccgtctgttgtgtgactctggtaactaga"
          "gatccctcagacccttttagtcagtgtggaaaatctctagca"
          "CTAATCATGACAAAGACTCGGCACTATGAAGCAATGAAGTTTCAGGCACACTAGCCTACACTTCCACAGCTGAAGTCTACCTGGCTAC")
PLUS_R = ("GGGCAGGTCCCAAACTATCATCAACTATTATAAAAGTCAACACTTCATTTGCGCCCAAAGTCACTCTTAATTTCATAACTCCTAAAAAC"
          "AATTTCTTAAAATTAATTTCAGAAAAATAAAAGATAATATTCTATCTAAT"
          "tggaagggctaattcactcccaaagaagacaagatctgctttttgcctgtactgggtctctctggttagaccagatctgagcctgggag"
          "ctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgctcaaagtagtgtgtgcccgtctgttgtgtga"
          "ctctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagtggcgcccgaa")

# PD37580 3:17443679-17443685 -- NOT an integration in that patient: present in all 84 colonies
# and already NOISE/ROOT to the joint step. Its left clip carries a U3-deleted LTR with vector
# sequence upstream of the terminus, its right clip an intact nef-overlapping U3 (and further
# along a CMV enhancer, GGGACTTTCC NF-kB sites -- an engineered construct). Both LTRs of one
# provirus share a structure, so the disagreement is the contamination signature.
CONFLICT_L = ("tcctttccttggctgctcgcctgtgttgccacctggattctgcgcgggacgtccttctgctacgtcccttcggccctcaatccagcgg"
              "accttccttcccgcggcctgctgccggctctgcggcctcttccgcgtcttcgccttcgccctcagacgagtcggatctccctttgggc"
              "cgcctccccgcctgagatcctttaagaccaatgacttacaaggcagctgtagatcttagccactttttaaaagaaaagggggga"
              "ctggaagggctaattcactcccaacgaagacaagatctgctttttgcttgtactgggtctctctggttagaccagatctgagcctggg"
              "agctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtg"
              "tgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagc"
              "AGTAGTAGTTCATGTCATCTTATTATTCAGTATTTATAACTTGCAAAGAAATGAATATCAGAGAGTGAGAGGCCC")
CONFLICT_R = ("ACCAAAATCATCCCAAACTTCCCACCCCATACCCTATTACCACTGCCAATTACCTGTGGTTTCATTTACTCTAAACCTGTGATTCCTC"
              "TGAATTATTTTCATTTTAAAGAAATTGTATTTGTTAAATATGTACTACAAACTTAGTAGT"
              "tggaagggctaattcactcccaaagaagacaagatatccttgatctgtggatctaccacacacaaggctacttccctgattagcagaa"
              "ctacacaccagggccaggggtcagatatccactgacctttggatggtgctacaagctagtaccagttgagccagataaggtagaagag"
              "gccaataaaggagagaacaccagcttgttacaccctgtgagcctgcatgggatggatgacccggagagagaagtgttagagtggaggt"
              "ttgacagccgcctagcatttcatcacgtggcccgagagctgcatccggagtacttcaagaactgctgatatcgagcttgctacaaggg"
              "actttccgctggggactttccagggaggcgtggcctgggcgggactggggagtggcgagccctcagatcctgcatataagcagctgct"
              "ttttgcctgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaa"
              "taaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctgg")

# a plain Alu insertion with a poly-A tail: the negative control the gate must not fire on
ALU_L = ("aaaaaaaaaaaaaaaaaaaaaaaaaagagacggagtctcgctctgtcgcccaggctggagtgcagtggcgcgatctcggctcactgcaa"
         "GCTTCTGCACAGCAAAAGAAACAATCAGCAGAATAGACAACCCAATGAGTGAGAGAAAATCTCACAATCTATACATCTGACAAAGGAC")
ALU_R = ("GAATGACGTGAACCCCAGGGGGCGGAGCCTGCAGTGAGCCGAGATTGCACCACTGCACTCCAGCCTGGG"
         "tgaaagagcgagattccgtctcaaaaaaaaaaaaaaaaaaaaa")

fails = []


def check(name, got, want):
    if got != want:
        fails.append(f"{name}:\n    got  {got!r}\n    want {want!r}")


def contains(name, got, needle):
    if needle not in (got or ""):
        fails.append(f"{name}: {needle!r} not in {got!r}")


def absent(name, got, needle):
    if needle in (got or ""):
        fails.append(f"{name}: {needle!r} unexpectedly in {got!r}")


# ---- 1. both orientations are detected, on both junctions
for label, (l, r), strand in (("minus-strand vector", (MINUS_L, MINUS_R), "-"),
                              ("plus-strand vector", (PLUS_L, PLUS_R), "+")):
    concl = Insertion(f"{label}", l, r)._provirus_conclusion()
    contains(f"{label}: detected", concl, "proviral LTR")
    contains(f"{label}: both junctions", concl, "both junctions")
    contains(f"{label}: orientation", concl, f"provirus {strand}")
    check(f"{label}: class", element_class(concl), "PROVIRUS_LTR")

# ---- 2. a SIN (U3-deleted) vector is distinguished from an intact provirus
contains("plus-strand vector: SIN", Insertion("p", PLUS_L, PLUS_R)._provirus_conclusion(),
         "lentiviral vector, U3-deleted (SIN)")

# ---- 3. the class beats the element call: 1:8545322 sits in a reference AluY, and the
#         upper-case flank alone annotates as ALU. The provirus must win.
check("vector-in-Alu does not class as ALU", element_class("AluY - (11-129|312)"), "ALU")
check("vector-in-Alu classes as PROVIRUS_LTR",
      element_class(Insertion("v", MINUS_L, MINUS_R)._provirus_conclusion()), "PROVIRUS_LTR")

# ---- 4. conflicting U3 structure between junctions is called out, not reported as a clean vector
cf = Insertion("contamination", CONFLICT_L, CONFLICT_R)._provirus_conclusion()
contains("conflict: flagged", cf, "likely contamination")
absent("conflict: not called a clean SIN vector", cf, "U3-deleted (SIN);")
check("conflict: still PROVIRUS_LTR", element_class(cf), "PROVIRUS_LTR")

# ---- 5. specificity: a plain Alu insertion, and junctions with nothing inserted
check("Alu insertion is not a provirus", Insertion("alu", ALU_L, ALU_R)._provirus_conclusion(), None)
check("empty junctions", Insertion("empty", "", "")._provirus_conclusion(), None)
check("insert shorter than the anchor",
      Insertion("short", "acgtacgtac" + "GGGGCCCC", "GGGGCCCC" + "acgtacgtac")._provirus_conclusion(), None)
# a terminus present but NOT at the junction must not fire (the anchoring requirement)
unanchored = ("tggaagggctaattcactcccaa" + "t" * 80).lower() + "ACGTACGTACGTACGTACGTACGT"
check("terminus far from the junction", Insertion("far", unanchored, "")._provirus_conclusion(), None)

# ---- 6. the additive [site:] / [SV:] suffixes are still stripped by element_class
check("suffix stripping", element_class("proviral LTR (U3 not covered; one junction: left U5+) [site: intron of FOO]"),
      "PROVIRUS_LTR")

# ---- 7. a single resolved junction is reported as such
one = Insertion("one-sided", PLUS_L, "")._provirus_conclusion()
contains("one junction", one, "one junction")

if fails:
    print(f"FAIL ({len(fails)})")
    for f in fails:
        print(" -", f)
    sys.exit(1)
print("test_annotate_provirus: all checks passed")
