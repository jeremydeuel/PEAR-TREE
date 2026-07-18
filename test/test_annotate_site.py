# PEAR-TREE - insertion-SITE annotation regression test
#
# Locks in the insertion-site extension of tools/annotate_v2.py (GeneModel + Insertion.site()).
# The element logic answers "what was inserted"; this answers "WHERE did it land": the annotator
# takes the junction's own reference locus (the title "contig:start-end", on the sample's
# discovery/BAM genome) and classifies it against a gene model into
#     exon / splice donor / splice acceptor / polypyrimidine tract / branch point (lariat) /
#     intron / promoter (core = disrupts, proximal = near) / intergenic
# strand-aware (donor vs acceptor and the TSS both depend on gene strand). The note is ADDITIVE:
# conclusion() appends "[site: ...]" without changing the element class, and element_class()
# strips it. Pure classifier logic, so no nhmmscan / bowtie2 / pysam runtime is needed.
#
# Run with:  python test/test_annotate_site.py     (from the repo root)
#        or:  pytest test/test_annotate_site.py

import os
import sys
import types
import importlib.util

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)

try:
    import pysam  # noqa: F401
except Exception:
    sys.modules["pysam"] = types.ModuleType("pysam")

_spec = importlib.util.spec_from_file_location(
    "annotate_v2", os.path.join(REPO, "tools", "annotate_v2.py"))
annotate_v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(annotate_v2)
Insertion = annotate_v2.Insertion
GeneModel = annotate_v2.GeneModel
VAC = annotate_v2.VariantAnnotationContainer

# A two-gene model on chr1, both with two exons and a wide (3.7 kb) intron so the splice/intron
# probes sit well clear of the promoter windows:
#   PLUS  (+): exon [1000,1300)  intron [1300,5000)  exon [5000,5300)   TSS = 1000
#   MINUS (-): exon [8000,8300)  intron [8300,12000) exon [12000,12300) TSS = 12300
# For '+' the donor is the intron's LEFT boundary (1300) and the acceptor its RIGHT (5000);
# for '-' it is mirrored (donor at 12000, acceptor at 8300).
_TRACK = "\n".join([
    "chr1\t1000\t1300\tPLUS\t+",
    "chr1\t5000\t5300\tPLUS\t+",
    "chr1\t8000\t8300\tMINUS\t-",
    "chr1\t12000\t12300\tMINUS\t-",
]) + "\n"

_track_path = os.path.join(
    os.environ.get("TMPDIR", "/tmp"), "peartree_test_gene_model.tsv")
with open(_track_path, "w") as _fh:
    _fh.write(_TRACK)
GM = GeneModel(_track_path)   # default windows (== the config defaults)


def _region(contig, lo, hi):
    return GM.annotate(contig, lo, hi)


def test_exon():
    r = _region("chr1", 1150, 1150)
    assert r[0] == "exon" and r[1] == "PLUS" and r[2] == "+", r
    assert "exonic of PLUS (+)" in r[3], r


def test_splice_donor_plus():
    r = _region("chr1", 1305, 1305)          # 6 nt into the intron from the + donor boundary
    assert r[0] == "splice_donor", r
    assert "splice donor of PLUS (+)" in r[3], r


def test_splice_acceptor_plus():
    r = _region("chr1", 4998, 4998)          # 2 nt before the + acceptor boundary
    assert r[0] == "splice_acceptor", r
    assert "splice acceptor of PLUS (+)" in r[3], r


def test_ppt_plus():
    r = _region("chr1", 4990, 4990)          # 10 nt upstream of the + acceptor
    assert r[0] == "polypyrimidine_tract", r


def test_branch_plus():
    r = _region("chr1", 4970, 4970)          # 30 nt upstream of the + acceptor
    assert r[0] == "branch_point", r
    assert "lariat" in r[3], r


def test_deep_intron_plus():
    r = _region("chr1", 3000, 3000)
    assert r[0] == "intron" and r[1] == "PLUS", r


def test_splice_donor_minus():
    r = _region("chr1", 11995, 11995)        # 5 nt into the intron from the - donor boundary (12000)
    assert r[0] == "splice_donor" and r[2] == "-", r


def test_splice_acceptor_minus():
    r = _region("chr1", 8301, 8301)          # 2 nt from the - acceptor boundary (8300)
    assert r[0] == "splice_acceptor" and r[2] == "-", r


def test_ppt_minus():
    r = _region("chr1", 8310, 8310)
    assert r[0] == "polypyrimidine_tract" and r[2] == "-", r


def test_branch_minus():
    r = _region("chr1", 8330, 8330)
    assert r[0] == "branch_point" and r[2] == "-", r


def test_promoter_core_plus():
    r = _region("chr1", 800, 800)            # 200 bp upstream of the + TSS (1000), not in the body
    assert r[0] == "promoter_core", r
    assert "disrupts promoter of PLUS (+)" in r[3], r


def test_promoter_proximal_plus():
    r = _region("chr1", 500, 500)            # 500 bp upstream of the + TSS
    assert r[0] == "promoter_proximal", r
    assert "near promoter of PLUS (+)" in r[3], r


def test_promoter_core_minus():
    r = _region("chr1", 12500, 12500)        # 200 bp upstream (sense) of the - TSS (12300)
    assert r[0] == "promoter_core", r
    assert "disrupts promoter of MINUS (-)" in r[3], r


def test_intergenic():
    r = _region("chr1", 50000, 50000)
    assert r[0] == "intergenic" and r[3] == "intergenic", r


def test_unknown_contig_is_intergenic():
    r = _region("chrZ", 1000, 1000)
    assert r[0] == "intergenic", r


def test_chr_prefix_tolerance():
    # a numeric-contig title (GRCh37 style) must still hit a chr-prefixed model
    r = _region("1", 1150, 1150)
    assert r[0] == "exon" and r[1] == "PLUS", r


def test_additive_and_class_strip():
    # site note is appended to the element conclusion without changing the class, and
    # element_class() strips the [site: ...] suffix. Use an artefact junction (no dfam/map).
    Insertion.gene_model = GM
    try:
        ins = Insertion("chr1:1150-1150", "acgtACGT", "ACGTacgt")   # -> 'artefact'
        c = ins.conclusion()
        assert c.startswith("artefact"), c
        assert "[site: exonic of PLUS (+)" in c, c
        assert VAC.element_class(c) == "artefact", c
    finally:
        Insertion.gene_model = None


def test_no_model_no_site_suffix():
    # with no gene model configured (the default) conclusion() appends nothing, so exact-string
    # assertions in the other suites are unaffected.
    Insertion.gene_model = None
    ins = Insertion("chr1:1150-1150", "acgtACGT", "ACGTacgt")
    assert ins.conclusion() == "artefact", ins.conclusion()
    assert ins.site() is None


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok   {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL {name}: {e}")
    sys.exit(1 if fails else 0)
