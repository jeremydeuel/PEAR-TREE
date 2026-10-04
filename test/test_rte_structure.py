# PEAR-TREE - tools/rte element / 5'-structure / tag tests on synthetic insertions.
#
# Every case is simulated with test/rte_sim.py from the stand-in library in
# test/fixtures/rte_library (L1HS consensus from L1Base, Dfam Alu/SVA, real hg38/hs1 flanks):
# reads tile the allele, junction strings + evidence rows follow the SPEC formats.
#
# Run:  pytest test/test_rte_structure.py
import os
import random
import sys
from functools import lru_cache

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rte_sim  # noqa: E402
from tools.rte.annotator import RteAnnotator  # noqa: E402
from tools.rte.genome import FastaGenome  # noqa: E402
from tools.rte.sequtil import rc  # noqa: E402

LIB = rte_sim.library()
L1 = LIB.consensus["L1HS"]
CEND = LIB.cons_end["L1HS"]
A = "A" * 25


def annotate(insert, ann_kw=None, inp_kw=None, **kw):
    inp, ev, genome, info = rte_sim.build(insert, **kw)
    for k, v in (inp_kw or {}).items():
        setattr(inp, k, v)
    ann = RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome, **(ann_kw or {}))
    return ann.annotate(inp, ev)


@lru_cache(None)
def full_length(strand=1):
    return annotate(L1[:CEND] + A, strand=strand)


@lru_cache(None)
def truncated():
    return annotate(L1[5500:CEND] + A)


@lru_cache(None)
def twin_priming():
    fs = CEND - 2000          # forward 2 kb
    p2 = fs - 14              # 14 bp deleted at the inversion point
    return annotate(rc(L1[p2 - 1000:p2]) + L1[fs:CEND] + A)


# ------------------------------------------------------------------------------------- L1
@pytest.mark.parametrize("strand", [1, -1])
def test_full_length_l1(strand):
    r = full_length(strand)
    assert r.element == "L1"
    assert r.structure == "FULL_LENGTH"
    assert r.covered_5p <= 10 and r.covered_3p >= CEND - 10
    assert r.strand == strand
    assert r.nearest_active != "." and r.element_identity >= 0.98


def test_truncated_at_5_5kb():
    r = truncated()
    assert r.element == "L1" and r.structure == "TRUNCATED_5P"
    assert abs(r.detail["j5"] - 5500) <= 5
    assert r.detail.get("j5_feature") == "ORF2"


def test_twin_priming_inversion_with_14bp_deletion():
    r = twin_priming()
    assert r.element == "L1"
    assert r.structure == "INVERTED_5P"
    assert r.detail["inv_junction"] == "del14"
    assert r.detail["inv_exact"] == 1
    assert abs(r.detail["fwd_start"] - (CEND - 2000)) <= 3
    # two covered blocks: the inverted 1 kb and the forward 2 kb
    assert len(r.covered_intervals) == 2
    assert "twin_priming_ge590" in r.tprt_points
    assert "FOLDBACK_INVDUP_5P" not in r.tags


def test_twin_priming_with_template_switch():
    ins = rc(L1[3000:3500]) + L1[3600:4000] + rc(L1[4100:4500]) + L1[4600:CEND] + A
    r = annotate(ins)
    assert r.structure == "INVERTED_5P_SWITCH"


def test_foldback_inverted_duplication():
    # inverted piece is the mirror image of the start of the forward piece (snap-back)
    fs = CEND - 1500
    r = annotate(rc(L1[fs:fs + 400]) + L1[fs:CEND] + A)
    assert r.structure.startswith("INVERTED_5P")
    assert "FOLDBACK_INVDUP_5P" in r.tags


def test_5p_unresolved_without_5p_junction_reads():
    ins = L1[4000:CEND] + A
    # drop every read and the junction consensus crossing the 5' junction
    inp, ev, genome, info = rte_sim.build(ins, exclude=(0,))
    inp.right_seq = ""
    ev.junctions.pop("RIGHT")
    r = RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome).annotate(inp, ev)
    assert r.element == "L1"
    assert r.structure == "5P_UNRESOLVED"


# ------------------------------------------------------------------------------------- Alu / SVA
def test_alu_minus_strand():
    a = LIB.consensus["ALU_YA5"][:LIB.cons_end["ALU_YA5"]]
    r = annotate(a + A, strand=-1, tsd=12)
    assert r.element == "ALU" and r.structure == "FULL_LENGTH"
    assert r.consensus.startswith("ALU_Y")
    assert r.tsd_len == 12
    assert not any(t.startswith("TD3P") for t in r.tags)


def test_sva_with_5p_transduction():
    f5 = list(rte_sim.fasta("flanks_5p_sva.fa").values())[0]
    s = LIB.consensus["SVA_F"][:LIB.cons_end["SVA_F"]]
    r = annotate(f5[-300:] + s + A)
    assert r.element == "SVA"
    assert "TD5P" in r.tags
    assert r.structure == "FULL_LENGTH"


# ------------------------------------------------------------------------------------- other types
def test_solitary_polya():
    r = annotate("A" * 30, tsd=12)
    assert r.element == "POLYA_ONLY"
    assert r.tprt_call in ("TPRT", "LIKELY_TPRT")


def test_orphan_transduction_element():
    fl = rte_sim.fasta("flanks_3p.fa")["TD_UID-5"]
    r = annotate(fl[500:900] + A)
    assert r.element == "ORPHAN_TD"
    assert "TD3P_SOURCE=TD_UID-5" in r.tags


def test_en_independent():
    # internal L1 piece, no tail, no TSD, no EN motif: both ends truncated
    r = annotate(L1[3000:5000], tsd=0, motif=None)
    assert r.element == "L1" and "EN_INDEPENDENT" in r.tags
    assert r.tprt_call not in ("TPRT", "LIKELY_TPRT")


def test_l1_mediated_deletion():
    r = annotate(L1[5000:CEND] + A, tsd=-300)
    assert r.tsd_len == -300
    assert "TSD_DELETION" in r.tags and "L1_MED_DELETION" in r.tags


def test_chimeric_ends_artefact():
    alu = LIB.consensus["ALU_YA5"]
    r = annotate(alu[:200] + L1[5700:CEND] + A)
    assert "CHIMERIC_ENDS" in r.tags
    assert r.tprt_call not in ("TPRT", "LIKELY_TPRT")


# ------------------------------------------------------------------------------------- pseudogene
def _pseudogene_setup():
    rng = random.Random(7)
    g = "".join(rng.choice("ACGT") for _ in range(5000))
    exons = [(500, 700), (1500, 1650), (2800, 3100)]
    remap = FastaGenome(records={"chrG": g})
    mrna = "".join(g[s:e] for s, e in exons)
    return remap, {"PGENE": [("chrG", s, e) for s, e in exons]}, mrna


def test_pseudogene_with_exon_junction_read():
    remap, exons, mrna = _pseudogene_setup()
    r = annotate(mrna + A, ann_kw={"remap_genome": remap, "exons_by_gene": exons},
                 inp_kw={"pseudogene_genes": ["PGENE"]})
    assert r.element == "PSEUDOGENE"
    assert "EXON_JUNCTION" in r.tags
    assert "exon_junction" in r.tprt_points


def test_pseudogene_without_exon_junction_read_is_candidate_only():
    remap, exons, mrna = _pseudogene_setup()
    # no read (and no junction consensus) crosses the exon1|exon2 and exon2|exon3 joins
    r = annotate(mrna + A, exclude=(200, 350), clip_len=100,
                 ann_kw={"remap_genome": remap, "exons_by_gene": exons},
                 inp_kw={"pseudogene_genes": ["PGENE"]})
    assert r.element != "PSEUDOGENE"
    assert "PSEUDOGENE_CANDIDATE" in r.tags
    assert "EXON_JUNCTION" not in r.tags


# ------------------------------------------------------------------------------------- templated / pre-mRNA
def test_templated_local_segment():
    _, _, _, info = rte_sim.build(L1[5000:CEND] + A)          # same seed -> same reference
    L = info["L"]
    templ = info["ref"][L + 120:L + 160]                      # 40 bp from ~120 bp downstream
    r = annotate(templ + L1[5000:CEND] + A)
    assert "TEMPLATED_LOCAL" in r.tags
    assert r.detail["templated_dist"] <= 250


class _GeneModel:
    """Duck-typed annotate_v2.GeneModel: one + strand gene with two exons."""
    def _candidates(self, contig, p):
        return [(10000, 30000, "HOSTG", "+", ((10000, 10500), (29000, 30000)))]

    def _genic_feature(self, exons, strand, p):
        return (5, "intron", "intron")


def test_premrna_coinsert():
    kw = dict(glen=40000, site=4000)
    _, _, _, info = rte_sim.build(L1[5000:CEND] + A, **kw)
    intron = info["ref"][20000:20120]                        # 16 kb away, inside HOSTG's intron
    r = annotate(L1[5000:CEND] + "A" * 6 + intron + A, ann_kw={"gene_model": _GeneModel()}, **kw)
    assert "PREMRNA_COINSERT" in r.tags
    assert "TD3P" not in r.tags                       # host-gene sequence, not a transduction
    assert r.detail["premrna"].startswith("HOSTG:intron")
