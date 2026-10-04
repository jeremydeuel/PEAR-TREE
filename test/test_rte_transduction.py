# PEAR-TREE - tools/rte 3' transduction tests: known source lookup in flanks_3p.fa (partnered and
# orphan), the transduction endpoint offset, and the SPEC novel-source rule on a real hs1 region
# (test/fixtures/rte_library/novel_source_region.fa: a reference L1HS, chrX:11289795-11295826 (+),
# plus 16 kb downstream; novel_source.rmsk.out: its RepeatMasker rows).
#
# Run:  pytest test/test_rte_transduction.py
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rte_sim  # noqa: E402
from tools.rte.annotator import RteAnnotator  # noqa: E402
from tools.rte.assembly import Segment  # noqa: E402
from tools.rte.genome import FastaGenome  # noqa: E402
from tools.rte.transduction import (known_source, L1Rmsk, MappyLocator,  # noqa: E402
                                    NovelSourceFinder)

LIB = rte_sim.library()
L1 = LIB.consensus["L1HS"]
CEND = LIB.cons_end["L1HS"]
A = "A" * 25
FLANK = rte_sim.fasta("flanks_3p.fa")["TD_UID-5"]
REGION_FA = os.path.join(rte_sim.FIX, "novel_source_region.fa")
RMSK = os.path.join(rte_sim.FIX, "novel_source.rmsk.out")
REGION = FastaGenome(REGION_FA)
SRC_END = 11295826            # reference L1HS (+) ends here (0-based exclusive) on hs1 chrX


def run(insert, ann_kw=None, **kw):
    inp, ev, genome, info = rte_sim.build(insert, **kw)
    return RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome, **(ann_kw or {})).annotate(inp, ev)


def test_partnered_td3p_from_real_flank():
    # 5'-truncated L1 + its own short tail + 400 bp of the TD_UID-5 downstream flank + poly-A
    r = run(L1[5600:CEND] + "A" * 8 + FLANK[20:420] + A)
    assert r.element == "L1"
    assert "TD3P" in r.tags and "TD3P_SOURCE=TD_UID-5" in r.tags
    # transduction endpoint = offset of the tag's 3' end inside the source flank
    assert abs(r.detail["td_end"] - 420) <= 5
    assert "NOVEL_SOURCE" not in r.tags


def test_partnered_td3p_minus_strand():
    r = run(L1[5600:CEND] + "A" * 8 + FLANK[20:420] + A, strand=-1)
    assert "TD3P_SOURCE=TD_UID-5" in r.tags and r.strand == -1


def test_td3p_source_matches_5p_element():
    # the 5' part comes from exactly the source element UID-5 -> concordance bonus
    uid5 = rte_sim.fasta("l1_intact.fa")["UID-5"]
    r = run(uid5[:len(uid5) - 40] + "A" * 8 + FLANK[20:420] + A)
    assert "TD3P_SOURCE=TD_UID-5" in r.tags
    assert r.nearest_active == "UID-5" and r.element_identity >= 0.999
    assert "td_source_matches_5p" in r.tprt_points


def test_orphan_transduction():
    r = run(FLANK[500:900] + A)
    assert r.element == "ORPHAN_TD"
    assert "TD3P_SOURCE=TD_UID-5" in r.tags
    assert abs(r.detail["td_end"] - 900) <= 5


def test_known_source_unit():
    segs = [Segment(0, 100, "FLANK3P", "TD_UID-8", 200, 300, 1, 1.0, 100),
            Segment(0, 50, "FLANK3P", "TD_UID-5", 10, 60, 1, 1.0, 50),
            Segment(0, 60, "FLANK3P", "TD_UID-8", 300, 360, 1, 1.0, 60)]
    sc = known_source(segs, LIB)
    assert sc.source_id == "TD_UID-8" and sc.td_end == 360 and sc.td_start == 200
    assert known_source([Segment(0, 100, "FLANK3P", "TD_UID-8", 0, 100, -1, 1.0, 100)], LIB) is None


def test_l1_rmsk_upstream_lookup():
    rm = L1Rmsk(RMSK)
    hits = rm.upstream_of("chrX", SRC_END + 500, SRC_END + 900, "+", 15000)
    assert hits and hits[0][3] == "L1HS"
    assert rm.upstream_of("chrX", SRC_END + 500, SRC_END + 900, "-", 15000) == []
    assert rm.upstream_of("chrX", SRC_END + 20000, SRC_END + 20400, "+", 15000) == []


def _finder():
    return NovelSourceFinder(LIB, rmsk=RMSK, locator=MappyLocator(REGION_FA), remap_genome=REGION)


def test_novel_source_rule_accepts_downstream_of_young_l1():
    tag = REGION.fetch("chrX", SRC_END + 374, SRC_END + 774)     # unique, between L1MB8 and LTR50
    sc = _finder().find(tag)
    assert sc is not None and sc.novel
    assert sc.source_id == "novel:chrX:11289795-11295826"
    assert sc.identity >= 0.95
    assert sc.td_start == 374 and sc.td_end == 774


def test_novel_source_rule_rejects_wrong_strand_and_far():
    from tools.rte.sequtil import rc
    tag = REGION.fetch("chrX", SRC_END + 374, SRC_END + 774)
    assert _finder().find(rc(tag)) is None            # tag reads antisense to the L1
    far = REGION.fetch("chrX", 11311000 - 20, 11311300)   # > 15 kb downstream
    assert _finder().find(far) is None


def test_novel_source_end_to_end():
    tag = REGION.fetch("chrX", SRC_END + 374, SRC_END + 774)
    r = run(L1[5600:CEND] + "A" * 8 + tag + A,
            ann_kw={"locator": MappyLocator(REGION_FA), "rmsk": RMSK, "remap_genome": REGION})
    assert "TD3P" in r.tags and "NOVEL_SOURCE" in r.tags
    assert "TD3P_SOURCE=novel:chrX:11289795-11295826" in r.tags
    assert r.detail["source_identity"] >= 0.95
    assert "novel_source" in r.tprt_points


def test_unsourced_td3p_without_locator():
    tag = REGION.fetch("chrX", SRC_END + 374, SRC_END + 774)
    r = run(L1[5600:CEND] + "A" * 8 + tag + A)
    assert "TD3P" in r.tags
    assert not any(t.startswith("TD3P_SOURCE=") for t in r.tags)
