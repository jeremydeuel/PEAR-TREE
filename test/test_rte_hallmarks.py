# PEAR-TREE - tools/rte hallmark tests: TSD (genome and flank-only), target-site deletion, EN
# motif 7-mer + PCAWG mismatch bins, strand-consistent poly-A, beyond-poly-A support, poly-A
# slippage context, both-sided poly-A, fold-back clip.
#
# Run:  pytest test/test_rte_hallmarks.py
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rte_sim  # noqa: E402
from tools.rte import hallmarks as H  # noqa: E402
from tools.rte.annotator import RteAnnotator  # noqa: E402
from tools.rte.genome import FastaGenome  # noqa: E402
from tools.rte.sequtil import rc  # noqa: E402

LIB = rte_sim.library()
L1 = LIB.consensus["L1HS"]
CEND = LIB.cons_end["L1HS"]
INS = L1[5000:CEND] + "A" * 25


def run(insert=INS, genome_on=True, **kw):
    inp, ev, genome, info = rte_sim.build(insert, **kw)
    ann = RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome if genome_on else None)
    return ann.annotate(inp, ev), info


@pytest.mark.parametrize("strand", [1, -1])
def test_tsd_sequence_from_genome(strand):
    r, info = run(strand=strand, tsd=15)
    assert r.tsd_len == 15
    assert r.tsd_seq == info["ref"][info["L"]:info["R"]]
    assert "tsd_4_25" in r.tprt_points


def test_tsd_flank_only_fallback():
    r, info = run(tsd=15, genome_on=False)
    assert r.tsd_len == 15 and r.tsd_seq == info["ref"][info["L"]:info["R"]]
    assert r.en_motif == ""      # no genome -> no motif


def test_tsd_with_one_mismatch_in_reads():
    # a sample SNV in one copy of the duplication still verifies (<= 1 mismatch)
    inp, ev, genome, info = rte_sim.build(INS, tsd=15)
    rs = list(ev.junctions["RIGHT"].clip_consensus)
    i = next(k for k, c in enumerate(rs) if c.islower()) - 3   # inside the TSD copy on the right
    rs[i] = "A" if rs[i] != "A" else "C"
    _, lf, rf, _ = H.split_junction(ev.junctions["LEFT"].clip_consensus, "".join(rs))
    si = H.locate_site(inp.title, lf, rf, genome)
    H.target_site(si, lf, rf, genome)
    assert si.tsd_len == 15 and si.tsd_verified
    # two mismatches no longer verify the duplicated copy
    rf2 = rf[:-8] + ("A" if rf[-8] != "A" else "C") + rf[-7:]
    si2 = H.target_site(H.locate_site(inp.title, lf, rf, genome), lf, rf2, genome)
    assert not si2.tsd_verified


def test_target_site_deletion():
    r, info = run(tsd=-8)
    assert r.tsd_len == -8
    assert "TSD_DELETION" in r.tags
    assert "tsd_deletion_le20" in r.tprt_points


def test_tsd_over_50_is_artefact_indicator():
    r, info = run(tsd=60)
    assert r.tsd_len == 60
    assert "tsd_gt50" in r.tprt_points
    assert r.tprt_call not in ("TPRT", "LIKELY_TPRT")


@pytest.mark.parametrize("motif,mm", [("TTTTTAA", 0), ("TTTTTGA", 0), ("CTTTTAA", 0),
                                      ("TCTTTAA", 1), ("TTCTCAA", 2), ("GGCGCCC", 5)])
@pytest.mark.parametrize("strand", [1, -1])
def test_en_motif_strand_corrected(motif, mm, strand):
    r, info = run(motif=motif, strand=strand)
    assert r.en_motif == f"{motif[:5]}/{motif[5:]}"
    assert r.en_mismatches == mm


def test_en_bins():
    assert [H.en_bin(m) for m in (0, 1, 2, 3, 4, 5, None)] == ["0-1", "0-1", "2", "3", "4-5", "4-5", "."]


def test_polya_length_on_strand_consistent_side():
    # the L1 3' end is A-rich (...TATAAT), so the tolerant edge run may add a few bases
    r, _ = run(L1[5000:CEND] + "A" * 31, strand=-1, polya_len=0)
    assert 31 <= r.polya_len <= 35 and r.strand == -1
    r, _ = run(L1[5000:CEND] + "A" * 31, strand=-1, polya_len=40)   # pooled evidence median wins
    assert r.polya_len == 40


def test_beyond_polya_needs_two_fragments_for_extra_points():
    r, _ = run()
    assert r.beyond_polya and r.beyond_polya_support >= 2
    assert "beyond_polya:+2" in r.tprt_points and "beyond_polya_2frag" in r.tprt_points
    # keep only the junction consensus (no reads): sequence beyond the tail, but < 2 fragments
    inp, ev, genome, info = rte_sim.build(INS)
    ev.reads = []
    r = RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome).annotate(inp, ev)
    assert r.beyond_polya and "beyond_polya:+2" in r.tprt_points
    assert "beyond_polya_2frag" not in r.tprt_points


def test_polya_slippage_artefact():
    r, _ = run("A" * 14, tsd=0, ref_homopolymer=15)
    assert r.detail.get("slippage", "").startswith("ref A")
    assert "slippage_context" in r.tprt_points
    assert "en_0_1" not in r.tprt_points          # motif void inside a homopolymer
    assert r.tprt_call == "ARTEFACT_LIKE"


def test_both_sided_polya_flag():
    pa = H.polya_info("GATTACA" + "A" * 15, "T" * 14 + "GATTACA")
    assert pa.both_sided
    pa = H.polya_info("GATTACA" + "A" * 15, "A" * 15)        # solitary poly-A: consistent
    assert not pa.both_sided and pa.strand == 1


def test_edge_run_tolerates_one_substitution():
    assert H.edge_run("CCCG" + "AAAAAAGAAAAAAA", at_end=True) == ("A", 14)
    assert H.edge_run("TTTTCTTTTTTTGGC", at_end=False) == ("T", 12)


def test_foldback_clip():
    flank = "ACGTTGCAAGGCTTACCGATTGCAGGT"
    right_insert = rc(flank)[:20] + "GGGCCC"
    assert H.foldback("", "", right_insert, flank) >= 15
    left_flank = "TTGACCAGGATCCGATGCCATAGGCATT"
    left_insert = "CCCAAA" + rc(left_flank[:22])
    assert H.foldback(left_insert, left_flank, "", "") >= 15
    assert H.foldback("ACGTACGTAGCTAGCTAGCATCGA", left_flank, "TTTTTTTTT", flank) == 0
