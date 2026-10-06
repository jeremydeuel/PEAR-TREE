# PEAR-TREE - tools/rte: which 3'/5' flank hits may name a transduction source.
#
# PD37590 had 13 spurious TD3P_SOURCE calls (ALU/SVA/L1 inserts "sourced" 1-14 kb into some L1
# flank): every one was a 20-55 bp, 90-95 % identity hit inside a SOFT-MASKED (RepeatMasker)
# stretch of the 15 kb flank -- the insert's own Alu/L1/SVA sequence matching a repeat copy
# in that flank. The true orphan transduction chr7:131206236 -> L1_chrX_11289796_f hits
# unmasked flank at 100 %. Rules locked in here (structure.py):
#   * a FLANK3P/FLANK5P hit on >= td_max_masked_frac (0.5) soft-masked flank is no source;
#   * the source class must match the insert class (no L1 source behind an Alu/SVA);
#   * a short (< td_min_flank_bp) tail-position hit on an A/T-rich flank piece (the source's own
#     poly-A remnant at the flank start; chr10:115802225 td_end=32) is no source.
#
# Run:  pytest test/test_rte_td_sources.py
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rte_sim  # noqa: E402
from tools.rte.annotator import RteAnnotator  # noqa: E402

LIB = rte_sim.library()
L1 = LIB.consensus["L1HS"]
CEND = LIB.cons_end["L1HS"]
ALU = LIB.consensus["ALU_Y"][:LIB.cons_end["ALU_Y"]]
A = "A" * 25
FLANK = rte_sim.fasta("flanks_3p.fa")["TD_UID-5"]


@pytest.fixture(scope="module")
def masked_lib(tmp_path_factory):
    """The fixture library with TD_UID-5 flank 500-900 soft-masked (a repeat in the flank) and
    flank 0-40 replaced by an A-rich poly-A remnant (as real flanks start)."""
    d = tmp_path_factory.mktemp("rte_lib")
    lib = os.path.join(d, "rte_library")
    shutil.copytree(rte_sim.FIX, lib)
    fl = rte_sim.fasta("flanks_3p.fa")
    remnant = "AAAAAATAACAATAAAATGAGATAAAATCTACTAAAAATC"
    fl["TD_UID-5"] = remnant + FLANK[40:500] + FLANK[500:900].lower() + FLANK[900:]
    with open(os.path.join(lib, "flanks_3p.fa"), "w") as fh:
        for k, v in fl.items():
            fh.write(f">{k}\n{v}\n")
    return lib, fl["TD_UID-5"]


def run(insert, lib=rte_sim.FIX, **kw):
    inp, ev, genome, _ = rte_sim.build(insert, **kw)
    return RteAnnotator({"rte_library": lib}, genome=genome).annotate(inp, ev)


def sources(r):
    return [t for t in r.tags if t.startswith("TD3P_SOURCE=")]


@pytest.mark.parametrize("strand", [1, -1])
def test_orphan_from_masked_flank_is_not_sourced(masked_lib, strand):
    lib, flank = masked_lib
    r = run(flank[500:900] + A, lib=lib, strand=strand)
    assert sources(r) == []
    assert r.element != "ORPHAN_TD"


def test_orphan_from_unmasked_flank_still_sourced(masked_lib):
    lib, flank = masked_lib
    r = run(flank[1200:1600] + A, lib=lib)
    assert r.element == "ORPHAN_TD" and sources(r) == ["TD3P_SOURCE=TD_UID-5"]
    assert abs(r.detail["td_end"] - 1600) <= 5


def test_partnered_td_crossing_into_masked_part_keeps_source(masked_lib):
    """A real transduction may run into a flank repeat; the unmasked part names the source."""
    lib, flank = masked_lib
    r = run(L1[5600:CEND] + "A" * 8 + flank[100:600] + A, lib=lib)
    assert sources(r) == ["TD3P_SOURCE=TD_UID-5"]


def test_alu_insert_gets_no_l1_source():
    r = run(ALU + "A" * 8 + FLANK[20:420] + A)
    assert r.element == "ALU"
    assert sources(r) == []


def test_short_tail_hit_on_polya_remnant_is_not_a_source(masked_lib):
    """L1 3' end + tail + 26 bp of the A-rich flank start + tail (chr10:115802225 shape)."""
    lib, flank = masked_lib
    r = run(L1[5000:CEND] + "A" * 20 + flank[6:32] + "A" * 30, lib=lib)
    assert r.element == "L1"
    assert sources(r) == []
