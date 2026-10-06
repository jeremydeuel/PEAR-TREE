# PEAR-TREE - final locus class (annotate_v2 table `class`): a confident tools/rte verdict
# overrides the legacy fallback classes (unknown / artefact / templated_insertion); precedence
# in tools/rte/locus_class.py.
#
# Real PD37590 loci, replayed through RteAnnotator with the REAL resources/rte_library and the
# hg38 reads/junction strings combine kept (test/fixtures/rte_regression: reads, junction
# strings, +-2.5 kb hg38 windows around each site):
#   chr7:131206236  legacy templated_insertion "(source chrX:11296797)" (hs1) -- an orphan 3'
#                   transduction from L1_chrX_11289796_f (Xp22.2-1, UID-50): the left clip is
#                   flank offset 971-1287 at 100 % (unmasked) -> LINE1
#   chr7:87178413   legacy unknown "(unmapped complex insertion)" -- L1HS 5'-inverted, TPRT -> LINE1
#   chr2:99279093   legacy templated_insertion, ALU with a (spurious, repeat-in-flank) L1 source
#                   -> no source any more, low score -> stays templated_insertion
#
# Run:  pytest test/test_rte_locus_class.py
import csv
import gzip
import os
import sys
from collections import defaultdict
from functools import lru_cache

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rte_sim  # noqa: E402,F401  (puts the repo on sys.path)
from tools.rte.annotator import RteAnnotator, InsertionInput  # noqa: E402
from tools.rte.genome import FastaGenome  # noqa: E402
from tools.rte.inputs import InsertionEvidence, EvidenceRead  # noqa: E402
from tools.rte.locus_class import locus_class, flank_source_at  # noqa: E402
from tools.rte.record import RteRecord  # noqa: E402

FIX = os.path.join(HERE, "fixtures", "rte_regression")
LIB_REAL = os.path.join(rte_sim.REPO, "resources", "rte_library")


def _tsv(name):
    with gzip.open(os.path.join(FIX, name), "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@lru_cache(None)
def _annotator():
    genome = FastaGenome(os.path.join(FIX, "pd37590_hg38_windows.fa.gz"))
    return RteAnnotator({"rte_library": LIB_REAL, "rte_wide_window": 2000}, genome=genome)


@lru_cache(None)
def _loci():
    reads = defaultdict(list)
    for r in _tsv("pd37590_reads.tsv.gz"):
        reads[r["locus"]].append(EvidenceRead(r["side"], r["role"], r["colony"], r["fragment"],
                                              r["read"], r["sequence"]))
    return {r["locus"]: (r, reads[r["locus"]]) for r in _tsv("pd37590_loci.tsv.gz")}


@lru_cache(None)
def replay(locus):
    row, reads = _loci()[locus]
    ann = _annotator()
    rec = ann.annotate(InsertionInput(locus, row["L_junction"], row["R_junction"]),
                       InsertionEvidence(locus, reads=list(reads)))
    return row, rec


def final(locus):
    import re
    row, rec = replay(locus)
    m = re.search(r"\(source ([^:\s()]+):(\d+)\)", row["conclusion"])
    tsrc = (m.group(1), int(m.group(2))) if m else None
    return locus_class(row["legacy_class"], rec, _annotator().lib, tsrc)


def test_chr7_131206236_orphan_transduction_is_line1():
    row, rec = replay("chr7:131206236-131206252")
    assert row["legacy_class"] == "templated_insertion"
    assert rec.element == "ORPHAN_TD"
    assert "TD3P_SOURCE=L1_chrX_11289796_f" in rec.tags
    assert abs(int(rec.detail["td_end"]) - 1308) <= 10
    cls, note = final("chr7:131206236-131206252")
    assert cls == "LINE1"
    assert "3' transduction from L1_chrX_11289796_f (Xp22.2-1, UID-50)" in note


def test_chr7_87178413_inverted_l1_is_line1():
    row, rec = replay("chr7:87178413-87178423")
    assert row["legacy_class"] == "unknown"
    assert rec.element == "L1" and rec.structure == "INVERTED_5P"
    cls, note = final("chr7:87178413-87178423")
    assert cls == "LINE1"
    assert note.startswith("L1HS, 5'-inverted (twin priming), inv 5155-")


def test_chr2_99279093_spurious_source_not_promoted():
    row, rec = replay("chr2:99279093-99279092")
    assert not any(t.startswith("TD3P_SOURCE=") for t in rec.tags)
    assert final("chr2:99279093-99279092")[0] == "templated_insertion"


def test_chr10_115802225_no_source_from_32bp_of_flank():
    row, rec = replay("chr10:115802225-115802229")
    assert rec.element == "L1"
    assert not any(t.startswith("TD3P_SOURCE=") for t in rec.tags)
    assert final("chr10:115802225-115802229")[0] == "LINE1"     # legacy class kept


@pytest.mark.parametrize("locus,structure", [("chr20:14211933-14211949", "INVERTED_5P"),
                                             ("chr2:126577464-126577479", "INVERTED_5P")])
def test_inverted_loci_with_hg38_windows(locus, structure):
    row, rec = replay(locus)
    assert rec.element == "L1" and rec.structure == structure


# ------------------------------------------------------------------ precedence (unit level)
def _rec(element, call="TPRT", tags=(), **detail):
    r = RteRecord("x", element=element, tprt_call=call, tags=list(tags), consensus="L1HS")
    r.detail.update(detail)
    return r


def test_positive_legacy_classes_are_kept():
    for legacy in ("ALU", "SVA", "LINE1", "processed_pseudogene", "RTE_other", "non_RTE_SV",
                   "microsatellite"):
        assert locus_class(legacy, _rec("L1"), None)[0] == legacy


def test_unconfident_rte_does_not_override():
    for call in ("UNCERTAIN", "ARTEFACT_LIKE"):
        assert locus_class("artefact", _rec("ALU", call), None)[0] == "artefact"


def test_orphan_td_needs_credible_source():
    lib = _annotator().lib
    known = _rec("ORPHAN_TD", tags=("TD3P", "TD3P_SOURCE=L1_chrX_11289796_f"), td_end=1308)
    assert locus_class("templated_insertion", known, lib)[0] == "LINE1"
    novel_b = _rec("ORPHAN_TD", tags=("TD3P", "TD3P_SOURCE=novel:chr2:70428272", "NOVEL_SOURCE"),
                   novel_tier="B")
    assert locus_class("templated_insertion", novel_b, lib)[0] == "templated_insertion"
    novel_a = _rec("ORPHAN_TD", tags=("TD3P", "TD3P_SOURCE=novel:chrX:1-2", "NOVEL_SOURCE"),
                   novel_tier="A")
    assert locus_class("unknown", novel_a, lib)[0] == "LINE1"
    unsourced = _rec("ORPHAN_TD", tags=("TD3P",))
    assert locus_class("unknown", unsourced, lib)[0] == "unknown"


def test_pseudogene_needs_exon_junction():
    assert locus_class("unknown", _rec("PSEUDOGENE", tags=("EXON_JUNCTION",)), None)[0] == \
        "processed_pseudogene"
    assert locus_class("unknown", _rec("PSEUDOGENE"), None)[0] == "unknown"


def test_templated_source_in_unmasked_l1_flank_without_rte_record():
    lib = _annotator().lib
    # hs1 chrX:11296797 = 971 bp into the UID-50 3' flank (source 3' end hs1 chrX:11295826, +)
    assert flank_source_at(lib, "chrX", 11296797) == ("L1_chrX_11289796_f", 971)
    cls, note = locus_class("templated_insertion", None, lib, ("chrX", 11296797))
    assert cls == "LINE1" and "L1_chrX_11289796_f (Xp22.2-1, UID-50)" in note
    # upstream of the source / beyond 15 kb: no transduction
    assert flank_source_at(lib, "chrX", 11289000) is None
    assert flank_source_at(lib, "chrX", 11295826 + 15500) is None


# ------------------------------------------------------------------ ALU vs SVA vs reference
def test_chr11_127257706_reference_alu_is_not_the_insert():
    """Legacy unknown; the RTE record said ALU (LIKELY_TPRT) and was suspected SVA. The only
    'element' pieces (ALU_Y 198-276 at 80-88 %) are 100 % identical to the reference Alu right
    behind the T23 at the breakpoint (reads through a poly-T of different length lose their
    whole-read REF hit); the clip after the poly-T is degraded and closer to that local Alu
    than to any Alu (<= 0.69) or SVA (<= 0.64) consensus, with no hexamer / VNTR / SINE-R. The
    honest output is no element (not a confident ALU, not SVA) and no class promotion."""
    row, rec = replay("chr11:127257706-127257729")
    assert rec.element not in ("ALU", "SVA")
    assert final("chr11:127257706-127257729")[0] == "unknown"


def test_element_piece_identical_to_site_window_is_ref():
    import random
    from tools.rte.assembly import Assembler, SiteContext
    from tools.rte.sequtil import rc
    lib = rte_sim.library()
    alu = lib.consensus["ALU_Y"][:lib.cons_end["ALU_Y"]]
    rng = random.Random(3)
    old = "".join(ch if rng.random() > 0.12 else rng.choice("ACGT".replace(ch, "")) for ch in alu)
    left = rte_sim.rnd(400, rng)
    window = left + "T" * 23 + rc(old) + rte_sim.rnd(400, rng)
    ctx = SiteContext("chrS:400-423", "chrS", 400, 423, window_start=0, window_seq=window)
    asm = Assembler(lib)
    # read through a poly-T expanded from 23 to 35: the (low-quality, 25 % errors) flank before
    # the run + T35 + the reference Alu -- the whole-read REF hit falls below min_ref_identity
    noisy = "".join(ch if rng.random() > 0.25 else rng.choice("ACGT".replace(ch, ""))
                    for ch in left[-40:])
    read = noisy + "T" * 35 + rc(old)[:76]
    lay = asm.layout(read, ctx, side="LEFT", role="CLIP")
    assert not any(s.kind == "ELEMENT" for s in lay.segments)
    assert lay.segments[-1].kind == "REF"
    # a young inserted Alu (consensus-identical) next to that old reference copy stays ELEMENT
    read2 = left[-60:] + alu[150:] + "A" * 20
    lay2 = asm.layout(read2, ctx, side="RIGHT", role="CLIP")
    assert any(s.kind == "ELEMENT" for s in lay2.segments)
