# PEAR-TREE - tools/rte: 5'-inverted L1 (twin priming) whose only sequenced element piece is the
# INVERTED 5' part, plus a poly-A that the junction strings do not carry.
#
# Regression: PD37590 chr20:14211933-14211949 (GRCh38; TSD AAATTTTCCTCAATTT). Reference-forward
# the insert is  rc(L1HS 5478-5650) ... [unsequenced] ... A15 : the poly-A sits at the LEFT
# junction (element on +), so the 5' junction (RIGHT record) is ANTI-sense = INVERTED_5P. The
# combined.txt.gz LEFT junction string had no clip (no poly-A), so the strand fell back to the
# element pieces next to REF -- the inverted piece -- and flipped the frame: the "5' junction"
# became REF | poly-T and the call 5P_UNRESOLVED. Only a clip read showed A15 | REF.
#
# Run:  pytest test/test_rte_inverted5p.py
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rte_sim  # noqa: E402
from tools.rte.annotator import RteAnnotator, InsertionInput  # noqa: E402
from tools.rte.assembly import AssemblyResult, ReadLayout, Segment  # noqa: E402
from tools.rte.inputs import InsertionEvidence, JunctionEvidence, EvidenceRead  # noqa: E402
from tools.rte.sequtil import rc  # noqa: E402

LIB = rte_sim.library()
L1 = LIB.consensus["L1HS"]
CEND = LIB.cons_end["L1HS"]

TITLE = "chr20:14211933-14211949"
# combined.txt.gz junction strings as annotated (REF upper | clip lower): LEFT has NO clip
LEFT_JUNCTION = "AAATTTTCCTCAATTTATCCTAAGTTGACCATGTACCTACTTACTTCTCAGACCCTAAAAGATGCATAATACTGTTTAATTTATGAAATATTTACTAAATTACTGTCTGAATACATACACATGCAGAGGCACAAAC"
RIGHT_JUNCTION = (
    "CATTATGACCCAGAGGTGGTAATAGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTT"
    "tcttaatccagtctatcattgttggacatttgggttggttccaagtctttgctattgtgaatagtgccgcaataaacatacgtgtgcatgtgtctttatagcagcatgatttatactcatttgggtatatacccagtaatgggatggctgggtcaaatggtatttctagttc")
# evidence clip consensus of the LEFT junction (carries the 15 bp poly-A)
LEFT_CLIP_CONSENSUS = "aaaaaaaaaaaaaaaAAATTTTCCTCAATTTATCCTAAGTTGACCATGTACCTACTTACTTCTCAGACCCTAAAAGATGCATAATACTGTTTAATTTATGAAATATTTACTAAATTACTGTCTGAATACATACACATGCAGAGGCACAAAC"
# every evidence read combine kept (colony PD37590b_lo0028), allele-forward
READS = [
    ("LEFT", "CLIP", "7fae098f527ea5c6", "2",
     "AAAAAAAAAAAAAAAAAATTTTCCTCAATTTATCCTAAGTTGACCATGTACCTACTTACTTCTCAGACCCTAAAAGATGCATAATACTGTTTAATTTATGAAATATTTACTAAATTACTGTCTGAATACATACACATGCAGAGGCACAAAC"),
    ("LEFT", "DISC", "5faf3a1bced39348", "1",
     "ATTATATTAGCATAGGAATTTTACTGTATTTTATAAATTATTACTGTTCAAAAAATTATATTAATGTTCAGTTAAACACAGAGCATTTAAATGGCCTACAGGTAACACCACCAACATCTCAGTTCTAAAAATATGAAAAATGTCATTGAGC"),
    ("LEFT", "MATE", "5faf3a1bced39348", "2",
     "AAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCACCCAAAATAACCCCAAATTCAACAACAACACATACACAACAAACCCCAAAAACAAGAAAAAATCATTTTTATTTTTTTAAAATTTA"),
    ("LEFT", "MATE", "7fae098f527ea5c6", "1",
     "AAAAAGCAGCAGCAAGATTAGACCTTAATTCCACTTGATCACTTGTTCATCAACAAATGGTCTTCTCACAAATTCCAAGACGCCTTAAGAAGAGATGCCCAGCCGGTAGCATTAACATAACTTTATCTAAATTATTTCGCCTTCTCTCAGT"),
    ("RIGHT", "CLIP", "d5e5ef671214ffd4", "1",
     "CATTATGACCCAGAGGTGGTAATAGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGG"),
    ("RIGHT", "CLIP", "b44ad24686737ae7", "2",
     "TTATGACCCAGAGGTGGTAATAGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTT"),
    ("RIGHT", "CLIP", "77e4f38835071ceb", "1",
     "ATGACCCAGAGGTGGTAATAGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTGTACATTTGGGTTGG"),
    ("RIGHT", "CLIP", "978a2c7623c2a294", "1",
     "TGACCCAGAGGTGGTAATAGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGT"),
    ("RIGHT", "CLIP", "8e845770f1902b4b", "1",
     "CCCAGAGGTGGTAATAGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCC"),
    ("RIGHT", "CLIP", "baba56f77b8cd2d4", "2",
     "AGAGGTGGTAATAGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAANAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAG"),
    ("RIGHT", "CLIP", "8b76a888001246a8", "1",
     "AGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTG"),
    ("RIGHT", "CLIP", "8bfcbf0a13eec460", "2",
     "GTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGC"),
    ("RIGHT", "CLIP", "d45f6df9d077f1f2", "1",
     "GCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATACAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCG"),
    ("RIGHT", "CLIP", "76e1287e786b194d", "1",
     "GCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATACAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCG"),
    ("RIGHT", "DISC", "c01b91f1bf583b7e", "1",
     "CTCCAGAATCTACTCAGGGAAAATGCTAGATAGCTTTCATTTATACCTCCAGATTCACTGTCCACCTTTCTCCATCTAGCTATGTGTCCAGGAGTCTGATAGTATAGATTGTAAGTATGCAGTCTGTGTTTACAGGCTCTCTGCTATCTGC"),
    ("RIGHT", "DISC", "90a3da2fed1e00bd", "1",
     "CTCCAGAATCTACTCAGGGAAAATGCTAGATAGCTTTCATTTATACCTCCAGATTCACTGTCCACCTTTCTCCATCTAGCTATGTGTCCAGGAGTCTGATAGTATAGATTGTAAGTATGCAGTCTGTGTTTACAGGCTCTCTGCTATCTGC"),
    ("RIGHT", "DISC", "f2e4ef81e8966e71", "2",
     "GTGAGAGAGGTATTTAATTCACCTAGTTCTCTTCATGAGTGGTCACCGTAGTGTGTCTCTAGACTGAAGGTCTCAGCATGTGAAAGGTACCACCCCTCACAAGGCCTCTTAGTCTTCAGCTTTCAGTAACAGCTCCCTTCCTTTACCCATT"),
    ("RIGHT", "DISC", "eb3074cd030b565c", "1",
     "ATGAGTGGTCACCGTAGTGTGTCTCTAGACTGAAGGTCTCAGCATGTGAAAGGTACCACCCCTCACAAGGCCTCTTAGTCTTCAGCTTTCAGTAACAGCTCCCTTCCTTTACCCATTATGACCCAGAGGTGGTAATAGTCACTACAGTGCT"),
    ("RIGHT", "DISC", "d1e89c35a7165805", "1",
     "CCCCTCACAAGGCCTCTTAGTCTTCAGCTTTCAGTAACAGCTCCCTTCCTTTACCCATTATGACCCAGAGGTGGTAATAGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCT"),
    ("RIGHT", "MATE", "90a3da2fed1e00bd", "2",
     "TATCAAGCTGTCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCNAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACG"),
    ("RIGHT", "MATE", "f2e4ef81e8966e71", "1",
     "TTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTGTCTTTATAGCAGCATGATTTATACTCATTTGGGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATTGCCACACTGA"),
    ("RIGHT", "MATE", "c01b91f1bf583b7e", "2",
     "TATCAAGCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCNAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACG"),
    ("RIGHT", "MATE", "d1e89c35a7165805", "2",
     "CAATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTGTCTTTATAGCAGCATGATTTATACTCATTTGGGTATATACCCAGTAATGGGAT"),
    ("RIGHT", "MATE", "eb3074cd030b565c", "2",
     "TCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTGTCTTTATAGCAGCATGATTTATACTCATTTGGGTATATACCCAGTAATGGGATGGCTGGGTCAAA"),
    ("RIGHT", "MATE", "d45f6df9d077f1f2", "2",
     "TTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTGTCTTTATAGCAGCATGATTTATACTCATTTGGGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTTCTAGTTC"),
    ("RIGHT", "MATE", "76e1287e786b194d", "2",
     "TTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTGTCTTTATAGCAGCATGATTTATACTCATTTGGGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTTCTAGTTC"),
    ("RIGHT", "MATE", "978a2c7623c2a294", "2",
     "TGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTNGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTGTCTTT"),
    ("RIGHT", "MATE", "8b76a888001246a8", "2",
     "GCTATCCTTTGTGGTTTCTCTACACCTTAGCAAATAGTGCCTTTATTAAATTTTCCTCAATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCA"),
    ("RIGHT", "MATE", "77e4f38835071ceb", "2",
     "CTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACANACGTGTGCATGTGTCTTTATAGCAGCATGATTTATACTCATTTGGGTATATACCCAGTAATGGGATGGCTGGG"),
    ("RIGHT", "MATE", "8e845770f1902b4b", "2",
     "TCTCCATCTAGCTATGTGTCCAGGAGTCTGATAGTATAGATTGTAAGTATGCAGTCTGTGTTTACAGGCTCTCNGCTATCTGCTTTCCCATTAGCTCTGGCCGATGGAAGTCACTGGCGAGATATTGGAGTGAGGCAGCAGAGTGAGAGAG"),
    ("RIGHT", "MATE", "d5e5ef671214ffd4", "2",
     "CAGCAGAGTGAGAGAGGTATTTAATTCACCTAGTTCTCTTCATGAGTGGTCACCGTAGTGTGTCTCTAGACTGNAGGTCTCAGCATGTGAAAGGTACCACCCCTCACAAGGCCTCTTAGTCTTCAGCTTTCAGTAACAGCTCCCTTCCTTT"),
    ("RIGHT", "MATE", "b44ad24686737ae7", "1",
     "GTTATTTAATTCACCTAGTTCTCTTCATGAGTGGTCACCGTAGTGTGTCTCTAGACTGAAGGTCTCAGCATGTGAAAGGTACCACCCCTCACAAGGCCTCTGAGTCTTCAGCTTTCAGTAACAGCTCCCTTCCTTTACCCATTATGACCCA"),
    ("RIGHT", "MATE", "8bfcbf0a13eec460", "1",
     "GTGGTCACCGTAGTGTGTCTCTAGACTGAAGGTCTCAGCATGTGAAAGGTACCACCCCTCACAAGGCCTCTTAGTCTTCAGCTTTCAGTAACAGCTCCCTTCCTTTACCCATTATGACCCAGAGGTGGTAATAGTCACTACAGTGCTAGCT"),
    ("RIGHT", "MATE", "baba56f77b8cd2d4", "1",
     "CTCTAGACTGAAGGTCTCAGCATGTGAAAGGTACCACCCCTCACAAGGCCTCTTAGTCTTCAGCTTTCAGTAACAGCTCCCTTCCTTTACCCATTATGACCCAGAGGTGGTAATAGTCACTACAGTGCTAGCTTCAGGTTATCAAGCTATC"),]


# ---- PD37590 chr2:126577464-126577479 (all reads PD37590b_lo0117): minimal twin priming.
# RIGHT record: REF | ATTATTAT T10 rc(L1HS 5993-6019) = - element, 3' terminus + tail;
# LEFT record: A19 | REF = the inverted copy of the tail at the 5' junction. Was POLYA_ONLY:
# 26 element bp < min_element_bp and the strand came from the bare tails.
CHR2_TITLE = "chr2:126577464-126577479"
CHR2_LEFT_JUNCTION = "AAAAATTTCAATGTTTATCCTAAAATGTTAAAATAAATGTTAGATTTAATTTTTAAAAGAATTGGATACTATTTAAAAAAAAAGAAGAAGAGATAATAAGAACAAGTTGTTGGGAAAAGAAGAAAAAGATTG"
CHR2_RIGHT_JUNCTION = "TCAATGGTTGTCTCCCATTCACCTTGAATCAAACCCTGCCAGAAACAGGGGCAATGATACACAAACAACCTTCCTTGCCAGGAAAAAGACTTATCTACCATCAGGGGACAAAAGAGAACTCAAAAATTTCAATGTTattattattttttttttattatactctaag"
CHR2_LEFT_CLIP = "aaaaaaaaaaaaaaaaaaaAAAAATTTCAATGTTTATCCTAAAATGTTAAAATAAATGTTAGATTTAATTTTTAAAAGAATTGGATACTATTTAAAAAAAAAGAAGAAGAGATAATAAGAACAAGTTGTTGGGAAAAGAAGAAAAAGATTG"
CHR2_RIGHT_CLIP = "TCAATGGTTGTCTCCCATTCACCTTGAATCAAACCCTGCCAGAAACAGGGGCAATGATACACAAACAACCTTCCTTGCCAGGAAAAAGACTTATCTACCATCAGGGGACAAAAGAGAACTCAAAAATTTCAATGTTattattattttttttttattatactctaagttttagggtacat"
CHR2_READS = [
    ("LEFT", "CLIP", "e7b4601cd8e0d53a", "1",
     "AAAAAAAAAAAAAAAAAAAAAAAATTTCAATGTTTATCCTAAAATGTTAAAATAAATGTTAGATTTAATTTTTAAAAGAATTGGATACTATTTAAAAAAAAAGAAGAAGAGATAATAAGAACAAGTTGTTGGGAAAAGAAGAAAAAGATTG"),
    ("LEFT", "MATE", "e7b4601cd8e0d53a", "2",
     "AGAGAAAGAAAAAAGAAACTCCTTAGAAACTCCTGAGAATCAGACTGATACCATATTTCTCTTTAGCAACACTTATGTAAGACAGAATGGTGCAACAGCACCAAGTTCTAAGGGAAATTATTTCAATCCTAGAGTTGATACCAAGCCAAAA"),
    ("RIGHT", "CLIP", "e9fccf786a6c6d89", "2",
     "TCAATGGTTGTCTCCCATTCACCTTGAATCAAACCCTGCCAGAAACAGGGGCAATGATACACAAACAACCTTCCTTGCCAGGAAAAAGACTTATCTACCATCAGGGGACAAAAGAGAACTCAAAAATTTCAATGTTATTATTATTTTTTTT"),
    ("RIGHT", "CLIP", "13eec6213774460a", "2",
     "CATTCACCTTGAATCAAACCCTGCCAGAAACAGGGGCAATGATACACAAACAACCTTCCTTGCCAGGAAAAAGACTTATCTACCATCAGGGGACAAAAGAGAACTCAAAAATTTCAATGTTATTATTATTTTTTTTTTATTATACTCTAAG"),
    ("RIGHT", "CLIP", "2478c262e62b21c4", "2",
     "TCAAACCCTGCCAGAAACAGGGGCAATGATACACAAACAACCTTCCTTGCCAGGAAAAAGACTTATCTACCATCAGGGGACAAAAGAGAACTCAAAAATTTCAATGTTATTATTATTTTTTTTTTATTATACTCTAAGTTTTAGGGTACAT"),
    ("RIGHT", "MATE", "13eec6213774460a", "1",
     "CAGGTAAAAACACACATATGAACAGCAGAAAGGAGCTGAGGCTTCTAGGGTGATTGCTGGAGCCATGGATTATAGCACAGACAACCTTGGCAATGGGGAAGGGGCTGCAGTGGATGAGAGAAGGTGATGAGTCAATGGTTGTCTCCCATTC"),
    ("RIGHT", "MATE", "e9fccf786a6c6d89", "1",
     "CAGCAGAAAGGAGCTGAGGCTTCTAGGGTGATTGCTGGAGCCATGGATTATAGCACAGACAACCTTGGCAATGGGGAAGGGGCTGCAGTGGATGAGAGAAGGTGATGAGTCAATGGTTGTCTCCCATTCACCTTGAATCAAACCCTGCCAG"),
    ("RIGHT", "MATE", "2478c262e62b21c4", "1",
     "GGCTGCAGTGGATGAGAGAAGGTGATGAGTCAATGGTTGTCTCCCATTCACCTTGAATCAAACCCTGCCAGAAACAGGGGCAATGATACACAAACAACCTTCCTTGCCAGGAAAAAGACTTATCTACCATCAGGGGACAAAAGAGAACTCA"),
]

def _locus(with_polya_clip):
    reads = [EvidenceRead(side, role, "PD37590b_lo0028", frag, r12, seq)
             for side, role, frag, r12, seq in READS]
    ev = InsertionEvidence(TITLE, reads=reads)
    if with_polya_clip:
        ev.junctions["LEFT"] = JunctionEvidence("LEFT", clip_consensus=LEFT_CLIP_CONSENSUS,
                                                supported=1, n_samples=1)
    ann = RteAnnotator({"rte_library": rte_sim.FIX})
    return ann.annotate(InsertionInput(TITLE, LEFT_JUNCTION, RIGHT_JUNCTION), ev)


@pytest.mark.parametrize("with_polya_clip", [False, True])
def test_pd37590_chr20_inverted_5p(with_polya_clip):
    r = _locus(with_polya_clip)
    assert r.element == "L1"
    assert r.strand == 1                  # poly-A at the LEFT junction: + element
    assert r.structure == "INVERTED_5P"
    lo, p2 = (int(x) for x in r.detail["inv"].split("-"))
    assert abs(p2 - 5650) <= 3            # inverted piece ends at L1HS ~5650 at the flank
    # the forward part was never sequenced: no inversion point, no twin-priming points
    assert r.detail["inv_junction"] == "unresolved"
    assert "fwd_start" not in r.detail and r.score_input.inv_p1 is None
    if with_polya_clip:
        assert r.polya_len >= 10
    else:                                 # tail only in the clip read: reported, not scored
        assert r.polya_len == 0 and r.detail["polya_reads"] >= 10


def _strip_clip(inp, ev, side):
    """Drop the clip (lower case) of one junction string, as combine delivered the regression
    locus: the poly-A is then only in the clip reads."""
    s = getattr(inp, "left_seq" if side == "LEFT" else "right_seq")
    s = "".join(ch for ch in s if ch.isupper())
    setattr(inp, "left_seq" if side == "LEFT" else "right_seq", s)
    ev.junctions[side].clip_consensus = s
    ev.junctions[side].polya_len_median = 0.0


def _sim(insert, strand, strip_polya):
    inp, ev, genome, _ = rte_sim.build(insert, strand=strand)
    if strip_polya:
        _strip_clip(inp, ev, "LEFT" if strand > 0 else "RIGHT")
    ann = RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome)
    return ann.annotate(inp, ev)


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("strip_polya", [False, True])
@pytest.mark.parametrize("forward", [False, True])
def test_inverted_piece_only(strand, strip_polya, forward):
    """Twin priming whose forward part is not sequenced. forward=False: the insert is only
    rc(L1[p2-400:p2]) + poly-A (the regression locus' shape); forward=True: a 1.5 kb forward
    part whose reads are dropped except those at the poly-A junction."""
    p2 = CEND - 1500
    fwd = L1[p2 + 10:CEND] if forward else ""
    insert = rc(L1[p2 - 400:p2]) + fwd + "A" * 25
    excl = range(420, len(insert) - 25 - 60, 20) if forward else ()
    inp, ev, genome, _ = rte_sim.build(insert, strand=strand, exclude=excl)
    if strip_polya:
        _strip_clip(inp, ev, "LEFT" if strand > 0 else "RIGHT")
    r = RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome).annotate(inp, ev)
    assert r.element == "L1"
    assert r.strand == strand
    assert r.structure == "INVERTED_5P"
    assert abs(int(r.detail["inv"].split("-")[1]) - p2) <= 5
    if not forward:
        assert r.detail["inv_junction"] == "unresolved"


@pytest.mark.parametrize("strand", [1, -1])
def test_truncated_without_junction_polya_unchanged(strand):
    r = _sim(L1[5500:CEND] + "A" * 25, strand, strip_polya=True)
    assert r.element == "L1" and r.strand == strand
    assert r.structure == "TRUNCATED_5P"
    assert abs(r.detail["j5"] - 5500) <= 5


def _lay(*segs, frag=("S", "f")):
    return ReadLayout("r", "", "CLIP", frag, "N" * 200, [Segment(*s) for s in segs])


def test_polya_strand_vote_is_unanimous():
    a_ref = _lay((0, 15, "POLYA", "A", -1, -1, 1), (15, 150, "REF", "site"), frag=("S", "1"))
    ref_t = _lay((0, 100, "REF", "site"), (100, 120, "POLYA", "T", -1, -1, -1), frag=("S", "2"))
    ref_a = _lay((0, 100, "REF", "site"), (100, 120, "POLYA", "A", -1, -1, 1), frag=("S", "3"))
    assert AssemblyResult._strand_from_polya([a_ref])[0] == 1
    assert AssemblyResult._strand_from_polya([ref_t])[0] == -1
    # tails on both sides (slippage / ligation): no poly-A strand
    assert AssemblyResult._strand_from_polya([a_ref, ref_t])[0] == 0
    # an A-run AFTER the reference (REF | A) is not a tail of a + element at the LEFT junction
    assert AssemblyResult._strand_from_polya([ref_a])[0] == 0


def _replay(title, lj, rj, reads, lclip=None, rclip=None):
    ev = InsertionEvidence(title, reads=[EvidenceRead(sd, ro, "S", fr, r12, seq)
                                         for sd, ro, fr, r12, seq in reads])
    for side, c in (("LEFT", lclip), ("RIGHT", rclip)):
        if c:
            ev.junctions[side] = JunctionEvidence(side, clip_consensus=c, supported=1, n_samples=1)
    return RteAnnotator({"rte_library": rte_sim.FIX}).annotate(InsertionInput(title, lj, rj), ev)


@pytest.mark.parametrize("with_clips", [False, True])
def test_pd37590_chr2_minimal_inversion(with_clips):
    r = _replay(CHR2_TITLE, CHR2_LEFT_JUNCTION, CHR2_RIGHT_JUNCTION, CHR2_READS,
                CHR2_LEFT_CLIP if with_clips else None, CHR2_RIGHT_CLIP if with_clips else None)
    assert r.element == "L1"              # 26 bp at the consensus terminus next to the tail
    assert r.strand == -1                 # the terminus + tail is at the RIGHT record
    assert r.structure == "INVERTED_5P"
    assert r.detail["inv"] == "polyA" and r.detail["inv_junction"] == "unresolved"
    assert r.detail["j3"] >= CEND - 5


@pytest.mark.parametrize("strand", [1, -1])
def test_minimal_inversion_synthetic(strand):
    """element sense: REF | T19 (inverted tail) | [300 bp never sequenced] | L1 last 40 bp |
    A12 | REF -- reads see either the 5' REF | T19 or the 3' terminus + tail, never both."""
    import random
    mid = rte_sim.rnd(300, random.Random(7))
    inp, ev, genome, _ = rte_sim.build("T" * 19 + mid + L1[CEND - 40:CEND] + "A" * 12,
                                       strand=strand, exclude=range(45, 300, 10), clip_len=30,
                                       step=10)
    r = RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome).annotate(inp, ev)
    assert r.element == "L1" and r.strand == strand
    assert r.structure == "INVERTED_5P" and r.detail["inv"] == "polyA"
    assert r.detail["j3"] >= CEND - 5


@pytest.mark.parametrize("strand", [1, -1])
def test_tails_on_both_sides_without_terminus_stay_polya_only(strand):
    """Slippage-like: poly-T | junk | poly-A with no element terminus -> no element."""
    r = _sim("T" * 19 + "GCATGCTAGCATCG" + "A" * 20, strand, strip_polya=False)
    assert r.element in ("POLYA_ONLY", "UNKNOWN")
    assert r.structure == "5P_UNRESOLVED"


def test_non_terminal_piece_is_not_a_minimal_inversion():
    """26 bp of L1 from the middle of the consensus before the tail is no 3' terminus: the
    minimal-inversion rule (inv=polyA) must not fire."""
    r = _sim("T" * 19 + L1[3000:3026] + "A" * 12, 1, strip_polya=False)
    assert r.detail.get("inv") != "polyA"
