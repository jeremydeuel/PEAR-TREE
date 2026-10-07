# PEAR-TREE - structural-variant annotation regression test
#
# Locks in the SV extension of tools/annotate_v2.py (Insertion._sv_subtype / _sv_partner_maps /
# link_reciprocal_translocations). A rearrangement presents as a junction whose clipped side is
# reference sequence from the partner breakpoint rather than an inserted element, so annotate
# compares the junction's own locus (from the title, "contig:start-end") to where the clip
# uniquely remaps:
#   * other contig                       -> translocation junction
#   * same contig, distal, reverse clip  -> inversion
#   * same contig, distal, forward clip  -> intrachromosomal SV (deletion/duplication)
# Crucially SV and MEI are NOT mutually exclusive: the SV note is ADDITIVE to the Alu/L1/SVA
# call, never a substitute, and the class stays the element's. Specificity comes from requiring
# the partner clip to map UNIQUELY (MAPQ >= sv_min_mapq) so dispersed repeats (MAPQ 0) do not
# read as translocations, plus a poly-A transduction guard.
#
# Coordinates are the real ones seen in analysis/mei9x10 (the legacy non_RTE_SV flag rows), so
# the distance/contig/strand arithmetic runs on genuine values. Pure conclusion()-level logic,
# so no nhmmscan / bowtie2 / pysam runtime is needed.
#
# Run with:  python test/test_annotate_sv.py     (from the repo root)
#        or:  pytest test/test_annotate_sv.py

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
Dfam_Annotation = annotate_v2.Dfam_Annotation
RepeatMasker_Annotation = annotate_v2.RepeatMasker_Annotation
VAC = annotate_v2.VariantAnnotationContainer


def _dfam(model, bits, strand="-"):
    line = [model, "ACC", "target", str(bits), "1e-10", "0.0",
            "1", "100", strand, "1", "100", "1", "100", "1400", "desc"]
    return Dfam_Annotation(line)


def _rmsk_rte(name="AluY", cls="SINE/Alu"):
    """A RepeatMasker_Annotation whose repClass is an RTE (so a partner carrying it is a
    repeat-mediated SV partner rather than a clean unique locus)."""
    line = ["100", "0", "0", "0", "chr1", "1", "300", "300", "+", name, cls, "1", "300", "1"]
    return RepeatMasker_Annotation(line)


def _map(pos, mapq, strand, rmsks=()):
    return (pos, mapq, list(rmsks), strand)


# Default junction sequences carry a long, high-complexity clip on each side so that the SV
# clip-map trust guard (_clip_trustworthy_for_sv: len>=25, entropy>=1.9, AT<0.72) is satisfied —
# a real translocation/inversion partner comes from a clip that maps uniquely on merit. The two
# clips and the flank are mutually unrelated, so the reference-free duplication/STR detector does
# NOT fire on these fixtures (they exercise the distal-map SV path).
_LCLIP = "gcattgacctaggtcaagcttgaccatgca"     # 30 bp, complex, AT~0.5
_RCLIP = "ttcaggcatactgacatctgagctacatgc"     # 30 bp, complex, unrelated to _LCLIP/flank
_UFLANK = "TGCAGTCAGTCTGACAGTCA"              # 20 bp reference flank (upper)


def _make(title, left_seq=_LCLIP + _UFLANK, right_seq=_UFLANK + _RCLIP,
          left_dfams=(), right_dfams=(), left_maps=(), right_maps=()):
    ins = Insertion(title, left_seq, right_seq)
    ins.left_dfams = list(left_dfams)
    ins.right_dfams = list(right_dfams)
    ins.left_maps = list(left_maps)
    ins.right_maps = list(right_maps)
    return ins


# poly-A seqs (>= polya_min_len A's flush to the aligned part) for the transduction cases.
LEFT_POLYA = "a" * 20 + "ACGTACGT"       # a{20}[ACGT] -> has_left_polyA


# --- pure structural variants (no element call) -------------------------------------------

# chr8:14978384 -> both clips map uniquely to chr6 (a different contig): translocation junction.
TRANSLOCATION = _make("chr8:14978384-14978386",
                      left_maps=[_map("chr6:43826342+", 60, "+")],
                      right_maps=[_map("chr6:43826352+", 60, "+")])

# chr11:24263179 -> chr11:24390378(+) and chr11:24390496(-): same contig, distal, a reverse
# clip => inversion (the reverse-strand junction is the inversion signature).
INVERSION = _make("chr11:24263179-24263186",
                  left_maps=[_map("chr11:24390378+", 60, "+")],
                  right_maps=[_map("chr11:24390496-", 60, "-")])

# chr8:102238952 -> chr8:103364844(+)/103364878(+): same contig, distal, co-oriented =>
# deletion/duplication (an intrachromosomal SV that is not an inversion).
DELDUP = _make("chr8:102238952-102238971",
               left_maps=[_map("chr8:103364844+", 60, "+")],
               right_maps=[_map("chr8:103364878+", 60, "+")])

# a unique non-RTE partner only ~190 bp away is local context, not a rearrangement: the legacy
# flag still fires (nothing lost) but no specific subtype is claimed.
LOCAL = _make("chr5:1000-1010", left_maps=[_map("chr5:1200+", 60, "+")])

# a dispersed repeat clip maps at MAPQ 0 (multi-maps): must NOT read as a translocation.
LOWMAPQ = _make("chr1:5000-5010", left_maps=[_map("chr2:99999+", 0, "+")])


# --- SV and MEI are not mutually exclusive: additive note on an element call ---------------

# An Alu call (right AluY + left poly-A) that ALSO has a unique reverse-strand same-contig
# partner: the Alu class is kept and an inversion SV note is appended.
ALU_PLUS_INVERSION = _make("chr1:5000-5010",
                           left_seq=LEFT_POLYA,
                           right_dfams=[_dfam("AluY", 60.0, "+")],
                           right_maps=[_map("chr1:200000-", 60, "-")])

# Same Alu call but the extra partner is co-oriented (forward) and there is a poly-A hallmark:
# the transduction guard suppresses the co-oriented SV note (a 3' transduction look-alike).
ALU_TRANSDUCTION_GUARD = _make("chr1:5000-5010",
                               left_seq=LEFT_POLYA,
                               right_dfams=[_dfam("AluY", 60.0, "+")],
                               right_maps=[_map("chr1:200000+", 60, "+")])

# A partner clip that maps only to a paralogous repeat copy (RTE-annotated) is NOT reported on
# an element call without reciprocal confirmation: it is far more likely the Alu itself aligning
# to a family member than a real rearrangement, so it must not decorate an ordinary Alu call.
ALU_REPEAT_PARTNER = _make("chr1:5000-5010",
                           left_seq=LEFT_POLYA,
                           right_dfams=[_dfam("AluY", 60.0, "+")],
                           right_maps=[_map("chr1:200000-", 60, "-", rmsks=[_rmsk_rte()])])


def test_translocation():
    c = TRANSLOCATION.conclusion()
    assert "translocation" in c, c
    assert VAC.element_class(c) == "non_RTE_SV", c


def test_inversion():
    c = INVERSION.conclusion()
    assert "inversion" in c, c
    assert VAC.element_class(c) == "non_RTE_SV", c


def test_deletion_duplication():
    c = DELDUP.conclusion()
    assert "deletion/duplication" in c, c
    assert "inversion" not in c, c
    assert VAC.element_class(c) == "non_RTE_SV", c


def test_local_partner_is_not_subtyped():
    # legacy flag preserved (still non_RTE_SV) but no specific rearrangement subtype claimed.
    c = LOCAL.conclusion()
    assert VAC.element_class(c) == "non_RTE_SV", c
    assert "inversion" not in c and "translocation junction" not in c, c


def test_lowmapq_repeat_is_not_translocation():
    c = LOWMAPQ.conclusion()
    assert "translocation" not in c and "inversion" not in c, c
    assert VAC.element_class(c) == "unknown", c


def test_sv_and_mei_not_mutually_exclusive():
    c = ALU_PLUS_INVERSION.conclusion()
    assert "[SV: inversion" in c, c            # additive note present
    assert VAC.element_class(c) == "ALU", c    # class stays the element's, not SV


def test_transduction_guard_suppresses_cooriented_note():
    c = ALU_TRANSDUCTION_GUARD.conclusion()
    assert "[SV:" not in c, c                   # co-oriented partner + poly-A -> no note
    assert VAC.element_class(c) == "ALU", c


def test_repeat_paralog_partner_not_reported_on_element():
    # a lone repeat-paralog partner (no reciprocity) must NOT be annotated on the Alu call.
    c = ALU_REPEAT_PARTNER.conclusion()
    assert "[SV:" not in c, c
    assert VAC.element_class(c) == "ALU", c


# --- artefact vs uncharacterised complex insertion (all-evidence-empty split) -------------

# both clips carry a long, high-complexity inserted sequence but nothing identified them (no
# Dfam, no remap, no poly-A): a real but uncharacterised non-MEI/complex insertion, NOT noise.
_HICX_L = "gcgtctgcacaggggcccagctgaggcagcccggc"     # 35 bp, high-entropy lower-case insert...
_HICX_R = "cggcaggcagcctcctcgtcagctgcttg"           # 29 bp insert (real chr8:142834932 clips)
COMPLEX_INS = _make("chr8:142834932-142834940",
                    left_seq=_HICX_L + "ACGTACGTACGTACGT",       # lower insert then ref (upper)
                    right_seq="ACGTACGTACGTACGT" + _HICX_R)       # ref (upper) then lower insert

# short / low-complexity clips with nothing identified: a genuine artefact (unchanged).
SHORT_ARTEFACT = _make("chr9:5000-5010",
                       left_seq="acg" + "ACGTACGTACGT",           # 3 bp insert -> too short
                       right_seq="ACGTACGTACGT" + "aca")          # 3 bp insert
LOWCOMPLEX_ARTEFACT = _make("chr9:6000-6010",
                            left_seq="a" * 30 + "ACGTACGTACGT",   # 30 bp but homopolymer
                            right_seq="ACGTACGTACGT" + "t" * 30)


def test_complex_insertion_not_artefact():
    c = COMPLEX_INS.conclusion()
    assert "unmapped complex insertion" in c, c
    assert VAC.element_class(c) == "unknown", c


def test_short_clips_stay_artefact():
    c = SHORT_ARTEFACT.conclusion()
    assert c == "artefact", c
    assert VAC.element_class(c) == "artefact", c


def test_lowcomplexity_clips_stay_artefact():
    # 30 bp but poly-A/poly-T homopolymer: low entropy -> stays artefact (this is the slippage
    # mechanism the artefact bin is meant to hold).
    c = LOWCOMPLEX_ARTEFACT.conclusion()
    assert c == "artefact", c


# --- reference-free local-duplication / microsatellite subtypes (workers on TP "unknown") ------

def _rc(s):
    return s.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]

# Clean synthetic tandem duplication: the duplicated unit is (Uend|core|Dstart). The LEFT clip
# reads Uend+core (its first 20 bp = the END of the RIGHT junction's flank); the RIGHT clip reads
# core+Dstart (its last 20 bp = the START of the LEFT junction's flank). Reciprocal flank copy.
_UEND = "CATGACTGCATTGGCATCGA"
_CORE = "GGTTCAAGTCCATCGATCCA"
_DSTART = "TTGCACCGATAGCATGGACA"
_MOREU = "ACGTGACCTTAGCTAGCTGA"
_MORED = "GCTTAGGCATCGATCGATGC"
TANDEM_DUP = _make("chr7:5000-5003",
                   left_seq=(_UEND + _CORE).lower() + (_DSTART + _MORED),
                   right_seq=(_MOREU + _UEND) + (_CORE + _DSTART).lower())

# Inverted duplication: each clip is the REVERSE-COMPLEMENT copy of the opposite flank.
INVERTED_DUP = _make("chr7:8000-8003",
                     left_seq=(_UEND + _CORE).lower() + (_DSTART + _MORED),
                     right_seq=(_MOREU + _rc(_UEND)) + (_CORE + _rc(_DSTART)).lower())

# Real microsatellite (CA)n length change from chr3:190742584 (worker 9): both clips are the
# ordinary flanks on the far side of a (CA)n tract flush to both breakpoints.
MICROSAT = _make(
    "chr3:190742584-190742595",
    left_seq=("gattattttaatgtgtagcctaaagtaaa"
              "ACACACACACACACACACACACACAAATAGGCACAATTAGTGAGACTCAATAGAGAAGAGAAACCCATAC"
              "AAGTCAGACCTGCAATTATATCAGGTACTGAAATACAAAATATATGACAGAATATTTAGAATATATGTTAAT"),
    right_seq=("GTCACCTACCAAAGCAAACACTAATTCTCTGAAACAAATACAATTTCAACCTAGACCTTAAAACTTCCCACAA"
               "TTTCTTAGAATATGAATATTCTTAGAATATGAATCTGTAGCCTAAAGTAAAACACACACACACACACACACACA"
               "aataggcacatagaga"))


def test_tandem_duplication():
    c = TANDEM_DUP.conclusion()
    assert "tandem/segmental duplication" in c, c
    assert VAC.element_class(c) == "SV_DUPLICATION", c


# --- plain local SVs: both clips are the reference across the partner breakpoint (PD51635) ------
# Real junction consensus strings from the PD51635 run (2026-10-07).

# chr2:99670833-99670814: a 19 bp deletion (hg38 99670815-99670833 = TTCCTTTTATTGCAAAAAT), typed
# TSD_DELETION / "unknown" before; microhomology AG.
SV_DEL = _make(
    "chr2:99670833-99670814",
    left_seq="ataaaaattattcattctgtcAGGAAGGTATATTCACATAGGAATTAAATGCAAGAATAAAGGATCTGACCTGTACATTCAT",
    right_seq="GGTTATAATTATAGGCCAAAGCAGGACAACTTAATGTATAAAAATTATTCATTCTGTCAGgaaggtatattcaca")

# chrX:141289677-141289690: a 13 bp tandem duplication with a 2 bp untemplated seam (TA), which
# the TPRT score read as a TSD (LIKELY_TPRT).
SV_DUP_SEAM = _make(
    "chrX:141289677-141289690",
    left_seq="ttcctagtcctcagtaTCCTAGTCCTCAGAACTACTGAAGACTGATACAAACTGAATCTTTCACC",
    right_seq="CCCTAGAACCTACTTTGTCTGTTCCTAGTCCTCAGtatcctagtcctcagaactact")

# A TSD insertion: the clips are an Alu 3' end / poly-A, never the partner flank -> no SV label.
TSD_ALU = _make(
    "chr5:1000-1015",
    left_seq="ggccgggcgcggtggctcacgcctgtaatcccagcaCTTAGGCATCGATCGATGCTTGCACCGATAGC",
    right_seq="ACGTGACCTTAGCTAGCTGACATGACTGCATTGGCATCGAaaaaaaaaaaaaaaaaaaaaaaa")

# poly-A clips that happen to sit against an A-run in the partner flank: low complexity -> no SV.
POLYA_CHANCE = _make(
    "chr5:2000-1990",
    left_seq="aaaaaaaaaaaaaaaaCTTAGGCATCGATCGATGCTTGCACCG",
    right_seq="ACGTGACCTTAGCTAGCTGAAAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaa")


# chr6:70651185-70651180: the 5 bp "seam" AATAG replaces hg38 70651181-70651185 = CTATT, its
# reverse complement -> a micro-inversion inside a palindrome, not a deletion.
MICRO_INV = _make(
    "chr6:70651185-70651180",
    left_seq="cacccttatcaatagGATAAGGGTGGAGCCCTCGTGACTTAATTACTGTCCCAAGG",
    right_seq="GGCCAATTTTATTACTACCCTTATCaataggataagggtag")


class _StubGenome:
    """hg38 chr6:70651120-70651250 (0-based start), enough for the micro-inversion fixture."""
    OFF, SEQ = 70651120, ("TCCTCATATGGTCGAAGAGGCTCTCTGCTCCCTTAGGCCAATTTTATTACTACCCTTATCCTATTGATAAGGGTGG"
                          "AGCCCTCGTGACTTAATTACTGTCCCAAGGGCCCAACTCTATATCACTTCAATG")

    def fetch(self, contig, start, end):
        assert contig == "chr6"
        return self.SEQ[start - self.OFF:end - self.OFF]


def test_micro_inversion_needs_the_reference():
    saved = Insertion._GENOME
    try:
        Insertion._GENOME = None                   # no genome: seam == deleted length -> no call
        assert MICRO_INV._breakpoint_sv_subtype() is None
        Insertion._GENOME = _StubGenome()
        c = MICRO_INV.conclusion()
        assert c.startswith("SV_INVERSION: 5 bp micro-inversion (CTATT -> AATAG"), c
        assert VAC.element_class(c) == "SV_INVERSION", c
    finally:
        Insertion._GENOME = saved


def test_plain_deletion_is_sv_deletion():
    c = SV_DEL.conclusion()
    assert c.startswith("SV_DELETION: 19 bp deletion"), c
    assert "microhomology 2 bp" in c, c
    assert VAC.element_class(c) == "SV_DELETION", c


def test_tandem_dup_with_seam_is_sv_duplication():
    c = SV_DUP_SEAM.conclusion()
    assert c.startswith("SV_DUPLICATION: 13 bp tandem duplication"), c
    assert "seam TA" in c, c
    assert VAC.element_class(c) == "SV_DUPLICATION", c


def test_tsd_insertion_is_not_sv():
    assert TSD_ALU._breakpoint_sv_subtype() is None


def test_polya_chance_match_is_not_sv():
    assert POLYA_CHANCE._breakpoint_sv_subtype() is None


def test_inverted_duplication():
    c = INVERTED_DUP.conclusion()
    assert "inverted duplication" in c, c
    assert VAC.element_class(c) == "non_RTE_SV", c


def test_microsatellite_length_change():
    c = MICROSAT.conclusion()
    assert "microsatellite" in c, c
    assert VAC.element_class(c) == "microsatellite", c


def test_flank_leak_dfam_demoted():
    # A short novel insert whose soft-clip drags in genomic flank that carries an (old-repeat)
    # Dfam hit: the hit sits on the flank-copied portion of the clip, so demote_flank_leak_dfams
    # removes it and no false element is called. The clip region under the hit is an exact
    # copy (>= flankleak_min = 25 bp) of this insertion's own left flank.
    FLANK30 = "TTGCACCGATAGCATGGACAGCTTAGGCAT"        # 30 bp left flank
    left_seq = "acgtacgt" + FLANK30                     # left uppercase flank == FLANK30
    novel = "GGTTCAAGTCCATCGATCCA"                      # 20 bp novel
    right_ins = novel + FLANK30                          # 50 bp clip = novel + flank copy
    right_seq = "AAAACCCCGGGGTTTTACGTAC" + right_ins.lower()
    ins = _make("chr9:100-103", left_seq=left_seq, right_seq=right_seq,
                right_dfams=[_dfam("L1MC4_3end", 40.0, "+")])
    ins.right_dfams[0].ali_start = 21   # hit covers clip[20:50] = the 30 bp flank copy
    ins.right_dfams[0].ali_end = 50
    cont = VAC.__new__(VAC)
    cont.insertions = {"x": ins}
    cont.demote_flank_leak_dfams()
    assert ins.right_dfams == [], "flank-leak L1 hit should be demoted"


def test_sv_clip_trust_guard_blocks_short_lowcomplex_map():
    # A short, AT-rich clip that 'uniquely' maps to another chromosome must NOT be reported as a
    # translocation partner (the map is a chance/paralogous placement).
    ins = _make("chr6:48672036-48672049",
                left_seq="attat" + "GCATCGATCGATCGATCGAT",     # 5 bp AT-rich clip -> untrusted
                left_maps=[_map("chr12:33095027+", 60, "+")])
    c = ins.conclusion()
    assert "translocation" not in c, c


def test_reciprocal_balanced_translocation():
    # Build two junctions that point back at each other; the container pairing pass must mark
    # both balanced. Bypass __init__ (which reads files) with __new__.
    a = _make("chrA:1000-1010", left_maps=[_map("chrB:5000+", 60, "+")])
    b = _make("chrB:4990-5010", left_maps=[_map("chrA:1005+", 60, "+")])
    cont = VAC.__new__(VAC)
    cont.insertions = {"a": a, "b": b}
    cont.link_reciprocal_translocations()
    assert a.reciprocal_partner is not None and b.reciprocal_partner is not None
    ca = a.conclusion()
    assert "balanced translocation" in ca, ca
    assert VAC.element_class(ca) == "non_RTE_SV", ca


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
