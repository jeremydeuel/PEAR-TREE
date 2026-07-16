# PEAR-TREE - genotyping unit tests
#
# Deterministic tests for the allele-fraction (VAF) genotype model and the
# ±1 bp qscore register selection. These lock in the behaviour that motivated
# the genotyping code review (see plans/improve_discovery): a stray read must
# not flip a homozygote to het, low-level contamination must read as wild-type,
# both-ends-alt reads are chimeric artefacts, and low-evidence calls degrade to
# the explicit uncertain classes rather than fabricating confident calls.
#
# Run with:  python test/test_genotyping.py     (from the repo root)
#        or:  pytest test/test_genotyping.py

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from genotyping_insertion import (
    Insertion,
    GT_WILDTYPE, GT_WILDTYPE_UNCERTAIN,
    GT_HETEROZYGOUS, GT_HOMOZYGOUS,
    GT_INSERTION_UNCERTAIN, GT_ARTEFACT,
    GT_NO_COVERAGE,
)
from genotype_qscore import qscore, LEFT_TO_RIGHT, RIGHT_TO_LEFT
from quality_seq import QualitySeq


class FakeRead:
    """Minimal stand-in for EvidenceRead: only left/right genotype tuples matter."""
    def __init__(self, left, right):
        self.left_genotype = left
        self.right_genotype = right


# Canonical per-side score tuples (ref, alt, art):
REF_SIDE = (300.0, -100.0, 0.0)   # matches genomic reference -> 'ref'
ALT_SIDE = (-100.0, 300.0, 0.0)   # matches inserted element  -> 'alt'
ART_SIDE = (-50.0, -50.0, 90.0)   # matches neither           -> 'art'
NONE_SIDE = (0.0, 0.0, 0.0)       # breakpoint not covered


def call(reads):
    ins = Insertion("1:100-108")
    ins.evidence_reads = reads
    return ins.summarise_evidence()


def wildtype_read():
    return FakeRead(REF_SIDE, REF_SIDE)


def insertion_read():
    # a real junction read: alt on the crossed junction, genomic ref on the other
    return FakeRead(ALT_SIDE, REF_SIDE)


CASES = []


def case(fn):
    CASES.append(fn)
    return fn


@case
def test_no_coverage():
    assert call([])[0] == GT_NO_COVERAGE


@case
def test_clean_heterozygous():
    gt, sg, sa = call([wildtype_read()] * 5 + [insertion_read()] * 5)
    assert gt == GT_HETEROZYGOUS, gt          # VAF 0.5


@case
def test_homozygous_not_flipped_by_one_stray_wildtype_read():
    # G1 scenario (a): a single mismapped wild-type read must NOT turn a
    # homozygote into a heterozygote (the old q_hom<0 logic did exactly that).
    gt, sg, sa = call([insertion_read()] * 10 + [wildtype_read()])
    assert gt == GT_HOMOZYGOUS, gt            # VAF 10/11 = 0.909


@case
def test_low_level_contamination_reads_as_wildtype():
    # G1 scenario (b): a handful of alt reads over many wild-type reads is
    # contamination / index hopping, not a clonal het.
    gt, sg, sa = call([wildtype_read()] * 40 + [insertion_read()] * 3)
    assert gt == GT_WILDTYPE, gt              # VAF 3/43 = 0.07 < 0.10


@case
def test_double_alt_is_artefact():
    # G12: one short read matching the inserted element on BOTH junctions is
    # geometrically impossible for a real long insertion -> artefact.
    gt, sg, sa = call([FakeRead(ALT_SIDE, ALT_SIDE)] * 5)
    assert gt == GT_ARTEFACT, gt


@case
def test_artefact_dominated_locus():
    gt, sg, sa = call([FakeRead(ART_SIDE, NONE_SIDE)] * 3 + [wildtype_read()])
    assert gt == GT_ARTEFACT, gt              # 3/4 art, >= fraction 0.5


@case
def test_uncertain_insertion_below_support():
    # G4/G13: VAF looks like an insertion (0.33) but only one supporting read,
    # below min_supporting_reads -> explicit 'insertion?' not a confident het.
    gt, sg, sa = call([wildtype_read()] * 2 + [insertion_read()])
    assert gt == GT_INSERTION_UNCERTAIN, gt


@case
def test_uncertain_wildtype_between_bands():
    # VAF in (vaf_wildtype_max, vaf_het_min) -> 'wild-type?'
    gt, sg, sa = call([wildtype_read()] * 5 + [insertion_read()])
    assert gt == GT_WILDTYPE_UNCERTAIN, gt    # VAF 1/6 = 0.167


@case
def test_qscore_prefers_matching_hypothesis():
    # sequence identical to ref: ref score positive, alt negative; and vice versa.
    seq = QualitySeq("ACGTACGT", [30] * 8)
    ref_pos, alt_pos, art = qscore(seq, "ACGTACGT", "TGCATGCA", LEFT_TO_RIGHT)
    assert ref_pos > 0 and alt_pos < 0, (ref_pos, alt_pos, art)
    ref_neg, alt_hi, art2 = qscore(seq, "TGCATGCA", "ACGTACGT", LEFT_TO_RIGHT)
    assert alt_hi > 0 and ref_neg < 0, (ref_neg, alt_hi, art2)


@case
def test_qscore_reports_artefact_for_neither():
    # a read matching neither ref nor alt must surface a positive artefact score
    seq = QualitySeq("GGGGGGGG", [30] * 8)
    ref, alt, art = qscore(seq, "ACGTACGT", "TTTTTTTT", LEFT_TO_RIGHT)
    assert art > 0 and ref < 0 and alt < 0, (ref, alt, art)


def main():
    failed = 0
    for fn in CASES:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except AssertionError as exc:
            failed += 1
            print(f"FAIL  {fn.__name__}: {exc}")
    print(f"\n{len(CASES) - failed}/{len(CASES)} passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
