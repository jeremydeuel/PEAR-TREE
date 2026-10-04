"""Indel-aware clip consensus (src/indel_consensus.py).

The motivating failure: Illumina SBS poly-A length jitter. Clips =
20 bp element + poly-A of 18 +- 3 (per read) + 30 bp beyond-poly-A sequence.
The legacy column-wise majority (consensus.find_consensus) left-aligns the clips and
stops at the first ambiguous column -- inside/just after the poly-A -- so the beyond
sequence is lost. The RLE star alignment must recover it exactly.

Run:  pytest test/test_indel_consensus.py
"""
import random

import pytest

from _combine_shim import consensus  # legacy find_consensus (needs the config shim)
from indel_consensus import ClipRead, indel_aware_consensus, revcomp, rle
from quality_seq import QualitySeq

ELEM = "GGCTCACGCCTGTAATCCCG"            # 20 bp, does not end in A (keeps poly-A length exact)
BEYOND = "GTCAGGATCCTGACCGTTAGCCAGTCTTGC"  # 30 bp
POLYA = 18


def _mutate(rng, s, err=0.01):
    s = list(s)
    for k in range(len(s)):
        if rng.random() < err:
            s[k] = rng.choice([b for b in "ACGT" if b != s[k]])
    return s


def synth_reads(seed, n=14, jitter=3, err=0.01, n_hp_indel=3):
    """n clips; per-read poly-A jitter; ~1% substitutions (never on the last 2 bases of
    the two full-length reads); a few reads with a homopolymer indel in the element
    (CCC -> CCCC / CC, GG -> G)."""
    rng = random.Random(seed)
    reads = []
    pa_true = []
    for i in range(n):
        pa = POLYA + rng.randint(-jitter, jitter)
        elem = ELEM
        if i < n_hp_indel:
            elem = [ELEM.replace("CCCG", "CCCCG"), ELEM.replace("CCCG", "CCG"),
                    ELEM.replace("GG", "G", 1)][i % 3]
        full = elem + "A" * pa + BEYOND
        s = _mutate(rng, full[:-2], err) + list(full[-2:]) if i < 2 else _mutate(rng, full, err)
        pa_true.append(pa)
        if i >= 2:  # the first two reads are full length; the rest end randomly
            s = s[:rng.randint(len(s) - 12, len(s))]
        q = [rng.randint(25, 38) for _ in s]
        reads.append(ClipRead("".join(s), q, group=i))
    synth_reads.true_polya = pa_true
    return reads


@pytest.mark.parametrize("seed", range(25))
def test_recovers_beyond_polya_exactly(seed):
    reads = synth_reads(seed)
    true = sorted(synth_reads.true_polya)
    true_median = (true[len(true) // 2 - 1] + true[len(true) // 2]) / 2
    res = indel_aware_consensus(reads)
    assert res.seq.startswith(ELEM)
    assert res.polya_base == "A"
    # the median of the reads' true lengths (sampled from 18 +- 3), within 1 bp
    assert abs(res.polya_len_median - true_median) <= 1
    assert abs(res.polya_len_median - POLYA) <= 2
    assert res.polya_len_min >= POLYA - 3 and res.polya_len_max <= POLYA + 3
    assert res.beyond_polya == BEYOND
    assert res.beyond_polya_support >= 2
    assert len(res.depth) == len(res.seq)


def test_legacy_column_majority_loses_beyond_polya():
    """Show the old method fails on the same data (and the new one does not)."""
    n_legacy_ok = n_new_ok = 0
    for seed in range(25):
        reads = synth_reads(seed)
        legacy = consensus.find_consensus([QualitySeq(r.seq, list(r.qual)) for r in reads])
        new = indel_aware_consensus(reads)
        n_legacy_ok += BEYOND in str(legacy)
        n_new_ok += new.beyond_polya == BEYOND
        # legacy stops inside the poly-A / right after it
        assert len(legacy) < len(ELEM) + POLYA + 3 + 5
    assert n_legacy_ok == 0
    assert n_new_ok == 25


def test_stops_when_beyond_has_single_fragment():
    rng = random.Random(7)
    reads = []
    for i in range(6):
        s = ELEM + "A" * (POLYA + rng.randint(-2, 2)) + (BEYOND if i == 0 else BEYOND[:3])
        reads.append(ClipRead(s, [30] * len(s), group=i))
    res = indel_aware_consensus(reads)
    assert res.stop_reason == "depth"
    assert res.beyond_polya == BEYOND[:3]


def test_same_fragment_counts_once_for_depth():
    """Two reads of ONE fragment (read + overlapping mate, or a PCR family) carry the
    beyond sequence: depth must be 1 -> not reported."""
    reads = [ClipRead(ELEM + "A" * 18 + BEYOND, [30] * 68, group="f1"),
             ClipRead(ELEM + "A" * 19 + BEYOND, [30] * 69, group="f1"),
             ClipRead(ELEM + "A" * 17, [30] * 37, group="f2"),
             ClipRead(ELEM + "A" * 16, [30] * 36, group="f3")]
    res = indel_aware_consensus(reads)
    assert res.beyond_polya == ""
    assert res.seq.startswith(ELEM + "A")


def test_quality_weighting_overrides_low_quality_majority():
    hi = ELEM + "T" + "GATTACA"
    lo = ELEM + "C" + "GATTACA"
    reads = [ClipRead(hi, [40] * len(hi), 1), ClipRead(hi, [40] * len(hi), 2),
             ClipRead(lo, [30] * 20 + [3] + [30] * 7, 3), ClipRead(lo, [30] * 20 + [3] + [30] * 7, 4),
             ClipRead(lo, [30] * 20 + [3] + [30] * 7, 5)]
    res = indel_aware_consensus(reads)
    assert res.seq == hi


def test_genuine_disagreement_stops():
    a = ELEM + "TTTTGACCA"
    b = ELEM + "CCGAGGTCA"
    reads = [ClipRead(a, [30] * len(a), 1), ClipRead(a, [30] * len(a), 2),
             ClipRead(b, [30] * len(b), 3), ClipRead(b, [30] * len(b), 4)]
    res = indel_aware_consensus(reads)
    assert res.seq == ELEM
    assert res.stop_reason == "disagreement"


def test_mates_extend_consensus_where_they_overlap():
    """Junction reads cover ELEM+polyA+10 of BEYOND; two mates of different fragments
    (given in the opposite orientation, as an unmapped mate may be) overlap the clip
    and carry the rest of BEYOND."""
    clip = ELEM + "A" * 18 + BEYOND[:10]
    reads = [ClipRead(clip, [30] * len(clip), g) for g in ("a", "b", "c")]
    mate = BEYOND[:]
    reads += [ClipRead(revcomp(ELEM[-5:] + "A" * 18 + mate), [30] * (23 + len(mate)), g, anchored=False)
              for g in ("a", "b")]
    # an unrelated mate must not be placed
    junk = "ACGTTGCAAGGCTTAACCGGTA"
    reads.append(ClipRead(junk, [30] * len(junk), "c", anchored=False))
    res = indel_aware_consensus(reads)
    assert res.beyond_polya == BEYOND
    assert res.beyond_polya_support == 2


def test_polyT_at_junction_has_empty_beyond():
    """3' junction seen from the flank: outward clip starts with the poly-T."""
    rng = random.Random(3)
    reads = []
    for i in range(5):
        s = "T" * (15 + rng.randint(-2, 2)) + revcomp(ELEM)
        reads.append(ClipRead(s, [30] * len(s), i))
    res = indel_aware_consensus(reads)
    assert res.polya_base == "T" and res.polya_start == 0
    assert abs(res.polya_len_median - 15) <= 1
    assert res.beyond_polya == ""
    assert res.seq.endswith(revcomp(ELEM))


def test_rle():
    assert rle("AAACG", [10, 20, 30, 40, 40]) == ("ACG", [3, 1, 1], [20.0, 40.0, 40.0])


def test_empty_and_single():
    assert indel_aware_consensus([]).seq == ""
    r = indel_aware_consensus([ClipRead("ACGTACGTAC", [30] * 10, 1)])
    assert r.seq == "" and r.stop_reason == "depth"
