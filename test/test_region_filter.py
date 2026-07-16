# PEAR-TREE - dense-region filter equivalence tests
#
# The dense-region filter in combine_insertions used to rescan the whole insertion
# list once per hot bin (O(hot_regions * insertions)). It now makes a single pass
# against a precomputed set of hot bins. That rewrite is only worth anything if it
# is exactly output-neutral, so these tests pin it against the original loop,
# copied verbatim below as an oracle.
#
# Run with:  python test/test_region_filter.py     (from the repo root)
#        or:  pytest test/test_region_filter.py

import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from combine_insertions_region_filter import filter_dense_regions


class FakeInsertion:
    """Minimal stand-in: the filter only reads reference_name/left_pos/right_pos."""

    def __init__(self, reference_name, left_pos, right_pos, uid):
        self.reference_name = reference_name
        self.left_pos = left_pos
        self.right_pos = right_pos
        self.uid = uid  # identity tag, so we can compare surviving sets exactly

    def __repr__(self):
        return f"Ins({self.reference_name}, {self.left_pos}, {self.right_pos}, #{self.uid})"


def reference_filter(insertions, bin_range=100, ins_cutoff=4):
    """The ORIGINAL implementation, copied verbatim from combine_insertions.py
    (commit 249eee2, lines 94-113) with only the print removed. Do not 'clean up':
    its job is to be the thing we did not write."""
    bins = []
    for i in insertions:
        bins.append((i.reference_name, int(i.right_pos / bin_range) * bin_range if i.right_pos is not None else int(
            i.left_pos / bin_range) * bin_range))
    c = Counter(bins)
    removed = 0
    regions = 0
    for (reference_name, bin), n in c.items():
        clean_ins = []
        if n >= ins_cutoff:
            regions += 1
            for i in insertions:
                pos = i.right_pos if i.right_pos is not None else i.left_pos
                if i.reference_name == reference_name and pos > bin - bin_range / 2 and pos < bin + bin_range * 1.5:
                    removed += 1
                else:
                    clean_ins.append(i)
            insertions = clean_ins
    return insertions, removed, regions


def assert_agree(insertions, bin_range=100, ins_cutoff=4, label=""):
    exp_kept, exp_removed, exp_regions = reference_filter(list(insertions), bin_range, ins_cutoff)
    got_kept, got_removed, got_regions = filter_dense_regions(list(insertions), bin_range, ins_cutoff)
    assert [i.uid for i in got_kept] == [i.uid for i in exp_kept], (
        f"{label}: surviving set/order differs\n  expected {exp_kept}\n  got      {got_kept}")
    assert got_removed == exp_removed, f"{label}: removed count {got_removed} != {exp_removed}"
    assert got_regions == exp_regions, f"{label}: region count {got_regions} != {exp_regions}"


def test_randomised_agreement():
    """Dense random positions over few contigs: guarantees plenty of hot bins, and
    positions clustered tightly enough to exercise the overlapping-window edges."""
    rng = random.Random(20260716)
    for trial in range(300):
        n = rng.randint(0, 60)
        insertions = []
        for uid in range(n):
            # Narrow span => many collisions per bin => many hot bins.
            pos = rng.randint(0, 600)
            contig = rng.choice(["chr1", "chr2", "1"])
            if rng.random() < 0.3:
                insertions.append(FakeInsertion(contig, pos, None, uid))  # right_pos None => left_pos used
            else:
                insertions.append(FakeInsertion(contig, rng.randint(0, 600), pos, uid))
        assert_agree(insertions, label=f"random trial {trial}")


def test_boundary_positions():
    """Every position in a window around the bin edges, at the cutoff exactly.

    The coverage predicate is strict (`>` / `<`) on bin-50 and bin+150, and bins are
    multiples of 100 -- so positions landing exactly on a boundary are where an
    off-by-one in the candidate enumeration would surface.
    """
    for pos in range(0, 400):
        # ins_cutoff=1 makes every occupied bin hot, maximising the number of hot bins
        # whose windows can reach `pos` and so the chance of a missed candidate.
        insertions = [FakeInsertion("chr1", pos, pos, 0)]
        assert_agree(insertions, ins_cutoff=1, label=f"boundary pos={pos}")
        for other in range(0, 400, 25):
            pair = [FakeInsertion("chr1", pos, pos, 0), FakeInsertion("chr1", other, other, 1)]
            assert_agree(pair, ins_cutoff=1, label=f"boundary pos={pos} other={other}")


def test_low_positions_near_zero():
    """Positions below bin_range*1.5 make the candidate base offset go negative."""
    for pos in range(0, 200):
        insertions = [FakeInsertion("chr1", pos, pos, u) for u in range(4)]
        assert_agree(insertions, ins_cutoff=1, label=f"near-zero pos={pos}")


def test_contig_isolation():
    """A hot bin on one contig must not reach an insertion on another."""
    insertions = [FakeInsertion("chr1", 100, 100, u) for u in range(5)]
    insertions += [FakeInsertion("chr2", 100, 100, 99)]
    kept, removed, regions = filter_dense_regions(insertions)
    assert [i.uid for i in kept] == [99], kept
    assert (removed, regions) == (5, 1)
    assert_agree(insertions, label="contig isolation")


def test_below_cutoff_survives():
    """Three insertions in a bin is under the cutoff of four -- nothing is dropped."""
    insertions = [FakeInsertion("chr1", 100, 100, u) for u in range(3)]
    kept, removed, regions = filter_dense_regions(insertions)
    assert len(kept) == 3 and removed == 0 and regions == 0
    assert_agree(insertions, label="below cutoff")


def test_neighbouring_bin_reach():
    """A bin's window is twice its width, so a hot bin drops insertions in the
    adjacent bin too -- the behaviour a naive same-bin-only rewrite would lose."""
    insertions = [FakeInsertion("chr1", 100, 100, u) for u in range(4)]     # bin 100, hot
    insertions.append(FakeInsertion("chr1", 200, 200, 50))                  # bin 200, cold, but 100<200<250
    kept, removed, regions = filter_dense_regions(insertions)
    assert [i.uid for i in kept] == [], f"hot bin at 100 must also reach pos 200: {kept}"
    assert_agree(insertions, label="neighbour reach")


def test_candidate_enumeration_is_exhaustive():
    """The rewrite's load-bearing claim, checked directly rather than by sampling.

    The old code tested every hot bin against every insertion. The new code only probes
    bins in `range(base, pos + bin_range/2 + bin_range, bin_range)`. That is only safe if
    NO bin outside that range can ever satisfy the coverage predicate. Here we brute-force
    every multiple of bin_range within a wide margin either side of pos and assert the
    candidate list is a superset of the genuinely-covering ones -- so the enumeration can
    never skip a hot bin the old rescan would have found.
    """
    bin_range = 100
    for pos in range(0, 3000):
        base = ((pos - int(bin_range * 1.5)) // bin_range) * bin_range
        candidates = set(range(base, pos + int(bin_range / 2) + bin_range, bin_range))
        covering = {
            b for b in range(-1000, 5000, bin_range)
            if pos > b - bin_range / 2 and pos < b + bin_range * 1.5
        }
        assert covering <= candidates, (
            f"pos={pos}: bins {sorted(covering - candidates)} cover it but are never probed")


def test_empty():
    assert filter_dense_regions([]) == ([], 0, 0)


def _run_all():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"\n{len(tests)} tests passed")


if __name__ == "__main__":
    _run_all()
