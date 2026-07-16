# PEAR-TREE - paired ends of aberrant retrotransposons in phylogenetic trees
#
# Copyright (C) 2025 Jeremy Deuel <jeremy.deuel@usz.ch>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Drop insertions that fall in regions carrying an implausibly high insertion rate.

Kept as a pure function (no pysam/py2bit/CONFIG imports) so it can be exercised
directly against the reference implementation it replaced -- see
test/test_region_filter.py.
"""

from collections import Counter


def _binned_pos(i):
    """The insertion's position, and the bin it lands in."""
    pos = i.right_pos if i.right_pos is not None else i.left_pos
    return pos


def filter_dense_regions(insertions, bin_range=100, ins_cutoff=4):
    """Remove every insertion covered by a bin holding >= ins_cutoff insertions.

    A bin at offset `b` covers `b - bin_range/2 < pos < b + bin_range*1.5`, i.e. a
    window twice the bin width, so neighbouring bins' windows overlap and an insertion
    can be covered by a bin other than its own.

    Returns (kept, n_removed, n_regions).
    """
    bins = [(i.reference_name, int(_binned_pos(i) / bin_range) * bin_range) for i in insertions]
    c = Counter(bins)
    hot_bins = {name for name, n in c.items() if n >= ins_cutoff}

    # The hot set is fixed before any filtering, so an insertion is dropped iff SOME hot
    # bin covers it -- no need to rescan every insertion once per hot bin. Only bins
    # strictly inside (pos - bin_range*1.5, pos + bin_range/2) can satisfy the coverage
    # predicate, so probing the handful of candidates around pos is exhaustive; the
    # predicate applied to each candidate is unchanged from the per-region rescan.
    kept = []
    removed = 0
    for i in insertions:
        pos = _binned_pos(i)
        base = ((pos - int(bin_range * 1.5)) // bin_range) * bin_range
        for bin in range(base, pos + int(bin_range / 2) + bin_range, bin_range):
            if (i.reference_name, bin) in hot_bins and pos > bin - bin_range / 2 and pos < bin + bin_range * 1.5:
                removed += 1
                break
        else:
            kept.append(i)
    return kept, removed, len(hot_bins)
