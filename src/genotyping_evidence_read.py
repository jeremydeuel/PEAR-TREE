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


import pysam
import os
from typing import Tuple
from genotype_qscore import qscore, RIGHT_TO_LEFT, LEFT_TO_RIGHT
from quality_seq import QualitySeq
# define constants used to identify matching conditions
MATCHES_WT = 2
MATCHES_CLIP = 1
MATCHES_NEITHER = -1
NO_COVERAGE_OF_BREAKPOINT = 0
HIGH_COVERAGE = -2

DEBUG = os.environ.get("PEARTREE_DEBUG", "") not in ("", "0", "false", "False")

# CIGAR op codes (BAM spec) grouped by what they consume. A "match" op aligns a query
# base to a reference base (M/=/X); those are exactly the pairs pysam's
# get_aligned_pairs(matches_only=True) returns.
_CONSUMES_QUERY = frozenset((0, 1, 4, 7, 8))   # M I S = X
_CONSUMES_REF = frozenset((0, 2, 3, 7, 8))     # M D N = X
_MATCH_OPS = frozenset((0, 7, 8))              # M = X


def query_index_at_ref(cigar, reference_start, target_ref):
    """Query offset aligned to reference position ``target_ref``.

    Returns the query index q such that (q, target_ref) is an aligned match pair, or
    None if target_ref is not aligned to a query base (it falls in a deletion / skip,
    or outside the read). This is exactly what searching
    ``get_aligned_pairs(matches_only=True)`` for ``ref == target_ref`` yields, but in
    O(#cigar ops) with no per-read list materialisation -- get_aligned_pairs was the
    genotyping hot spot (built a ~read-length tuple list twice per read).
    """
    if cigar is None:
        return None
    qpos = 0
    rpos = reference_start
    for op, length in cigar:
        if op in _MATCH_OPS and rpos <= target_ref < rpos + length:
            return qpos + (target_ref - rpos)
        if op in _CONSUMES_QUERY:
            qpos += length
        if op in _CONSUMES_REF:
            rpos += length
    return None


class EvidenceRead:
    def __init__(self, read: pysam.AlignedRead):
        self.read = read
        self.cigar = self.read.cigartuples
        self.left_genotype: Tuple[float, float, float] = 0,0,0 #ref, alt, art
        self.right_genotype: Tuple[float, float, float] = 0,0,0


    def has_original_sequence_direction(self):
        """
        Returns True if the read has the original sequencing direction, else False
        Original sequencing direction:

            R1 =======>     <========= R2
        """
        return self.read.is_read1 == self.read.is_forward

    def qleft(self, breakpoint: int, ref: str, alt: str):
        read = self.read
        if breakpoint < read.reference_start or breakpoint > read.reference_end:
            self.left_genotype = 0, 0, 0
            return
        # query base aligned at the breakpoint; score the read portion 5' of it (RIGHT_TO_LEFT)
        breakpoint_query = query_index_at_ref(self.cigar, read.reference_start, breakpoint)
        if breakpoint_query is None:
            self.left_genotype = 0, 0, 0
            return
        self.left_genotype = qscore(
            QualitySeq(read.query_sequence[:breakpoint_query], read.query_qualities[:breakpoint_query]),
            ref, alt, RIGHT_TO_LEFT)

    def qright(self, breakpoint: int, ref: str, alt: str):
        read = self.read
        if breakpoint < read.reference_start or breakpoint > read.reference_end:
            self.right_genotype = 0, 0, 0
            return
        # query base aligned just 5' of the breakpoint (ref == breakpoint-1); score the read
        # portion 3' of it (query index +1, LEFT_TO_RIGHT), matching the original semantics.
        q = query_index_at_ref(self.cigar, read.reference_start, breakpoint - 1)
        if q is None:
            self.right_genotype = 0, 0, 0
            return
        breakpoint_query = q + 1
        self.right_genotype = qscore(
            QualitySeq(read.query_sequence[breakpoint_query:], read.query_qualities[breakpoint_query:]),
            ref, alt, LEFT_TO_RIGHT)

    def __str__(self):
        return f"read {self.read.query_name}: right_gt={self.right_genotype}, left_gt={self.left_genotype}"