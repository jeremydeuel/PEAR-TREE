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
from typing import Tuple
from genotype_qscore import qscore, ARTEFACT, ALT_MATCH, REF_MATCH, RIGHT_TO_LEFT, LEFT_TO_RIGHT
from quality_seq import QualitySeq
# define constants used to identify matching conditions
MATCHES_WT = 2
MATCHES_CLIP = 1
MATCHES_NEITHER = -1
NO_COVERAGE_OF_BREAKPOINT = 0
HIGH_COVERAGE = -2

DEBUG = True

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
        if breakpoint >= self.read.reference_start and breakpoint <= self.read.reference_end:
            breakpoint_query = None
            for query, ref_pos in self.read.get_aligned_pairs(True):
                if ref_pos == breakpoint:
                    breakpoint_query = query
                    break
            if breakpoint_query is None:
                self.left_genotype = 0,0,0
            self.left_genotype = qscore(QualitySeq(self.read.query_sequence[:breakpoint_query],self.read.query_qualities[:breakpoint_query]),
                          ref, alt, RIGHT_TO_LEFT)

        else:
            self.left_genotype = 0,0,0

    def qright(self, breakpoint: int, ref: str, alt: str):
        if breakpoint >= self.read.reference_start and breakpoint <= self.read.reference_end:
            breakpoint_query = None
            for query, ref_pos in reversed(self.read.get_aligned_pairs(True)):
                if ref_pos == breakpoint-1:
                    breakpoint_query = query+1
                    break
            if breakpoint_query is None:
                self.right_genotype = 0,0,0
            self.right_genotype = qscore(QualitySeq(self.read.query_sequence[breakpoint_query:],self.read.query_qualities[breakpoint_query:]),
                          ref, alt, LEFT_TO_RIGHT)
        else:
            self.right_genotype = 0,0,0

    def __str__(self):
        return f"read {self.read.query_name}: right_gt={self.right_genotype}, left_gt={self.left_genotype}"