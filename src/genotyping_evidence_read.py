# Geminate - Genomic events of mutation through insertional alterations by transposable elements
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
from config import CONFIG
from revcomp import revcomp

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
        self.left_match = None
        self.right_match = None


    def has_original_sequence_direction(self):
        """
        Returns True if the read has the original sequencing direction, else False
        Original sequencing direction:

            R1 =======>     <========= R2
        """
        return self.read.is_read1 == self.read.is_forward

    def classify_breakpoint(self, score_expected: int, score_clipped: int) -> int:
        if abs(score_clipped - score_expected) < 6 and max(score_expected, score_clipped) >= 2:
            return NO_COVERAGE_OF_BREAKPOINT
        if score_clipped > score_expected + 5 and score_clipped >= 3:
            return MATCHES_CLIP
        if score_expected > score_clipped + 5 and score_expected >= 3:
            return MATCHES_WT
        return None

    def check_left(self, breakpoint: int, clipped: str, expected: str) -> int:
        if breakpoint >= self.read.reference_start and breakpoint  < self.read.reference_end:
            breakpoint_query = None
            for query, ref in self.read.get_aligned_pairs(True):
                if ref==breakpoint:
                    breakpoint_query = query
                    break
            if not breakpoint_query:
                self.left_match = NO_COVERAGE_OF_BREAKPOINT
                return NO_COVERAGE_OF_BREAKPOINT

            test_seq = self.read.seq[:breakpoint_query]
            # remove adapter at read end
            test_seq_revcomp = revcomp(test_seq)
            for a in CONFIG['adapters']:
                if test_seq in a:
                    self.left_match = NO_COVERAGE_OF_BREAKPOINT
                    return NO_COVERAGE_OF_BREAKPOINT
                a_n = 0
                ts_n = 0
                score = 0
                while ts_n < min(len(test_seq_revcomp),16) and a_n < len(a):
                    if test_seq_revcomp[ts_n] == a[a_n]:
                        score += 1
                        a_n += 1
                        ts_n += 1
                        if score > 7:
                            break
                        continue
                    else:
                        if score > 4:
                            score -= 2
                            a_n += 1
                            ts_n += 1
                            continue
                        if a_n > 0:
                            a_n = 0
                            score = 0
                            continue
                        score = 0
                        a_n += 1
                        ts_n += 1
                if score > 3:
                    DEBUG and print(f"removed adapter LEFT: {test_seq} -> {test_seq[-(ts_n-a_n):]}")
                    test_seq = test_seq[-(ts_n-a_n):]
                    break
            #remove polyG at end (inverted since this is a left clipped read)
            for n in range(len(test_seq)):
                if test_seq[n] == 'C':
                    continue
                else:
                    if n>2:
                        #print(f"removed polyG LEFT: {test_seq} -> {test_seq[n:]}")
                        test_seq = test_seq[n:]
                    break
            min_len = min([len(x) for x in (test_seq, clipped, expected)])
            if min_len < 4:
                self.left_match = NO_COVERAGE_OF_BREAKPOINT
                return NO_COVERAGE_OF_BREAKPOINT
            score_clipped = 0
            score_expected = 0
            for i in range(1,min_len+1):
                if test_seq[-i] == clipped[-i]:
                    score_clipped += 1
                else:
                    score_clipped -= 2
                if test_seq[-i] == expected[-i]:
                    score_expected += 1
                else:
                    score_expected -= 2
            self.left_match = self.classify_breakpoint(score_expected, score_clipped)
            if self.left_match is not None:
                return self.left_match
            DEBUG and print(f"detected wrong match side=L score_wt={score_expected}, score_clip={score_clipped}")
            max_len = max([len(x) for x in (test_seq, clipped, expected)])
            DEBUG and print(f"measured seq: {' '*(max_len-len(test_seq)) + test_seq}")
            DEBUG and print(f"clipped seq:  {' '*(max_len-len(clipped)) + clipped}")
            DEBUG and print(f"wt seq:       {' '*(max_len-len(expected)) + expected}")
            self.left_match = MATCHES_NEITHER
            return MATCHES_NEITHER

        else:
            self.left_match = NO_COVERAGE_OF_BREAKPOINT
            return NO_COVERAGE_OF_BREAKPOINT

    def check_right(self, breakpoint: int, clipped: str, expected: str) -> int:
        if breakpoint > self.read.reference_start and breakpoint  <= self.read.reference_end:
            breakpoint_query = None
            for query, ref in reversed(self.read.get_aligned_pairs(True)):
                if ref == breakpoint-1:
                    breakpoint_query = query+1
                    break
            if not breakpoint_query:
                self.right_match = NO_COVERAGE_OF_BREAKPOINT
                return NO_COVERAGE_OF_BREAKPOINT
            test_seq = self.read.seq[breakpoint_query:]

            for a in CONFIG['adapters']:
                if test_seq in a:
                    self.right_match = NO_COVERAGE_OF_BREAKPOINT
                    return NO_COVERAGE_OF_BREAKPOINT
                a_n = 0
                ts_n = 0
                score = 0
                while ts_n < min(len(test_seq),16) and a_n < len(a):
                    if test_seq[ts_n] == a[a_n]:
                        score += 1
                        a_n += 1
                        ts_n += 1
                        if score > 7:
                            break
                        continue
                    else:
                        if score > 4:
                            score -= 2
                            a_n += 1
                            ts_n += 1
                            continue
                        if a_n > 0:
                            a_n = 0
                            score = 0
                            continue
                        score = 0
                        a_n += 1
                        ts_n += 1
                if score > 3:
                    DEBUG and print(f"clipped adapter RIGHT: {test_seq} -> {test_seq[:ts_n-a_n]}")
                    test_seq = test_seq[:ts_n-a_n]
                    break
            #remove polyG at end (inverted since this is a left clipped read)
            for n in range(1,len(test_seq)):
                if test_seq[-n] == 'G':
                    continue
                else:
                    if n > 3:
                        #print(f"clipped polyG RIGHT: {test_seq} -> {test_seq[:-n+1]}")
                        test_seq = test_seq[:-n]
                    break

            min_len = min([len(x) for x in (test_seq, clipped, expected)])
            if min_len < 4:
                self.right_match = NO_COVERAGE_OF_BREAKPOINT
                return NO_COVERAGE_OF_BREAKPOINT
            score_clipped = 0
            score_expected = 0
            for i in range(min_len):
                if test_seq[i] == clipped[i]:
                    score_clipped += 1
                else:
                    score_clipped -= 2
                if test_seq[i] == expected[i]:
                    score_expected += 1
                else:
                    score_expected -= 2
            # decision time
            self.right_match = self.classify_breakpoint(score_expected, score_clipped)
            if self.right_match is not None:
                return self.right_match
            DEBUG and print(f"detected wrong match side=R score_wt={score_expected}, score_clip={score_clipped}")
            DEBUG and print(f"measured seq: {test_seq}")
            DEBUG and print(f"clipped seq:  {clipped}")
            DEBUG and print(f"wt seq:       {expected}")
            self.right_match = MATCHES_NEITHER
            return MATCHES_NEITHER

        else:
            self.right_match = NO_COVERAGE_OF_BREAKPOINT
            return NO_COVERAGE_OF_BREAKPOINT

    def genotype_str(self):
        if self.left_match == MATCHES_WT and self.right_match == MATCHES_WT:
            return 'wt' #wild type
        if self.left_match == NO_COVERAGE_OF_BREAKPOINT and self.right_match == NO_COVERAGE_OF_BREAKPOINT:
            return 'nc' #no coverage
        if self.left_match == MATCHES_WT and self.right_match == NO_COVERAGE_OF_BREAKPOINT:
            return 'nc'
        if self.left_match == NO_COVERAGE_OF_BREAKPOINT and self.right_match == MATCHES_WT:
            return 'nc'
        if self.left_match == MATCHES_CLIP and self.right_match == MATCHES_CLIP:
            return 'id' #double insertion
        if self.left_match == MATCHES_CLIP:
            return 'il' #insertion-left
        if self.right_match == MATCHES_CLIP:
            return 'ir' #insertion-right
        if self.right_match == MATCHES_NEITHER:
            if self.left_match == MATCHES_NEITHER:
                return 'ad' #artefact: matches neither on either side
            else:
                return 'ar' #artefact right
        if self.left_match == MATCHES_NEITHER:
            return 'al' #artefact left
        DEBUG and print(f"found other genotype: r={self.right_match}, l={self.left_match}")
        return '??'

    def __str__(self):
        return f"read {self.read.query_name}: genotype={self.genotype_str()}"