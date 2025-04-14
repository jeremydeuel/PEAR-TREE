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


import gzip
import pysam
from genotyping_evidence_read import EvidenceRead
from typing import Iterator, List
from genotype_qscore import REF_MATCH, ALT_MATCH, ARTEFACT
from math import log10

GT_ARTEFACT = 'artefact'
GT_WILDTYPE = 'wild-type'
GT_HETEROZYGOUS = 'heterozygous'
GT_HOMOZYGOUS = 'homozygous'
class Insertion:
    def __init__(self, name: str):
        self.name = name
        self.chr, pos = name.split(":", maxsplit=2)
        self.left_pos, self.right_pos = pos.split("-", maxsplit=2)
        self.right_pos, self.left_pos = int(self.right_pos), int(self.left_pos)
        self.left_clipped = None
        self.right_clipped = None
        self.left_ref = None
        self.right_ref = None
        self.evidence_reads: List[EvidenceRead] = []
    def __str__(self):
        return self.name

    @staticmethod
    def import_file(path: str) -> Iterator['Insertion']:
        insertion = None
        with gzip.open(path, 'rt') as f:
            for line in f:
                line = line.strip()
                if not len(line): continue
                if line[0] == ">":
                    if insertion: yield insertion
                    insertion = Insertion(line[1:])
                    status = None
                    continue
                if line[0] == "@":
                    status = line[1:]
                    continue
                if status == 'RIGHT_INSERTION':
                    insertion.right_clipped = line
                elif status == 'LEFT_INSERTION':
                    insertion.left_clipped = line
                elif status == 'RIGHT_REFERENCE':
                    insertion.right_ref = line
                elif status == 'LEFT_REFERENCE':
                    insertion.left_ref = line
        if insertion: yield insertion
    def genotype(self, bam: pysam.AlignmentFile):
        start = min(self.left_pos, self.right_pos)-1
        end = max(self.left_pos, self.right_pos)+1
        # remove qname duplicates
        for read in bam.fetch(self.chr, start, end):
            if read.mapq < 40: continue
            if read.is_secondary: continue
            if read.is_unmapped: continue
            if read.is_qcfail: continue
            if read.reference_name != self.chr: continue
            evi_read = EvidenceRead(read)
            if self.left_pos: evi_read.qleft(self.left_pos, self.left_ref, self.left_clipped)
            if self.right_pos: evi_read.qright(self.right_pos, self.right_ref, self.right_clipped)
            self.evidence_reads.append(evi_read)
        return

    def summarise_evidence(self):
        if not len(self.evidence_reads):
            print(f"omitting {self.name}: no reads found")
            return GT_WILDTYPE, 0, 0
        # calculate best call using likelihood ratios
        q_art = 0
        q_alt = 0
        q_ref = 0
        q_hom = 0
        #print(f"evidence for {self.name}:")
        for er in self.evidence_reads:
            left_q_ref, left_q_alt, left_q_art = er.left_genotype
            right_q_ref, right_q_alt, right_q_art = er.right_genotype
            if max(left_q_ref,left_q_alt) <= left_q_art:
                q_art += left_q_art-max(left_q_ref, left_q_alt)
                continue
            if max(right_q_ref, right_q_alt) <= right_q_art:
                q_art += right_q_art-max(right_q_ref, right_q_alt)
                continue
            if right_q_art > 0 and left_q_art > 0:
                if right_q_ref >= right_q_alt and right_q_ref > 0:
                    if left_q_ref >= left_q_alt and left_q_ref > 0:
                        #double wt
                        q_hom -= 1
                        q_ref += right_q_ref-right_q_alt + left_q_ref-right_q_alt
                    elif left_q_alt > 0:
                        q_alt += left_q_alt-left_q_ref
                        q_ref += right_q_ref-right_q_alt
                elif right_q_alt > 0:
                    if left_q_ref >= left_q_alt and left_q_ref > 0:
                        q_alt += right_q_alt-right_q_ref
                        q_ref += left_q_ref-left_q_alt
                    elif left_q_alt > 0:
                        #double alt, this is an artefact
                        q_art += left_q_alt-left_q_ref + right_q_alt - right_q_ref
            elif right_q_art > 0:
                if right_q_ref >= right_q_alt and right_q_ref > 0:
                    q_ref += right_q_ref-right_q_alt
                elif right_q_alt > 0:
                    q_alt += right_q_alt-right_q_ref
            elif left_q_art > 0:
                if left_q_ref >= left_q_alt and left_q_ref > 0:
                    q_ref += left_q_ref-left_q_alt
                elif left_q_alt > 0:
                    q_alt += left_q_alt-left_q_ref

        if q_art >= max(q_ref,q_alt):
            return GT_ARTEFACT, q_art, max(q_ref, q_alt)
        if q_alt > q_art:
            if q_hom < 0:
                if q_ref >= q_art:
                    return GT_HETEROZYGOUS, q_alt, q_ref
                else:
                    return GT_ARTEFACT, q_art, max(q_alt, q_ref)
            if q_alt >= q_ref:
                return GT_HOMOZYGOUS, q_alt, q_ref
            else:
                return GT_HETEROZYGOUS, q_alt, q_ref
        if q_ref > q_art:
            return GT_WILDTYPE, q_ref, q_alt
        return GT_ARTEFACT, 0, 0
