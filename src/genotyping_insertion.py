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
from typing import Iterator

class Insertion:
    def __init__(self, name: str):
        self.name = name
        self.chr, pos = name.split(":", maxsplit=2)
        self.right_pos, self.left_pos = pos.split("-", maxsplit=2)
        self.right_pos, self.left_pos = int(self.right_pos), int(self.left_pos)
        self.left_clipped = None
        self.right_clipped = None
        self.left_ref = None
        self.right_ref = None
        self.evidence = {'wt':0,'id':0,'il':0,'ir':0,'ad':0,'al':0,'ar':0, '??': 0, 'hc': 0, 'nc': 0}
        self.evidence_reads = []
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
            if read.is_supplementary: continue
            if read.is_qcfail: continue
            if read.reference_name != self.chr: continue
            evi_read = EvidenceRead(read)
            evi_read.check_left(self.left_pos, self.left_clipped, self.left_ref)
            evi_read.check_right(self.right_pos, self.right_clipped, self.right_ref)
            self.evidence_reads.append(evi_read)
        return

    def summarise_evidence(self):
        genotypes = {}
        for er in self.evidence_reads:
            gt = er.genotype_str()
            if gt not in genotypes.keys():
                genotypes[gt] = set()
            genotypes[gt].add(er.read.query_name)
        for gt, queries in genotypes.items():
            self.evidence[gt] += len(queries)


    @property
    def genotype_string(self):
        if self.evidence['hc']:
            return 'NA'
        good_reads = self.evidence['wt']+self.evidence['ir']+self.evidence['il']
        artefact_reads = self.evidence['ar']+self.evidence['al']+self.evidence['ad']
        if artefact_reads > 1+int(good_reads/10): #number of artefact reads is more than 10% of the number of good reads.
            return 'artefact'
        if self.evidence['id']>0:
            return 'artefact'
        if self.evidence['wt'] and not self.evidence['ir'] and not self.evidence['il']:
            return 'wild-type'
        if self.evidence['ir'] and self.evidence['il'] and not self.evidence['wt']:
            return 'homozygous'
        if self.evidence['ir'] and self.evidence['il'] and self.evidence['wt']:
            return 'heterozygous'
        if self.evidence['il'] or self.evidence['ir']:
            #figure out a few edge cases
            insertion_sum = self.evidence['il'] + self.evidence['ir']
            if insertion_sum > max(8,8* self.evidence['wt']): #8 times more insertion reads, but not both ends covered, thats suspicious
                return 'artefact'
            if self.evidence['wt'] > max(8,insertion_sum * 8): #8 times more wt reads than insertion reads, thats suspicious
                return 'artefact'
            return 'insertion?'
        return 'NA'