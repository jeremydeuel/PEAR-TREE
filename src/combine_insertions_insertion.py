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


from revcomp import revcomp
from typing import List, Iterator, Dict
from quality_seq import QualitySeq
import gzip
import os

TYPE_FULL_INFO = 3
TYPE_RIGHT_POLYA = 1
TYPE_LEFT_POLYA = 2
class Insertion:
    def __init__(self, reference_name: str, start: str, end: str, data: Dict, file: str):
        assert reference_name is not None and len(reference_name), f"reference_name not given."
        self.name = f"{reference_name}:{start}-{end}"
        self.type = None
        if 'RIGHT:ALIGNED' in data.keys() and 'RIGHT:CLIPPED' in data.keys():
            self.right_aligned = data['RIGHT:ALIGNED'].revcomp()
            self.right_clipped = data['RIGHT:CLIPPED']
            self.right_pos = int(end)
        else:
            assert end[:6] == 'polyA_', f"false right poly A detected in {reference_name}:{start}-{end} -> {data}"
            self.right_clipped = data['RIGHT:CLIPPED_POLYA']
            self.right_aligned = None
            self.right_pos = None
            self.type = TYPE_RIGHT_POLYA
        if 'LEFT:ALIGNED' in data.keys() and 'LEFT:CLIPPED' in data.keys():
            self.left_clipped = data['LEFT:CLIPPED'].revcomp()
            self.left_aligned =  data['LEFT:ALIGNED']
            self.left_pos = int(start)
        else:
            assert start[:6] == 'polyA_' , f"false left poly A detected in {reference_name}:{start}-{end} -> {data}"
            self.left_clipped = data['LEFT:CLIPPED_POLYA'].revcomp().lower()
            self.left_aligned = None
            self.left_pos = None
            self.type = TYPE_LEFT_POLYA
        if self.type is None:
            self.type = TYPE_FULL_INFO
        self.reference_name = reference_name
        self.right_mates = data['RIGHT:MATE']
        self.left_mates = data['LEFT:MATE']
        self.files = [os.path.basename(file)]

    @property
    def left_consensus(self) -> str:
        assert self.type is not TYPE_LEFT_POLYA, f"can not extract consensus from left polyA type"
        return self.left_clipped.revcomp().lower() + self.left_aligned

    @property
    def right_consensus(self) -> str:
        assert self.type is not TYPE_RIGHT_POLYA, f"can not extract consensus from right polyA type"
        return self.right_aligned.revcomp() + self.right_clipped.lower()

    def __str__(self) -> str:
        output = self.right_clipped.fastq(f"{self.name}:RIGHT:CLIPPED")
        if type is not TYPE_RIGHT_POLYA:
            output += self.right_consensus.fastq(f"{self.name}:RIGHT:ALIGNED")
        output += self.left_clipped.fastq(f"{self.name}:LEFT:CLIPPED")
        if type is not TYPE_LEFT_POLYA:
            output += self.left_consensus.fastq(f"{self.name}:LEFT:ALIGNED")
        i=0
        for m in self.left_mates:
            i += 1
            output += m.fastq(f"{self.name}:LEFT:MATE{i}")
        i=0
        for m in self.right_mates:
            i+= 1
            output += m.fastq(f"{self.name}:RIGHT:MATE{i}")
        return output

    def __iadd__(self, other: 'Insertion') -> 'Insertion':
        if self.type is TYPE_LEFT_POLYA and other.type is not TYPE_LEFT_POLYA:
            self.left_clipped = other.left_clipped
            self.left_aligned = other.left_aligned
            self.left_pos = other.left_pos
            self.name = f"{self.reference_name}:{self.right_pos}-{self.left_pos}"
            self.type = TYPE_FULL_INFO
            self.left_mates.append(QualitySeq('tttttttttttt',[0]*12) + self.left_clipped)
        elif self.type is not TYPE_LEFT_POLYA and other.type is TYPE_LEFT_POLYA:
            self.left_mates.append(QualitySeq('tttttttttttt',[0]*12) + other.left_clipped)
        else:
            if len(self.left_clipped) < len(other.left_clipped):
                self.left_clipped = other.left_clipped
            if self.type is not TYPE_LEFT_POLYA and len(self.left_aligned) < len(other.left_aligned):
                self.left_aligned = other.left_aligned
        if self.type is TYPE_RIGHT_POLYA and other.type is not TYPE_RIGHT_POLYA:
            self.right_clipped = other.right_clipped
            self.right_aligned = other.right_aligned
            self.right_pos = other.right_pos
            self.name = f"{self.reference_name}:{self.right_pos}-{self.left_pos}"
            self.right_mates.append(QualitySeq('tttttttttttt',[0]*12) + self.right_clipped)
            self.type = TYPE_FULL_INFO
        elif self.type is not TYPE_RIGHT_POLYA and other.type is TYPE_RIGHT_POLYA:
            self.right_mates.append(QualitySeq('tttttttttttt',[0]*12) + other.right_clipped)
        else:
            if len(self.right_clipped) < len(other.right_clipped):
                self.right_clipped = other.right_clipped
            if self.type is not TYPE_RIGHT_POLYA and len(self.right_aligned) < len(other.right_aligned):
                self.right_aligned = other.right_aligned
        self.right_mates += other.right_mates
        self.left_mates += other.left_mates
        self.files += other.files
        return self

    @staticmethod
    def parseFile(path) -> Iterator['Insertion']:
        currentInsertion = None
        currentData = {'LEFT:MATE': [], 'RIGHT:MATE': []}
        with gzip.open(path, 'rt') as f:
            while l := f.readline():
                l = l.strip()
                if not len(l): continue
                reference_name, positions, side, field = l[1:].split(":")
                start, end = positions.split("-")
                insertion = (reference_name, start, end)
                if insertion != currentInsertion:
                    if currentInsertion is not None:
                        reference_name, start, end = currentInsertion
                        yield Insertion(reference_name, start, end, currentData, os.path.basename(path))
                    currentInsertion = insertion
                    currentData = {'LEFT:MATE':[],'RIGHT:MATE':[]}
                assert l[0] == '@', f"Expected label field, got {l}"
                seq = f.readline().strip()
                plus = f.readline().strip()
                assert plus == '+', f"Expected + separator, got {plus}"
                qual = [ord(i)-33 for i in f.readline().strip()]
                qs = QualitySeq(seq, qual)
                if "MATE" in field:
                    currentData[f"{side}:MATE"].append(qs)
                else:
                    currentData[f"{side}:{field}"] = qs
            if currentInsertion is not None:
                reference_name, start, end = currentInsertion
                yield Insertion(reference_name, start, end, currentData, os.path.basename(path))
