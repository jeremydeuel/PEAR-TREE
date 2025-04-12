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
from typing import List, Iterator
import gzip
import os
class Insertion:
    def __init__(self, name: str, right_consensus: str, left_consensus: str, right_mates: List[str], left_mates: List[str], file: str):
        self.name = name
        self.chr, pos = name.split(":", maxsplit=2)
        self.right_pos, self.left_pos = pos.split("-", maxsplit=2)
        self.right_pos, self.left_pos = int(self.right_pos), int(self.left_pos)
        self.right_consensus = right_consensus
        self.left_consensus = left_consensus
        self.right_mates = right_mates
        self.left_mates = left_mates
        self.files = [os.path.basename(file)]

    @property
    def left_clipped(self) -> str:
        for i in range(len(self.left_consensus)):
            if self.left_consensus[i] in ("a","t","g","c","n"):
                continue
            return revcomp(self.left_consensus[:i])
        return ''


    @property
    def right_clipped(self) -> str:
        for i in range(len(self.right_consensus)):
            if self.right_consensus[i] in ("A","T","G","C","N"):
                continue
            return self.right_consensus[i:]
        return ''

    def __str__(self) -> str:
        output = f'>{self.name}\n@RIGHT_CONSENSUS\n{self.right_consensus}\n@RIGHT_MATES\n'
        output += '\n'.join(self.right_mates)
        output += f'\n@LEFT_CONSENSUS\n{self.left_consensus}\n@LEFT_MATES\n'
        output += '\n'.join(self.left_mates)
        output += f'\n@FILES\n'
        output += '\n'.join(self.files)
        output += '\n\n'
        return output

    def __iadd__(self, other: 'Insertion') -> 'Insertion':
        if len(self.left_clipped) < len(other.left_clipped):
            self.left_consensus = other.left_consensus
        if len(self.right_clipped) < len(other.right_clipped):
            self.right_consensus = other.right_consensus
        self.right_mates += other.right_mates
        self.left_mates += other.left_mates
        self.files += other.files
        return self

    @staticmethod
    def parseFile(path) -> Iterator['Insertion']:
        output = None
        state = None
        with gzip.open(path, 'rt') as f:
            for s in f:
                s = s.strip()
                if not len(s): continue
                if s[0] == '>':
                    if output: yield output
                    output = Insertion(s[1:], '', '', [], [], path)
                    state = None
                if s[0] == '@':
                    state = s[1:]
                    if state == 'FILES':
                        #reset files if such a field is provided
                        output.files = []
                    continue
                else:
                    if state is None:
                        continue
                    elif state == 'LEFT_CONSENSUS':
                        output.left_consensus = s
                    elif state == 'RIGHT_CONSENSUS':
                        output.right_consensus = s
                    elif state == 'LEFT_MATES':
                        output.left_mates.append(s)
                    elif state == 'RIGHT_MATES':
                        output.right_mates.append(s)
                    elif state == 'FILES':
                        output.files.append(s)
                    else:
                        raise ValueError(f"Unknown state {state} in file {f}")