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


def revcomp(seq: str):
    seq = [b for b in seq]
    seq.reverse()
    for i in range(len(seq)):
        if seq[i] == "A":
            seq[i] = "T"
        elif seq[i] == "T":
            seq[i] = "A"
        elif seq[i] == "G":
            seq[i] = "C"
        elif seq[i] == "C":
            seq[i] = "G"
        elif seq[i] == 'a':
            seq[i] = 't'
        elif seq[i] == 't':
            seq[i] = 'a'
        elif seq[i] == 'g':
            seq[i] = 'c'
        elif seq[i] == 'c':
            seq[i] = 'g'
    return "".join(seq)