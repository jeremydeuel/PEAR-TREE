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

# helper function to find adapters

from config import CONFIG
from revcomp import revcomp


def is_adapter(seq: str) -> bool:
    """
    Checks if seq has an adapter sequence. seq has to have exactly MIN_CLIP_LEN bases
    left-cipped input has to be reverse complemented!
    """
    seq = seq.upper()
    if len(seq) < 2: return False
    # remove all homopolymers except polyT (which corresponds to a polyA tail)
    for b in ('G', 'C', 'A'):
        if b * CONFIG['discovery']['min_clip_len'] == seq:
            return True
    if seq[0] != seq[1]: #if the first two bases dont match:
        if seq[0:2] * int(CONFIG['discovery']['min_clip_len'] / 2) == seq[:(int(CONFIG['discovery']['min_clip_len'] / 2) * 2)]:
            return True #repetitive dimer
    # remove adapter sequences
    revseq = revcomp(seq)
    for a in CONFIG['adapters']:
        if seq in a or revseq in a:
            return True
    return False


def clean_clipped_seq(seq: str) -> str:
    # remove poly G at tail of read
    if len(seq) < 2:
        return ''
    i = len(seq)
    while i > 1 and seq[i-1] == 'G':
        i -= 1
    seq = seq[:i]
    #remove adapter sequence, use the first of the REMOVE_ADAPTERS sequence for this.
    i = 0
    j = 0
    while i < len(seq):
        if seq[i] == CONFIG['adapters'][0][j]:
            j += 1
        else:
            if j > 0: #give j=0 another chance in case the adapter starts exactly here
                j = 0
                continue
        i += 1
        #if more than 9 bases match, abort this and cut of.
        if j > 9:
            break
    if j >= CONFIG['discovery']['min_adapterlen_for_clip']:
        seq = seq[:i-j]
    return seq

