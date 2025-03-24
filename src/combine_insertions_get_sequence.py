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


import py2bit
from config import CONFIG

GENOME = py2bit.open(CONFIG['combine_insertions']['genome_2bit'])

def get_sequence(seqname, start, end):
    """
    Wrapper for py2bit, renaming seqname if UCSC / NCBI mismatch
    :param seqname: seqname, can be either "chr1" or "1"
    :param start: start position, 0-based
    :param end: end-position, 0-based
    :return: sequence as string.
    """
    seqnames = GENOME.chroms().keys()
    if seqname not in seqnames:
        if seqname == 'MT':
            seqname = 'chrM'
        else:
            new_seqname = f"chr{seqname}"
            if new_seqname not in seqnames:
                if len(seqname)>2:
                    if seqname[-2:] == '.1':
                        seqname = seqname[:-2]
                    for s in seqnames:
                        if seqname in s:
                            seqname = s
                            break
            else:
                seqname = new_seqname
    if seqname not in seqnames:
        print(f"Seqname {seqname} not found in file {CONFIG['combine_insertions']['genome_2bit']}")
        return ''
    if start>=end: return ''
    try:
        return GENOME.sequence(seqname, start, end)
    except Exception as e:
        print(f'Failure to extract {seqname}:{start}-{end} -> {e}')
        return ''
