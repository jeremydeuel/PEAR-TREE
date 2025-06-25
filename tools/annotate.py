# PEAR-TREE - Genomic events of mutation through insertional alterations by transposable elements
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

import pandas as pd
import gzip
import os
import pysam
import pyliftover


#HUMAN
EXCEL_IN = "../../spar/human/insertions.colon.xlsx"
EXCEL_OUT = "../../spar/human/insertions.colon.annot.xlsx"
INSERTIONS_FILE = '../../spar/human/insertions.colon.combined.txt.gz'
TMP_FASTA = '../../spar/human/insertions.colon.fa.gz'
DFAM_SCRIPT = '../../genomes/dfam/dfamscan.pl'
DFAM_HMM = '../../genomes/dfam/Dfam_hs.hmm'
TMP_DFAM = '../../spar/human/insertions.colon.dfam'
TMP_SAM = '../../spar/human/insertions.colon.sam'
BOWTIE_SCRIPT = '/opt/homebrew/bin/bowtie2'
BOWTIE_REF = '../../genomes/bowtie2_indices/hs1'
CHAINFILE = '../../genomes/hs1.hg38.all.chain.gz'
#MOUSE
EXCEL_IN = "../../spar/mouse/germline.insertions.xlsx"
EXCEL_OUT = "../../spar/mouse/germline.insertions.annot.xlsx"
INSERTIONS_FILE = '../../spar/mouse/discovery.combined.txt.gz'
TMP_FASTA = '../../spar/mouse/annotation.germline.fa.gz'
DFAM_SCRIPT = '../../genomes/dfam/dfamscan.pl'
DFAM_HMM = '../../genomes/dfam/Dfam_mm.hmm'
TMP_DFAM = '../../spar/mouse/annotation.germline.dfam'
TMP_SAM = '../../spar/mouse/annotation.germline.bam'
BOWTIE_SCRIPT = '/opt/homebrew/bin/bowtie2'
BOWTIE_REF = '../../genomes/bowtie2_indices/mm39'
CHAINFILE = None#


if __name__ == '__main__':
    if CHAINFILE is not None:
        print(f"reading chainfile {CHAINFILE}, this might take a while...")
        lo = pyliftover.LiftOver(CHAINFILE)
    else:
        lo = None
    print("done with reading chainfile.")
    d = pd.read_excel(EXCEL_IN)
    print(d)
    seqs = {i: {'L': None, 'R': None} for i in list(d['insertion'])}
    with gzip.open(INSERTIONS_FILE, 'rt') as ifh:
        with gzip.open(TMP_FASTA, 'wt') as tmpfah:
            state = None
            insertion = None
            for line in ifh:
                line = line.strip()
                if not line: continue
                if line[0] == ">":
                    if line[1:] in seqs.keys():
                        insertion = line[1:]
                        state = None
                    else:
                        insertion = None
                if insertion is None:
                    continue
                if line[0] == "@":
                    state = line[1:]
                    continue
                if state == 'RIGHT_CONSENSUS':
                    insertion_start = min([line.find(b) for b in 'acgt' if b in line])
                    tmpfah.write(f">{insertion}:L\n{line[insertion_start:]}\n\n")
                    seqs[insertion]['L'] = line[insertion_start:]
                if state == 'LEFT_CONSENSUS':
                    insertion_end = min([line.find(b) for b in 'ACGT' if b in line])
                    tmpfah.write(f">{insertion}:R\n{line[:insertion_end]}\n\n")
                    seqs[insertion]['R'] = line[:insertion_end]
    if not os.path.exists(TMP_DFAM):
        print("running dfam")
        print(os.system(f"{DFAM_SCRIPT} --fastafile {TMP_FASTA} --hmmfile {DFAM_HMM} --cpu {os.cpu_count()} --dfam_outfile {TMP_DFAM}"))
    if not os.path.exists(TMP_SAM):
        print('running bowtie2')
        print(os.system(f"{BOWTIE_SCRIPT} -x {BOWTIE_REF} --end-to-end -f {TMP_FASTA} > {TMP_SAM}"))
    dfams = {i: {'L': None, 'R': None} for i in list(d['insertion'])}

    #annotate dfam
    with open(TMP_DFAM, 'r') as dfam:
        for line in dfam:
            line = line.strip()
            if not line: continue
            if line[0] == '#': continue
            line = [s for s in line.split(" ") if len(s)]
            insertion = line[2][:-2]
            if dfams[insertion][line[2][-1]] is None:
                dfams[insertion][line[2][-1]] = line[0], float(line[3])
            else:
                if float(line[3]) > dfams[insertion][line[2][-1]][1]:
                    dfams[insertion][line[2][-1]] = line[0], float(line[3])
    for insertion, dfam_sides in dfams.items():
        for side, dfam in dfam_sides.items():
            if dfam is None:
                if side == 'L':
                    if seqs[insertion][side][:6] == 'tttttt':
                        dfams[insertion][side] = ('polyA', 0)
                if side == 'R':
                    if seqs[insertion][side][-6:] == 'aaaaaa':
                        dfams[insertion][side] = ('polyA', 0)
    mappings = {i: {'L': None, 'R': None} for i in list(d['insertion'])}
    hg38 = {i: {'L': None, 'R': None} for i in list(d['insertion'])}
    with pysam.AlignmentFile(TMP_SAM) as sam:
        for read in sam:
            if read.is_qcfail: continue
            if read.is_secondary: continue
            if read.is_supplementary: continue
            if read.is_unmapped: continue
            if read.mapq >= 40:
                mappings[read.query_name[:-2]][read.query_name[
                    -1]] = f"{read.reference_name}:{read.reference_start}-{read.reference_end}{'+' if read.is_forward else '-'}"
                if ( read.query_name[-1] == "R" ) ^ read.is_forward:
                    if lo is not None:
                        co = lo.convert_coordinate(read.reference_name, read.reference_start, '+' if read.is_forward else '-')
                    else:
                        co = [(read.reference_name, read.reference_start, '+' if read.is_forward else '-')]
                else:
                    if lo is not None:
                        co = lo.convert_coordinate(read.reference_name, read.reference_end, '+' if read.is_forward else '-')
                    else:
                        co = [(read.reference_name, read.reference_end, '+' if read.is_forward else '-')]
                if co:
                    hg38[read.query_name[:-2]][read.query_name[
                    -1]] = f"{co[0][0]}:{co[0][1]}{co[0][2]}"

    for i in range(d.shape[0]):
        d.loc[i,'left_seq'] = seqs[d['insertion'][i]]['L']
        d.loc[i, 'right_seq'] = seqs[d['insertion'][i]]['R']
        d.loc[i, 'left_dfam'] = dfams[d['insertion'][i]]['L'][0] if dfams[d['insertion'][i]]['L'] is not None else pd.NA
        d.loc[i, 'right_dfam'] = dfams[d['insertion'][i]]['R'][0] if dfams[d['insertion'][i]]['R'] is not None else pd.NA
        d.loc[i, 'left_align'] = mappings[d['insertion'][i]]['L'] if mappings[d['insertion'][i]][
                                                                         'R'] is not None else pd.NA
        d.loc[i, 'right_align'] = mappings[d['insertion'][i]]['R'] if mappings[d['insertion'][i]][
                                                                         'R'] is not None else pd.NA
        d.loc[i, 'hg38_ins_bp_left'] = hg38[d['insertion'][i]]['L'] if hg38[d['insertion'][i]][
                                                                         'L'] is not None else pd.NA
        d.loc[i, 'hg38_ins_bp_right'] = hg38[d['insertion'][i]]['R'] if hg38[d['insertion'][i]][
                                                                          'R'] is not None else pd.NA
    print(d)
    d.to_excel(EXCEL_OUT)



