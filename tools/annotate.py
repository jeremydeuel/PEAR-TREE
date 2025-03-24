import pandas as pd
import gzip, os
import pysam

EXCEL_IN = "all.insertions.xlsx"
EXCEL_OUT = "all.insertions.annot.xlsx"
INSERTIONS_FILE = 'all_insertions.txt.gz'
TMP_FASTA = 'all.insertions.fa.gz'
DFAM_SCRIPT = '/Users/jeremy/Desktop/isa_all/dfam/dfamscan.pl'
DFAM_HMM = '/Users/jeremy/Desktop/isa_all/dfam/Dfam_hs.hmm'
TMP_DFAM = 'all.insertions.dfam'
TMP_SAM = 'all.insertions.sam'
BOWTIE_SCRIPT = '/opt/homebrew/bin/bowtie2'
BOWTIE_REF = '/Users/jeremy/Desktop/MSC/hs1/hs1'


if __name__ == '__main__':
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

    with pysam.AlignmentFile(TMP_SAM) as sam:
        for read in sam:
            if read.is_qcfail: continue
            if read.is_secondary: continue
            if read.is_supplementary: continue
            if read.is_unmapped: continue
            if read.mapq > 40:
                mappings[read.query_name[:-2]][read.query_name[-1]] = f"{read.reference_name}:{read.reference_start}-{read.reference_end}{'+' if read.is_forward else '-'}"


    print(seqs)
    print(dfams)
    print(mappings)

    for i in range(d.shape[0]):
        d.loc[i,'left_seq'] = seqs[d['insertion'][i]]['L']
        d.loc[i, 'right_seq'] = seqs[d['insertion'][i]]['R']
        d.loc[i, 'left_dfam'] = dfams[d['insertion'][i]]['L'][0] if dfams[d['insertion'][i]]['L'] is not None else pd.NA
        d.loc[i, 'right_dfam'] = dfams[d['insertion'][i]]['R'][0] if dfams[d['insertion'][i]]['R'] is not None else pd.NA
        d.loc[i, 'left_align'] = mappings[d['insertion'][i]]['L'] if mappings[d['insertion'][i]][
                                                                         'R'] is not None else pd.NA
        d.loc[i, 'right_align'] = mappings[d['insertion'][i]]['R'] if mappings[d['insertion'][i]][
                                                                         'R'] is not None else pd.NA
    print(d)
    d.to_excel(EXCEL_OUT)



