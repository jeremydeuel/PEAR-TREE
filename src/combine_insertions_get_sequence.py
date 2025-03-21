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
