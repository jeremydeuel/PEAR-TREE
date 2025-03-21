from config import CONFIG
from typing import List
def is_low_complexity(seq: str) -> bool:
    """
    checks if a sequence is low complexity
    this is done by counting the number of different bases. If only two or less of the four possible bases
    are found, returns False, else True
    """
    bases = {'A':0,'T':0,'G':0,'C':0,'N':0}
    for s in seq:
        bases[s] += 1
    del bases['N']
    missing_bases = 0
    summed_bases = 0
    for base, n in bases.items():
        if n == 0:
            missing_bases += 1
        summed_bases += n
    if missing_bases > 1:
        return True
    return False


def has_well_defined_breakpoint(left_seq: str, right_seq: str) -> bool:
    """
    This function looks at the sequences of the clipped and unclipped sequence and determines if they are so repetitive
    that the breakpoint can not be determined safely.
    """
    if is_low_complexity(left_seq[-10:]) and is_low_complexity(right_seq[10:]):
        return False
    # check for repeating n-mers
    full_seq = left_seq[-10:] + right_seq[:10]
    for n in range(2,6):
        repetitive_sequence = full_seq[:n]*max(2,1+int(11/n))
        if repetitive_sequence in full_seq:
            return False
        repetitive_sequence = full_seq[-n:] * max(2, 1 + int(11 / n))
        if repetitive_sequence in full_seq:
            return False
    return True


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

def sequence_matching_score(seqs: List[str]) -> int:
    min_len = min(32, min([len(s) for s in seqs]))
    score = 0
    for i in range(min_len):
        bases = {
            'a': 0,
            't': 0,
            'g': 0,
            'c': 0,
            'n': 0
        }
        for s in seqs:
            bases[s[i]] += 1
        bases['n'] = 0
        if max(bases.values()) >= sum(bases.values())*0.51:
            score += 1
        else:
            score -= 2
    return score/min_len