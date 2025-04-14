from quality_seq import QualitySeq
from math import log10
LEFT_TO_RIGHT = 1
RIGHT_TO_LEFT = 2

REF_MATCH = 1
ALT_MATCH = 2
ARTEFACT = 3

K = 10.0 ** -0.1
def qscore(seq: QualitySeq, ref: str, alt: str, direction=LEFT_TO_RIGHT):


    #set prior probabilities
    q_ref = 0
    q_alt = 0
    q_art = 0
    if direction == RIGHT_TO_LEFT:
        r = range(-min(len(ref),len(alt), len(seq)), 0)
    else:
        r = range(min(len(ref),len(alt),len(seq)))
    for i in r:
        base, qual = seq.getPos(i)
        if base == "N": continue
        if ref[i] == alt[i]:
            if base == ref[i]:
                q_ref += qual
                q_alt += qual
                q_art += qual /4
            else:
                q_ref -= qual
                q_alt -= qual
                q_art += qual/4
        else:
            if base == ref[i]:
                q_ref += qual
                q_alt -= qual
                q_art += qual/4
            elif base == alt[i]:
                q_alt += qual
                q_ref -= qual
                q_art += qual/4
            else:
                q_ref -= qual
                q_alt -= qual
                q_art += qual/4
    return q_ref, q_alt, q_art
if __name__ == '__main__':
    measured = QualitySeq('ACGGTTTTTTTTTTTTT',[20,20,20,20,20,20,3,20,3,20,3,20,3,3,20,3,20])
    alt = 'TTTTTTTTTTTT'
    ref = 'TATGTCTCGTCT'
    #ref = 'TTATCTTTTTGT'
    print(qscore(measured, ref, alt, RIGHT_TO_LEFT))