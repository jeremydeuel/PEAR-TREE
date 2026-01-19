from quality_seq import QualitySeq
from math import log10
LEFT_TO_RIGHT = 1
RIGHT_TO_LEFT = 2

REF_MATCH = 1
ALT_MATCH = 2
ARTEFACT = 3

K = 10.0 ** -0.1

def qscore(seq: QualitySeq, ref: str, alt: str, direction:int=LEFT_TO_RIGHT):
    # create two zipped objects for comparison
    # z is without offset
    # z2 is with an offset of 1
    if direction is LEFT_TO_RIGHT:
        z = zip(seq._sequence, seq._quality, ref, alt)
        z2 = zip(seq._sequence[1:], seq._quality[1: ], ref, alt)
        z3 = zip(seq._sequence, seq._quality, ref[1:], alt[1:])
    else:
        z = zip(reversed(seq._sequence), reversed(seq._quality), reversed(ref), reversed(alt))
        z2 = zip(reversed(seq._sequence[:-1]), reversed(seq._quality[:-1]), reversed(ref), reversed(alt))
        z3 = zip(reversed(seq._sequence), reversed(seq._quality), reversed(ref[:-1]), reversed(alt[:-1]))
    ref_score, alt_score, art_score = 0, 0, 0
    for s, q, r, a in z:
        if s == 'N': continue
        if s == r:
            ref_score += q
        else:
            ref_score -= q
        if s == a:
            alt_score += q
        else:
            alt_score -= q
        if s != a and s != r:
            art_score += q
    s1 = ref_score, alt_score, art_score
    ref_score, alt_score, art_score = 0, 0, 0
    for s, q, r, a in z2:
        if s == 'N': continue
        if s == r:
            ref_score += q
        else:
            ref_score -= q
        if s == a:
            alt_score += q
        else:
            alt_score -= q
        if s != a and s != r:
            art_score += q
    s2 = ref_score, alt_score, art_score
    ref_score, alt_score, art_score = 0, 0, 0
    for s, q, r, a in z3:
        if s == 'N': continue
        if s == r:
            ref_score += q
        else:
            ref_score -= q
        if s == a:
            alt_score += q
        else:
            alt_score -= q
        if s != a and s != r:
            art_score += q
    s3 = ref_score, alt_score, art_score
    if s2[0]>s1[0]  and s2[0] > s2[2]:
        #print(f"swapped {s2} for {s1} for better ref read (-1)")
        return s2
    if s2[1]>s1[1]  and s2[1] > s2[2]:
        #print(f"swapped {s2} for {s1} for better alt read (-1)")
        return s2
    if s3[0]>s1[0]  and s3[0] > s3[2]:
        #print(f"swapped {s3} for {s1} for better ref read (+1)")
        return s3
    if s3[1]>s1[1]  and s3[1] > s3[2]:
        #print(f"swapped {s3} for {s1} for better alt read (+1)")
        return s3
    return s1


def qscore_previous(seq: QualitySeq, ref: str, alt: str, direction=LEFT_TO_RIGHT):


    #set prior probabilities
    q_ref = 0
    q_alt = 0
    q_art = 0
    if direction == RIGHT_TO_LEFT:
        r = range(-min([len(ref),len(alt), len(seq)]), 0)
    else:
        r = range(min([len(ref),len(alt),len(seq)]))
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