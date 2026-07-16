from quality_seq import QualitySeq
LEFT_TO_RIGHT = 1
RIGHT_TO_LEFT = 2


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
    # ±1 bp register search to tolerate breakpoint imprecision. Pick the register
    # that best explains the read under EITHER hypothesis, scored as max(ref,alt),
    # and report THAT register's (ref, alt, art) jointly. The previous ladder
    # tested ref before alt and returned the first register that improved either,
    # which (a) resolved ties toward ref and (b) reported an artefact score from a
    # register chosen to maximise ref/alt — systematically suppressing artefact
    # evidence. Ties prefer the no-shift register s1 so a shift only wins when it
    # genuinely improves the alignment.
    best = s1
    best_key = max(s1[0], s1[1])
    for cand in (s2, s3):
        cand_key = max(cand[0], cand[1])
        if cand_key > best_key:
            best, best_key = cand, cand_key
    return best