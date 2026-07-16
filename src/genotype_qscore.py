from quality_seq import QualitySeq
LEFT_TO_RIGHT = 1
RIGHT_TO_LEFT = 2


def _pass(z):
    """Score one register: returns (ref_score, alt_score, art_score, total_quality)."""
    ref_score = alt_score = art_score = total_q = 0
    for s, q, r, a in z:
        if s == 'N':
            continue
        total_q += q
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
    return ref_score, alt_score, art_score, total_q


def qscore(seq: QualitySeq, ref: str, alt: str, direction: int = LEFT_TO_RIGHT):
    # z is the no-shift register; z2/z3 are the +/-1 bp shifts (built lazily so a clean
    # read never pays for them). z2 slides the read by one, z3 slides the consensus.
    S, Q = seq._sequence, seq._quality
    if direction is LEFT_TO_RIGHT:
        z1 = zip(S, Q, ref, alt)
        mk2 = lambda: zip(S[1:], Q[1:], ref, alt)
        mk3 = lambda: zip(S, Q, ref[1:], alt[1:])
    else:
        z1 = zip(reversed(S), reversed(Q), reversed(ref), reversed(alt))
        mk2 = lambda: zip(reversed(S[:-1]), reversed(Q[:-1]), reversed(ref), reversed(alt))
        mk3 = lambda: zip(reversed(S), reversed(Q), reversed(ref[:-1]), reversed(alt[:-1]))

    r1, a1, art1, total_q = _pass(z1)
    s1 = (r1, a1, art1)
    # Fast path: if the no-shift register already matches one hypothesis at EVERY scored
    # base (ref_score or alt_score == the summed quality), the alignment is perfect and
    # optimal for this register; the +/-1 bp shifts exist only to rescue an imperfect
    # alignment and cannot improve on a perfect one here (validated: golden call and the
    # genotyping unit tests are byte-identical). Skips 2 of 3 passes on clean reads.
    if r1 == total_q or a1 == total_q:
        return s1

    # ±1 bp register search to tolerate breakpoint imprecision. Pick the register that
    # best explains the read under EITHER hypothesis, scored as max(ref,alt), and report
    # THAT register's (ref, alt, art) jointly. Ties prefer the no-shift register s1 so a
    # shift only wins when it genuinely improves the alignment.
    s2 = _pass(mk2())[:3]
    s3 = _pass(mk3())[:3]
    best = s1
    best_key = max(r1, a1)
    for cand in (s2, s3):
        cand_key = max(cand[0], cand[1])
        if cand_key > best_key:
            best, best_key = cand, cand_key
    return best