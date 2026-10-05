"""Generate consensus_cases.json.gz: python indel_aware_consensus results for consensus.rs.

  VENV=/private/tmp/claude-501/-Users-jeremy-Documents-PEAR-TREE/fde0700f-e325-4651-8daf-0cdd52bd072b/scratchpad/venv/bin/python
  $VENV rust/peartree-combine/tests/fixtures/gen_consensus_fixture.py N_RANDOM SEED OUT.json.gz

Cases = the deterministic read sets of test/test_indel_consensus.py + N_RANDOM randomised
sets (poly-A/T jitter, substitutions, indels, N's, lowercase, constant qualities -> weight
ties, two-haplotype ties, floating mates in both orientations, junk mates, shared fragment
groups with 1/len(cluster) weights, empty reads, > 40 anchored reads for probe sampling).
Each read is [seq, qual as chr(q+33), group, weight, anchored].
"""
import gzip
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "test"))

from indel_consensus import ClipRead, indel_aware_consensus, revcomp  # noqa: E402


def dump_result(r):
    return {"seq": r.seq, "depth": r.depth, "score": r.score, "stop_reason": r.stop_reason,
            "polya_base": r.polya_base, "polya_start": r.polya_start,
            "polya_len_median": r.polya_len_median, "polya_len_min": r.polya_len_min,
            "polya_len_max": r.polya_len_max, "beyond_polya": r.beyond_polya,
            "beyond_start": r.beyond_start, "beyond_polya_support": r.beyond_polya_support}


def case(name, reads, min_depth=2, polya_min_len=8):
    gmap = {}
    for r in reads:
        gmap.setdefault(r.group, len(gmap) * 7 + 3)   # non-dense ids on purpose
    reads = [ClipRead(r.seq, list(r.qual), gmap[r.group], r.weight, r.anchored) for r in reads]
    res = indel_aware_consensus(reads, min_depth=min_depth, polya_min_len=polya_min_len)
    return {"name": name, "min_depth": min_depth, "polya_min_len": polya_min_len,
            "reads": [[r.seq, "".join(chr(q + 33) for q in r.qual), r.group, r.weight, r.anchored]
                      for r in reads],
            "result": dump_result(res)}


def deterministic_cases():
    import test_indel_consensus as T
    out = []
    for seed in range(25):
        out.append(case(f"synth_reads_{seed}", T.synth_reads(seed)))
    rng = random.Random(7)
    reads = []
    for i in range(6):
        s = T.ELEM + "A" * (T.POLYA + rng.randint(-2, 2)) + (T.BEYOND if i == 0 else T.BEYOND[:3])
        reads.append(ClipRead(s, [30] * len(s), group=i))
    out.append(case("stops_single_fragment", reads))
    E, B = T.ELEM, T.BEYOND
    out.append(case("same_fragment_once", [
        ClipRead(E + "A" * 18 + B, [30] * 68, group="f1"), ClipRead(E + "A" * 19 + B, [30] * 69, group="f1"),
        ClipRead(E + "A" * 17, [30] * 37, group="f2"), ClipRead(E + "A" * 16, [30] * 36, group="f3")]))
    hi, lo = E + "T" + "GATTACA", E + "C" + "GATTACA"
    lq = [30] * 20 + [3] + [30] * 7
    out.append(case("quality_weighting", [ClipRead(hi, [40] * len(hi), 1), ClipRead(hi, [40] * len(hi), 2),
                                          ClipRead(lo, lq, 3), ClipRead(lo, lq, 4), ClipRead(lo, lq, 5)]))
    a, b = E + "TTTTGACCA", E + "CCGAGGTCA"
    out.append(case("disagreement", [ClipRead(a, [30] * len(a), 1), ClipRead(a, [30] * len(a), 2),
                                     ClipRead(b, [30] * len(b), 3), ClipRead(b, [30] * len(b), 4)]))
    clip = E + "A" * 18 + B[:10]
    reads = [ClipRead(clip, [30] * len(clip), g) for g in ("a", "b", "c")]
    reads += [ClipRead(revcomp(E[-5:] + "A" * 18 + B), [30] * (23 + len(B)), g, anchored=False) for g in ("a", "b")]
    junk = "ACGTTGCAAGGCTTAACCGGTA"
    reads.append(ClipRead(junk, [30] * len(junk), "c", anchored=False))
    out.append(case("mates_extend", reads))
    rng = random.Random(3)
    reads = []
    for i in range(5):
        s = "T" * (15 + rng.randint(-2, 2)) + revcomp(E)
        reads.append(ClipRead(s, [30] * len(s), i))
    out.append(case("polyT_at_junction", reads))
    out.append(case("empty", []))
    out.append(case("single", [ClipRead("ACGTACGTAC", [30] * 10, 1)]))
    return out


def rand_template(rng):
    alpha = rng.choice(["ACGT", "ACGT", "AACGTT", "AAAACGT", "ACGTTTT", "AC", "AG"])
    parts = []
    if rng.random() < 0.25:
        parts.append("T" * rng.randint(4, 22))
    parts.append("".join(rng.choice(alpha) for _ in range(rng.randint(5, 60))))
    if rng.random() < 0.6:
        parts.append("A" * rng.randint(4, 25))
    if rng.random() < 0.8:
        parts.append("".join(rng.choice(alpha) for _ in range(rng.randint(0, 50))))
    s = "".join(parts)
    # sprinkle homopolymers
    s = list(s)
    for _ in range(rng.randint(0, 4)):
        k = rng.randrange(len(s))
        s[k:k + 1] = [s[k]] * rng.randint(2, 6)
    return "".join(s)


def jitter(rng, s, j):
    """vary every homopolymer run >= 4 by +-j"""
    out, k = [], 0
    while k < len(s):
        e = k
        while e < len(s) and s[e] == s[k]:
            e += 1
        n = e - k
        if n >= 4 and j:
            n = max(1, n + rng.randint(-j, j))
        out.append(s[k] * n)
        k = e
    return "".join(out)


def mutate(rng, s, sub, indel, nrate, lower):
    out = []
    for c in s:
        x = rng.random()
        if x < sub:
            out.append(rng.choice([b for b in "ACGT" if b != c]))
        elif x < sub + indel / 2:
            continue
        elif x < sub + indel:
            out.append(c)
            out.append(rng.choice("ACGT"))
        else:
            out.append(c)
        if out and rng.random() < nrate:
            out[-1] = "N"
    s = "".join(out)
    if lower and s:
        a = rng.randrange(len(s))
        b = rng.randint(a, len(s))
        s = s[:a] + s[a:b].lower() + s[b:]
    return s


def rand_case(rng, i):
    tmpl = rand_template(rng)
    alt = None
    if rng.random() < 0.35:   # second haplotype -> close/tied votes
        t = list(tmpl)
        for _ in range(rng.randint(1, 3)):
            k = rng.randrange(len(t))
            t[k] = rng.choice([b for b in "ACGT" if b != t[k]])
        if rng.random() < 0.3:
            k = rng.randrange(len(t))
            t.insert(k, rng.choice("ACGT") * rng.randint(1, 3))
        alt = "".join(t)
    big = rng.random() < 0.06
    n_reads = rng.randint(41, 90) if big else rng.randint(1, 24)
    sub = rng.choice([0, 0, 0.005, 0.01, 0.03, 0.08])
    indel = rng.choice([0, 0, 0.005, 0.02, 0.05])
    nrate = rng.choice([0, 0, 0, 0.01, 0.05])
    jit = rng.choice([0, 1, 2, 3])
    constq = rng.random() < 0.4
    q0 = rng.choice([2, 20, 30, 40])
    lower = rng.random() < 0.05
    # fragment clusters
    n_clusters = max(1, rng.randint(max(1, n_reads // 4), n_reads))
    reads = []
    for k in range(n_reads):
        src = alt if (alt is not None and rng.random() < 0.5) else tmpl
        s = mutate(rng, jitter(rng, src, jit), sub, indel, nrate, lower)
        if rng.random() < 0.7 and s:
            s = s[:rng.randint(max(1, len(s) // 3), len(s))]
        if rng.random() < 0.03:
            s = ""
        anchored = True
        r = rng.random()
        if r < 0.15 and len(tmpl) > 15:   # floating mate: a later window of the template
            a = rng.randint(0, len(tmpl) // 2)
            s = mutate(rng, jitter(rng, src[a:a + rng.randint(15, 80)], jit), sub, indel, nrate, lower)
            if rng.random() < 0.5:
                s = revcomp(s)
            anchored = False
        elif r < 0.2:                       # junk mate
            s = "".join(rng.choice("ACGT") for _ in range(rng.randint(10, 60)))
            anchored = False
        if constq:
            q = [q0] * len(s)
        else:
            q = [rng.randint(2, 41) for _ in s]
        if s and rng.random() < 0.01:       # qual shorter than seq (zip truncation)
            q = q[:rng.randint(0, len(q))]
        reads.append([s, q, rng.randrange(n_clusters), anchored])
    if rng.random() < 0.8:                  # junction.rs order: by cluster
        reads.sort(key=lambda x: x[2])
    sizes = {}
    for x in reads:
        sizes[x[2]] = sizes.get(x[2], 0) + 1
    unit = rng.random() < 0.2
    cr = [ClipRead(s, q, g, 1.0 if unit else 1.0 / sizes[g], a) for s, q, g, a in reads]
    return case(f"random_{i}", cr, min_depth=rng.choice([1, 2, 2, 3]),
                polya_min_len=rng.choice([5, 8, 8, 10]))


def main():
    n = int(sys.argv[1])
    seed = int(sys.argv[2])
    out = sys.argv[3]
    cases = deterministic_cases()
    rng = random.Random(seed)
    n_raise = 0
    for i in range(n):
        try:
            cases.append(rand_case(rng, i))
        except TypeError:   # _pick_seed with only empty-RLE anchored reads (qual == []): python
            n_raise += 1    # raises; consensus.rs panics (tests::empty_rle_seed_panics)
    print("skipped (python raises):", n_raise)
    with gzip.open(out, "wt") as fh:
        json.dump(cases, fh, separators=(",", ":"))
    from collections import Counter
    print(len(cases), "cases;", Counter(c["result"]["stop_reason"] for c in cases),
          "polyA:", sum(c["result"]["polya_base"] is not None for c in cases))


if __name__ == "__main__":
    main()
