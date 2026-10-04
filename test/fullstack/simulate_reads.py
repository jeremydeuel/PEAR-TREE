#!/usr/bin/env python3
"""Paired-end reads for one sample of a donor_types.py donor (replaces wgsim there).

Why not wgsim: the TPRT-hallmark caller needs (a) per-read poly-A / homopolymer length
jitter (Illumina SBS slippage) and post-homopolymer phasing loss, (b) PCR duplicates that are
NOT flagged and carry independent sequencing noise, (c) exact per-junction fragment counts for
the >=2-independent-fragments rule. All come from test/simlib/reads.py.

Inputs: the donor dir (haps.tsv, S<i>.hap<k>.fa, ref.hap.fa, S<i>.molecules.fa,
junctions.tsv, slippage.tsv). Outputs: <prefix>_R1.fq, <prefix>_R2.fq and
<prefix>.support.tsv (event_id, side, fragments, reads incl. duplicates).
"""
import argparse
import os
import random
import sys
from bisect import bisect_left, bisect_right

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from simlib.reads import FragmentSampler, Hap, entries_span  # noqa: E402
from simlib.seqs import read_fasta, revcomp  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--donor-dir", required=True)
    p.add_argument("--sample", type=int, required=True, help="1-based sample index")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--depth", type=float, default=15.0)
    p.add_argument("--read-len", type=int, default=151)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--polya-jitter", type=float, default=1.0)
    p.add_argument("--phasing-burst", type=float, default=1.0)
    p.add_argument("--pcr-dup-unflagged-frac", type=float, default=0.05)
    p.add_argument("--error-rate", type=float, default=0.002)
    a = p.parse_args()

    rng = random.Random(a.seed * 1000 + a.sample)
    d = a.donor_dir
    S = a.sample
    sampler = FragmentSampler(rng, read_len=a.read_len, jitter_scale=a.polya_jitter,
                              error_rate=a.error_rate, pcr_dup_frac=a.pcr_dup_unflagged_frac,
                              burst_scale=a.phasing_burst)
    haps = []
    with open(os.path.join(d, "haps.tsv")) as f:
        next(f)
        for line in f:
            s, h, w = line.rstrip("\n").split("\t")
            if int(s) == S:
                haps.append((h, float(w)))
    junc = {}      # (hap, record) -> sorted [(pos, event_id, side)]
    with open(os.path.join(d, "junctions.tsv")) as f:
        next(f)
        for line in f:
            h, eid, side, rec, pos = line.rstrip("\n").split("\t")
            junc.setdefault((h, rec), []).append((int(pos), eid, side))
    for k in junc:
        junc[k].sort()
    slip = {}
    with open(os.path.join(d, "slippage.tsv")) as f:
        next(f)
        for line in f:
            h, rec, s0, s1 = line.rstrip("\n").split("\t")
            if h.startswith(f"S{S}:"):
                h = h.split(":", 1)[1]
            elif ":" in h:
                continue
            slip.setdefault((h, rec), []).append((int(s0), int(s1)))

    support = {}
    o1 = open(a.out_prefix + "_R1.fq", "w")
    o2 = open(a.out_prefix + "_R2.fq", "w")
    n_pairs = 0

    def emit(name, pair, r1_fwd):
        nonlocal n_pairs
        (ef, sf, qf, _), (er, sr, qr, _) = pair
        fwd = (sf, qf)
        rev = (revcomp(sr), qr[::-1])          # as sequenced
        r1, r2 = (fwd, rev) if r1_fwd else (rev, fwd)
        o1.write(f"@{name}/1\n{r1[0]}\n+\n{r1[1]}\n")
        o2.write(f"@{name}/2\n{r2[0]}\n+\n{r2[1]}\n")
        n_pairs += 1

    for hfile, w in haps:
        recs = read_fasta(os.path.join(d, hfile))
        hname = hfile.replace(".fa", "")
        for rec, seq in recs.items():
            tracts = slip.get((hfile, rec), [])
            boost = {t0: 3.0 for t0, _ in tracts}
            hap = Hap(seq, None, hname, boost)
            js = junc.get((hfile, rec), [])
            jpos = [x[0] for x in js]
            n = sampler.n_fragments(len(seq), a.depth, w)
            for i in range(n):
                fl = sampler.flen()
                st = rng.randint(0, max(0, len(seq) - fl))
                en = min(len(seq), st + fl)
                in_tract = any(st < t1 and en > t0 for t0, t1 in tracts)
                if in_tract:
                    old = (sampler.slip_min_delta, sampler.burst_scale)
                    sampler.slip_min_delta, sampler.burst_scale = 3, max(2.5, 2.5 * sampler.burst_scale)
                ncopy = sampler.copies()
                r1_fwd = rng.random() < 0.5
                lo, hi = bisect_left(jpos, st + 10), bisect_right(jpos, en - 10)
                for c in range(ncopy + 1):
                    pair = sampler.pair(hap, st, en)
                    emit(f"S{S}:{hname}:{rec}:{i}" + (f"_d{c}" if c else ""), pair, r1_fwd)
                    for (pos, eid, side) in js[lo:hi]:
                        if any(x[0] <= pos - 10 and x[1] >= pos + 10
                               for x in (entries_span(r[0]) for r in pair)):
                            k = (eid, side)
                            fr, rd = support.get(k, (0, 0))
                            support[k] = (fr + (c == 0), rd + 1)
                if in_tract:
                    sampler.slip_min_delta, sampler.burst_scale = old
    # single-molecule library artefacts (+ unflagged PCR copies)
    mpath = os.path.join(d, f"S{S}.molecules.fa")
    if os.path.exists(mpath):
        for name, mol in read_fasta(mpath).items():
            copies = int(name.split("copies=")[1]) if "copies=" in name else 0
            hap = Hap(mol, None, "mol")
            r1_fwd = rng.random() < 0.5
            for c in range(copies + 1):
                emit(f"S{S}:mol:{name.split('|')[0]}" + (f"_d{c}" if c else ""),
                     sampler.pair(hap, 0, len(mol)), r1_fwd)
    o1.close(); o2.close()
    with open(a.out_prefix + ".support.tsv", "w") as f:
        f.write("event_id\tside\tfragments\treads\n")
        for (eid, side), (fr, rd) in sorted(support.items(), key=lambda x: (int(x[0][0]), x[0][1])):
            f.write(f"{eid}\t{side}\t{fr}\t{rd}\n")
    print(f"sample {S}: {n_pairs} pairs -> {a.out_prefix}_R1/_R2.fq", file=sys.stderr)


if __name__ == "__main__":
    main()
