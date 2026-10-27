#!/usr/bin/env python3
"""Synthetic coordinate-sorted BAM that stresses the `disc_agree_second_fragment` PENDING
path of discovery (analysis/discovery_mem/REPORT.md). No reference file is needed (BAM
input only): each contig gets a private random reference string from which every aligned
read part is taken, so clusters have a real unclipped consensus.

Per contig it writes
  * proper-pair background at --cov x (151 bp, MAPQ 60), 3 % of which carry a 1..11 bp
    soft clip (the junk `temporary_breakpoints` a real BAM has);
  * one SINGLE-MOLECULE junction read every --junction-every bp, alternating LEFT/RIGHT
    and alternating `has_mate` orientation, with a 40..70 bp random (complex) clip ->
    `join()` holds it PENDING;
  * for --anchor-frac of them, 1..3 discordant anchors in the 500 bp flank pointing at
    the junction (mate on the other contig, MAPQ 0); --agree-frac of those mates carry
    the insert, so the junction is promoted in the mate pass, the rest are rejected;
  * one real two-sided insertion (3 fragments per side, 15 bp TSD) every 20 kb so the
    .txt.gz / evidence sidecar are non-empty.

Usage: gen_pending_stress.py --out stress.bam [--contigs 2 --len 8000000 --cov 25
       --junction-every 150 --anchor-frac 0.9 --agree-frac 0.15 --seed 1]
"""
import argparse, random, os, sys
import pysam

RL = 151
QUAL = [40] * RL


def dna(rng, n):
    return "".join(rng.choices("ACGT", k=n))


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class W:
    def __init__(self, path, header):
        self.f = pysam.AlignmentFile(path, "wb", header=header)
        self.h = header
        self.n = 0

    def rec(self, name, flag, tid, pos, cigar, seq, mapq=60, mtid=-1, mpos=-1, tlen=0):
        a = pysam.AlignedSegment(self.h)
        a.query_name = name
        a.flag = flag
        a.reference_id = tid
        a.reference_start = pos
        a.mapping_quality = mapq
        a.cigartuples = cigar
        a.query_sequence = seq
        a.query_qualities = QUAL[: len(seq)]
        a.next_reference_id = mtid
        a.next_reference_start = mpos
        a.template_length = tlen
        self.f.write(a)
        self.n += 1

    def pair(self, name, tid, p1, p2, cig1, cig2, s1, s2, r1_fwd=True, q=60):
        """proper pair: forward read at p1, reverse read at p2 (p1 < p2), FR orientation;
        r1_fwd says which of the two is read1."""
        end2 = p2 + sum(l for op, l in cig2 if op in (0, 2, 3))
        tlen = end2 - p1
        if r1_fwd:
            f1, f2 = 0x1 | 0x2 | 0x20 | 0x40, 0x1 | 0x2 | 0x10 | 0x80
        else:
            f1, f2 = 0x1 | 0x2 | 0x20 | 0x80, 0x1 | 0x2 | 0x10 | 0x40
        self.rec(name, f1, tid, p1, cig1, s1, q, tid, p2, tlen)
        self.rec(name, f2, tid, p2, cig2, s2, q, tid, p1, -tlen)

    def close(self):
        self.f.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--contigs", type=int, default=2)
    ap.add_argument("--len", type=int, default=8_000_000)
    ap.add_argument("--cov", type=float, default=25.0)
    ap.add_argument("--junction-every", type=int, default=150)
    ap.add_argument("--anchor-frac", type=float, default=0.9)
    ap.add_argument("--agree-frac", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    names = [str(i + 1) for i in range(a.contigs)]
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": n, "LN": a.len} for n in names]}
    )
    tmp = a.out + ".unsorted.bam"
    w = W(tmp, header)
    k = 0
    for tid in range(a.contigs):
        other = (tid + 1) % a.contigs
        ref = dna(rng, a.len)
        # background proper pairs
        n_pairs = int(a.len * a.cov / (2 * RL))
        for _ in range(n_pairs):
            ins = max(2 * RL + 10, int(rng.gauss(400, 50)))
            p1 = rng.randrange(0, a.len - ins)
            p2 = p1 + ins - RL
            c1 = [(0, RL)]
            c2 = [(0, RL)]
            s1 = ref[p1 : p1 + RL]
            s2 = ref[p2 : p2 + RL]
            if rng.random() < 0.03:  # junk short clip
                s = rng.randrange(1, 12)
                if rng.random() < 0.5:
                    c1 = [(4, s), (0, RL - s)]
                    s1 = dna(rng, s) + ref[p1 : p1 + RL - s]
                else:
                    c2 = [(0, RL - s), (4, s)]
                    s2 = ref[p2 : p2 + RL - s] + dna(rng, s)
            k += 1
            w.pair(f"bg{k}", tid, p1, p2, c1, c2, s1, s2, rng.random() < 0.5)

        def left_read(b, clen, insert):
            """LEFT junction at b: clip to the left, aligned part starts at b."""
            return [(4, clen), (0, RL - clen)], insert + ref[b : b + RL - clen]

        def right_read(b, clen, insert):
            """RIGHT junction at b: aligned part ends at b, clip to the right."""
            return [(0, RL - clen), (4, clen)], ref[b - (RL - clen) : b] + insert

        def bg(p):
            return ref[p : p + RL]

        # single-molecule pending junctions with discordant anchors
        b = 2000
        j = 0
        while b < a.len - 3000:
            j += 1
            left = j % 2 == 0
            hm = j % 4 < 2  # has_mate orientation (LEFT reverse / RIGHT forward)
            clen = rng.randrange(40, 71)
            insert = dna(rng, clen)
            k += 1
            name = f"jx{k}"
            if left:
                cig, seq = left_read(b, clen, insert)
                if hm:
                    p1 = b - rng.randrange(150, 260) - RL
                    w.pair(name, tid, p1, b, [(0, RL)], cig, bg(p1), seq, True)
                else:
                    p2 = b + rng.randrange(150, 260)
                    w.pair(name, tid, b, p2, cig, [(0, RL)], seq, bg(p2), True)
            else:
                cig, seq = right_read(b, clen, insert)
                p0 = b - (RL - clen)
                if hm:
                    p2 = b + rng.randrange(150, 260)
                    w.pair(name, tid, p0, p2, cig, [(0, RL)], seq, bg(p2), True)
                else:
                    p1 = p0 - rng.randrange(150, 260)
                    w.pair(name, tid, p1, p0, [(0, RL)], cig, bg(p1), seq, True)
            if rng.random() < a.anchor_frac:
                agree = rng.random() < a.agree_frac
                for _ in range(rng.randrange(1, 4)):
                    k += 1
                    an = f"an{k}"
                    mp = rng.randrange(1000, a.len - 1000)
                    if left:
                        st = b + rng.randrange(0, 300)  # reverse anchor at/after b
                        aflag = 0x1 | 0x10 | 0x40
                        mflag = 0x1 | 0x20 | 0x80
                    else:
                        st = b - RL - rng.randrange(0, 300)  # forward anchor ending before b
                        aflag = 0x1 | 0x40
                        mflag = 0x1 | 0x80
                    w.rec(an, aflag, tid, st, [(0, RL)], bg(st), 60, other, mp, 0)
                    if agree:
                        pad = RL - 2 * clen
                        mseq = insert + dna(rng, max(pad, 0)) + rc(insert)
                        mseq = mseq[:RL] if len(mseq) >= RL else mseq + dna(rng, RL - len(mseq))
                    else:
                        mseq = dna(rng, RL)
                    w.rec(an, mflag, other, mp, [(0, RL)], mseq, 0, tid, st, 0)
            b += a.junction_every
        # real two-sided insertions (3 fragments per side, 15 bp TSD: the LEFT junction
        # sits at the TSD's right end, the RIGHT junction 15 bp upstream of it)
        for b in range(10_000, a.len - 10_000, 20_000):
            ins_l = dna(rng, 90)
            ins_r = dna(rng, 90)
            e = b - 15
            for i, clen in enumerate((40, 50, 60)):
                k += 1
                cig, seq = left_read(b, clen, ins_l[-clen:])
                if i == 0:  # one has_mate (reverse) read per side
                    p1 = b - 400
                    w.pair(f"tl{k}", tid, p1, b, [(0, RL)], cig, bg(p1), seq, True)
                else:
                    p2 = b + 160 + 20 * i
                    w.pair(f"tl{k}", tid, b, p2, cig, [(0, RL)], seq, bg(p2), True)
                k += 1
                cig, seq = right_read(e, clen, ins_r[:clen])
                p0 = e - (RL - clen)
                if i == 0:
                    p2 = e + 300
                    w.pair(f"tr{k}", tid, p0, p2, cig, [(0, RL)], seq, bg(p2), True)
                else:
                    p1 = p0 - 170 - 20 * i
                    w.pair(f"tr{k}", tid, p1, p0, [(0, RL)], cig, bg(p1), seq, True)
    w.close()
    pysam.sort("-o", a.out, tmp)
    os.unlink(tmp)
    pysam.index(a.out)
    print(f"wrote {w.n} records -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
