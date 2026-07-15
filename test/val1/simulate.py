#!/usr/bin/env python3
"""VAL-1 synthetic spike-in simulator (Phase 3).

Emits a coordinate-sorted BAM of labelled retrotransposon insertions plus a truth
TSV, so discovery recall/precision can be measured across element classes, TSD
lengths and variant-allele fractions (VAF).

Each insertion is rendered as the read signals discovery actually keys on — no
reference genome or aligner required, exactly as the differential generators do:

  * a LEFT breakpoint  = left-soft-clipped reads  (CIGAR `<clip>S<anchor>M`),
    breakpoint = reference_start = L
  * a RIGHT breakpoint = right-soft-clipped reads (CIGAR `<anchor>M<clip>S`),
    breakpoint = reference_end   = R = L + TSD
  * output pairs the two when 2 <= R - L <= 40, emitting `contig:L-R`

VAF is modelled by the alt-read count at each junction versus the count of
fully-aligned reference reads spanning the locus.

*** SYNTHETIC-ONLY LIMITATION ***
This cannot prove specificity against real EN/ERV/twin-priming artefacts: clips
are clean, junctions exact, no chimeric/PCR/mapping artefacts. Treat green VAL-1
as necessary, not sufficient (see discovery_plan.md VAL-1). Orthogonal (IGV /
long-read / PCR) confirmation of real calls remains required.
"""
import argparse
import random
import sys

import pysam

# a fixed "LTR consensus" shared by all ERV insertions — models the biology that
# ERV clips are the element LTR consensus (so many insertions share one clip).
LTR_CONSENSUS = "TGCTAGGCAACCGTATTCAGGTACCGATTCAGGCATTGACCGATTCAGGTACCGATTCAGG"
ANCHOR_POOL = "ACGTTGCACCGATTACGGATCCGTTAAGGCACTTGACCGTATCGGATCACGTTAGCCGATA"

BASES = "ACGT"


def rnd_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def make_read(hdr, tid, name, seq, start, cigar, mapq, flag=0):
    a = pysam.AlignedSegment(hdr)
    a.query_name = name
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.flag = flag
    a.reference_id = tid
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigarstring = cigar
    return a


def simulate(args):
    rng = random.Random(args.seed)
    contigs = ["1", "2", "3"]
    hdr = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": 20_000_000} for c in contigs],
    })
    tid = {c: i for i, c in enumerate(contigs)}

    anchor_m = 90    # aligned bases per junction read (>40 so unclipped_cons passes)
    clip_s = 60      # clipped bases (>12 min_clip_len)
    ref_m = 150      # fully-aligned reference read length

    records = []
    truth = []
    # space insertions far apart so their breakpoints never cross-cluster
    step = 100_000
    slot = [10_000, 10_000, 10_000]  # next free position per contig

    plan = [("L1", args.n_l1), ("ERV", args.n_erv)]
    idx = 0
    for cls, n in plan:
        for _ in range(n):
            contig = contigs[idx % len(contigs)]
            ti = tid[contig]
            L = slot[ti]
            slot[ti] += step
            tsd = rng.randint(args.tsd_min, args.tsd_max)
            R = L + tsd
            alt = rng.randint(args.alt_min, args.alt_max)
            ref = rng.randint(args.ref_min, args.ref_max)

            # element clip sequences: ERV shares the LTR consensus; L1 is unique.
            if cls == "ERV":
                clip5 = LTR_CONSENSUS[:clip_s].ljust(clip_s, "A")
                clip3 = LTR_CONSENSUS[:clip_s].ljust(clip_s, "A")
            else:
                clip5 = rnd_seq(rng, clip_s)
                clip3 = rnd_seq(rng, clip_s - 15) + "A" * 15  # L1 3' poly-A-ish tail
            anchor_l = (ANCHOR_POOL * 2)[:anchor_m]  # shared anchor => consensus forms
            anchor_r = (ANCHOR_POOL[::-1] * 2)[:anchor_m]

            vaf = round(alt / (alt + ref), 3) if (alt + ref) else 0.0
            tag = f"{cls}_{contig}_{L}"

            # RIGHT breakpoint at R: right-soft-clipped, reference_end = R
            for k in range(alt):
                seq = anchor_r + clip5
                records.append(make_read(hdr, ti, f"{tag}_R{k}", seq,
                                         R - anchor_m, f"{anchor_m}M{clip_s}S", 60,
                                         flag=0x1 | 0x40))
            # LEFT breakpoint at L: left-soft-clipped, reference_start = L
            for k in range(alt):
                seq = clip3 + anchor_l
                records.append(make_read(hdr, ti, f"{tag}_L{k}", seq,
                                         L, f"{clip_s}S{anchor_m}M", 60,
                                         flag=0x1 | 0x40))
            # fully-aligned reference reads spanning the locus (set coverage / VAF)
            for k in range(ref):
                st = L - ref_m // 2 + (k % 7)
                records.append(make_read(hdr, ti, f"{tag}_ref{k}", rnd_seq(rng, ref_m),
                                         max(st, 0), f"{ref_m}M", 60, flag=0x1 | 0x40))

            truth.append((contig, L, R, cls, tsd, alt, ref, vaf))
            idx += 1

    erv5 = LTR_CONSENSUS[:clip_s].ljust(clip_s, "A")
    anchor_l = (ANCHOR_POOL * 2)[:anchor_m]
    anchor_r = (ANCHOR_POOL[::-1] * 2)[:anchor_m]

    def coverage_reads(tag, ti, L):
        for k in range(rng.randint(args.ref_min, args.ref_max)):
            records.append(make_read(hdr, ti, f"{tag}_ref{k}", rnd_seq(rng, ref_m),
                                     max(L - ref_m // 2 + (k % 7), 0), f"{ref_m}M", 60, flag=0x1 | 0x40))

    # SENS-1 targets: ERV insertions whose junction breakpoints wobble over a few bp,
    # so the exact modal count is 1 (missed) but the windowed count clears the floor.
    # Each read is one continuous genome+element sequence with the clip boundary (and
    # so the breakpoint) shifted by `off` — the realistic form the consensus can
    # delta-align back together.
    t_right = anchor_r + erv5   # [anchor | element], boundary slides for RIGHT clips
    t_left = erv5 + anchor_l    # [element | anchor], boundary slides for LEFT clips
    for _ in range(args.n_wobble):
        contig = contigs[idx % len(contigs)]; ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(args.tsd_min, args.tsd_max); R = L + tsd
        alt = max(args.alt_min, 3)
        tag = f"WOB_{contig}_{L}"
        for k in range(alt):
            off = k - alt // 2  # spread breakpoints across +/- a couple bp
            # RIGHT bp = reference_end = R + off: fixed start, slide the M/S boundary
            records.append(make_read(hdr, ti, f"{tag}_R{k}", t_right,
                                     R - anchor_m, f"{anchor_m + off}M{clip_s - off}S", 60, flag=0x1 | 0x40))
            # LEFT bp = reference_start = L + off: slide start and the S/M boundary
            records.append(make_read(hdr, ti, f"{tag}_L{k}", t_left,
                                     L + off, f"{clip_s + off}S{anchor_m - off}M", 60, flag=0x1 | 0x40))
        coverage_reads(tag, ti, L)
        truth.append((contig, L, R, "WOBBLE", tsd, alt, 0, 0.0))
        idx += 1

    # SENS-2 targets: ERV insertions with one high-MAPQ + one low-MAPQ read per
    # junction. At the default floor only the high-MAPQ read survives (n=1, missed);
    # a lowered min_mapq recovers n=2, and the low-MAPQ fraction (0.5) stays under
    # the guard so it is not dropped.
    for _ in range(args.n_lowmapq):
        contig = contigs[idx % len(contigs)]; ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(args.tsd_min, args.tsd_max); R = L + tsd
        tag = f"LMQ_{contig}_{L}"
        for mq in (60, 30):
            records.append(make_read(hdr, ti, f"{tag}_R{mq}", anchor_r + erv5,
                                     R - anchor_m, f"{anchor_m}M{clip_s}S", mq, flag=0x1 | 0x40))
            records.append(make_read(hdr, ti, f"{tag}_L{mq}", erv5 + anchor_l,
                                     L, f"{clip_s}S{anchor_m}M", mq, flag=0x1 | 0x40))
        coverage_reads(tag, ti, L)
        truth.append((contig, L, R, "LOWMAPQ", tsd, 2, 0, 0.0))
        idx += 1

    # pileup artefacts: high local coverage + a spurious clipped pair, far from real
    # insertions and NOT in truth. These are the false positives SPEC-3 (coverage
    # mask) and SPEC-4 (adaptive evidence floor) must remove.
    for a in range(args.n_artefacts):
        contig = contigs[a % len(contigs)]
        ti = tid[contig]
        A = 5_000_000 + a * 100_000
        pileup = rng.randint(args.artefact_cov_min, args.artefact_cov_max)
        for k in range(pileup):
            st = A + (k % 400)  # starts within one 500 bp bin, aligned with the breakpoints below
            records.append(make_read(hdr, ti, f"art_{contig}_{A}_p{k}", rnd_seq(rng, ref_m),
                                     st, f"{ref_m}M", 60, flag=0x1 | 0x40))
        tsd = rng.randint(args.tsd_min, args.tsd_max)
        # each junction's reads must share seq so a (false) consensus forms
        seq_r = rnd_seq(rng, anchor_m) + rnd_seq(rng, clip_s)
        seq_l = rnd_seq(rng, clip_s) + rnd_seq(rng, anchor_m)
        for k in range(2):  # spurious RIGHT breakpoint at A+tsd
            records.append(make_read(hdr, ti, f"art_{contig}_{A}_R{k}", seq_r,
                                     (A + tsd) - anchor_m, f"{anchor_m}M{clip_s}S", args.artefact_mapq, flag=0x1 | 0x40))
        for k in range(2):  # spurious LEFT breakpoint at A
            records.append(make_read(hdr, ti, f"art_{contig}_{A}_L{k}", seq_l,
                                     A, f"{clip_s}S{anchor_m}M", args.artefact_mapq, flag=0x1 | 0x40))

    records.sort(key=lambda a: (a.reference_id, a.reference_start))
    with pysam.AlignmentFile(args.out_bam, "wb", header=hdr) as out:
        for a in records:
            out.write(a)
    pysam.index(args.out_bam)

    with open(args.out_truth, "w") as f:
        f.write("contig\tleft\tright\tclass\ttsd\talt_reads\tref_reads\tvaf\n")
        for row in truth:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"wrote {len(records)} reads, {len(truth)} insertions to {args.out_bam}; truth -> {args.out_truth}",
          file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-bam", required=True)
    p.add_argument("--out-truth", required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--n-l1", type=int, default=30)
    p.add_argument("--n-erv", type=int, default=20)
    p.add_argument("--tsd-min", type=int, default=2)
    p.add_argument("--tsd-max", type=int, default=40)
    p.add_argument("--alt-min", type=int, default=2, help="min alt (junction) reads per side => low VAF")
    p.add_argument("--alt-max", type=int, default=8)
    p.add_argument("--ref-min", type=int, default=5)
    p.add_argument("--ref-max", type=int, default=40)
    p.add_argument("--n-wobble", type=int, default=0, help="SENS-1: insertions with wobbled breakpoints")
    p.add_argument("--n-lowmapq", type=int, default=0, help="SENS-2: insertions with mixed-MAPQ junction reads")
    p.add_argument("--n-artefacts", type=int, default=0, help="pileup false-positive regions (not in truth)")
    p.add_argument("--artefact-cov-min", type=int, default=150)
    p.add_argument("--artefact-cov-max", type=int, default=300)
    p.add_argument("--artefact-mapq", type=int, default=60, help="MAPQ of artefact clipped reads (30 => low-MAPQ pileup)")
    simulate(p.parse_args())


if __name__ == "__main__":
    main()
