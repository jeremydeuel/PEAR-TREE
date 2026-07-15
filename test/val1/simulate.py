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


def make_read(hdr, tid, name, seq, start, cigar, mapq, flag=0,
              next_tid=-1, next_start=-1, tlen=0):
    a = pysam.AlignedSegment(hdr)
    a.query_name = name
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.flag = flag
    a.reference_id = tid
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigarstring = cigar
    # mate/discordant fields (Feature A/B): PNEXT/RNEXT + template length. A read is
    # "discordant" when it is paired (0x1) but not proper (no 0x2); its mate is mapped
    # (no 0x8) either on another contig or far away on the same one.
    a.next_reference_id = next_tid
    a.next_reference_start = next_start
    a.template_length = tlen
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

    # SENS-8 targets: an insertion whose RIGHT junction clip is a short pure poly-A
    # tail (>= min_good_bases but <= the 12 bp floor), so baseline rejects that side
    # (no pair); SENS-8 accepts it because it is a pure poly-A terminus.
    short = "A" * args.shortpolya_len
    for _ in range(args.n_shortpolya):
        contig = contigs[idx % len(contigs)]; ti = tid[contig]
        L = slot[ti]; slot[ti] += step
        tsd = rng.randint(args.tsd_min, min(args.tsd_max, 20)); R = L + tsd
        alt = max(args.alt_min, 2)
        tag = f"SPA_{contig}_{L}"
        for k in range(alt):
            records.append(make_read(hdr, ti, f"{tag}_R{k}", anchor_r + short,
                                     R - anchor_m, f"{anchor_m}M{len(short)}S", 60, flag=0x1 | 0x40))
            records.append(make_read(hdr, ti, f"{tag}_L{k}", erv5 + anchor_l,
                                     L, f"{clip_s}S{anchor_m}M", 60, flag=0x1 | 0x40))
        coverage_reads(tag, ti, L)
        truth.append((contig, L, R, "SHORTPOLYA", tsd, alt, 0, 0.0))
        idx += 1

    # ------------------------------------------------------------------ Feature A
    # DISCORDANT targets (A1): a one-sided junction — a real LEFT soft-clip at L, no
    # RIGHT soft-clip — rescued by a cluster of reverse discordant anchor reads that
    # start at ~R (reference_start = R => RIGHT-role partner) whose mates map into an
    # element locus. `discordant` variants send mates to the RTE band (contig 3,
    # covered by the emitted rmsk track => RTE-origin true); `disc_artefact` variants
    # send mates to a non-RTE band. With `discordant_anchor` on (no gate) BOTH are
    # rescued; under `discordant_rte_only` only the RTE-origin ones survive.
    elem_tid = tid[contigs[2]]  # contig "3" holds the synthetic element / gene loci
    RTE_BAND = 15_000_000       # mates here fall inside the emitted rmsk track (young)
    NONRTE_BAND = 16_000_000    # mates here are random genome (no rmsk entry)
    rmsk_rows = []  # (contig, begin, end, div) young RepeatMasker element intervals

    def discordant_block(n, mate_band, cls, in_truth):
        nonlocal idx
        for j in range(n):
            contig = contigs[idx % len(contigs)]
            if contig == contigs[2]:
                contig = contigs[0]  # keep the insertion off the element contig
            ti = tid[contig]
            L = slot[ti]; slot[ti] += step
            tsd = rng.randint(args.tsd_min, min(args.tsd_max, 20)); R = L + tsd
            alt = max(args.alt_min, 2)
            ndisc = max(args.discordant_min_reads, 3)
            tag = f"{cls}_{contig}_{L}"
            # real LEFT breakpoint at L
            for k in range(alt):
                records.append(make_read(hdr, ti, f"{tag}_L{k}", erv5 + anchor_l,
                                         L, f"{clip_s}S{anchor_m}M", 60, flag=0x1 | 0x40))
            # reverse discordant anchors at R, mate -> element band on contig 3
            for k in range(ndisc):
                mpos = mate_band + j * 2000 + k * 10
                records.append(make_read(hdr, ti, f"{tag}_D{k}", rnd_seq(rng, ref_m),
                                         R, f"{ref_m}M", 60,
                                         flag=0x1 | 0x10 | 0x80,  # paired, reverse, read2, NOT proper
                                         next_tid=elem_tid, next_start=mpos))
            coverage_reads(tag, ti, L)
            if in_truth:
                truth.append((contig, L, R, cls.upper(), tsd, alt, 0, 0.0))
            idx += 1

    discordant_block(args.n_discordant, RTE_BAND, "discordant", in_truth=True)
    discordant_block(args.n_disc_artefact, NONRTE_BAND, "disc_artefact", in_truth=False)
    if args.n_discordant or args.n_disc_artefact:
        # one young-element interval spanning the whole RTE band the mates map into
        rmsk_rows.append((contigs[2], RTE_BAND, RTE_BAND + 1_000_000, 3.0))

    # ------------------------------------------------------------------ Feature B
    # PSEUDOGENE targets (B3): a normal insertion whose RIGHT-junction reads are read1
    # of pairs whose mates map into >=2 exons of a single synthetic gene (introns
    # skipped) — the processed-pseudogene signature. A single-exon variant must NOT be
    # flagged. Exons are emitted to the companion annotation file.
    GENE_START = 17_000_000
    EXON_LEN = 400
    EXON_GAP = 10_000  # intron between exons
    exon_rows = []  # (contig, begin, end, gene_id)
    n_gene_exons = 3
    for e in range(n_gene_exons):
        b = GENE_START + e * EXON_GAP
        exon_rows.append((contigs[2], b, b + EXON_LEN, "G1"))

    def pseudogene_block(n, n_exons_hit, cls, in_truth):
        nonlocal idx
        for _ in range(n):
            contig = contigs[idx % len(contigs)]
            if contig == contigs[2]:
                contig = contigs[0]
            ti = tid[contig]
            L = slot[ti]; slot[ti] += step
            tsd = rng.randint(args.tsd_min, min(args.tsd_max, 20)); R = L + tsd
            alt = max(args.alt_min, n_exons_hit)
            tag = f"{cls}_{contig}_{L}"
            # LEFT breakpoint (plain) at L
            for k in range(alt):
                records.append(make_read(hdr, ti, f"{tag}_L{k}", erv5 + anchor_l,
                                         L, f"{clip_s}S{anchor_m}M", 60, flag=0x1 | 0x40))
            # RIGHT breakpoint at R: forward read1 => has_mate; mate maps into exon (k mod n_exons_hit)
            for k in range(alt):
                exon = k % n_exons_hit
                mpos = GENE_START + exon * EXON_GAP + 50
                nm = f"{tag}_R{k}"
                records.append(make_read(hdr, ti, nm, anchor_r + erv5,
                                         R - anchor_m, f"{anchor_m}M{clip_s}S", 60,
                                         flag=0x1 | 0x40,  # paired, read1, forward
                                         next_tid=elem_tid, next_start=mpos))
                # the mate as a real record inside the exon (read2)
                records.append(make_read(hdr, elem_tid, nm, rnd_seq(rng, ref_m),
                                         mpos, f"{ref_m}M", 60,
                                         flag=0x1 | 0x80, next_tid=ti, next_start=R - anchor_m))
            coverage_reads(tag, ti, L)
            if in_truth:
                truth.append((contig, L, R, cls.upper(), tsd, alt, 0, 0.0))
            idx += 1

    pseudogene_block(args.n_pseudogene, n_gene_exons, "pseudogene", in_truth=True)
    pseudogene_block(args.n_single_exon, 1, "single_exon", in_truth=True)

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

    # Feature A: emit a matching RepeatMasker .out track (div-gated parser format:
    # col1=%div, col4=contig, col5=begin 1-based, col6=end) covering the RTE band the
    # discordant mates map into, so the D3 config can point `discordant_rte_track` at it.
    if args.out_rmsk and rmsk_rows:
        with open(args.out_rmsk, "w") as f:
            f.write("   SW  perc perc perc  query  begin  end (left) strand repeat class ...\n")
            for contig, begin0, end0, div in rmsk_rows:
                # 1-based inclusive begin; parser stores [begin-1, end)
                f.write(f"1000 {div} 0.0 0.0 {contig} {begin0 + 1} {end0} (0) + ELEM LINE/L1 1 100 (0) 1\n")

    # Feature B: emit the exon annotation (contig, begin, end, gene_id; 0-based
    # half-open) so the D5 config can point `exon_annotation` at it.
    if args.out_exons and exon_rows:
        with open(args.out_exons, "w") as f:
            for contig, begin, end, gene in exon_rows:
                f.write(f"{contig}\t{begin}\t{end}\t{gene}\n")

    print(f"wrote {len(records)} reads, {len(truth)} insertions to {args.out_bam}; truth -> {args.out_truth}",
          file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-bam", required=True)
    p.add_argument("--out-truth", required=True)
    p.add_argument("--out-rmsk", help="Feature A: write a RepeatMasker .out track for the RTE band")
    p.add_argument("--out-exons", help="Feature B: write the exon annotation (contig,begin,end,gene_id)")
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
    p.add_argument("--n-shortpolya", type=int, default=0, help="SENS-8: insertions with a short poly-A clip")
    p.add_argument("--shortpolya-len", type=int, default=11, help="length of the short poly-A clip (<= 12 floor)")
    p.add_argument("--n-discordant", type=int, default=0, help="A1: one-sided junctions rescued by RTE-origin discordant mates")
    p.add_argument("--n-disc-artefact", type=int, default=0, help="A2: one-sided junctions whose discordant mates are non-RTE (not in truth)")
    p.add_argument("--discordant-min-reads", type=int, default=3, help="discordant anchor reads per one-sided junction")
    p.add_argument("--n-pseudogene", type=int, default=0, help="B3: insertions whose mates span >=2 exons of one gene")
    p.add_argument("--n-single-exon", type=int, default=0, help="B3 negative control: mates hit a single exon (must not flag)")
    p.add_argument("--n-artefacts", type=int, default=0, help="pileup false-positive regions (not in truth)")
    p.add_argument("--artefact-cov-min", type=int, default=150)
    p.add_argument("--artefact-cov-max", type=int, default=300)
    p.add_argument("--artefact-mapq", type=int, default=60, help="MAPQ of artefact clipped reads (30 => low-MAPQ pileup)")
    simulate(p.parse_args())


if __name__ == "__main__":
    main()
