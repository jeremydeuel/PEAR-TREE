#!/usr/bin/env python3
# PEAR-TREE - build the SMALL stand-in RTE library used by the tools/rte tests.
#
# The real library (resources/rte_library/, built by tools/rte_library/build.py) is produced by
# a separate work package. The tests must not depend on it, so this script writes a tiny
# library with the SAME layout (see plans/tprt_hallmarks/SPEC.md) from:
#   * L1Base hsflil1_8438.fa (146 intact human L1, GRCh38, 8 kb windows incl. ~1 kb flanks)
#   * Dfam consensus FASTA (AluY*, AluSx, AluJb, SVA_E/F, L1HS_5end/3end, L1PA7_3end)
#   * hg38.2bit (3' flanks of three L1Base sources) and hs1.2bit + hs1 RepeatMasker (one
#     reference L1HS + 15 kb downstream for the novel-source test, one SVA_F 5' flank).
#
# It is committed together with its outputs so the tests run without the genomes. Rebuild:
#   python test/fixtures/rte_library/build_fixture.py \
#       --l1base /path/hsflil1_8438.fa --dfam /path/dfam_consensus.fa \
#       --hg38 /path/hg38.2bit --hs1 /path/hs1.2bit --hs1-rmsk /path/hs1.repeatMasker.out.gz
import argparse
import gzip
import os
import sys
from collections import Counter

import mappy
import py2bit

OUT = os.path.dirname(os.path.abspath(__file__))
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(_COMP)[::-1]


def read_fa(path):
    return [(n, s.upper()) for n, s, _ in mappy.fastx_read(path)]


def write_fa(path, recs, width=80):
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def trim_l1(seq, five, three):
    """Return (start, end) of the L1 body inside an L1Base window (already sense)."""
    a5 = mappy.Aligner(seq=five, preset="sr")
    a3 = mappy.Aligner(seq=three, preset="sr")
    st = en = None
    # align the window's pieces onto the dfam 5'/3' end models
    for h in mappy.Aligner(seq=seq, preset="map-ont").map(five):
        if h.strand == 1 and h.q_st < 30:
            st = h.r_st - h.q_st
            break
    for h in mappy.Aligner(seq=seq, preset="map-ont").map(three):
        if h.strand == 1 and len(three) - h.q_en < 30:
            en = h.r_en + (len(three) - h.q_en)
            break
    del a5, a3
    return st, en


def majority_consensus(ref, others):
    """Simple pileup consensus of `others` placed on `ref` (substitutions only + indels by vote)."""
    al = mappy.Aligner(seq=ref, preset="map-ont")
    votes = [Counter({ref[i]: 1}) for i in range(len(ref))]
    for s in others:
        for h in al.map(s):
            if h.strand != 1 or not h.is_primary:
                continue
            r, q = h.r_st, h.q_st
            for ln, op in h.cigar:
                if op in (0, 7, 8):
                    for k in range(ln):
                        votes[r + k][s[q + k]] += 1
                    r += ln; q += ln
                elif op == 1:
                    q += ln
                elif op == 2:
                    for k in range(ln):
                        votes[r + k]["-"] += 1
                    r += ln
            break
    return "".join(v.most_common(1)[0][0] for v in votes if v.most_common(1)[0][0] != "-")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--l1base", required=True)
    ap.add_argument("--dfam", required=True)
    ap.add_argument("--hg38", required=True)
    ap.add_argument("--hs1", required=True)
    ap.add_argument("--hs1-rmsk", required=True)
    a = ap.parse_args()

    dfam = {n.split("|")[0]: s for n, s in read_fa(a.dfam)}
    l1b = read_fa(a.l1base)
    elements = []          # (id, hg38 contig, body_start, body_end, strand, seq)
    for name, seq in l1b[:12]:
        uid, loc = name.split("|")
        contig, rest = loc.split(":")
        span, strand = rest[:-3], rest[-2]
        ws, we = (int(x) for x in span.split("-"))
        sense = rc(seq) if strand == "-" else seq
        st, en = trim_l1(sense, dfam["L1HS_5end"], dfam["L1HS_3end"])
        if st is None or en is None or not (5900 <= en - st <= 6200):
            continue
        body = sense[max(0, st):en]
        # genome coords of the body (window is [ws, we) on hg38)
        if strand == "+":
            gs, ge = ws + st, ws + en
        else:
            gs, ge = we - en, we - st
        elements.append((uid, contig, gs, ge, strand, body))
        if len(elements) == 8:
            break
    assert len(elements) >= 5, "could not trim enough L1Base elements"
    cons_l1 = majority_consensus(elements[0][5], [e[5] for e in elements[1:]])

    write_fa(os.path.join(OUT, "l1_intact.fa"), [(e[0], e[5]) for e in elements])
    with open(os.path.join(OUT, "l1_intact.tsv"), "w") as fh:
        fh.write("id\tclass\tsubfamily\thg38\tstrand\tlength\n")
        for e in elements:
            fh.write(f"{e[0]}\tL1\tL1HS\t{e[1]}:{e[2]}-{e[3]}\t{e[4]}\t{len(e[5])}\n")

    alus = [("ALU_Y", "AluY"), ("ALU_YA5", "AluYa5"), ("ALU_YB8", "AluYb8")]
    old = [("ALU_SX", "AluSx"), ("ALU_JB", "AluJb"), ("L1PA7", "L1PA7_3end")]
    svas = [("SVA_E", "SVA_E"), ("SVA_F", "SVA_F")]
    write_fa(os.path.join(OUT, "alu_y_intact.fa"), [(f"{n}_ref", dfam[d]) for n, d in alus])
    write_fa(os.path.join(OUT, "sva_intact.fa"), [(f"{n}_ref", dfam[d]) for n, d in svas])
    cons = [("L1HS", cons_l1)] + [(n, dfam[d]) for n, d in alus + svas + old]
    write_fa(os.path.join(OUT, "consensus.fa"), cons)

    L = len(cons_l1)
    with open(os.path.join(OUT, "consensus_landmarks.tsv"), "w") as fh:
        fh.write("consensus\tfeature\tstart\tend\n")
        # L1.3-like anatomy (approximate). Intervals below are 0-based half-open and are
        # WRITTEN 1-based inclusive (start + 1, end), the convention of resources/rte_library.
        for f, s, e in (("5UTR", 0, 910), ("ORF1", 910, 1927), ("ORF2", 1990, 5817), ("3UTR", 5817, L)):
            fh.write(f"L1HS\t{f}\t{s + 1}\t{min(e, L)}\n")
        for n, d in alus:
            ln = len(dfam[d])
            for f, s, e in (("left_monomer", 0, 120), ("A_linker", 120, 140), ("right_monomer", 140, ln)):
                fh.write(f"{n}\t{f}\t{s + 1}\t{e}\n")
        for n, d in svas:
            ln = len(dfam[d])
            for f, s, e in (("hexamer", 0, 60), ("alu_like", 60, 430), ("VNTR", 430, 860), ("SINE_R", 860, ln)):
                fh.write(f"{n}\t{f}\t{s + 1}\t{e}\n")

    with open(os.path.join(OUT, "active.tsv"), "w") as fh:
        fh.write("id\tclass\tconsensus\thot\n")
        for e in elements[:4]:
            fh.write(f"{e[0]}\tL1\tL1HS\t1\n")
        for n, _ in alus[1:]:
            fh.write(f"{n}_ref\tALU\t{n}\t1\n")
        for n, _ in svas:
            fh.write(f"{n}_ref\tSVA\t{n}\t1\n")

    # 3' flanks (5 kb downstream, element sense) of the first three elements, from hg38
    tb = py2bit.open(a.hg38)
    flanks = []
    with open(os.path.join(OUT, "transduction_sources.tsv"), "w") as fh:
        fh.write("id\tclass\ths1\thg38\tstrand\tstatus\tevidence\thot\tintact_id\n")
        for e in elements[:3]:
            uid, contig, gs, ge, strand, _ = e
            if strand == "+":
                fl = tb.sequence(contig, ge, ge + 5000)
            else:
                fl = rc(tb.sequence(contig, gs - 5000, gs))
            flanks.append((f"TD_{uid}", fl))
            fh.write(f"TD_{uid}\tL1\t.\t{contig}:{gs}-{ge}\t{strand}\treference\tL1Base\t1\t{uid}\n")
    write_fa(os.path.join(OUT, "flanks_3p.fa"), flanks)

    # hs1: SVA_F 5' flank + reference L1HS region for the novel-source test
    hs = py2bit.open(a.hs1)
    l1hs = sva = None
    rmsk_lines = []
    with gzip.open(a.hs1_rmsk, "rt") as fh:
        for i, line in enumerate(fh):
            f = line.split()
            if len(f) < 15 or not f[0].isdigit():
                continue
            if f[4] != "chrX":
                continue
            ln = int(f[6]) - int(f[5])
            if l1hs is None and f[9] == "L1HS" and f[8] == "+" and ln > 6000:
                l1hs = (f[4], int(f[5]) - 1, int(f[6]))
            if sva is None and f[9] == "SVA_F" and ln > 1300 and float(f[1]) < 3:
                sva = (f[4], int(f[5]) - 1, int(f[6]), f[8])
            if l1hs and sva:
                break
    c, s, e = l1hs
    rs, re_ = s - 2000, e + 16000
    write_fa(os.path.join(OUT, "novel_source_region.fa"), [(f"{c}:{rs}-{re_}", hs.sequence(c, rs, re_).upper())])
    with gzip.open(a.hs1_rmsk, "rt") as fh, open(os.path.join(OUT, "novel_source.rmsk.out"), "w") as out:
        for i, line in enumerate(fh):
            if i < 3:
                out.write(line)
                continue
            f = line.split()
            if len(f) < 15 or f[4] != c:
                continue
            if int(f[6]) >= rs and int(f[5]) <= re_:
                out.write(line)
            elif int(f[5]) > re_ + 100000:
                break
    c, s, e, st = sva
    up = rc(hs.sequence(c, e, e + 2000)) if st == "C" else hs.sequence(c, s - 2000, s)
    write_fa(os.path.join(OUT, "flanks_5p_sva.fa"), [(f"TD5_SVA_F_{c}_{s}", up.upper())])
    print(f"wrote fixture library to {OUT}: {len(elements)} L1, L1HS consensus {L} bp; "
          f"novel-source L1HS {l1hs}; SVA {sva}")


if __name__ == "__main__":
    sys.exit(main())
