#!/usr/bin/env python3
"""Scaled full-stack donor builder (10k implants + centromere/telomere FP compartment).

Same design as test/fullstack/build_donor.py but:
 * implants N_IMPLANTS (~10k) real young elements across several euchromatin TP windows
   at TTAAAA (L1 EN) sites, spaced >= --min-site-gap, each with TSD + poly-A + variant;
 * emits a large FP compartment = hs1 centromere/satellite + telomere windows (from
   fp_windows.bed) verbatim, so mapping hs1->hg38 yields assembly-discordance FPs.

Outputs donor.fa, truth_hs1.tsv, flanks.fa (as build_donor.py).
"""
import argparse
import random
import pysam

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]

# real young elements to implant (hs1 coords, 1-based inclusive, strand) -- from build_donor.py
ELEMENT_LOCI = {
    "L1HS":  ("chrX", 11289796, 11295826, "+"),
    "AluYa5": ("chrX", 2759615, 2759925, "+"),
    "HERVK": ("chr10", 5039043, 5047200, "-"),
    "SVA_E": ("chr3",  148851419, 148853373, "+"),
    "SVA_F": ("chr18",  54545509,  54547673, "+"),
}

# family/variant sampling weights (relative frequency of implanted classes)
PLAN_WEIGHTS = [
    ("L1HS", "full", 3), ("L1HS", "5p_truncated", 5), ("L1HS", "5p_inversion", 3),
    ("L1HS", "3p_transduction", 2),
    ("AluYa5", "full", 8), ("AluYa5", "5p_truncated", 3),
    ("HERVK", "provirus", 1),
    ("SVA_F", "full", 2), ("SVA_F", "5p_truncated", 1), ("SVA_E", "full", 1),
]
POLYA_MIN, POLYA_MAX = 15, 40
TSD_MIN, TSD_MAX = 8, 20
FLANK = 300


def fetch(fa, contig, s1, e1, strand="+"):
    seq = fa.fetch(contig, s1 - 1, e1).upper()
    return rc(seq) if strand == "-" else seq


def build_element(fa, family, variant, rng, cache):
    if family not in cache:
        c, s, e, st = ELEMENT_LOCI[family]
        cache[family] = fetch(fa, c, s, e, st)
    elem = cache[family]
    polya = "A" * rng.randint(POLYA_MIN, POLYA_MAX)
    if family == "HERVK":
        return elem, {"tsd": rng.randint(5, 6)}
    if variant == "full":
        body = elem
    elif variant == "5p_truncated":
        keep = min(len(elem), rng.randint(1200, 2000))
        body = elem[-keep:]
    elif variant == "5p_inversion":
        inv = min(len(elem) - 100, rng.randint(800, 2500))
        body = rc(elem[:inv]) + elem[inv:]
    elif variant == "3p_transduction":
        tag = "".join(rng.choice("ACGT") for _ in range(rng.randint(150, 400)))
        body = elem + tag
    else:
        body = elem
    return body + polya, {"tsd": rng.randint(TSD_MIN, TSD_MAX)}


def find_motif_sites(seq, rng, motif="TTAAAA", margin=5000):
    sites = []
    i = seq.find(motif, margin)
    while i != -1 and i < len(seq) - margin:
        sites.append(i)
        i = seq.find(motif, i + 1)
    rng.shuffle(sites)
    return sites


def read_bed(path):
    out = []
    with open(path) as f:
        for line in f:
            p = line.split()
            out.append((p[0], int(p[1]), int(p[2]), p[3] if len(p) > 3 else "."))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hs1", default="/Users/jeremy/Downloads/hs1.fa")
    p.add_argument("--tp-bed", required=True)
    p.add_argument("--fp-bed", required=True)
    p.add_argument("--out-donor", required=True)
    p.add_argument("--out-truth", required=True)
    p.add_argument("--out-flanks", required=True)
    p.add_argument("--n-implants", type=int, default=10000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--min-site-gap", type=int, default=2500)
    args = p.parse_args()

    rng = random.Random(args.seed)
    fa = pysam.FastaFile(args.hs1)
    cache = {}
    weighted = [(f, v) for (f, v, w) in PLAN_WEIGHTS for _ in range(w)]

    donor = open(args.out_donor, "w")
    flanks = open(args.out_flanks, "w")
    truth = open(args.out_truth, "w")
    truth.write("id\tdonor_contig\tdonor_pos\ths1_contig\ths1_pos\tfamily\tvariant\t"
                "tsd\tstrand\telem_len\trole\n")

    tp = read_bed(args.tp_bed)
    # distribute implants proportionally to window length
    tp_len = [e - s for _, s, e, _ in tp]
    total_len = sum(tp_len)
    tid = 0
    for wi, (contig, rs, re, _) in enumerate(tp):
        want = round(args.n_implants * tp_len[wi] / total_len)
        if wi == len(tp) - 1:
            want = args.n_implants - tid  # remainder into last window
        region = fetch(fa, contig, rs + 1, re)
        dname = f"{contig}_{rs}_{re}_TP"
        sites = find_motif_sites(region, rng)
        chosen, last = [], -10**9
        for s in sorted(sites):
            if len(chosen) >= want:
                break
            if s - last >= args.min_site_gap:
                chosen.append(s); last = s
        if len(chosen) < want:
            print(f"WARN {dname}: only {len(chosen)} sites for {want} implants")

        records = []
        for site in chosen:
            family, variant = rng.choice(weighted)
            elem_seq, meta = build_element(fa, family, variant, rng, cache)
            tsd = meta["tsd"]
            strand = rng.choice(["+", "-"])
            ins = elem_seq if strand == "+" else rc(elem_seq)
            insert_point = site + tsd
            tsd_seq = region[site:insert_point]
            left_flank = region[max(0, insert_point - FLANK):insert_point]
            right_flank = region[insert_point:insert_point + FLANK]
            records.append((tid, dname, insert_point, contig, rs + insert_point,
                            family, variant, tsd, strand, len(elem_seq),
                            tsd_seq, ins, insert_point, left_flank, right_flank))
            tid += 1

        out = []
        prev = 0
        for rec in sorted(records, key=lambda r: r[12]):
            ip = rec[12]; tsd_seq = rec[10]; ins = rec[11]
            out.append(region[prev:ip]); out.append(ins); out.append(tsd_seq)
            prev = ip
        out.append(region[prev:])
        donor.write(f">{dname}\n{''.join(out)}\n")

        for rec in records:
            (rid, dnm, dpos, hc, hpos, fam, var, tsd, strand, elen,
             tsd_seq, ins, ip, lflank, rflank) = rec
            truth.write(f"{rid}\t{dnm}\t{dpos}\t{hc}\t{hpos}\t{fam}\t{var}\t{tsd}\t"
                        f"{strand}\t{elen}\tTP\n")
            flanks.write(f">{rid}_L\n{lflank}\n>{rid}_R\n{rflank}\n")
    donor_tp = tid
    print(f"implanted {tid} MEIs across {len(tp)} TP windows")

    # FP compartment: emit hs1 satellite/telomere windows verbatim (no implants)
    fp = read_bed(args.fp_bed)
    fp_bp = 0
    for contig, s, e, lab in fp:
        seq = fetch(fa, contig, s + 1, e)
        donor.write(f">{contig}_{s}_{e}_FP_{lab}\n{seq}\n")
        fp_bp += len(seq)
    print(f"FP compartment: {len(fp)} windows, {fp_bp/1e6:.1f} Mb")

    donor.close(); flanks.close(); truth.close()


if __name__ == "__main__":
    main()
