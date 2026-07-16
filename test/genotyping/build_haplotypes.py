#!/usr/bin/env python3
"""Two-haplotype donor builder for the GENOTYPING panel (test/genotyping).

Unlike test/fullstack/build_donor.py (homozygous implants -> VAF~1), genotyping needs
*het* loci (VAF~0.5). We build TWO donors from the SAME hs1 windows and the SAME implant
plan, differing only at the implant sites:

  * donor_ins.fa  -- TP windows WITH the implanted elements (+ TSD) + FP compartment
  * donor_ref.fa  -- the SAME TP windows verbatim (pre-insertion allele) + FP compartment

run_genotyping.sh then mixes wgsim reads from the two donors per sample: a het sample
draws depth/2 from each (VAF 0.5), a wild-type sample draws only from donor_ref (VAF 0).
Both donors share the FP compartment byte-for-byte, so the hs1->hg38 assembly-discordance
false positives are identical in every sample -- the property combine_genotypes' 3-negative
discriminator is meant to reject.

Element classes: the young retrotransposons of build_donor_10k.py (L1HS/AluYa5/HERVK/
SVA_E/SVA_F) PLUS a **processed-pseudogene** class whose inserted sequence is a real mature
mRNA (test/genotyping/pseudogenes.fa: DUX4/MALAT1/HNRNPA1/CASP12/RPL21) + poly-A + TSD +
EN-motif target. Because that mRNA is identical to the parent gene's transcript, exonic
spanning reads can bwa-mismap to the parent locus in hg38 -> alt reads lost -> depressed VAF:
the genotyping failure mode no pure-retrotransposon class exercises.

Outputs: donor_ins.fa, donor_ref.fa, truth_hs1.tsv, flanks.fa (truth/flanks as build_donor_10k;
lift_truth.py maps them into hg38 for scoring).
"""
import argparse
import random
import pysam

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]

# real young retrotransposons to implant (hs1 coords, 1-based inclusive, strand)
ELEMENT_LOCI = {
    "L1HS":  ("chrX", 11289796, 11295826, "+"),
    "AluYa5": ("chrX", 2759615, 2759925, "+"),
    "HERVK": ("chr10", 5039043, 5047200, "-"),
    "SVA_E": ("chr3",  148851419, 148853373, "+"),
    "SVA_F": ("chr18",  54545509,  54547673, "+"),
}

# retrotransposon family/variant sampling weights (relative frequency)
RTE_WEIGHTS = [
    ("L1HS", "full", 3), ("L1HS", "5p_truncated", 5), ("L1HS", "5p_inversion", 3),
    ("L1HS", "3p_transduction", 2),
    ("AluYa5", "full", 8), ("AluYa5", "5p_truncated", 3),
    ("HERVK", "provirus", 1),
    ("SVA_F", "full", 2), ("SVA_F", "5p_truncated", 1), ("SVA_E", "full", 1),
]
# processed-pseudogene variants (parent-gene symbols filled from --pseudogene-fasta)
PG_VARIANTS = ["full", "5p_truncated"]

POLYA_MIN, POLYA_MAX = 15, 40
TSD_MIN, TSD_MAX = 8, 20
FLANK = 300


def fetch(fa, contig, s1, e1, strand="+"):
    seq = fa.fetch(contig, s1 - 1, e1).upper()
    return rc(seq) if strand == "-" else seq


def read_fasta(path):
    recs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs[name] = "".join(buf).upper()
    return recs


def read_bed(path):
    out = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            out.append((p[0], int(p[1]), int(p[2]), p[3] if len(p) > 3 else "."))
    return out


def build_element(fa, family, variant, rng, cache, pseudogenes):
    """Return (inserted_sequence, {tsd}). inserted_sequence is 5'->3' on the + strand."""
    if family in pseudogenes:
        # processed pseudogene: the mature mRNA (+ poly-A), optionally 5' truncated.
        elem = pseudogenes[family]
        polya = "A" * rng.randint(POLYA_MIN, POLYA_MAX)
        if variant == "5p_truncated" and len(elem) > 500:
            keep = rng.randint(400, len(elem))
            body = elem[-keep:]
        else:
            body = elem
        return body + polya, {"tsd": rng.randint(TSD_MIN, TSD_MAX)}

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


def build_assignment_plan(n_implants, n_pseudogene, pseudogenes, rng):
    """A length-n_implants list of (family, variant): n_pseudogene retrocopies (spread
    evenly over the parent genes) + the rest from the retrotransposon weight table."""
    plan = []
    pg_names = sorted(pseudogenes)
    if pg_names:
        for k in range(n_pseudogene):
            fam = pg_names[k % len(pg_names)]
            var = "full" if k % 3 else "5p_truncated"     # ~1/3 truncated
            plan.append((fam, var))
    weighted = [(f, v) for (f, v, w) in RTE_WEIGHTS for _ in range(w)]
    while len(plan) < n_implants:
        plan.append(rng.choice(weighted))
    rng.shuffle(plan)
    return plan[:n_implants]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hs1", default="/Users/jeremy/Downloads/hs1.fa")
    p.add_argument("--tp-bed", required=True)
    p.add_argument("--fp-bed", required=True)
    p.add_argument("--pseudogene-fasta", required=True,
                   help="mature parent transcripts (record id = gene symbol)")
    p.add_argument("--out-donor-ins", required=True)
    p.add_argument("--out-donor-ref", required=True)
    p.add_argument("--out-truth", required=True)
    p.add_argument("--out-flanks", required=True)
    p.add_argument("--n-implants", type=int, default=1000)
    p.add_argument("--n-pseudogene", type=int, default=100)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--min-site-gap", type=int, default=10000,
                   help="site spacing; must exceed max element + 2*flank so flanks stay clean")
    args = p.parse_args()

    rng = random.Random(args.seed)
    fa = pysam.FastaFile(args.hs1)
    pseudogenes = read_fasta(args.pseudogene_fasta)
    cache = {}

    d_ins = open(args.out_donor_ins, "w")
    d_ref = open(args.out_donor_ref, "w")
    flanks = open(args.out_flanks, "w")
    truth = open(args.out_truth, "w")
    truth.write("id\tdonor_contig\tdonor_pos\ths1_contig\ths1_pos\tfamily\tvariant\t"
                "tsd\tstrand\telem_len\trole\n")

    tp = read_bed(args.tp_bed)
    tp_len = [e - s for _, s, e, _ in tp]
    total_len = sum(tp_len)
    plan = build_assignment_plan(args.n_implants, args.n_pseudogene, pseudogenes, rng)
    plan_i = 0

    tid = 0
    for wi, (contig, rs, re, _) in enumerate(tp):
        want = round(args.n_implants * tp_len[wi] / total_len)
        if wi == len(tp) - 1:
            want = args.n_implants - tid            # remainder into the last window
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
            print(f"WARN {dname}: only {len(chosen)} sites for {want} implants "
                  f"(widen the TP window or lower --min-site-gap)")

        records = []
        for site in chosen:
            family, variant = plan[plan_i]; plan_i += 1
            elem_seq, meta = build_element(fa, family, variant, rng, cache, pseudogenes)
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

        # donor_ins: region with elements spliced in (+ duplicated target site = TSD)
        out = []
        prev = 0
        for rec in sorted(records, key=lambda r: r[12]):
            ip = rec[12]; tsd_seq = rec[10]; ins = rec[11]
            out.append(region[prev:ip]); out.append(ins); out.append(tsd_seq)
            prev = ip
        out.append(region[prev:])
        d_ins.write(f">{dname}\n{''.join(out)}\n")
        # donor_ref: the SAME window verbatim (pre-insertion allele) -> same contig name
        d_ref.write(f">{dname}\n{region}\n")

        for rec in records:
            (rid, dnm, dpos, hc, hpos, fam, var, tsd, strand, elen,
             tsd_seq, ins, ip, lflank, rflank) = rec
            truth.write(f"{rid}\t{dnm}\t{dpos}\t{hc}\t{hpos}\t{fam}\t{var}\t{tsd}\t"
                        f"{strand}\t{elen}\tTP\n")
            flanks.write(f">{rid}_L\n{lflank}\n>{rid}_R\n{rflank}\n")
    print(f"implanted {tid} MEIs across {len(tp)} TP windows "
          f"({args.n_pseudogene} processed pseudogenes)")

    # FP compartment: emit hs1 satellite/telomere windows verbatim into BOTH donors,
    # identically, so the assembly-discordance FPs are the same in every sample.
    fp = read_bed(args.fp_bed)
    fp_bp = 0
    for contig, s, e, lab in fp:
        seq = fetch(fa, contig, s + 1, e)
        rec = f">{contig}_{s}_{e}_FP_{lab}\n{seq}\n"
        d_ins.write(rec); d_ref.write(rec)
        fp_bp += len(seq)
    print(f"FP compartment: {len(fp)} windows, {fp_bp/1e6:.1f} Mb (identical in both donors)")

    d_ins.close(); d_ref.close(); flanks.close(); truth.close()


if __name__ == "__main__":
    main()
