#!/usr/bin/env python3
"""FP-stress donor builder: 10k true insertions (all classes) + a large non-MEI decoy
population tuned to make raw *discovery* throw >=10 000 emergent false positives.

Design (window-splice, same proven geometry as build_donor_10k.py so coverage stays
uniform and `coverage_mask` behaves as it would on real WGS -- a *sparse per-cassette*
donor breaks the median-of-populated-bins estimate and masks the real loci):

  * TRUE insertions  -- canonical TPRT MEIs (element + poly-A + TSD) spliced at TTAAAA
    (L1 endonuclease) sites into euchromatin TP windows. Spans ALL classes: young
    (L1HS/AluYa5/HERVK/SVA_E/SVA_F), older (L1PA2/AluSx/AluJb) and a processed-pseudogene
    class (spliced mRNAs). role=TP.

  * DECOY false positives -- non-MEI structural differences spliced into a separate set of
    euchromatin FP windows (so the flanks anchor uniquely and the junction is a clean,
    unmasked clip cluster). The inserted material is real, genome-wide hs1 repeat sequence
    so the FP population *is* "all low-complexity + young/old retrotransposon sites", but
    every decoy deliberately OMITS the poly-A + TSD hallmarks, i.e. it is a copy-number /
    non-canonical difference a good MEI caller must reject:
      - lcr       : extra tandem copies of a real Low_complexity/Simple_repeat/Satellite
                    tract (a length polymorphism, not a retrotransposition)
      - young_rte : a divergent fragment of a real young element (AluY/L1P/SVA/ERVK)
      - old_rte   : a divergent fragment of a real old element (AluS/AluJ/MIR/L2/CR1/ERVL)
    role=FP (bookkeeping; the scorer counts any call away from a true insertion as an
    emergent FP regardless of whether it lands on a labelled decoy).

Outputs donor.fa, truth_hs1.tsv (id..role, incl. a `class` column), flanks.fa.
"""
import argparse
import random
import pysam

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]

# ---- TRUE-insertion element sources (hs1 coords, 1-based inclusive, strand) --------------
ELEMENT_LOCI = {
    "L1HS":   ("chrX", 11289796, 11295826, "+"),
    "AluYa5": ("chrX",  2759615,  2759925, "+"),
    "HERVK":  ("chr10", 5039043,  5047200, "-"),
    "SVA_E":  ("chr3", 148851419, 148853373, "+"),
    "SVA_F":  ("chr18", 54545509, 54547673, "+"),
    "L1PA2":  ("chrX",  1500493,  1506508, "-"),   # older LINE-1
    "AluSx":  ("chrX",   108031,   108331, "-"),   # older Alu (AluS)
    "AluJb":  ("chrX",   129238,   129536, "-"),   # oldest Alu (AluJ)
}
YOUNG_FAMS = {"L1HS", "AluYa5", "HERVK", "SVA_E", "SVA_F"}

PLAN_WEIGHTS = [
    ("L1HS", "full", 3), ("L1HS", "5p_truncated", 5), ("L1HS", "5p_inversion", 3),
    ("L1HS", "3p_transduction", 2),
    ("AluYa5", "full", 8), ("AluYa5", "5p_truncated", 3),
    ("HERVK", "provirus", 1),
    ("SVA_F", "full", 2), ("SVA_F", "5p_truncated", 1), ("SVA_E", "full", 1),
    ("L1PA2", "full", 1), ("L1PA2", "5p_truncated", 2),
    ("AluSx", "full", 2), ("AluJb", "full", 2),
    ("PSEUDOGENE", "spliced", 3),
]
POLYA_MIN, POLYA_MAX = 15, 40
TSD_MIN, TSD_MAX = 8, 20
FLANK = 300           # hs1 flank emitted to flanks.fa for the hg38 lift
POOL_PER_BED = 4000   # decoy source sequences sampled per repeat category


def fetch(fa, contig, s1, e1, strand="+"):
    seq = fa.fetch(contig, s1 - 1, e1).upper()
    return rc(seq) if strand == "-" else seq


# ---- TRUE-insertion element construction -------------------------------------------------
def build_true_element(fa, family, variant, rng, cache, pseudo):
    polya = "A" * rng.randint(POLYA_MIN, POLYA_MAX)
    if family == "PSEUDOGENE":
        name, seq = rng.choice(pseudo)
        if rng.random() < 0.4 and len(seq) > 800:
            seq = seq[-rng.randint(500, len(seq)):]
        return seq + polya, {"tsd": rng.randint(TSD_MIN, TSD_MAX), "sub": name}
    if family not in cache:
        c, s, e, st = ELEMENT_LOCI[family]
        cache[family] = fetch(fa, c, s, e, st)
    elem = cache[family]
    if family == "HERVK":
        return elem, {"tsd": rng.randint(5, 6), "sub": family}
    if variant == "full":
        body = elem
    elif variant == "5p_truncated":
        keep = min(len(elem), rng.randint(1200, 2000) if len(elem) > 2000 else max(150, len(elem) // 2))
        body = elem[-keep:]
    elif variant == "5p_inversion":
        inv = min(len(elem) - 100, rng.randint(800, 2500))
        body = rc(elem[:inv]) + elem[inv:]
    elif variant == "3p_transduction":
        tag = "".join(rng.choice("ACGT") for _ in range(rng.randint(150, 400)))
        body = elem + tag
    else:
        body = elem
    return body + polya, {"tsd": rng.randint(TSD_MIN, TSD_MAX), "sub": family}


# ---- non-MEI decoy INSERT construction (no poly-A, no TSD) --------------------------------
def build_lcr_insert(pool, rng):
    tract = rng.choice(pool)
    want = rng.randint(100, 400)
    unit = tract if len(tract) <= want else tract[:want]
    return (unit * (want // max(1, len(unit)) + 1))[:want]


def build_frag_insert(pool, rng):
    src = rng.choice(pool)
    L = rng.randint(150, 600)
    frag = src if len(src) <= L else src[rng.randint(0, len(src) - L):][:L]
    return rc(frag) if rng.random() < 0.5 else frag


# ---- helpers -----------------------------------------------------------------------------
def find_motif_sites(seq, rng, motif="TTAAAA", margin=5000):
    sites, i = [], seq.find(motif, margin)
    while i != -1 and i < len(seq) - margin:
        sites.append(i); i = seq.find(motif, i + 1)
    rng.shuffle(sites)
    return sites


def spaced(sorted_sites, want, gap):
    out, last = [], -10**9
    for s in sorted_sites:
        if len(out) >= want:
            break
        if s - last >= gap:
            out.append(s); last = s
    return out


def read_bed(path):
    out = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 3:
                out.append((p[0], int(p[1]), int(p[2]), p[3] if len(p) > 3 else "."))
    return out


def load_pool(fa, bed_path, rng, n, min_len=40, max_fetch=4000):
    """Sample up to n repeat intervals from a RepeatMasker BED and fetch their hs1 seq."""
    rows = read_bed(bed_path)
    rng.shuffle(rows)
    pool = []
    for c, s, e, _ in rows:
        if len(pool) >= n:
            break
        if e - s < min_len:
            continue
        seq = fetch(fa, c, s + 1, min(e, s + max_fetch))
        if seq and set(seq) <= set("ACGTN") and seq.count("N") < len(seq) // 10:
            pool.append(seq.replace("N", "A"))
    return pool


def load_pseudogenes(path):
    recs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    recs.append((name, "".join(buf)))
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip().upper())
    if name:
        recs.append((name, "".join(buf)))
    return [(n, s) for n, s in recs if s and set(s) <= set("ACGT")]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hs1", default="/Users/jeremy/Downloads/hs1.fa")
    p.add_argument("--tp-bed", required=True, help="euchromatin windows for TRUE implants")
    p.add_argument("--fp-bed", required=True, help="euchromatin windows for DECOY implants")
    p.add_argument("--lcr-bed", required=True)
    p.add_argument("--young-bed", required=True)
    p.add_argument("--old-bed", required=True)
    p.add_argument("--pseudogenes", required=True)
    p.add_argument("--out-donor", required=True)
    p.add_argument("--out-truth", required=True)
    p.add_argument("--out-flanks", required=True)
    p.add_argument("--n-implants", type=int, default=10000)
    p.add_argument("--n-decoys", type=int, default=16000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--tp-gap", type=int, default=2500)
    p.add_argument("--fp-gap", type=int, default=1500)
    args = p.parse_args()

    rng = random.Random(args.seed)
    fa = pysam.FastaFile(args.hs1)
    cache = {}
    pseudo = load_pseudogenes(args.pseudogenes)
    weighted = [(f, v) for (f, v, w) in PLAN_WEIGHTS for _ in range(w)]
    print("loading decoy source pools (genome-wide LCR / young / old RTE) ...")
    lcr_pool = load_pool(fa, args.lcr_bed, rng, POOL_PER_BED, min_len=40, max_fetch=400)
    young_pool = load_pool(fa, args.young_bed, rng, POOL_PER_BED, min_len=200)
    old_pool = load_pool(fa, args.old_bed, rng, POOL_PER_BED, min_len=200)
    print(f"  pools: lcr={len(lcr_pool)} young={len(young_pool)} old={len(old_pool)}")

    donor = open(args.out_donor, "w")
    flanks = open(args.out_flanks, "w")
    truth = open(args.out_truth, "w")
    truth.write("id\tdonor_contig\tdonor_pos\ths1_contig\ths1_pos\tfamily\tvariant\t"
                "tsd\tstrand\telem_len\tclass\trole\n")
    tid = 0

    def emit_window(contig, rs, re, dtag, records):
        """records: list of dicts with keys ip(insert point in region), ins, tsd_seq,
        family, variant, tsd, strand, elem_len, class, role. Splices them into the region."""
        nonlocal tid
        region = fetch(fa, contig, rs + 1, re)
        dname = f"{contig}_{rs}_{re}_{dtag}"
        out, prev = [], 0
        for r in sorted(records, key=lambda r: r["ip"]):
            ip = r["ip"]
            out.append(region[prev:ip]); out.append(r["ins"]); out.append(r["tsd_seq"])
            prev = ip
        out.append(region[prev:])
        donor.write(f">{dname}\n{''.join(out)}\n")
        for r in records:
            ip = r["ip"]
            truth.write(f"{tid}\t{dname}\t{ip}\t{contig}\t{rs + ip}\t{r['family']}\t"
                        f"{r['variant']}\t{r['tsd']}\t{r['strand']}\t{r['elem_len']}\t"
                        f"{r['class']}\t{r['role']}\n")
            lflank = region[max(0, ip - FLANK):ip]
            rflank = region[ip:ip + FLANK]
            flanks.write(f">{tid}_L\n{lflank}\n>{tid}_R\n{rflank}\n")
            tid += 1

    # ---- TRUE implants into TP windows ----------------------------------------------------
    tp = read_bed(args.tp_bed)
    tp_len = [e - s for _, s, e, _ in tp]
    tot = sum(tp_len) or 1
    n_tp = 0
    for wi, (contig, rs, re, _) in enumerate(tp):
        want = (args.n_implants - n_tp) if wi == len(tp) - 1 else round(args.n_implants * tp_len[wi] / tot)
        region = fetch(fa, contig, rs + 1, re)
        sites = spaced(sorted(find_motif_sites(region, rng)), want, args.tp_gap)
        if len(sites) < want:
            # TTAAAA motif supply exhausted -> top up with generic positions (junctions only
            # need to be distinct, ~fill_gap bp, not EN-motif spaced) so the window meets its
            # TRUE-implant quota. The EN motif is realism, not required for the test.
            fill_gap = 800
            occ = sorted(sites)
            cand = list(range(5000, len(region) - 5000, fill_gap))
            rng.shuffle(cand)
            import bisect
            for pos in cand:
                if len(sites) >= want:
                    break
                j = bisect.bisect_left(occ, pos)
                near = min([abs(pos - occ[k]) for k in (j - 1, j) if 0 <= k < len(occ)] or [fill_gap])
                if near >= fill_gap:
                    sites.append(pos); occ.insert(j, pos)
            sites = sorted(sites)
            print(f"NOTE TP {contig}:{rs}-{re}: TTAAAA short, topped up to {len(sites)}/{want}")
        recs = []
        for site in sites:
            family, variant = rng.choice(weighted)
            elem_seq, meta = build_true_element(fa, family, variant, rng, cache, pseudo)
            strand = rng.choice(["+", "-"])
            ins = elem_seq if strand == "+" else rc(elem_seq)
            ip = site + meta["tsd"]
            klass = "pseudogene" if family == "PSEUDOGENE" else ("young" if family in YOUNG_FAMS else "old")
            recs.append({"ip": ip, "ins": ins, "tsd_seq": region[site:ip],
                         "family": meta["sub"], "variant": variant, "tsd": meta["tsd"],
                         "strand": strand, "elem_len": len(elem_seq), "class": klass, "role": "TP"})
        emit_window(contig, rs, re, "TP", recs)
        n_tp += len(recs)
    print(f"implanted {n_tp} TRUE insertions (all classes)")

    # ---- DECOY implants into FP windows ---------------------------------------------------
    fp = read_bed(args.fp_bed)
    fp_len = [e - s for _, s, e, _ in fp]
    tot_fp = sum(fp_len) or 1
    n_fp = 0
    for wi, (contig, rs, re, _) in enumerate(fp):
        want = (args.n_decoys - n_fp) if wi == len(fp) - 1 else round(args.n_decoys * fp_len[wi] / tot_fp)
        region = fetch(fa, contig, rs + 1, re)
        cand = list(range(5000, len(region) - 5000))
        rng.shuffle(cand)
        pts = spaced(sorted(cand[:want * 4]), want, args.fp_gap)
        if len(pts) < want:
            print(f"WARN FP {contig}:{rs}-{re}: {len(pts)} slots < {want} wanted")
        recs = []
        for ip in pts:
            kind = rng.random()
            if kind < 0.45 and lcr_pool:
                ins, klass, var = build_lcr_insert(lcr_pool, rng), "lcr", "expansion"
            elif kind < 0.72 and young_pool:
                ins, klass, var = build_frag_insert(young_pool, rng), "young_rte", "reffrag"
            elif old_pool:
                ins, klass, var = build_frag_insert(old_pool, rng), "old_rte", "reffrag"
            else:
                ins, klass, var = build_lcr_insert(lcr_pool or ["ACGT"], rng), "lcr", "expansion"
            recs.append({"ip": ip, "ins": ins, "tsd_seq": "", "family": klass,
                         "variant": var, "tsd": 0, "strand": "+", "elem_len": len(ins),
                         "class": klass, "role": "FP"})
        emit_window(contig, rs, re, "FP", recs)
        n_fp += len(recs)
    print(f"planted {n_fp} non-MEI decoy loci (lcr/young_rte/old_rte)")

    donor.close(); flanks.close(); truth.close()
    print(f"total: {tid} insertions  (TP={n_tp}, FP-decoy={n_fp})")


if __name__ == "__main__":
    main()
