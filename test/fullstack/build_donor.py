#!/usr/bin/env python3
"""Full-stack donor builder.

Constructs a DONOR genome from curated hs1 (T2T-CHM13v2.0) regions and implants
retrotransposon insertions (true positives), then the pipeline (run_fullstack.sh)
simulates reads from it with wgsim and maps them to **GRCh38** with bwa-mem. Mapping
hs1-derived reads to GRCh38 is deliberate: every hs1<->hg38 assembly difference
(T2T centromeres/telomeres, resolved segdups, novel sequence) becomes a mapping-driven
false positive that PEAR-TREE must reject, alongside the implanted true positives.

Curated regions (see REGIONS):
  * TP euchromatin windows (unique, good hg38 homology) -> implanted MEIs live here.
  * FP windows: acrocentric p-arm satellite/rDNA, pericentromeric alpha-satellite, and a
    telomeric end -> pure assembly-discordance / repeat artefact sources (no truth).

Preferential retrotransposition: L1 endonuclease nicks at the degenerate 5'-TTAAAA motif,
so implant sites are chosen at TTAAAA occurrences and flanked by a target-site duplication
(TSD) — the TPRT scar. Elements are REAL young copies extracted from hs1 (so clips carry
genuine Alu/L1/HERV-K sequence a repeat annotator would recognise), implanted full-length
and as the structural variants (5' truncation, 5' inversion, 3' transduction; HERV-K as a
no-poly-A LTR provirus).

Outputs: donor.fa (regions, TP windows implanted), truth_hs1.tsv (one row per implant,
in donor+hs1 coords), flanks.fa (300 bp hs1 flanks per implant; run_fullstack maps these to
hg38 to lift the truth into hg38 coordinates for scoring).
"""
import argparse
import random
import pysam

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]

# --- real young elements to implant (hs1 coords, 1-based inclusive, strand) ----------
# SVA loci picked from the hs1 (= chm13v2.0) RepeatMasker track (hs1.repeatMasker.out.gz,
# the local equivalent of UCSC chm13v2.0_rmsk.bb): for each of SVA_A..SVA_F the youngest
# full-length copy (1.3-2.2 kb, lowest % divergence) on an autosome. SVA_E/SVA_F are the
# youngest, human-specific, retrotranspositionally-active subfamilies (implanted below);
# SVA_A..D are registered as positive controls a caller should still detect if mobilised.
ELEMENT_LOCI = {
    "L1HS":  ("chrX", 11289796, 11295826, "+"),   # 6030 bp, div 0.3
    "AluYa5": ("chrX", 2759615, 2759925, "+"),     # 310 bp,  div 0.0
    "HERVK": ("chr10", 5039043, 5047200, "-"),     # 8157 bp HERVK13-int, div 0.4
    "SVA_A": ("chr7",   93267620,  93269590, "+"),  # 1970 bp, div 4.5
    "SVA_B": ("chr9",  130871189, 130872953, "-"),  # 1764 bp, div 3.1
    "SVA_C": ("chr11", 116360451, 116362140, "+"),  # 1689 bp, div 2.2
    "SVA_D": ("chr1",    5830038,   5831880, "-"),  # 1842 bp, div 1.8
    "SVA_E": ("chr3",  148851419, 148853373, "+"),  # 1954 bp, div 2.0
    "SVA_F": ("chr18",  54545509,  54547673, "+"),  # 2164 bp, div 1.1 (youngest)
}

# curated donor regions (hs1). role: "TP" implanted; "FP" artefact source only.
REGIONS = [
    ("chr21", 20_000_000, 26_000_000, "TP"),   # euchromatin
    ("chr22", 25_000_000, 31_000_000, "TP"),   # euchromatin
    ("chr21",  3_000_000,  6_000_000, "FP"),   # acrocentric p-arm satellite / rDNA
    ("chr21", 10_000_000, 13_000_000, "FP"),   # pericentromeric alpha-satellite
    ("chr22", 50_800_000, 51_324_926, "FP"),   # telomeric end of chr22
]

# per-family variant plan (implant counts are per TP region)
PLAN = [
    ("L1HS", "full", 3), ("L1HS", "5p_truncated", 3), ("L1HS", "5p_inversion", 2),
    ("L1HS", "3p_transduction", 2),
    ("AluYa5", "full", 4), ("AluYa5", "5p_truncated", 2),
    ("HERVK", "provirus", 1),
    # SVA_E/F: youngest, active subfamilies. TPRT hallmarks (poly-A tail + TSD) and the
    # common 5' truncation; the extracted hs1 sequence carries the real SVA hexamer/VNTR/
    # SINE-R structure a repeat annotator recognises.
    ("SVA_F", "full", 2), ("SVA_F", "5p_truncated", 1), ("SVA_E", "full", 1),
]
POLYA_MIN, POLYA_MAX = 15, 40
TSD_MIN, TSD_MAX = 8, 20
FLANK = 300


def fetch(fa, contig, s1, e1, strand="+"):
    seq = fa.fetch(contig, s1 - 1, e1).upper()
    return rc(seq) if strand == "-" else seq


def build_element(fa, family, variant, rng):
    """Return the inserted sequence (5'->3', to be placed on the + strand of the donor)."""
    c, s, e, st = ELEMENT_LOCI[family]
    elem = fetch(fa, c, s, e, st)
    polya = "A" * rng.randint(POLYA_MIN, POLYA_MAX)
    if family == "HERVK":                       # LTR provirus: no poly-A
        return elem, {"tsd": rng.randint(5, 6)}
    if variant == "full":
        body = elem
    elif variant == "5p_truncated":             # keep the 3' ~1.2-2 kb (5' loss)
        keep = min(len(elem), rng.randint(1200, 2000))
        body = elem[-keep:]
    elif variant == "5p_inversion":             # twin priming: invert a 5' segment
        inv = min(len(elem) - 100, rng.randint(800, 2500))
        body = rc(elem[:inv]) + elem[inv:]
    elif variant == "3p_transduction":          # co-mobilised unique downstream tag
        tag = "".join(rng.choice("ACGT") for _ in range(rng.randint(150, 400)))
        body = elem + tag
    else:
        body = elem
    return body + polya, {"tsd": rng.randint(TSD_MIN, TSD_MAX)}


def find_motif_sites(seq, motif="TTAAAA", n=None, rng=None, margin=5000):
    sites = []
    i = seq.find(motif, margin)
    while i != -1 and i < len(seq) - margin:
        sites.append(i)
        i = seq.find(motif, i + 1)
    if rng is not None:
        rng.shuffle(sites)
    return sites[:n] if n else sites


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hs1", default="/Users/jeremy/Downloads/hs1.fa")
    p.add_argument("--out-donor", required=True)
    p.add_argument("--out-truth", required=True)
    p.add_argument("--out-flanks", required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--min-site-gap", type=int, default=60_000,
                   help="minimum spacing between implant sites in a TP region")
    args = p.parse_args()

    rng = random.Random(args.seed)
    fa = pysam.FastaFile(args.hs1)

    donor = open(args.out_donor, "w")
    flanks = open(args.out_flanks, "w")
    truth = open(args.out_truth, "w")
    truth.write("id\tdonor_contig\tdonor_pos\ths1_contig\ths1_pos\tfamily\tvariant\t"
                "tsd\tstrand\telem_len\trole\n")

    tid = 0
    for (contig, rs, re, role) in REGIONS:
        region = fetch(fa, contig, rs + 1, re)     # 0-based region -> 1-based fetch
        dname = f"{contig}_{rs}_{re}_{role}"
        if role != "TP":
            donor.write(f">{dname}\n{region}\n")   # FP region: emit as-is (assembly-discordance source)
            continue

        # choose EN-motif implant sites, spaced apart
        sites = find_motif_sites(region, n=None, rng=rng)
        chosen, last = [], -1e9
        plan = [(f, v) for (f, v, k) in PLAN for _ in range(k)]
        rng.shuffle(plan)
        for s in sorted(sites):
            if s - last >= args.min_site_gap and plan:
                chosen.append(s); last = s
        chosen = chosen[:len(plan)]

        # implant from the 3' end backwards so earlier positions stay valid
        pieces = []
        cursor = len(region)
        implants = sorted(zip(chosen, plan), key=lambda x: x[0], reverse=True)
        records = []
        for site, (family, variant) in implants:
            elem_seq, meta = build_element(fa, family, variant, rng)
            tsd = meta["tsd"]
            strand = rng.choice(["+", "-"])
            ins = elem_seq if strand == "+" else rc(elem_seq)
            tsd_seq = region[site:site + tsd]
            # donor = ... region[:site+tsd] + INS + tsd_seq + region[site+tsd:] ...
            insert_point = site + tsd
            left_flank = region[max(0, insert_point - FLANK):insert_point]
            right_flank = region[insert_point:insert_point + FLANK]
            records.append((tid, dname, insert_point, contig, rs + insert_point,
                            family, variant, tsd, strand, len(elem_seq),
                            tsd_seq, ins, insert_point, left_flank, right_flank))
            tid += 1

        # assemble donor sequence with inserts (process ascending, build with offsets)
        out = []
        prev = 0
        for rec in sorted(records, key=lambda r: r[12]):
            ip = rec[12]; tsd_seq = rec[10]; ins = rec[11]
            out.append(region[prev:ip])          # up to and including the first TSD copy
            out.append(ins)                       # inserted element (+ poly-A)
            out.append(tsd_seq)                   # duplicated target site (TSD)
            prev = ip
        out.append(region[prev:])
        donor.write(f">{dname}\n{''.join(out)}\n")

        for rec in records:
            (rid, dnm, dpos, hc, hpos, fam, var, tsd, strand, elen,
             tsd_seq, ins, ip, lflank, rflank) = rec
            truth.write(f"{rid}\t{dnm}\t{dpos}\t{hc}\t{hpos}\t{fam}\t{var}\t{tsd}\t"
                        f"{strand}\t{elen}\t{role}\n")
            flanks.write(f">{rid}_L\n{lflank}\n>{rid}_R\n{rflank}\n")

    donor.close(); flanks.close(); truth.close()
    print(f"implanted {tid} MEIs across TP regions; donor -> {args.out_donor}")


if __name__ == "__main__":
    main()
