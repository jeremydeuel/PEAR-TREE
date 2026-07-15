#!/usr/bin/env python3
"""Lift the implanted-MEI truth from hs1/donor coordinates into hg38.

Two methods, in order of robustness:

1. **Chain lift (preferred, --chain):** convert the implant's hs1 insertion point
   (`hs1_contig`:`hs1_pos`) straight to hg38 with the UCSC hs1->hg38 liftOver chain. This
   is a coordinate map, so it works even where the surrounding sequence is repeat-dense —
   exactly where method 2 fails. The hg38 breakpoint brackets the lifted point by the TSD.

2. **Flank mapping (fallback):** each implant's 300 bp hs1 left/right flanks were mapped to
   hg38 (flanks.bam); the junction sits between the left flank's 3' end and the right
   flank's 5' start. This fails in LINE/segdup-dense regions (flanks map at low MAPQ), which
   is why chain lifting is preferred when a chain is supplied.

An implant is `scoreable` if either method places it; otherwise `unmappable` (genuinely no
hg38 homolog — e.g. a T2T-only region) and scoring ignores it. The `lift_method` column
records which method placed each row. Without --chain the behaviour is flank-only (legacy).
"""
import argparse
import pysam

try:
    import pyliftover
except ImportError:
    pyliftover = None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flanks-bam", required=True)
    p.add_argument("--truth-hs1", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--chain", default=None,
                   help="hs1->hg38 liftOver chain (e.g. hs1ToHg38.over.chain.gz). When set, "
                        "the insertion point is chain-lifted first; flank mapping is the fallback.")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--max-gap", type=int, default=5000)
    args = p.parse_args()

    # collect best primary alignment per flank id/side (method 2, fallback)
    aln = {}
    bam = pysam.AlignmentFile(args.flanks_bam)
    for r in bam:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        rid, side = r.query_name.rsplit("_", 1)
        aln[(rid, side)] = (r.reference_name, r.reference_start, r.reference_end, r.mapping_quality)

    lo = None
    if args.chain:
        if pyliftover is None:
            raise SystemExit("--chain given but pyliftover is not installed (pip install pyliftover)")
        lo = pyliftover.LiftOver(args.chain)

    rows = []
    with open(args.truth_hs1) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))

    n_chain = n_flank = n_bad = 0
    with open(args.out, "w") as out:
        out.write("id\tfamily\tvariant\ttsd\tstrand\telem_len\t"
                  "hg38_contig\thg38_left\thg38_right\tflank_mapq\tstatus\tlift_method\n")
        for r in rows:
            rid = r["id"]
            status = "unmappable"; method = "none"
            contig = left = right = mq = "."

            # method 1: chain-lift the insertion point (robust in repeat-dense regions)
            if lo is not None and r.get("hs1_contig") and r.get("hs1_pos"):
                hits = lo.convert_coordinate(r["hs1_contig"], int(r["hs1_pos"]))
                if hits:
                    contig, pos = hits[0][0], hits[0][1]
                    tsd = int(r["tsd"]) if r.get("tsd") else 0
                    left, right = pos, pos + tsd    # bracket the TSD like the L-R breakpoint pair
                    status, method, mq = "scoreable", "chain", "."

            # method 2 (fallback): flank mapping
            if status == "unmappable":
                l = aln.get((rid, "L")); rr = aln.get((rid, "R"))
                if l and rr and l[0] == rr[0] and min(l[3], rr[3]) >= args.min_mapq:
                    glo, ghi = sorted([l[2], rr[1]])
                    if ghi - glo <= args.max_gap:
                        status, method = "scoreable", "flank"
                        contig, left, right, mq = l[0], glo, ghi, min(l[3], rr[3])

            if status == "scoreable":
                n_chain += method == "chain"; n_flank += method == "flank"
            else:
                n_bad += 1
            out.write(f"{rid}\t{r['family']}\t{r['variant']}\t{r['tsd']}\t{r['strand']}\t"
                      f"{r['elem_len']}\t{contig}\t{left}\t{right}\t{mq}\t{status}\t{method}\n")
    print(f"lifted truth: {n_chain + n_flank} scoreable ({n_chain} chain, {n_flank} flank), "
          f"{n_bad} unmappable (no hg38 homolog)")


if __name__ == "__main__":
    main()
