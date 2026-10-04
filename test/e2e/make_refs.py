#!/usr/bin/env python3
"""E2E helper: side files combine/annotate need for the reduced discovery reference.

* `<out>/identity.chain.gz`  UCSC chain mapping every contig of the reduced reference onto
  itself. In the E2E the clip re-map genome (bowtie2_index2) IS the discovery reference, so
  combine's "clip maps near the breakpoint" liftover is the identity.
* `<out>/reduced.rmsk.out.gz`  RepeatMasker .out rows (the format annotate_v2.read_rmsk parses)
  for the reduced reference's chr contigs, converted from UCSC hg38 rmsk.txt.gz.
"""
import argparse
import gzip
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="reduced reference FASTA (with .fai)")
    ap.add_argument("--hg38-rmsk", required=True, help="UCSC hg38 rmsk.txt.gz")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    contigs = []
    with open(a.ref + ".fai") as fh:
        for line in fh:
            f = line.split("\t")
            contigs.append((f[0], int(f[1])))
    chain = os.path.join(a.out_dir, "identity.chain.gz")
    with gzip.open(chain, "wt") as out:
        for i, (c, n) in enumerate(contigs, 1):
            out.write(f"chain 1000 {c} {n} + 0 {n} {c} {n} + 0 {n} {i}\n{n}\n\n")
    names = {c for c, _ in contigs}
    rmsk = os.path.join(a.out_dir, "reduced.rmsk.out.gz")
    if not os.path.exists(rmsk):
        n = 0
        with gzip.open(a.hg38_rmsk, "rt") as fh, gzip.open(rmsk, "wt") as out:
            out.write("   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat\n")
            out.write("score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID\n")
            out.write("\n")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                # bin swScore milliDiv milliDel milliIns genoName genoStart genoEnd genoLeft strand
                # repName repClass repFamily repStart repEnd repLeft id
                if f[5] not in names:
                    continue
                strand = "+" if f[9] == "+" else "C"
                fam = f"{f[11]}/{f[12]}"
                out.write(" ".join([f[1], f"{int(f[2]) / 10:.1f}", f"{int(f[3]) / 10:.1f}",
                                    f"{int(f[4]) / 10:.1f}", f[5], str(int(f[6]) + 1), f[7],
                                    f"({f[8].lstrip('-')})", strand, f[10], fam, f[13], f[14],
                                    f"({f[15].lstrip('-')})", f[16]]) + "\n")
                n += 1
        print(f"wrote {n} rmsk rows -> {rmsk}")
    print(f"wrote identity chain for {len(contigs)} contigs -> {chain}")


if __name__ == "__main__":
    main()
