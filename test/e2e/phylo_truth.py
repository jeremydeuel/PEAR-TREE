#!/usr/bin/env python3
"""Truth table for tools/phylo/calibration.py from a tree-mode E2E run (test/e2e/run_phylo_e2e.sh).

Maps every genotyped locus (names in the per-colony genotype files) to the simulated event whose
lifted breakpoints (results_by_event.tsv: ref_contig / ref_lo / ref_hi) lie within --window bp of
its breakpoints (order-free; a one-sided locus needs its one real breakpoint within the window of
either event breakpoint), and writes `locus  truth  carriers  event_id  variant  placement`:

  ROOT placement           -> germline   (somatic in every colony; purity-diluted)
  GERMLINE placement       -> germline   (germline het: VAF 0.5 in every colony)
  branch, 1 carrier        -> private
  branch, >= 2 carriers    -> clade
  NONCLADE                 -> nonclade   (real insertion sequence, carriers not a clade)
  NONCLADE_ARTEFACT        -> nonclade   (systematic artefact: same tract slips in many colonies)
  ARTEFACT                 -> artefact   (per-library, one colony)
  no event                 -> none       (e.g. hs1-vs-GRCh38 germline differences)
"""
import argparse
import csv
import gzip
import os
import re

NAME_RE = re.compile(r"^([^:]+):(polyA_|disc_|oneside_)?(\d+)-(polyA_|disc_|oneside_)?(\d+)$")


def parse(name):
    m = NAME_RE.match(name)
    if not m:
        return None
    a, b = int(m.group(3)), int(m.group(5))
    if m.group(2):          # oneside_R-R: real = right
        return m.group(1), [b]
    if m.group(4):          # L-oneside_L: real = left
        return m.group(1), [a]
    return m.group(1), [a, b]


def loci_from_dir(d):
    out = set()
    for f in os.listdir(d):
        if f.endswith(".txt.gz"):
            with gzip.open(os.path.join(d, f), "rt") as fh:
                next(fh)
                for line in fh:
                    out.add(line.split("\t", 1)[0])
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--e2e-dir", required=True)
    ap.add_argument("--genotype-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=30)
    a = ap.parse_args()
    O = a.e2e_dir
    pl = {}
    with open(os.path.join(O, "donor", "phylo_truth.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pl[r["id"]] = (r["placement"], [f"S{x}" for x in r["carriers"].split(",") if x])
    events = []
    with open(os.path.join(O, "results_by_event.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["ref_contig"] in ("", ".", "None") or r["id"] not in pl:
                continue
            events.append((r["id"], r["variant"], r["ref_contig"], int(r["ref_lo"]), int(r["ref_hi"])))
    by_contig = {}
    for e in events:
        by_contig.setdefault(e[2], []).append(e)
    n_match = 0
    with open(a.out, "w") as out:
        out.write("locus\ttruth\tcarriers\tevent_id\tvariant\tplacement\n")
        for name in loci_from_dir(a.genotype_dir):
            p = parse(name)
            best = None
            if p:
                c, bps = p
                for (eid, var, ec, lo, hi) in by_contig.get(c, []):
                    if len(bps) == 2:
                        d = min(max(abs(bps[0] - lo), abs(bps[1] - hi)), max(abs(bps[0] - hi), abs(bps[1] - lo)))
                    else:
                        d = min(abs(bps[0] - lo), abs(bps[0] - hi))
                    if d <= a.window and (best is None or d < best[0]):
                        best = (d, eid, var)
            if best is None:
                out.write(f"{name}\tnone\t\t.\t.\t.\n")
                continue
            n_match += 1
            placement, car = pl[best[1]]
            if placement in ("ROOT", "GERMLINE"):
                t = "germline"
            elif placement.startswith("NONCLADE"):
                t = "nonclade"
            elif placement == "ARTEFACT":
                t = "artefact"
            else:
                t = "clade" if len(car) > 1 else "private"
            out.write(f"{name}\t{t}\t{','.join(car)}\t{best[1]}\t{best[2]}\t{placement}\n")
    print(f"{n_match} genotyped loci matched to {len(events)} lifted events -> {a.out}")


if __name__ == "__main__":
    main()
