#!/usr/bin/env python3
"""E2E helper: the exon tracks annotate needs so that its processed-pseudogene proof uses the SAME
gene models the simulator built its pseudogene parents from (truth and annotate must agree).

Source of the gene models, in order:
  1. `<donor>/genes_hs1.tsv` written by test/fullstack/donor_types.py (every parent gene, hs1,
     tools/build_gene_model.py track format: contig start end gene strand, 0-based half-open);
  2. else recovered from the truth: every `EXON<n>[gene]` part of a PSEUDOGENE / PSEUDOGENE_DECOY
     event (truth_types.ins.fa) is located on hs1 near the event's `gene_loc` (exact match).

Writes
  --out-hs1   exon track on hs1 (= the remap genome `remap_2bit`): tools/rte builds the exon-exon
              junction cores from it (CONFIG['annotate']['rte_exon_annotation'])
  --out-ref   the same exons placed on the E2E clip-remap reference (`--ref`, the reduced GRCh38;
              placed by mapping the hs1 exon sequence with minimap2): annotate_v2's clip-exon
              candidate test (CONFIG['annotate']['exon_annotation']).
In a real run both are ONE hs1 track (the clip remap genome IS hs1): see E2E_REPORT.md."""
import argparse
import csv
import os
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, REPO)
from tools.rte.genome import open_genome  # noqa: E402
from tools.rte.sequtil import rc  # noqa: E402

PART = re.compile(r"^(EXON\d+):(\d+)-(\d+)\[([^\]]+)\]$")


def read_fa(path):
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split("|")[0]
                out[name] = []
            elif name is not None:
                out[name].append(line)
    return {k: "".join(v).upper() for k, v in out.items()}


def recover(donor, hs1):
    truth = list(csv.DictReader(open(os.path.join(donor, "truth_types_hs1.tsv")), delimiter="\t"))
    ins = read_fa(os.path.join(donor, "truth_types.ins.fa"))
    loc = {}
    for r in truth:
        m = re.search(r"gene_loc=([^:;]+):(\d+)-(\d+):([+-])", r["info"])
        if m:
            g = re.search(r"gene=([^;]+)", r["info"]).group(1)
            loc[g] = (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
    exons = {}
    for r in truth:
        if r["variant"] not in ("PSEUDOGENE", "PSEUDOGENE_DECOY") or r["id"] not in ins:
            continue
        seq = ins[r["id"]]
        for p in r["parts"].split(";"):
            m = PART.match(p)
            if not m:
                continue
            g = m.group(4)
            gl = loc.get(g)
            if gl is None:
                continue
            c, s, e, st = gl
            piece = seq[int(m.group(2)):int(m.group(3))]
            if len(piece) < 20:
                continue
            win_lo = max(0, s - 5000)
            win = hs1.fetch(c, win_lo, e + 5000)
            q = piece if st == "+" else rc(piece)
            i = win.find(q)
            if i < 0:
                continue
            exons.setdefault((g, c, st), set()).add((win_lo + i, win_lo + i + len(q)))
    rows = []
    for (g, c, st), ivs in exons.items():
        for a, b in _merge(sorted(ivs)):
            rows.append((c, a, b, g, st))
    return rows


def _merge(ivs):
    out = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def place_on_ref(rows, hs1, ref_fa):
    import mappy
    al = mappy.Aligner(ref_fa, preset="sr")
    out = []
    for c, a, b, g, st in rows:
        seq = hs1.fetch(c, a, b)
        best = None
        for h in al.map(seq):
            if h.q_en - h.q_st < 0.9 * len(seq) or h.mlen < 0.9 * h.blen:
                continue
            if best is None or h.mapq > best.mapq:
                best = h
        if best is not None:
            out.append((best.ctg, best.r_st, best.r_en, g, st))
    return out


def write(rows, path):
    with open(path, "w") as fh:
        for r in sorted(rows):
            fh.write("\t".join(map(str, r)) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor-dir", required=True)
    ap.add_argument("--hs1-2bit", required=True)
    ap.add_argument("--ref", required=True, help="E2E clip-remap reference FASTA (reduced.fa)")
    ap.add_argument("--out-hs1", required=True)
    ap.add_argument("--out-ref", required=True)
    a = ap.parse_args()
    hs1 = open_genome(a.hs1_2bit)
    gpath = os.path.join(a.donor_dir, "genes_hs1.tsv")
    if os.path.exists(gpath):
        rows = []
        with open(gpath) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5:
                    rows.append((f[0], int(f[1]), int(f[2]), f[3], f[4]))
        how = "genes_hs1.tsv"
    else:
        rows = recover(a.donor_dir, hs1)
        how = "recovered from truth EXON parts"
    write(rows, a.out_hs1)
    ref_rows = place_on_ref(rows, hs1, a.ref)
    write(ref_rows, a.out_ref)
    print(f"gene track ({how}): {len(rows)} exons of {len({r[3] for r in rows})} genes on hs1, "
          f"{len(ref_rows)} placed on {os.path.basename(a.ref)}")


if __name__ == "__main__":
    main()
