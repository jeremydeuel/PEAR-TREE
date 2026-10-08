"""Generates the WP-TD unit-test fixtures + python expectations (run from the repo root):

    PYTHONPATH=. PY rust/peartree-rte/golden/make_td_expect.py

needs src/config.py (gitignored; copy it in temporarily) because tools/annotate_v2 imports it.
Writes rust/peartree-rte/golden/data/td/{track.tsv, rmsk.out, rmsk.out.gz, rmsk_ucsc.txt,
expect.json}: tests/td_python_parity.rs replays them against the Rust port.
"""
import gzip
import json
import os
import sys

sys.path.insert(0, "tools")
from annotate_v2 import GeneModel                      # noqa: E402
from rte.pseudogene import load_exons_by_gene, load_gene_strands   # noqa: E402
from rte.transduction import L1Rmsk                    # noqa: E402

OUT = "rust/peartree-rte/golden/data/td"
os.makedirs(OUT, exist_ok=True)

# ---- exon / gene track (contig start end gene strand), deliberately messy
TRACK = """# comment
chr1\t1000\t1200\tGA\t+
chr1\t2000\t2100\tGA\t+
chr1\t2050\t2300\tGA\t+
chr1\t5000\t5400\tGA\t+
chr1\t100\t180\tGP\t+
chr1\t8000\t8200\tGM\t-
chr1\t9000\t9100\tGM\t-
chr1\t10000\t10200\tGM\t-
chr1\t7000\t30000\tLONG\t+
chr1\t7000\t7050\tLONG\t+
chr1\t29900\t30000\tLONG\t+
chr1\t8100\t8150\tOVER\t-
chr1\t8100\t8150\tOVER\t+
chr1\t40000\t40500\tSINGLE\t-

chr1\t50000\t50010\tSHORTLINE
chr1\t60000\t60100\tDUP\t+
chr1\t60000\t60100\tDUP\t+
chr1\t6000\t6100\tSAMESTART\t+
chr1\t6000\t6100\tSAMESTA2\t-
2\t3000\t3100\tNUM\t+
2\t3500\t3600\tNUM\t+
2\t3500\t3600\tNUM\t+
2\t5000\t6000\tNUM2\t-
2\t5100\t5300\tNUM2\t-
2\t5400\t5500\tNUM2\t-
chrX\t100\t200\tXG\t-
"""
open(f"{OUT}/track.tsv", "w").write(TRACK)

gm = GeneModel(f"{OUT}/track.tsv", {"splice_donor_window": 6, "splice_acceptor_window": 3,
                                    "splice_ppt_window": 17, "splice_branch_window": 45,
                                    "promoter_up": 2000})
cfg2 = {"splice_donor_window": 9, "splice_acceptor_window": 5, "splice_ppt_window": 30,
        "splice_branch_window": 80, "promoter_up": 300}
gm2 = GeneModel(f"{OUT}/track.tsv", cfg2)

exp = {"cfg2": cfg2}
exp["genes"] = {c: [[g[0], g[1], g[2], g[3], [list(x) for x in g[4]]] for g in v] for c, v in gm.genes.items()}
exp["maxspan"] = gm._maxspan

# candidates + genic_feature over a position grid, several contig spellings
contigs = ["chr1", "1", "chr2", "2", "chrX", "X", "chrY", "MT"]
cand, feat, feat2 = [], [], []
pts = list(range(-50, 600, 37)) + list(range(900, 10500, 53)) + list(range(29800, 32500, 97)) \
    + list(range(39000, 42800, 211)) + list(range(2800, 6200, 41)) + list(range(60000 - 2500, 60100 + 2500, 211)) \
    + list(range(49000, 53000, 777))
for c in contigs:
    for p in pts:
        r = gm._candidates(c, p)
        cand.append([c, p, None if r is None else [[g[0], g[1], g[2], g[3]] for g in r]])
for c in contigs:
    for p in pts:
        r = gm._candidates(c, p) or []
        for g in r:
            if g[0] <= p < g[1]:
                feat.append([g[2], g[3], p, list(gm._genic_feature(g[4], g[3], p))])
                feat2.append([g[2], g[3], p, list(gm2._genic_feature(g[4], g[3], p))])
exp["candidates"] = cand
exp["features"] = feat
exp["features2"] = feat2
exp["resolve"] = [[c, gm._resolve(c)] for c in contigs]
# dense splice-class grid
exp["splice"] = [[dd, da, list(gm._splice_class(dd, da)), list(gm2._splice_class(dd, da))]
                 for dd in (1, 3, 6, 7, 9, 10, 100) for da in (1, 3, 4, 5, 6, 17, 18, 30, 31, 45, 46, 80, 81, 100)]
# denser intron walk of every multi-exon gene
walk = []
for c, v in gm.genes.items():
    for g in v:
        for p in range(g[0] - 1, g[1] + 2):
            if g[0] <= p < g[1] and len(g[4]) > 1 and p % 3 == 0:
                walk.append([g[2], g[3], p, list(gm._genic_feature(g[4], g[3], p)), list(gm2._genic_feature(g[4], g[3], p))])
exp["walk"] = walk

# ---- pseudogene loaders
exp["exons_by_gene"] = {g: [list(x) for x in v] for g, v in load_exons_by_gene(f"{OUT}/track.tsv").items()}
exp["gene_strands"] = load_gene_strands(f"{OUT}/track.tsv")

# ---- L1Rmsk: RepeatMasker .out (15/16 fields), gz, and UCSC rmsk.txt layouts
def out_line(score, div, contig, b, e, strand, name, fam, extra=""):
    return (f"  {score:5d} {div:4.1f}  1.0  2.0 {contig:>8s} {b:9d} {e:9d} (1000) {strand} {name} {fam}"
            f"  1 6000 (0) {7}{extra}")

rows = [
    ("chr1", 101, 6200, "+", "L1HS", "LINE/L1"),
    ("chr1", 20001, 26100, "C", "L1PA2", "LINE/L1"),
    ("chr1", 30001, 30500, "+", "L1HS", "LINE/L1"),          # too short
    ("chr1", 40001, 46000, "+", "AluY", "SINE/Alu"),         # not L1
    ("chr1", 50001, 56000, "+", "L1PA3", "LINE/L1"),
    ("chr1", 60001, 66000, "C", "L1PA3", "LINE/L1"),
    ("chr1", 60001, 66000, "C", "L1PA2", "LINE/L1"),         # same span, name tiebreak
    ("chr1", 70001, 76000, "+", "L1ME", "LINE/L1-dep"),
    ("2", 1001, 7000, "+", "L1HS", "LINE/L1"),
    ("2", 9001, 15000, "C", "L1HS", "LINE/L1"),
]
lines = ["   SW  perc perc perc  query      position in query           matching       repeat              class/family",
         "score  div. del. ins.  sequence    begin     end    (left)    repeat         position in repeat",
         ""]
for i, (c, b, e, s, n, f) in enumerate(rows):
    lines.append(out_line(1000 + i, 1.5 + i * 0.7, c, b, e, s, n, f, "  *" if i % 3 == 0 else ""))
txt = "\n".join(lines) + "\n"
open(f"{OUT}/rmsk.out", "w").write(txt)
with gzip.open(f"{OUT}/rmsk.out.gz", "wt") as fh:
    fh.write(txt)
ucsc = []
for i, (c, b, e, s, n, f) in enumerate(rows):
    cls, fam = f.split("/", 1)
    ucsc.append("\t".join(map(str, [585, 1000 + i, 15 + i * 7, 5, 4, c, b - 1, e, -1000, "-" if s == "C" else s, n,
                                    cls, fam, 1, 6000, 0, i])))
open(f"{OUT}/rmsk_ucsc.txt", "w").write("\n".join(ucsc) + "\n")


def dump(r):
    return {c: [[x[0], x[1], x[2], x[3], x[4]] for x in v] for c, v in r.by_contig.items()}


exp["rmsk"] = {}
queries = []
for c in ["chr1", "1", "chr2", "2", "chrZ"]:
    for s in range(0, 80000, 1700):
        for strand in "+-":
            queries.append([c, s, s + 400, strand])
for name, path in [("out", "rmsk.out"), ("gz", "rmsk.out.gz"), ("ucsc", "rmsk_ucsc.txt")]:
    r = L1Rmsk(f"{OUT}/{path}", 5500)
    res = []
    for c, s, e, strand in queries:
        res.append(r.upstream_of(c, s, e, strand, 15000))
    exp["rmsk"][name] = {"rows": dump(r), "queries": queries,
                         "up": [[list(x) for x in q] for q in res],
                         "up_small": [[list(x) for x in r.upstream_of(c, s, e, strand, 3000)] for c, s, e, strand in queries]}

json.dump(exp, open(f"{OUT}/expect.json", "w"), indent=None, separators=(",", ":"))
print("candidates", len(cand), "features", len(feat), "walk", len(walk), "queries", len(queries))
