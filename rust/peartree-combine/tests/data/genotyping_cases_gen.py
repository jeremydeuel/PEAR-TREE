import sys, types, random, json, re, os
SRC = "/Users/jeremy/Documents/PEAR_TREE-rc-p6/src"
sys.path.insert(0, SRC)
class _D(dict):
    def __missing__(self, k): return 10
cfg = types.ModuleType("config"); cfg.CONFIG = _D({"genotyping": {"max_bases": 30}, "discovery": _D()})
sys.modules["config"] = cfg
from quality_seq import QualitySeq
from revcomp import revcomp
from sequence_checks import sequence_matching_score
TYPE_FULL_INFO, TYPE_RIGHT_POLYA, TYPE_LEFT_POLYA, TYPE_RIGHT_DISC, TYPE_LEFT_DISC = 3, 1, 2, 4, 5
MAXB = 30
CONFIG = cfg.CONFIG
random.seed(7)
genome = {}
for c in ["chr1", "chr2", "chrM_long_name_x"]:
    s = "".join(random.choice("ACGT") for _ in range(2000))
    if c == "chr2":
        s = s[:500] + "N" * 40 + s[540:]
    genome[c] = s
genome["chr3"] = "".join(random.choice("ACGT") for _ in range(300))
def get_sequence(name, start, end):
    s = genome.get(name)
    if s is None or start >= end: return ""
    end = min(end, len(s))
    if start < 0 or start >= end: return ""
    return s[start:end]
def mut(s, p):
    return "".join(random.choice("ACGT") if random.random() < p else c for c in s)
def rnd(n): return "".join(random.choice("ACGTN" if random.random()<0.01 else "ACGT") for _ in range(n))
class Ins:  pass
cases = []
for k in range(400):
    c = random.choice(["chr1", "chr1", "chr2", "chr3", "chrM_long_name_x", "chrZ"])
    L = len(genome.get(c, "A"*100))
    ty = random.choice([3, 3, 3, 1, 2, 4, 5])
    lp = random.randint(-5, L+10) if random.random() < 0.15 else random.randint(40, max(41, L-60))
    rp = lp + random.randint(0, 12)
    if c == "chr2" and random.random() < 0.3: lp = rp = random.randint(480, 560)
    ins = Ins()
    ins.reference_name = c
    ins.name = f"{c}:{lp}-{rp}"
    ins.type = ty
    ins.left_pos, ins.right_pos = lp, rp
    # clipped sequences; sometimes derived from the reference (to trigger sms > 0 exclusion)
    def clip(side):
        n = random.choice([5, 12, 25, 40, 60])
        mode = random.random()
        r = genome.get(c, "")
        if mode < 0.3 and r:
            if side == "R":
                base = r[max(0, rp):max(0, rp)+n]
            else:
                base = r[max(0, lp-n):max(0, lp)]
            base = mut(base, 0.05) if base else rnd(n)
            if not base: base = rnd(n)
            return base
        return rnd(n)
    rc_ = clip("R"); lc_ = clip("L")
    if random.random() < 0.3: rc_ = rc_.lower()
    if random.random() < 0.3: lc_ = revcomp(lc_)
    ins.right_clipped = QualitySeq(rc_, [random.randint(2, 40) for _ in rc_])
    ins.left_clipped = QualitySeq(lc_, [random.randint(2, 40) for _ in lc_])
    cases.append(ins)
# identical-sequence case: left_ref == right_ref  (rp == lp -> left_ref = ref[lp-30:lp], right_ref = ref[lp:lp+30]; not equal) -> use a repeat genome region
X = genome["chr1"][100:130]
genome["chr1"] = genome["chr1"][:900] + X + X + revcomp(X) + genome["chr1"][990:]
for lp in (930, 960):
    ins = Ins(); ins.reference_name = "chr1"; ins.name = f"chr1:{lp}-{lp}"; ins.type = 3; ins.left_pos = ins.right_pos = lp
    ins.right_clipped = QualitySeq(rnd(20), [30]*20); ins.left_clipped = QualitySeq(rnd(20), [30]*20)
    cases.append(ins)
filter_names = {cases[3].name, cases[10].name}
src = open(SRC + "/combine_insertions.py").read()
a = src.index("    n_excluded = 0\n    n_included = 0\n    with gzip.open(insertions_genotyping_file")
body = src[a:].replace("with gzip.open(insertions_genotyping_file, 'wt') as f:", "with open(insertions_genotyping_file, 'w') as f:")
ns = dict(TYPE_FULL_INFO=3, TYPE_RIGHT_POLYA=1, TYPE_LEFT_POLYA=2, TYPE_RIGHT_DISC=4, TYPE_LEFT_DISC=5,
          CONFIG=CONFIG, get_sequence=get_sequence, revcomp=revcomp, sequence_matching_score=sequence_matching_score)
exec("def run(insertions, filter_reads, insertions_genotyping_file):\n" + body, ns)
import io, contextlib
outp = sys.argv[1]
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ns["run"](cases, filter_names, outp + ".expected.txt")
import collections
cnt = collections.Counter()
for l in buf.getvalue().splitlines():
    m = re.search(r"(due to [a-zA-Z ]+|since the breakpoint|wrote)", l)
    cnt[m.group(1) if m else l[:30]] += 1
print(cnt)
data = {
    "max_bases": MAXB,
    "genome": genome,
    "filter": sorted(filter_names),
    "insertions": [{"name": i.name, "contig": i.reference_name, "type": i.type, "left_pos": i.left_pos, "right_pos": i.right_pos,
                    "left_clipped": str(i.left_clipped), "right_clipped": str(i.right_clipped)} for i in cases],
    "expected": open(outp + ".expected.txt").read(),
}
json.dump(data, open(outp, "w"))
