import random, gzip, sys, json
import pyliftover
random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 1)
out = sys.argv[1]
srcs = {"chr1": 100000, "chr2": 5000, "chrX": 777, "chrEmpty": 50, "1": 3000}
tgts = {"chr1": 90000, "chrA": 20000, "chrB": 1234, "chr2": 6000}
lines = []
cid = 0
for _ in range(400):
    s = random.choice(list(srcs)); t = random.choice(list(tgts))
    ssz, tsz = srcs[s], tgts[t]
    nb = random.randint(1, 12)
    sizes = [random.randint(0 if random.random()<0.05 else 1, 300) for _ in range(nb)]
    sgaps = [random.randint(0, 200) for _ in range(nb-1)]
    tgaps = [random.randint(0, 200) for _ in range(nb-1)]
    tot_s = sum(sizes)+sum(sgaps); tot_t = sum(sizes)+sum(tgaps)
    if tot_s > ssz or tot_t > tsz: continue
    sstart = random.randint(0, ssz-tot_s); tstart = random.randint(0, tsz-tot_t)
    strand = random.choice("+-")
    score = random.choice([100, 200, 200, 500, 500, 500, 1000])
    cid += 1
    # chain target strand '-' : tstart coordinates are on the reverse strand per python code (no conversion at parse)
    lines.append(f"chain {score} {s} {ssz} + {sstart} {sstart+tot_s} {t} {tsz} {strand} {tstart} {tstart+tot_t} {cid}")
    for i in range(nb-1):
        lines.append(f"{sizes[i]} {sgaps[i]} {tgaps[i]}")
    lines.append(f"{sizes[-1]}")
    lines.append("")
txt = "\n".join(lines)+"\n"
with gzip.open(out, "wt") as f: f.write(txt)
lo = pyliftover.LiftOver(out)
qs = []
for c, sz in list(srcs.items()) + [("chrNone", 100)]:
    for _ in range(3000):
        qs.append((c, random.randint(-5, sz+5)))
with open(out + ".expected", "w") as f:
    for c, p in qs:
        r = lo.convert_coordinate(c, p)
        f.write(json.dumps([c, p, None if r is None else [list(x) for x in r]]) + "\n")
print(len(lines), "lines", len(qs), "queries")
