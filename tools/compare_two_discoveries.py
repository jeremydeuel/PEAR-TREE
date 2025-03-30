import gzip

d1 = "../../spar/human/discovery/PD45517.insertions.combined.old.txt.gz"
d2 = "../../spar/human/discovery/PD45517.insertions.combined.txt.gz"

i1 = {}
i2 = {}

with gzip.open(d1, 'rt') as f:
    insertion = None
    step = None
    for line in f:
        line = line.strip()
        if not line: continue
        if line[0] == ">":
            insertion = line[1:]
            i1[insertion] = []
            step = "none"
            continue
        if line[0] == "@":
            step = line[1:]
            continue
        if step == "FILES":
            i1[insertion].append(line)

with gzip.open(d2, 'rt') as f:
    insertion = None
    step = None
    for line in f:
        line = line.strip()
        if not line: continue
        if line[0] == ">":
            insertion = line[1:]
            i2[insertion] = []
            step = "none"
            continue
        if line[0] == "@":
            step = line[1:]
            continue
        if step == "FILES":
            i2[insertion].append(line)

ins1 = set(i1.keys())
ins2 = set(i2.keys())
both = ins1.intersection(ins2)
print(f"d1={len(ins1)}, d2={len(ins2)}, shared={len(both)}")

s1 = [i for i in i1.keys() if len(i1[i])==1]
s2 = [i for i in i2.keys() if len(i2[i])==1]
print(f"singles in f1={len(s1)}, singles in f2={len(s2)}")

s1 = [i for i in both if len(i1[i])==1]
s2 = [i for i in both if len(i2[i])==1]
print(f"shared singles in f1={len(s1)}, singles in f2={len(s2)}")

s1 = [i for i in ins1.difference(both) if len(i1[i])==1]
s2 = [i for i in ins2.difference(both) if len(i2[i])==1]
print(f"new singles in f1={len(s1)}, singles in f2={len(s2)}")

s1 = [i for i in ins1.difference(both) if len(i1[i])>1 and len(i1[i])<10]
s2 = [i for i in ins2.difference(both) if len(i2[i])>1 and len(i2[i])<10]
print(f"new multi-insertions in f1={len(s1)}, singles in f2={len(s2)}")
