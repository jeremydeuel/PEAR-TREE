#!/usr/bin/env python3
import gzip, glob, os, random, statistics
from collections import Counter

DIR = "/Users/jeremy/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents - STEM_Green_Lab/HNRNPA1/JAK2_HNRNPA1_paper/previous analyses bin/trees_mpn/discovery"
files = sorted(glob.glob(os.path.join(DIR, "*.discovery.fq.gz")))
random.seed(7)
sample = random.sample(files, 25)

def revcomp(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

# Alu / L1 diagnostic k-mers (canonical consensus fragments)
ALU = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGG"
L1  = "GGAGGAGCCAAGATGGCCGAATAGGAACAGCTCCGGTCTACAGCTCCCAGCGTGAGCGA"  # L1HS 5' region
def looks_repeat(s):
    for k in range(0,len(ALU)-15,6):
        if ALU[k:k+15] in s: return 'Alu'
    for k in range(0,len(L1)-15,6):
        if L1[k:k+15] in s: return 'L1'
    return None

clip=[]; loci_side_polya=Counter()
by_side_field=Counter()
palindrome_examples=[]; polyg_examples=[]
for fp in sample:
    with gzip.open(fp,'rt') as fh:
        cur=None; data={}
        def flush(locus,data):
            for fk,seq in data.items():
                if 'CLIPPED' in fk:
                    clip.append(seq.upper())
                    side=fk.split(':')[0]
                    aln=data.get(f'{side}:ALIGNED')
                    u=seq.upper()
                    if aln and len(u)>=14 and u[:14] in revcomp(aln.upper()) and len(palindrome_examples)<6:
                        palindrome_examples.append((locus,side,u[:50],aln.upper()[:40]))
                    if u.endswith('GGGGGGGG') and len(polyg_examples)<6:
                        polyg_examples.append((locus,fk,u[:50]))
        while True:
            h=fh.readline()
            if not h: break
            h=h.strip()
            if not (h.startswith('@') or h.startswith('>')): continue
            seq=fh.readline().strip(); fh.readline(); fh.readline()
            p=h[1:].split(':')
            if len(p)<4: continue
            locus=f'{p[0]}:{p[1]}'; fk=f'{p[2]}:{p[3]}'
            if locus!=cur:
                if data: flush(cur,data)
                cur=locus; data={}
            data[fk]=seq
        if data: flush(cur,data)

tot=len(clip)
# recurrence
c=Counter(s for s in clip if len(s)>=20)
n_unique=len(c)
n_recur_seqs=sum(1 for s,n in c.items() if n>=5)
n_recur_instances=sum(n for s,n in c.items() if n>=5)
# homopolymer classification
def maxhomo(s):
    b=cur=1
    for i in range(1,len(s)):
        if s[i]==s[i-1]: cur+=1; b=max(b,cur)
        else: cur=1
    return b if s else 0
polyAT=sum(1 for s in clip if s and maxhomo(s)>=12 and (s.count('A')+s.count('T'))/len(s)>0.7)
# repeat-derived among recurrent (>=5) non-homopolymer
rep=Counter()
for s,n in c.items():
    if n<5: continue
    if maxhomo(s)/len(s)>=0.5: continue
    r=looks_repeat(s)
    if r: rep[r]+=n

print(f"Sampled {len(sample)} files; clipped seqs={tot}")
print(f"clipped >=20bp: {sum(1 for s in clip if len(s)>=20)}  unique={n_unique}")
print(f"recurrent seqs (identical, appearing >=5x): {n_recur_seqs} distinct, {n_recur_instances} instances "
      f"({100*n_recur_instances/max(1,tot):.1f}% of all clips)")
print(f"poly-A/T dominated clips (maxhomo>=12 & >70% A/T): {polyAT} ({100*polyAT/max(1,tot):.1f}%)")
print(f"repeat-consensus instances among recurrent non-homopolymer clips: {dict(rep)}")
print("\n-- palindrome self-fold examples (clip[:14] found in revcomp of same-side ALIGNED flank) --")
for locus,side,u,aln in palindrome_examples:
    print(f"  {locus} {side}\n     clip   {u}\n     flankRC{revcomp(aln)}")
print("\n-- poly-G tail examples --")
for locus,fk,u in polyg_examples:
    print(f"  {locus} {fk}: {u}")
