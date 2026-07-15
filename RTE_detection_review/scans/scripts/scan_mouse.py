#!/usr/bin/env python3
import gzip, glob, os, sys, re
from collections import Counter, defaultdict

DIR = "/Users/jeremy/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents - STEM_Green_Lab/HNRNPA1/spar_2ndrev/mouse/pear-tree/data/discovery"
files = sorted(glob.glob(os.path.join(DIR, "*.txt.gz")))

def revcomp(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def maxhomo(s):
    if not s: return 0
    b=cur=1
    for i in range(1,len(s)):
        if s[i]==s[i-1]: cur+=1; b=b if b>=cur else cur
        else: cur=1
    return b
def indiv(fn):
    m=re.match(r'(M[DX]\d+)', os.path.basename(fn))
    return m.group(1) if m else os.path.basename(fn)[:7]
ADAPTERS=['AGATCGGAAGAGC','GCTCTTCCGATCT']

# tallies
total_clip=0; n_loci=0; polyA_loci=0
low_complex=homop_dom=polyAT=poly_g_any=adapter_ct=short_clip=palindrome=has_polyA_run=0
clip_counter=Counter()          # exact clips >=20 bp, global recurrence
clip_indiv=defaultdict(set)      # clip>=20bp -> set of individuals (only for clips seen >=2x within an indiv)
indiv_files=Counter(); indiv_clip=Counter()
per_indiv_local=defaultdict(Counter)  # reset per file group? no—accumulate per individual, flush at indiv change

def process_locus(locus,data,cur_indiv):
    global total_clip,n_loci,polyA_loci,low_complex,homop_dom,polyAT,poly_g_any
    global adapter_ct,short_clip,palindrome,has_polyA_run
    if not data: return
    n_loci+=1
    if 'polyA' in locus: polyA_loci+=1
    for fk,seq in data.items():
        if 'CLIPPED' not in fk: continue
        total_clip+=1; indiv_clip[cur_indiv]+=1
        u=seq.upper(); L=len(u)
        if L<15: short_clip+=1
        if 'GGGGGGGGG' in u: poly_g_any+=1
        for a in ADAPTERS:
            if a in u or revcomp(a) in u: adapter_ct+=1; break
        nb=set(u); nb.discard('N')
        if len(nb)<=2: low_complex+=1
        if L:
            mh=maxhomo(u)
            if mh/L>=0.5: homop_dom+=1
            if mh>=12 and (u.count('A')+u.count('T'))/L>0.7: polyAT+=1
        if ('A'*12 in u) or ('T'*12 in u): has_polyA_run+=1
        if L>=20:
            clip_counter[u]+=1
            per_indiv_local[cur_indiv][u]+=1
        side=fk.split(':')[0]; aln=data.get(f'{side}:ALIGNED')
        if aln and L>=14 and u[:14] in revcomp(aln.upper()): palindrome+=1

cur=None; data={}; cur_indiv=None
for fp in files:
    fi=indiv(fp); indiv_files[fi]+=1
    if cur_indiv is None: cur_indiv=fi
    try:
        with gzip.open(fp,'rt') as fh:
            floc=None; fdata={}
            while True:
                h=fh.readline()
                if not h: break
                h=h.strip()
                if not (h.startswith('@') or h.startswith('>')): continue
                seq=fh.readline().rstrip('\n'); fh.readline(); fh.readline()
                p=h[1:].split(':')
                if len(p)<4: continue
                locus=f'{p[0]}:{p[1]}'; fk=f'{p[2]}:{p[3]}'
                if locus!=floc:
                    process_locus(floc,fdata,fi); floc=locus; fdata={}
                fdata[fk]=seq
            process_locus(floc,fdata,fi)
    except Exception as e:
        print(f"ERR {os.path.basename(fp)}: {e}", file=sys.stderr)

# build clip_indiv from per_indiv_local: clip present >=2x within an individual -> that individual
for iv,cc in per_indiv_local.items():
    for s,n in cc.items():
        if n>=2: clip_indiv[s].add(iv)
per_indiv_local.clear()

def pct(x): return f"{x} ({100*x/max(1,total_clip):.2f}%)"
n_indiv=len(indiv_files)
print(f"=== MOUSE FULL SCAN: {len(files)} files, {n_indiv} individuals ===")
print(f"Loci: {n_loci}   polyA-rescued loci: {polyA_loci} ({100*polyA_loci/max(1,n_loci):.2f}%)")
print(f"Clipped sequences: {total_clip}")
print(f"clip contains polyA/T run>=12:  {pct(has_polyA_run)}   <-- L1Md/SINE-like; ERV/IAP should LACK this")
print(f"poly-A/T dominated:             {pct(polyAT)}")
print(f"low-complexity (<=2 bases):     {pct(low_complex)}")
print(f"homopolymer-dominated (>=50%):  {pct(homop_dom)}")
print(f"poly-G run>=9:                  {pct(poly_g_any)}")
print(f"short clip (<15bp):             {pct(short_clip)}")
print(f"palindrome self-fold:           {pct(palindrome)}")
print(f"adapter substring:              {pct(adapter_ct)}")

clips20=sum(clip_counter.values()); uniq=len(clip_counter)
recur=[(s,n) for s,n in clip_counter.items() if n>=5]
recur_inst=sum(n for _,n in recur)
print(f"\nclips>=20bp: {clips20}  unique: {uniq}")
print(f"recurrent (>=5x): {len(recur)} distinct, {recur_inst} instances ({100*recur_inst/max(1,clips20):.2f}% of >=20bp)")

# cross-individual recurrence: clips present in many individuals (>=2x within each) = reference-repeat/artefact
cross = Counter({s:len(ivs) for s,ivs in clip_indiv.items()})
multi_indiv = sum(1 for s,k in cross.items() if k>=int(0.5*n_indiv))
print(f"clips (>=2x within an indiv) present in >=50% of individuals: {multi_indiv}  (reference-repeat / recurrent artefact)")
one_indiv = sum(1 for s,k in cross.items() if k==1)
print(f"clips present in exactly ONE individual (>=2x): {one_indiv}  (candidate germline / private)")

print("\n--- top 40 recurrent clipped sequences (>=20bp) with #individuals ---")
for s,n in clip_counter.most_common(40):
    tag=[]
    if 'GGGGGGGGG' in s: tag.append('polyG')
    if (s.count('A')+s.count('T'))/len(s)>0.7: tag.append('polyA/T')
    if len(set(s))<=2: tag.append('lowcplx')
    if ('A'*12 in s) or ('T'*12 in s): tag.append('hasPolyA')
    ki=len(clip_indiv.get(s,()))
    print(f"{n:6d} [ind={ki:2d}] {s[:56]:56s} {'/'.join(dict.fromkeys(tag))}")
