#!/usr/bin/env python3
import gzip, glob, os, sys
from collections import Counter, defaultdict

ROOT = "/Users/jeremy/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents - STEM_Green_Lab/HNRNPA1/spar_2ndrev/human/pear-tree/data/discovery"
indiv_dirs = sorted([d for d in glob.glob(os.path.join(ROOT,'*')) if os.path.isdir(d) and os.path.basename(d)!='genotyping'])

def revcomp(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def maxhomo(s):
    if not s: return 0
    b=cur=1
    for i in range(1,len(s)):
        if s[i]==s[i-1]: cur+=1; b=b if b>=cur else cur
        else: cur=1
    return b
ADAPTERS=['AGATCGGAAGAGC','GCTCTTCCGATCT']
ALU="GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGG"
ALU_K=[ALU[k:k+15] for k in range(0,len(ALU)-15,6)]

total_clip=0; n_loci=0; polyA_loci=0
low_complex=homop_dom=polyAT=poly_g_any=adapter_ct=short_clip=palindrome=0
clip_counter=Counter()
clip_indiv=defaultdict(int)        # clip>=20bp (>=2x within indiv) -> #individuals
indiv_clipcount=Counter(); indiv_filecount=Counter()

def process_locus(locus,data,local):
    global total_clip,n_loci,polyA_loci,low_complex,homop_dom,polyAT,poly_g_any
    global adapter_ct,short_clip,palindrome
    if not data: return
    n_loci+=1
    if 'polyA' in locus: polyA_loci+=1
    for fk,seq in data.items():
        if 'CLIPPED' not in fk: continue
        total_clip+=1
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
        if L>=20:
            clip_counter[u]+=1; local[u]+=1
        side=fk.split(':')[0]; aln=data.get(f'{side}:ALIGNED')
        if aln and L>=14 and u[:14] in revcomp(aln.upper()): palindrome+=1

for d in indiv_dirs:
    iv=os.path.basename(d)
    files=sorted(glob.glob(os.path.join(d,'*.txt.gz')))
    indiv_filecount[iv]=len(files)
    local=Counter()
    before=total_clip
    for fp in files:
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
                        process_locus(floc,fdata,local); floc=locus; fdata={}
                    fdata[fk]=seq
                process_locus(floc,fdata,local)
        except Exception as e:
            print(f"ERR {os.path.basename(fp)}: {e}", file=sys.stderr)
    indiv_clipcount[iv]=total_clip-before
    for s,n in local.items():
        if n>=2: clip_indiv[s]+=1
    print(f"...{iv}: {len(files)} files, {indiv_clipcount[iv]} clips", file=sys.stderr, flush=True)

def pct(x): return f"{x} ({100*x/max(1,total_clip):.2f}%)"
n_indiv=len(indiv_dirs)
print(f"=== HUMAN FULL SCAN (folder=individual): {n_indiv} individuals, {sum(indiv_filecount.values())} files ===")
for iv in sorted(indiv_filecount): print(f"   {iv}: {indiv_filecount[iv]} files, {indiv_clipcount[iv]} clips")
print(f"Loci: {n_loci}   polyA-rescued loci: {polyA_loci} ({100*polyA_loci/max(1,n_loci):.2f}%)")
print(f"Clipped sequences: {total_clip}")
print(f"poly-A/T dominated:             {pct(polyAT)}")
print(f"low-complexity (<=2 bases):     {pct(low_complex)}")
print(f"homopolymer-dominated (>=50%):  {pct(homop_dom)}")
print(f"poly-G run>=9:                  {pct(poly_g_any)}")
print(f"short clip (<15bp):             {pct(short_clip)}")
print(f"palindrome self-fold:           {pct(palindrome)}")
print(f"adapter substring:              {pct(adapter_ct)}")

clips20=sum(clip_counter.values()); uniq=len(clip_counter)
recur=[(s,n) for s,n in clip_counter.items() if n>=5]; recur_inst=sum(n for _,n in recur)
print(f"\nclips>=20bp: {clips20}  unique: {uniq}")
print(f"recurrent (>=5x): {len(recur)} distinct, {recur_inst} instances ({100*recur_inst/max(1,clips20):.2f}% of >=20bp)")
# cross-individual framing (the germline logic)
multi = sum(1 for s,k in clip_indiv.items() if k>=int(0.5*n_indiv))
allind= sum(1 for s,k in clip_indiv.items() if k==n_indiv)
one   = sum(1 for s,k in clip_indiv.items() if k==1)
tot_ci= len(clip_indiv)
print(f"clips (>=2x within an indiv), tracked: {tot_ci}")
print(f"  present in ALL {n_indiv} individuals: {allind}   (reference-repeat / recurrent artefact)")
print(f"  present in >=50% individuals:        {multi}")
print(f"  present in exactly ONE individual:   {one}  ({100*one/max(1,tot_ci):.1f}%)  (candidate germline/private)")
# Alu among cross-individual recurrent
alu_multi=sum(1 for s,k in clip_indiv.items() if k>=int(0.5*n_indiv) and any(km in s for km in ALU_K))
print(f"  of the >=50%-individual clips, Alu-consensus: {alu_multi}")

print("\n--- top 35 recurrent clipped sequences (>=20bp) with #individuals present ---")
for s,n in clip_counter.most_common(35):
    tag=[]
    if 'GGGGGGGGG' in s: tag.append('polyG')
    if (s.count('A')+s.count('T'))/len(s)>0.7: tag.append('polyA/T')
    if len(set(s))<=2: tag.append('lowcplx')
    if any(km in s for km in ALU_K): tag.append('Alu')
    ki=clip_indiv.get(s,0)
    print(f"{n:6d} [ind={ki:2d}] {s[:56]:56s} {'/'.join(dict.fromkeys(tag))}")
