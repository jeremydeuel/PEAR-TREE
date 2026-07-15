#!/usr/bin/env python3
import gzip, glob, os, sys, statistics
from collections import Counter

DIR = "/Users/jeremy/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents - STEM_Green_Lab/HNRNPA1/JAK2_HNRNPA1_paper/previous analyses bin/trees_mpn/discovery"
files = sorted(glob.glob(os.path.join(DIR, "*.discovery.fq.gz")))

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
L1 ="GGAGGAGCCAAGATGGCCGAATAGGAACAGCTCCGGTCTACAGCTCCCAGCGTGAGCGA"
ALU_K=[ALU[k:k+15] for k in range(0,len(ALU)-15,6)]
L1_K =[L1[k:k+15]  for k in range(0,len(L1)-15,6)]
def repeat_kind(s):
    for k in ALU_K:
        if k in s: return 'Alu'
    for k in L1_K:
        if k in s: return 'L1'
    return None

n_files=0; empty_files=0
n_loci=0; polyA_loci=0
total_clip=0
poly_g_tail=poly_g_any=adapter_ct=low_complex=homop_dom=short_clip=palindrome=polyAT=0
gc_sum=0.0; gc_n=0
clip_counter=Counter()   # exact clips >=20bp (for recurrence)
side_field=Counter()

def process_locus(locus, data):
    global n_loci, polyA_loci, total_clip, poly_g_tail, poly_g_any, adapter_ct
    global low_complex, homop_dom, short_clip, palindrome, polyAT, gc_sum, gc_n
    if not data: return
    n_loci+=1
    if 'polyA' in locus: polyA_loci+=1
    for fk,seq in data.items():
        if 'CLIPPED' not in fk: continue
        total_clip+=1
        u=seq.upper()
        L=len(u)
        if L<15: short_clip+=1
        if u.endswith('GGGGGG'): poly_g_tail+=1
        if 'GGGGGGGGG' in u: poly_g_any+=1
        for a in ADAPTERS:
            if a in u or revcomp(a) in u: adapter_ct+=1; break
        nb=set(u); nb.discard('N')
        if len(nb)<=2: low_complex+=1
        if L:
            mh=maxhomo(u)
            if mh/L>=0.5: homop_dom+=1
            if mh>=12 and (u.count('A')+u.count('T'))/L>0.7: polyAT+=1
            gc_sum+=(u.count('G')+u.count('C'))/L; gc_n+=1
        if L>=20: clip_counter[u]+=1
        # palindrome self-fold vs same-side aligned flank
        side=fk.split(':')[0]
        aln=data.get(f'{side}:ALIGNED')
        if aln and L>=14 and u[:14] in revcomp(aln.upper()):
            palindrome+=1

for fp in files:
    n_files+=1
    try:
        with gzip.open(fp,'rt') as fh:
            cur=None; data={}
            while True:
                h=fh.readline()
                if not h: break
                h=h.strip()
                if not (h.startswith('@') or h.startswith('>')): continue
                seq=fh.readline().rstrip('\n'); fh.readline(); fh.readline()
                p=h[1:].split(':')
                if len(p)<4: continue
                locus=f'{p[0]}:{p[1]}'; fk=f'{p[2]}:{p[3]}'
                if locus!=cur:
                    process_locus(cur,data); cur=locus; data={}
                data[fk]=seq
            process_locus(cur,data)
    except Exception as e:
        print(f"ERR {os.path.basename(fp)}: {e}", file=sys.stderr)
    if n_files % 50 == 0:
        print(f"...{n_files}/{len(files)} files, {total_clip} clips so far", file=sys.stderr, flush=True)

def pct(x): return f"{x} ({100*x/max(1,total_clip):.2f}%)"
print(f"=== FULL SCAN: {n_files} files ===")
print(f"Loci: {n_loci}   polyA-rescued loci: {polyA_loci} ({100*polyA_loci/max(1,n_loci):.2f}%)")
print(f"Clipped sequences: {total_clip}")
print(f"poly-G tail (ends GGGGGG):      {pct(poly_g_tail)}")
print(f"poly-G run >=9 anywhere:        {pct(poly_g_any)}")
print(f"adapter substring:             {pct(adapter_ct)}")
print(f"low-complexity (<=2 bases):     {pct(low_complex)}")
print(f"homopolymer-dominated (>=50%):  {pct(homop_dom)}")
print(f"poly-A/T dominated (>=12,>70%): {pct(polyAT)}")
print(f"short clip (<15bp):             {pct(short_clip)}")
print(f"palindrome self-fold:           {pct(palindrome)}")
print(f"clip GC mean: {gc_sum/max(1,gc_n):.3f}")

# recurrence
clips20 = sum(v for v in clip_counter.values())
uniq = len(clip_counter)
recur_seqs = [(s,n) for s,n in clip_counter.items() if n>=5]
recur_inst = sum(n for _,n in recur_seqs)
print(f"\nclips>=20bp: {clips20}  unique: {uniq}")
print(f"recurrent (>=5x): {len(recur_seqs)} distinct seqs, {recur_inst} instances "
      f"({100*recur_inst/max(1,clips20):.2f}% of >=20bp clips, {100*recur_inst/max(1,total_clip):.2f}% of all clips)")
# repeat classification of recurrent non-homopolymer
rep=Counter()
for s,n in recur_seqs:
    if maxhomo(s)/len(s)>=0.5: continue
    k=repeat_kind(s)
    if k: rep[k]+=n
print(f"repeat-consensus among recurrent non-homopolymer clips: {dict(rep)}")
print("\n--- top 30 recurrent clipped sequences (>=20bp) ---")
for s,n in clip_counter.most_common(30):
    tag=[]
    if 'GGGGGGGGG' in s: tag.append('polyG')
    if (s.count('A')+s.count('T'))/len(s)>0.7: tag.append('polyA/T')
    if repeat_kind(s): tag.append(repeat_kind(s))
    if len(set(s))<=2: tag.append('lowcplx')
    print(f"{n:6d}  {s[:58]:58s} {'/'.join(dict.fromkeys(tag))}")
