#!/usr/bin/env python3
import gzip, glob, os, re, sys, random
from collections import Counter

DIR = "/Users/jeremy/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents - STEM_Green_Lab/HNRNPA1/JAK2_HNRNPA1_paper/previous analyses bin/trees_mpn/discovery"
files = sorted(glob.glob(os.path.join(DIR, "*.discovery.fq.gz")))
random.seed(42)
sample = random.sample(files, min(40, len(files)))

def revcomp(s):
    return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

def homopolymer_frac(s):
    if not s: return 0
    # longest homopolymer run
    best=cur=1
    for i in range(1,len(s)):
        if s[i]==s[i-1]: cur+=1; best=max(best,cur)
        else: cur=1
    return best/len(s)

ADAPTERS = ['AGATCGGAAGAGC','GCTCTTCCGATCT']

clip_seqs = []            # all CLIPPED / CLIPPED_POLYA sequences
records = {}              # locus -> {field: seq}
n_records = 0
n_loci = 0
polyA_loci = 0
total_clip = 0
poly_g_tail = 0
poly_g_any = 0
adapter_ct = 0
low_complex = 0          # <=2 distinct bases
homop_dom = 0            # longest homopolymer >=50% of clip
short_clip = 0
palindrome_selffold = 0  # clip is revcomp-ish of aligned flank
gc_clip = []

def flush(locus, data):
    global n_loci, polyA_loci, total_clip, poly_g_tail, poly_g_any, adapter_ct
    global low_complex, homop_dom, short_clip, palindrome_selffold
    if not data: return
    n_loci += 1
    if 'polyA' in locus: polyA_loci += 1
    for field, seq in data.items():
        if 'CLIPPED' not in field: continue
        total_clip += 1
        clip_seqs.append(seq)
        u = seq.upper()
        if len(u) < 15: short_clip += 1
        if u.endswith('GGGGGG'): poly_g_tail += 1
        if 'GGGGGGGGG' in u: poly_g_any += 1
        for a in ADAPTERS:
            if a in u or revcomp(a) in u:
                adapter_ct += 1; break
        if len(set(u.replace('N','')))<=2: low_complex += 1
        if len(u) and homopolymer_frac(u) >= 0.5: homop_dom += 1
        if u: gc_clip.append((u.count('G')+u.count('C'))/len(u))
        # palindrome self-fold: clipped matches revcomp of the aligned flank of same side
        side = field.split(':')[0]  # LEFT/RIGHT
        aln = data.get(f'{side}:ALIGNED')
        if aln and len(u)>=12:
            probe = u[:12]
            if probe in revcomp(aln.upper()):
                palindrome_selffold += 1

for fp in sample:
    try:
        with gzip.open(fp,'rt') as fh:
            cur_locus=None; data={}
            while True:
                h=fh.readline()
                if not h: break
                h=h.strip()
                if not h.startswith('@') and not h.startswith('>'): continue
                seq=fh.readline().strip(); plus=fh.readline(); qual=fh.readline().strip()
                n_records_local=1
                # header: @seqname:R-L:SIDE:FIELD
                parts=h[1:].split(':')
                if len(parts)<4: continue
                seqname, pos, side, field = parts[0], parts[1], parts[2], parts[3]
                locus=f'{seqname}:{pos}'
                fkey=f'{side}:{field}'
                if locus!=cur_locus:
                    flush(cur_locus, data)
                    cur_locus=locus; data={}
                data[fkey]=seq
                global_n=0
            flush(cur_locus, data)
    except Exception as e:
        print(f"ERR {os.path.basename(fp)}: {e}", file=sys.stderr)

print(f"Sampled {len(sample)} of {len(files)} files")
print(f"Loci: {n_loci}   polyA-rescued loci: {polyA_loci} ({100*polyA_loci/max(1,n_loci):.1f}%)")
print(f"Clipped sequences examined: {total_clip}")
print(f"--- artefact-pattern prevalence among clipped seqs ---")
def pct(x): return f"{x} ({100*x/max(1,total_clip):.1f}%)"
print(f"poly-G tail (ends GGGGGG):     {pct(poly_g_tail)}")
print(f"poly-G run >=9 anywhere:       {pct(poly_g_any)}")
print(f"adapter substring present:     {pct(adapter_ct)}")
print(f"low-complexity (<=2 bases):    {pct(low_complex)}")
print(f"homopolymer-dominated (>=50%): {pct(homop_dom)}")
print(f"short clip (<15 bp):           {pct(short_clip)}")
print(f"palindrome self-fold vs flank: {pct(palindrome_selffold)}")
import statistics
if gc_clip:
    print(f"clip GC content mean={statistics.mean(gc_clip):.2f} median={statistics.median(gc_clip):.2f}")

# recurrence: most common exact clipped sequences (>=20bp) across the sample
c = Counter(s.upper() for s in clip_seqs if len(s)>=20)
print("\n--- top 25 recurrent clipped sequences (>=20bp), a hallmark of recurrent artefacts ---")
for seq,n in c.most_common(25):
    tag=[]
    if 'GGGGGGGGG' in seq: tag.append('polyG')
    if seq.count('A')/len(seq)>0.7 or seq.count('T')/len(seq)>0.7: tag.append('polyA/T')
    for a in ADAPTERS:
        if a in seq or revcomp(a) in seq: tag.append('adapter')
    if len(set(seq))<=2: tag.append('lowcomplex')
    print(f"{n:5d}  {seq[:60]}  {'/'.join(sorted(set(tag)))}")
