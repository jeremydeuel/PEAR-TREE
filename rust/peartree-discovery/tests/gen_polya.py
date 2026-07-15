import pysam, sys, random
random.seed(42)
def rnd(n): return ''.join(random.choice('ACGT') for _ in range(n))
targets=['1']
hdr=pysam.AlignmentHeader.from_dict({'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':t,'LN':300_000_000} for t in targets]})
recs=[]
def mk(name,seq,qual_ch,flag,start,cigar,mapq,mref=None,mstart=0):
    a=pysam.AlignedSegment(hdr); a.query_name=name; a.query_sequence=seq
    a.query_qualities=pysam.qualitystring_to_array(qual_ch*len(seq)); a.flag=flag
    a.reference_id=0 if start is not None else -1
    a.reference_start=start if start is not None else (mstart or 0)
    a.mapping_quality=mapq; a.cigarstring=cigar
    if mref is not None:
        a.next_reference_id=0; a.next_reference_start=mstart
    recs.append(a)

# --- LEFT real breakpoint at L=900000: forward reads, cigar 60S91M ---
CL=rnd(60); RUL=rnd(91)
for k in range(3):
    mk(f"leftbp{k}", CL+RUL, 'I', 0x1|0x40, 900000, "60S91M", 60, mref=0, mstart=900500)
# --- RIGHT real breakpoint at R=1000000: reverse reads, cigar 90M61S ---
RUR=rnd(90); E5=rnd(61)
for k in range(3):
    mk(f"rightbp{k}", RUR+E5, 'I', 0x1|0x80|0x10, 999910, "90M61S", 60, mref=0, mstart=999400)
# --- polyA read1 (mapq0, fwd, not proper, mate mapped) + mate read2 (mapq60, reverse) ---
E3=rnd(40)
mk("pa", E3+"A"*25, 'I', 0x1|0x40, 999950, "65M", 0, mref=0, mstart=999860)   # low mapq -> polyA candidate
mk("pa", rnd(90), 'I', 0x1|0x80|0x10, 999860, "90M", 60, mref=0, mstart=999950) # mate: reference_end=999950
recs.sort(key=lambda x:(x.reference_id,x.reference_start))
with pysam.AlignmentFile(sys.argv[1],'wb',header=hdr) as o:
    for a in recs: o.write(a)
print("wrote",len(recs))
