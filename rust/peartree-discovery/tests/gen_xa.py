import pysam, sys, random
random.seed(7)
def rnd(n): return ''.join(random.choice('ACGT') for _ in range(n))
hdr=pysam.AlignmentHeader.from_dict({'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'1','LN':300_000_000}]})
recs=[]
def mk(name,seq,flag,start,cigar,mapq,xa=None,mref=None,mstart=0):
    a=pysam.AlignedSegment(hdr); a.query_name=name; a.query_sequence=seq
    a.query_qualities=pysam.qualitystring_to_array('I'*len(seq)); a.flag=flag
    a.reference_id=0; a.reference_start=start; a.mapping_quality=mapq; a.cigarstring=cigar
    if mref is not None: a.next_reference_id=0; a.next_reference_start=mstart
    if xa is not None: a.set_tag('XA', xa)
    recs.append(a)

def insertion(pos, tag, xa=None):
    # left breakpoint at pos (forward reads, 60S91M) + right breakpoint at pos+8 (reverse, 90M61S)
    CL=rnd(60); RUL=rnd(91); RUR=rnd(90); E5=rnd(61)
    for k in range(3):
        mk(f"{tag}_L{k}", CL+RUL, 0x1|0x40, pos, "60S91M", 60, xa=xa)
    for k in range(3):
        mk(f"{tag}_R{k}", RUR+E5, 0x1|0x80|0x10, pos+8-90, "90M61S", 60, xa=xa)

# insertion A at 1,000,000: clean (no XA) -> should be called
insertion(1_000_000, "clean", xa=None)
# insertion B at 2,000,000: every supporting read carries a full-length XA -> should be filtered
insertion(2_000_000, "fullmap", xa="1,+123456,151M,0;")

recs.sort(key=lambda x:(x.reference_id,x.reference_start))
with pysam.AlignmentFile(sys.argv[1],'wb',header=hdr) as o:
    for a in recs: o.write(a)
print("wrote", len(recs))
