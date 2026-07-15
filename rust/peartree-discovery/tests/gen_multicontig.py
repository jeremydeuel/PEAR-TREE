import pysam, sys
src = pysam.AlignmentFile('test_data/test.bam')
# new header: several contigs incl MT and a >5-char name (both should be skipped for breakpoints)
targets = ['1','2','7','MT','chrLongName']   # order deliberately non-sorted-friendly
newhdr = pysam.AlignmentHeader.from_dict({
    'HD': {'VN':'1.6','SO':'coordinate'},
    'SQ': [{'SN':t,'LN':300_000_000} for t in targets],
})
tid = {t:i for i,t in enumerate(targets)}
reads = list(src)
out_records = []
# place a full copy of the read set onto each of contigs 1,2,7,MT,chrLongName
for cp in targets:
    for r in reads:
        a = pysam.AlignedSegment(newhdr)
        a.query_name = f"{r.query_name}_{cp}"
        a.query_sequence = r.query_sequence
        a.query_qualities = r.query_qualities
        a.flag = r.flag
        a.reference_id = tid[cp]
        a.reference_start = r.reference_start
        a.mapping_quality = r.mapping_quality
        a.cigartuples = r.cigartuples
        if r.is_paired:
            a.next_reference_id = tid[cp]
            a.next_reference_start = r.next_reference_start
            a.template_length = r.template_length
        for t in ('SA',):
            if r.has_tag(t):
                a.set_tag(t, r.get_tag(t))
        out_records.append(a)
out_records.sort(key=lambda x:(x.reference_id, x.reference_start))
with pysam.AlignmentFile(f"{sys.argv[1]}", 'wb', header=newhdr) as out:
    for a in out_records:
        out.write(a)
print("wrote", len(out_records), "records across", targets)
