import gzip
peartree_file = 'MD7180s.peartree.fq.gz'
searchstrings_file = 'searchstrings.txt'


if __name__ == '__main__':
    searchstrings = set()
    with open(searchstrings_file, 'r') as f:
        for l in f:
            l = l.strip()
            if len(l):
                seqname, pos = l.split(":")
                start, end = pos.split("-")
                searchstrings.add((seqname, int(end), int(start)))
    ident = set()
    with gzip.open(peartree_file, 'rt') as f:
        for l in f:
            l = l.strip()
            if len(l):
                if l[0] == '@':
                    seqname, pos, side, clip = l[1:].split(":")
                    start, end = pos.split("-")
                    try:
                        start = int(start)
                    except:
                        start = None
                        pass
                    try:
                        end = int(end)
                    except:
                        end = None
                        pass
                    ident.add((seqname, start, end))
    print(f"imported {len(ident)} insertions.")
    searchstrings = set()
    with open(searchstrings_file, 'r') as f:
        for l in f:
            l = l.strip()
            if len(l):
                seqname, pos = l.split(":")
                start, end = pos.split("-")
                searchstrings.add((seqname, int(end), int(start)))
                if (seqname, int(start), int(end)) in ident:
                    print(f"found {seqname}:{end}-{start}")
                elif (seqname, None, int(end)) in ident:
                    print(f"found {seqname}:{end}-{start} as right-polyA")
                elif (seqname,  int(start), None) in ident:
                    print(f"found {seqname}:{end}-{start} as left-polyA")
                else:
                    print(f"did not find {seqname}:{end}-{start}")

