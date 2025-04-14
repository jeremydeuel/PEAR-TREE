import gzip, os
peartree_file = 'MD7180s.peartree.fq.gz'
searchstrings_file = 'searchstrings.txt'


if __name__ == '__main__':
    found_stats = {}
    for peartree_file in os.listdir("."):
        if ".gz" in peartree_file:
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
                            try:
                                seqname, pos, side, clip = l[1:].split(":")
                            except Exception as e:
                                continue #probably a quality string
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
            searchstrings = set()
            found = 0
            not_found = 0
            with open(searchstrings_file, 'r') as f:
                for l in f:
                    l = l.strip()
                    if len(l):
                        if not l in found_stats.keys():
                            found_stats[l] = [0,0,0]
                        seqname, pos = l.split(":")
                        start, end = pos.split("-")
                        searchstrings.add((seqname, int(end), int(start)))
                        if (seqname, int(start), int(end)) in ident:
                            #print(f"  found {seqname}:{end}-{start}")
                            found += 1
                            found_stats[l][1] += 1
                        elif (seqname,  int(start), None) in ident:
                            print(f"  found {seqname}:{end}-{start} as right-polyA")
                            found += 1
                            found_stats[l][1] += 1
                            found_stats[l][2] += 1
                        elif (seqname,  None, int(end)) in ident:
                            print(f"  found {seqname}:{end}-{start} as left-polyA")
                            found += 1
                            found_stats[l][1] += 1
                        else:
                            #print(f"  did not find {seqname}:{end}-{start}")
                            not_found += 1
                            found_stats[l][0] += 1
            print(f"{peartree_file} found {100*found/(not_found+found):.0f}%, total {len(ident)} insertions")
    print(f"\nsensitivity stats:")
    for insertion, (not_found, found, polyA) in found_stats.items():
        print(f"{insertion}: {found}/{not_found+found} = {100*found/(not_found+found):.0f}%, thereof {100*polyA/found if found>0 else 0:.0f}% as polyA rescues.")
