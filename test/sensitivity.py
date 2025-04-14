import gzip

gt_file = "MD7180s.genotype.txt.gz"
true_pos_file = "true_positives.txt"

if __name__ == '__main__':
    true_pos_list = set()
    with open(true_pos_file, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line):
                s = line.split(":")
                if len(s)!= 2:
                    print(f"error: not recognised {line}")
                    continue
                reference_name, pos = line.split(":")
                end, start = pos.split("-")
                true_pos_list.add((reference_name, start, end))
    found = set()
    with gzip.open(gt_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if len(line):
                name, genotype, score_gt, score_alt = line.split("\t")
                if name == "insertion":
                    continue
                reference_name, pos = name.split(":")
                start, end = pos.split("-")
                if (reference_name, start, end) in true_pos_list:
                    print(f"found {reference_name}:{start}-{end} with genotype {genotype}")
                    found.add((reference_name, start, end))
    print(f"found {len(found)} of {len(true_pos_list)}")
    not_found = true_pos_list ^ found
    for reference_name, start, end in not_found:
        print(f" - {reference_name}:{end}-{start}")