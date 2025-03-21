def revcomp(seq: str):
    seq = [b for b in seq]
    seq.reverse()
    for i in range(len(seq)):
        if seq[i] == "A":
            seq[i] = "T"
        elif seq[i] == "T":
            seq[i] = "A"
        elif seq[i] == "G":
            seq[i] = "C"
        elif seq[i] == "C":
            seq[i] = "G"
        elif seq[i] == 'a':
            seq[i] = 't'
        elif seq[i] == 't':
            seq[i] = 'a'
        elif seq[i] == 'g':
            seq[i] = 'c'
        elif seq[i] == 'c':
            seq[i] = 'g'
    return "".join(seq)