from typing import Tuple, List, Iterable
from config import CONFIG


def find_consensus(seqs: Iterable[str]) -> Tuple[List[int], List[float], str]:
    """
    Find consensus sequence of an iterable of sequences, that must be left-aligned
    returns a tuple of
        - list of integers: number of sequences matching consensus base at each position
        - list of floats: Fraction of sequences matching consensus sequences at each position
        - string: Consensus sequence
    """
    lengths = [len(s) for s in seqs]
    seqs = [s.upper() for s in seqs] #ensure upper case
    consensus_seq = ""
    consensus_score = []
    consensus_n = []
    for position in range(max(lengths)):
        base_stat = {
            'A': 0, 'T': 0, 'G': 0, 'C': 0
        }
        for sequence in range(len(lengths)):
            if lengths[sequence] <= position: continue
            base = seqs[sequence][position]
            if not base in base_stat.keys():
                continue
            base_stat[base] += 1
        n_sequences = sum(base_stat.values())
        n_sorted = sorted(base_stat.items(), key=lambda x: x[1], reverse=True)
        delta_first_two_hits = n_sorted[0][1]-n_sorted[1][1]
        if (n_sequences < 3 and delta_first_two_hits == n_sequences) or (delta_first_two_hits > 2):
            consensus_score.append(delta_first_two_hits)
            consensus_n.append(n_sorted[0][1])
            consensus_seq += n_sorted[0][0]
        else:
            break #only extract seq to the first ambigous base.
            #consensus_score.append(0)
            #consensus_n.append(0)
            #consensus_seq += "N"
    return consensus_n, consensus_score, consensus_seq


def is_good_consensus(consensus_n: List[int], consensus_score: List[float]) -> bool:
    """
    Check if a consensus is "good", thus passes filtering criteria
    """
    for i in range(len(consensus_n)):
        if consensus_n[i] >= CONFIG['discovery']['min_evidence_reads_per_breakpoint'] and \
                consensus_score[i] >= CONFIG['discovery']['min_consensus_score_for_good_base']:
            pass
        else:
            break
        if i >= CONFIG['discovery']['min_good_bases']:
            return True
    return False
