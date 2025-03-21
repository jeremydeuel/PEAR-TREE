from typing import Dict, List
from sequence_checks import sequence_matching_score
from combine_insertions_insertion import Insertion
from collections import Counter

def intersect_insertions(insertions: Dict[str, Dict[str, Insertion]]) -> List[Insertion]:
    names = []
    for file_insertions in insertions.values():
        names += file_insertions.keys()
    print(f"found a total of {len(names)} unconsolidated insertions, reducing...")
    names = Counter(names)
    print(f"reduced to {len(names)} unconsolidated insertions, continuing...")
    single = []
    multi = []
    n_excluded = 0
    for name, n in names.items():
        if n == 1:
            for file_insertions in insertions.values():
                if name in file_insertions.keys():
                    single.append(file_insertions[name])
                    continue
        else:
            insertions_with_name = []
            for file_insertions in insertions.values():
                if name in file_insertions.keys():
                    insertions_with_name.append(file_insertions[name])
            if sequence_matching_score([s.right_clipped for s in insertions_with_name]) < 0.6 or sequence_matching_score([s.left_clipped for s in insertions_with_name]) < 0.6:
                #print(f"Excluding {name}")
                #print('L: '+'\nL: '.join([s.left_clipped for s in insertions_with_name]))
                #print("R: "+'\nR: '.join([s.right_clipped for s in insertions_with_name]))
                #print('')
                n_excluded += 1
            else:
                new_multi = insertions_with_name[0]
                for i in insertions_with_name[1:]:
                    new_multi += i
                multi.append(new_multi)
    print(f"found a total of {len(single)+len(multi)+n_excluded} insertions, thereof {len(single)} a single time and {len(multi)} multiple times, {n_excluded} where excluded due to different sequences in different files.")
    return single + multi