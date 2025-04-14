# PEAR-TREE - paired ends of aberrant retrotransposons in phylogenetic trees
#
# Copyright (C) 2025 Jeremy Deuel <jeremy.deuel@usz.ch>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


from typing import Dict, List
from sequence_checks import sequence_matching_score
from combine_insertions_insertion import Insertion, TYPE_LEFT_POLYA, TYPE_RIGHT_POLYA, TYPE_FULL_INFO
from collections import Counter

def intersect_insertions(insertions: List[Insertion]) -> List[Insertion]:
    full_insertions = {}
    polyA = {}
    # sort everything neatly by name and type
    full_insertion_assoc = {}
    for i in insertions:
        name = (i.reference_name, i.left_pos, i.right_pos)
        if i.type is TYPE_FULL_INFO:
            if name not in full_insertions.keys():
                full_insertions[name] = [i]
            else:
                full_insertions[name].append(i)
            left_polyA_name = name[0], None, name[2]
            right_polyA_name = name[0], name[1], None
            if left_polyA_name in full_insertion_assoc.keys():
                full_insertion_assoc[left_polyA_name].append(name)
            else:
                full_insertion_assoc[left_polyA_name] = [name]
            if right_polyA_name in full_insertion_assoc.keys():
                full_insertion_assoc[right_polyA_name].append(name)
            else:
                full_insertion_assoc[right_polyA_name] = [name]
        elif i.type is TYPE_RIGHT_POLYA or i.type is TYPE_LEFT_POLYA:
            if name not in polyA.keys():
                polyA[name] = [i]
            else:
                polyA[name].append(i)
        else:
            raise ValueError(f"Unknown Insertion Type: {i.type}")

    print(f"imported {len(insertions)}  insertions, divided up in {len(full_insertions)} unique full-information insertions and {len(polyA)}  polyA insertions")
    for name, insertions in full_insertions.items():
        if len(insertions) == 1:
            full_insertions[name] = insertions[0]
        else:
            #combining
            if (score := sequence_matching_score([s.right_clipped for s in insertions])) < 0.6:
                print(f"ignoring {name} with {len(insertions)} entries, failed right_clip alignment, score is {score} with sequences {','.join([str(s.right_clipped) for s in insertions])}")
                full_insertions[name] = None
                continue
            if (score := sequence_matching_score([s.left_clipped for s in insertions])) < 0.6:
                print(f"ignoring {name} with {len(insertions)} entries, failed left_clip alignment, score is {score} with sequences {','.join([str(s.left_clipped) for s in insertions])}")
                full_insertions[name] = None
                continue
            if (score := sequence_matching_score([s.left_aligned for s in insertions])) < 0.6:
                print(f"ignoring {name} with {len(insertions)} entries, failed left_aligned alignment, score is {score} with sequences {','.join([str(s.left_aligned) for s in insertions])}")
                full_insertions[name] = None
                continue
            if (score := sequence_matching_score([s.right_aligned for s in insertions])) < 0.6:
                print(f"ignoring {name} with {len(insertions)} entries, failed right_aligned alignment, score is {score} with sequences {','.join([str(s.right_aligned) for s in insertions])}")
                full_insertions[name] = None
                continue
            combined = None
            for h in insertions:
                if combined is None: combined = h
                else: combined += h
            full_insertions[name] = combined

    print(f"processed {len(full_insertions)} full insertions, now polyA insertions.")
    for name, hits in polyA.items():
        continue #dont process these for now. This need a heavy filter
        full_hit = None
        #find full hit
        if name in full_insertion_assoc.keys():
            full_hit = [full_insertions[fa] for fa in full_insertion_assoc[name] if full_insertions[fa] is not None]
            if len(full_hit)!=1:
                full_hit = None
        if full_hit is not None:
            hits = hits + full_hit #combine final hit if it exists.
        if len(hits)==1:
            full_insertions[name] = hits[0]
            continue
        side = hits[0].type
        if len(hits)>1:
            if side is not TYPE_RIGHT_POLYA and (score := sequence_matching_score([s.right_clipped for s in hits if s.right_clipped is not None])) < 0.6:
                print(
                    f"ignoring polyA {name} with {len(hits)} entries, failed right_clip alignment, score is {score} with sequences {','.join([str(s.right_clipped) for s in hits])}")
                continue
            if side is not TYPE_RIGHT_POLYA and (score := sequence_matching_score([s.right_aligned for s in hits if s.right_aligned is not None])) < 0.6:
                print(
                    f"ignoring polyA {name} with {len(hits)} entries, failed right_aligned alignment, score is {score} with sequences {','.join([str(s.right_aligned) for s in hits])}")
                continue
            if side is not TYPE_LEFT_POLYA and (score := sequence_matching_score([s.left_clipped for s in hits if s.left_clipped is not None])) < 0.6:
                print(
                    f"ignoring polyA {name} with {len(hits)} entries, failed left_clipped alignment, score is {score} with sequences {','.join([str(s.left_clipped) for s in hits])}")
                continue
            if side is not TYPE_LEFT_POLYA and (score := sequence_matching_score([s.left_aligned for s in hits if s.left_aligned is not None])) < 0.6:
                print(
                    f"ignoring polyA {name} with {len(hits)} entries, failed left_aligned alignment, score is {score} with sequences {','.join([str(s.left_aligned) for s in hits])}")
                continue
            combined = None
            for h in hits:
                if combined is None:
                    combined = h
                else:
                    combined += h
            if full_hit is not None:
                #print(f"achieved full info for insertion {combined.name} by integrating a total of {len(hits)-1} polyAs")
                full_insertions[full_hit[0].name] = combined
            else:
                full_insertions[combined.name] = combined
    print(f"found a total of {len(full_insertions)} insertions")
    return [i for i in full_insertions.values() if i is not None]