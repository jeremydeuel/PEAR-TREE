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
from combine_insertions_insertion import (
    Insertion, TYPE_LEFT_POLYA, TYPE_RIGHT_POLYA, TYPE_FULL_INFO,
    TYPE_LEFT_DISC, TYPE_RIGHT_DISC,
)
from collections import Counter

def _polya_to_one_sided(i: Insertion) -> Insertion:
    """TPRT: turn a legacy Bp+polyA record (`contig:L-polyA_P` / `contig:polyA_P-R`) into a
    one-sided locus on its real junction: the poly-A end keeps its coordinate (P, from the
    poly-A read's mate) but no sequence, so combine remaps / filters / writes only the real
    side, exactly like a discovery `oneside_` locus. The name is unchanged (sidecar link)."""
    pos = i.name.rsplit(':', 1)[1]
    start, end = pos.split('-')
    if i.type is TYPE_RIGHT_POLYA:
        i.right_clipped, i.right_aligned = None, None
        i.right_pos = int(end[len('polyA_'):])
        i.type, i.open_side = TYPE_RIGHT_DISC, 'RIGHT'
    else:
        i.left_clipped, i.left_aligned = None, None
        i.left_pos = int(start[len('polyA_'):])
        i.type, i.open_side = TYPE_LEFT_DISC, 'LEFT'
    return i


def _real_clip(i: Insertion):
    return i.right_clipped if i.type is TYPE_LEFT_DISC else i.left_clipped


def intersect_insertions(insertions: List[Insertion], keep_polya_one_sided: bool = None) -> List[Insertion]:
    """Pool per-sample discovery records by locus name. Full (two-sided) records are merged
    after a clip/flank agreement check. One-sided records — Feature-A discordant calls and
    TPRT `oneside_` loci — are kept: one representative per name (longest real-side clip);
    for `oneside_` loci the mates and source files of every sample are pooled onto it.
    Legacy Bp+polyA records stay parked unless `keep_polya_one_sided` (default:
    CONFIG['combine_insertions']['keep_polya_one_sided'], False), which keeps them as
    one-sided loci on their real junction."""
    if keep_polya_one_sided is None:
        try:
            from config import CONFIG
            keep_polya_one_sided = bool(CONFIG.get('combine_insertions', {}).get('keep_polya_one_sided', False))
        except ImportError:
            keep_polya_one_sided = False
    full_insertions = {}
    polyA = {}
    disc = {}  # Feature A: discordant-anchored one-sided junctions
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
        elif i.type is TYPE_RIGHT_DISC or i.type is TYPE_LEFT_DISC:
            # Feature A: parked like polyA until a dedicated discordant filter exists.
            disc.setdefault(name, []).append(i)
        else:
            raise ValueError(f"Unknown Insertion Type: {i.type}")

    print(f"imported {len(insertions)}  insertions, divided up in {len(full_insertions)} unique full-information insertions, {len(polyA)}  polyA insertions and {len(disc)} discordant-anchored insertions")
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
    if keep_polya_one_sided:
        # TPRT: one-sided loci on the real junction (genome-aware filters still apply)
        for name, hits in polyA.items():
            for h in hits:
                disc.setdefault(name, []).append(_polya_to_one_sided(h))
        polyA = {}
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
    # Feature A: surface discordant-anchored one-sided calls (parked above). Each is a
    # singleton within a discovery file; keep the longest-consensus representative so the
    # downstream clean-remap / clipped-remap filters can validate the one real side.
    n_one_sided = 0
    for name, hits in disc.items():
        rep = hits[0]
        for h in hits[1:]:
            real_side = h.right_clipped if h.type is TYPE_LEFT_DISC else h.left_clipped
            rep_side = rep.right_clipped if rep.type is TYPE_LEFT_DISC else rep.left_clipped
            if real_side is not None and (rep_side is None or len(real_side) > len(rep_side)):
                rep = h
        if getattr(rep, 'open_side', None) is not None:
            # TPRT one-sided: pool every sample's mates and files onto the representative
            # (and the longest aligned flank of the real side).
            n_one_sided += 1
            for h in hits:
                if h is rep:
                    continue
                rep.left_mates = rep.left_mates + h.left_mates
                rep.right_mates = rep.right_mates + h.right_mates
                rep.files = rep.files + [f for f in h.files if f not in rep.files]
                if rep.type is TYPE_RIGHT_DISC and h.left_aligned is not None and len(h.left_aligned) > len(rep.left_aligned):
                    rep.left_aligned = h.left_aligned
                if rep.type is TYPE_LEFT_DISC and h.right_aligned is not None and len(h.right_aligned) > len(rep.right_aligned):
                    rep.right_aligned = h.right_aligned
        full_insertions[rep.name] = rep
    if n_one_sided:
        print(f"kept {n_one_sided} one-sided (oneside_/poly-A) loci")
    print(f"found a total of {len(full_insertions)} insertions")
    return [i for i in full_insertions.values() if i is not None]