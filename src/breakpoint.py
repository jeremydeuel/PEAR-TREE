
CLIP_RIGHT = 1
CLIP_LEFT = 2
from quality_seq import QualitySeq
from typing import Iterable, Tuple, List
from sequence_checks import clean_clipped_seq
from collections import Counter
from config import CONFIG
from consensus import find_consensus
class Breakpoint:
    def __init__(self, side: int, reference_name: str, breakpoint: int, query_name: str,
                 clipped: QualitySeq, unclipped: QualitySeq, is_read1: bool, is_forward: bool, exclude: bool):
        self.exclude: bool = exclude
        self.is_forward: bool = is_forward
        self.is_read1: bool = is_read1
        self.query_name: str = query_name
        self.breakpoint: int = breakpoint
        self.clipped: QualitySeq = clipped
        self.unclipped: QualitySeq = unclipped
        self.reference_name: str = reference_name
        self.side: int = side
        self.has_mate:bool = False
        self.bp_precise: bool = True #marker if breakpoint is precise or not (discordant reads)
        self.mates: List[Tuple[bool, str]] | None = []
        self.mate_seqs: List = []
        self.n_reads: int = 1

    def __str__(self):
        return f"{self.reference_name}:{self.breakpoint}:{'L' if self.side is CLIP_LEFT else 'R'}"


    stats = {
        CLIP_LEFT: {
            'too_few': 0,
            'rescued_pA': 0,
            'excluded': 0,
            'too_few_after_filter':0,
            'clipped_failed': 0,
            'unclipped_failed': 0,
            'polymer': 0,
            'passed': 0

        },
        CLIP_RIGHT: {
            'too_few': 0,
            'rescued_pA': 0,
            'excluded': 0,
            'too_few_after_filter':0,
            'clipped_failed': 0,
            'unclipped_failed': 0,
            'polymer': 0,
            'passed': 0

        },
    }
    @staticmethod
    def DEBUG_check_if_breakpoint_of_interest(breakpoint: 'Breakpoint') -> bool:
        if breakpoint.reference_name == '14' and abs(breakpoint.breakpoint-35038214)<40: return True
        if breakpoint.reference_name == '6' and abs(breakpoint.breakpoint-143309573)<40: return True
        if breakpoint.reference_name == '16' and abs(breakpoint.breakpoint-51440715)<40: return True
        if breakpoint.reference_name == '11' and abs(breakpoint.breakpoint-16580607)<40: return True
        if breakpoint.reference_name == '8' and abs(breakpoint.breakpoint-12146960)<40: return True
        if breakpoint.reference_name == 'X' and abs(breakpoint.breakpoint-57020400)<40: return True
        if breakpoint.reference_name == '9' and abs(breakpoint.breakpoint-70209274)<40: return True
        if breakpoint.reference_name == '5' and abs(breakpoint.breakpoint-61887543)<40: return True
        if breakpoint.reference_name == '18' and abs(breakpoint.breakpoint-83460413)<40: return True
        if breakpoint.reference_name == '7' and abs(breakpoint.breakpoint-40606189)<40: return True
        if breakpoint.reference_name == '7' and abs(breakpoint.breakpoint-62347358)<40: return True
        if breakpoint.reference_name == '3' and abs(breakpoint.breakpoint-108708290)<40: return True
        return False
    @staticmethod
    def join(breakpoints: Iterable['Breakpoint']) -> 'Breakpoint':
        #print(f"joining {', '.join([str(b) for b in breakpoints])}")
        if len(breakpoints)<2:
            # rescue solo polyA breakpoints
            if breakpoints[0].side is CLIP_LEFT and str(breakpoints[0].clipped[-8:]) == 'AAAAAAAA':
                bp = breakpoints[0]
                bp.clipped = clean_clipped_seq(bp.clipped.revcomp())
                if len(bp.clipped) < CONFIG['discovery']['min_clip_len']: return None
                if not bp.is_forward:
                    bp.has_mate = True
                    bp.mates = [(not bp.is_read1, bp.query_name)]
                if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]):
                    print(
                    f"rescued {'L' if breakpoints[0].side is CLIP_LEFT else 'R'} breakpoint {breakpoints[0].reference_name}:{breakpoints[0].breakpoint} with polyA: {str(breakpoints[0].clipped)}")
                Breakpoint.stats[breakpoints[0].side]['rescued_pA'] += 1
                return bp
            elif breakpoints[0].side is CLIP_RIGHT and str(breakpoints[0].clipped[:8]) == 'TTTTTTTT':
                bp = breakpoints[0]
                bp.clipped = clean_clipped_seq(bp.clipped)
                if len(bp.clipped) < CONFIG['discovery']['min_clip_len']: return None
                if bp.is_forward:
                    bp.has_mate = True
                    bp.mates = [(not bp.is_read1, bp.query_name)]
                if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]):
                    print(
                    f"rescued {'L' if breakpoints[0].side is CLIP_LEFT else 'R'} breakpoint {breakpoints[0].reference_name}:{breakpoints[0].breakpoint} with polyA: {str(breakpoints[0].clipped)}")
                Breakpoint.stats[breakpoints[0].side]['rescued_pA'] += 1
                return bp
            if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]): print(f"omitting {'L' if breakpoints[0].side is CLIP_LEFT else 'R'} breakpoint {breakpoints[0].reference_name}:{breakpoints[0].breakpoint}: only a single one.")
            Breakpoint.stats[breakpoints[0].side]['too_few'] += 1
            return None
        for bp in breakpoints:
            if bp.exclude:
                if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]): print(
                    f"omitting {'L' if breakpoints[0].side is CLIP_LEFT else 'R'} {breakpoints[0].breakpoint} exclusion flag set.")
                Breakpoint.stats[breakpoints[0].side]['excluded'] += 1
                return None #if any has the exclude flag set, "poison" all other breakpoints in that area. exclusion means that the clipped part maps somewhere else.
        side = breakpoints[0].side
        reference_name = breakpoints[0].reference_name
        # cleanup by removing adapters and marking if there is a mate
        bps = []
        for bp in breakpoints:
            assert side == bp.side
            assert reference_name == bp.reference_name
            if side is CLIP_LEFT:
                bp.clipped = clean_clipped_seq(bp.clipped.revcomp()) #left clipped is reverse complemented from now on!
            else:
                bp.clipped = clean_clipped_seq(bp.clipped)
            if side is CLIP_LEFT and not bp.is_forward:
                bp.has_mate = True
            if side is CLIP_RIGHT and bp.is_forward:
                bp.has_mate = True
            if bp.bp_precise and len(bp.clipped) >= CONFIG['discovery']['min_good_bases']:
                bps.append(bp.breakpoint)
        if not len(bps):
            Breakpoint.stats[side]['too_few_after_filter'] += 1
            if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]): print(f"removing {'L' if side is CLIP_LEFT else 'R'} breakpoint because no reads pass qc")
            return None#no bp found
        best_bp, n = Counter(bps).most_common()[0]
        clipped = []
        unclipped = []
        if n < CONFIG['discovery']['min_evidence_reads_per_breakpoint']:
            for bp in breakpoints:
                if bp.side is CLIP_LEFT and str(bp.clipped[-8:]) == 'AAAAAAAA':
                    bp.clipped = clean_clipped_seq(bp.clipped.revcomp())
                    if len(bp.clipped) < CONFIG['discovery']['min_clip_len']: return None
                    if not bp.is_forward:
                        bp.has_mate = True
                        bp.mates = [(not bp.is_read1, bp.query_name)]
                    if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]):
                        print(
                        f"rescued {'L' if breakpoints[0].side is CLIP_LEFT else 'R'} breakpoint {breakpoints[0].reference_name}:{breakpoints[0].breakpoint} with polyA: {str(breakpoints[0].clipped)}")
                    Breakpoint.stats[breakpoints[0].side]['rescued_pA'] += 1
                    return bp
                elif bp.side is CLIP_RIGHT and str(bp.clipped[:8]) == 'TTTTTTTT':
                    bp.clipped = clean_clipped_seq(bp.clipped)
                    if len(bp.clipped) < CONFIG['discovery']['min_clip_len']: return None
                    if bp.is_forward:
                        bp.has_mate = True
                        bp.mates = [(not bp.is_read1, bp.query_name)]
                    if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]):
                        print(
                        f"rescued {'L' if breakpoints[0].side is CLIP_LEFT else 'R'} breakpoint {breakpoints[0].reference_name}:{breakpoints[0].breakpoint} with polyA: {str(breakpoints[0].clipped)}")
                    Breakpoint.stats[breakpoints[0].side]['rescued_pA'] += 1
                    return bp
            Breakpoint.stats[side]['too_few_after_filter'] += 1
            if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]): print(f"removing {'L' if side is CLIP_LEFT else 'R'} breakpoint because bps don't converge: {Counter(bps)}")
            return None
        mates = []
        deltas = []
        for bp in breakpoints:
            if bp.has_mate:
                mates.append((not bp.is_read1, bp.query_name))
            if not bp.bp_precise: continue #disjoint mates, useless for exact bp determination
            delta_bp = bp.breakpoint - best_bp
            deltas.append(delta_bp)
            if side is CLIP_RIGHT:
                if delta_bp == 0:
                    clipped.append(bp.clipped)
                    unclipped.append(bp.unclipped.revcomp())
                elif delta_bp > 0 and delta_bp<len(bp.clipped):
                    clipped.append(bp.unclipped[-delta_bp:] + bp.clipped)
                    unclipped.append(bp.unclipped[:-delta_bp].revcomp())
                elif delta_bp < 0 and delta_bp>-len(bp.unclipped):
                    clipped.append(bp.clipped[-delta_bp:])
                    unclipped.append(bp.clipped[:-delta_bp].revcomp() + bp.unclipped.revcomp())
            elif side is CLIP_LEFT:
                clipped.append(bp.clipped)  # clipped already revcomped
                unclipped.append(bp.unclipped)
                if delta_bp == 0:
                   clipped.append(bp.clipped)  # clipped already revcomped
                   unclipped.append(bp.unclipped)
                elif delta_bp > 0 and delta_bp<len(bp.unclipped):
                    clipped.append(bp.clipped[delta_bp:])
                    unclipped.append(bp.clipped[:delta_bp].revcomp() + bp.unclipped)
                elif delta_bp < 0 and delta_bp>-len(bp.clipped):
                    clipped.append(bp.unclipped[:-delta_bp].revcomp() + bp.clipped)
                    unclipped.append(bp.unclipped[-delta_bp:])
        clipped_cons = find_consensus(clipped)
        unclipped_cons = find_consensus(unclipped)
        if len(clipped_cons) <= CONFIG['discovery']['min_clip_len']:
            Breakpoint.stats[side]['clipped_failed'] += 1
            if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]): print(f"removing {'L' if side is CLIP_LEFT else 'R'} breakpoint because clipped breakpoint did not converge, breakpoints={Counter(bps)}, {','.join([str(bp.clipped) for bp in breakpoints])}")
            return None
        if len(unclipped_cons) <= 40:
            Breakpoint.stats[side]['unclipped_failed'] += 1
            if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]): print(f"removing {'L' if side is CLIP_LEFT else 'R'} breakpoint because unclipped consensus breakpoint did not converge, breakpoints={Counter(bps)} {','.join([str(bp.unclipped) for bp in breakpoints])}")

            return None
        # check if unclipped sequence is n-polymer
        for n in range(1,5):
            query = unclipped_cons[:n] * int(24/n)
            if unclipped_cons[:24] == query:
                Breakpoint.stats[side]['polymer'] += 1
                if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]): print(f"removing {'L' if side is CLIP_LEFT else 'R'} {best_bp}  breakpoint due to {n}-polymer at breakpoint: {str(query)} in {str(unclipped_cons)}")
                return None
        Breakpoint.stats[side]['passed'] += 1
        if Breakpoint.DEBUG_check_if_breakpoint_of_interest(breakpoints[0]): print(f"passed {'L' if side is CLIP_LEFT else 'R'} {best_bp} with clip={clipped_cons} and unclipped={unclipped_cons} of {len(clipped)} reads.")
        b = Breakpoint(side, reference_name, best_bp, None, clipped_cons, unclipped_cons, None, None, False)
        b.mates = mates
        b.n_reads = len(clipped)
        return b







