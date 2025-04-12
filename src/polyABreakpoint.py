import pysam
from quality_seq import QualitySeq
from breakpoint import CLIP_RIGHT, CLIP_LEFT
from sequence_checks import clean_clipped_seq
POLYA_CUTOFF = 12
from config import CONFIG
POLYA_SEQ = 'A' * POLYA_CUTOFF
POLYT_SEQ = 'T' * POLYA_CUTOFF
class PolyABreakpoint:
    """
    In contrast to the other breakpoints, this starts with a polyA containing read
    """


    def __init__(self, read: pysam.AlignedSegment, polyA=True):
        self.seq: QualitySeq = QualitySeq(read.seq, read.query_qualities)
        self.polyA: bool= polyA
        self.qname: str = read.query_name
        self.is_forward: bool = read.is_forward
        self.is_read1:bool = read.is_read1
        self.polyA_len: int = 0
        self.firstInside: int|None= None #first base inside the RE
        self.firstOutside: int|None = None #first base outside the RE
        self.find_parts()
        assert self.polyA_len >= POLYA_CUTOFF, f"poly{'A' if self.polyA else 'T'} is {self.polyA_len}, should be at least {POLYA_CUTOFF}, seq = {read.query_sequence}"

        self.reference_name: str |None = None
        self.breakpoint: int|None = None
        self.clipped: QualitySeq|None = None
        self.clip: int = CLIP_RIGHT if self.is_forward ^ self.polyA else CLIP_LEFT
        if self.polyA:
            #clean_clipped_seq has to executed in the direction of the read
            if self.firstInside is not None and self.firstInside > 6:
                if read.is_read1 ^ read.is_forward:
                    self.clipped = clean_clipped_seq(self.seq[:self.firstInside].revcomp())
                else:
                    self.clipped = clean_clipped_seq(self.seq[:self.firstInside]).revcomp()
        else:
            if read.is_read1 ^ read.is_forward:
                self.clipped = clean_clipped_seq(self.seq[(self.firstOutside+1):].revcomp()).revcomp()
            else:
                self.clipped = clean_clipped_seq(self.seq[(self.firstOutside + 1):])

    def setMate(self, mateRead: pysam.AlignedSegment):
        self.reference_name = mateRead.reference_name
        if mateRead.mapq < CONFIG['discovery']['min_mapq']: return
        if self.is_forward ^ self.polyA ^ mateRead.is_forward:
            self.breakpoint = mateRead.reference_start
        else:
            self.breakpoint = mateRead.reference_end
        if self.breakpoint is None:
            self.reference_name = None

    def find_parts(self):
        current_start = 0
        current_end = 0
        in_A = False
        if self.polyA:
            for i, b in enumerate(self.seq._sequence):
                if b == 'A':
                    if in_A:
                        current_end = i
                    else:
                        in_A = True
                        current_start = i
                        current_end = i
                else:
                    if current_end - current_start + 1> self.polyA_len:
                        self.firstInside = current_start
                        self.firstOutside = current_end
                        self.polyA_len = current_end - current_start + 1
                    current_end = 0
                    current_start = 0
                    in_A = False
            if current_end - current_start + 1 > self.polyA_len:
                self.firstInside = current_start
                self.firstOutside = current_end
                self.polyA_len = current_end - current_start + 1
        else:
            for i, b in enumerate(self.seq._sequence):
                if b == 'T':
                    if in_A:
                        current_end = i
                    else:
                        in_A = True
                        current_start = i
                        current_end = i
                else:
                    if current_end - current_start + 1 > self.polyA_len:
                        self.firstInside = current_end
                        self.firstOutside = current_start
                        self.polyA_len = current_end - current_start + 1
                    current_end = 0
                    current_start = 0
                    in_A = False
            if current_end - current_start + 1 > self.polyA_len:
                self.firstInside = current_end
                self.firstOutside = current_start
                self.polyA_len = current_end - current_start + 1

    @staticmethod
    def findPolyA(read: pysam.AlignedSegment) -> 'PolyABreakpoint':
        if read.is_proper_pair: return None
        if not read.mate_is_mapped: return None #read needs to be mapped for this to work
        b = None
        if POLYA_SEQ in read.query_sequence:
            if read.is_reverse:
                if read.is_read2:
                   b = PolyABreakpoint(read, polyA=True)
            else:
                if read.is_read1:
                    b = PolyABreakpoint(read, polyA=True)
        elif POLYT_SEQ in read.query_sequence:
            if read.is_reverse:
                if read.is_read1:
                    b = PolyABreakpoint(read, polyA=False)
            else:
                if read.is_read2:
                    b = PolyABreakpoint(read, polyA = False)
        if b is None: return None
        if not b.clipped: return None
        return b


