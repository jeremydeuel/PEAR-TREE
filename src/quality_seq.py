from typing import List, Tuple, Iterator
from revcomp import revcomp


class QualitySeq:
    """
    A DNA sequence paired with a per-base quality/score track.

    The quality track holds phred qualities for raw reads but *consensus scores*
    for consensus sequences, and those can exceed the 0-255 phred range, so it is
    kept as a plain Python ``list`` of ints (an ``array('B')`` would overflow and,
    measured, an ``array('i')`` is actually slower to build at read length).

    The track is treated as immutable: no method mutates it in place, so it can
    be shared between instances. ``__init__`` therefore stores an existing list
    by reference (avoiding a copy on every slice/upper/revcomp result) and only
    copies when handed a non-list (e.g. pysam's ``array('B')`` query_qualities).
    """
    __slots__ = ('_sequence', '_quality')

    def __init__(self, seq: str, quality):
        self._sequence = seq
        self._quality = quality if type(quality) is list else list(quality)
        assert len(seq) == len(self._quality)

    def split(self, breakpoint) -> Tuple['QualitySeq', 'QualitySeq']:
        return QualitySeq(self._sequence[:breakpoint], self._quality[:breakpoint]), \
               QualitySeq(self._sequence[breakpoint:], self._quality[breakpoint:])

    def __len__(self) -> int:
        return len(self._sequence)

    def __iter__(self) -> Iterator[Tuple[str, int]]:
        return zip(self._sequence, self._quality)

    def __str__(self) -> str:
        return self._sequence

    def revcomp(self) -> 'QualitySeq':
        return QualitySeq(revcomp(self._sequence), self._quality[::-1])

    def upper(self) -> 'QualitySeq':
        upper_seq = self._sequence.upper()
        if upper_seq == self._sequence:
            return self  # already uppercase (the common case) — skip the copy
        return QualitySeq(upper_seq, self._quality)

    def lower(self) -> 'QualitySeq':
        lower_seq = self._sequence.lower()
        if lower_seq == self._sequence:
            return self
        return QualitySeq(lower_seq, self._quality)

    def __getitem__(self, item):
        return QualitySeq(self._sequence[item], self._quality[item])

    def getPos(self, position: int) -> Tuple[str, int]:
        return self._sequence[position], self._quality[position]

    def qual(self, item = slice(None)):
        return self._quality[item]

    def seq(self, item = slice(None)):
        return self._sequence[item]

    def __add__(self, other: 'QualitySeq'):
        return QualitySeq(self._sequence + other._sequence, self._quality + other._quality)

    def __mul__(self, other: int) -> str:
        return self._sequence * other

    def __eq__(self, other: str) -> bool:
        return self._sequence == other

    def phred(self, base=33):
        return "".join([chr(i + base if i + base <= 126 else 126) for i in self._quality])

    def fastq(self, title, phred_base=33):
        return "@" + title + "\n" + self._sequence + "\n+\n" + self.phred(phred_base) + "\n"
