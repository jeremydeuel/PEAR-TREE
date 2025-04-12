from typing import List, Tuple, Iterator
from ctypes import c_ubyte
from revcomp import revcomp
class QualitySeq:
    def __init__(self, seq: str, quality: List[c_ubyte]):
        assert len(seq) == len(quality)
        self._sequence = seq
        self._quality = list(quality)

    def split(self, breakpoint) -> Tuple['QualitySeq', 'QualitySeq']:
        return QualitySeq(self._sequence[:breakpoint], self._quality[:breakpoint]), \
               QualitySeq(self._sequence[breakpoint:], self._quality[breakpoint:])

    def __len__(self) -> int:
        return len(self._sequence)

    def __iter__(self) -> Iterator[Tuple[str, int]]:
        for i in range(len(self._sequence)):
            yield self._sequence[i], self._quality[i]

    def __str__(self) -> str:
        return self._sequence

    def revcomp(self) -> 'QualitySeq':
        return QualitySeq(revcomp(self._sequence), list(reversed(self._quality)))

    def upper(self) -> 'QualitySeq':
        return QualitySeq(self._sequence.upper(), self._quality)

    def lower(self) -> 'QualitySeq':
        return QualitySeq(self._sequence.lower(), self._quality)

    def __getitem__(self, item):
        return QualitySeq(self._sequence[item],self._quality[item])

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

