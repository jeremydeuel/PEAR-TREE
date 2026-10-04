"""Genome sequence access for tools/rte.

Two back-ends with one interface (`fetch(contig, start, end) -> str`, 0-based half-open,
upper-case, '' when out of range / unknown contig):

* TwoBitGenome  -- py2bit on a .2bit (the discovery genome `genome_2bit`, GRCh38/hg19/hs1, or
                   the remap genome `remap_2bit`, hs1).
* FastaGenome   -- small FASTA files (tests, fixtures). A record named `contig:start-end` is a
                   REGION of `contig` starting at `start`, so fixtures can carry real coordinates
                   without the whole chromosome.

Contig lookup is chr-prefix tolerant (chr8 <-> 8), mirroring annotate_v2.GeneModel.
"""
from __future__ import annotations

import re

_REGION = re.compile(r"^(.+):(\d+)-(\d+)$")


def _alt(contig: str) -> str:
    return contig[3:] if contig.startswith("chr") else "chr" + contig


class TwoBitGenome:
    def __init__(self, path: str):
        import py2bit
        self.path = path
        self._tb = py2bit.open(path)
        self._chroms = self._tb.chroms()

    def _resolve(self, contig):
        if contig in self._chroms:
            return contig
        a = _alt(contig)
        return a if a in self._chroms else None

    def length(self, contig):
        c = self._resolve(contig)
        return self._chroms[c] if c else 0

    def fetch(self, contig: str, start: int, end: int) -> str:
        c = self._resolve(contig)
        if c is None:
            return ""
        n = self._chroms[c]
        start, end = max(0, start), min(n, end)
        if end <= start:
            return ""
        return self._tb.sequence(c, start, end).upper()


class FastaGenome:
    def __init__(self, *paths: str, records: dict | None = None):
        # contig -> list of (offset, seq)
        self._regions: dict[str, list[tuple[int, str]]] = {}
        for p in paths:
            import mappy
            for name, seq, _ in mappy.fastx_read(p):
                self.add(name, seq)
        for name, seq in (records or {}).items():
            self.add(name, seq)

    def add(self, name: str, seq: str):
        m = _REGION.match(name)
        if m:
            contig, off = m.group(1), int(m.group(2))
        else:
            contig, off = name, 0
        self._regions.setdefault(contig, []).append((off, seq.upper()))

    def _resolve(self, contig):
        if contig in self._regions:
            return contig
        a = _alt(contig)
        return a if a in self._regions else None

    def length(self, contig):
        c = self._resolve(contig)
        return max(o + len(s) for o, s in self._regions[c]) if c else 0

    def fetch(self, contig: str, start: int, end: int) -> str:
        c = self._resolve(contig)
        if c is None:
            return ""
        for off, seq in self._regions[c]:
            if start >= off and start < off + len(seq):
                return seq[max(0, start - off):max(0, min(len(seq), end - off))]
        return ""

    def regions(self):
        for c, lst in self._regions.items():
            for off, seq in lst:
                yield c, off, seq


def open_genome(spec):
    """None -> None; a genome object -> itself; a path -> TwoBitGenome/FastaGenome by suffix."""
    if spec is None or spec == "":
        return None
    if hasattr(spec, "fetch"):
        return spec
    s = str(spec)
    if s.endswith(".2bit"):
        return TwoBitGenome(s)
    return FastaGenome(s)
