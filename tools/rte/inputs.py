"""Readers for the combine -> annotate sidecars (plans/tprt_hallmarks/SPEC.md).

`<patient>.insertions.evidence.tsv.gz`  one row per insertion junction (LEFT / RIGHT; a
    separate poly-A row may use side POLYA). Columns by header name: insertion_id, side, n_reads,
    n_fragments, n_independent, n_samples, n_mates, supported, clip_consensus, consensus_depth,
    polya_len_median, polya_len_range, beyond_polya, beyond_polya_support (+ optional
    cross_sample_identical). Missing columns default to 0 / ''.

`<patient>.insertions.reads.fa.gz`  FASTA, name `insertion_id|side|role|sample|frag|r12`,
    sequence in reference-forward orientation (mates already oriented like the allele).
"""
from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass, field


def _int(v, d=0):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return d


def _float(v, d=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


@dataclass
class JunctionEvidence:
    side: str
    n_reads: int = 0
    n_fragments: int = 0
    n_independent: int = 0
    n_samples: int = 0
    n_mates: int = 0
    supported: int = 0
    clip_consensus: str = ""
    polya_len_median: float = 0.0
    polya_len_range: str = ""
    beyond_polya: str = ""
    beyond_polya_support: int = 0
    cross_sample_identical: int = 0
    n_short_used: int = 0           # SHORT overhang fragments counted in n_independent
    n_short_mate_inside: int = 0    # ... of which the mate lies inside the element

    @classmethod
    def from_row(cls, r: dict) -> "JunctionEvidence":
        return cls(side=(r.get("side") or "").upper(),
                   n_reads=_int(r.get("n_reads")), n_fragments=_int(r.get("n_fragments")),
                   n_independent=_int(r.get("n_independent")), n_samples=_int(r.get("n_samples")),
                   n_mates=_int(r.get("n_mates")), supported=_int(r.get("supported")),
                   clip_consensus=(r.get("clip_consensus") or "").strip(".") ,
                   polya_len_median=_float(r.get("polya_len_median")),
                   polya_len_range=r.get("polya_len_range") or "",
                   beyond_polya=(r.get("beyond_polya") or "").strip("."),
                   beyond_polya_support=_int(r.get("beyond_polya_support")),
                   cross_sample_identical=_int(r.get("cross_sample_identical")),
                   n_short_used=_int(r.get("n_short_used")),
                   n_short_mate_inside=_int(r.get("n_short_mate_inside")))


@dataclass
class EvidenceRead:
    side: str
    role: str
    sample: str
    frag: str
    r12: str
    seq: str

    @property
    def fragment_key(self):
        return (self.sample, self.frag)


@dataclass
class InsertionEvidence:
    insertion_id: str
    junctions: dict = field(default_factory=dict)   # side -> JunctionEvidence
    reads: list = field(default_factory=list)       # [EvidenceRead]


def read_evidence_tsv(path: str, store: dict | None = None) -> dict:
    store = {} if store is None else store
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        rdr = csv.DictReader((l for l in fh if l.strip()), delimiter="\t")
        for r in rdr:
            iid = r.get("insertion_id")
            if not iid:
                continue
            ev = store.setdefault(iid, InsertionEvidence(iid))
            j = JunctionEvidence.from_row(r)
            ev.junctions[j.side] = j
    return store


def read_reads_fa(path: str, store: dict | None = None, wanted=None) -> dict:
    """Parse the pooled reads FASTA. `wanted` (a set of insertion ids) skips the rest."""
    store = {} if store is None else store
    op = gzip.open if str(path).endswith(".gz") else open
    name = None
    chunks = []

    def flush():
        if name is None:
            return
        parts = name.split("|")
        if len(parts) < 6:
            parts += [""] * (6 - len(parts))
        iid = "|".join(parts[:-5]) if len(parts) > 6 else parts[0]
        side, role, sample, frag, r12 = parts[-5:] if len(parts) > 6 else parts[1:6]
        if wanted is not None and iid not in wanted:
            return
        ev = store.setdefault(iid, InsertionEvidence(iid))
        ev.reads.append(EvidenceRead(side.upper(), role.upper(), sample, frag, r12, "".join(chunks)))

    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0] if line[1:].strip() else ""
                chunks = []
            elif line:
                chunks.append(line.strip())
    flush()
    return store
