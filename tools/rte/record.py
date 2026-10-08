"""Structured per-insertion result of the tools/rte annotation (replaces string parsing for
the new fields). `RteRecord.COLUMNS` is the column order appended to annotate_v2's table
(plans/tprt_hallmarks/SPEC.md "annotate output -- new columns", plus `rte_detail`)."""
from __future__ import annotations

from dataclasses import dataclass, field


def _f(v, nd=4):
    if v is None or v == "":
        return "."
    if isinstance(v, float):
        return f"{v:.{nd}g}" if abs(v) < 1 else f"{round(v, 2):g}"
    return str(v)


@dataclass
class RteRecord:
    insertion_id: str
    element: str = "UNKNOWN"
    structure: str = "5P_UNRESOLVED"
    tags: list = field(default_factory=list)
    covered_5p: int | None = None
    covered_3p: int | None = None
    covered_intervals: list = field(default_factory=list)
    consensus: str = ""
    element_identity: float | None = None
    nearest_active: str = "."
    tsd_seq: str = ""
    tsd_len: int | None = None
    en_motif: str = ""
    en_mismatches: int | None = None
    polya_len: float | None = None
    beyond_polya: str = ""
    beyond_polya_support: int = 0
    strand: int = 0
    tprt_score: float = 0.0
    tprt_points: str = "."
    tprt_call: str = "UNCERTAIN"
    detail: dict = field(default_factory=dict)
    score_input: object = None      # score.ScoreInput (kept for re-scoring passes)
    covered_seqs: list = field(default_factory=list)
    site: tuple = (None, None, None)   # (contig, L, R) on the discovery genome
    # genotype2 extra-pass reads (GT_* roles, <P>.insertions.genotype_reads.fa.gz) used in the
    # assembly, and the calls they changed vs the combine reads alone (annotator.gt_changes);
    # written as the optional gt_reads / gt_changed columns, never part of COLUMNS
    gt_reads: int = 0
    gt_changed: str = ""

    COLUMNS = ["element", "structure", "tags", "covered_5p", "covered_3p", "element_identity",
               "nearest_active", "tsd_seq", "tsd_len", "en_motif", "en_mismatches", "polya_len",
               "beyond_polya", "tprt_score", "tprt_points", "tprt_call", "rte_detail"]

    def detail_string(self):
        d = dict(self.detail)
        if self.consensus:
            d.setdefault("consensus", self.consensus)
        if self.covered_intervals:
            d["covered"] = ",".join(f"{a}-{b}" for a, b in self.covered_intervals)
        if self.strand:
            d["strand"] = "+" if self.strand > 0 else "-"
        if self.beyond_polya_support:
            d["beyond_polya_support"] = self.beyond_polya_support
        return ";".join(f"{k}={v}" for k, v in d.items()) or "."

    def row(self):
        vals = [self.element, self.structure, ",".join(self.tags) or ".",
                _f(self.covered_5p), _f(self.covered_3p), _f(self.element_identity),
                self.nearest_active or ".", self.tsd_seq or ".", _f(self.tsd_len),
                self.en_motif or ".", _f(self.en_mismatches), _f(self.polya_len),
                self.beyond_polya or ".", _f(self.tprt_score), self.tprt_points or ".",
                self.tprt_call, self.detail_string()]
        return [str(v).replace("\t", " ").replace("\n", " ") for v in vals]

    @classmethod
    def empty_row(cls):
        return ["."] * len(cls.COLUMNS)
