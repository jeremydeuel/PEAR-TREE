"""annotate_v2 <-> rust/peartree-rte file contract (rust/peartree-rte/SPEC.md "Interface").

Not wired into annotate_v2 yet (the integration stage does that). Three helpers:

    write_inputs(inputs, path)     {key: InsertionInput} -> JSONL, one object per insertion, in order
    write_config(cfg, path)        CONFIG['annotate'] -> JSON (callables dropped; rte_library made
                                   absolute so the binary needs no repository-relative lookup)
    read_records(path)             the binary's TSV -> {locus: RustRecord} (row() / COLUMNS /
                                   element / tags / detail ... as annotate_v2 and locus_class read
                                   them from an RteRecord)
"""
from __future__ import annotations

import csv
import gzip
import json
import os

from .record import RteRecord


def _open(path, mode):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def write_inputs(inputs, path):
    with _open(path, "wt") as fh:
        for key, i in inputs.items():
            fh.write(json.dumps({"locus": key, "left_seq": i.left_seq or "", "right_seq": i.right_seq or "",
                                 "pseudogene_genes": list(i.pseudogene_genes or []),
                                 "legacy_class": i.legacy_class,
                                 "sv": list(i.sv) if i.sv is not None else None}) + "\n")


def _plain(v):
    if callable(v):
        return None
    if isinstance(v, dict):
        return {str(k): _plain(x) for k, x in v.items() if not callable(x)}
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    return str(v)


def write_config(cfg, path):
    from .library import resolve_library_path
    out = {k: _plain(v) for k, v in cfg.items() if not callable(v)}
    if out.get("rte_library"):
        out["rte_library"] = os.path.abspath(resolve_library_path(out["rte_library"]))
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)


class RustRecord:
    """The fields of an RteRecord that annotate_v2 / locus_class read, from one output row."""
    COLUMNS = RteRecord.COLUMNS

    def __init__(self, r):
        self._row = [r[c] for c in RteRecord.COLUMNS]
        self.insertion_id = r["locus"]
        self.element = r["element"]
        self.structure = r["structure"]
        self.tags = [] if r["tags"] == "." else r["tags"].split(",")
        self.tprt_call = r["tprt_call"]
        self.tprt_score = r["tprt_score"]
        self.tprt_points = r["tprt_points"]
        self.consensus = "" if r["consensus"] == "." else r["consensus"]
        self.detail = json.loads(r["rte_detail_json"])
        self.gt_reads = int(r["gt_reads"])
        self.gt_changed = "" if r["gt_changed"] == "." else r["gt_changed"]
        self.strand = int(r["strand"])
        dot = lambda v: None if v == "." else v      # noqa: E731
        self.site = (dot(r["site_contig"]), None if r["site_L"] == "." else int(r["site_L"]),
                     None if r["site_R"] == "." else int(r["site_R"]))

    def row(self):
        return list(self._row)


def read_records(path):
    with _open(path, "rt") as fh:
        return {r["locus"]: RustRecord(r) for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)}
