"""Truth-TSV schema shared by val1 and fullstack (documented in TRUTH_SCHEMA.md)."""

# label columns appended to every truth table (val1 after its 8 legacy columns; fullstack
# after its coordinate columns)
LABEL_COLUMNS = [
    "role", "type_id", "variant", "element", "structure", "tags", "strand", "ins_len",
    "polya_len", "polya_side", "tsd_seq", "en_motif", "en_mismatches", "mh_seq",
    "element_id", "source_id", "subfamily", "samples", "vaf_by_sample",
    "frags_R", "frags_L", "reads_R", "reads_L", "parts", "info",
]

# legacy val1 classes -> SPEC element (rows emitted by the pre-existing val1 generators)
_LEGACY_ELEMENT = {"L1HS": "L1", "ALUY": "ALU", "SVA": "SVA", "PSEUDOGENE": "PSEUDOGENE"}


def legacy_labels(cls):
    el = "UNKNOWN"
    for k, v in _LEGACY_ELEMENT.items():
        if cls.upper().startswith(k):
            el = v
    d = {c: "." for c in LABEL_COLUMNS}
    d.update(role="TP", type_id="0", variant=f"legacy:{cls}", element=el, samples="1")
    return d


def fmt(v):
    if v is None or v == "":
        return "."
    if isinstance(v, float):
        return f"{v:.3f}"
    if isinstance(v, (list, tuple)):
        return ",".join(fmt(x) for x in v) if v else "."
    if isinstance(v, dict):
        return ";".join(f"{k}={v[k]}" for k in sorted(v)) if v else "."
    return str(v)


def event_labels(ev, tr, samples, vafs, counts_by_sample):
    """Build the label-column dict for one catalogue event."""
    from .models import parts_str
    d = {c: "." for c in LABEL_COLUMNS}
    d.update(role=ev.role, type_id=ev.type_id, variant=ev.key, element=ev.element,
             structure=ev.structure, tags=",".join(ev.tags) if ev.tags else ".",
             strand=tr.get("strand", "."), ins_len=tr.get("ins_len", 0), polya_len=ev.polya_len,
             polya_side=tr.get("polya_side", "."), tsd_seq=tr.get("tsd_seq") or ".",
             en_motif=tr.get("en_motif") or ".", en_mismatches=tr.get("en_mm", "."),
             mh_seq=tr.get("mh_seq") or ".", element_id=ev.element_id, source_id=ev.source_id,
             subfamily=ev.subfamily, samples=",".join(str(s) for s in samples),
             vaf_by_sample=",".join(f"{v:.3f}" for v in vafs),
             frags_R=",".join(str(c["R_frags"]) for c in counts_by_sample),
             frags_L=",".join(str(c["L_frags"]) for c in counts_by_sample),
             reads_R=",".join(str(c["R_reads"]) for c in counts_by_sample),
             reads_L=",".join(str(c["L_reads"]) for c in counts_by_sample),
             parts=parts_str(tr.get("parts") or []) or ".", info=ev.info)
    return {k: fmt(v) for k, v in d.items()}
