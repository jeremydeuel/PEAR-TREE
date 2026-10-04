"""Insertion-type catalogue: build an `Event` and apply it to a reference window.

Every event is modelled in the ELEMENT-SENSE ("+") orientation and then mapped back:

    TSD  (t>0):  alt = ref[:s+t] + X + ref[s:]        RIGHT bp = s+t, LEFT bp = s
    DEL  (d>0):  alt = ref[:s-d] + X + ref[s:]        RIGHT bp = s-d, LEFT bp = s  (TSD_DELETION)
    NONE      :  alt = ref[:s]   + X + ref[s:]        RIGHT bp = LEFT bp = s
    L1DEL (D) :  alt = ref[:s-D] + X[mh:] + ref[s:]   5' microhomology X[:mh] == ref[s-D-mh:s-D]

`s` is the L1 endonuclease nick: top strand 5'-TT|AAAA-3' (bottom 5'-TTTT/AA-3'), so the
poly-A of a sense-strand insertion abuts the second TSD copy (LEFT breakpoint) — exactly the
PEAR-TREE convention (LEFT = 3'/poly-A junction for '+' insertions, RIGHT = 5' junction). A
'-' insertion is the reverse complement: poly-T at the RIGHT breakpoint, 5' end at LEFT.

`X` = [ref-derived 5' extras: fold-back inverted duplication | templated local copy |
co-inserted local pre-mRNA] + the element-derived insert (`Event.ins`, ends in poly-A).

Truth labels follow plans/tprt_hallmarks/SPEC.md (`element`, `structure`, `tags`) plus
`role` (TP | ARTEFACT) and `type_id` (literature catalogue numbering, see TYPE_IDS).
"""
import random
from dataclasses import dataclass, field

from .seqs import (EN_NICK_OFFSET, degenerate_en_motif, en_mismatches, lognormal_int,
                   make_polya, revcomp, rnd_seq, sample_deletion_len, sample_en_mismatches,
                   sample_l1_truncated_len, sample_polya_len, sample_tsd_len, sample_twin_priming)

# literature catalogue numbering (docs/insertion_types.html). 9-11 are the out-of-scope
# translocation bridge / chimeric bridge / complex reciprocal inversion types; 12 is used
# here for EN-independent insertion. Artefacts carry type_id 0.
TYPE_IDS = {
    1: "solo L1 (full-length / 5'-truncated / 5'-inverted twin priming / TSD deletion)",
    2: "partnered 3' transduction",
    3: "orphan 3' transduction",
    4: "Alu / SVA insertion (incl. SVA 5' and 3' transduction)",
    5: "processed pseudogene (+ no-junction decoy)",
    6: "solitary poly(A/T) insertion",
    7: "L1-mediated deletion",
    8: "L1-mediated tandem duplication",
    9: "translocation bridge (OUT OF SCOPE)",
    10: "chimeric bridge (OUT OF SCOPE)",
    11: "complex reciprocal inversion (OUT OF SCOPE)",
    12: "EN-independent insertion",
    13: "twin priming + 5' switching",
    14: "templated local insertion",
    15: "co-inserted local pre-mRNA",
    16: "fold-back inverted duplication 5' of the target site",
}


@dataclass
class Event:
    key: str
    type_id: int
    role: str                       # TP | ARTEFACT
    element: str
    structure: str
    tags: list
    ins: str                        # element-derived insert, element sense, incl. poly-A
    parts: list = field(default_factory=list)   # [(label, start, end, info)] in `ins` coords
    target: str = "TSD"             # TSD | DEL | NONE | L1DEL
    target_len: int = 0
    mh: int = 0
    extras5: list = field(default_factory=list)  # ref-derived 5' additions (see apply_event)
    polya_len: int = 0
    plant_en: bool = True           # whether the nick carries an EN motif (val1 plants it)
    en_mm: int = 0                  # planned EN-motif mismatch count
    element_id: str = "."
    source_id: str = "."
    subfamily: str = "."
    render: str = "normal"          # normal | single_fragment | slippage | mismap | foldback_reads
    info: dict = field(default_factory=dict)


# ------------------------------------------------------------------------------ helpers
def _l1_structure(rng, seq, structure):
    """5' structure of an L1 body (no poly-A). Returns (body, parts, info)."""
    n = len(seq)
    if structure == "FULL_LENGTH":
        return seq, [("L1", 0, n, f"0-{n}+")], {}
    if structure == "TRUNCATED_5P":
        k = sample_l1_truncated_len(rng, n)
        return seq[n - k:], [("L1", 0, k, f"{n - k}-{n}+")], {"trunc_at": n - k}
    tp = sample_twin_priming(rng, n)
    a, b, c = tp["a"], tp["b"], tp["c"]
    inv = revcomp(seq[a:b])
    fwd = seq[c:]
    info = {"inv_breakpoint": c, "inv_a": a, "inv_b": b, "inv_junction": tp["junction"],
            "inv_jlen": tp["jlen"], "fwd_len": n - c, "inv_len": b - a}
    if structure == "INVERTED_5P":
        body = inv + fwd
        parts = [("L1_INV", 0, len(inv), f"{a}-{b}-"), ("L1", len(inv), len(body), f"{c}-{n}+")]
        return body, parts, info
    # INVERTED_5P_SWITCH: the internally primed cDNA switches to an upstream (more 5')
    # region of the RNA -> an extra segment, forward relative to the element, 5' of the
    # inverted part (Zumalave 2026 Fig. 2c "twin priming & 5' switching").
    if a < 120:
        a2, b2 = 0, max(40, a)
    else:
        gap = rng.randint(0, min(300, a - 100))
        b2 = a - gap
        a2 = max(0, b2 - rng.randint(80, 600))
    sw = seq[a2:b2]
    body = sw + inv + fwd
    parts = [("L1_SWITCH", 0, len(sw), f"{a2}-{b2}+"),
             ("L1_INV", len(sw), len(sw) + len(inv), f"{a}-{b}-"),
             ("L1", len(sw) + len(inv), len(body), f"{c}-{n}+")]
    info.update({"switch_a": a2, "switch_b": b2})
    return body, parts, info


def _with_polya(rng, body, parts, cls, scale):
    pa = sample_polya_len(rng, cls, scale)
    tail = make_polya(rng, pa)
    parts = parts + [("POLYA", len(body), len(body) + pa, str(pa))]
    return body + tail, parts, pa


def _tsd_target(rng):
    return "TSD", sample_tsd_len(rng)


# ------------------------------------------------------------------------------ catalogue
class Ctx:
    """Per-run context for event builders."""

    def __init__(self, lib, genes=None, polya_scale=1.0, max_del=20000, max_dup=5000):
        self.lib = lib
        self.genes = genes or []
        self.polya_scale = polya_scale
        self.max_del = max_del
        self.max_dup = max_dup


def _l1_event(rng, ctx, key, type_id, structure, tags=(), elem=None):
    el = elem or ctx.lib.element(rng, "L1")
    body, parts, info = _l1_structure(rng, el.seq, structure)
    ins, parts, pa = _with_polya(rng, body, parts, "L1", ctx.polya_scale)
    target, tl = _tsd_target(rng)
    return Event(key, type_id, "TP", "L1", structure, list(tags), ins, parts, target, tl,
                 polya_len=pa, element_id=el.id, subfamily=el.subfamily, info=info,
                 en_mm=sample_en_mismatches(rng))


def b_l1_full(rng, ctx):
    return _l1_event(rng, ctx, "L1_FULL", 1, "FULL_LENGTH")


def b_l1_trunc(rng, ctx):
    return _l1_event(rng, ctx, "L1_TRUNC", 1, "TRUNCATED_5P")


def b_l1_inv(rng, ctx):
    return _l1_event(rng, ctx, "L1_INV", 1, "INVERTED_5P")


def b_l1_inv_switch(rng, ctx):
    return _l1_event(rng, ctx, "L1_INV_SWITCH", 13, "INVERTED_5P_SWITCH")


def b_l1_tsd_deletion(rng, ctx):
    ev = _l1_event(rng, ctx, "L1_TSD_DELETION", 1, "TRUNCATED_5P", ["TSD_DELETION"])
    ev.target, ev.target_len = "DEL", lognormal_int(rng, 8, 0.7, 1, 30)   # Nam 2023 negative TSDs
    return ev


def _td3p_tag(rng, src):
    """Real downstream flank up to one of the source's fixed termination offsets (+-3 bp)."""
    end = rng.choice(src.endpoints) + rng.randint(-3, 3)
    end = min(max(end, 20), len(src.flank3))
    return src.flank3[:end], end


def b_l1_td3p(rng, ctx):
    src = ctx.lib.source(rng, "L1")
    u = rng.random()
    structure = "TRUNCATED_5P" if u < 0.7 else ("INVERTED_5P" if u < 0.85 else "FULL_LENGTH")
    body, parts, info = _l1_structure(rng, src.element_seq, structure)
    tag, end = _td3p_tag(rng, src)
    pre = src.tail                   # the source copy's own genomic A-rich tail
    p0 = len(body)
    parts = parts + ([("SRC_TAIL", p0, p0 + len(pre), "")] if pre else []) + \
        [("TD3P", p0 + len(pre), p0 + len(pre) + len(tag), f"{src.id}:0-{end}")]
    ins, parts, pa = _with_polya(rng, body + pre + tag, parts, "L1", ctx.polya_scale)
    target, tl = _tsd_target(rng)
    info.update({"td_len": len(tag), "td_end": end})
    return Event("L1_TD3P", 2, "TP", "L1", structure, ["TD3P", f"TD3P_SOURCE={src.id}"], ins, parts,
                 target, tl, polya_len=pa, element_id=src.id, source_id=src.id,
                 subfamily=src.subfamily, info=info, en_mm=sample_en_mismatches(rng))


def b_orphan_td3p(rng, ctx):
    src = ctx.lib.source(rng, "L1")
    tag, end = _td3p_tag(rng, src)
    # orphan: the element part is lost (severe 5' truncation into the flank); keep a random
    # 5' cut inside the tag so tags start at varied flank offsets
    cut = rng.randint(0, max(0, len(tag) - 40)) if len(tag) > 80 and rng.random() < 0.5 else 0
    body = tag[cut:]
    parts = [("TD3P", 0, len(body), f"{src.id}:{cut}-{end}")]
    ins, parts, pa = _with_polya(rng, body, parts, "ORPHAN_TD", ctx.polya_scale)
    target, tl = _tsd_target(rng)
    return Event("ORPHAN_TD3P", 3, "TP", "ORPHAN_TD", "5P_UNRESOLVED",
                 ["TD3P", f"TD3P_SOURCE={src.id}"], ins, parts, target, tl, polya_len=pa,
                 source_id=src.id, subfamily=src.subfamily,
                 info={"td_len": len(body), "td_start": cut, "td_end": end},
                 en_mm=sample_en_mismatches(rng))


def _alu(rng, ctx, key, sub):
    el = ctx.lib.element(rng, "ALU", sub)
    if rng.random() < 0.1:                       # rare 5'-truncated Alu
        k = rng.randint(120, len(el.seq) - 20)
        body, structure = el.seq[-k:], "TRUNCATED_5P"
        parts = [("ALU", 0, k, f"{len(el.seq) - k}-{len(el.seq)}+")]
    else:
        body, structure = el.seq, "FULL_LENGTH"
        parts = [("ALU", 0, len(body), f"0-{len(body)}+")]
    ins, parts, pa = _with_polya(rng, body, parts, "ALU", ctx.polya_scale)
    target, tl = _tsd_target(rng)
    return Event(key, 4, "TP", "ALU", structure, [], ins, parts, target, tl, polya_len=pa,
                 element_id=el.id, subfamily=el.subfamily, en_mm=sample_en_mismatches(rng))


def b_alu_ya5(rng, ctx):
    return _alu(rng, ctx, "ALU_YA5", "AluYa5")


def b_alu_yb8(rng, ctx):
    return _alu(rng, ctx, "ALU_YB8", "AluYb8")


def _sva_body(rng, seq):
    if rng.random() < 0.6:
        return seq, "FULL_LENGTH", [("SVA", 0, len(seq), f"0-{len(seq)}+")]
    k = rng.randint(300, len(seq) - 50)
    return seq[-k:], "TRUNCATED_5P", [("SVA", 0, k, f"{len(seq) - k}-{len(seq)}+")]


def _sva(rng, ctx, key, sub):
    el = ctx.lib.element(rng, "SVA", sub)
    body, structure, parts = _sva_body(rng, el.seq)
    ins, parts, pa = _with_polya(rng, body, parts, "SVA", ctx.polya_scale)
    target, tl = _tsd_target(rng)
    return Event(key, 4, "TP", "SVA", structure, [], ins, parts, target, tl, polya_len=pa,
                 element_id=el.id, subfamily=el.subfamily, en_mm=sample_en_mismatches(rng))


def b_sva_e(rng, ctx):
    return _sva(rng, ctx, "SVA_E", "SVA_E")


def b_sva_f(rng, ctx):
    return _sva(rng, ctx, "SVA_F", "SVA_F")


def b_sva_td5p(rng, ctx):
    """SVA 5' transduction: transcription starts upstream of the source SVA, so the insert
    carries the source's real upstream flank 5' of a full-length SVA (Damert 2009)."""
    src = ctx.lib.source(rng, "SVA", need_flank5=True)
    n5 = lognormal_int(rng, 400, 0.7, 40, min(3000, len(src.flank5)))
    up = src.flank5[len(src.flank5) - n5:]
    body = up + src.element_seq
    parts = [("TD5P", 0, n5, f"{src.id}:-{n5}-0"), ("SVA", n5, len(body), f"0-{len(src.element_seq)}+")]
    ins, parts, pa = _with_polya(rng, body, parts, "SVA", ctx.polya_scale)
    target, tl = _tsd_target(rng)
    return Event("SVA_TD5P", 4, "TP", "SVA", "FULL_LENGTH", ["TD5P", f"TD5P_SOURCE={src.id}"], ins,
                 parts, target, tl, polya_len=pa, element_id=src.id, source_id=src.id,
                 subfamily=src.subfamily, info={"td5_len": n5}, en_mm=sample_en_mismatches(rng))


def b_sva_td3p(rng, ctx):
    src = ctx.lib.source(rng, "SVA")
    body, structure, parts = _sva_body(rng, src.element_seq)
    tag, end = _td3p_tag(rng, src)
    pre = src.tail
    p0 = len(body)
    parts = parts + ([("SRC_TAIL", p0, p0 + len(pre), "")] if pre else []) + \
        [("TD3P", p0 + len(pre), p0 + len(pre) + len(tag), f"{src.id}:0-{end}")]
    ins, parts, pa = _with_polya(rng, body + pre + tag, parts, "SVA", ctx.polya_scale)
    target, tl = _tsd_target(rng)
    return Event("SVA_TD3P", 4, "TP", "SVA", structure, ["TD3P", f"TD3P_SOURCE={src.id}"], ins,
                 parts, target, tl, polya_len=pa, element_id=src.id, source_id=src.id,
                 subfamily=src.subfamily, info={"td_len": len(tag), "td_end": end},
                 en_mm=sample_en_mismatches(rng))


def _pick_gene(rng, ctx):
    pool = [g for g in ctx.genes if len(g.exon_seqs) >= 2]
    if not pool:
        raise RuntimeError("no gene models available for pseudogene types")
    return rng.choice(pool)


def b_pseudogene(rng, ctx):
    """Processed pseudogene: the spliced mRNA (>= 2 exons, exon-exon junctions inside the
    insert), possibly 5'-truncated, + poly-A, TPRT scar (Ewing 2013)."""
    g = _pick_gene(rng, ctx)
    ex = g.exon_seqs
    first = rng.randint(0, len(ex) - 2) if rng.random() < 0.4 else 0
    body, parts, pos, junctions = "", [], 0, []
    for i in range(first, len(ex)):
        seg = ex[i]
        if i == first and first > 0:              # 5' truncation inside the first kept exon
            seg = seg[rng.randint(0, max(0, len(seg) - 30)):]
        if i > first:
            junctions.append(pos)
        parts.append((f"EXON{i + 1}", pos, pos + len(seg), g.id))
        body += seg
        pos += len(seg)
    ins, parts, pa = _with_polya(rng, body, parts, "PSEUDOGENE", ctx.polya_scale)
    target, tl = _tsd_target(rng)
    structure = "FULL_LENGTH" if first == 0 else "TRUNCATED_5P"
    return Event("PSEUDOGENE", 5, "TP", "PSEUDOGENE", structure, ["EXON_JUNCTION"], ins, parts,
                 target, tl, polya_len=pa, element_id=g.id,
                 info={"gene": g.id, "exon_junctions": ",".join(map(str, junctions)),
                       "gene_loc": f"{g.contig}:{g.exons[0][0]}-{g.exons[-1][1]}:{g.strand}"},
                 en_mm=sample_en_mismatches(rng))


def b_pseudogene_decoy(rng, ctx):
    """Decoy: exonic sequence but NO exon-exon junction (one exon + its downstream intron
    start, i.e. an unspliced chunk). Reads hit exons, but it must not be called PSEUDOGENE."""
    g = _pick_gene(rng, ctx)
    i = rng.randint(0, len(g.exon_seqs) - 2)
    # locate exon i inside the pre-mRNA (sense) and take it plus some intron
    pre = g.premrna
    start = pre.find(g.exon_seqs[i])
    if start < 0:
        start = 0
    elen = len(g.exon_seqs[i])
    extra = rng.randint(40, 300)
    body = pre[start:start + elen + extra]
    parts = [(f"EXON{i + 1}", 0, elen, g.id), ("INTRON", elen, len(body), g.id)]
    ins, parts, pa = _with_polya(rng, body, parts, "PSEUDOGENE", ctx.polya_scale)
    target, tl = _tsd_target(rng)
    return Event("PSEUDOGENE_DECOY", 5, "TP", "UNKNOWN", "5P_UNRESOLVED", [], ins, parts, target,
                 tl, polya_len=pa, element_id=g.id, info={"gene": g.id, "decoy": "no_exon_junction"},
                 en_mm=sample_en_mismatches(rng))


def b_polya_only(rng, ctx):
    pa = sample_polya_len(rng, "POLYA_ONLY", ctx.polya_scale)
    ins = make_polya(rng, pa)
    target, tl = _tsd_target(rng)
    return Event("POLYA_ONLY", 6, "TP", "POLYA_ONLY", "5P_UNRESOLVED", [], ins,
                 [("POLYA", 0, pa, str(pa))], target, tl, polya_len=pa,
                 en_mm=sample_en_mismatches(rng))


def b_l1_med_deletion(rng, ctx):
    """L1-mediated deletion: 3' end by TPRT (poly-A, EN nick), 5' end joined upstream with
    1-5 bp microhomology, the segment between deleted, no TSD (Gilbert 2002;
    Rodriguez-Martin 2020)."""
    ev = _l1_event(rng, ctx, "L1_MED_DELETION", 7, "TRUNCATED_5P", ["L1_MED_DELETION"])
    ev.target, ev.target_len = "L1DEL", sample_deletion_len(rng, 100, ctx.max_del)
    ev.mh = rng.randint(1, 5)
    return ev


def b_l1_med_duplication(rng, ctx):
    ev = _l1_event(rng, ctx, "L1_MED_DUPLICATION", 8, "TRUNCATED_5P", ["L1_MED_DUPLICATION"])
    ev.target, ev.target_len = "TSD", sample_deletion_len(rng, 50, ctx.max_dup)
    return ev


def b_en_independent(rng, ctx):
    """EN-independent insertion: no TSD, no poly-A, both ends truncated (internal L1 piece)."""
    el = ctx.lib.element(rng, "L1")
    n = len(el.seq)
    k = lognormal_int(rng, 600, 0.7, 100, n - 200)
    b = rng.randint(n - 1500 if n > 1600 else k, n - 30)   # 3' truncated: ends before the 3' end
    a = max(0, b - k)
    ins = el.seq[a:b]
    ev = Event("EN_INDEPENDENT", 12, "TP", "L1", "TRUNCATED_5P", ["EN_INDEPENDENT"], ins,
               [("L1", 0, len(ins), f"{a}-{b}+")], "NONE", 0, polya_len=0, plant_en=False,
               element_id=el.id, subfamily=el.subfamily, info={"l1_a": a, "l1_b": b})
    if rng.random() < 0.4:
        ev.target, ev.target_len = "DEL", rng.randint(1, 10)
    return ev


def b_templated_local(rng, ctx):
    """<250 bp copied from within 15 bp of the target site, embedded 5' of the L1 body
    (Zumalave 2026 'templated insertions'; Nam 2023)."""
    ev = _l1_event(rng, ctx, "TEMPLATED_LOCAL", 14, "TRUNCATED_5P", ["TEMPLATED_LOCAL"])
    n = lognormal_int(rng, 60, 0.6, 20, 249)
    off = rng.randint(-15, 15)               # template start relative to the 5'-junction
    ev.extras5.append(("REFCOPY", "TEMPLATED", off, n, rng.choice("+-")))
    ev.info.update({"templ_len": n, "templ_off": off})
    return ev


def b_premrna(rng, ctx):
    """Co-inserted local pre-mRNA: a template switch to an unspliced transcript of a nearby
    gene, inserted 5' of the L1 body (Nam 2023 Fig. 4g)."""
    ev = _l1_event(rng, ctx, "PREMRNA_COINSERT", 15, "TRUNCATED_5P", ["PREMRNA_COINSERT"])
    n = lognormal_int(rng, 300, 0.5, 120, 900)
    off = rng.choice([-1, 1]) * rng.randint(400, 2400)
    ev.extras5.append(("REFCOPY", "PREMRNA", off, n, "+"))
    ev.info.update({"premrna_len": n, "premrna_off": off})
    return ev


def b_foldback(rng, ctx):
    """Fold-back inverted duplication of the sequence immediately 5' (upstream) of the
    target site, between the flank and the element's 5' end (Nam 2023: 0.3%)."""
    ev = _l1_event(rng, ctx, "FOLDBACK_INVDUP_5P", 16, "TRUNCATED_5P", ["FOLDBACK_INVDUP_5P"])
    f = lognormal_int(rng, 60, 0.6, 15, 400)
    g = rng.randint(0, 10)
    ev.extras5.append(("FOLDBACK", f, g))
    ev.info.update({"foldback_len": f, "foldback_gap": g})
    return ev


# ----------------------------------------------------------------- artefacts (role=ARTEFACT)
def _art(ev, kind):
    ev.role, ev.type_id = "ARTEFACT", 0
    ev.info["artefact"] = kind
    return ev


def b_art_ligation(rng, ctx, copies=0, key="ART_LIGATION_CHIMERA"):
    """Chimeric ligation: an L1 3' end (+ poly-A) ligated to an unrelated locus. ONE
    molecule -> one independent fragment (PCR copies optional)."""
    el = ctx.lib.element(rng, "L1")
    k = rng.randint(150, 600)
    body = el.seq[-k:]
    ins, parts, pa = _with_polya(rng, body, [("L1", 0, k, f"{len(el.seq) - k}-{len(el.seq)}+")],
                                 "L1", ctx.polya_scale)
    ev = Event(key, 0, "ARTEFACT", "L1", "TRUNCATED_5P", [], ins, parts, "NONE", 0, polya_len=pa,
               plant_en=False, element_id=el.id, subfamily=el.subfamily, render="single_fragment")
    ev.info.update({"sides": "POLYA", "pcr_copies": copies})
    return _art(ev, "ligation_chimera")


def b_art_ligation_pcr(rng, ctx):
    return b_art_ligation(rng, ctx, copies=rng.randint(2, 5), key="ART_LIGATION_PCR")


def b_art_long_tsd(rng, ctx):
    """Two chimeric molecules overlapping by > 50 bp: looks like an insertion with a
    TSD > 50 bp; one fragment per junction."""
    ev = _l1_event(rng, ctx, "ART_LONG_TSD", 0, "TRUNCATED_5P")
    ev.target, ev.target_len = "TSD", rng.randint(51, 150)
    ev.render, ev.plant_en = "single_fragment", False
    ev.info.update({"sides": "BOTH", "pcr_copies": rng.randint(0, 2)})
    return _art(ev, "chimera_long_tsd")


def b_art_chimeric_ends(rng, ctx):
    """5' end from one element class, 3' end from another (Alu/SVA 5' + L1 3'): PCR /
    ligation chimera between two element-derived molecules."""
    other = ctx.lib.element(rng, rng.choice(["ALU", "SVA"]))
    l1 = ctx.lib.element(rng, "L1")
    n5 = min(len(other.seq), rng.randint(120, 400))
    n3 = rng.randint(150, 600)
    body = other.seq[:n5] + l1.seq[-n3:]
    parts = [(other.cls, 0, n5, f"0-{n5}+"), ("L1", n5, n5 + n3, f"{len(l1.seq) - n3}-{len(l1.seq)}+")]
    ins, parts, pa = _with_polya(rng, body, parts, "L1", ctx.polya_scale)
    ev = Event("ART_CHIMERIC_ENDS", 0, "ARTEFACT", "L1", "TRUNCATED_5P", ["CHIMERIC_ENDS"], ins,
               parts, "TSD", sample_tsd_len(rng), polya_len=pa, plant_en=False,
               element_id=f"{other.id}+{l1.id}", subfamily=f"{other.subfamily}+{l1.subfamily}",
               render="single_fragment")
    ev.info.update({"sides": "BOTH", "pcr_copies": rng.randint(0, 2)})
    return _art(ev, "chimeric_ends")


def b_art_polya_slippage(rng, ctx):
    """Poly-A slippage at a reference A-tract (e.g. an Alu tail): no insertion; reads
    crossing the tract carry exaggerated homopolymer slippage and are soft-clipped there."""
    n = rng.randint(18, 40)
    ev = Event("ART_POLYA_SLIPPAGE", 0, "ARTEFACT", "POLYA_ONLY", "5P_UNRESOLVED", [], "", [],
               "NONE", 0, plant_en=False, render="slippage")
    ev.info.update({"tract_len": n})
    return _art(ev, "polya_slippage")


def b_art_subfamily_mismap(rng, ctx):
    """Mismapping-like clips: reads from a paralogous young L1 (same subfamily) extend into
    THAT copy's unique 3' flank but align onto a reference L1 copy here -> clip = foreign
    flank at the reference element's 3' end; mostly MAPQ 0, a minority confident."""
    src = ctx.lib.source(rng, "L1")
    ev = Event("ART_SUBFAMILY_MISMAP", 0, "ARTEFACT", "L1", "5P_UNRESOLVED", [], "", [], "NONE", 0,
               plant_en=False, render="mismap", source_id=src.id, subfamily=src.subfamily)
    ev.info.update({"paralog": src.id, "ref_elem_len": rng.randint(400, 900)})
    return _art(ev, "subfamily_mismap")


def b_art_foldback(rng, ctx):
    """Library fold-back palindrome: read = genomic + reverse complement of the adjacent
    genomic sequence (a hairpin molecule). 1-3 molecules."""
    ev = Event("ART_FOLDBACK_PALINDROME", 0, "ARTEFACT", "UNKNOWN", "5P_UNRESOLVED", [], "", [],
               "NONE", 0, plant_en=False, render="foldback_reads")
    ev.info.update({"n_molecules": rng.randint(1, 3), "pcr_copies": rng.randint(0, 3)})
    return _art(ev, "foldback_palindrome")


CATALOGUE = {
    "L1_FULL": b_l1_full, "L1_TRUNC": b_l1_trunc, "L1_INV": b_l1_inv,
    "L1_INV_SWITCH": b_l1_inv_switch, "L1_TSD_DELETION": b_l1_tsd_deletion,
    "L1_TD3P": b_l1_td3p, "ORPHAN_TD3P": b_orphan_td3p,
    "ALU_YA5": b_alu_ya5, "ALU_YB8": b_alu_yb8, "SVA_E": b_sva_e, "SVA_F": b_sva_f,
    "SVA_TD5P": b_sva_td5p, "SVA_TD3P": b_sva_td3p,
    "PSEUDOGENE": b_pseudogene, "PSEUDOGENE_DECOY": b_pseudogene_decoy,
    "POLYA_ONLY": b_polya_only, "L1_MED_DELETION": b_l1_med_deletion,
    "L1_MED_DUPLICATION": b_l1_med_duplication, "EN_INDEPENDENT": b_en_independent,
    "TEMPLATED_LOCAL": b_templated_local, "PREMRNA_COINSERT": b_premrna,
    "FOLDBACK_INVDUP_5P": b_foldback,
    "ART_LIGATION_CHIMERA": b_art_ligation, "ART_LIGATION_PCR": b_art_ligation_pcr,
    "ART_LONG_TSD": b_art_long_tsd, "ART_CHIMERIC_ENDS": b_art_chimeric_ends,
    "ART_POLYA_SLIPPAGE": b_art_polya_slippage, "ART_SUBFAMILY_MISMAP": b_art_subfamily_mismap,
    "ART_FOLDBACK_PALINDROME": b_art_foldback,
}
TP_TYPES = [k for k in CATALOGUE if not k.startswith("ART_")]
ART_TYPES = [k for k in CATALOGUE if k.startswith("ART_")]
_KEY_TYPE_ID = {"L1_FULL": 1, "L1_TRUNC": 1, "L1_INV": 1, "L1_INV_SWITCH": 13,
                "L1_TSD_DELETION": 1, "L1_TD3P": 2, "ORPHAN_TD3P": 3, "ALU_YA5": 4,
                "ALU_YB8": 4, "SVA_E": 4, "SVA_F": 4, "SVA_TD5P": 4, "SVA_TD3P": 4,
                "PSEUDOGENE": 5, "PSEUDOGENE_DECOY": 5, "POLYA_ONLY": 6, "L1_MED_DELETION": 7,
                "L1_MED_DUPLICATION": 8, "EN_INDEPENDENT": 12, "TEMPLATED_LOCAL": 14,
                "PREMRNA_COINSERT": 15, "FOLDBACK_INVDUP_5P": 16}


def parse_types(spec):
    """'all' | 'tp' | 'artefact' | comma list of keys and/or numeric type ids."""
    if not spec:
        return []
    out = []
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        low = tok.lower()
        if low == "all":
            out += TP_TYPES + ART_TYPES
        elif low in ("tp", "tps"):
            out += TP_TYPES
        elif low in ("art", "artefact", "artefacts", "artifact", "artifacts"):
            out += ART_TYPES
        elif tok.isdigit():
            out += [k for k, v in _KEY_TYPE_ID.items() if v == int(tok)]
        elif tok.upper() in CATALOGUE:
            out.append(tok.upper())
        else:
            raise SystemExit(f"unknown insertion type '{tok}'; known: {', '.join(CATALOGUE)}")
    seen, uniq = set(), []
    for k in out:
        if k not in seen:
            seen.add(k)
            uniq.append(k)
    return uniq


def build_event(key, rng, ctx):
    return CATALOGUE[key](rng, ctx)


# ------------------------------------------------------------------------------ apply
@dataclass
class Seg:
    a0: int
    a1: int
    kind: str          # REF | INS
    r0: int = -1       # REF: window coords [r0, r1)
    r1: int = -1
    strand: str = "+"  # REF orientation of alt relative to ref
    label: str = ""

    def refpos(self, a):
        if self.strand == "+":
            return self.r0 + (a - self.a0)
        return self.r1 - 1 - (a - self.a0)


@dataclass
class Alt:
    alt: str
    segs: list
    left: int          # LEFT-clip breakpoint (window coords, original orientation)
    right: int         # RIGHT-clip breakpoint
    x_span: tuple      # alt coords [x0, x1) of the inserted sequence (original orientation)
    x_seq: str         # inserted sequence in ELEMENT sense
    tsd_seq: str
    en_motif: str      # 7-mer around the nick, element orientation (3 bp 5' of nick + 4 3')
    en_mm: int         # mismatches of the 6-mer TT|AAAA
    polya_side: str    # LEFT | RIGHT | -
    mh_seq: str = ""
    parts: list = field(default_factory=list)   # (label, x-start, x-end, info) in x_seq coords


def plant_en_motif(ref, nick, strand, motif):
    """Write the 6-mer `motif` (top-strand form TT|AAAA for a '+' insertion) at `nick`."""
    if strand == "+":
        p = nick - EN_NICK_OFFSET
        return ref[:p] + motif + ref[p + 6:]
    rcm = revcomp(motif)            # '-' insertion: top strand TTTT|AA around the nick
    p = nick - 4
    return ref[:p] + rcm + ref[p + 6:]


def find_en_sites(seq, max_mm=3, margin=3000):
    """All candidate nicks with <= max_mm mismatches to TT|AAAA: [(nick, strand, mm)]."""
    from .seqs import EN_CONSENSUS
    out = []
    rc_cons = revcomp(EN_CONSENSUS)
    n = len(seq)
    for p in range(margin, n - margin - 6):
        w = seq[p:p + 6]
        if "N" in w:
            continue
        mp = sum(1 for a, b in zip(w, EN_CONSENSUS) if a != b)
        if mp <= max_mm:
            out.append((p + EN_NICK_OFFSET, "+", mp))
        mm = sum(1 for a, b in zip(w, rc_cons) if a != b)
        if mm <= max_mm:
            out.append((p + 4, "-", mm))
    return out


def resolve_extra(ref_p, x5_bp, extra):
    """Resolve one ref-derived 5' extra against the '+'-orientation window. x5_bp is the
    window position of the 5' junction (where X starts). Returns ('REF', r0, r1, strand, label)."""
    kind = extra[0]
    if kind == "FOLDBACK":
        f, g = extra[1], extra[2]
        r1 = x5_bp - g
        return ("REF", max(0, r1 - f), r1, "-", "FOLDBACK")
    if kind == "REFCOPY":
        _, label, off, n, orient = extra
        r0 = x5_bp + off if off >= 0 else x5_bp + off - n
        r0 = max(0, min(r0, len(ref_p) - n))
        return ("REF", r0, r0 + n, orient, label)
    raise ValueError(kind)


def apply_event(ref, nick, strand, ev):
    """Apply `ev` at `nick` (window coords, original orientation) with insertion strand.
    Returns an Alt (window coords)."""
    N = len(ref)
    ref_p = ref if strand == "+" else revcomp(ref)
    s = nick if strand == "+" else N - nick
    t = ev.target_len
    if ev.target == "TSD":
        left_part_end, right_part_start, x5 = s + t, s, s + t
    elif ev.target == "DEL":
        left_part_end, right_part_start, x5 = s - t, s, s - t
    elif ev.target == "NONE":
        left_part_end, right_part_start, x5 = s, s, s
    elif ev.target == "L1DEL":
        left_part_end, right_part_start, x5 = s - t, s, s - t
    else:
        raise ValueError(ev.target)

    # assemble X pieces (5'->3')
    pieces = []
    for extra in ev.extras5:
        pieces.append(resolve_extra(ref_p, x5, extra))
    ins = ev.ins
    mh_seq = ""
    if ev.target == "L1DEL" and ev.mh:
        mh_seq = ref_p[left_part_end - ev.mh:left_part_end]
        ins = ins[ev.mh:]               # the first mh bases of the element ARE the reference mh
    pieces.append(("INS", ins, ev.key))

    segs, alt_chunks, a = [], [], 0
    segs.append(Seg(0, left_part_end, "REF", 0, left_part_end, "+", "flank5"))
    alt_chunks.append(ref_p[:left_part_end]); a = left_part_end
    x0 = a
    x_seq_parts, parts = [], []
    for p in pieces:
        if p[0] == "REF":
            _, r0, r1, st, label = p
            sq = ref_p[r0:r1] if st == "+" else revcomp(ref_p[r0:r1])
            segs.append(Seg(a, a + len(sq), "REF", r0, r1, st, label))
            parts.append((label, a - x0, a - x0 + len(sq), f"ref:{r0}-{r1}{st}"))
        else:
            sq = p[1]
            segs.append(Seg(a, a + len(sq), "INS", label=p[2]))
            off = a - x0 - (ev.mh if ev.target == "L1DEL" else 0)
            for (lab, ps, pe, inf) in ev.parts:
                ps2, pe2 = max(ps + off, a - x0), pe + off
                if pe2 > ps2:
                    parts.append((lab, ps2, pe2, inf))
        alt_chunks.append(sq); x_seq_parts.append(sq); a += len(sq)
    x1 = a
    segs.append(Seg(a, a + N - right_part_start, "REF", right_part_start, N, "+", "flank3"))
    alt_chunks.append(ref_p[right_part_start:])
    alt_p = "".join(alt_chunks)
    x_seq = "".join(x_seq_parts)

    # truth breakpoints in '+' orientation
    R_p = left_part_end              # RIGHT-clip reads end here (ref | X)
    L_p = right_part_start           # LEFT-clip reads start here (X | ref)
    tsd_seq = ref_p[s:s + t] if ev.target == "TSD" else ""
    en7 = ref_p[s - 3:s + 4] if s >= 3 else ""
    en6 = ref_p[s - 2:s + 4]
    polya_side = "LEFT" if ev.polya_len > 0 else "-"

    if strand == "+":
        return Alt(alt_p, segs, L_p, R_p, (x0, x1), x_seq, tsd_seq, en7, en_mismatches(en6),
                   polya_side, mh_seq, parts)
    # map back to original orientation
    M = len(alt_p)
    segs_o = []
    for sg in reversed(segs):
        if sg.kind == "REF":
            # alt and ref are both reverse-complemented, so a segment's orientation
            # relative to the reference is unchanged
            segs_o.append(Seg(M - sg.a1, M - sg.a0, "REF", N - sg.r1, N - sg.r0,
                              sg.strand, sg.label))
        else:
            segs_o.append(Seg(M - sg.a1, M - sg.a0, "INS", label=sg.label))
    return Alt(revcomp(alt_p), segs_o, N - R_p, N - L_p, (M - x1, M - x0), x_seq,
               ref[N - R_p:N - L_p] if ev.target == "TSD" else "",
               en7, en_mismatches(en6), "RIGHT" if ev.polya_len > 0 else "-", mh_seq, parts)


def parts_str(parts):
    return ";".join(f"{lab}:{a}-{b}" + (f"[{inf}]" if inf else "") for lab, a, b, inf in parts)
