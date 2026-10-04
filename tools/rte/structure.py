"""Element / 5'-structure / tag decisions from the element-sense read layouts.

Frame: element sense (assembly.py). The 5' junction reads REF | (extra) | element..., the 3'
junction reads ...element | (tail) | POLYA | REF.

structure (needs a read/consensus crossing the 5' junction, else 5P_UNRESOLVED):
  FULL_LENGTH          first element segment after REF is sense and starts within
                       `full_length_tolerance[class]` bp of the consensus 5' end (SVA: or inside
                       the hexamer landmark)
  TRUNCATED_5P         sense, starts further in (5' truncation point = covered_5p)
  INVERTED_5P          first element segment after REF is ANTI-sense (twin priming, Ostertag &
                       Kazazian 2001). p2 = its consensus end at the flank; p1 = start of the
                       forward part (from a read crossing the inversion point, else the lowest
                       sense coordinate); junction = p1 - p2 (>0 deletion, <0 duplication,
                       Zumalave 2026 i-del / i-dup)
  INVERTED_5P_SWITCH   as INVERTED_5P with >= 2 distinct orientation-switch points (twin priming
                       + template switch)

tags: see plans/tprt_hallmarks/SPEC.md vocabulary. Decided here: TD3P, TD3P_SOURCE=, NOVEL_SOURCE,
TD5P, TEMPLATED_LOCAL, PREMRNA_COINSERT, FOLDBACK_INVDUP_5P, CHIMERIC_ENDS, EXON_JUNCTION,
PSEUDOGENE_CANDIDATE, POLYA_ONLY/ORPHAN_TD (as element). Site-level tags (TSD_DELETION,
EN_INDEPENDENT, L1_MED_*) are added by the annotator from hallmarks.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from .transduction import known_source

DEFAULTS = {
    "full_length_tolerance": {"L1": 60, "ALU": 20, "SVA": 60},
    "truncated_3p_tolerance": 60,   # element 3' end this far before the consensus end = 3'-truncated
    "switch_cluster_bp": 50,        # orientation switches closer than this are one switch point
    "foldback_tolerance": 30,       # |b - p1| at the inversion point => fold-back inverted dup
    "min_unknown_bp": 20,           # unexplained segment length that counts as a tag candidate
    "templated_max_dist": 250,      # LOCAL template within this distance of the site
    "premrna_min_bp": 30,
}


@dataclass
class StructureCall:
    element: str = "UNKNOWN"
    structure: str = "5P_UNRESOLVED"
    tags: list = field(default_factory=list)
    detail: dict = field(default_factory=dict)
    source: object = None          # transduction.SourceCall
    j5_class: str = ""
    j3_class: str = ""
    j5_pos: int | None = None
    inv_p1: int | None = None
    three_prime_truncated: bool = False
    has_polya_3p: bool = False
    td3p_seq: str = ""

    def add(self, tag):
        if tag not in self.tags:
            self.tags.append(tag)


def _first_element_after_ref(lay):
    """(element segment, [segments between REF and it]) for a 5' layout, or (None, extras)."""
    extras = []
    for s in lay.segments[1:]:
        if s.kind == "ELEMENT":
            return s, extras
        extras.append(s)
    return None, extras


def _last_element_before_ref(lay):
    tail = []
    for s in reversed(lay.segments[:-1]):
        if s.kind == "ELEMENT":
            return s, list(reversed(tail))
        tail.append(s)
    return None, list(reversed(tail))


def classify(res, lib, ctx=None, cfg=None, novel_finder=None, premrna=None,
             pseudogene=None, legacy_class=None):
    """res: assembly.AssemblyResult. premrna: callable(seq) -> label or None. pseudogene:
    (candidate_genes, exon_junction_hits). legacy_class: annotate_v2 element_class() string."""
    c = dict(DEFAULTS)
    c.update(cfg or {})
    call = StructureCall()
    cls_ = res.element_class
    cons = res.consensus
    cend = lib.cons_end.get(cons, 0) if cons else 0
    sense_layouts = res.layouts

    five = [l for l in sense_layouts if l.ref_left() and len(l.segments) > 1]
    three = [l for l in sense_layouts if l.ref_right() and len(l.segments) > 1]
    # junction consensus first, then reads
    five.sort(key=lambda l: l.role != "JUNCTION")
    three.sort(key=lambda l: l.role != "JUNCTION")

    # ---------------------------------------------------------------- 5' junction
    j5 = None
    j5_extras = []
    votes = Counter()
    cands = []
    for lay in five:
        seg, extras = _first_element_after_ref(lay)
        if seg is None:
            continue
        if lib.cons_class.get(seg.target) != cls_ and cls_:
            # a different element class at the 5' junction (chimera candidate)
            cands.append((lay, seg, extras))
            continue
        cands.append((lay, seg, extras))
        votes[seg.strand > 0] += 3 if lay.role == "JUNCTION" else 1
    if not cands:
        # the 5' piece next to REF is a transduced/templated segment longer than a read: chain
        # through another layout that shows the same piece (same source) joined to the element
        for lay in five:
            extras = [x for x in lay.segments[1:] if x.kind in ("FLANK5P", "LOCAL", "FLANK3P")]
            if not extras:
                continue
            seg = _chain(sense_layouts, extras[0], after=True, cls_=cls_, lib=lib)
            if seg is not None:
                cands.append((lay, seg, extras))
                votes[seg.strand > 0] += 1
                break
    if cands:
        if votes:
            want = max(votes, key=votes.get)
            pick = [x for x in cands if (x[1].strand > 0) == want
                    and (not cls_ or lib.cons_class.get(x[1].target) == cls_)]
        else:
            pick = cands
        lay, j5, j5_extras = pick[0] if pick else cands[0]
        call.j5_class = lib.cons_class.get(j5.target, "")

    # ---------------------------------------------------------------- 3' junction
    j3 = None
    tail = []
    for lay in three:
        seg, t = _last_element_before_ref(lay)
        if seg is not None:
            j3, tail = seg, t
            break
    if j3 is None and three:
        tail = list(three[0].segments[:-1])
        fl = [x for x in tail if x.kind == "FLANK3P"]
        if fl:
            j3 = _chain(sense_layouts, fl[0], after=False, cls_=cls_, lib=lib)
    call.j3_class = lib.cons_class.get(j3.target, "") if j3 is not None else ""
    call.has_polya_3p = any(l.segments[-2].kind == "POLYA" for l in three if len(l.segments) >= 2)
    if j3 is not None and cls_ and lib.cons_class.get(j3.target) == cls_:
        call.three_prime_truncated = j3.strand > 0 and j3.t_en < cend - c["truncated_3p_tolerance"]
        call.detail["j3"] = j3.t_en

    # ---------------------------------------------------------------- structure
    if cls_ and j5 is not None and lib.cons_class.get(j5.target) == cls_:
        tol = c["full_length_tolerance"].get(cls_, 50)
        if j5.strand > 0:
            call.j5_pos = j5.t_st
            in_hex = cls_ == "SVA" and lib.landmark_at(j5.target, j5.t_st).lower() == "hexamer"
            call.structure = "FULL_LENGTH" if (j5.t_st <= tol or in_hex) else "TRUNCATED_5P"
            call.detail["j5"] = j5.t_st
        else:
            call.j5_pos = j5.t_en
            call.structure = "INVERTED_5P"
            _inversion(call, res, j5, lib, c)
        lm = lib.landmark_at(cons, call.j5_pos)
        if lm != ".":
            call.detail["j5_feature"] = lm
    elif cls_:
        call.structure = "5P_UNRESOLVED"

    # ---------------------------------------------------------------- segments by kind
    flank3 = []
    flank5 = []
    unknown = []
    local = []
    for lay in sense_layouts:
        for s in lay.segments:
            if s.kind == "FLANK3P":
                flank3.append(s)
            elif s.kind == "FLANK5P":
                flank5.append(s)
            elif s.kind == "UNKNOWN" and s.qlen >= c["min_unknown_bp"]:
                unknown.append((lay, s))
            elif s.kind == "LOCAL":
                local.append(s)

    # ---------------------------------------------------------------- 3' transduction
    src = known_source(flank3, lib) if flank3 else None
    td_seq = ""
    # an unexplained segment sitting between the element 3' end and the poly-A
    for lay in three:
        seg, t = _last_element_before_ref(lay)
        cand = [s for s in t if s.kind in ("UNKNOWN", "FLANK3P") and s.qlen >= c["min_unknown_bp"]]
        if cand:
            td_seq = lay.seq[cand[0].q_st:cand[-1].q_en]
            break
        if seg is None:
            # read entirely in the tag: X | POLYA | REF
            cand = [s for s in lay.segments[:-1] if s.kind in ("UNKNOWN", "FLANK3P")
                    and s.qlen >= c["min_unknown_bp"]]
            if cand:
                td_seq = lay.seq[cand[0].q_st:cand[-1].q_en]
                break
    if src is None and td_seq and novel_finder is not None:
        src = novel_finder.find(td_seq)
    premrna_lab = None
    if src is None and td_seq and premrna is not None and len(td_seq) >= c["premrna_min_bp"]:
        # host-gene sequence between element and poly-A is a pre-mRNA co-insertion, not a
        # transduction (a source flank match / credible novel source wins over this)
        premrna_lab = premrna(td_seq)
        if premrna_lab:
            call.add("PREMRNA_COINSERT")
            call.detail["premrna"] = premrna_lab
            td_seq = ""
    has_td = src is not None or (bool(td_seq) and bool(cls_))
    if has_td and cls_:
        call.add("TD3P")
    if src is not None:
        call.source = src
        call.add(f"TD3P_SOURCE={src.source_id}")
        call.detail["td_end"] = src.td_end
        if src.novel:
            call.add("NOVEL_SOURCE")
            call.detail["source_identity"] = src.identity
            if src.tier:
                call.detail["novel_tier"] = src.tier
        elif src.detail:
            call.detail["td_flank"] = src.detail
    call.td3p_seq = td_seq

    # ---------------------------------------------------------------- SVA 5' transduction
    if cls_ == "SVA":
        if any(s.strand > 0 for s in flank5) or any(
                e.kind in ("UNKNOWN", "FLANK5P") and e.qlen >= 30 for e in j5_extras):
            call.add("TD5P")
            f5 = [s for s in flank5 if s.strand > 0]
            if f5:
                call.detail["td5_source"] = f5[0].target

    # ---------------------------------------------------------------- templated / pre-mRNA
    for s in local:
        d = _local_distance(s, ctx)
        if d is None or d <= c["templated_max_dist"]:
            call.add("TEMPLATED_LOCAL")
            if d is not None:
                call.detail["templated_dist"] = d
            break
    if premrna is not None and not premrna_lab:
        for lay, s in unknown:
            if s.qlen < c["premrna_min_bp"]:
                continue
            seq = lay.seq[s.q_st:s.q_en]
            if td_seq and seq in td_seq:
                continue
            lab = premrna(seq)
            if lab:
                call.add("PREMRNA_COINSERT")
                call.detail["premrna"] = lab
                break

    # ---------------------------------------------------------------- chimeric ends
    if call.j5_class and call.j3_class and call.j5_class != call.j3_class:
        call.add("CHIMERIC_ENDS")
        call.detail["ends"] = f"{call.j5_class}/{call.j3_class}"

    # ---------------------------------------------------------------- element
    pg_genes, pg_hits = pseudogene if pseudogene else ([], [])
    nonpolya_unknown = sum(s.qlen for _, s in unknown)
    if cls_:
        call.element = {"L1": "L1", "ALU": "ALU", "SVA": "SVA"}.get(cls_, "UNKNOWN")
    elif pg_genes:
        if pg_hits:
            call.element = "PSEUDOGENE"
            call.add("EXON_JUNCTION")
            call.detail["exon_junction"] = pg_hits[0][0]
        else:
            call.add("PSEUDOGENE_CANDIDATE")
    if call.element == "UNKNOWN" and not cls_:
        if src is not None and call.has_polya_3p:
            call.element = "ORPHAN_TD"
        elif call.has_polya_3p and nonpolya_unknown < c["min_unknown_bp"] and not flank3 and not pg_genes:
            call.element = "POLYA_ONLY"
        elif legacy_class in ("non_RTE_SV", "microsatellite", "templated_insertion"):
            call.element = "NON_TPRT"
    if call.element in ("POLYA_ONLY", "ORPHAN_TD", "PSEUDOGENE", "NON_TPRT", "UNKNOWN"):
        if not cls_:
            call.structure = "5P_UNRESOLVED"
    return call


def _chain(layouts, piece, after, cls_, lib):
    """Find, in any layout, a segment of the same kind/target as `piece` directly joined
    (<= 10 bp) to an element segment of class `cls_`; return that element segment. after=True:
    piece then element (5' side); after=False: element then piece (3' side)."""
    for lay in layouts:
        segs = lay.segments
        for i, s in enumerate(segs):
            if s.kind != piece.kind or s.target != piece.target or s.strand != piece.strand:
                continue
            d = 1 if after else -1
            j = i + d
            # the element's own (short) poly-A may sit between element and transduced flank
            while 0 <= j < len(segs) and segs[j].kind == "POLYA" and segs[j].qlen <= 30:
                j += d
            if not (0 <= j < len(segs)):
                continue
            e = segs[j]
            near = segs[j - d]
            gap = (near.q_st - e.q_en) if not after else (e.q_st - near.q_en)
            if e.kind == "ELEMENT" and gap <= 10 and (not cls_ or lib.cons_class.get(e.target) == cls_):
                return e
    return None


def _local_distance(seg, ctx):
    if ctx is None or not ctx.window_seq or seg.t_st < 0:
        return None
    p = ctx.window_start + (seg.t_st + seg.t_en) // 2
    bps = [x for x in (ctx.left_bp, ctx.right_bp) if x is not None]
    if not bps:
        return None
    return min(abs(p - b) for b in bps)


def _inversion(call, res, j5, lib, c):
    """Twin-priming geometry from switch points observed inside single layouts."""
    cls_ = res.element_class
    p2 = j5.t_en
    switches = []         # (position on consensus, kind)
    inner = None          # (b, p1) of an anti->sense switch
    for lay in res.layouts:
        el = [s for s in lay.segments if s.kind == "ELEMENT" and lib.cons_class.get(s.target) == cls_]
        for a, b in zip(el, el[1:]):
            if (a.strand > 0) == (b.strand > 0):
                continue
            if b.q_st - a.q_en > 10:
                continue      # not adjacent in the read
            if a.strand < 0 and b.strand > 0:
                switches.append((b.t_st, "anti>sense"))
                if inner is None:
                    inner = (a.t_st, b.t_st)
            else:
                switches.append((a.t_en, "sense>anti"))
    sense_starts = [s for s, e, sense, _ in res.segments_on_cons if sense]
    p1 = inner[1] if inner else (min(sense_starts) if sense_starts else None)
    # cluster switch positions
    pts = sorted(set(round(p) for p, _ in switches))
    clusters = []
    for p in pts:
        if not clusters or p - clusters[-1] > c["switch_cluster_bp"]:
            clusters.append(p)
    if len(clusters) >= 2:
        call.structure = "INVERTED_5P_SWITCH"
    anti_lo = min((s for s, e, sense, _ in res.segments_on_cons if not sense), default=None)
    call.detail["inv"] = f"{anti_lo if anti_lo is not None else '?'}-{p2}"
    if p1 is not None:
        call.inv_p1 = p1
        call.detail["fwd_start"] = p1
        d = p1 - p2
        call.detail["inv_junction"] = (f"del{d}" if d > 0 else (f"dup{-d}" if d < 0 else "blunt"))
        call.detail["inv_exact"] = int(inner is not None)
    if inner is not None and abs(inner[0] - inner[1]) <= c["foldback_tolerance"]:
        call.add("FOLDBACK_INVDUP_5P")
