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
                       Zumalave 2026 i-del / i-dup). Only the inverted piece sequenced (no sense
                       piece at all): p1 unknown, detail inv_junction=unresolved. Minimal form:
                       3' end = consensus terminus + poly-A and the 5' junction reads REF |
                       poly-T only (inverted copy of the tail): inv=polyA, p2 = consensus end
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

from .assembly import counts_as_fragment
from .sequtil import low_complexity
from .transduction import known_source

DEFAULTS = {
    "full_length_tolerance": {"L1": 60, "ALU": 20, "SVA": 60},
    "truncated_3p_tolerance": 60,   # element 3' end this far before the consensus end = 3'-truncated
    "en_independent_3p_tolerance": 20,  # ... for EN_INDEPENDENT (no poly-A): any 3' end >= 20 bp short
    "switch_cluster_bp": 50,        # orientation switches closer than this are one switch point
    "foldback_tolerance": 30,       # |b - p1| at the inversion point => fold-back inverted dup
    "min_unknown_bp": 20,           # unexplained segment length that counts as a tag candidate
    "templated_max_dist": 250,      # LOCAL template's near end within this distance of the site
    "premrna_min_bp": 30,
    # annotate round 2 (E2E over-tagging): a tag needs sequence that is long, complex, in the
    # right place and seen in >= 2 independent fragments (or a library hit)
    "td_min_bp": 30,                # unexplained piece between element 3' end and poly-A (TD3P)
    "td_min_flank_bp": 30,          # FLANK3P hit length that counts as a source-flank hit
    "td_short_flank_identity": 0.95,  # ... or a shorter (>= 20 bp) hit at >= this identity sitting
                                      # directly before the poly-A (the tail position)
    "td_max_masked_frac": 0.5,      # FLANK3P hit on >= this soft-masked (repeat) flank: no source
    "td_short_flank_max_at": 0.5,   # short tail-position hit whose flank piece is this A- (or
                                    # T-) rich is the source's own poly-A remnant region, not a tag
    "td_min_fragments": 2,          # fragments showing an UNEXPLAINED tag (a flank hit needs 1)
    "templated_min_bp": 20,
    "templated_min_identity": 0.90,
    "templated_min_fragments": 2,   # fragments (junction consensus = 1) showing the template
    "switch_min_seg": 20,           # each side of a sense>anti switch read (L1_INV_SWITCH)
    "pseudogene_full_length_tol": 15,  # 5' insert starts this close to the transcript start
    "terminal_inv_tolerance": 5,    # 3' end this close to the consensus end = terminus
    "inverted_tail_min": 10,        # poly-T run at the 5' junction (element sense) = inverted tail
    "inverted_tail_max_gap": 12,    # ... at most this far from REF
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
    three_prime_short: bool = False
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
        if not _ref_at_breakpoint(lay.segments[0], ctx):
            continue      # REF | element read from elsewhere in the window (a reference copy)
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
        call.three_prime_short = j3.strand > 0 and j3.t_en < cend - c["en_independent_3p_tolerance"]
        call.detail["j3"] = j3.t_en

    # ---------------------------------------------------------------- structure
    if cls_ and j5 is not None and lib.cons_class.get(j5.target) == cls_:
        tol = c["full_length_tolerance"].get(cls_, 50)
        if j5.strand > 0:
            call.j5_pos = j5.t_st
            in_hex = cls_ == "SVA" and lib.landmark_at(j5.target, j5.t_st).lower() == "hexamer"
            call.structure = "FULL_LENGTH" if (j5.t_st <= tol or in_hex) else "TRUNCATED_5P"
            call.detail["j5"] = j5.t_st
            sw = _sense_switch(res, j5, lib, cls_, c)
            if sw is not None:
                # twin priming + 5' switching (catalogue 13): the 5'-most piece is sense, joined
                # (in a read) to an INVERTED piece further 3' on the consensus
                call.structure = "INVERTED_5P_SWITCH"
                call.detail["switch"] = f"{sw[0]}>{sw[1]}"
        else:
            call.j5_pos = j5.t_en
            call.structure = "INVERTED_5P"
            _inversion(call, res, j5, lib, c)
        lm = lib.landmark_at(cons, call.j5_pos)
        if lm != ".":
            call.detail["j5_feature"] = lm
    elif cls_:
        call.structure = "5P_UNRESOLVED"
        if (j3 is not None and j3.strand > 0 and call.has_polya_3p
                and j3.t_en >= cend - c["terminal_inv_tolerance"] and _inverted_tail_5p(five, c)):
            # minimal twin priming: the 3' end is the consensus terminus + poly-A, and the 5'
            # junction reads REF | poly-T, i.e. an inverted copy of the tail (and possibly the
            # terminal end) -- no element piece inside the inversion was sequenced. The
            # terminal piece is required, so reference A-run slippage (POLYA_ONLY) can't get here
            call.structure = "INVERTED_5P"
            call.j5_pos = cend
            call.detail["inv"] = "polyA"
            call.detail["inv_junction"] = "unresolved"
    # SVA 5' transduction: transcription started upstream of the SVA, so the SVA itself is
    # complete at its 5' end; a 5' junction reading REF | source 5' flank (sense) is therefore a
    # FULL_LENGTH SVA even when no read joins the (long) flank to the hexamer
    if cls_ == "SVA" and call.structure == "5P_UNRESOLVED" and _td5_at_junction(five):
        call.structure = "FULL_LENGTH"
        call.detail["j5"] = "td5p"

    # ---------------------------------------------------------------- segments by kind
    flank3 = []
    flank5 = []
    unknown = []
    for lay in sense_layouts:
        for s in lay.segments:
            if s.kind == "FLANK3P":
                flank3.append(s)
            elif s.kind == "FLANK5P":
                flank5.append(s)
            elif s.kind == "UNKNOWN" and s.qlen >= c["min_unknown_bp"]:
                unknown.append((lay, s))

    # ---------------------------------------------------------------- 3' transduction
    # a source-flank hit counts when it is long enough and sits on the 3' side (no element of
    # the class after it in the read): short (20 bp) chance hits of an element end / tail junk
    # to a 15 kb flank library, and genomic pieces at the 5' end (SVA 5' transduction, local
    # templates) were the E2E TD3P_SOURCE false positives
    side3 = "LEFT" if res.strand >= 0 else "RIGHT"
    flank3_ok = []
    side3_sources = set()
    long_sources = set()
    for lay in sense_layouts:
        for i, s in enumerate(lay.segments):
            if s.kind == "FLANK3P" and not _element_after(lay, i, cls_, lib):
                if s.qlen < 60 and _explained_by_consensus(lay.seq[s.q_st:s.q_en], lib):
                    continue      # an element end that also sits (as a repeat) in some flank
                if _flank_masked_frac(lib, s) >= c["td_max_masked_frac"]:
                    continue      # the hit lies in a repeat of the flank (soft-masked): an
                                  # Alu/L1/SVA/MER copy matches many flanks and loci, it does
                                  # not identify a source (PD37590: 13/15 spurious sources)
                sid = lib.source_for_flank(s.target)
                if not _source_class_ok(lib, sid, cls_):
                    continue      # an L1 source cannot transduce behind an Alu/SVA (and v.v.)
                flank3_ok.append(s)
                nxt = lay.segments[i + 1] if i + 1 < len(lay.segments) else None
                tail_pos = nxt is not None and nxt.kind == "POLYA" and nxt.q_st - s.q_en <= 3
                short_ok = (tail_pos and s.identity >= c["td_short_flank_identity"]
                            and not _at_rich(lib.flanks3.get(s.target, "")[s.t_st:s.t_en],
                                             c["td_short_flank_max_at"]))
                if s.qlen >= c["td_min_flank_bp"] or short_ok:
                    long_sources.add(sid)
                if lay.side in (side3, "") or j5 is None:
                    side3_sources.add(sid)
    # a source needs one hit >= td_min_flank_bp (its shorter pieces -- a flank split by an
    # internal A-run -- then still count for td_end), and must be seen from the 3' junction
    # side (or the insertion has no element 5' end): mates of the element's 5' junction that
    # land in some flank are outside the insert
    flank3_ok = [s for s in flank3_ok if lib.source_for_flank(s.target) in side3_sources & long_sources]
    src = known_source(flank3_ok, lib) if flank3_ok else None
    td_seq = _unexplained_tail(sense_layouts, cls_, lib, c, call)
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
    if has_td and (cls_ or src is not None):
        call.add("TD3P")          # ORPHAN_TD (no element) is a 3' transduction too
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
        # a hit in a soft-masked part of an SVA 5' flank (the source's own hexamer / an Alu /
        # an L1 nearby) is a repeat match, not the source's unique upstream sequence
        f5 = [s for s in flank5 if s.strand > 0
              and _flank_masked_frac(lib, s, lib.flanks5) < c["td_max_masked_frac"]]
        if f5 or any(e.kind == "UNKNOWN" and e.qlen >= 30 for e in j5_extras) or any(
                e.kind == "FLANK5P" and e.qlen >= 30 and e in f5 for e in j5_extras):
            call.add("TD5P")
            if f5:
                best = Counter()
                for s5 in f5:
                    best[lib.source_for_flank(s5.target)] += s5.matches
                sid = max(best, key=best.get)
                call.detail["td5_source"] = sid
                call.add(f"TD5P_SOURCE={sid}")

    # ---------------------------------------------------------------- templated / pre-mRNA
    fb_segs = _foldback_5p(five, c)
    if fb_segs:
        call.add("FOLDBACK_INVDUP_5P")
        call.detail["foldback_frags"] = len({k for k, _ in fb_segs})
    tmpl = _local_templates(sense_layouts, ctx, c, exclude={id(x) for _, x in fb_segs})
    if tmpl.get("templated"):
        d, n, iv = tmpl["templated"]
        call.add("TEMPLATED_LOCAL")
        if d is not None:
            call.detail["templated_dist"] = d
        call.detail["templated_frags"] = n
    if tmpl.get("distal") and not premrna_lab:
        # a local template further than templated_max_dist (up to rte_wide_window): a template
        # switch onto a nearby transcript, i.e. co-inserted local pre-mRNA (Nam 2023 Fig. 4g)
        d, n, iv = tmpl["distal"]
        premrna_lab = f"local:{ctx.contig}:{iv[0]}-{iv[1]}(d={d})"
        call.add("PREMRNA_COINSERT")
        call.detail["premrna"] = premrna_lab
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
    # SVA carries an Alu-like domain right after its 5' hexamer, so an SVA whose 5' junction
    # piece lands there reads ALU/SVA: compatible, not a chimera (E2E: 4/8 SVA_TD5P/TD3P TPs
    # were tagged CHIMERIC_ENDS before this)
    compatible = {("ALU", "SVA")}
    if (call.j5_class and call.j3_class and call.j5_class != call.j3_class
            and (call.j5_class, call.j3_class) not in compatible):
        call.add("CHIMERIC_ENDS")
        call.detail["ends"] = f"{call.j5_class}/{call.j3_class}"

    # ---------------------------------------------------------------- element
    pg = tuple(pseudogene) if pseudogene else ([], [])
    pg_genes, pg_hits = pg[0], pg[1]
    pg_structure = pg[2] if len(pg) > 2 else None
    nonpolya_unknown = sum(s.qlen for _, s in unknown)
    if pg_genes and pg_hits:
        # an exon-exon junction read of a candidate gene is proof of a spliced mRNA: it wins
        # over an element class, which then comes from an Alu/L1 piece inside the mRNA (UTR
        # Alus; the E2E's parent genes were cut from repeat-rich chr22 sequence)
        call.element = "PSEUDOGENE"
        call.add("EXON_JUNCTION")
        call.detail["exon_junction"] = pg_hits[0][0]
        if cls_:
            call.detail["rte_in_mrna"] = cls_
        # RTE-structure tags read from pieces of the mRNA are meaningless here; a parent gene
        # near the site also makes the mRNA look like a local pre-mRNA template
        for t in list(call.tags):
            if t.split("=")[0] in ("TD3P", "TD3P_SOURCE", "NOVEL_SOURCE", "TD5P", "TD5P_SOURCE",
                                   "CHIMERIC_ENDS", "FOLDBACK_INVDUP_5P", "PREMRNA_COINSERT",
                                   "TEMPLATED_LOCAL"):
                call.tags.remove(t)
        call.detail.pop("premrna", None)
        call.source = None
        st = pg_structure(sense_layouts) if pg_structure is not None else None
        call.structure = st or "5P_UNRESOLVED"
        return call
    if cls_:
        call.element = {"L1": "L1", "ALU": "ALU", "SVA": "SVA"}.get(cls_, "UNKNOWN")
    elif pg_genes:
        call.add("PSEUDOGENE_CANDIDATE")
    if call.element == "UNKNOWN" and not cls_:
        if src is not None and call.has_polya_3p:
            call.element = "ORPHAN_TD"
        elif call.has_polya_3p and nonpolya_unknown < c["min_unknown_bp"] and not flank3 and not pg_genes:
            call.element = "POLYA_ONLY"
        elif legacy_class in ("non_RTE_SV", "SV_DELETION", "SV_DUPLICATION", "SV_INVERSION",
                              "microsatellite",
                              "templated_insertion"):
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
            if s.kind == "LOCAL" and (s.t_en < piece.t_st - 20 or s.t_st > piece.t_en + 20):
                continue      # a different local template
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
    else:
        # only the inverted piece (+ poly-A) was sequenced: the forward part and hence the
        # inversion point / i-del size are unknown -- say so instead of guessing
        call.detail["inv_junction"] = "unresolved"
    if inner is not None and abs(inner[0] - inner[1]) <= c["foldback_tolerance"]:
        call.add("FOLDBACK_INVDUP_5P")


def _element_after(lay, i, cls_, lib):
    """True when a class element segment follows segment i in this (element-sense) layout."""
    return any(t.kind == "ELEMENT" and (not cls_ or lib.cons_class.get(t.target) == cls_)
               for t in lay.segments[i + 1:])


def _unexplained_tail(layouts, cls_, lib, c, call):
    """The unexplained (UNKNOWN, complex, >= td_min_bp) piece between the element 3' end and the
    poly-A: in element sense ELEMENT ... X | POLYA (X directly before the poly-A), or a read lying
    entirely in the tag (X | POLYA with no element). Needs >= td_min_fragments distinct fragments
    (the junction consensus counts as one). Returns the longest such sequence or ''."""
    frags = set()
    best = ""
    best_j = ""
    for lay in layouts:
        if not counts_as_fragment(lay):
            continue
        segs = lay.segments
        for i in range(len(segs) - 1):
            x, nxt = segs[i], segs[i + 1]
            if x.kind != "UNKNOWN" or nxt.kind != "POLYA" or nxt.q_st - x.q_en > 5:
                continue
            after = segs[i + 2] if i + 2 < len(segs) else None
            if after is not None and after.kind != "REF" and not (after.kind == "UNKNOWN" and i + 3 == len(segs)):
                continue      # an A-run inside genomic sequence, not the insertion's poly-A tail
            if x.qlen < c["td_min_bp"]:
                continue
            piece = lay.seq[x.q_st:x.q_en]
            if low_complexity(piece):
                continue
            prev = [t for t in segs[:i] if t.kind == "ELEMENT"]
            if prev and cls_ and lib.cons_class.get(prev[-1].target) != cls_:
                continue
            frags.add(lay.frag_key)
            if lay.role == "JUNCTION":
                best_j = piece if len(piece) > len(best_j) else best_j
            elif len(piece) > len(best):
                best = piece
    if frags:
        call.detail["td_frags"] = len(frags)
    if len(frags) < c["td_min_fragments"]:
        return ""
    return best_j or best


def _template_abs(seg, ctx):
    """Absolute genome interval of a LOCAL segment (target 'wide' already absolute; 'site' is
    relative to the +-local_window window). None without genome coordinates."""
    if seg.t_st < 0:
        return None
    if seg.target == "wide":
        return seg.t_st, seg.t_en
    if ctx is None or not ctx.window_seq:
        return None
    return ctx.window_start + seg.t_st, ctx.window_start + seg.t_en


def _foldback_5p(five_layouts, c, max_gap=12, min_len=15):
    """Fold-back inverted duplication 5' of the site (catalogue 16): a 5' junction read reading
    REF | rc(the flank right next to the junction) | ... -- a LOCAL piece on the opposite strand
    to the flank whose template ends within `max_gap` bp of the junction, on the flank side.
    Needs >= templated_min_fragments fragments. Returns [(frag_key, segment)]."""
    hits = []
    for lay in five_layouts:
        if not counts_as_fragment(lay):
            continue
        segs = lay.segments
        if len(segs) < 2 or segs[0].kind != "REF" or segs[0].strand == 0 or segs[0].t_st < 0:
            continue
        ref = segs[0]
        x = segs[1]
        if x.kind != "LOCAL" or x.target != "site" or x.qlen < min_len or x.q_st - ref.q_en > max_gap:
            continue
        if x.strand == ref.strand or x.t_st < 0:
            continue
        if ref.strand > 0:            # junction at the REF end
            J = ref.t_en
            ok = J - max_gap <= x.t_en <= J + 2
        else:
            J = ref.t_st
            ok = J - 2 <= x.t_st <= J + max_gap
        if ok:
            hits.append((lay.frag_key, x))
    if len({k for k, _ in hits}) < c["templated_min_fragments"]:
        return []
    return hits


def _local_templates(layouts, ctx, c, exclude=()):
    """Qualifying LOCAL templates: >= templated_min_bp, >= templated_min_identity, not low
    complexity, not explained by the target site (template inside the TSD +-5 bp), seen in >=
    templated_min_fragments fragments (overlapping template intervals). Split by the distance
    of the template's NEAR end to the closest breakpoint: <= templated_max_dist -> 'templated',
    beyond -> 'distal' (pre-mRNA). Each value: (distance or None, n_fragments, interval or None)."""
    bps = [x for x in ((ctx.left_bp, ctx.right_bp) if ctx is not None else ()) if x is not None]
    cands = []
    for lay in layouts:
        if not counts_as_fragment(lay):
            continue
        for s in lay.segments:
            if s.kind != "LOCAL" or s.qlen < c["templated_min_bp"] or id(s) in exclude:
                continue
            if s.identity and s.identity < c["templated_min_identity"]:
                continue
            if low_complexity(lay.seq[s.q_st:s.q_en]):
                continue
            iv = _template_abs(s, ctx)
            d = None
            if iv is not None and bps:
                lo, hi = min(bps), max(bps)
                if iv[0] >= lo - 5 and iv[1] <= hi + 5:
                    continue          # the TSD itself (a read through the duplication)
                d = min(0 if iv[0] <= b <= iv[1] else min(abs(iv[0] - b), abs(iv[1] - b)) for b in bps)
            cands.append((iv, d, lay.frag_key))
    out = {}
    for iv, d, _ in cands:
        if iv is None:
            n = len({f for v, _, f in cands if v is None})
        else:
            n = len({f for v, _, f in cands if v is not None and v[0] < iv[1] + 20 and iv[0] < v[1] + 20})
        if n < c["templated_min_fragments"]:
            continue
        key = "templated" if d is None or d <= c["templated_max_dist"] else "distal"
        if key not in out or n > out[key][1]:
            out[key] = (d, n, iv)
    return out


def _sense_switch(res, j5, lib, cls_, c):
    """A read joining a SENSE piece to an ANTI-SENSE piece that lies further 3' on the
    consensus (sense>anti switch, <= 10 bp apart in the read), with the sense piece at or after
    the 5' junction piece: twin priming whose 5'-most part switched back to sense. Returns
    (sense end, anti end) on the consensus or None."""
    m = c["switch_min_seg"]
    for lay in res.layouts:
        el = [s for s in lay.segments if s.kind == "ELEMENT" and lib.cons_class.get(s.target) == cls_]
        for a, b in zip(el, el[1:]):
            if not (a.strand > 0 and b.strand < 0) or b.q_st - a.q_en > 10:
                continue
            if a.qlen < m or b.qlen < m:
                continue
            if a.t_st >= j5.t_st - 50 and b.t_en > a.t_en + 20:
                return a.t_en, b.t_en
    return None


def _inverted_tail_5p(five_layouts, c):
    """A 5' junction (element sense) reading REF | poly-T with no element segment after it."""
    for lay in five_layouts:
        segs = lay.segments
        if len(segs) < 2 or segs[0].kind != "REF":
            continue
        t = segs[1]
        if (t.kind == "POLYA" and t.strand < 0 and t.qlen >= c["inverted_tail_min"]
                and t.q_st - segs[0].q_en <= c["inverted_tail_max_gap"]
                and not any(x.kind == "ELEMENT" for x in segs[2:])):
            return True
    return False


def _flank_masked_frac(lib, seg, flanks=None):
    """Soft-masked (lower-case = RepeatMasker) fraction of a FLANK3P/FLANK5P hit's flank
    interval."""
    fs = (lib.flanks3 if flanks is None else flanks).get(seg.target, "")
    piece = fs[max(0, seg.t_st):max(0, seg.t_en)]
    if not piece:
        return 0.0
    return sum(ch.islower() for ch in piece) / len(piece)


def _source_class_ok(lib, sid, cls_):
    """A transduction source must be of the inserted element's class (L1 flanks behind an L1,
    SVA flanks behind an SVA; no Alu sources exist). An element-less insert (orphan
    transduction) takes any source."""
    if not cls_:
        return True
    src = lib.sources.get(sid) if hasattr(lib, "sources") else None
    sc = ((src.get("element_class") or src.get("class") or "") if isinstance(src, dict) else "")
    return not sc or sc == cls_


def _at_rich(piece, max_frac):
    piece = piece.upper()
    if not piece:
        return False
    return max(piece.count("A"), piece.count("T")) / len(piece) > max_frac


def _td5_at_junction(five_layouts):
    for lay in five_layouts:
        segs = lay.segments
        if len(segs) >= 2 and segs[0].kind == "REF":
            nxt = next((x for x in segs[1:] if x.kind != "POLYA"), None)
            if nxt is not None and nxt.kind == "FLANK5P" and nxt.strand > 0:
                return True
    return False


def _explained_by_consensus(piece, lib, max_frac=0.15):
    """A short flank hit whose sequence is also an element-consensus match (an L1 3' end / Alu
    copy inside a 15 kb source flank) says nothing about the source."""
    from .sequtil import edlib_best
    for name, seq in lib.consensus.items():
        if edlib_best(piece, seq[:lib.cons_end.get(name, len(seq))], max_frac=max_frac) is not None:
            return True
    return False


def _ref_at_breakpoint(ref, ctx, tol=30):
    """True when a read's REF piece ends at an insertion breakpoint (or carries no window
    coordinates, e.g. the junction consensus)."""
    if ref.t_st < 0 or ctx is None:
        return True
    if not ctx.window_seq:
        # genome-free local reference = right_flank + N*30 + left_flank (SiteContext)
        n = len(ctx.right_flank)
        bps = ([n] if ctx.right_flank else []) + ([n + 30] if ctx.left_flank else [])
        a, b = ref.t_st, ref.t_en
    else:
        bps = [b for b in (ctx.left_bp, ctx.right_bp) if b is not None]
        a, b = ctx.window_start + ref.t_st, ctx.window_start + ref.t_en
    if not bps:
        return True
    return any(min(abs(a - x), abs(b - x)) <= tol for x in bps)
