"""Transparent additive TPRT point system (plans/tprt_hallmarks/SPEC.md "TPRT point system").

Every feature that fires contributes `feature:+n` / `feature:-n` to `tprt_points`; the sum is
`tprt_score`, thresholded into `tprt_call` (TPRT / LIKELY_TPRT / UNCERTAIN / ARTEFACT_LIKE).
Weights and thresholds live in WEIGHTS / THRESHOLDS and can be overridden from
CONFIG['annotate']['rte_score'] = {'weights': {...}, 'thresholds': {...}}. The defaults are
literature-anchored priors, meant to be re-fit with tools/rte/calibrate.py on simulator truth.
"""
from __future__ import annotations

from dataclasses import dataclass, field

WEIGHTS = {
    # --- target-site duplication -------------------------------------------------------------
    # TPRT leaves a TSD from the staggered EN / second-strand nicks; L1 TSDs peak at ~15 bp and
    # are mostly 7-20 bp (Szak 2002, 4-25 bp window here). Verified identity (<=1 mismatch) of the
    # duplicated copy on both flanks is required for the full points.
    "tsd_4_25": 3.0,
    "tsd_other_1_50": 1.0,          # 1-3 or 26-50 bp: still a duplication, weaker evidence
    # target-site deletions accompany a minority of L1 insertions (Gilbert 2002; PCAWG
    # Rodriguez-Martin 2020 L1-mediated deletions); a small one is weakly TPRT-compatible.
    "tsd_deletion_le20": 1.0,
    # --- poly-A ------------------------------------------------------------------------------
    # poly-A >= 10 on the strand-consistent side (the side the element's 3' end points to). A
    # true TPRT tail is 15-40 bp (annotate_v2 polya_min_len notes); >= 10 keeps margin for jitter.
    "polya_ge10": 2.0,
    "polya_5_9": 0.5,
    # --- beyond-poly-A (SPEC: extra points) ----------------------------------------------------
    # sequence recovered on the far side of the poly-A proves the read really spans the tail
    # (a slippage soft clip ends inside the homopolymer); >= 2 independent fragments doubles it.
    "beyond_polya": 2.0,
    "beyond_polya_2frag": 1.5,      # added on top of beyond_polya
    # --- L1 endonuclease motif ----------------------------------------------------------------
    # SoL1Rs are 190-fold enriched in EN-motif bins (Nam 2023, 95% CI 78.8-459); PCAWG bins by
    # mismatches to TTTT|R (Rodriguez-Martin 2020). The preference is graded, not absolute: the
    # best single site 5'-TTTTT/AA-3' accounts for only 9.6% of insertions (Flasch 2019), so a
    # 2-3 mismatch site still earns points.
    "en_0_1": 2.0,
    "en_2": 1.0,
    "en_3": 0.5,
    # --- element identity / structure ---------------------------------------------------------
    "ends_concordant": 1.0,          # 5' and 3' junction pieces from the same class
    "td_source_matches_5p": 2.0,     # 3' tag is the transduced flank of exactly the 5'-end element
    "active_identity_ge98": 1.0,     # nearest active element >= 98% identical (young source)
    # twin-priming inversion points are never in the first 590 bp of L1 (Gardner 2017 / MELT,
    # 0/298 in 1000 Genomes); one deeper in is consistent, one inside 590 bp is suspicious.
    "twin_priming_ge590": 1.0,
    "twin_priming_lt590": -1.5,
    "exon_junction": 3.0,            # pseudogene: exon-exon junction read (intron removed)
    "novel_source": 1.0,             # credible (tier A, >= 98 % to L1HS) novel transduction source
    "novel_source_tier_b": 0.5,      # tier B source (95-98 %, L1PA2 / young L1PA3): weaker
    # --- support (pooled after combine, SPEC independence rule) -------------------------------
    "junction_supported": 0.5,       # per junction with >= 2 independent fragments (max 3)
    "multi_colony": 1.0,             # seen in >= 2 colonies/samples (independent libraries)
    # --- artefact indicators ------------------------------------------------------------------
    "polya_both_sides": -4.0,        # poly-A LEFT and poly-T RIGHT: orientation-incompatible ends
    # reference homopolymer at the 3' breakpoint: the dominant PEAR-TREE discovery FP mode (bwa
    # soft-clips at reference poly-A tails; a reference-context gate flagged 88% of that FP band
    # with 0/17 TP survivors on PD44579). The EN-motif points are void in that context (the
    # homopolymer itself reads as a perfect TTTT/A site).
    "slippage_context": -4.0,
    # chimeric 5'/3' ends from different elements = library chimera (Evrony 2016 / Zumalave 2026)
    "chimeric_ends": -6.0,
    # TSD > 50 bp: Evrony 2016 checked 10 RC-seq candidates with >50 bp TSDs, all were chimeras.
    "tsd_gt50": -6.0,
    "inactive_only": -2.0,           # only old/inactive subfamily sequence (no active match)
    "foldback_clip": -3.0,           # clip = reverse complement of the adjacent reference
    "recurrence": -2.0,              # same junction signature at many unrelated loci
    "cross_sample_identical": -2.0,  # single sample + fragments identical across samples
    "no_element_no_tail": -1.0,      # nothing inserted that TPRT could explain
    # EN-independent insertion (no TSD, no tail, both ends truncated; Morrish 2002): a real
    # retrotransposition but NOT a TPRT event, so it must not read as TPRT.
    "en_independent": -2.0,
}

THRESHOLDS = {"TPRT": 7.0, "LIKELY_TPRT": 4.0, "UNCERTAIN": 0.0}


@dataclass
class ScoreInput:
    element: str = "UNKNOWN"
    structure: str = "5P_UNRESOLVED"
    tags: list = field(default_factory=list)
    tsd_len: int | None = None
    tsd_verified: bool = False
    polya_len: float = 0.0
    polya_both_sides: bool = False
    slippage: bool = False
    beyond_polya_len: int = 0
    beyond_polya_support: int = 0
    en_mismatches: int | None = None
    ends_concordant: bool = False
    td_source_matches_5p: bool = False
    element_identity: float = 0.0
    inactive_only: bool = False
    inv_p1: int | None = None
    junctions_supported: int = 0
    n_samples: int = 0
    foldback: bool = False
    recurrent: bool = False
    cross_sample_identical: bool = False
    novel_tier: str = ""


def score(si: ScoreInput, weights=None, thresholds=None):
    """Returns (score, points_string, call)."""
    w = dict(WEIGHTS)
    w.update(weights or {})
    th = dict(THRESHOLDS)
    th.update(thresholds or {})
    pts = []

    def add(name, n=None):
        v = w[name] if n is None else n
        if v:
            pts.append((name, v))

    t = si.tsd_len
    # an L1-mediated duplication (TPRT 3' end: element + strand-consistent poly-A + EN motif)
    # legitimately duplicates > 50 bp; a long-TSD chimera has no EN-nick context and one
    # fragment per junction, so it keeps the penalty
    # (the chimera signature is a single-colony event of modest size: below 150 bp, one sample,
    # the penalty stays -- simulated long-TSD chimeras are 51-150 bp, single-sample)
    l1dup = ("L1_MED_DUPLICATION" in si.tags and si.en_mismatches is not None
             and si.en_mismatches <= 2 and not si.slippage
             and t is not None and (t > 150 or si.n_samples >= 2))
    if t is not None:
        if t > 50 and not l1dup:
            add("tsd_gt50")
        elif t > 50:
            pass
        elif 4 <= t <= 25 and si.tsd_verified:
            add("tsd_4_25")
        elif t > 0 and si.tsd_verified:
            add("tsd_other_1_50")
        elif -20 <= t < 0:
            add("tsd_deletion_le20")
    if si.polya_len >= 10:
        add("polya_ge10")
    elif si.polya_len >= 5:
        add("polya_5_9")
    if si.beyond_polya_len >= 10 and si.polya_len >= 5:
        add("beyond_polya")
        if si.beyond_polya_support >= 2:
            add("beyond_polya_2frag")
    mm = si.en_mismatches
    if mm is not None and not si.slippage:
        if mm <= 1:
            add("en_0_1")
        elif mm == 2:
            add("en_2")
        elif mm == 3:
            add("en_3")
    if si.td_source_matches_5p:
        add("td_source_matches_5p")
    elif si.ends_concordant:
        add("ends_concordant")
    if si.element_identity >= 0.98:
        add("active_identity_ge98")
    if si.element == "L1" and si.structure.startswith("INVERTED_5P") and si.inv_p1 is not None:
        add("twin_priming_ge590" if si.inv_p1 >= 590 else "twin_priming_lt590")
    if "EXON_JUNCTION" in si.tags:
        add("exon_junction")
    if "EN_INDEPENDENT" in si.tags:
        add("en_independent")
    if "NOVEL_SOURCE" in si.tags:
        add("novel_source_tier_b" if si.novel_tier == "B" else "novel_source")
    if si.junctions_supported:
        add("junction_supported", w["junction_supported"] * min(3, si.junctions_supported))
    if si.n_samples >= 2:
        add("multi_colony")
    # artefacts
    if si.polya_both_sides:
        add("polya_both_sides")
    if si.slippage:
        add("slippage_context")
    if "CHIMERIC_ENDS" in si.tags:
        add("chimeric_ends")
    if si.inactive_only:
        add("inactive_only")
    if si.foldback:
        add("foldback_clip")
    if si.recurrent:
        add("recurrence")
    if si.cross_sample_identical and si.n_samples <= 1:
        add("cross_sample_identical")
    if si.element in ("UNKNOWN", "NON_TPRT") and si.polya_len < 5:
        add("no_element_no_tail")
    total = round(sum(v for _, v in pts), 2)
    if total >= th["TPRT"]:
        call = "TPRT"
    elif total >= th["LIKELY_TPRT"]:
        call = "LIKELY_TPRT"
    elif total >= th["UNCERTAIN"]:
        call = "UNCERTAIN"
    else:
        call = "ARTEFACT_LIKE"
    # hard artefact signatures cap the call regardless of the other points
    if call in ("TPRT", "LIKELY_TPRT") and any(
            n in ("chimeric_ends", "tsd_gt50", "polya_both_sides") for n, _ in pts):
        call = "UNCERTAIN"
    points = ";".join(f"{n}:{v:+g}" for n, v in pts) or "."
    return total, points, call
