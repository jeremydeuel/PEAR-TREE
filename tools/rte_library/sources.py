"""Transduction-source catalogue for the RTE library (called from build.py).

Published source lists (all GRCh37/hg19 as published) are lifted hg19 -> hg38 -> hs1, merged per
locus, classified reference / non-reference against RepeatMasker full-length young L1s, and
complemented by seed lists (all full-length L1HS on hg38, L1Base intact L1PA2/3, hs1-only
full-length L1HS, near-full-length SVA_E/F). Every source gets a strand-aware downstream flank
(element sense) extracted from hs1 and soft-masked with the hs1 RepeatMasker track.
"""
import collections
import re

import rmsk as rm
from common import revcomp, identity, band_of, cons_identity
from build import (ta_status, IntervalIndex, overlap, young_fl_l1, xlsx_rows, log,
                   FLANK_3P_L1, FLANK_3P_SVA, FLANK_5P_SVA)

SOURCE_COLUMNS = [
    "id", "element_class", "subfamily", "ta_status", "band", "band_published", "reference",
    "hs1_status", "hs1_chrom", "hs1_start", "hs1_end", "hs1_strand", "hg38_chrom", "hg38_start",
    "hg38_end", "strand", "strand_source", "hg19_published", "l1base_id", "orf_intact",
    "identity_L1HS", "canonical_pas_3p", "evidence", "seed", "n_daughters", "daughters_by_study",
    "hotness", "alt_id", "cell_culture_activity", "flank_3p", "flank_3p_len", "flank_3p_genome",
    "pas_hexamers_3p", "flank_5p", "notes", "origin"]
# `origin` (appended last, so readers keyed by column name are unaffected): germline, or somatic
# for a source that is itself a somatic insertion in one tumour (Tubio 2014 S5 status Somatic).

YOUNG_L1 = ["L1HS"] + ["L1PA%d" % i for i in range(2, 9)]
MATCH_MIN_SPAN = 4000
# published non-reference positions of one insertion differ between studies by up to ~230 bp
# (MELT single-point calls vs TraFiC / long-read breakpoints, e.g. 6p22.1: Gardner 29,919,990 vs
# Tubio 29,920,101-277 vs Nam 29,920,212, hg19); two distinct non-reference L1 insertions within
# 300 bp of each other are not expected, so one window as for reference matching
NONREF_MERGE = 300


# ----------------------------------------------------------------------------- parsing
def _td_dir(s, e, ts, te):
    """Strand of a source from where its transduced segment lies (downstream = 3')."""
    if ts >= e - 50:
        return "+", te - e
    if te <= s + 50:
        return "-", s - ts
    return None, None


def parse_published(P):
    """Returns (entries, transduction_geometry). Entries keep hg19 coordinates as published."""
    ents = []
    # --- Nam et al. 2023 (Nature 617:540) Supplementary Table 4: 276 sources, GRCh37
    rows = xlsx_rows(P["nam_st4"])
    hdr = rows[3]
    for r in rows[4:]:
        if not r[0]:
            continue
        c, s, e = str(r[0]).split("_")
        n = sum(int(r[i]) for i in range(12, len(hdr)) if isinstance(r[i], (int, float)))
        ents.append(dict(study="Nam2023", chrom=c, start=int(s), end=int(e), band=r[1],
                         ref={"REF": "yes", "NonREF": "no"}.get(r[2], "?"), strand=r[3],
                         strand_source="Nam2023",
                         subfamily=r[10] if r[10] not in ("NA", None) else "", daughters=n,
                         note=("truncating:" + str(r[11])) if r[11] not in ("NA", "not_found", None) else ""))
    # --- Rodriguez-Martin et al. 2020 (Nat Genet 52:306) Supplementary Table 5: 124, hg19
    rows = xlsx_rows(P["rm_supp"], "suppl_table5")
    for r in rows[5:]:
        if not r[1]:
            continue
        c, se = str(r[1]).split(":")
        s, e = se.split("-")
        n = sum(int(x) for x in r[2:] if isinstance(x, (int, float)))
        ents.append(dict(study="RodriguezMartin2020", chrom=c, start=int(s), end=int(e),
                         band=r[0], ref="?", strand="?", subfamily="", daughters=n, note=""))
    # strand + transduction geometry from RM Supplementary Table 2 (every partnered/orphan TD)
    rows2 = xlsx_rows(P["rm_supp"], "suppl_table2")
    ix = {k: i for i, k in enumerate(rows2[4])}
    rm_strand = collections.defaultdict(collections.Counter)
    geom = []
    for r in rows2[5:]:
        if r[ix["ins_type"]] not in ("Partnered", "Orphan"):
            continue
        sc, tc = r[ix["source_coord"]], r[ix["td_coord"]]
        try:
            c, s, e = str(sc).split("_")
            tcc, ts, te = str(tc).split("_")
            s, e, ts, te = int(s), int(e), int(ts), int(te)
        except ValueError:
            continue
        if tcc != c:
            continue
        st, dist = _td_dir(s, e, ts, te)
        if st is None:
            continue
        rm_strand[(c, s, e)][st] += 1
        geom.append(dict(study="RodriguezMartin2020", type=r[ix["ins_type"]], td_len=te - ts, distal=dist))
    for e_ in ents:
        if e_["study"] != "RodriguezMartin2020":
            continue
        for (c, s, e), cnt in rm_strand.items():
            if c == e_["chrom"] and abs(s - e_["start"]) < 50 and abs(e - e_["end"]) < 50:
                e_["strand"] = cnt.most_common(1)[0][0]
                e_["strand_source"] = "RM2020_td_geometry"
                break
    # Nam ST2 transduction geometry (stats only)
    rows = xlsx_rows(P["nam_st2"])
    ix = {k: i for i, k in enumerate(rows[3])}
    for r in rows[4:]:
        if r[ix["MEI_TYPE"]] not in ("Partnered_transduction", "Orphan_transduction"):
            continue
        try:
            c, s, e = str(r[ix["source_coordinate"]]).split("_")[:3]
            tcc, ts, te = str(r[ix["TD_coordinate"]]).split("_")[:3]
            s, e, ts, te = int(s), int(e), int(ts), int(te)
        except ValueError:
            continue
        st, dist = _td_dir(s, e, ts, te)
        if st is None or tcc != c:
            continue
        geom.append(dict(study="Nam2023", type=r[ix["MEI_TYPE"]], td_len=te - ts, distal=dist))
    # --- Gardner et al. 2017 (MELT; Genome Res 27:1916) Supplemental Table S9 B: 38 sources
    rows = xlsx_rows(P["melt_s9"], "B")
    for r in rows[5:]:
        if r[0] is None or r[1] is None:
            continue
        tub = r[5] if isinstance(r[5], (int, float)) else 0
        ents.append(dict(study="Gardner2017", chrom=str(r[0]), start=int(r[1]), end=int(r[1]) + 1,
                         band="", ref={"REF": "yes", "NONREF": "no"}.get(str(r[2]).upper(), "?"),
                         strand="?", subfamily=r[9] or "", daughters=int(r[4] or 0),
                         alt_id=r[3] or "", tubio=int(tub),
                         cell_activity=("%s" % r[6]) if r[6] not in (None, "") else "",
                         note=""))
    # --- Gardner 2017 S9 C: literature-assessed FL-L1 sources (Tubio, Beck, Brouha, Scott ...)
    rows = xlsx_rows(P["melt_s9"], "C")
    hdr = rows[4]
    for r in rows[5:]:
        if r[0] is None or r[1] is None:
            continue
        studies = []
        for i in range(6, 21):
            if r[i] in (None, ""):
                continue
            h = str(hdr[i])
            m = re.match(r"(\w[\w\-]*) et al \((\d{4})\)", h)
            studies.append("%s%s" % (m.group(1), m.group(2)) if m else ("Gardner2017" if h == "This Study" else h))
        ents.append(dict(study="Gardner2017_S9C", chrom=str(r[0]), start=int(r[1]),
                         end=int(r[1]) + 1, band="", ref="?", strand="?", subfamily="",
                         daughters=0, alt_id=r[2] or "", lit=studies,
                         note=re.sub(r"\s+", " ", str(r[21]))[:120] if len(r) > 21 and r[21] else ""))
    # --- Tubio et al. 2014 (Science 345:1251343) Table S5: 89 sources (72 germline, 17 somatic)
    if P.get("tubio"):
        ents += parse_tubio_s5(P["tubio"])
    return ents, geom


def _tubio_activity(brouha, beck):
    """Cell-culture retrotransposition activity as tabulated by Tubio S5 (fraction of the L1RP /
    L1.3 control activity); '.' = not tested, ND = tested, not detected."""
    out = []
    for lab, v, ctrl in (("Brouha2003", brouha, "L1RP"), ("Beck2010", beck, "L1.3")):
        if v in (None, ".", ""):
            continue
        out.append("%s:%s" % (lab, ("%gx%s" % (v, ctrl)) if isinstance(v, (int, float)) else v))
    return ";".join(out)


def parse_tubio_s5(path):
    """Tubio 2014 Table S5. Coordinates hg19 (verified: 22q12.1/TTC28 at chr22:29,059,272-
    29,065,303, the Table S8 bisulfite primers sit at 29,058,931-29,059,399 = its 5' end). For
    non-reference and somatic sources Position1/Position2 are the two insertion breakpoints and
    are not ordered (Position1 > Position2 for ~1/3 of rows), so they are sorted here. Total
    derived transductions is an Excel SUM in the sheet: recomputed from the per-sample columns."""
    rows = xlsx_rows(path, "TableS5")
    hdr = next(i for i, r in enumerate(rows) if r[0] == "Chromosome" and r[1] == "Position1")
    ents = []
    for r in rows[hdr + 1:]:
        if r[0] in (None, "unk") or not isinstance(r[1], (int, float)):
            continue
        p1, p2 = int(r[1]), int(r[2])
        n = sum(int(x) for x in r[10:] if isinstance(x, (int, float)))
        somatic = str(r[6]).lower() == "somatic"
        ents.append(dict(study="Tubio2014", chrom=str(r[0]), start=min(p1, p2), end=max(p1, p2),
                         band=r[4] or "", ref="?", strand={"plus": "+", "minus": "-"}.get(r[3], "?"),
                         strand_source="Tubio2014", subfamily="", daughters=n,
                         cell_activity=_tubio_activity(r[7], r[8]),
                         origin="somatic" if somatic else "germline",
                         note="Tubio2014_somatic_source" if somatic else ""))
    return ents


def parse_tubio_s3(path):
    """Tubio 2014 Table S3: every partnered / orphan transduction with its source (hg19).
    Returns dicts: type, source (chrom, start, end, strand as in S3), segment (ts, te: 0-based
    half-open hg19), td_len, distal (= distance of the segment's distal end from the source's 3'
    end). Partnered rows give one breakpoint (the distal end; the segment starts at the source's
    3' end); orphans give both. 'partnered2' (2 rows) counts as partnered."""
    rows = xlsx_rows(path, "TableS3")
    hdr = next(i for i, r in enumerate(rows) if r[0] == "Sample" and r[1] == "Chromosome")
    out = []
    for r in rows[hdr + 1:]:
        typ = {"partnered": "partnered", "partnered2": "partnered", "orphan": "orphan"}.get(r[6])
        if typ is None:
            continue
        try:
            sc, s1, s2 = str(r[7]), int(r[8]), int(r[9])
        except (TypeError, ValueError):
            continue                      # source not identified in the paper
        strand = {"plus": "+", "minus": "-"}.get(r[10])
        bps = [int(x) for x in (r[11], r[12]) if isinstance(x, (int, float))]
        if strand is None or not bps:
            continue
        s, e = min(s1, s2), max(s1, s2)
        if typ == "partnered":
            b = bps[0]
            ts, te = (e, b) if strand == "+" else (b - 1, s - 1)
        else:
            ts, te = min(bps) - 1, max(bps)
        if te <= ts:
            continue
        distal = (te - e) if strand == "+" else (s - 1 - ts)
        out.append(dict(type=typ, sample=r[0], chrom=sc, start=s, end=e, strand=strand,
                        ts=ts, te=te, td_len=te - ts, distal=distal))
    return out


def _lift_seg(lifter, c, ts, te):
    """(chrom, start, end, end_only) of a lifted segment: the interval lift, else the span of
    whichever ends lift (same chromosome) — segments next to an element absent from the target
    assembly often sit in a chain gap. end_only = True when only one end could be placed."""
    li = lifter.interval(c, ts, te, max_len_ratio=3)
    if li:
        return li[0], li[1], li[2], False
    pts = [p for p in (lifter.point(c, ts), lifter.point(c, te - 1)) if p]
    if not pts or len({p[0] for p in pts}) > 1:
        return None
    lo, hi = min(p[1] for p in pts), max(p[1] for p in pts)
    if hi - lo > 3 * (te - ts) + 50:
        return None
    return pts[0][0], lo, hi + 1, len(pts) == 1


def _seq_hit(hg38, li, strand, names, seqs, max_err=0.1):
    """Sequence fallback when the coordinate check fails (chain artefacts around an element
    that is present in only one assembly): the hg38 segment, in source sense, aligned (edlib
    infix) to the source's flank. Returns ('hit_by_sequence', offset) or None."""
    import edlib
    if li[2] - li[1] < 30:
        return None
    g = hg38.seq(li[0], li[1], li[2])
    q = g if strand == "+" else revcomp(g)
    for nm in names:
        r = edlib.align(q, seqs[nm].upper(), mode="HW", task="locations",
                        k=int(max_err * len(q)))
        if r["editDistance"] != -1:
            return "hit_by_sequence", r["locations"][0][0] + 1
    return None


def validate_tubio_s3(tds, ents, rows, fl3, lift19, lift_hs1, hg38=None, slack=100):
    """Does every Tubio S3 transduced segment fall inside flanks_3p of its source, on the
    downstream side in element sense? Segments are lifted hg19 -> hg38 -> hs1 and compared with
    the flank window (FASTA description) of the source row the matching S5 entry was merged into.
    `slack` bp are allowed at the proximal (source) end of the window: published element ends
    differ from RepeatMasker's by a few bp, and partnered segments start exactly there.
    When the coordinate check fails, the hg38 segment sequence (source sense) is aligned to the
    flank instead (`_seq_hit`). Annotates each td with `flank_hit` (hit / hit_unknown_strand /
    hit_end_only / hit_by_sequence / miss_outside_window / miss_wrong_chrom / unlifted /
    no_source) and `flank_offset` (1-based start in the flank)."""
    t5 = [e for e in ents if e["study"] == "Tubio2014" and e.get("source_id")]
    by_id = {r["id"]: r for r in rows}
    win, seqs = {}, {}
    for nm, sq, desc in fl3:
        seqs[nm] = sq
        m = re.match(r"(hs1|hg38):(\S+):(\d+)-(\d+)\(([+-])\)", desc)
        win[nm] = (m.group(1), m.group(2), int(m.group(3)) - 1, int(m.group(4)), m.group(5))
    for t in tds:
        c = t["chrom"]
        src = next((e for e in t5 if e["chrom"] == c and abs(e["start"] - t["start"]) <= 5
                    and abs(e["end"] - t["end"]) <= 5), None)
        if src is None or src["source_id"] not in by_id:
            t["flank_hit"], t["source_id"] = "no_source", "."
            continue
        row = by_id[src["source_id"]]
        t["source_id"] = row["id"]
        li = _lift_seg(lift19, "chr" + c, t["ts"], t["te"])
        if not li:
            t["flank_hit"] = "unlifted"
            continue
        fails = []
        res = None
        for nm in row["flank_3p"].split(","):
            g, wc, f0, f1, fst = win[nm]
            if g == "hs1":
                lh = _lift_seg(lift_hs1, li[0], li[1], li[2])
                if not lh:
                    fails.append("unlifted")
                    continue
                sc, a, b = lh[0], lh[1], lh[2]
            else:
                sc, a, b = li[0], li[1], li[2]
            if sc != wc:
                fails.append("miss_wrong_chrom")
                continue
            lo, hi = (f0 - slack, f1) if fst == "+" else (f0, f1 + slack)
            if lo <= a and b <= hi:
                off = max(1, (a - f0 + 1) if fst == "+" else (f1 - b + 1))
                res = ("hit" if "," not in row["flank_3p"] else "hit_unknown_strand", off)
                if li[3] or (g == "hs1" and lh[3]):
                    res = ("hit_end_only", off)
                break
            fails.append("miss_outside_window")
        if res is None and hg38 is not None:
            res = _seq_hit(hg38, li, t["strand"], row["flank_3p"].split(","), seqs)
        if res is None:
            res = (sorted(fails, key=["miss_outside_window", "miss_wrong_chrom", "unlifted"].index)[0], None)
        t["flank_hit"], t["flank_offset"] = res
    return tds


STATS_COLUMNS = ["study", "metric", "n", "median", "p90", "p95", "p99", "max", "frac_le_1kb",
                 "frac_le_5kb", "frac_le_10kb", "frac_le_15kb", "flank_hits", "flank_hit_rate"]


def transduction_stats(geom, tubio_tds=None):
    """Length / distal-end distributions per study (+ 'all'); for Tubio 2014 also per type and
    the flanks_3p validation hit rate (appended columns flank_hits = hits/tested with a source
    in the library, flank_hit_rate)."""
    import numpy as np
    out = []
    groups = [("RodriguezMartin2020", lambda t: t["study"] == "RodriguezMartin2020"),
              ("Nam2023", lambda t: t["study"] == "Nam2023")]
    if tubio_tds:
        geom = geom + [dict(t, study="Tubio2014") for t in tubio_tds]
        groups += [("Tubio2014", lambda t: t["study"] == "Tubio2014"),
                   ("Tubio2014_partnered", lambda t: t["study"] == "Tubio2014" and t["type"] == "partnered"),
                   ("Tubio2014_orphan", lambda t: t["study"] == "Tubio2014" and t["type"] == "orphan")]
    groups.append(("all", lambda t: True))
    for study, sel in groups:
        g = [t for t in geom if sel(t)]
        hits = "."
        rate = "."
        tested = [t for t in g if t.get("flank_hit") not in (None, "no_source")]
        if tested:
            k = sum(1 for t in tested if t["flank_hit"].startswith("hit"))
            hits, rate = "%d/%d" % (k, len(tested)), "%.4f" % (k / len(tested))
        for key in ("td_len", "distal"):
            v = np.array([t[key] for t in g if t[key] is not None and t[key] >= 0])
            if not len(v):
                continue
            out.append(dict(study=study, metric=key, n=len(v), median=int(np.median(v)),
                            p90=int(np.percentile(v, 90)), p95=int(np.percentile(v, 95)),
                            p99=int(np.percentile(v, 99)), max=int(v.max()),
                            frac_le_1kb="%.4f" % (v <= 1000).mean(),
                            frac_le_5kb="%.4f" % (v <= 5000).mean(),
                            frac_le_10kb="%.4f" % (v <= 10000).mean(),
                            frac_le_15kb="%.4f" % (v <= 15000).mean(),
                            flank_hits=hits, flank_hit_rate=rate))
    return out


def pas_hexamers(flank, n=5):
    res = []
    up = flank.upper()
    for hexa in ("AATAAA", "ATTAAA"):
        pos = [m.start() + 1 for m in re.finditer("(?=%s)" % hexa, up)][:n]
        res.append("%s:%s" % (hexa, ",".join(map(str, pos)) or "-"))
    return ";".join(res)


def hotness(n):
    if n >= 20:
        return "hot"
    if n >= 5:
        return "strong"
    if n >= 1:
        return "active"
    return "none_reported"


def _canon_pas(tb, c, s, e, strand):
    g = tb.seq(c, s, e)
    if strand == "+":
        end = g[-30:] + tb.seq(c, e, e + 10)
    else:
        end = revcomp(g)[-30:] + revcomp(tb.seq(c, s - 10, s))
    return "yes" if "AATAAA" in end else "no"


# ----------------------------------------------------------------------------- catalogue
def build_sources(ents, hg38, hs1, lift19, lift_hs1, young38, young_hs1, mask_hs1,
                  l1rows, l1cons, sva_all, bands38):
    # seeds use the strict full-length definition (>= 5.9 kb); reference-status matching of
    # published sources uses a looser one (>= 4 kb: some published reference sources are
    # 5.4-5.5 kb, e.g. Nam 2023 6p12.3 / 2q21.1)
    fl38 = young_fl_l1(young38, YOUNG_L1)
    fl_hs1 = young_fl_l1(young_hs1, YOUNG_L1)
    idx38 = IntervalIndex(young_fl_l1(young38, YOUNG_L1, min_span=MATCH_MIN_SPAN))
    idx_hs1 = IntervalIndex(young_fl_l1(young_hs1, YOUNG_L1, min_span=MATCH_MIN_SPAN))
    l1b_idx = IntervalIndex([dict(chrom=r["hg38_chrom"], start=r["hg38_start"] - 1,
                                  end=r["hg38_end"], row=r) for r in l1rows])
    clusters = {}
    unlifted = []
    strand_check = collections.Counter()

    for e_ in ents:
        c = e_["chrom"] if e_["chrom"].startswith("chr") else "chr" + e_["chrom"]
        s0, e0 = e_["start"] - 1, max(e_["end"], e_["start"])
        li = lift19.interval(c, s0, e0, max_len_ratio=3)
        if not li:
            unlifted.append(e_)
            continue
        hc, h0, h1, _ = li
        hits = idx38.query(hc, h0 - 300, h1 + 300)
        if hits:
            el = max(hits, key=lambda x: (overlap(h0 - 300, h1 + 300, x["start"], x["end"]),
                                          x["end"] - x["start"]))
            k = ("ref", el["chrom"], el["start"], el["end"])
        else:
            el = None
            pos = (h0 + h1) // 2
            k = next((kk for kk in clusters if kk[0] == "nonref" and kk[1] == hc
                      and abs(kk[2] - pos) <= NONREF_MERGE), ("nonref", hc, pos, pos))
        cl = clusters.setdefault(k, dict(element=el, ents=[], seed=set()))
        cl["ents"].append(e_)
        if el is not None and e_.get("strand") in ("+", "-"):
            strand_check["%s/%s/%s" % (e_["study"], e_.get("strand_source", "table"),
                                       "agree" if e_["strand"] == el["strand"] else "DISAGREE")] += 1
        if el is not None and e_.get("ref") == "no":
            strand_check["%s/listed_nonref_but_FL_L1_in_hg38" % e_["study"]] += 1
        if el is None and e_.get("ref") == "yes":
            strand_check["%s/listed_ref_but_no_FL_L1_in_hg38" % e_["study"]] += 1
    log("published entries: %d, unlifted hg19->hg38: %d, loci: %d"
        % (len(ents), len(unlifted), len(clusters)))
    for kk, v in sorted(strand_check.items()):
        log("  check %s: %d" % (kk, v))
    for e_ in unlifted:
        log("  unlifted: %s %s:%s-%s" % (e_["study"], e_["chrom"], e_["start"], e_["end"]))

    # seeds: every full-length L1HS in hg38 + L1Base intact non-L1HS
    for f in fl38:
        if f["name"] == "L1HS":
            k = ("ref", f["chrom"], f["start"], f["end"])
            clusters.setdefault(k, dict(element=f, ents=[], seed=set()))["seed"].add("young_fl_ref")
    for r in l1rows:
        if r["subfamily"] == "L1HS":
            continue
        hits = idx38.query(r["hg38_chrom"], r["hg38_start"] - 1, r["hg38_end"])
        if hits:
            f = max(hits, key=lambda x: x["end"] - x["start"])
            k = ("ref", f["chrom"], f["start"], f["end"])
            clusters.setdefault(k, dict(element=f, ents=[], seed=set()))["seed"].add("l1base_intact")

    rows = []
    hs1_used = []
    for k in sorted(clusters, key=lambda kk: (kk[1], kk[2])):
        cl = clusters[k]
        ents_, el = cl["ents"], cl["element"]
        studies = collections.OrderedDict()
        lit, subf, alt, cell, notes, bandp, hg19 = set(), set(), set(), set(), set(), set(), set()
        daughters = collections.OrderedDict()
        pub_strand = collections.Counter()
        for e_ in ents_:
            studies[e_["study"].replace("_S9C", "")] = 1
            if e_.get("daughters") and e_["study"] == "Tubio2014":
                # the same count also reaches us through Gardner S9 B ("Tubio et al. Activity",
                # a subset): never sum the two
                daughters["Tubio2014"] = max(daughters.get("Tubio2014", 0), e_["daughters"])
            elif e_.get("daughters"):
                daughters[e_["study"]] = daughters.get(e_["study"], 0) + e_["daughters"]
            if e_.get("tubio"):
                daughters["Tubio2014"] = max(daughters.get("Tubio2014", 0), e_["tubio"])
            lit.update(x for x in e_.get("lit", []) if x != "Gardner2017")
            for key, bag in (("subfamily", subf), ("alt_id", alt), ("cell_activity", cell),
                             ("note", notes), ("band", bandp)):
                if e_.get(key):
                    bag.add(str(e_[key]))
            if e_.get("strand") in ("+", "-"):
                pub_strand[e_["strand"]] += 1
            hg19.add("%s:%d-%d" % (e_["chrom"].replace("chr", ""), e_["start"], e_["end"]))
        evidence = list(studies) + sorted(l for l in lit if l not in studies)
        origin = "somatic" if any(e_.get("origin") == "somatic" for e_ in ents_) else "germline"
        if el is not None:
            reference = "yes"
            c38, s38, e38, strand = el["chrom"], el["start"], el["end"], el["strand"]
            strand_src = "rmsk_hg38"
            subfamily = el["name"]
            if pub_strand and pub_strand.most_common(1)[0][0] != strand:
                notes.add("published_strand_disagrees_with_rmsk")
            pubsub = sorted(x for x in subf if x)
            if pubsub:
                subfamily += "(" + "/".join(pubsub) + ")"
        else:
            reference = "no"
            c38, s38, e38 = k[1], k[2], k[3] + 1
            if pub_strand:
                strand, strand_src = pub_strand.most_common(1)[0][0], "published"
            else:
                strand, strand_src = ".", "unknown"
            subfamily = "/".join(sorted(subf)) or "."
        l1b = l1b_idx.query(c38, s38, e38) if reference == "yes" else []
        l1b = l1b[0]["row"] if l1b else None
        ident = canon = "."
        ta = l1b["ta_status"] if l1b else "."
        if reference == "yes":
            g = hg38.seq(c38, s38, e38)
            gs = g if strand == "+" else revcomp(g)
            ident = "%.4f" % cons_identity(gs, l1cons)
            canon = _canon_pas(hg38, c38, s38, e38, strand)
            if ta == ".":
                t = ta_status(gs)
                ta = "Ta" if t == "Ta" else ("preTa" if el["name"] == "L1HS" and t == "nonTa" else t)
        # hs1 placement
        if reference == "yes":
            li = lift_hs1.interval(c38, s38, e38)
            hits = idx_hs1.query(li[0], li[1], li[2]) if li else []
            if hits:
                hb = max(hits, key=lambda x: overlap(li[1], li[2], x["start"], x["end"]))
                hc, h0, h1, hstrand = hb["chrom"], hb["start"], hb["end"], hb["strand"]
                hs1_status = "present"
                hs1_used.append(dict(chrom=hc, start=h0, end=h1))
            else:
                # element not in hs1 (polymorphic, absent from the CHM13 haplotype): place its
                # 3' junction via the first liftable downstream base; the flank is still there
                j = lift_hs1.junction(c38, e38, +1) if strand == "+" else lift_hs1.junction(c38, s38 - 1, -1)
                if j:
                    hc, hp, flip = j[0], j[1] if strand == "+" else j[1] + 1, j[2] == "-"
                    if flip:
                        hp = j[1] + 1 if strand == "+" else j[1]
                    hstrand = strand if not flip else ("-" if strand == "+" else "+")
                    h0 = h1 = hp
                    hs1_status = "insertion_point"
                    notes.add("reference_in_hg38_absent_in_hs1")
                else:
                    hc, h0, h1, hstrand, hs1_status = ".", -1, -1, ".", "unlifted"
        else:
            p = lift_hs1.point(c38, s38)
            if p:
                hc, hp, flip = p[0], p[1], p[2] == "-"
                hstrand = strand if (strand == "." or not flip) else ("-" if strand == "+" else "+")
                hits = idx_hs1.query(hc, hp - 150, hp + 150)
                if hits:
                    hb = max(hits, key=lambda x: x["end"] - x["start"])
                    hc, h0, h1 = hb["chrom"], hb["start"], hb["end"]
                    if hstrand == ".":
                        hstrand, strand_src = hb["strand"], "rmsk_hs1"
                    elif hstrand != hb["strand"]:
                        notes.add("hs1_copy_strand_disagrees_with_published")
                    hs1_status = "present_in_hs1"
                    hs1_used.append(dict(chrom=hc, start=h0, end=h1))
                else:
                    h0 = h1 = hp
                    hs1_status = "insertion_point"
            else:
                hc, h0, h1, hstrand, hs1_status = ".", -1, -1, ".", "unlifted"
        nd = sum(daughters.values())
        sl = {"+": "f", "-": "r"}.get(hstrand if hc != "." else strand, "u")
        if hc != ".":
            sid = "L1_%s_%d_%s" % (hc, (h0 + 1) if h1 > h0 else h0, sl)
        else:
            sid = "L1_hg38_%s_%d_%s" % (c38, s38 + 1, sl)
        for e_ in ents_:
            e_["source_id"] = sid
        rows.append(dict(
            id=sid, element_class="L1", subfamily=subfamily, ta_status=ta,
            band=band_of(bands38, c38, s38), band_published="/".join(sorted(bandp)) or ".",
            reference=reference, hs1_status=hs1_status, hs1_chrom=hc,
            hs1_start=(h0 + 1) if h1 > h0 else h0, hs1_end=h1, hs1_strand=hstrand,
            hg38_chrom=c38, hg38_start=s38 + 1, hg38_end=e38, strand=strand,
            strand_source=strand_src, hg19_published=";".join(sorted(hg19)) or ".",
            l1base_id=l1b["id"] if l1b else ".", orf_intact="yes" if l1b else ".",
            identity_L1HS=ident, canonical_pas_3p=canon,
            evidence=";".join(evidence) or ".",
            seed=";".join(sorted(cl["seed"]) + (["published"] if ents_ else [])),
            n_daughters=nd,
            daughters_by_study=";".join("%s:%d" % kv for kv in daughters.items()) or ".",
            hotness=hotness(nd) if ents_ else "candidate",
            alt_id=";".join(sorted(alt)) or ".", cell_culture_activity=";".join(sorted(cell)) or ".",
            notes=";".join(sorted((n if len(n) <= 200 else n[:197] + "...") for n in notes)) or ".",
            origin=origin))

    # hs1-only full-length L1HS (T2T-resolved regions or present only in the CHM13 haplotype)
    lifted = []
    for f in fl38:
        li = lift_hs1.interval(f["chrom"], f["start"], f["end"])
        if li:
            lifted.append(dict(chrom=li[0], start=li[1], end=li[2]))
    lidx = IntervalIndex(lifted + hs1_used)
    n_hs1_only = 0
    for f in fl_hs1:
        if f["name"] != "L1HS" or lidx.query(f["chrom"], f["start"], f["end"]):
            continue
        n_hs1_only += 1
        g = hs1.seq(f["chrom"], f["start"], f["end"])
        gs = g if f["strand"] == "+" else revcomp(g)
        t = ta_status(gs)
        rows.append(dict(
            id="L1_%s_%d_%s" % (f["chrom"], f["start"] + 1, "f" if f["strand"] == "+" else "r"),
            element_class="L1", subfamily="L1HS", ta_status="Ta" if t == "Ta" else ("preTa" if t == "nonTa" else t),
            band=".", band_published=".", reference="hs1_only", hs1_status="present",
            hs1_chrom=f["chrom"], hs1_start=f["start"] + 1, hs1_end=f["end"], hs1_strand=f["strand"],
            hg38_chrom=".", hg38_start=-1, hg38_end=-1, strand=f["strand"], strand_source="rmsk_hs1",
            hg19_published=".", l1base_id=".", orf_intact=".",
            identity_L1HS="%.4f" % cons_identity(gs, l1cons),
            canonical_pas_3p=_canon_pas(hs1, f["chrom"], f["start"], f["end"], f["strand"]),
            evidence=".", seed="hs1_only_young_fl", n_daughters=0, daughters_by_study=".",
            hotness="candidate", alt_id=".", cell_culture_activity=".", notes=".", origin="germline"))
    log("hs1-only full-length L1HS added: %d" % n_hs1_only)

    # SVA sources: near-full-length SVA_E / SVA_F (human-specific, youngest; Wang et al. 2005)
    for rid, s, r in sva_all:
        if r["subfamily"] not in ("SVA_E", "SVA_F"):
            continue
        on_hs1 = r["hs1_chrom"] != "."
        sl = {"+": "f", "-": "r"}.get(r["hs1_strand"] if on_hs1 else r["strand"], "u")
        rows.append(dict(
            id=("SVA_%s_%s_%s" % (r["hs1_chrom"], r["hs1_start"], sl)) if on_hs1
            else "SVA_hg38_%s_%d_%s" % (r["hg38_chrom"], r["hg38_start"], sl),
            element_class="SVA", subfamily=r["subfamily"], ta_status=".",
            band=band_of(bands38, r["hg38_chrom"], r["hg38_start"]), band_published=".",
            reference="yes", hs1_status="present" if on_hs1 else "unlifted",
            hs1_chrom=r["hs1_chrom"], hs1_start=r["hs1_start"], hs1_end=r["hs1_end"],
            hs1_strand=r["hs1_strand"], hg38_chrom=r["hg38_chrom"], hg38_start=r["hg38_start"],
            hg38_end=r["hg38_end"], strand=r["strand"], strand_source="rmsk_hg38",
            hg19_published=".", l1base_id=".", orf_intact=".", identity_L1HS=".",
            canonical_pas_3p=".", evidence=".", seed="young_fl_ref_sva", n_daughters=0,
            daughters_by_study=".", hotness="candidate", alt_id=".", cell_culture_activity=".",
            notes=".", origin="germline"))

    seen, uniq = set(), []
    for r in rows:
        if r["id"] in seen:
            log("  duplicate source id dropped: %s (%s:%s)" % (r["id"], r["hg38_chrom"], r["hg38_start"]))
            continue
        seen.add(r["id"])
        uniq.append(r)
    rows = uniq
    fl3, fl5 = extract_flanks(rows, hg38, hs1, mask_hs1)
    return rows, fl3, fl5


def flank_window(start1, end, strand, length, side="3p", point=False):
    """0-based half-open genome window of the flank and the orientation to report it in.
    start1/end: 1-based inclusive element coords (or, for an insertion point, both = the
    0-based junction offset). 3p = downstream in element sense; 5p = upstream."""
    a, b = (start1, start1) if point else (start1 - 1, end)
    down = (side == "3p") == (strand == "+")
    if down:
        return max(0, b), b + length
    return max(0, a - length), a


def extract_flanks(rows, hg38, hs1, mask_hs1):
    fl3, fl5 = [], []
    for row in rows:
        L = FLANK_3P_L1 if row["element_class"] == "L1" else FLANK_3P_SVA
        if row["hs1_chrom"] != ".":
            gname, tb, c, st, mk = "hs1", hs1, row["hs1_chrom"], row["hs1_strand"], mask_hs1
            s1, e1 = row["hs1_start"], row["hs1_end"]
        else:
            gname, tb, c, st, mk = "hg38", hg38, row["hg38_chrom"], row["strand"], None
            s1, e1 = row["hg38_start"], row["hg38_end"]
        point = row["hs1_status"] == "insertion_point" and gname == "hs1"
        names, pas = [], []
        for strand_try in ([st] if st in ("+", "-") else ["+", "-"]):
            f0, f1 = flank_window(s1, e1, strand_try, L, "3p", point)
            g = tb.seq(c, f0, f1)
            if mk is not None:
                g = rm.softmask(g, c, f0, mk)
            seq = g if strand_try == "+" else revcomp(g)
            nm = row["id"] if st in ("+", "-") else "%s/%s" % (row["id"], strand_try)
            fl3.append((nm, seq, "%s:%s:%d-%d(%s) 3p_flank len=%d" % (gname, c, f0 + 1, f1, strand_try, len(seq))))
            names.append(nm)
            pas.append(pas_hexamers(seq))
        row["flank_3p"] = ",".join(names)
        row["flank_3p_len"] = L
        row["flank_3p_genome"] = gname
        row["pas_hexamers_3p"] = "|".join(pas)
        row["flank_5p"] = "."
        if row["element_class"] == "SVA" and st in ("+", "-"):
            f0, f1 = flank_window(s1, e1, st, FLANK_5P_SVA, "5p")
            g = tb.seq(c, f0, f1)
            if mk is not None:
                g = rm.softmask(g, c, f0, mk)
            seq = g if st == "+" else revcomp(g)
            fl5.append((row["id"], seq, "%s:%s:%d-%d(%s) 5p_flank len=%d sense; ends at the SVA 5' end"
                        % (gname, c, f0 + 1, f1, st, len(seq))))
            row["flank_5p"] = row["id"]
    return fl3, fl5


# ----------------------------------------------------------------------------- polymorphic L1s
POLY_COLUMNS = ["id", "hg19_chrom", "hg19_pos_positive", "hg19_pos_negative", "n_samples",
                "hg38_chrom", "hg38_pos", "hs1_chrom", "hs1_pos", "fl_l1_hg38", "fl_l1_hs1",
                "identity_L1HS", "source_id", "status"]


def polymorphic_candidates(path, lift19, lift_hs1, young38, young_hs1, rows, hs1, l1cons, win=300):
    """Tubio 2014 Table S7: 1,478 putative polymorphic (non-reference vs hg19) L1 insertions that
    TraFiC found in the 244 matched normals. TraFiC reports two breakpoints, no strand and no
    length, so these are *not* sources: they are kept as positions for the novel-source rule
    (tier B: an orphan/partnered segment whose upstream window holds one of these positions).
    Each is lifted hg19 -> hg38 -> hs1 (hs1_pos = 1-based position between the two breakpoints)
    and classified (first match wins):
      in_library            a non-reference or hs1-only library source lies within `win` bp
                            (source_id)
      reference_l1_at_site  a young (L1HS/L1PA2-8) >= 4 kb L1 is at the site in hg38: the call
                            sits on an L1 that is already in the reference (nested / misplaced
                            call), so it says nothing about a new full-length copy
      resolved_in_hs1       such an L1 is at the site only in hs1: the polymorphic insertion is
                            present in CHM13, so its length and identity_L1HS are known
      no_length_info        everything else (the large majority)."""
    rs = xlsx_rows(path, "TableS7")
    hdr = next(i for i, r in enumerate(rs) if r[0] == "Chromosome" and r[1] == "Positive breakpoint")
    idx38 = IntervalIndex(young_fl_l1(young38, YOUNG_L1, min_span=MATCH_MIN_SPAN))
    idxh = IntervalIndex(young_fl_l1(young_hs1, YOUNG_L1, min_span=MATCH_MIN_SPAN))
    nonref = [r for r in rows if r["element_class"] == "L1" and r["reference"] in ("no", "hs1_only")]
    src38 = IntervalIndex([dict(chrom=r["hg38_chrom"], start=int(r["hg38_start"]) - 1 - win,
                                end=int(r["hg38_end"]) + win, id=r["id"])
                           for r in nonref if r["hg38_chrom"] != "."])
    srch = IntervalIndex([dict(chrom=r["hs1_chrom"], start=int(r["hs1_start"]) - 1 - win,
                               end=int(r["hs1_end"]) + win, id=r["id"])
                          for r in nonref if r["hs1_chrom"] != "."])
    out = []
    for r in rs[hdr + 1:]:
        if r[0] is None or not isinstance(r[1], (int, float)):
            continue
        c, a, b = str(r[0]), int(r[1]), int(r[2])
        mid = (a + b) // 2
        p38 = lift19.point("chr" + c, mid - 1)
        c38, q38 = (p38[0], p38[1]) if p38 else (".", -1)
        ph = lift_hs1.point(c38, q38) if p38 else None
        ch, qh = (ph[0], ph[1]) if ph else (".", -1)
        fl38 = idx38.query(c38, q38 - win, q38 + win) if p38 else []
        flh = idxh.query(ch, qh - win, qh + win) if ph else []
        src = (src38.query(c38, q38, q38 + 1) if p38 else []) or (srch.query(ch, qh, qh + 1) if ph else [])
        ident = "."
        if flh:
            f = max(flh, key=lambda x: x["end"] - x["start"])
            g = hs1.seq(f["chrom"], f["start"], f["end"])
            ident = "%.4f" % cons_identity(g if f["strand"] == "+" else revcomp(g), l1cons)

        def fmt(lst):
            if not lst:
                return "."
            f = max(lst, key=lambda x: x["end"] - x["start"])
            return "%s:%d-%d(%s)%s:%d" % (f["chrom"], f["start"] + 1, f["end"], f["strand"], f["name"],
                                          f["end"] - f["start"])
        status = ("in_library" if src else "reference_l1_at_site" if fl38 else
                  "resolved_in_hs1" if flh else "no_length_info")
        out.append(dict(id="PL1_%s_%d" % (ch, qh + 1) if ph else "PL1_hg19_chr%s_%d" % (c, mid),
                        hg19_chrom=c, hg19_pos_positive=a, hg19_pos_negative=b, n_samples=r[4],
                        hg38_chrom=c38, hg38_pos=q38 + 1 if p38 else -1, hs1_chrom=ch,
                        hs1_pos=qh + 1 if ph else -1, fl_l1_hg38=fmt(fl38), fl_l1_hs1=fmt(flh),
                        identity_L1HS=ident, source_id=src[0]["id"] if src else ".", status=status))
    log("Tubio S7 polymorphic L1s: %d (%s)" % (len(out), dict(collections.Counter(o["status"] for o in out))))
    return out
