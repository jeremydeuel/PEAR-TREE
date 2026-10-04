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
    "pas_hexamers_3p", "flank_5p", "notes"]

YOUNG_L1 = ["L1HS"] + ["L1PA%d" % i for i in range(2, 9)]
MATCH_MIN_SPAN = 4000


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
    return ents, geom


def transduction_stats(geom):
    import numpy as np
    out = []
    for study in ("RodriguezMartin2020", "Nam2023", "all"):
        g = [t for t in geom if study == "all" or t["study"] == study]
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
                            frac_le_15kb="%.4f" % (v <= 15000).mean()))
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
                      and abs(kk[2] - pos) <= 100), ("nonref", hc, pos, pos))
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
            if e_.get("daughters"):
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
            notes=";".join(sorted((n if len(n) <= 200 else n[:197] + "...") for n in notes)) or "."))

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
            hotness="candidate", alt_id=".", cell_culture_activity=".", notes="."))
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
            notes="."))

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
