"""Final locus class for annotate_v2's table: the legacy element_class() of the conclusion
string, overridden by a confident tools/rte verdict where the legacy class is only a fallback.

Precedence (first match wins):
  1. legacy class is a positive call -- ALU, SVA, LINE1, processed_pseudogene, RTE_other,
     non_RTE_SV, microsatellite: KEPT (Dfam / remap evidence for the element or the SV).
  2. legacy class is a fallback -- unknown, artefact, templated_insertion -- and the RTE record
     is confident:
       a. element L1 / ALU / SVA with tprt_call TPRT or LIKELY_TPRT       -> LINE1 / ALU / SVA
       b. element PSEUDOGENE (exon-exon junction proven, tag EXON_JUNCTION) -> processed_pseudogene
       c. a 3' transduction (element ORPHAN_TD, or tag TD3P) from a CREDIBLE source -- a library
          source (TD3P_SOURCE=<id>) or a tier-A novel source; a tier-B novel source (Tubio
          polymorphic L1 position, no source sequence) is not enough -- with tprt_call TPRT or
          LIKELY_TPRT                                        -> the source's class (LINE1 / SVA)
  3. legacy class templated_insertion whose remap (hs1) source coordinate lies inside the
     UNMASKED part of a known source's 3' flank (0-15 kb downstream of the source, strand-aware)
     -> 3' transduction, the source's class, even without an RTE verdict (a templated insert
     sourced < 15 kb behind a hot L1 is its orphan transduction). The soft-masked guard is the
     same as structure.py's: a repeat inside the flank matches many loci and names no source.
  4. otherwise the legacy class.
Orphan transductions are L1-mediated (or SVA-mediated for SVA sources), so they take the
LINE1 / SVA class; there is no separate transduction class (consumers -- check_known,
compare_arms, the somatic table -- key on the existing vocabulary). The conclusion text keeps
the legacy wording and gets the RTE verdict appended.
"""
from __future__ import annotations

FALLBACK = ("unknown", "artefact", "templated_insertion")
CONFIDENT = ("TPRT", "LIKELY_TPRT")
CLASS_OF = {"L1": "LINE1", "ALU": "ALU", "SVA": "SVA"}
STRUCTURE_TEXT = {
    "FULL_LENGTH": "full-length",
    "TRUNCATED_5P": "5'-truncated",
    "INVERTED_5P": "5'-inverted (twin priming)",
    "INVERTED_5P_SWITCH": "5'-inverted with template switch",
    "5P_UNRESOLVED": "5' end unresolved",
}


def source_class(lib, sid):
    """LINE1 / SVA for a source id (library source or novel:...)."""
    if sid.startswith("novel:"):
        return "LINE1"            # the novel-source rule only accepts L1 sources
    src = lib.sources.get(sid) if lib is not None else None
    sc = (src.get("element_class") or src.get("class") or "") if src else ""
    if not sc:
        sc = "SVA" if sid.upper().startswith("SVA") else "L1"
    return CLASS_OF.get(sc, "LINE1")


def source_label(lib, sid):
    """'L1_chrX_11289796_f (Xp22.2-1, UID-50)' -- band and L1Base id when known."""
    src = lib.sources.get(sid) if lib is not None and not sid.startswith("novel:") else None
    if not src:
        return sid
    bits = []
    band = src.get("band_published") or src.get("band") or ""
    band = band.split("/")[-1] if band not in (".", "") else ""
    if band:
        bits.append(band)
    for k in ("l1base_id", "intact_id", "alt_id"):
        v = src.get(k) or ""
        if v not in ("", "."):
            bits.append(v)
            break
    return f"{sid} ({', '.join(bits)})" if bits else sid


def _credible_source(rec):
    """(source id, label kind) of a credible 3' transduction source on the record, else None."""
    sid = next((t.split("=", 1)[1] for t in rec.tags if t.startswith("TD3P_SOURCE=")), None)
    if sid is None:
        return None
    if sid.startswith("novel:"):
        tier = str(rec.detail.get("novel_tier", ""))
        if tier != "A":
            return None
    return sid


def describe(rec, lib):
    """Human verdict of an RTE record, e.g. "L1HS, 5'-inverted (twin priming), inv 5155-5499"."""
    parts = []
    if rec.element in CLASS_OF:
        parts.append(rec.consensus or rec.element)
        parts.append(STRUCTURE_TEXT.get(rec.structure, rec.structure))
        d = rec.detail
        if rec.structure.startswith("INVERTED_5P") and d.get("inv"):
            inv = f"inv {d['inv']}"
            if d.get("inv_junction") and d["inv_junction"] != "unresolved":
                inv += f" {d['inv_junction']}"
            parts.append(inv)
        elif rec.structure == "TRUNCATED_5P" and d.get("j5") not in (None, ""):
            parts.append(f"from {d['j5']}")
    elif rec.element == "PSEUDOGENE":
        parts.append(f"processed pseudogene {rec.detail.get('exon_junction', '')}".strip())
    elif rec.element == "ORPHAN_TD":
        parts.append("orphan 3' transduction")
    sid = _credible_source(rec)
    if sid is not None:
        txt = f"3' transduction from {source_label(lib, sid)}"
        if rec.detail.get("td_end") not in (None, ""):
            txt += f", ends {rec.detail['td_end']} bp into its flank"
        parts.append(txt)
    return ", ".join(p for p in parts if p)


def flank_source_at(lib, contig, pos, assembly="hs1", max_masked=0.5, pad=25):
    """Known source whose 3' flank (on `assembly`) contains contig:pos in an unmasked stretch.
    Returns (source id, offset in the flank) or None."""
    if lib is None or contig is None or pos is None:
        return None
    for sid, src in lib.sources.items():
        if (src.get("flank_3p_genome") or "hs1") != assembly:
            continue
        if src.get(f"{assembly}_chrom") != contig:
            continue
        try:
            s0, e0 = int(src[f"{assembly}_start"]), int(src[f"{assembly}_end"])
            flen = int(src.get("flank_3p_len") or 15000)
        except (KeyError, TypeError, ValueError):
            continue
        strand = src.get(f"{assembly}_strand") or src.get("strand") or ""
        if strand == "+":
            off = pos - e0
        elif strand == "-":
            off = s0 - pos
        else:
            continue
        if not 0 <= off < flen:
            continue
        flank = lib.flanks3.get(src.get("flank_3p") or sid, "")
        piece = flank[max(0, off - pad):off + pad]
        if piece and sum(ch.islower() for ch in piece) / len(piece) >= max_masked:
            continue          # inside a repeat of the flank: no source
        return sid, off
    return None


def locus_class(legacy_class, rec=None, lib=None, templated_source=None, assembly="hs1"):
    """(class, appended verdict text or '') following the module precedence.
    templated_source: (contig, pos) of the legacy distal-remap source (rule 3)."""
    if legacy_class not in FALLBACK:
        return legacy_class, ""
    if rec is not None:
        verdict = describe(rec, lib)
        if rec.element in CLASS_OF and rec.tprt_call in CONFIDENT:
            return CLASS_OF[rec.element], verdict
        if rec.element == "PSEUDOGENE" and "EXON_JUNCTION" in rec.tags:
            return "processed_pseudogene", verdict
        sid = _credible_source(rec)
        if sid is not None and rec.tprt_call in CONFIDENT and (
                rec.element == "ORPHAN_TD" or "TD3P" in rec.tags):
            return source_class(lib, sid), verdict
    if legacy_class == "templated_insertion" and templated_source is not None:
        hit = flank_source_at(lib, templated_source[0], templated_source[1], assembly)
        if hit is not None:
            sid, off = hit
            return source_class(lib, sid), (f"3' transduction from {source_label(lib, sid)}, "
                                            f"source {off} bp into its 3' flank")
    return legacy_class, ""
