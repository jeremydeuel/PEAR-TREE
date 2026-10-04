#!/usr/bin/env python3
"""Classify the unexplained combined loci of the TPRT E2E (test/e2e/run_e2e.sh) by cause, from
READ-LEVEL simulator truth.

Every junction read of an unexplained locus is traced back to its source: the simulator qname
(`S<k>:<hap>:<record>:<i>[_d<c>]`) names the haplotype it was drawn from, and the read is
re-aligned (mappy sr) to THAT haplotype (hs1-derived, carries the planted events):

  planted_junction   the read spans (or is split near) a planted event junction of its haplotype
  contiguous         aligns end to end with no indel (<= 3 bp end clip, NM <= 4): the read is a
                     faithful copy of hs1, so its GRCh38 soft clip is an hs1-vs-GRCh38 difference
  slipped            indel / end clip / partial alignment on its own source (sequencing
                     homopolymer slippage, post-homopolymer phasing junk)
  artefact_molecule  read of a simulated library-artefact molecule

Per locus (priority order):
  far_pair_tp_junction     L1DEL/L1DUP geometry, one junction made of planted-junction reads,
                           paired with an unrelated breakpoint
  planted_event_reads      a junction made of planted-junction reads (displaced / unlifted event)
  slippage                 >= half of the reads slipped on their source (sequencing artefact)
  germline_rte             >= half contiguous and a junction clip hits the RTE library
                           (hs1 carries an Alu/L1/SVA that GRCh38 lacks: a TRUE germline insertion)
  germline_str             >= half contiguous, a period-1..6 reference repeat >= 10 bp touches a
                           junction (STR / VNTR length difference between the assemblies)
  germline_other           >= half contiguous, anything else (indel / divergent block)
  far_pair_unrelated       L1DEL/L1DUP whose two junctions fall in different read classes
  other
plus geometry, number of colonies with evidence, the hs1-vs-GRCh38 alignment of an 800 bp
window (indels >= 5 bp), the reference repeat at each junction and the annotate call.

Writes <out>/unexplained_classes.tsv and prints markdown tables (cause x colonies, cause x
geometry, cause x tprt_call)."""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

import mappy
import py2bit
import pysam

NAME_RE = re.compile(r"^([^:]+):(oneside_)?(\d+)-(oneside_)?(\d+)$")


def fnv1a(s: bytes) -> int:
    h = 0xcbf29ce484222325
    for b in s:
        h ^= b
        h = (h * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))


def geometry(name):
    m = NAME_RE.match(name)
    a, b = int(m.group(3)), int(m.group(5))
    if m.group(2) or m.group(4):
        return "one-sided", 0
    gap = b - a
    if 2 <= gap <= 40:
        return "TSD", gap
    if -30 <= gap < 0:
        return "TSD-deletion", gap
    if gap in (0, 1):
        return "blunt", gap
    return ("L1DEL" if gap < -30 else "L1DUP"), gap


def rep_ctx(tb, c, b):
    """Longest period-1..6 repeat starting or ending within 2 bp of junction b (0-based)."""
    w = tb.sequence(c, max(0, b - 60), b + 60).upper()
    j = min(60, b)
    best = (0, "")
    for p in range(1, 7):
        for st in range(j - 2, j + 3):
            k = st
            while k + p < len(w) and w[k + p] == w[k]:
                k += 1
            L = k + p - st
            if L >= 2 * p and L > best[0]:
                best = (L, w[st:st + p])
            k = st
            while k - p - 1 >= 0 and w[k - p - 1] == w[k - 1]:
                k -= 1
            L = st - k + p
            if L >= 2 * p and L > best[0]:
                best = (L, w[st - p:st])
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samples", type=int, required=True)
    ap.add_argument("--hg38-2bit", required=True)
    ap.add_argument("--library", required=True, help="resources/rte_library/consensus.fa")
    ap.add_argument("--names", default=None, help="default <out>/unexplained.txt (score_e2e.py)")
    a = ap.parse_args()
    O = a.out_dir
    names = [l.strip() for l in open(a.names or os.path.join(O, "unexplained.txt")) if l.strip()]
    samples = [f"S{i + 1}" for i in range(a.samples)]

    ev = defaultdict(dict)
    with gzip.open(os.path.join(O, "combine", "P1.insertions.evidence.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ev[r["insertion_id"]][r["side"]] = r
    want = defaultdict(set)
    for n in names:
        ml = {n}
        for r in ev.get(n, {}).values():
            if r["member_loci"] not in (".", ""):
                ml |= set(r["member_loci"].split(","))
        for m in ml:
            want[m].add(n)
    rows = defaultdict(list)
    for s in samples:
        with gzip.open(os.path.join(O, f"{s}.discovery.txt.gz.evidence.tsv.gz"), "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["locus"] in want and r["role"] in ("CLIP", "POLYA"):
                    r["sample"] = s
                    for n in want[r["locus"]]:
                        rows[(n, r["side"])].append(r)
    qn = {}
    for s in samples:
        bam = pysam.AlignmentFile(os.path.join(O, f"{s}.bam"))
        need = {r["frag"] for v in rows.values() for r in v if r["sample"] == s}
        seen = set()
        for n in names:
            c = n.split(":")[0]
            nums = [int(x) for x in re.findall(r"\d+", n.split(":")[1])]
            for rd in bam.fetch(c, max(0, min(nums) - 800), max(nums) + 800):
                q = rd.query_name
                if q in seen:
                    continue
                seen.add(q)
                h = f"{fnv1a(q.encode()):016x}"
                if h in need:
                    qn[(s, h)] = q
    junc = defaultdict(list)
    with open(os.path.join(O, "donor", "junctions.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            junc[r["hap"].replace(".fa", "")].append(int(r["pos"]))
    aligners = {}

    def aligner(h):
        if h not in aligners:
            aligners[h] = mappy.Aligner(os.path.join(O, "donor", f"{h}.fa"), preset="sr")
        return aligners[h]

    def read_class(qname, seq):
        hap = qname.split(":")[1]
        if hap == "mol":
            return "artefact_molecule"
        hits = list(aligner(hap).map(seq))
        if not hits:
            return "slipped"
        h = max(hits, key=lambda x: x.mlen)
        cov = (h.q_en - h.q_st) / len(seq)
        for p in junc.get(hap, []):
            if h.r_st + 5 <= p <= h.r_en - 5 or (cov < 0.9 and (abs(p - h.r_st) < 200 or abs(p - h.r_en) < 200)):
                return "planted_junction"
        if cov < 0.9 or any(op in (1, 2) for _, op in h.cigar) or h.q_st > 3 or len(seq) - h.q_en > 3 or h.NM > 4:
            return "slipped"
        return "contiguous"

    tb = py2bit.open(a.hg38_2bit)
    lib = mappy.Aligner(a.library, k=11, w=3, min_cnt=1, min_chain_score=15, min_dp_score=20, best_n=3)
    refhap = aligner("ref.hap")
    ann = {}
    ap_ = os.path.join(O, "annot", "P1.annotated.tsv")
    if os.path.exists(ap_):
        with open(ap_) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ann[r["locus"]] = r

    def outward_clip(n, sd):
        r = ev.get(n, {}).get(sd)
        if not r:
            return ""
        low = "".join(ch for ch in r["clip_consensus"] if ch.islower()).upper()
        return rc(low) if sd == "LEFT" else low

    def lib_hit(s):
        hs = [h for h in lib.map(s) if h.mlen >= 18] if len(s) >= 18 else []
        if not hs:
            return "."
        h = max(hs, key=lambda x: x.mlen)
        return f"{h.ctg}{'+' if h.strand > 0 else '-'}{h.mlen}"

    def asm_diff(c, lo, hi):
        w = tb.sequence(c, max(0, lo - 400), hi + 400).upper()
        hs = list(refhap.map(w))
        if not hs:
            return "no_hs1_hit"
        h = max(hs, key=lambda x: x.mlen)
        big = [f"{l}{'ID'[op - 1]}" for l, op in h.cigar if op in (1, 2) and l >= 5]
        return f"aligned{h.q_st}-{h.q_en}/{len(w)};indels:{','.join(big) or '-'}"

    out_rows = []
    for n in names:
        geom, gap = geometry(n)
        m = NAME_RE.match(n)
        c, pa, pb = m.group(1), int(m.group(3)), int(m.group(5))
        per = {}
        cols = set()
        for sd in ("LEFT", "RIGHT"):
            cl = Counter()
            seen = set()
            for r in rows.get((n, sd), []):
                k = (r["sample"], r["frag"])
                if k in seen or k not in qn:
                    continue
                seen.add(k)
                cl[read_class(qn[k], r["seq"])] += 1
                cols.add(r["sample"])
            per[sd] = cl
        tot = per["LEFT"] + per["RIGHT"]
        nt = sum(tot.values()) or 1
        maj = {sd: (per[sd].most_common(1)[0][0] if per[sd] else None) for sd in per}
        lc, rcl = outward_clip(n, "LEFT"), outward_clip(n, "RIGHT")
        lh, rh = lib_hit(lc), lib_hit(rcl)
        rl, rr = rep_ctx(tb, c, pa), rep_ctx(tb, c, pb)
        far = geom in ("L1DEL", "L1DUP")
        tpj = [sd for sd in per if per[sd] and per[sd]["planted_junction"] * 2 >= sum(per[sd].values())]
        if tpj and far:
            cause = "far_pair_tp_junction"
        elif tpj:
            cause = "planted_event_reads"
        elif far and maj["LEFT"] and maj["RIGHT"] and maj["LEFT"] != maj["RIGHT"]:
            cause = "far_pair_unrelated"
        elif (tot["slipped"] + tot["artefact_molecule"]) * 2 >= nt:
            cause = "slippage"
        elif tot["contiguous"] * 2 >= nt:
            if lh != "." or rh != ".":
                cause = "germline_rte"
            elif max(rl[0], rr[0]) >= 10:
                cause = "germline_str"
            else:
                cause = "germline_other"
        else:
            cause = "other"
        if far and cause.startswith("germline"):
            cause = "far_pair_unrelated"
        an = ann.get(n, {})
        out_rows.append({
            "locus": n, "cause": cause, "geometry": geom, "gap": gap, "colonies": len(cols),
            "L_reads": ",".join(f"{k}={v}" for k, v in per["LEFT"].most_common()) or ".",
            "R_reads": ",".join(f"{k}={v}" for k, v in per["RIGHT"].most_common()) or ".",
            "L_library": lh, "R_library": rh, "L_ref_repeat": f"{rl[0]}{rl[1]}", "R_ref_repeat": f"{rr[0]}{rr[1]}",
            "hs1_vs_grch38": asm_diff(c, min(pa, pb), max(pa, pb)) if abs(gap) < 3000 else "far",
            "tprt_call": an.get("tprt_call", "."), "annot_class": an.get("class", "."),
            "L_clip": lc[:60] or ".", "R_clip": rcl[:60] or ".",
        })
    keys = list(out_rows[0].keys()) if out_rows else ["locus", "cause"]
    with open(os.path.join(O, "unexplained_classes.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r[k]) for k in keys) + "\n")
    causes = ["far_pair_tp_junction", "far_pair_unrelated", "planted_event_reads", "slippage", "germline_rte",
              "germline_str", "germline_other", "other"]
    print(f"Unexplained combined loci: {len(out_rows)}")
    for title, key, vals in (("colonies", "colonies", sorted({r['colonies'] for r in out_rows})),
                             ("geometry", "geometry", ["TSD", "TSD-deletion", "blunt", "L1DEL", "L1DUP", "one-sided"]),
                             ("tprt_call", "tprt_call", ["TPRT", "LIKELY_TPRT", "UNCERTAIN", "ARTEFACT_LIKE"])):
        cnt = Counter((r["cause"], r[key]) for r in out_rows)
        print()
        print(f"| cause | " + " | ".join(f"{title} {v}" for v in vals) + " | total |")
        print("|---|" + "---|" * (len(vals) + 1))
        for cs in causes:
            n_ = sum(cnt[(cs, v)] for v in vals)
            if n_:
                print(f"| {cs} | " + " | ".join(str(cnt[(cs, v)]) for v in vals) + f" | {n_} |")
        print(f"| **total** | " + " | ".join(str(sum(cnt[(cs, v)] for cs in causes)) for v in vals)
              + f" | {len(out_rows)} |")


if __name__ == "__main__":
    main()
