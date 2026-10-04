#!/usr/bin/env python3
"""Score the TPRT-hallmark E2E run (test/e2e/run_e2e.sh) per simulated insertion type.

Per variant (truth `variant`, role TP / ARTEFACT):
  disc      discovery recall: some colony's discovery call has both breakpoints within --window
  ev        the insertion reached combine's evidence evaluation (insertions.evidence.tsv.gz)
  pooled    ... and BOTH junctions have >= 2 independent fragments pooled over colonies
  comb      the insertion survives combine (combined.txt.gz)
  elem/struct/strand   annotate (tools/rte) element / structure / strand correct, of `comb`
  tagR      recall of the truth tags (names before '='), of `comb`; `src` exact TD3P_SOURCE id
  beyond    annotate recovered >= 10 bp beyond the poly-A (sequence across the tail), of `comb`
  calls     tprt_call histogram T/L/U/A (TPRT / LIKELY_TPRT / UNCERTAIN / ARTEFACT_LIKE)
Truth strand / poly-A side are converted to the discovery reference orientation (the hs1
event window can map reverse onto GRCh38).

Also: within-sample duplicate collapse (qname-truth: simulate_reads.py names PCR copies
`<name>_d<k>`; sidecar frag = FNV-1a(qname)), SHORT-only rescued junctions, unexplained calls,
and TPRT-score separation (AUC) on two halves (fit / held-out eval, split by event id parity;
unexplained loci split by a hash of the locus). Writes calibrate.py truth files
`<calib-truth>.fit.tsv` / `.eval.tsv`.
"""
import argparse
import csv
import gzip
import hashlib
import os
import re
import sys
from collections import Counter, defaultdict

import pysam

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "test", "fullstack"))
sys.path.insert(0, REPO)
from score_types import lift  # noqa: E402
from tools.rte.calibrate import auc  # noqa: E402

NAME_RE = re.compile(r"^([^:]+):(polyA_|disc_|oneside_)?(\d+)-(polyA_|disc_|oneside_)?(\d+)$")
EV_ROLES = ("CLIP", "POLYA", "DISC", "SPAN", "SHORT")


def fnv1a(s: bytes) -> int:
    h = 0xcbf29ce484222325
    for b in s:
        h ^= b
        h = (h * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h


def parse_name(n):
    """(contig, a, b, one_sided). One-sided loci (`oneside_`/`polyA_`/`disc_` token) carry one
    real breakpoint."""
    m = NAME_RE.match(n)
    if not m:
        return None
    return (m.group(1), int(m.group(3)), int(m.group(5)), bool(m.group(2) or m.group(4)))


def read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def disc_calls(path):
    out = set()
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0 and line.startswith("@"):
                n = line[1:].split(":LEFT:")[0].split(":RIGHT:")[0]
                p = parse_name(n)
                if p:
                    out.add(p)
    return out


def combined_names(path):
    out = set()
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0 and line.startswith("@"):
                out.add(line[1:].strip()[:-2])
    return out


def tagnames(s):
    if not s or s == ".":
        return set()
    return {t.split("=", 1)[0] for t in s.split(",") if t}


def detail(s):
    d = {}
    for kv in (s or "").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samples", type=int, required=True)
    ap.add_argument("--window", type=int, default=30)
    ap.add_argument("--report", required=True, help="markdown tables")
    ap.add_argument("--calib-truth", required=True, help="prefix for calibrate.py truth files")
    a = ap.parse_args()
    O = a.out_dir
    truth = read_tsv(os.path.join(O, "donor", "truth_types_hs1.tsv"))
    aln = lift(os.path.join(O, "flanks.bam"), 20)
    for r in truth:
        l, rr = aln.get((r["id"], "L")), aln.get((r["id"], "R"))
        span = abs(int(r["tsd"]))
        r["c"] = r["lo"] = r["hi"] = None
        r["flip"] = False
        if l and rr and l[0] == rr[0] and not l[3] and not rr[3]:
            r["c"], r["lo"], r["hi"] = l[0], l[2], rr[1]
        elif l and rr and l[0] == rr[0] and l[3] and rr[3]:
            r["c"], r["lo"], r["hi"], r["flip"] = l[0], rr[2], l[1], True
        elif l and not l[3]:
            r["c"], r["lo"], r["hi"] = l[0], l[2], l[2] + span
        elif rr and not rr[3]:
            r["c"], r["lo"], r["hi"] = rr[0], rr[1] - span, rr[1]
        st = r.get("strand", ".")
        if r["flip"] and st in "+-":
            st = "+" if st == "-" else "-"
        r["strand_ref"] = st

    def match(name_tuple, r):
        """True when both breakpoints agree; for a one-sided locus its one real breakpoint
        must agree with either truth breakpoint."""
        if name_tuple is None or r["c"] is None or name_tuple[0] != r["c"]:
            return False
        lo, hi = sorted((r["lo"], r["hi"]))
        if name_tuple[3]:
            return any(abs(x - y) <= a.window for x in name_tuple[1:3] for y in (lo, hi))
        cl, ch = sorted(name_tuple[1:3])
        return abs(cl - lo) <= a.window and abs(ch - hi) <= a.window

    samples = [f"S{i + 1}" for i in range(a.samples)]
    dcalls = {s: disc_calls(os.path.join(O, f"{s}.discovery.txt.gz")) for s in samples}
    ev_rows = read_tsv(os.path.join(O, "combine", "P1.insertions.evidence.tsv.gz"))
    ev = defaultdict(dict)
    for row in ev_rows:
        ev[row["insertion_id"]][row["side"]] = row
    comb = combined_names(os.path.join(O, "combine", "P1.combined.txt.gz"))
    annot = {row["locus"]: row for row in read_tsv(os.path.join(O, "annot", "P1.annotated.tsv"))}

    # ---------------------------------------------------------------- per event
    explained_comb = set()
    explained_disc = set()
    for r in truth:
        r["disc"] = any(match(c, r) for s in samples for c in dcalls[s])
        r["disc1s"] = r["disc"] and not any(match(c, r) for s in samples for c in dcalls[s] if not c[3])
        for s in samples:
            for c in dcalls[s]:
                if match(c, r):
                    explained_disc.add((s, c))
        # two-sided names first (an event may also have a one-sided duplicate)
        cand = sorted((n for n in ev if match(parse_name(n), r)), key=lambda n: (parse_name(n)[3], n))
        r["ev_name"] = cand[0] if cand else None
        # every junction the locus has (a one-sided locus has one real junction)
        r["pooled"] = bool(cand) and any(
            ev[n] and all(row.get("supported") == "1" for row in ev[n].values()) for n in cand)
        cc = sorted((n for n in comb if match(parse_name(n), r)), key=lambda n: (parse_name(n)[3], n))
        r["comb_name"] = cc[0] if cc else None
        explained_comb.update(cc)
        r["ann"] = annot.get(r["comb_name"]) if r["comb_name"] else None

    # ---------------------------------------------------------------- per variant table
    var = defaultdict(lambda: Counter())
    roles = {}
    calls_by_role = defaultdict(Counter)
    scores = defaultdict(lambda: {"TP": [], "ARTEFACT": []})
    for r in truth:
        v = var[r["variant"]]
        roles[r["variant"]] = r["role"]
        v["n"] += 1
        if r["c"] is None:
            continue
        v["lifted"] += 1
        v["disc"] += r["disc"]
        v["disc1s"] += r["disc1s"]
        v["comb1s"] += r["comb_name"] is not None and parse_name(r["comb_name"])[3]
        v["ev"] += r["ev_name"] is not None
        v["pooled"] += r["pooled"]
        v["comb"] += r["comb_name"] is not None
        an = r["ann"]
        if an is None:
            continue
        v["ann"] += 1
        v["elem"] += an["element"] == r["element"]
        v["struct"] += an["structure"] == r["structure"]
        d = detail(an.get("rte_detail"))
        if r["strand_ref"] in "+-":
            v["strand_n"] += 1
            v["strand"] += d.get("strand") == r["strand_ref"]
        tt = tagnames(r["tags"])
        ct = tagnames(an["tags"])
        v["tag_n"] += len(tt)
        v["tag_hit"] += len(tt & ct)
        v["tag_extra"] += len(ct - tt)
        src = [t for t in r["tags"].split(",") if t.startswith("TD3P_SOURCE=")]
        if src:
            v["src_n"] += 1
            v["src"] += src[0] in an["tags"].split(",")
        if float(r.get("polya_len") or 0) > 0:
            v["polya_n"] += 1
            bp = an.get("beyond_polya", ".")
            v["beyond"] += bp not in (".", "") and len(bp) >= 10
            v["beyond2"] += bp not in (".", "") and len(bp) >= 10 and int(d.get("beyond_polya_support", 0)) >= 2
        call = an["tprt_call"]
        v["call_" + call] += 1
        calls_by_role[r["role"]][call] += 1
        half = "fit" if int(r["id"]) % 2 == 0 else "eval"
        try:
            scores[half][r["role"]].append(float(an["tprt_score"]))
        except ValueError:
            pass
    # unexplained combined calls (organic FPs) count as artefacts for the score separation
    unexplained = sorted(n for n in comb if n not in explained_comb)
    with open(os.path.join(O, "unexplained.txt"), "w") as fh:
        fh.writelines(n + "\n" for n in unexplained)
    for n in unexplained:
        an = annot.get(n)
        if an is None:
            continue
        half = "fit" if int(hashlib.md5(n.encode()).hexdigest(), 16) % 2 == 0 else "eval"
        calls_by_role["UNEXPLAINED"][an["tprt_call"]] += 1
        try:
            scores[half + "_unexpl"]["ARTEFACT"].append(float(an["tprt_score"]))
        except ValueError:
            pass

    # ---------------------------------------------------------------- calibrate.py truth
    for half in ("fit", "eval"):
        with open(f"{a.calib_truth}.{half}.tsv", "w") as fh:
            fh.write("insertion_id\tvariant\telement\tstructure\ttags\trole\n")
            for r in truth:
                if r["comb_name"] and ("fit" if int(r["id"]) % 2 == 0 else "eval") == half:
                    fh.write(f"{r['comb_name']}\t{r['variant']}\t{r['element']}\t{r['structure']}\t{r['tags']}\t{r['role']}\n")
            for n in unexplained:
                if ("fit" if int(hashlib.md5(n.encode()).hexdigest(), 16) % 2 == 0 else "eval") == half:
                    fh.write(f"{n}\tUNEXPLAINED\tUNKNOWN\t5P_UNRESOLVED\t.\tARTEFACT\n")

    # ---------------------------------------------------------------- dedup truth
    dup = Counter()
    want = defaultdict(set)      # sample -> frag hashes
    fa_frags = defaultdict(lambda: defaultdict(set))   # (ins, side) -> sample -> {frag}
    with gzip.open(os.path.join(O, "combine", "P1.insertions.reads.fa.gz"), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                ins, side, role, smp, frag, r12 = line[1:].strip().rsplit("|", 5)
                if role in EV_ROLES:
                    fa_frags[(ins, side)][smp].add(frag)
                    want[smp].add(frag)
    qmap = {}
    for smp in samples:
        loci = {parse_name(k[0]) for k in fa_frags}
        bam = pysam.AlignmentFile(os.path.join(O, f"{smp}.bam"))
        sm = smp + ".discovery"
        need = want.get(sm, set())
        seen = set()
        for loc in loci:
            if loc is None:
                continue
            lo, hi = min(loc[1:3]), max(loc[1:3])
            for rd in bam.fetch(loc[0], max(0, lo - 1500), hi + 1500):
                q = rd.query_name
                if q in seen:
                    continue
                seen.add(q)
                h = f"{fnv1a(q.encode()):016x}"
                if h in need:
                    qmap[(sm, h)] = re.sub(r"_d\d+$", "", q)
    for key, per in fa_frags.items():
        row = ev.get(key[0], {}).get(key[1])
        if row is None or row.get("n_independent") in (None, "", "NA"):
            continue
        frags = sum(len(v) for v in per.values())
        mols = sum(len({qmap.get((s, f), f) for f in v}) for s, v in per.items())
        unresolved = sum(1 for s, v in per.items() for f in v if (s, f) not in qmap)
        n_ind = int(row["n_independent"])
        dup["junctions"] += 1
        dup["fragments"] += frags
        dup["true_dup_fragments"] += frags - mols
        dup["reported_duplicates"] += int(row.get("n_duplicates") or 0)
        dup["dup_coord"] += int(row.get("n_dup_coord") or 0)
        dup["dup_seq"] += int(row.get("n_dup_seq") or 0)
        dup["unresolved_frags"] += unresolved
        dup["exact"] += n_ind == mols
        dup["over_collapsed"] += n_ind < mols
        dup["under_collapsed"] += n_ind > mols
        dup["missed_dups"] += max(0, n_ind - mols)
        dup["false_merges"] += max(0, mols - n_ind)
    short_only = sum(1 for row in ev_rows if row.get("n_independent_no_short") not in (None, "")
                     and int(row["n_independent"] or 0) >= 2 > int(row["n_independent_no_short"] or 0))
    short_used = sum(int(row.get("n_short_used") or 0) for row in ev_rows)

    # ---------------------------------------------------------------- output
    lines = []
    P = lines.append
    P("| variant | role | n | lifted | disc (1-sided) | ev | pooled>=2 | comb (1-sided) | elem | struct | strand | tag recall | extra tags | src id | beyond (>=10bp) | beyond >=2 frag | T/L/U/A |")
    P("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    tot = {"TP": Counter(), "ARTEFACT": Counter()}
    for k in sorted(var, key=lambda x: (roles[x] != "TP", x)):
        v = var[k]
        tot[roles[k]].update(v)
        P(_row(k, roles[k], v))
    for role in ("TP", "ARTEFACT"):
        P(_row(f"**all {role}**", role, tot[role]))
    P("")
    P(f"TPRT call by role (combined, annotated): " + "; ".join(
        f"{role}: " + ", ".join(f"{c}={n}" for c, n in sorted(cnt.items())) for role, cnt in calls_by_role.items()))
    P("")
    P("| half | AUC TP vs simulated ARTEFACT | AUC TP vs ARTEFACT+unexplained | n TP | n ART | n unexplained |")
    P("|---|---|---|---|---|---|")
    for half in ("fit", "eval"):
        tp, art = scores[half]["TP"], scores[half]["ARTEFACT"]
        un = scores[half + "_unexpl"]["ARTEFACT"]
        P(f"| {half} | {_f(auc(tp, art))} | {_f(auc(tp, art + un))} | {len(tp)} | {len(art)} | {len(un)} |")
    P("")
    P(f"Within-sample duplicates (qname truth over {dup['junctions']} junctions, {dup['fragments']} evidence fragments): "
      f"{dup['true_dup_fragments']} true PCR-copy fragments; combine merged {dup['reported_duplicates']} "
      f"(coord {dup['dup_coord']}, seq {dup['dup_seq']}); n_independent == true molecules at {dup['exact']} junctions, "
      f"under-collapsed {dup['under_collapsed']} ({dup['missed_dups']} dups missed), over-collapsed "
      f"{dup['over_collapsed']} ({dup['false_merges']} false merges); {dup['unresolved_frags']} frags without qname.")
    P("")
    P(f"SHORT overhang fragments used: {short_used}; junctions reaching >= 2 independent fragments only thanks "
      f"to SHORT reads: {short_only}.")
    disc_all = set().union(*[{(s, c) for c in dcalls[s]} for s in samples])
    P(f"Discovery calls: {len(disc_all)} over {len(samples)} colonies, {len(disc_all - explained_disc)} unexplained; "
      f"combined insertions: {len(comb)}, unexplained (organic FP) {len(unexplained)}.")
    txt = "\n".join(lines)
    with open(a.report, "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    with open(os.path.join(O, "e2e_events.tsv"), "w") as fh:
        keys = ["id", "variant", "role", "element", "structure", "tags", "strand_ref", "polya_len", "c", "lo", "hi",
                "disc", "ev_name", "pooled", "comb_name"]
        fh.write("\t".join(keys + ["a_element", "a_structure", "a_tags", "a_call", "a_score", "a_points", "a_detail"]) + "\n")
        for r in truth:
            an = r["ann"] or {}
            fh.write("\t".join(str(r.get(k)) for k in keys) + "\t" + "\t".join(
                str(an.get(k, ".")) for k in ("element", "structure", "tags", "tprt_call", "tprt_score",
                                              "tprt_points", "rte_detail")) + "\n")


def _f(x):
    return "." if x is None else f"{x:.3f}"


def _pct(a, b):
    return f"{a}/{b}" if b else "."


def _row(k, role, v):
    calls = "/".join(str(v["call_" + c]) for c in ("TPRT", "LIKELY_TPRT", "UNCERTAIN", "ARTEFACT_LIKE"))
    return (f"| {k} | {role} | {v['n']} | {v['lifted']} | {v['disc']} ({v['disc1s']}) | {v['ev']} | {v['pooled']} | "
            f"{v['comb']} ({v['comb1s']}) | "
            f"{_pct(v['elem'], v['ann'])} | {_pct(v['struct'], v['ann'])} | {_pct(v['strand'], v['strand_n'])} | "
            f"{_pct(v['tag_hit'], v['tag_n'])} | {v['tag_extra']} | {_pct(v['src'], v['src_n'])} | "
            f"{_pct(v['beyond'], v['polya_n'])} | {_pct(v['beyond2'], v['polya_n'])} | {calls} |")


if __name__ == "__main__":
    main()
