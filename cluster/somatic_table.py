#!/usr/bin/env python3
"""Annotated Excel table of a patient's somatic insertion candidates (stdlib only, no openpyxl).

Candidates = loci the genotype2 joint step places on a tree branch below the root (clade or
private), plus the loci the legacy run's tree_fit calls private or phylo-consistent shared, plus
the known insertions. One row per locus with: the joint verdict and per-carrier reads (genotype2),
tree_fit's class on both genotypers, the annotate_v2 annotation (element class, TPRT call/score,
TSD, poly-A, insertion site), whether it is a known insertion, and a priority tier.

  python cluster/somatic_table.py --patient PD37590 --joint bias3.joint.tsv \
      --genotype-dir V2_refbias/genotypes --fit-v2 V2_refbias/fit/phylo_fit.tsv \
      --fit-legacy C_rust/eval/fit/phylo_fit.tsv --annotation C_rust/PD37590/PD37590.annotated.csv.gz \
      --known patients/colorectum/PD37590/known_insertions.tsv --out PD37590.somatic.xlsx

Sheets: README, somatic (the table, sorted by tier), known (the known insertions' rows),
refbias (the b estimates, when --refbias is given).
"""
import argparse
import collections
import csv
import gzip
import os
import re
import sys
import zipfile
from xml.sax.saxutils import escape

RTE_CLASSES = {"LINE1", "ALU", "SVA", "RTE_other"}
TPRT_CALLS = {"TPRT", "LIKELY_TPRT"}
TIERS = {
    "A": "known insertion",
    "B": "genotype2 clade/private, every carrier P >= 0.9, RTE class or TPRT/LIKELY_TPRT call",
    "C": "genotype2 clade/private, every carrier P >= 0.9",
    "D": "genotype2 clade/private with a carrier P < 0.9, or only the legacy tree_fit supports it",
}


# ------------------------------------------------------------------ input
def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        lines = [l for l in fh if not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def parse_locus(name):
    m = re.match(r"^([^:]+):(\d+)-(\d+)", name or "")
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def joint_class(r):
    b = r.get("best", "")
    if b in ("ROOT", "NOISE", "INDEP", ""):
        return b or "NA"
    try:
        n = int(r.get("n_carriers") or 0)
    except ValueError:
        n = 0
    return "clade" if n >= 2 else "private"


def fnum(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def match_known(locus, known, tol=30):
    p = parse_locus(locus)
    if not p:
        return None
    c, a, b = p
    lo, hi = min(a, b), max(a, b)
    for k in known:
        q = parse_locus(k["locus"])
        if q and q[0] == c and min(q[1], q[2]) - tol <= hi and lo <= max(q[1], q[2]) + tol:
            return k
    return None


# ------------------------------------------------------------------ xlsx (minimal SpreadsheetML)
def col_letter(i):
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def sheet_xml(rows, widths=None, freeze=True, autofilter=True):
    out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
           '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">']
    if freeze and rows:
        out.append('<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
                   'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>')
    if widths:
        out.append("<cols>" + "".join(f'<col min="{i + 1}" max="{i + 1}" width="{w}" customWidth="1"/>'
                                      for i, w in enumerate(widths)) + "</cols>")
    out.append("<sheetData>")
    for ri, row in enumerate(rows):
        cells = []
        for ci, v in enumerate(row):
            ref = f"{col_letter(ci)}{ri + 1}"
            style = ' s="1"' if ri == 0 and autofilter else ""
            if v is None or v == "":
                continue
            if isinstance(v, bool):
                v = "yes" if v else "no"
            if isinstance(v, (int, float)) and not (isinstance(v, float) and v != v):
                cells.append(f'<c r="{ref}"{style}><v>{v}</v></c>')
            else:
                t = escape(str(v))[:32000]
                cells.append(f'<c r="{ref}" t="inlineStr"{style}><is><t xml:space="preserve">{t}</t></is></c>')
        out.append(f'<row r="{ri + 1}">' + "".join(cells) + "</row>")
    out.append("</sheetData>")
    if autofilter and rows and rows[0]:
        out.append(f'<autoFilter ref="A1:{col_letter(len(rows[0]) - 1)}{max(1, len(rows))}"/>')
    out.append("</worksheet>")
    return "".join(out)


def write_xlsx(path, sheets):
    """sheets: list of (name, rows, widths, is_table)"""
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        n = len(sheets)
        z.writestr("[Content_Types].xml",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                   '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                   '<Default Extension="xml" ContentType="application/xml"/>'
                   '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
                   '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
                   + "".join(f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" '
                             'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                             for i in range(n)) + "</Types>")
        z.writestr("_rels/.rels",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
                   "</Relationships>")
        z.writestr("xl/workbook.xml",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                   'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
                   + "".join(f'<sheet name="{escape(s[0])}" sheetId="{i + 1}" r:id="rId{i + 1}"/>' for i, s in enumerate(sheets))
                   + "</sheets></workbook>")
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   + "".join(f'<Relationship Id="rId{i + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i + 1}.xml"/>'
                             for i in range(n))
                   + f'<Relationship Id="rId{n + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
                   "</Relationships>")
        z.writestr("xl/styles.xml",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                   '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
                   '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
                   '<fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill>'
                   '<fill><patternFill patternType="solid"><fgColor rgb="FFDDEBF7"/></patternFill></fill></fills>'
                   '<borders count="1"><border/></borders>'
                   '<cellStyleXfs count="1"><xf/></cellStyleXfs>'
                   '<cellXfs count="2"><xf xfId="0"/><xf xfId="0" fontId="1" fillId="2" applyFont="1" applyFill="1"/></cellXfs>'
                   '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
                   "</styleSheet>")
        for i, (name, rows, widths, is_table) in enumerate(sheets):
            z.writestr(f"xl/worksheets/sheet{i + 1}.xml", sheet_xml(rows, widths, freeze=is_table, autofilter=is_table))


# ------------------------------------------------------------------ combine outputs
def read_consensus(path, wanted):
    """combined.txt.gz (FASTQ, titles <locus>:L / <locus>:R) -> {(locus, 'L'|'R'): seq}"""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        while True:
            t = fh.readline()
            if not t:
                break
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            name = t.strip()[1:]
            loc, _, side = name.rpartition(":")
            if loc in wanted and side in ("L", "R"):
                out[(loc, side)] = seq
    return out


def read_evidence(path, wanted):
    """insertions.evidence.tsv.gz -> {(locus, 'L'|'R'): row} (side LEFT/RIGHT mapped to L/R)"""
    out = {}
    for r in read_tsv(path):
        loc = r.get("insertion_id", "")
        if loc in wanted:
            side = {"LEFT": "L", "RIGHT": "R"}.get(r.get("side", "").upper(), r.get("side", ""))
            out.setdefault((loc, side), r)
    return out


def read_reads(path, wanted):
    """insertions.reads.fa.gz (>locus|side|role|sample|frag|r12) -> {locus: [(side, role, sample, frag, r12, seq)]}"""
    out = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        head = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                head = line[1:].split("|")
            elif head is not None:
                if head[0] in wanted:
                    h = head + [""] * (6 - len(head))
                    out[head[0]].append((h[1], h[2], h[3], h[4], h[5], line))
                head = None
    return out


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--patient", required=True)
    ap.add_argument("--joint", required=True, help="genotype2 joint table (<P>.joint.tsv)")
    ap.add_argument("--genotype-dir", required=True, help="genotype2 per-colony files <colony>.txt.gz")
    ap.add_argument("--fit-v2", help="tree_fit phylo_fit.tsv on the genotype2 files")
    ap.add_argument("--fit-legacy", help="tree_fit phylo_fit.tsv on the legacy genotypes")
    ap.add_argument("--annotation", help="annotate_v2 table <P>.annotated.csv.gz")
    ap.add_argument("--known", help="patients/<organ>/<P>/known_insertions.tsv")
    ap.add_argument("--refbias", help="joint reference-bias table (<P>.joint.refbias.tsv)")
    ap.add_argument("--insertions-dir",
                    help="combine output dir: adds clip consensus (<P>.combined.txt.gz), per-side evidence "
                         "(<P>.insertions.evidence.tsv.gz) and a reads sheet with every clipped read and mate "
                         "(<P>.insertions.reads.fa.gz)")
    ap.add_argument("--min-p", type=float, default=0.9, help="carrier P(carrier) for tiers B/C (0.9 = annotate_v2's)")
    ap.add_argument("--out", required=True, help="output .xlsx")
    a = ap.parse_args()

    joint = {r["locus"]: r for r in read_tsv(a.joint)}
    if not joint:
        sys.exit(f"no rows in {a.joint}")
    colonies = [k[2:] for k in next(iter(joint.values())) if k.startswith("p_")]
    fit2 = {r["locus"]: r for r in read_tsv(a.fit_v2)}
    fitl = {r["locus"]: r for r in read_tsv(a.fit_legacy)}
    ann = {r["locus"]: r for r in read_tsv(a.annotation)}
    known = [k for k in read_tsv(a.known) if k.get("locus")]

    # candidate set
    cand = collections.OrderedDict()
    for loc, r in joint.items():
        if joint_class(r) in ("clade", "private"):
            cand[loc] = "genotype2"
    for loc, r in fitl.items():
        cls, lab = r.get("class", ""), r.get("label", "")
        if cls == "private" or (cls == "informative_shared" and lab == "phylo_consistent"):
            cand[loc] = cand.get(loc, "") and "both" or "legacy"
    known_hit = {}
    for loc in joint:
        k = match_known(loc, known)
        if k:
            known_hit[loc] = k
            cand.setdefault(loc, "known")

    # per-colony reads for the candidates
    gt = {}
    for c in colonies:
        p = os.path.join(a.genotype_dir, c + ".txt.gz")
        if os.path.exists(p):
            gt[c] = {g["locus"]: g for g in read_tsv(p) if g["locus"] in cand}

    cons, evid, reads_by = {}, {}, {}
    if a.insertions_dir:
        stem = os.path.join(a.insertions_dir, a.patient)
        cons = read_consensus(stem + ".combined.txt.gz", cand)
        evid = read_evidence(stem + ".insertions.evidence.tsv.gz", cand)
        reads_by = read_reads(stem + ".insertions.reads.fa.gz", cand)
        print(f"combine outputs: consensus for {len({k[0] for k in cons})} loci, evidence for "
              f"{len({k[0] for k in evid})}, reads for {len(reads_by)} ({sum(map(len, reads_by.values()))} reads)")

    header = ["tier", "locus", "chrom", "start", "end", "kind", "source", "known", "known_carriers",
              "joint_class", "joint_best", "n_carriers", "carriers", "min_P_carrier", "post_best", "log10_bf_tree",
              "carrier_reads (alt/ref/uninf)", "carrier_vaf_mean", "alt_reads_in_noncarriers", "n_noncarriers_with_alt",
              "treefit_v2_class", "treefit_v2_label", "treefit_v2_carriers",
              "treefit_legacy_class", "treefit_legacy_label", "treefit_legacy_carriers",
              "element_class", "element", "tprt_call", "tprt_score", "tsd_len", "tsd_seq", "polya_len",
              "left_polyA", "right_polyA", "en_motif", "structure", "element_identity", "nearest_active",
              "site_region", "site_gene", "site_strand", "conclusion",
              "L_n_reads", "L_n_fragments", "L_n_samples", "L_n_mates", "L_polya_len", "L_beyond_polya",
              "R_n_reads", "R_n_fragments", "R_n_samples", "R_n_mates", "R_polya_len", "R_beyond_polya",
              "reads_by_role", "L_junction (REF upper | clip lower)", "R_junction (REF upper | clip lower)",
              "L_clip_consensus", "R_clip_consensus"]
    rows = []
    for loc, src in cand.items():
        j = joint.get(loc, {})
        f2, fl, an = fit2.get(loc, {}), fitl.get(loc, {}), ann.get(loc, {})
        k = known_hit.get(loc)
        jc = joint_class(j) if j else "NA"
        carriers = [c for c in (j.get("carriers") or "").split(",") if c] if jc in ("clade", "private") else []
        ps = [fnum(j.get("p_" + c)) for c in carriers]
        min_p = min(ps) if ps else float("nan")
        reads, vafs = [], []
        for c in carriers:
            g = gt.get(c, {}).get(loc)
            if g:
                reads.append(f"{c}:{g['n_alt']}/{g['n_ref']}/{g['n_uninf']}")
                v = fnum(g.get("vaf"))
                if v == v:
                    vafs.append(v)
        nc_alt, nc_n = 0, 0
        for c in gt:
            if c in carriers:
                continue
            g = gt[c].get(loc)
            if g:
                n = int(fnum(g.get("n_alt"), 0))
                nc_alt += n
                nc_n += n > 0
        ecls = an.get("class", "")
        tcall = an.get("tprt_call", "")
        good_p = ps and min_p >= a.min_p
        if k:
            tier = "A"
        elif jc in ("clade", "private") and good_p and (ecls in RTE_CLASSES or tcall in TPRT_CALLS):
            tier = "B"
        elif jc in ("clade", "private") and good_p:
            tier = "C"
        else:
            tier = "D"
        p = parse_locus(loc) or ("", "", "")
        num = lambda x: (round(fnum(x), 4) if fnum(x) == fnum(x) else (x or ""))
        clean = lambda x: "" if x in (None, ".", "NA", "nan") else x
        rows.append([tier, loc, p[0], p[1], p[2], j.get("locus_kind") or f2.get("locus_kind", ""), src,
                     bool(k), k.get("carriers", "") if k else "",
                     jc, j.get("best", ""), len(carriers) if carriers else "", ",".join(carriers),
                     round(min_p, 4) if min_p == min_p else "", num(j.get("post_best")), num(j.get("log10_bf_tree")),
                     "; ".join(reads), round(sum(vafs) / len(vafs), 3) if vafs else "", nc_alt, nc_n,
                     f2.get("class", ""), f2.get("label", ""), f2.get("carriers", ""),
                     fl.get("class", ""), fl.get("label", ""), fl.get("carriers", ""),
                     clean(ecls) or ("(not annotated)" if not an else ""), clean(an.get("element", "")),
                     clean(tcall), num(clean(an.get("tprt_score", ""))), num(clean(an.get("tsd_len", ""))),
                     clean(an.get("tsd_seq", "")), num(clean(an.get("polya_len", ""))),
                     clean(an.get("left_polyA", "")), clean(an.get("right_polyA", "")),
                     clean(an.get("en_motif", "")), clean(an.get("structure", "")),
                     num(clean(an.get("element_identity", ""))), clean(an.get("nearest_active", "")),
                     clean(an.get("site_region", "")), clean(an.get("site_gene", "")), clean(an.get("site_strand", "")),
                     clean(an.get("conclusion", ""))] + side_cols(loc, evid, cons, reads_by))
    rank = {"A": 0, "B": 1, "C": 2, "D": 3}
    rows.sort(key=lambda r: (rank[r[0]], -(r[11] or 0), -(fnum(r[29], 0)), r[1]))
    widths = [5, 28, 7, 11, 11, 18, 10, 7, 30, 10, 16, 6, 40, 9, 9, 9, 50, 9, 9, 9,
              18, 16, 30, 18, 16, 30, 18, 18, 12, 8, 7, 14, 8, 6, 6, 9, 14, 9, 14, 14, 14, 6, 60,
              7, 7, 7, 7, 7, 9, 7, 7, 7, 7, 7, 9, 30, 60, 60, 50, 50]

    tier_n = collections.Counter(r[0] for r in rows)
    readme = [["PEAR-TREE somatic insertion candidates", a.patient],
              ["", ""],
              ["rows", len(rows)]] + [[f"tier {t}", f"{tier_n.get(t, 0)}  -  {d}"] for t, d in TIERS.items()] + [
              ["", ""],
              ["inputs", ""],
              ["genotype2 joint", os.path.abspath(a.joint)],
              ["genotype2 per-colony", os.path.abspath(a.genotype_dir)],
              ["tree_fit (genotype2)", os.path.abspath(a.fit_v2) if a.fit_v2 else "-"],
              ["tree_fit (legacy)", os.path.abspath(a.fit_legacy) if a.fit_legacy else "-"],
              ["annotation", os.path.abspath(a.annotation) if a.annotation else "-"],
              ["known insertions", os.path.abspath(a.known) if a.known else "-"],
              ["", ""],
              ["columns", ""],
              ["source", "genotype2 = joint step places it on a branch below the root; legacy = legacy tree_fit "
                         "calls it private / phylo-consistent shared; both; known = only via the known list"],
              ["joint_class", "clade (>= 2 carriers below one branch), private (one tip), ROOT, INDEP, NOISE"],
              ["min_P_carrier", "lowest per-colony P(carrier) among the joint carriers; annotate_v2 counts a carrier at >= 0.9. "
                                "Judge private calls on this, not on log10_bf_tree (a private event's BF saturates near 2)"],
              ["carrier_reads", "per carrier colony: alt / ref / uninformative reads from the genotype2 realignment"],
              ["alt_reads_in_noncarriers", "alt reads summed over all other colonies (background / missed carriers)"],
              ["treefit_*", "tools/phylo/tree_fit.py read-vote model on each genotyper: class (germline, informative_shared, "
                            "private, noise, uninformative_depth) and label (phylo_consistent / ambiguous / phylo_violating)"],
              ["element_class ... conclusion", "annotate_v2 on the legacy run (only loci with a legacy het/hom call are annotated; "
                                               "'(not annotated)' otherwise)"],
              ["tprt_call / tprt_score", "TPRT hallmark call (TPRT, LIKELY_TPRT, UNCERTAIN, ...) and score"],
              ["L_/R_ n_reads ... beyond_polya", "combine evidence per insertion end (L = left junction, R = right): reads, "
                                                 "distinct fragments, colonies, mates, median poly-A length, sequence beyond the poly-A"],
              ["L_/R_junction", "combine junction consensus: reference flank in UPPER case, clipped (inserted) sequence in lower case"],
              ["L_/R_clip_consensus", "the clip consensus written to <P>.combined.txt.gz"],
              ["reads sheet", "every read combine kept for the locus (<P>.insertions.reads.fa.gz): side, role (CLIP = split read at "
                              "the junction, POLYA, DISC = discordant pair, SPAN, SHORT, MATE = the mate of an evidence read), "
                              "colony, whether that colony is a joint carrier, fragment id, read 1/2, sequence in allele-forward "
                              "orientation. Filter by locus"]]
    sheets = [("README", readme, [26, 120], False), ("somatic", [header] + rows, widths, True)]
    krows = [r for r in rows if r[0] == "A"]
    sheets.append(("known", [header] + krows, widths, True))
    if reads_by:
        order = {r[1]: (i, r[0], r[12]) for i, r in enumerate(rows)}
        rr = [["tier", "locus", "side", "role", "colony", "joint_carrier", "fragment", "read", "length", "sequence"]]
        for loc in sorted(reads_by, key=lambda x: order.get(x, (10 ** 9, "", ""))[0]):
            _, tier, car = order.get(loc, (0, "", ""))
            cs = set(car.split(",")) if car else set()
            for side, role, sample, frag, r12, seq in sorted(reads_by[loc], key=lambda r: (r[0], r[1] != "CLIP", r[1], r[2])):
                rr.append([tier, loc, side, role, sample, "yes" if sample in cs else "no", frag, r12, len(seq), seq])
        sheets.append(("reads", rr, [5, 28, 7, 7, 18, 8, 18, 5, 7, 160], True))
    rb = read_tsv(a.refbias)
    if rb:
        h = list(rb[0].keys())
        sheets.append(("refbias", [h] + [[num_or(r[c]) for c in h] for r in rb], [10, 22] + [12] * (len(h) - 2), True))
    write_xlsx(a.out, sheets)
    print(f"wrote {a.out}: {len(rows)} loci ({', '.join(f'{t} {tier_n.get(t, 0)}' for t in TIERS)})")


def side_cols(loc, evid, cons, reads_by):
    out = []
    for side in ("L", "R"):
        e = evid.get((loc, side), {})
        out += [num_or(e.get("n_reads", "")), num_or(e.get("n_fragments", "")), num_or(e.get("n_samples", "")),
                num_or(e.get("n_mates", "")), num_or(e.get("polya_len_median", "")), e.get("beyond_polya", "")]
    roles = collections.Counter(f"{r[0]}:{r[1]}" for r in reads_by.get(loc, []))
    out.append(", ".join(f"{k} {v}" for k, v in sorted(roles.items())))
    out += [evid.get((loc, "L"), {}).get("clip_consensus", ""), evid.get((loc, "R"), {}).get("clip_consensus", ""),
            cons.get((loc, "L"), ""), cons.get((loc, "R"), "")]
    return out


def num_or(x):
    try:
        f = float(x)
        return int(f) if f.is_integer() and "." not in str(x) else f
    except (TypeError, ValueError):
        return x


if __name__ == "__main__":
    main()
