#!/usr/bin/env python3
"""Annotated Excel table of a patient's somatic insertion candidates (stdlib only, no openpyxl).

Candidates = loci the genotype2 joint step places on a tree branch below the root (clade or
private), plus the loci the legacy run's tree_fit calls private or phylo-consistent shared, plus
the known insertions. One row per locus with: the joint verdict and per-carrier reads (genotype2),
tree_fit's class on both genotypers, the annotate_v2 annotation (element class, TPRT call/score,
TSD, poly-A, insertion site), whether it is a known insertion, and a priority tier.

Artefact checks cap a call at tier D and say why in tier_note:
  * phylo_violating   tree_fit (genotype2) labels the locus phylo_violating
  * scattered alt     a private call with alt reads in more than --max-noncarrier-colonies other colonies
  * local origin      (--genome) every clip is reference sequence within --local-window bp and the
                      mates lie in the flanks: nothing was inserted (template switch between nearby
                      repeat copies, small del/dup) -- PD37580 lo0077's six L1 "insertions"

  python cluster/somatic_table.py --patient PD37590 --joint bias3.joint.tsv \
      --genotype-dir V2_refbias/genotypes --fit-v2 V2_refbias/fit/phylo_fit.tsv \
      --fit-legacy C_rust/eval/fit/phylo_fit.tsv --annotation C_rust/PD37590/PD37590.annotated.csv.gz \
      --known patients/colorectum/PD37590/known_insertions.tsv --out PD37590.somatic.xlsx

Hard rules (Jeremy, 2026-10-08) EXCLUDE a candidate (moved to the `excluded` sheet with the reason):
  * site gap          |TSD / target-site deletion| > 120 bp (MAX_SITE_GAP)
  * junction support  no colony with >= 2 INDEPENDENT fragments on BOTH ends: CLIP / POLYA reads count;
                      next to >= 1 such clip, a SHORT overhang or a DISC pair whose inside mate AGREES with
                      that clip's inserted sequence (>= 25 bp at >= 90 %) adds a fragment; a discordant
                      pair alone, or one that does not agree on the insert, never counts; templates whose
                      read AND mate agree within 5 bp (allele-forward shift, R1/R2 ignored) are one molecule
  Known insertions (tier A) are kept and only flagged in tier_note.
Sheets: README, somatic (the table, sorted by tier), excluded (the rule failures), known (the known
insertions' rows), refbias (the b estimates, when --refbias is given).
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
    "D": "genotype2 clade/private with a carrier P < 0.9, or only the legacy tree_fit supports it, "
         "or capped by an artefact check (tier_note)",
}
LOCAL_K_CLIP, LOCAL_K_MATE = 8, 15  # clip seeds short enough that one mismatch in 15 bp still seeds
LOCAL_MIN_CLIP = 15          # clip bp left after removing homopolymer runs >= 5
LOCAL_MATE_VOTES = 0.5       # fraction of a mate's 15-mers on one diagonal = it maps there (a paralog
                             # copy at ~90% identity keeps ~0.2-0.4; sequencing errors still leave >= 0.5)
LOCAL_MATE_FRAC = 0.7        # share of the locus' mates that must map inside the window


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


# ------------------------------------------------------------------ hard rules
MAX_SITE_GAP = 120            # bp, |TSD / target-site deletion| of a two-sided locus
MIN_JUNCTION_FRAGMENTS = 2    # independent junction-crossing fragments per end, in one colony
DUP_SHIFT = 5                 # bp: read and mate both within this shift = one molecule (= discovery/combine)
JUNCTION_ROLES = ("CLIP", "POLYA")
AGREE_MIN_BP = 25             # a DISC mate must share >= this many bp of a clip's insert ...
AGREE_MIN_ID = 0.9            # ... at >= this identity to count as a second fragment
AGREE_LC_WINDOW = 12          # low complexity: a base in any 12-bp window with <= 4 distinct 3-mers
AGREE_LC_MAX_TRIMERS = 4      # (poly-A/T, di-/tri-nucleotide repeats) never counts toward AGREE_MIN_BP


def site_gap(loc):
    p = parse_locus(loc)
    return None if p is None else p[2] - p[1]


def shift_match(a, b, tol=DUP_SHIFT, min_id=0.8):
    """True when b starts within tol bases of a (allele-forward, same orientation) -- the sequence-only
    stand-in for "same outer coordinate". Identity (>= min_id over the overlap) only LOCATES the offset:
    clipped tails of PCR duplicates carry many sequencing errors (low-quality poly-A/T), so it is lax."""
    if not a or not b:
        return False
    for d in range(-tol, tol + 1):
        x, y = (a[d:], b) if d >= 0 else (a, b[-d:])
        n = min(len(x), len(y))
        if n < 30:
            continue
        if sum(1 for i in range(n) if x[i] == y[i]) >= min_id * n:
            return True
    return False


def _molecules(rows):
    """single-linkage PCR-duplicate collapse of rows (kind, frag, read_seq, mate_seq): one component per
    molecule (read and mate both shift-match; a mate missing on either side: the read decides)"""
    n = len(rows)
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(n):
        for j in range(i + 1, n):
            mi, mj = rows[i][3], rows[j][3]
            # a missing mate is unknown, not evidence of a second molecule: the read alone decides
            if shift_match(rows[i][2], rows[j][2]) and (not (mi and mj) or shift_match(mi, mj)):
                parent[root(j)] = root(i)
    comp = collections.defaultdict(list)
    for i in range(n):
        comp[root(i)].append(rows[i])
    return list(comp.values())


def insert_part(read, side, record):
    """the inserted (clipped) part of an allele-forward junction read, located with the combined record's
    reference flank (`:L` = clip + FLANK, `:R` = FLANK + clip); falls back to the record's clip"""
    clip = "".join(c for c in record if c.islower()).upper()
    flank = "".join(c for c in record if c.isupper())
    if len(flank) >= 15:
        if side == "L":
            k = read.find(flank[:15])
            if k > 0:
                return read[:k]
        else:
            k = read.find(flank[-15:])
            if k >= 0:
                return read[k + 15:]
    return clip


def low_complexity_mask(s):
    """True per base inside low-complexity sequence (= discovery model.rs `low_complexity_mask`)"""
    w = min(AGREE_LC_WINDOW, len(s))
    if w < 3:
        return [True] * len(s)
    mask = [False] * len(s)
    for st in range(len(s) - w + 1):
        win = s[st:st + w]
        if len({win[i:i + 3] for i in range(w - 2)}) <= AGREE_LC_MAX_TRIMERS:
            mask[st:st + w] = [True] * w
    return mask


def agrees(mate, insert, k=12):
    """mate and insert share >= AGREE_MIN_BP on one diagonal at >= AGREE_MIN_ID (same orientation:
    both allele-forward), i.e. they agree on what was inserted; >= AGREE_MIN_BP of the matching
    bases must lie outside low-complexity sequence (a shared poly-A is no agreement)"""
    if not mate or len(insert) < AGREE_MIN_BP:
        return False
    lc = low_complexity_mask(insert)
    if lc.count(False) < AGREE_MIN_BP:
        return False
    seeds = collections.defaultdict(list)
    for i in range(len(insert) - k + 1):
        seeds[insert[i:i + k]].append(i)
    tried = set()
    for j in range(len(mate) - k + 1):
        for i in seeds.get(mate[j:j + k], ()):
            d = j - i
            if d in tried:
                continue
            tried.add(d)
            lo, hi = max(0, -d), min(len(insert), len(mate) - d)
            n = hi - lo
            if n < AGREE_MIN_BP:
                continue
            hit = [x for x in range(lo, hi) if insert[x] == mate[x + d]]
            if len(hit) >= AGREE_MIN_ID * n and sum(not lc[x] for x in hit) >= AGREE_MIN_BP:
                return True
    return False


def end_support(frags, side, record):
    """frags: {frag: (role, frag, read, mate)} of one colony at one end -> independent fragments:
    CLIP/POLYA molecules, plus (only next to >= 1 of them) SHORT-only molecules and DISC molecules whose
    inside mate agrees with a clip's insert"""
    rows = list(frags.values())
    clips = [r for r in rows if r[0] in JUNCTION_ROLES]
    if not clips:
        return 0
    inserts = [insert_part(r[2], side, record) for r in clips]
    usable = [r for r in rows if r[0] in JUNCTION_ROLES or r[0] == "SHORT"
              or (r[0] == "DISC" and any(agrees(r[3], ins) for ins in inserts))]
    return len(_molecules(usable))


def junction_support(reads, cons=None, loc=None):
    """reads: [(side, role, sample, frag, r12, seq)] of one locus -> {sample: (L, R)} independent
    fragments per end (reads.fa LEFT/RIGHT -> L/R); cons = combined records {(locus, 'L'|'R'): seq}"""
    cons = cons or {}
    mates = {(r[0], r[2], r[3]): r[5] for r in reads if r[1] == "MATE"}
    per = collections.defaultdict(lambda: collections.defaultdict(dict))
    for side, role, sample, frag, r12, seq in reads:
        if role in JUNCTION_ROLES or role in ("SHORT", "DISC"):
            s = {"LEFT": "L", "RIGHT": "R"}.get(side, side)
            per[sample][s].setdefault(frag, (role, frag, seq, mates.get((side, sample, frag), "")))
    return {smp: tuple(end_support(d.get(s, {}), s, cons.get((loc, s), "")) for s in ("L", "R"))
            for smp, d in per.items()}


def hard_rules(loc, reads_by, check_fragments, cons=None):
    """-> (violations, support text); support = the best colony's independent fragments per end"""
    bad = []
    g = site_gap(loc)
    if g is not None and abs(g) > MAX_SITE_GAP:
        bad.append(f"site gap {g} bp (|gap| > {MAX_SITE_GAP})")
    txt = ""
    if check_fragments:
        sup = junction_support(reads_by.get(loc, []), cons, loc)
        best = max(sup.items(), key=lambda kv: (min(kv[1]), sum(kv[1])), default=None)
        if best:
            txt = f"{best[0]}: L{best[1][0]}/R{best[1][1]}"
        if not best or min(best[1]) < MIN_JUNCTION_FRAGMENTS:
            bad.append(f"< {MIN_JUNCTION_FRAGMENTS} independent junction fragments on an end in every colony"
                       + (f" (best {txt})" if txt else " (no junction reads)"))
    return bad, txt


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


# ------------------------------------------------------------------ local-origin check
_RC = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(_RC)[::-1]


def kmer_index(ref, k):
    idx = collections.defaultdict(list)
    for i in range(len(ref) - k + 1):
        idx[ref[i:i + k]].append(i)
    return idx


def place(q, idx, k):
    """best ungapped placement of q (either strand) in the indexed window:
    (votes / n_kmers, strand, ref offset of q[0]) or None"""
    best = None
    for strand, s in (("+", q), ("-", revcomp(q))):
        n = len(s) - k + 1
        if n <= 0:
            continue
        diag = collections.Counter()
        for i in range(n):
            for p in idx.get(s[i:i + k], ()):
                diag[p - i] += 1
        if diag:
            d, v = diag.most_common(1)[0]
            if best is None or v / n > best[0]:
                best = (v / n, strand, d)
    return best


def usable_clip(clip):
    """the clip with homopolymer runs >= 5 removed (poly-A tails, slippage) -- None when too short
    or low-complexity to place on its own"""
    core = re.sub(r"(A{5,}|C{5,}|G{5,}|T{5,})", "", clip)
    if len(core) < LOCAL_MIN_CLIP:
        return None
    if len({core[i:i + 4] for i in range(len(core) - 3)}) < 0.5 * (len(core) - 3):
        return None
    return clip


def local_origin(loc, evid, cons, reads_by, genome, window):
    """'' unless the clips are (near-)exact reference within +-window of the locus (<= 1 mismatch per
    20 bp, gapless; homopolymer / low-complexity clips are skipped, an unplaceable clip < 20 bp too
    once the other end placed), no end carries a poly-A >= 10, AND >= LOCAL_MATE_FRAC of the mates
    map inside the window; else a description. Mates are not held to
    a side: around a small del/dup they fall on either."""
    p = parse_locus(loc)
    if not p or genome is None:
        return ""
    c, a, b = p
    lo, hi = min(a, b), max(a, b)
    start = max(1, lo - window)
    ref = genome.fetch(c, start - 1, hi + window)          # 1-based position = index + start
    if not ref:
        return ""
    clips = []
    for side in ("L", "R"):
        # combine's consensus (combined.txt.gz) carries the clip in lower case; the evidence table's
        # clip_consensus sometimes has none (PD37580 13:23156430 LEFT)
        cl = ""
        for src in (cons.get((loc, side), ""), evid.get((loc, side), {}).get("clip_consensus", "")):
            cl = "".join(re.findall("[a-z]+", src)).upper()
            if cl:
                break
        if usable_clip(cl):
            clips.append((side, cl))
    if not clips:
        return ""
    # a TPRT insertion keeps its poly-A: never call that local (a young Alu clip can match a
    # reference Alu next door)
    if any(fnum(evid.get((loc, sd), {}).get("polya_len_median"), 0) >= 10 for sd in ("L", "R")):
        return ""
    idx = kmer_index(ref, LOCAL_K_CLIP)
    where, missed = [], []
    for side, cl in clips:
        hit = place(cl, idx, LOCAL_K_CLIP)
        mm = None
        if hit:
            _, strand, d = hit
            s = cl if strand == "+" else revcomp(cl)
            if 0 <= d and d + len(s) <= len(ref):
                mm = sum(x != y for x, y in zip(s, ref[d:d + len(s)]))
        if mm is None or mm > max(1, len(cl) // 20):
            missed.append(cl)
            continue
        pos = d + start
        dist = 0 if lo <= pos <= hi else min(abs(pos - lo), abs(pos - hi), abs(pos + len(s) - 1 - lo), abs(pos + len(s) - 1 - hi))
        where.append(f"{side} clip {len(cl)}bp = ref {c}:{pos}({strand}, {mm} mm, {dist} bp away)")
    # a clip < 20 bp that does not place (one mismatch can leave it without a seed) is no evidence
    # either way once the other end placed; a longer one is foreign sequence
    if not where or any(len(cl) >= 20 for cl in missed):
        return ""
    mates = [(r[0], r[5]) for r in reads_by.get(loc, []) if r[1] == "MATE"]
    if len(mates) < 2:
        return ""
    midx = kmer_index(ref, LOCAL_K_MATE)
    ok = 0
    for _, seq in mates:
        hit = place(seq.upper(), midx, LOCAL_K_MATE)
        if not hit or hit[0] < LOCAL_MATE_VOTES:
            continue
        ok += 1
    if ok < LOCAL_MATE_FRAC * len(mates):
        return ""
    return "local origin: " + "; ".join(where) + f"; {ok}/{len(mates)} mates within {window} bp"


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
    ap.add_argument("--rte-library", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "resources", "rte_library"),
                    help="RTE library dir: active.tsv describes the nearest active element (default: the repo's resources/rte_library)")
    ap.add_argument("--min-p", type=float, default=0.9, help="carrier P(carrier) for tiers B/C (0.9 = annotate_v2's)")
    ap.add_argument("--max-noncarrier-colonies", type=int, default=2,
                    help="a private call with alt reads in more other colonies is capped at tier D (default 2)")
    ap.add_argument("--genome", help="reference of the discovery genome (.2bit or FASTA): enables the local-origin check")
    ap.add_argument("--local-window", type=int, default=5000, help="local-origin check: bp either side of the locus")
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
    active = {r["id"]: r for r in read_tsv(os.path.join(a.rte_library, "active.tsv")) if r.get("id")}
    td_sources = {r["id"]: r for r in read_tsv(os.path.join(a.rte_library, "transduction_sources.tsv")) if r.get("id")}
    genome = None
    if a.genome:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
        from tools.rte.genome import open_genome
        genome = open_genome(a.genome)

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

    header = ["tier", "tier_note", "locus", "chrom", "start", "end", "kind", "source", "known", "known_carriers",
              "joint_class", "joint_best", "n_carriers", "carriers", "min_P_carrier", "post_best", "log10_bf_tree",
              "carrier_reads (alt/ref/uninf)", "carrier_vaf_mean", "alt_reads_in_noncarriers", "n_noncarriers_with_alt",
              "treefit_v2_class", "treefit_v2_label", "treefit_v2_carriers",
              "treefit_legacy_class", "treefit_legacy_label", "treefit_legacy_carriers",
              "element_class", "element", "tprt_call", "tprt_score", "tsd_len", "tsd_seq", "polya_len",
              "left_polyA", "right_polyA", "en_motif", "structure", "element_identity", "nearest_active",
              "active_subfamily", "active_ta_status", "active_tier", "active_hotness", "active_n_daughters",
              "active_locus_hg38", "td3p_source", "td3p_source_band", "td3p_source_subfamily", "td3p_source_hotness",
              "td3p_source_n_daughters", "td3p_source_hg38", "td3p_end_in_flank", "covered_5p", "covered_3p", "rte_tags", "rte_detail",
              "site_region", "site_gene", "site_strand", "conclusion",
              "L_n_reads", "L_n_fragments", "L_n_samples", "L_n_mates", "L_polya_len", "L_beyond_polya",
              "R_n_reads", "R_n_fragments", "R_n_samples", "R_n_mates", "R_polya_len", "R_beyond_polya",
              "reads_by_role", "L_junction (REF upper | clip lower)", "R_junction (REF upper | clip lower)",
              "L_clip_consensus", "R_clip_consensus", "junction_fragments (best colony L/R)"]
    rows, excluded = [], []
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
        notes = []
        if f2.get("label") == "phylo_violating":
            notes.append("phylo_violating (tree_fit genotype2)")
        if jc == "private" and nc_n > a.max_noncarrier_colonies:
            notes.append(f"alt reads in {nc_n} non-carrier colonies")
        tprt_like = tcall in TPRT_CALLS or fnum(an.get("polya_len"), 0) >= 10
        lo_note = local_origin(loc, evid, cons, reads_by, genome, a.local_window) if not (k or tprt_like) else ""
        if lo_note:
            notes.append(lo_note)
        if notes and tier in ("B", "C"):
            notes.insert(0, f"tier {tier} -> D")
            tier = "D"
        violations, support = hard_rules(loc, reads_by, bool(a.insertions_dir), cons)
        if violations and k:
            notes.append("KNOWN but fails: " + "; ".join(violations))
        p = parse_locus(loc) or ("", "", "")
        num = lambda x: (round(fnum(x), 4) if fnum(x) == fnum(x) else (x or ""))
        clean = lambda x: "" if x in (None, ".", "NA", "nan") else x
        rows.append([tier, "; ".join(notes), loc, p[0], p[1], p[2], j.get("locus_kind") or f2.get("locus_kind", ""), src,
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
                     *active_cols(active.get(an.get("nearest_active", ""), {})),
                     *td_cols(an, td_sources),
                     num(clean(an.get("covered_5p", ""))), num(clean(an.get("covered_3p", ""))),
                     clean(an.get("tags", "")), clean(an.get("rte_detail", "")),
                     clean(an.get("site_region", "")), clean(an.get("site_gene", "")), clean(an.get("site_strand", "")),
                     clean(an.get("conclusion", ""))] + side_cols(loc, evid, cons, reads_by) + [support])
        if violations and not k:
            excluded.append(["; ".join(violations)] + rows.pop())
    rank = {"A": 0, "B": 1, "C": 2, "D": 3}
    H = {h: i for i, h in enumerate(header)}
    rows.sort(key=lambda r: (rank[r[0]], bool(r[1]), -(r[H["n_carriers"]] or 0), -(fnum(r[H["tprt_score"]], 0)), r[H["locus"]]))
    excluded.sort(key=lambda r: (rank[r[1]], r[1 + H["locus"]]))
    widths = [5, 40, 28, 7, 11, 11, 18, 10, 7, 30, 10, 16, 6, 40, 9, 9, 9, 50, 9, 9, 9,
              18, 16, 30, 18, 16, 30, 18, 18, 12, 8, 7, 14, 8, 6, 6, 9, 14, 9, 14, 9, 7, 16, 12, 8, 30, 20, 18, 18, 10, 8, 30, 9, 8, 8, 20, 60, 14, 14, 6, 60,
              7, 7, 7, 7, 7, 9, 7, 7, 7, 7, 7, 9, 30, 60, 60, 50, 50, 22]

    tier_n = collections.Counter(r[0] for r in rows)
    readme = [["PEAR-TREE somatic insertion candidates", a.patient],
              ["", ""],
              ["rows", len(rows)]] + [[f"tier {t}", f"{tier_n.get(t, 0)}  -  {d}"] for t, d in TIERS.items()] + [
              ["excluded", f"{len(excluded)}  -  hard rules: |TSD/deletion| <= {MAX_SITE_GAP} bp; >= {MIN_JUNCTION_FRAGMENTS} "
                           "independent fragments (CLIP/POLYA; beside a clip also SHORT, and DISC pairs whose inside mate "
                           f"agrees with the clip's insert over >= {AGREE_MIN_BP} bp at >= {int(AGREE_MIN_ID * 100)} % (poly-A / simple repeats never count); "
                           f"read+mate within {DUP_SHIFT} bp = one PCR molecule) on BOTH ends in one colony"
                           + ("" if a.insertions_dir else "  [fragment rule NOT checked: no --insertions-dir]")],
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
              ["tier_note", "why a call was capped at tier D: phylo_violating (tree_fit genotype2), alt reads in more than "
                            f"{a.max_noncarrier_colonies} non-carrier colonies (private calls), or local origin -- every clip is "
                            f"reference within {a.local_window} bp and the mates lie in the flanks, so nothing was inserted"
                            + ("" if genome else " (local-origin check OFF: no --genome)")],
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
              ["nearest_active / element_identity", "the active L1 (resources/rte_library/active.tsv) whose sequence best matches the "
                                                    "INSERTED sequence that was assembled, and the identity over that covered part only -- "
                                                    "short covered parts (see covered_5p/3p) cannot tell young L1s apart"],
              ["active_tier / hotness / n_daughters", "hot_source = published source element (hotness strong/hot/active from its "
                                                      "daughter count); L1HS_Ta_intact / L1HS_preTa_intact = intact young L1HS without "
                                                      "(candidate / none_reported) or with reported activity"],
              ["td3p_source ...", "3' transduction (TD3P / ORPHAN_TD): the source L1 (resources/rte_library/transduction_sources.tsv), "
                                  "its band, subfamily, hotness and published daughters, hg38 position, and how far into its 3' flank "
                                  "the transduced sequence reaches (td3p_end_in_flank, bp)"],
              ["covered_5p / covered_3p", "element consensus coordinates the assembled insert covers (5' truncation point / 3' end)"],
              ["L_/R_ n_reads ... beyond_polya", "combine evidence per insertion end (L = left junction, R = right): reads, "
                                                 "distinct fragments, colonies, mates, median poly-A length, sequence beyond the poly-A"],
              ["L_/R_junction", "combine junction consensus: reference flank in UPPER case, clipped (inserted) sequence in lower case"],
              ["L_/R_clip_consensus", "the clip consensus written to <P>.combined.txt.gz"],
              ["reads sheet", "every read combine kept for the locus (<P>.insertions.reads.fa.gz): side, role (CLIP = split read at "
                              "the junction, POLYA, DISC = discordant pair, SPAN, SHORT, MATE = the mate of an evidence read), "
                              "colony, whether that colony is a joint carrier, fragment id, read 1/2, sequence in allele-forward "
                              "orientation. Filter by locus"]]
    sheets = [("README", readme, [26, 120], False), ("somatic", [header] + rows, widths, True),
              ("excluded", [["excluded_by"] + header] + excluded, [60] + widths, True)]
    krows = [r for r in rows if r[0] == "A"]
    sheets.append(("known", [header] + krows, widths, True))
    if reads_by:
        order = {r[H["locus"]]: (i, r[0], r[H["carriers"]]) for i, r in enumerate(rows)}
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
    print(f"wrote {a.out}: {len(rows)} loci ({', '.join(f'{t} {tier_n.get(t, 0)}' for t in TIERS)}), "
          f"{len(excluded)} excluded by the hard rules")


def active_cols(r):
    """active.tsv metadata of the nearest active element: subfamily, Ta status, tier (hot_source = a published
    source element; L1HS_Ta/preTa_intact = intact young L1HS), hotness, daughters, hg38 position"""
    if not r:
        return ["", "", "", "", "", ""]
    loc = f"{r.get('hg38_chrom', '')}:{r.get('hg38_start', '')}-{r.get('hg38_end', '')}({r.get('strand', '')})"
    return [r.get("subfamily", ""), r.get("ta_status", ""), r.get("tier", ""), r.get("hotness", ""),
            num_or(r.get("n_daughters", "")), loc]


def td_cols(an, td_sources):
    """3' transduction source (tag TD3P_SOURCE=<id>) described from transduction_sources.tsv, and how far
    into its 3' flank the transduced sequence reaches (rte_detail td_end)"""
    tags = (an.get("tags") or "").split(",")
    sid = next((t.split("=", 1)[1] for t in tags if t.startswith("TD3P_SOURCE=")), "")
    m = re.search(r"(?:^|;)td_end=([^;]+)", an.get("rte_detail") or "")
    td_end = num_or(m.group(1)) if m else ""
    r = td_sources.get(sid, {})
    if not sid:
        return ["", "", "", "", "", "", td_end]
    loc = f"{r.get('hg38_chrom', '')}:{r.get('hg38_start', '')}-{r.get('hg38_end', '')}({r.get('strand', '')})" if r else ""
    return [sid, r.get("band_published") or r.get("band", ""), r.get("subfamily", ""), r.get("hotness", ""),
            num_or(r.get("n_daughters", "")), loc, td_end]


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
