#!/usr/bin/env python3
"""Build the PEAR-TREE RTE reference library (resources/rte_library/).

Produces (see plans/tprt_hallmarks/SPEC.md, "Reference libraries"):
  l1_intact.fa/.tsv, alu_y_intact.fa/.tsv, sva_intact.fa/.tsv, consensus.fa,
  consensus_landmarks.tsv, consensus_crosscheck.tsv, dfam_young.fa, active.tsv,
  transduction_sources.tsv, flanks_3p.fa.gz (+.fai/.gzi), flanks_5p_sva.fa.gz (+.fai/.gzi),
  transduction_stats.tsv, manifest.tsv

All inputs are explicit paths (download recipe: tools/rte_library/fetch_inputs.sh). Large
intermediates (rmsk caches, MSAs) go to --work. Deterministic for fixed inputs.

Example (Jeremy's scratchpad layout):
  python tools/rte_library/build.py --inputs $SP --work $SP/work --out resources/rte_library
where $SP contains genomes/, libs/, supp/, ncbi/ as written by fetch_inputs.sh.
"""
import argparse
import collections
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from common import (TwoBit, Lifter, revcomp, read_fasta, write_fasta, identity, column_map,
                    load_cytobands, band_of, ungapped_identity, cons_identity)
import rmsk as rm

PRIMARY = set(["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"])

# ----------------------------------------------------------------------------- thresholds
# (rationale in resources/rte_library/README.md)
ALU_MIN_CONS_COV = 280        # bp of the ~311 bp AluY consensus covered by the rmsk alignment
ALU_MAX_MILLIDIV = 30         # per mille divergence from the subfamily consensus (rmsk milliDiv)
ALU_MAX_5P_MISSING = 15       # consensus bases allowed missing at the 5' end (A box must survive)
ALU_GENOMIC_LEN = (270, 360)  # genomic span (rmsk excludes the poly-A tail)
ALU_CAP_PER_SUBFAMILY = 150   # youngest N per subfamily kept in alu_y_intact.fa
SVA_MIN_CONS_COV = 1300       # bp of the ~1,380 bp (VNTR-collapsed) Dfam SVA consensus covered
SVA_MAX_5P_MISSING = 60        # Dfam SVA 5' end = (CCCTCT)n hexamer, often annotated as simple repeat
SVA_CAP_PER_SUBFAMILY = 60
L1_FL_MIN_SPAN = 5900         # "full-length young reference L1" for the transduction seed list
L1_YOUNG_FL = ("L1HS",)       # seeded with 3' flanks; L1PA2/3 FL only if intact ORFs or published
FLANK_3P_L1 = 15000           # bp downstream of an L1 source (Tubio 2014: up to 12 kb observed)
FLANK_3P_SVA = 5000
FLANK_5P_SVA = 5000
MSA_N_L1 = 40
MSA_N_ALU = 120
MSA_N_SVA = 40
# Ta/pre-Ta diagnostic sites: (31-bp L1.3 context, index of the diagnostic base). L1.3 5931
# ('ACA'/'ACG') and 5712 (A/G); Ta carries A at both.
TA_SITES = (("AATGCTAGATGACACATTAGTGGGTGCAGCG", 15), ("AAATGATGAGTTCATATCCTTTGTAGGGACA", 15))


def log(*a):
    print("[rte_library]", *a, file=sys.stderr, flush=True)


# ============================================================================ inputs
def default_paths(sp):
    g = os.path.join(sp, "genomes")
    return dict(
        hg38_2bit=os.path.join(g, "hg38.2bit"), hs1_2bit=os.path.join(g, "hs1.2bit"),
        hg38_rmsk=os.path.join(g, "hg38.rmsk.txt.gz"),
        hs1_rmsk=os.path.join(g, "hs1.repeatMasker.out.gz"),
        chain_hg38_hs1=os.path.join(g, "hg38ToHs1.over.chain.gz"),
        chain_hg19_hg38=os.path.join(g, "hg19ToHg38.over.chain.gz"),
        cyto_hg38=os.path.join(g, "hg38.cytoBand.txt.gz"),
        l1base_fa=os.path.join(sp, "libs", "l1base", "hsflil1_8438.fa"),
        dfam=os.path.join(sp, "libs", "dfam_consensus.fa"),
        l13_gb=os.path.join(sp, "ncbi", "L19088.gb"),
        l12_gb=os.path.join(sp, "ncbi", "M80343.gb"),
        l1rp_gb=os.path.join(sp, "ncbi", "AF148856.gb"),
        nam_st4=os.path.join(sp, "supp", "nam_MOESM7.xlsx"),
        nam_st2=os.path.join(sp, "supp", "nam_MOESM5.xlsx"),
        rm_supp=os.path.join(sp, "supp", "41588_2019_562_MOESM3_ESM.xlsx"),
        melt_s9=os.path.join(sp, "supp", "supp_gr.218032.116_Supplemental_Table_S9.xlsx"),
        # manual download (fetch_inputs.sh); optional — without it Tubio enters via Gardner S9
        tubio=os.path.join(sp, "supp", "tubio2014_tables.xlsx"),
    )


OPTIONAL_INPUTS = ("tubio",)


def read_genbank_seq(path):
    """Sequence of the (single-record) GenBank flat file: the ORIGIN block."""
    seq, on = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith("ORIGIN"):
                on = True
            elif line.startswith("//"):
                break
            elif on:
                seq.append("".join(c for c in line if c.isalpha()))
    return "".join(seq).upper()


def xlsx_rows(path, sheet=None):
    import warnings
    import openpyxl
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws = wb[sheet] if sheet else wb.worksheets[0]
        return [r for r in ws.iter_rows(values_only=True)]


# ============================================================================ sequence utils
def strip_polya(s, min_run=5):
    m = re.search(r"A{%d,}[ACGTN]{0,3}A*$" % min_run, s)
    return s[:m.start()] if m else s


def polya_tail_len(seq_after):
    """Length of an A-rich tail at the start of seq_after (sense), tolerating 1 non-A per 10."""
    best = 0
    na = bad = 0
    for i, c in enumerate(seq_after[:120]):
        if c == "A":
            na += 1
        else:
            bad += 1
        if bad * 10 > (i + 1) + 5:
            break
        if c == "A" and na >= 0.8 * (i + 1):
            best = i + 1
    return best


def orfs(seq, min_aa=100):
    """All ATG..stop ORFs on the sense strand: list of (start0, end_excl_incl_stop, frame)."""
    stops = {"TAA", "TAG", "TGA"}
    out = []
    for f in range(3):
        start = None
        for i in range(f, len(seq) - 2, 3):
            cod = seq[i:i + 3]
            if start is None and cod == "ATG":
                start = i
            elif start is not None and cod in stops:
                if (i + 3 - start) // 3 >= min_aa:
                    out.append((start, i + 3, f))
                start = None
    return out


def l1_orfs(seq):
    """(orf1, orf2) as (start0, end) or None. ORF1: longest ORF of 900-1300 nt starting in the
    first 1.5 kb; ORF2: longest ORF >= 3000 nt downstream of ORF1."""
    o = orfs(seq)
    o1 = [x for x in o if x[0] < 1500 and 900 <= x[1] - x[0] <= 1300]
    o1 = max(o1, key=lambda x: x[1] - x[0]) if o1 else None
    o2 = [x for x in o if x[1] - x[0] >= 3000 and (o1 is None or x[0] > o1[0])]
    o2 = max(o2, key=lambda x: x[1] - x[0]) if o2 else None
    return o1, o2


def _site_base(window, center, seq):
    """Base of `seq` aligned to window[center] (window located by infix alignment, <=6 edits)."""
    import edlib
    r = edlib.align(window, seq, mode="HW", task="path", k=6)
    if r["editDistance"] == -1:
        return None
    nice = edlib.getNiceAlignment(r, window, seq)
    qi = 0
    for qa, ta in zip(nice["query_aligned"], nice["target_aligned"]):
        if qa != "-":
            if qi == center:
                return ta
            qi += 1
    return None


def ta_status(seq):
    """Ta vs pre-Ta/older from two 3'UTR diagnostic sites (L1.3 numbering):
    5931 ('ACA' Ta vs 'ACG' pre-Ta: the Boissinot et al. 2000 'ACA at 5930-5932' diagnostic)
    and 5712 (A in Ta, G in pre-Ta and L1PA2). Both sites were confirmed to separate the Ta-0/
    Ta-1 vs pre-Ta reference sources typed by long reads in Nam et al. 2023 (Supp. Table 4:
    50/50 Ta carry A/A, 28/29 pre-Ta carry G/G). Returns Ta, nonTa, ambiguous or unknown."""
    tail = seq[-700:]
    votes = []
    for win, c in TA_SITES:
        b = _site_base(win, c, tail)
        if b in ("A", "G"):
            votes.append(b)
    if not votes:
        return "unknown"
    if all(v == "A" for v in votes):
        return "Ta"
    if all(v == "G" for v in votes):
        return "nonTa"
    return "ambiguous"


def mafft_consensus(named_seqs, work, tag, min_frac=0.5, threads=4):
    """MSA (mafft --auto) then majority-rule consensus. Gap-majority columns are dropped
    (indel-aware: an insertion present in a minority of copies does not enter the consensus).
    Returns (consensus, n_seqs, msa_path)."""
    os.makedirs(work, exist_ok=True)
    inp = os.path.join(work, "msa_%s.in.fa" % tag)
    out = os.path.join(work, "msa_%s.aln.fa" % tag)
    with open(inp, "w") as fh:
        for n, s in named_seqs:
            write_fasta(fh, n, s)
    # cache keyed on the input's md5 only (inp is rewritten every run, so an mtime test always
    # re-ran mafft, whose multithreaded output is not bit-stable for the SVA sets)
    if not (os.path.exists(out) and _same_input(inp, out + ".md5")):
        with open(out, "w") as fh:
            subprocess.run(["mafft", "--auto", "--thread", str(threads), "--quiet", inp],
                           stdout=fh, check=True)
        with open(out + ".md5", "w") as fh:
            fh.write(_md5(inp))
    aln = [s.upper() for _, _, s in read_fasta(out)]
    L = len(aln[0])
    cons = []
    for j in range(L):
        col = collections.Counter(a[j] for a in aln)
        gaps = col.pop("-", 0)
        if gaps / len(aln) >= min_frac:
            continue
        col.pop("N", None)
        if not col:
            continue
        # deterministic tie-break: count desc, then base
        b = sorted(col.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        cons.append(b)
    return "".join(cons), len(aln), out


def _md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def _same_input(inp, md5file):
    return os.path.exists(md5file) and open(md5file).read().strip() == _md5(inp)


def evenly(lst, n):
    if len(lst) <= n:
        return list(lst)
    step = len(lst) / n
    return [lst[int(i * step)] for i in range(n)]


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


class IntervalIndex:
    def __init__(self, items, key=lambda x: (x["chrom"], x["start"], x["end"])):
        self.by = collections.defaultdict(list)
        for it in items:
            c, s, e = key(it)
            self.by[c].append((s, e, it))
        for c in self.by:
            self.by[c].sort(key=lambda t: t[0])
        self.maxlen = {c: max((e - s) for s, e, _ in v) for c, v in self.by.items()}

    def query(self, c, s, e):
        import bisect
        lst = self.by.get(c, [])
        if not lst:
            return []
        lo = bisect.bisect_left(lst, (s - self.maxlen[c] - 1,))
        out = []
        for i in range(lo, len(lst)):
            a, b, it = lst[i]
            if a >= e:
                break
            if b > s:
                out.append(it)
        return out


def trim_to_element(query, target, win=20, min_match=16):
    """Element boundaries (0-based, inclusive) of `query` (a consensus) inside `target`.
    edlib's infix (HW) alignment forces the whole query to align, so a 5'-truncated element
    drags flank into the alignment; the ends are therefore pulled in to the first/last window of
    `win` alignment columns with >= `min_match` identities."""
    import edlib
    r = edlib.align(query, target, mode="HW", task="path")
    ts, te = r["locations"][0]
    nice = edlib.getNiceAlignment(r, query, target)
    q, t = nice["query_aligned"], nice["target_aligned"]
    match = [1 if a == b else 0 for a, b in zip(q, t)]
    tpos, p = [], ts                          # target index of each column (or of the next base)
    for b in t:
        tpos.append(p)
        if b != "-":
            p += 1
    n = len(match)
    i = next((k for k in range(0, n - win + 1) if sum(match[k:k + win]) >= min_match), 0)
    while i < n and not match[i]:
        i += 1
    j = next((k for k in range(n, win - 1, -1) if sum(match[k - win:k]) >= min_match), n)
    while j > 0 and not match[j - 1]:
        j -= 1
    return tpos[i], tpos[j - 1]


# ============================================================================ L1 intact
def build_l1_intact(P, hg38, lift_hs1, young38, young_hs1, l13):
    """L1Base hsflil1_8438. The exported FASTA holds the BED interval [start-1, end) (L1Base BED
    starts are 1-based) *already in element sense* (reverse-complemented for '-'), and includes
    ~1 kb of genomic flank on both sides. Verified against hg38 here, then trimmed to the L1 by an
    infix alignment of L1.3 (poly-A stripped)."""
    l13_core = strip_polya(l13)
    rep_idx = IntervalIndex([r._asdict() for r in young38 if r.name.startswith("L1")])
    rep_idx_hs1 = IntervalIndex([r._asdict() for r in young_hs1 if r.name.startswith("L1")])
    rows = []
    seqs = []
    n_orient_ok = 0
    for name, _, s in read_fasta(P["l1base_fa"]):
        s = s.upper()
        uid, loc = name.split("|")
        chrom, rest = loc.split(":")
        a, b = rest.split("(")[0].split("-")
        strand = rest[-2]
        a, b = int(a), int(b)
        g = hg38.seq(chrom, a - 1, b)
        g_sense = g if strand == "+" else revcomp(g)
        if g_sense != s:
            raise SystemExit("L1Base %s: exported sequence does not match hg38 %s:%d-%d(%s)"
                             % (uid, chrom, a, b, strand))
        n_orient_ok += 1
        ts, te = trim_to_element(l13_core, s)
        el = s[ts:te + 1]
        tail = polya_tail_len(s[te + 1:])
        # genome coords of the trimmed element (0-based half-open)
        if strand == "+":
            g0, g1 = a - 1 + ts, a - 1 + te + 1
        else:
            g1 = b - ts
            g0 = b - (te + 1)
        assert (hg38.seq(chrom, g0, g1) if strand == "+" else revcomp(hg38.seq(chrom, g0, g1))) == el
        # subfamily by rmsk overlap (largest overlapping bp among L1 records)
        ov = collections.Counter()
        for r in rep_idx.query(chrom, g0, g1):
            ov[r["name"]] += overlap(g0, g1, r["start"], r["end"])
        subfam = ov.most_common(1)[0][0] if ov else "unknown"
        o1, o2 = l1_orfs(el)
        ta = ta_status(el)
        li = lift_hs1.interval(chrom, g0, g1)
        if li:
            hc, h0, h1, flip = li
            hstrand = strand if not flip else ("-" if strand == "+" else "+")
            hov = collections.Counter()
            for r in rep_idx_hs1.query(hc, h0, h1):
                hov[r["name"]] += overlap(h0, h1, r["start"], r["end"])
            hs1_status = "present" if sum(hov.values()) >= 0.9 * (h1 - h0) and h1 - h0 > 5500 else "partial_or_absent"
        else:
            # most often a polymorphic L1 absent from the CHM13 haplotype: report the hs1
            # position of its 3' junction (hs1_start = hs1_end = 0-based junction offset)
            j = lift_hs1.junction(chrom, g1, +1) if strand == "+" else lift_hs1.junction(chrom, g0 - 1, -1)
            if j:
                flip = j[2] == "-"
                hp = (j[1] if strand == "+" else j[1] + 1) if not flip else (j[1] + 1 if strand == "+" else j[1])
                hc, h0, h1 = j[0], hp - 1, hp
                hstrand = strand if not flip else ("-" if strand == "+" else "+")
                hs1_status = "absent_in_hs1"
            else:
                hc, h0, h1, hstrand, hs1_status = ".", -1, -1, ".", "unlifted"
        young = subfam in ("L1HS", "L1PA2", "L1PA3")
        rows.append(dict(
            id=uid, subfamily=subfam,
            ta_status=("preTa" if subfam == "L1HS" else "nonTa") if ta == "nonTa" else ta,
            subfamily_call="L1HS" if (subfam == "L1HS" or ta == "Ta") else subfam,
            young="yes" if young else "no",
            hg38_chrom=chrom, hg38_start=g0 + 1, hg38_end=g1, strand=strand,
            hs1_chrom=hc, hs1_start=(h0 + 1) if h0 >= 0 else -1, hs1_end=h1, hs1_strand=hstrand,
            hs1_status=hs1_status, length=len(el),
            orf1_nt=(o1[1] - o1[0]) if o1 else 0, orf2_nt=(o2[1] - o2[0]) if o2 else 0,
            polya_tail=tail, identity_L1_3="%.4f" % identity(el, l13_core)[0]))
        seqs.append((uid, el, rows[-1]))
    log("L1Base: %d/%d exported sequences match hg38 in element-sense orientation"
        % (n_orient_ok, len(rows)))
    return rows, seqs


# ============================================================================ Alu / SVA intact
def _lift_hs1(lift_hs1, chrom, g0, g1, strand):
    li = lift_hs1.interval(chrom, g0, g1)
    if not li:
        return ".", -1, -1, "."
    hc, h0, h1, flip = li
    return hc, h0 + 1, h1, (strand if not flip else ("-" if strand == "+" else "+"))


def build_alu(young38, hg38, lift_hs1):
    cand = collections.defaultdict(list)
    n_total = collections.Counter()
    for r in young38:
        if not r.name.startswith("AluY") or r.chrom not in PRIMARY:
            continue
        n_total[r.name] += 1
        cov = r.cons_end - r.cons_begin + 1
        glen = r.end - r.start
        if (cov >= ALU_MIN_CONS_COV and r.cons_begin - 1 <= ALU_MAX_5P_MISSING
                and r.milli_div <= ALU_MAX_MILLIDIV
                and ALU_GENOMIC_LEN[0] <= glen <= ALU_GENOMIC_LEN[1]):
            cand[r.name].append(r)
    rows, seqs, passed = [], [], {}
    for name in sorted(cand):
        lst = sorted(cand[name], key=lambda r: (r.milli_div, r.chrom, r.start))
        passed[name] = len(lst)
        kept = 0
        for r in lst:
            if kept >= ALU_CAP_PER_SUBFAMILY:
                break
            g = hg38.seq(r.chrom, r.start, r.end)
            if "N" in g:
                continue
            s = g if r.strand == "+" else revcomp(g)
            after = hg38.seq(r.chrom, r.end, r.end + 80) if r.strand == "+" else \
                revcomp(hg38.seq(r.chrom, r.start - 80, r.start))
            hc, h0, h1, hs = _lift_hs1(lift_hs1, r.chrom, r.start, r.end, r.strand)
            rid = "%s_%s_%d" % (name, r.chrom, r.start + 1)
            row = dict(id=rid, subfamily=name, hg38_chrom=r.chrom, hg38_start=r.start + 1,
                       hg38_end=r.end, strand=r.strand, hs1_chrom=hc, hs1_start=h0, hs1_end=h1,
                       hs1_strand=hs, length=len(s), cons_begin=r.cons_begin, cons_end=r.cons_end,
                       cons_len=r.cons_end + r.cons_left, milli_div=r.milli_div,
                       polya_tail=polya_tail_len(after))
            rows.append(row)
            seqs.append((rid, s, row))
            kept += 1
    stats = {k: dict(rmsk_records=n_total[k], passing=passed.get(k, 0)) for k in sorted(n_total)}
    return rows, seqs, stats


def build_sva(young38, hg38, lift_hs1):
    names = ["SVA_%s" % x for x in "ABCDEF"]
    rows, seqs, stats, all_pass = [], [], {}, []
    recs = [r for r in young38 if r.chrom in PRIMARY]
    for name in names:
        chained = rm.chain_fragments(recs, {name}, max_gap=150, max_overlap=1000)
        ok = [c for c in chained
              if c["cons_end"] - c["cons_begin"] + 1 >= SVA_MIN_CONS_COV
              and c["cons_begin"] - 1 <= SVA_MAX_5P_MISSING]
        ok.sort(key=lambda c: (c["milli_div"], c["chrom"], c["start"]))
        stats[name] = dict(elements=len(chained), passing=len(ok))
        kept = 0
        for c in ok:
            g = hg38.seq(c["chrom"], c["start"], c["end"])
            if "N" in g:
                continue
            s = g if c["strand"] == "+" else revcomp(g)
            hc, h0, h1, hs = _lift_hs1(lift_hs1, c["chrom"], c["start"], c["end"], c["strand"])
            rid = "%s_%s_%d" % (name, c["chrom"], c["start"] + 1)
            row = dict(id=rid, subfamily=name, hg38_chrom=c["chrom"], hg38_start=c["start"] + 1,
                       hg38_end=c["end"], strand=c["strand"], hs1_chrom=hc, hs1_start=h0,
                       hs1_end=h1, hs1_strand=hs, length=len(s), cons_begin=c["cons_begin"],
                       cons_end=c["cons_end"], cons_len=c["cons_len"],
                       milli_div=int(round(c["milli_div"])), n_rmsk_fragments=c["n_frag"])
            all_pass.append((rid, s, row))
            if kept < SVA_CAP_PER_SUBFAMILY:
                rows.append(row)
                seqs.append((rid, s, row))
                kept += 1
    return rows, seqs, stats, all_pass


def young_fl_l1(young, names, primary_only=True, min_span=None):
    recs = [r for r in young if (r.chrom in PRIMARY or not primary_only)]
    ch = rm.chain_fragments(recs, set(names), max_gap=150)
    return [c for c in ch if c["end"] - c["start"] >= (min_span or L1_FL_MIN_SPAN)]


# ============================================================================ landmarks
def landmarks_l1(name, cons, ref=None):
    """ref = (ref_name, ref_consensus) used to transfer ORF1 when the consensus has no intact
    ORF1 (majority-rule consensus of an older subfamily can carry a frameshift)."""
    out = []
    o1, o2 = l1_orfs(cons)
    note1 = "longest 900-1300 nt ORF in first 1.5 kb (incl. stop)"
    if o1 is None and ref is not None:
        r1, _ = l1_orfs(ref[1])
        if r1:
            q2t = column_map(ref[1], cons)
            a = next((q2t[i] for i in range(r1[0], r1[1]) if q2t[i] is not None), None)
            b = next((q2t[i] for i in range(r1[1] - 1, r1[0], -1) if q2t[i] is not None), None)
            if a is not None and b is not None:
                o1 = (a, b + 1, None)
                note1 = "no intact ORF1 in this consensus; ORF1 span transferred from %s by alignment" % ref[0]
    if o1:
        out.append((name, "5UTR", 1, o1[0], "consensus start .. base before ORF1 ATG"))
        out.append((name, "ORF1", o1[0] + 1, o1[1], note1))
    if o2:
        if o1:
            out.append((name, "INTER_ORF", o1[1] + 1, o2[0], "ORF1 stop .. ORF2 ATG"))
        out.append((name, "ORF2", o2[0] + 1, o2[1], "longest >=3 kb ORF (incl. stop)"))
        out.append((name, "ORF2_EN", o2[0] + 1, o2[0] + 239 * 3,
                    "approx.: ORF2p aa 1-239 endonuclease domain (Feng et al. 1996; Weichenrieder et al. 2004)"))
        out.append((name, "ORF2_RT", o2[0] + 498 * 3 - 2, o2[0] + 773 * 3,
                    "approx.: ORF2p aa 498-773 reverse-transcriptase domain (Mathias et al. 1991)"))
        out.append((name, "3UTR", o2[1] + 1, len(cons), "ORF2 stop .. consensus end (= poly-A start)"))
    # canonical pA signal: AATAAA; at the L1 3' end it overlaps the poly-A ('AAT|AAA')
    tail = cons[-300:] + "AAA"
    i = tail.rfind("AATAAA")
    if i >= 0:
        s = len(cons) - 300 + i + 1
        out.append((name, "POLYA_SIGNAL", s, min(s + 5, len(cons)),
                    "AATAAA spanning the element/poly-A boundary (last bases are the poly-A itself)"
                    if s + 5 > len(cons) else "AATAAA"))
    import edlib
    for (win, c), lab in zip(TA_SITES, ("L1.3 5931", "L1.3 5712")):
        r = edlib.align(win, cons[-700:], mode="HW", task="locations", k=6)
        if r["locations"]:
            p = len(cons) - 700 + r["locations"][0][0] + c + 1
            out.append((name, "TA_DIAGNOSTIC", p, p,
                        "Ta/pre-Ta site (%s): A = Ta, G = pre-Ta/older; this consensus carries %s"
                        % (lab, cons[p - 1])))
    return out


def landmarks_alu(name, cons):
    out = []
    m = re.search(r"GGCTCACGCC|GGCTCATGCC", cons[:40])
    if m:
        out.append((name, "A_BOX", m.start() + 1, m.end(), "Pol III promoter A box (left monomer)"))
    m = re.search(r"GAG[AT]TCGAGAC", cons[50:110])
    if m:
        out.append((name, "B_BOX", 50 + m.start() + 1, 50 + m.end(), "Pol III promoter B box"))
    m = re.search(r"TA{4,6}TACA{4,7}", cons[100:160])
    if m:
        s = 100 + m.start() + 1
        e = 100 + m.end()
        out.append((name, "LEFT_MONOMER", 1, s - 1, "left (FLAM-C derived) monomer"))
        out.append((name, "A_RICH_LINKER", s, e, "middle A-rich linker"))
        out.append((name, "RIGHT_MONOMER", e + 1, len(cons), "right monomer (ends at the poly-A)"))
    return out


def landmarks_sva(name, cons):
    out = []
    j = cons.find("CCCTCT")
    if 0 <= j < 100:
        # (CCCTCT)n region: extend while the next 12-mer is >= 10/12 C/T (tolerates the
        # CCCCCT / CTCCCT variants of the hexamer)
        k = j
        while k + 12 <= len(cons) and sum(ch in "CT" for ch in cons[k:k + 12]) >= 10:
            k += 1
        out.append((name, "HEXAMER", j + 1, k + 11, "(CCCTCT)n hexamer repeat at the SVA 5' end (C/T-rich run)"))
    i = (cons[-80:] + "AAA").rfind("AATAAA")
    if i >= 0:
        s = len(cons) - 80 + i + 1
        out.append((name, "POLYA_SIGNAL", s, min(s + 5, len(cons)),
                    "AATAAA in SINE-R" + (" (spans the element/poly-A boundary)" if s + 5 > len(cons) else "")))
    out.append((name, "SINE_R", max(1, len(cons) - 490 + 1), len(cons),
                "approx.: last ~490 bp = HERV-K10 LTR-derived SINE-R (Wang et al. 2005)"))
    return out


# ============================================================================ writers
def write_tsv(path, rows, cols=None):
    cols = cols or list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, ".")) for c in cols) + "\n")


def write_fa(path, recs):
    with open(path, "w") as fh:
        for rec in recs:
            name, seq = rec[0], rec[1]
            desc = rec[2] if len(rec) > 2 and isinstance(rec[2], str) else ""
            write_fasta(fh, name, seq, desc)


def bgzip_index(path):
    import pysam
    pysam.tabix_compress(path, path + ".gz", force=True)
    os.remove(path)
    pysam.faidx(path + ".gz")


def manifest(out):
    import gzip
    rows = []
    for fn in sorted(os.listdir(out)):
        p = os.path.join(out, fn)
        if not os.path.isfile(p) or fn in ("manifest.tsv", "README.md"):
            continue
        n = "."
        if fn.endswith(".fa"):
            n = sum(1 for l in open(p) if l.startswith(">"))
        elif fn.endswith(".fa.gz"):
            n = sum(1 for l in gzip.open(p, "rt") if l.startswith(">"))
        elif fn.endswith(".tsv"):
            n = sum(1 for _ in open(p)) - 1
        rows.append(dict(file=fn, records=n, bytes=os.path.getsize(p), md5=_md5(p)))
    write_tsv(os.path.join(out, "manifest.tsv"), rows)
    return rows


# ============================================================================ main
def main(argv=None):
    from sources import (parse_published, transduction_stats, build_sources, parse_tubio_s3,
                         validate_tubio_s3, polymorphic_candidates, STATS_COLUMNS, POLY_COLUMNS)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", required=True, help="directory laid out by fetch_inputs.sh")
    ap.add_argument("--work", required=True, help="scratch directory for caches / MSAs")
    ap.add_argument("--out", required=True, help="output directory (resources/rte_library)")
    for k in default_paths("X"):
        ap.add_argument("--" + k.replace("_", "-"), dest=k, help="override input path for %s" % k)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    P = default_paths(a.inputs)
    for k in P:
        v = getattr(a, k)
        if v:
            P[k] = v
    for k in OPTIONAL_INPUTS:
        if not os.path.exists(P[k]):
            log("optional input %s not found (%s): skipped" % (k, P[k]))
            P[k] = None
    missing = [k for k, v in P.items() if v and not os.path.exists(v)]
    if missing:
        raise SystemExit("missing inputs: %s" % ", ".join("%s=%s" % (k, P[k]) for k in missing))
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)

    hg38 = TwoBit(P["hg38_2bit"])
    hs1 = TwoBit(P["hs1_2bit"])
    lift_hs1 = Lifter(P["chain_hg38_hs1"])
    lift19 = Lifter(P["chain_hg19_hg38"])
    bands38 = load_cytobands(P["cyto_hg38"])
    log("loading RepeatMasker annotations")
    young38, _ = rm.load(P["hg38_rmsk"], a.work, "hg38")
    young_hs1, mask_hs1 = rm.load(P["hs1_rmsk"], a.work, "hs1")
    l13 = read_genbank_seq(P["l13_gb"])
    l12 = read_genbank_seq(P["l12_gb"])
    l1rp = read_genbank_seq(P["l1rp_gb"])
    dfam = {n.split("|")[0]: s.upper() for n, _, s in read_fasta(P["dfam"])}

    log("L1 intact (L1Base hsflil1_8438)")
    l1rows, l1seqs = build_l1_intact(P, hg38, lift_hs1, young38, young_hs1, l13)
    log("AluY intact")
    alurows, aluseqs, alustats = build_alu(young38, hg38, lift_hs1)
    log("SVA intact")
    svarows, svaseqs, svastats, sva_all = build_sva(young38, hg38, lift_hs1)

    # ---- consensus
    log("consensus MSAs (mafft)")
    cons = collections.OrderedDict()
    nseq = {}
    l1hs = sorted([(u, s) for u, s, r in l1seqs if r["subfamily"] == "L1HS" and r["ta_status"] == "Ta"],
                  key=lambda x: int(x[0].split("-")[1]))
    cons["L1HS"], nseq["L1HS"], _ = mafft_consensus(evenly(l1hs, MSA_N_L1), a.work, "L1HS", threads=a.threads)
    for fam in ("L1PA2", "L1PA3"):
        fl = [f for f in young_fl_l1(young38, [fam]) if f["cons_begin"] <= 20]
        fl.sort(key=lambda f: (f["milli_div"], f["chrom"], f["start"]))
        seqs = []
        for f in fl[:MSA_N_L1]:
            g = hg38.seq(f["chrom"], f["start"], f["end"])
            seqs.append(("%s_%s_%d" % (fam, f["chrom"], f["start"] + 1), g if f["strand"] == "+" else revcomp(g)))
        cons[fam], nseq[fam], _ = mafft_consensus(seqs, a.work, fam, threads=a.threads)
    for cname, sub in (("ALU_Y", "AluY"), ("ALU_YA5", "AluYa5"), ("ALU_YB8", "AluYb8")):
        seqs = [(u, s) for u, s, r in aluseqs if r["subfamily"] == sub][:MSA_N_ALU]
        cons[cname], nseq[cname], _ = mafft_consensus(seqs, a.work, cname, threads=a.threads)
    for cname in ("SVA_D", "SVA_E", "SVA_F"):
        seqs = [(u, s) for u, s, r in svaseqs if r["subfamily"] == cname][:MSA_N_SVA]
        if len(seqs) >= 10:
            cons[cname], nseq[cname], _ = mafft_consensus(seqs, a.work, cname, threads=a.threads)

    for k in cons:
        cons[k] = strip_polya(cons[k])

    # ---- cross-check
    xc = []
    pairs = [("L1HS", "L1.3 (L19088)", strip_polya(l13), "NW"),
             ("L1HS", "L1.2 (M80343)", strip_polya(l12), "NW"),
             ("L1HS", "L1RP (AF148856)", strip_polya(l1rp), "NW"),
             ("L1HS", "Dfam L1HS_5end", dfam["L1HS_5end"], "HW"),
             ("L1HS", "Dfam L1HS_3end", strip_polya(dfam["L1HS_3end"]), "HW"),
             ("L1PA2", "Dfam L1PA2_3end", strip_polya(dfam["L1PA2_3end"]), "HW"),
             ("L1PA3", "Dfam L1PA3_3end", strip_polya(dfam["L1PA3_3end"]), "HW"),
             ("L1PA2", "L1.3 (L19088)", strip_polya(l13), "NW"),
             ("ALU_Y", "Dfam AluY", dfam["AluY"], "NW"),
             ("ALU_YA5", "Dfam AluYa5", dfam["AluYa5"], "NW"),
             ("ALU_YB8", "Dfam AluYb8", dfam["AluYb8"], "NW")]
    for cname in ("SVA_D", "SVA_E", "SVA_F"):
        if cname in cons:
            pairs.append((cname, "Dfam " + cname, dfam[cname], "NW"))
    for cname, ref, rseq, mode in pairs:
        if cname not in cons:
            continue
        idt, ed, loc = identity(strip_polya(rseq), cons[cname], mode=mode)
        ung = ungapped_identity(strip_polya(rseq), cons[cname], mode=mode)
        xc.append(dict(consensus=cname, reference=ref, mode=mode, identity="%.4f" % idt,
                       identity_ungapped="%.4f" % ung,
                       edit_distance=ed, ref_len=len(rseq), consensus_len=len(cons[cname]),
                       consensus_span=("%d-%d" % (loc[0] + 1, loc[1] + 1)) if loc[0] is not None else "."))
    for u, s, r in l1seqs:
        r["identity_L1HS_consensus"] = "%.4f" % cons_identity(s, cons["L1HS"])
        _, _, (c0, c1) = identity(s, cons["L1HS"], mode="HW")
        r["l1hs_cons_start"], r["l1hs_cons_end"] = c0 + 1, c1 + 1

    lm = []
    for cname, cs in cons.items():
        if cname.startswith("L1"):
            lm += landmarks_l1(cname, cs, ref=("L1HS", cons["L1HS"]) if cname != "L1HS" else None)
        elif cname.startswith("ALU"):
            lm += landmarks_alu(cname, cs)
        elif cname.startswith("SVA"):
            lm += landmarks_sva(cname, cs)

    # ---- sources
    log("transduction sources")
    ents, td_geom = parse_published(P)
    src_rows, fl3, fl5 = build_sources(ents, hg38, hs1, lift19, lift_hs1, young38, young_hs1,
                                       mask_hs1, l1rows, cons["L1HS"], sva_all, bands38)
    tubio_tds, poly = None, None
    if P["tubio"]:
        tubio_tds = validate_tubio_s3(parse_tubio_s3(P["tubio"]), ents, src_rows, fl3, lift19, lift_hs1,
                                      hg38=hg38)
        log("Tubio S3 transductions vs flanks_3p: %s" % dict(collections.Counter(
            (t["type"], t["flank_hit"]) for t in tubio_tds)))
        poly = polymorphic_candidates(P["tubio"], lift19, lift_hs1, young38, young_hs1, src_rows,
                                      hs1, cons["L1HS"])
    tstats = transduction_stats(td_geom, tubio_tds)

    # ---- active
    act = []
    srcs_by_l1b = {r["l1base_id"]: r for r in src_rows if r["l1base_id"] != "."}
    for r in l1rows:
        if r["subfamily_call"] != "L1HS":
            continue
        s = srcs_by_l1b.get(r["id"])
        tier = "hot_source" if s and s["hotness"] in ("hot", "strong") else (
            "L1HS_Ta_intact" if r["ta_status"] == "Ta" else "L1HS_preTa_intact")
        act.append(dict(id=r["id"], element_class="L1", subfamily=r["subfamily"],
                        ta_status=r["ta_status"], tier=tier, reference="yes",
                        hg38_chrom=r["hg38_chrom"], hg38_start=r["hg38_start"], hg38_end=r["hg38_end"],
                        strand=r["strand"], hs1_chrom=r["hs1_chrom"], hs1_start=r["hs1_start"],
                        hs1_end=r["hs1_end"], hs1_strand=r["hs1_strand"], orf_intact="yes",
                        consensus="L1HS", identity_to_consensus=r["identity_L1HS_consensus"],
                        source_id=s["id"] if s else ".", n_daughters=s["n_daughters"] if s else 0,
                        hotness=s["hotness"] if s else "none_reported", in_fasta="l1_intact.fa"))
    seen = {x["source_id"] for x in act}
    for s in src_rows:
        if s["id"] in seen or s["element_class"] != "L1" or s["hotness"] not in ("hot", "strong"):
            continue
        act.append(dict(id=s["id"], element_class="L1", subfamily=s["subfamily"], ta_status=s["ta_status"],
                        tier="hot_source", reference=s["reference"], hg38_chrom=s["hg38_chrom"],
                        hg38_start=s["hg38_start"], hg38_end=s["hg38_end"], strand=s["strand"],
                        hs1_chrom=s["hs1_chrom"], hs1_start=s["hs1_start"], hs1_end=s["hs1_end"],
                        hs1_strand=s["hs1_strand"], orf_intact=s["orf_intact"], consensus="L1HS",
                        identity_to_consensus=s["identity_L1HS"], source_id=s["id"],
                        n_daughters=s["n_daughters"], hotness=s["hotness"], in_fasta="."))

    # ---- write
    out = a.out
    l1cols = ["id", "subfamily", "subfamily_call", "ta_status", "young", "hg38_chrom", "hg38_start", "hg38_end", "strand",
              "hs1_chrom", "hs1_start", "hs1_end", "hs1_strand", "hs1_status", "length", "orf1_nt",
              "orf2_nt", "polya_tail", "identity_L1_3", "identity_L1HS_consensus", "l1hs_cons_start",
              "l1hs_cons_end"]
    write_tsv(os.path.join(out, "l1_intact.tsv"), l1rows, l1cols)
    write_fa(os.path.join(out, "l1_intact.fa"),
             [(u, s, "%s %s hg38=%s:%d-%d(%s) hs1=%s:%s-%s(%s)" % (
                 r["subfamily"], r["ta_status"], r["hg38_chrom"], r["hg38_start"], r["hg38_end"],
                 r["strand"], r["hs1_chrom"], r["hs1_start"], r["hs1_end"], r["hs1_strand"]))
              for u, s, r in l1seqs])
    write_tsv(os.path.join(out, "alu_y_intact.tsv"), alurows)
    write_fa(os.path.join(out, "alu_y_intact.fa"),
             [(u, s, "%s hg38=%s:%d-%d(%s) milliDiv=%d" % (
                 r["subfamily"], r["hg38_chrom"], r["hg38_start"], r["hg38_end"], r["strand"], r["milli_div"]))
              for u, s, r in aluseqs])
    write_tsv(os.path.join(out, "sva_intact.tsv"), svarows)
    write_fa(os.path.join(out, "sva_intact.fa"),
             [(u, s, "%s hg38=%s:%d-%d(%s) milliDiv=%d" % (
                 r["subfamily"], r["hg38_chrom"], r["hg38_start"], r["hg38_end"], r["strand"], r["milli_div"]))
              for u, s, r in svaseqs])
    write_fa(os.path.join(out, "consensus.fa"),
             [(k, v, "majority-rule consensus of %d sense-oriented intact copies (mafft --auto); no poly-A" % nseq[k])
              for k, v in cons.items()])
    write_tsv(os.path.join(out, "consensus_crosscheck.tsv"), xc)
    write_tsv(os.path.join(out, "consensus_landmarks.tsv"),
              [dict(consensus=x[0], feature=x[1], start=x[2], end=x[3], note=x[4]) for x in lm])
    young_dfam = [k for k in dfam if re.match(r"^(Alu|L1HS|L1PA|SVA)", k)]
    write_fa(os.path.join(out, "dfam_young.fa"),
             [(k, dfam[k], "Dfam consensus (CC0); cross-check only") for k in young_dfam])
    write_tsv(os.path.join(out, "active.tsv"), act)
    from sources import SOURCE_COLUMNS
    write_tsv(os.path.join(out, "transduction_sources.tsv"), src_rows, SOURCE_COLUMNS)
    write_fa(os.path.join(out, "flanks_3p.fa"), fl3)
    bgzip_index(os.path.join(out, "flanks_3p.fa"))
    write_fa(os.path.join(out, "flanks_5p_sva.fa"), fl5)
    bgzip_index(os.path.join(out, "flanks_5p_sva.fa"))
    write_tsv(os.path.join(out, "transduction_stats.tsv"), tstats, STATS_COLUMNS)
    if poly is not None:
        write_tsv(os.path.join(out, "polymorphic_l1_candidates.tsv"), poly, POLY_COLUMNS)
    if tubio_tds is not None:
        write_tsv(os.path.join(a.work, "tubio2014_s3_validation.tsv"), tubio_tds,
                  ["type", "sample", "chrom", "start", "end", "strand", "ts", "te", "td_len",
                   "distal", "source_id", "flank_hit", "flank_offset"])
    mf = manifest(out)
    summary = dict(l1_intact=len(l1rows),
                   l1_subfamily=collections.Counter((r["subfamily"], r["ta_status"]) for r in l1rows).most_common(),
                   alu=alustats, sva=svastats, alu_kept=len(alurows), sva_kept=len(svarows),
                   consensus={k: len(v) for k, v in cons.items()}, msa_n=nseq,
                   sources=collections.Counter((r["element_class"], r["reference"]) for r in src_rows).most_common(),
                   hotness=collections.Counter(r["hotness"] for r in src_rows).most_common(),
                   active=len(act), total_bytes=sum(r["bytes"] for r in mf))
    log(json.dumps(summary, default=str, indent=1))
    with open(os.path.join(a.work, "build_summary.json"), "w") as fh:
        json.dump(summary, fh, default=str, indent=1)


if __name__ == "__main__":
    main()
