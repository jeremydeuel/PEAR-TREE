"""Real element sequences and real source flanks for the insertion simulators.

Two sources, same in-memory model:

1. `--rte-library DIR` — the committed `resources/rte_library/` (SPEC.md): `l1_intact.fa`,
   `alu_y_intact.fa`, `sva_intact.fa`, `transduction_sources.tsv`, `flanks_3p.fa`,
   `flanks_5p_sva.fa`. Parsed leniently (record name = first token; `id|...`/`id::...`
   suffixes ignored; TSV columns looked up by common names).
2. Fallback: extract young full-length copies from hs1 + its RepeatMasker `.out` track
   (L1HS >= 5.9 kb starting at Dfam consensus pos <= 150, fragments merged by rmsk ID; AluYa5/AluYb8 280-320 bp; SVA_E/SVA_F
   >= 1.3 kb; lowest divergence first) together with their REAL downstream (3') flanks
   (15 kb, element sense) and, for SVA, upstream (5') flanks (3 kb). The extraction is
   cached in exactly the rte_library file layout, so the cache dir can be passed back as
   `--rte-library`.

Transduction termination: each source gets 1-3 fixed termination offsets in its 3' flank,
placed ~10-30 bp after a canonical poly-A signal (AATAAA/ATTAAA), weighted towards short
tags (most < 1 kb, rare up to 12 kb) — Zumalave 2026 (per-source clustering at e.g. 234 /
64 bp), Tubio 2014 (length distribution). Deterministic per source id.
"""
import csv
import gzip
import hashlib
import math
import os
import random
import re
from dataclasses import dataclass, field

from .seqs import Genome, read_fasta, revcomp, write_fasta


def split_polya_tail(seq, win=10, min_a=7):
    """(body, tail): strip the genomic A-rich 3' tail RepeatMasker includes in young
    elements (e.g. ...AATAAAAAAAAAAGAAAA) so the simulator can add its own poly-A."""
    i = len(seq)
    while i >= win and seq[i - win:i].count("A") >= min_a:
        i -= 1
    # move back to the start of that A-rich stretch
    while i < len(seq) and seq[i] != "A":
        i += 1
    if len(seq) - i < 5:
        return seq, ""
    return seq[:i], seq[i:]


@dataclass
class Element:
    id: str
    cls: str              # L1 | ALU | SVA
    subfamily: str
    seq: str              # sense orientation, upper case, genomic poly-A tail stripped
    loc: str = "."        # hs1 contig:start-end:strand (documentation)
    tail: str = ""        # the stripped genomic A-rich tail of the reference copy

    def __post_init__(self):
        if not self.tail:
            self.seq, self.tail = split_polya_tail(self.seq)


@dataclass
class Source:
    id: str
    cls: str              # L1 | SVA
    subfamily: str
    element_seq: str      # sense
    flank3: str           # downstream flank, sense, upper case
    flank5: str = ""      # upstream flank, sense (SVA 5' transduction)
    loc: str = "."
    endpoints: list = field(default_factory=list)
    tail: str = ""        # source's own genomic A-rich tail (between element and flank3)

    def __post_init__(self):
        if self.element_seq and not self.tail:
            self.element_seq, self.tail = split_polya_tail(self.element_seq)


def _key(name):
    return re.split(r"[|:\s]", name, maxsplit=1)[0]


def _poly_a_signal_endpoints(sid, flank, max_len=12000):
    """1-3 deterministic termination offsets in `flank` for source `sid`."""
    rng = random.Random(int(hashlib.md5(sid.encode()).hexdigest()[:12], 16))
    f = flank[:max_len].upper()
    sig = [m.start() for m in re.finditer(r"(?=(AATAAA|ATTAAA))", f) if m.start() >= 20]
    cands = []
    for p in sig:
        end = p + 6 + rng.randint(10, 30)
        if end < len(f):
            # weight short tags: most transductions < 1 kb (Tubio 2014)
            cands.append((end, math.exp(-end / 700.0) + 0.002))
    if not cands:
        cands = [(min(len(f), x), 1.0) for x in (rng.randint(60, 400), rng.randint(400, 1200))]
    k = rng.choice([1, 1, 2, 2, 3])
    chosen = []
    pool = list(cands)
    for _ in range(min(k, len(pool))):
        tot = sum(w for _, w in pool)
        u, acc = rng.random() * tot, 0.0
        for i, (e, w) in enumerate(pool):
            acc += w
            if u <= acc:
                chosen.append(e)
                pool.pop(i)
                break
    return sorted(chosen)


class RteLibrary:
    def __init__(self):
        self.elements = {"L1": [], "ALU": [], "SVA": []}
        self.sources = {"L1": [], "SVA": []}
        self.origin = "?"

    # ---------------------------------------------------------------- selection
    def element(self, rng, cls, subfamily=None):
        pool = self.elements[cls]
        if subfamily:
            sub = [e for e in pool if e.subfamily.upper().startswith(subfamily.upper())]
            pool = sub or pool
        if not pool:
            raise RuntimeError(f"rte library has no {cls} elements")
        return rng.choice(pool)

    def source(self, rng, cls="L1", need_flank5=False, need_element=False):
        pool = [s for s in self.sources[cls] if s.flank3 and (s.flank5 or not need_flank5)
                and (s.element_seq or not need_element)]
        if not pool:
            raise RuntimeError(f"rte library has no {cls} transduction sources")
        return rng.choice(pool)

    def summary(self):
        return (f"rte library ({self.origin}): " +
                ", ".join(f"{k}={len(v)}" for k, v in self.elements.items()) +
                f"; sources L1={len(self.sources['L1'])} SVA={len(self.sources['SVA'])}")

    # ---------------------------------------------------------------- loaders
    @classmethod
    def from_dir(cls, d):
        lib = cls()
        lib.origin = d
        meta = {}
        for tsvname in ("l1_intact.tsv", "alu_y_intact.tsv", "sva_intact.tsv", "transduction_sources.tsv"):
            p = os.path.join(d, tsvname)
            if os.path.exists(p):
                with open(p) as f:
                    rows = [r for r in csv.DictReader((l for l in f if not l.startswith("##")), delimiter="\t")]
                for r in rows:
                    rid = r.get("id") or r.get("source_id") or r.get("name") or next(iter(r.values()))
                    meta[_key(str(rid).lstrip("#"))] = r
        for fname, klass in (("l1_intact.fa", "L1"), ("alu_y_intact.fa", "ALU"), ("sva_intact.fa", "SVA")):
            p = os.path.join(d, fname)
            if not os.path.exists(p):
                continue
            for name, seq in read_fasta(p).items():
                k = _key(name)
                m = meta.get(k, {})
                sub = m.get("subfamily") or m.get("repName") or m.get("family") or name
                lib.elements[klass].append(Element(k, klass, str(sub), seq.upper(),
                                                   m.get("hs1", m.get("hs1_coords", "."))))
        def _fa(name):
            for cand in (name, name + ".gz"):      # resources/rte_library ships bgzip .fa.gz
                pth = os.path.join(d, cand)
                if os.path.exists(pth):
                    return pth
            return None
        f3, f5 = _fa("flanks_3p.fa"), _fa("flanks_5p_sva.fa")
        # strandless sources ship two candidate flanks `<id>/+`, `<id>/-`: the simulator needs
        # an unambiguous source orientation, so those are not used as simulated sources
        fl3 = {_key(n): s.upper() for n, s in read_fasta(f3).items() if "/" not in n} if f3 else {}
        fl5 = {_key(n): s.upper() for n, s in read_fasta(f5).items() if "/" not in n} if f5 else {}
        by_id = {e.id: e for v in lib.elements.values() for e in v}
        # hs1 intervals of the intact elements (real library: hs1_chrom/hs1_start/hs1_end,
        # 1-based) -> a reference source without an explicit element link finds its copy
        iv = {}
        for e in by_id.values():
            m = meta.get(e.id, {})
            try:
                c, a, b = m["hs1_chrom"], int(m["hs1_start"]), int(m["hs1_end"])
            except (KeyError, ValueError):
                continue
            if c not in (".", "") and b > a:
                iv.setdefault(c, []).append((a, b, e))

        def _overlap_element(m):
            try:
                c, a, b = m["hs1_chrom"], int(m["hs1_start"]), int(m["hs1_end"])
            except (KeyError, ValueError):
                return None
            best = None
            for x, y, e in iv.get(c, ()):
                ov = min(b, y) - max(a, x)
                if ov > 0.5 * (b - a) and (best is None or ov > best[0]):
                    best = (ov, e)
            return best[1] if best else None

        for sid, flank in fl3.items():
            m = meta.get(sid, {})
            klass = (m.get("class") or m.get("element_class") or m.get("cls") or "").upper()
            el = by_id.get(sid)
            for k in ("intact_id", "l1base_id"):      # source -> its intact element (real library)
                if el is None and m.get(k) not in (None, "", "."):
                    el = by_id.get(m[k])
            if el is None and m:
                el = _overlap_element(m)
            if not klass:
                klass = el.cls if el else ("SVA" if "SVA" in sid.upper() else "L1")
            klass = "SVA" if "SVA" in klass else ("L1" if "L1" in klass or "LINE" in klass else klass)
            if klass not in lib.sources:
                continue
            src = Source(sid, klass, el.subfamily if el else str(m.get("subfamily", klass)),
                         (el.seq + el.tail) if el else "", flank, fl5.get(sid, ""),
                         str(m.get("hs1", m.get("hs1_coords", "."))))
            src.endpoints = _poly_a_signal_endpoints(sid, flank)
            lib.sources[klass].append(src)
        return lib

    @classmethod
    def from_hs1(cls, genome_path, rmsk_path, cache_dir=None, n_l1=40, n_alu=40, n_sva=24,
                 flank3_len=15000, flank5_len=3000):
        if cache_dir and os.path.exists(os.path.join(cache_dir, "l1_intact.fa")):
            lib = cls.from_dir(cache_dir)
            lib.origin = f"hs1 cache {cache_dir}"
            return lib
        g = Genome(genome_path)
        want = {"L1HS": [], "AluYa5": [], "AluYb8": [], "SVA_E": [], "SVA_F": []}
        l1_frag = {}
        with gzip.open(rmsk_path, "rt") as f:
            for line in f:
                t = line.split()
                if len(t) < 15 or t[9] not in want:
                    continue
                try:
                    div = float(t[1]); s = int(t[5]) - 1; e = int(t[6])
                except ValueError:
                    continue
                strand = "+" if t[8] == "+" else "-"
                L = e - s
                name = t[9]
                if t[4] not in g.lengths or "_" in t[4] or t[4] in ("chrM", "chrY"):
                    continue
                rbeg = int(t[11]) if strand == "+" else int(t[13])
                if name == "L1HS":                     # merge fragments of one element by rmsk ID
                    k = (t[4], t[14], strand)
                    d = l1_frag.setdefault(k, [s, e, rbeg, div])
                    d[0] = min(d[0], s); d[1] = max(d[1], e); d[2] = min(d[2], rbeg); d[3] = min(d[3], div)
                    continue
                if name.startswith("Alu") and not (280 <= L <= 320):
                    continue
                if name.startswith("SVA") and L < 1300:
                    continue
                want[name].append((div, t[4], s, e, strand))
        # full-length L1HS: >= 5.9 kb and starting at the 5' end of the (124-based) Dfam
        # L1HS consensus as RepeatMasker reports it
        for (c, _id, strand), (s, e, rbeg, div) in l1_frag.items():
            if e - s >= 5900 and rbeg <= 150:
                want["L1HS"].append((div, c, s, e, strand))
        lib = cls()
        lib.origin = f"hs1 rmsk {rmsk_path}"
        fa_el = {"L1": {}, "ALU": {}, "SVA": {}}
        tsv_el = {"L1": [], "ALU": [], "SVA": []}
        fl3, fl5, src_rows = {}, {}, []
        plan = [("L1HS", "L1", n_l1), ("AluYa5", "ALU", n_alu // 2), ("AluYb8", "ALU", n_alu // 2),
                ("SVA_E", "SVA", n_sva // 2), ("SVA_F", "SVA", n_sva // 2)]
        for name, klass, n in plan:
            for div, c, s, e, strand in sorted(want[name])[:n]:
                eid = f"{name}_{c}_{s}"
                seq = g.fetch(c, s, e, strand)
                loc = f"{c}:{s}-{e}:{strand}"
                fa_el[klass][eid] = seq
                tsv_el[klass].append({"id": eid, "subfamily": name, "class": klass, "hs1": loc, "div": div})
                if klass in ("L1", "SVA"):
                    if strand == "+":
                        f3 = g.fetch(c, e, e + flank3_len, "+")
                        f5 = g.fetch(c, s - flank5_len, s, "+")
                    else:
                        f3 = g.fetch(c, s - flank3_len, s, "-")
                        f5 = g.fetch(c, e, e + flank5_len, "-")
                    fl3[eid] = f3
                    if klass == "SVA":
                        fl5[eid] = f5
                    src_rows.append({"id": eid, "class": klass, "subfamily": name, "hs1": loc,
                                     "strand": strand, "evidence": "hs1_rmsk_fallback"})
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
            for fname, klass in (("l1_intact", "L1"), ("alu_y_intact", "ALU"), ("sva_intact", "SVA")):
                write_fasta(os.path.join(cache_dir, fname + ".fa"), fa_el[klass])
                _write_tsv(os.path.join(cache_dir, fname + ".tsv"), tsv_el[klass])
            write_fasta(os.path.join(cache_dir, "flanks_3p.fa"), fl3)
            write_fasta(os.path.join(cache_dir, "flanks_5p_sva.fa"), fl5)
            _write_tsv(os.path.join(cache_dir, "transduction_sources.tsv"), src_rows)
            lib = cls.from_dir(cache_dir)
            lib.origin = f"hs1 rmsk (cached {cache_dir})"
            return lib
        # no cache: assemble in memory
        for klass in fa_el:
            for row in tsv_el[klass]:
                lib.elements[klass].append(Element(row["id"], klass, row["subfamily"],
                                                   fa_el[klass][row["id"]], row["hs1"]))
        for row in src_rows:
            el = fa_el[row["class"]][row["id"]]
            src = Source(row["id"], row["class"], row["subfamily"], el, fl3[row["id"]],
                         fl5.get(row["id"], ""), row["hs1"])
            src.endpoints = _poly_a_signal_endpoints(src.id, src.flank3)
            lib.sources[row["class"]].append(src)
        return lib

    @classmethod
    def synthetic(cls, seed=0xE1E5):
        """Last-resort stand-in (no genome available): random sequences with the right
        lengths. Only for unit tests / genome-free runs without hs1."""
        rng = random.Random(seed)
        from .seqs import rnd_seq
        lib = cls()
        lib.origin = "synthetic"
        for i in range(4):
            lib.elements["L1"].append(Element(f"synL1_{i}", "L1", "L1HS", rnd_seq(rng, 6030)))
            lib.elements["ALU"].append(Element(f"synAluYa5_{i}", "ALU", "AluYa5", rnd_seq(rng, 300)))
            lib.elements["ALU"].append(Element(f"synAluYb8_{i}", "ALU", "AluYb8", rnd_seq(rng, 310)))
            lib.elements["SVA"].append(Element(f"synSVA_E_{i}", "SVA", "SVA_E", rnd_seq(rng, 1950)))
            lib.elements["SVA"].append(Element(f"synSVA_F_{i}", "SVA", "SVA_F", rnd_seq(rng, 2100)))
        for klass in ("L1", "SVA"):
            for e in lib.elements[klass][:4]:
                fl = list(rnd_seq(rng, 15000))
                for p in (180, 900, 4000):           # plant poly-A signals
                    fl[p:p + 6] = "AATAAA"
                src = Source(e.id, klass, e.subfamily, e.seq, "".join(fl), rnd_seq(rng, 3000))
                src.endpoints = _poly_a_signal_endpoints(src.id, src.flank3)
                lib.sources[klass].append(src)
        return lib


def _write_tsv(path, rows):
    if not rows:
        open(path, "w").close()
        return
    keys = list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in keys) + "\n")


def load_library(rte_library=None, genome=None, rmsk=None, cache_dir=None):
    """--rte-library DIR if given and populated; else hs1+rmsk extraction (cached);
    else the synthetic stand-in."""
    if rte_library and os.path.exists(os.path.join(rte_library, "l1_intact.fa")):
        return RteLibrary.from_dir(rte_library)
    if genome and rmsk and os.path.exists(genome) and os.path.exists(rmsk):
        if cache_dir is None:
            cache_dir = os.path.join(os.path.dirname(os.path.abspath(genome)), "simlib_rte_cache")
        return RteLibrary.from_hs1(genome, rmsk, cache_dir)
    return RteLibrary.synthetic()


# ------------------------------------------------------------------------------------
# gene models (processed pseudogene parents / co-inserted pre-mRNA)
# ------------------------------------------------------------------------------------
@dataclass
class Gene:
    id: str
    contig: str
    strand: str
    exons: list            # [(start, end)] 0-based half-open, genomic order
    exon_seqs: list        # sense-oriented exon sequences, transcript order
    premrna: str = ""      # sense-oriented unspliced sequence (first..last exon)

    def mrna(self):
        return "".join(self.exon_seqs)


def gene_from_sequence(gid, contig, seq, offset, rng, n_exons=None, strand="+"):
    """Build a gene model on real sequence `seq` (genome slice starting at `offset`)
    with canonical splice signals: every exon is followed by an intron starting `GT` and
    ending `AG`. Used when no annotation is available (fallback; documented)."""
    n_exons = n_exons or rng.randint(3, 6)
    pos, exons = rng.randint(100, 500), []
    for i in range(n_exons):
        elen = rng.randint(90, 260)
        end = pos + elen
        if i < n_exons - 1:
            q = seq.find("GT", end)
            if q < 0 or q - end > 60:
                q = end
            end = q
            exons.append((pos, end))
            ilen = rng.randint(300, 3000)
            r = seq.find("AG", end + ilen)
            if r < 0 or r + 2 + 300 > len(seq):
                break
            pos = r + 2
        else:
            exons.append((pos, min(end, len(seq))))
    if len(exons) < 2:
        return None
    exon_seqs = [seq[a:b] for a, b in exons]
    pre = seq[exons[0][0]:exons[-1][1]]
    gexons = [(offset + a, offset + b) for a, b in exons]
    if strand == "-":
        exon_seqs = [revcomp(s) for s in exon_seqs[::-1]]
        pre = revcomp(pre)
    return Gene(gid, contig, strand, gexons, exon_seqs, pre)


def load_gene_model(path, genome, min_exons=2, max_genes=2000):
    """Genes from a `contig start end gene strand` merged-exon TSV (tools/build_gene_model.py
    output, 0-based half-open) or a GTF (exon features). Contig names are matched to the
    genome with/without a chr prefix."""
    by_gene = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) >= 9 and t[2] == "exon":
                m = re.search(r'gene_name "([^"]+)"', t[8]) or re.search(r'gene_id "([^"]+)"', t[8])
                gid = m.group(1) if m else t[8][:30]
                c, s, e, st = t[0], int(t[3]) - 1, int(t[4]), t[6]
            elif len(t) >= 5:
                c, s, e, gid, st = t[0], int(t[1]), int(t[2]), t[3], t[4]
            else:
                continue
            by_gene.setdefault((gid, c, st), set()).add((s, e))
    genes = []
    for (gid, c, st), ex in by_gene.items():
        cc = c if c in genome.lengths else ("chr" + c if "chr" + c in genome.lengths else c.replace("chr", ""))
        if cc not in genome.lengths:
            continue
        ex = sorted(ex)
        if len(ex) < min_exons or ex[-1][1] - ex[0][0] > 300_000:
            continue
        seqs = [genome.fetch(cc, a, b) for a, b in ex]
        pre = genome.fetch(cc, ex[0][0], ex[-1][1])
        if st == "-":
            seqs = [revcomp(s) for s in seqs[::-1]]
            pre = revcomp(pre)
        genes.append(Gene(gid, cc, st, ex, seqs, pre))
        if len(genes) >= max_genes:
            break
    return genes
