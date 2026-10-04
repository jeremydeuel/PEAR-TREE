"""Reference RTE library loader + minimap2 (mappy) indices.

Layout (plans/tprt_hallmarks/SPEC.md, built by tools/rte_library/build.py into
resources/rte_library/; a stand-in lives in test/fixtures/rte_library/):

    consensus.fa              per-class consensus (L1HS, L1PA2, ALU_Y, ALU_YA5, SVA_E, ...)
    consensus_landmarks.tsv   consensus, feature, start, end   (1-based inclusive; optional)
    l1_intact.fa/.tsv, alu_y_intact.fa/.tsv, sva_intact.fa/.tsv   intact elements, sense
    active.tsv                id, element_class|class, source_id, ...   (active/hot subset)
    transduction_sources.tsv  id, element_class, hs1_*, strand, l1base_id, flank_3p, ...
                              (+ optional intact_id; l1base_id / active.source_id link a source
                              to its intact element)
    flanks_3p.fa(.gz)         0-15 kb downstream of each source, element sense, repeats soft-masked;
                              strandless sources ship two records `<id>/+` and `<id>/-`
    flanks_5p_sva.fa(.gz)     upstream flanks of SVA sources (resources/: bgzip + .fai/.gzi)

Every TSV is read by header name; unknown columns are kept, missing optional files are empty.
The class of a consensus / intact element is derived from its name (L1*, ALU*/Alu*, SVA*)
unless a `class` column says otherwise. Which consensus counts as YOUNG (active subfamily) is a
regex (config `young_consensus_regex`), because consensus.fa also carries old subfamilies used as
negative controls (AluSx/AluJb/L1PA7...): an insertion that only matches those is
"inactive-only" (an artefact indicator in score.py).
"""
from __future__ import annotations

import csv
import gzip
import os
import re

from .sequtil import trailing_polya_start

DEFAULT_YOUNG = r"^(L1HS|L1PA[23](?![0-9])|ALU_?Y|SVA)"


def element_class(name: str) -> str:
    u = name.upper()
    if u.startswith("L1") or u.startswith("LINE"):
        return "L1"
    if u.startswith("ALU"):
        return "ALU"
    if u.startswith("SVA"):
        return "SVA"
    return "OTHER"


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def resolve_library_path(root):
    """CONFIG['annotate']['rte_library'] may be relative (e.g. 'resources/rte_library'): resolve
    it against the current directory if it exists there, else against the repository root."""
    if not root or os.path.isabs(str(root)) or os.path.isdir(str(root)):
        return root
    cand = os.path.join(REPO_ROOT, str(root))
    return cand if os.path.isdir(cand) else root


def _find(root, name):
    for cand in (name, name + ".gz"):
        p = os.path.join(root, cand)
        if os.path.exists(p):
            return p
    return None


def _read_fa(path):
    if path is None:
        return {}
    import mappy
    return {n: s for n, s, _ in mappy.fastx_read(path, read_comment=False)}


def _read_tsv(path):
    if path is None:
        return []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("##")]
    if not lines:
        return []
    lines[0] = lines[0].lstrip("#")
    return list(csv.DictReader(lines, delimiter="\t"))


class RteLibrary:
    # mappy parameters per index. Reads are 100-150 bp and element pieces at a junction can be
    # 30 bp, so a short k-mer / small window is used (short pieces missed by mappy are rescued by
    # an edlib pass in assembly.py).
    MAP_OPTS = dict(k=11, w=3, min_chain_score=18, min_dp_score=25, best_n=8)

    def __init__(self, root: str, cfg: dict | None = None):
        cfg = cfg or {}
        root = resolve_library_path(root)
        if not root or not os.path.isdir(root):
            raise FileNotFoundError(f"RTE library directory {root!r} not found")
        self.root = root
        self.paths = {k: _find(root, f) for k, f in (
            ("consensus", "consensus.fa"), ("landmarks", "consensus_landmarks.tsv"),
            ("l1", "l1_intact.fa"), ("alu", "alu_y_intact.fa"), ("sva", "sva_intact.fa"),
            ("l1_tsv", "l1_intact.tsv"), ("alu_tsv", "alu_y_intact.tsv"), ("sva_tsv", "sva_intact.tsv"),
            ("active", "active.tsv"), ("sources", "transduction_sources.tsv"),
            ("flanks3", "flanks_3p.fa"), ("flanks5", "flanks_5p_sva.fa"))}
        if self.paths["consensus"] is None:
            raise FileNotFoundError(f"{root}/consensus.fa missing")
        self.consensus = {n: s.upper() for n, s in _read_fa(self.paths["consensus"]).items()}
        self.cons_class = {n: element_class(n) for n in self.consensus}
        # element 3' end on the consensus = start of its trailing poly-A (Dfam/intact-derived
        # consensus sequences carry an A-tail; it must not be counted as element sequence)
        self.cons_end = {n: trailing_polya_start(s) for n, s in self.consensus.items()}
        young = re.compile(cfg.get("young_consensus_regex", DEFAULT_YOUNG), re.I)
        self.young = {n for n in self.consensus if young.search(n)}
        # consensus_landmarks.tsv is 1-based INCLUSIVE (resources/rte_library/README.md); stored
        # here 0-based half-open. Feature names are kept as written; compare case-insensitively
        # (the real library writes HEXAMER, older fixtures hexamer).
        self.landmarks = {}
        for r in _read_tsv(self.paths["landmarks"]):
            try:
                self.landmarks.setdefault(r["consensus"], []).append(
                    (r["feature"], int(r["start"]) - 1, int(r["end"])))
            except (KeyError, ValueError):
                continue
        self.intact = {}
        self.intact_class = {}
        self.intact_meta = {}
        for fa, tsv, cls in (("l1", "l1_tsv", "L1"), ("alu", "alu_tsv", "ALU"), ("sva", "sva_tsv", "SVA")):
            for n, s in _read_fa(self.paths[fa]).items():
                self.intact[n] = s.upper()
                self.intact_class[n] = cls
            for r in _read_tsv(self.paths[tsv]):
                rid = r.get("id") or next(iter(r.values()), None)
                if rid:
                    self.intact_meta[rid] = r
        self.active = {}
        # source id -> intact element id (active.tsv `source_id`, real library)
        self._source_to_intact = {}
        for r in _read_tsv(self.paths["active"]):
            rid = r.get("id") or next(iter(r.values()), None)
            if rid:
                self.active[rid] = r
                sid = r.get("source_id")
                if sid and sid not in (".", "") and rid in self.intact:
                    self._source_to_intact[sid] = rid
        # classes with at least one active.tsv entry. The real active.tsv lists L1 only; the
        # Alu/SVA intact sets are by construction the youngest copies of the active subfamilies,
        # so for a class without active.tsv rows every intact element of that class counts.
        self.active_classes = {self.intact_class.get(k) or element_class(
            (r.get("element_class") or r.get("class") or k)) for k, r in self.active.items()}
        self.sources = {}
        for r in _read_tsv(self.paths["sources"]):
            rid = r.get("id") or next(iter(r.values()), None)
            if rid:
                self.sources[rid] = r
                for k in ("l1base_id", "intact_id", "element_id"):
                    v = r.get(k)
                    if v and v not in (".", "") and v in self.intact:
                        self._source_to_intact.setdefault(rid, v)
        # Tubio 2014 S7 polymorphic L1 positions on hs1 (novel-source tier B only; optional)
        self.polymorphic_l1 = []
        for r in _read_tsv(_find(root, "polymorphic_l1_candidates.tsv")):
            try:
                c, pos = r["hs1_chrom"], int(r["hs1_pos"])
            except (KeyError, ValueError, TypeError):
                continue
            if c not in (".", "") and pos >= 0:
                self.polymorphic_l1.append((c, pos, r.get("id", ".")))
        self.flanks3 = _read_fa(self.paths["flanks3"])     # case kept: lower = soft-masked
        self.flanks5 = _read_fa(self.paths["flanks5"])
        self._aligners = {}

    # ------------------------------------------------------------------ helpers
    def is_young(self, cons_name: str) -> bool:
        return cons_name in self.young

    def is_active_element(self, intact_id: str) -> bool:
        if intact_id in self.active:
            return True
        cls_ = self.intact_class.get(intact_id)
        return cls_ is not None and cls_ not in self.active_classes

    @staticmethod
    def _split_flank_name(flank_name: str):
        base = re.split(r"[|\s]", flank_name)[0]
        m = re.match(r"^(.*)/([+-])$", base)
        if m:                      # strandless source: both candidate flanks `<id>/+`, `<id>/-`
            return m.group(1), m.group(2)
        return base, ""

    def source_for_flank(self, flank_name: str) -> str:
        """Flank record -> source id. Records are named by source id, optionally with a suffix
        after '|' or whitespace; sources without a published strand ship two flanks named
        `<id>/+` and `<id>/-` (transduction_sources.tsv `flank_3p` lists both)."""
        if flank_name in self.sources:
            return flank_name
        return self._split_flank_name(flank_name)[0]

    def flank_strand(self, flank_name: str) -> str:
        """'+'/'-' for a `<id>/<strand>` flank of a strandless source, '' otherwise."""
        return self._split_flank_name(flank_name)[1]

    def source_element(self, source_id: str) -> str | None:
        """Intact element id of a transduction source: column intact_id / element_id /
        l1base_id (real library: the L1Base UID), or active.tsv `source_id`, else the source id
        itself when it is also an intact element id."""
        r = self.sources.get(source_id, {})
        for k in ("intact_id", "element_id", "element"):
            if r.get(k) and r[k] not in (".", ""):
                return r[k]
        if source_id in self._source_to_intact:
            return self._source_to_intact[source_id]
        return source_id if source_id in self.intact else None

    def landmark_at(self, cons_name: str, pos: int) -> str:
        for f, s, e in self.landmarks.get(cons_name, []):
            if s <= pos < e:
                return f
        return "."

    # ------------------------------------------------------------------ aligners
    def aligner(self, kind: str):
        """Lazy mappy index: 'consensus' | 'intact' | 'flanks3' | 'flanks5'. None if empty."""
        if kind in self._aligners:
            return self._aligners[kind]
        import mappy
        al = None
        if kind == "consensus":
            al = mappy.Aligner(self.paths["consensus"], **self.MAP_OPTS)
        elif kind == "intact":
            paths = [self.paths[k] for k in ("l1", "alu", "sva") if self.paths[k]]
            if paths:
                al = _multi_fasta_aligner(paths, self.MAP_OPTS)
        elif kind in ("flanks3", "flanks5") and self.paths[kind]:
            al = mappy.Aligner(self.paths[kind], **self.MAP_OPTS)
        if al is not None and not al:
            al = None
        self._aligners[kind] = al
        return al


def _multi_fasta_aligner(paths, opts):
    """mappy can index one FASTA path; concatenate several into a temp file when needed."""
    import mappy
    import tempfile
    if len(paths) == 1:
        return mappy.Aligner(paths[0], **opts)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    with tmp:
        for p in paths:
            for n, s, _ in mappy.fastx_read(p, read_comment=False):
                tmp.write(f">{n}\n{s}\n")
    al = mappy.Aligner(tmp.name, **opts)
    os.unlink(tmp.name)
    return al
