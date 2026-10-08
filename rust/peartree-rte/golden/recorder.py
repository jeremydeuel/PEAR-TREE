"""Record tools/rte module-boundary calls as JSON lines (golden data for the Rust port).

`install(path, case_fn)` monkeypatches the python modules IN PLACE (the defining module and every
module that imported the name), so any code path -- the pytest suite, annotate_v2, a script --
emits one JSON event per boundary call while behaving exactly as before (wrappers return the
original result untouched). Event schema: rust/peartree-rte/SPEC.md "Golden harness"; the Rust
side parses it in src/golden.rs.

Recorded kinds
  genome            genome definitions referenced by id from other events
  hallmark.<fn>     split_junction, parse_locus, edge_run, polya_info, locate_site, target_site,
                    tsd_from_flanks, en_motif, slippage_context, foldback, en_bin
  assemble          Assembler.assemble: ctx, junction_seqs, reads, strand_hint, cfg -> AssemblyResult
  classify          structure.classify: assembly (or `assembly_ref` = seq of the assemble event that
                    produced it), ctx, cfg, legacy_class, pseudogene genes/hits + every callback
                    answer (novel_answers, premrna_answers, pg_structure_answers) -> StructureCall
  known_source      transduction.known_source
  novel_find        NovelSourceFinder.find (+ the locator's answers, cohort_l1, rmsk path)
  cons_identity     transduction.cons_identity
  exon_find / exon_structure   ExonJunctionIndex.find / .structure (+ the candidate genes' exons)
  premrna           the annotator's pre-mRNA callback (seq, site) -> label
  score             score.score
  annotate          RteAnnotator.annotate (InsertionInput, evidence, capped read names) -> record
  annotate_key      RteAnnotator._annotate_key -> record (gt_reads / gt_changed)
  final             annotate_all's records after the cohort + recurrence passes, with row()
"""
from __future__ import annotations

import functools
import gzip
import itertools
import json
import os

_STATE = {"libdir": None, "fh": None, "seq": itertools.count(1), "case": lambda: "", "stack": [], "genomes": {},
          "asm_ids": {}, "full_reads": True, "depth": 0}


# ----------------------------------------------------------------------------- serializers
def _jsonable(v):
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items() if not callable(x)}
    if isinstance(v, set):
        return sorted(_jsonable(x) for x in v)
    return repr(v)


def cfg_json(cfg):
    out = {k: _jsonable(v) for k, v in (cfg or {}).items() if not callable(v)}
    if isinstance(out.get("rte_library"), str):
        out["rte_library"] = _snapshot_lib(out["rte_library"])
    return out


def seg_json(s):
    return {"q": [s.q_st, s.q_en], "kind": s.kind, "target": s.target, "t": [s.t_st, s.t_en],
            "strand": s.strand, "identity": s.identity, "matches": s.matches}


def layout_json(l):
    return {"name": l.name, "side": l.side, "role": l.role, "frag": [str(x) for x in l.frag_key],
            "seq": l.seq, "sense": l.sense, "segs": [seg_json(s) for s in l.segments]}


def ctx_json(c):
    if c is None:
        return None
    return {"title": c.title, "contig": c.contig, "left_bp": c.left_bp, "right_bp": c.right_bp,
            "window_start": c.window_start, "window_seq": c.window_seq, "left_flank": c.left_flank,
            "right_flank": c.right_flank, "wide_start": c.wide_start, "wide_seq": c.wide_seq}


def asm_json(a):
    return {"strand": a.strand, "strand_source": a.strand_source, "element_class": a.element_class,
            "consensus": a.consensus, "element_bp": a.element_bp,
            "layouts": [layout_json(l) for l in a.layouts],
            "raw_layouts": [layout_json(l) for l in a.raw_layouts],
            "covered": [list(x) for x in a.covered], "covered_seqs": list(a.covered_seqs),
            "segments_on_cons": [[s, e, bool(sense), li] for s, e, sense, li in a.segments_on_cons],
            "consensus_identity": a.consensus_identity, "nearest_intact": a.nearest_intact,
            "nearest_intact_identity": a.nearest_intact_identity, "nearest_active": a.nearest_active,
            "element_identity": a.element_identity, "class_bp": dict(a.class_bp)}


def src_json(s):
    if s is None:
        return None
    return {"source_id": s.source_id, "td_end": s.td_end, "td_start": s.td_start,
            "n_segments": s.n_segments, "novel": s.novel, "identity": s.identity, "detail": s.detail,
            "tier": s.tier}


def call_json(c):
    return {"element": c.element, "structure": c.structure, "tags": list(c.tags),
            "detail": _jsonable(c.detail), "source": src_json(c.source), "j5_class": c.j5_class,
            "j3_class": c.j3_class, "j5_pos": c.j5_pos, "inv_p1": c.inv_p1,
            "three_prime_truncated": c.three_prime_truncated,
            "three_prime_short": c.three_prime_short, "has_polya_3p": c.has_polya_3p,
            "td3p_seq": c.td3p_seq}


def si_json(si):
    if si is None:
        return None
    d = dict(vars(si))
    d["tags"] = list(si.tags)
    return _jsonable(d)


def site_json(si):
    return _jsonable(dict(vars(si)))


def pa_json(p):
    return {"strand": p.strand, "source": p.source, "left_run": list(p.left_run),
            "right_run": list(p.right_run), "both_sided": p.both_sided, "length": p.length}


def rec_json(r):
    return {"insertion_id": r.insertion_id, "element": r.element, "structure": r.structure,
            "tags": list(r.tags), "covered_5p": r.covered_5p, "covered_3p": r.covered_3p,
            "covered_intervals": [list(x) for x in r.covered_intervals], "consensus": r.consensus,
            "element_identity": r.element_identity, "nearest_active": r.nearest_active,
            "tsd_seq": r.tsd_seq, "tsd_len": r.tsd_len, "en_motif": r.en_motif,
            "en_mismatches": r.en_mismatches, "polya_len": r.polya_len,
            "beyond_polya": r.beyond_polya, "beyond_polya_support": r.beyond_polya_support,
            "strand": r.strand, "tprt_score": r.tprt_score, "tprt_points": r.tprt_points,
            "tprt_call": r.tprt_call, "detail": _jsonable(r.detail),
            "score_input": si_json(r.score_input), "site": list(r.site),
            "gt_reads": r.gt_reads, "gt_changed": r.gt_changed, "row": r.row()}


def je_json(j):
    return _jsonable(dict(vars(j)))


def read_json(r):
    return [r.side, r.role, r.sample, r.frag, r.r12, r.seq]


def inp_json(i):
    return {"locus": i.title, "left_seq": i.left_seq, "right_seq": i.right_seq,
            "pseudogene_genes": list(i.pseudogene_genes or []), "legacy_class": i.legacy_class,
            "sv": _jsonable(i.sv)}


# ----------------------------------------------------------------------------- plumbing
def emit(kind, **kw):
    fh = _STATE["fh"]
    if fh is None:
        return None
    seq = next(_STATE["seq"])
    ev = {"kind": kind, "seq": seq, "case": _STATE["case"]()}
    if _STATE["stack"]:
        ev["parent"] = _STATE["stack"][-1][0]
        ev["locus"] = _STATE["stack"][-1][1]
    ev.update(kw)
    fh.write(json.dumps(ev, separators=(",", ":")) + "\n")
    return seq


def genome_id(g):
    """Id of a genome object (emitting its definition the first time)."""
    if g is None:
        return None
    key = id(g)
    if key in _STATE["genomes"]:
        return _STATE["genomes"][key][0]
    gid = f"g{len(_STATE['genomes']) + 1}"
    _STATE["genomes"][key] = (gid, g)      # keep a reference: ids must not be recycled
    paths = getattr(g, "_golden_paths", None)
    if hasattr(g, "path") and str(getattr(g, "path", "")).endswith(".2bit"):
        emit("genome", id=gid, path=rel(g.path))
    elif paths:
        emit("genome", id=gid, path=rel(paths[0]))
    else:
        emit("genome", id=gid, regions=[[c, off, seq] for c, off, seq in g.regions()])
    return gid


REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))


def rel(p):
    """A path inside the repository is recorded relative to its root (src/golden.rs resolves it
    against the checkout running the tests); anything else stays absolute."""
    p = os.path.abspath(str(p))
    return os.path.relpath(p, REPO) if p.startswith(REPO + os.sep) else p


def _snapshot_lib(root):
    """A library directory outside the repository (a pytest tmp_path stand-in) is copied next
    to the golden file (`libs/<n>_<name>`) so the events stay replayable."""
    root = os.path.abspath(str(root))
    if root.startswith(REPO + os.sep) or _STATE["libdir"] is None:
        return rel(root)
    snaps = _STATE.setdefault("snaps", {})
    if root not in snaps:
        import shutil
        dst = os.path.join(_STATE["libdir"], f"{len(snaps) + 1}_{os.path.basename(root)}")
        if os.path.exists(dst):
            shutil.rmtree(dst)
        shutil.copytree(root, dst)
        snaps[root] = rel(dst)
    return snaps[root]


def _lib_root(lib):
    """Library directory; a test stand-in object without `root` -> {"mock": {...}} with the
    attributes the callee reads (polymorphic_l1; consensus = the fixture library's)."""
    if lib is None:
        return None
    root = getattr(lib, "root", None)
    if root is None:
        return {"mock": {"polymorphic_l1": _jsonable(getattr(lib, "polymorphic_l1", None)),
                         "consensus": sorted(getattr(lib, "consensus", {}) or {})}}
    return _snapshot_lib(root)


def _patch(module, name, wrapper_factory, also=()):
    orig = getattr(module, name)
    w = wrapper_factory(orig)
    setattr(module, name, w)
    for m in also:
        if getattr(m, name, None) is orig:
            setattr(m, name, w)
    return orig


def install(path, case_fn=lambda: os.environ.get("PYTEST_CURRENT_TEST", "").rsplit(" (", 1)[0],
            full_reads=True):
    """Start recording into `path` (.jsonl or .jsonl.gz)."""
    from tools.rte import annotator as A, assembly as AS, hallmarks as H, score as SC, \
        structure as ST, transduction as TD, pseudogene as PG, genome as GN
    _STATE["fh"] = gzip.open(path, "wt") if str(path).endswith(".gz") else open(path, "w")
    _STATE["libdir"] = os.path.join(os.path.dirname(os.path.abspath(str(path))), "libs")
    _STATE["case"] = case_fn
    _STATE["full_reads"] = full_reads

    # FastaGenome built from files: remember the paths (record a path, not 50 MB of regions)
    orig_fg_init = GN.FastaGenome.__init__

    def fg_init(self, *paths, records=None):
        orig_fg_init(self, *paths, records=records)
        self._golden_paths = list(paths) if paths and not records else None
    GN.FastaGenome.__init__ = fg_init

    # ---- hallmarks (pure functions; genome args by id; in-place SiteInfo before/after)
    def pure(name, genome_arg=None, mutates=False, enc=None):
        def factory(orig):
            @functools.wraps(orig)
            def w(*a, **kw):
                args = list(a)
                rec_args = []
                for k, x in enumerate(args):
                    if genome_arg is not None and k == genome_arg:
                        rec_args.append({"genome": genome_id(x)})
                    elif isinstance(x, H.SiteInfo):
                        rec_args.append({"site": site_json(x)})
                    elif hasattr(x, "polya_len_median"):
                        rec_args.append({"junction": je_json(x)})
                    else:
                        rec_args.append(_jsonable(x))
                rkw = {k: ({"genome": genome_id(v)} if k == "genome" else
                           {"junction": je_json(v)} if hasattr(v, "polya_len_median") else _jsonable(v))
                       for k, v in kw.items()}
                out = orig(*a, **kw)
                if isinstance(out, H.SiteInfo):
                    o = site_json(out)
                elif isinstance(out, H.PolyAInfo):
                    o = pa_json(out)
                else:
                    o = _jsonable(out)
                emit(f"hallmark.{name}", args=rec_args, kwargs=rkw, out=o)
                return out
            return w
        _patch(H, name, factory, also=(A,))

    pure("split_junction")
    pure("parse_locus")
    pure("edge_run")
    pure("polya_info")
    pure("locate_site", genome_arg=3)
    pure("tsd_from_flanks")
    pure("target_site", genome_arg=3)
    pure("en_motif", genome_arg=2)
    pure("en_bin")
    pure("slippage_context", genome_arg=2)
    pure("foldback")

    # ---- score
    def score_factory(orig):
        @functools.wraps(orig)
        def w(si, weights=None, thresholds=None):
            out = orig(si, weights, thresholds)
            emit("score", **{"in": {"score_input": si_json(si), "weights": _jsonable(weights),
                                    "thresholds": _jsonable(thresholds)}, "out": list(out)})
            return out
        return w
    _patch(SC, "score", score_factory, also=(A,))

    # ---- assembly
    orig_assemble = AS.Assembler.assemble

    def assemble(self, ctx, junction_seqs, reads, strand_hint=None):
        cj = ctx_json(ctx)
        js = [[k, v[0], list(v[1])] for k, v in junction_seqs.items()]
        rd = [read_json(r) for r in reads]
        res = orig_assemble(self, ctx, junction_seqs, reads, strand_hint)
        seq = emit("assemble", lib=_lib_root(self.lib), cfg=cfg_json(self.cfg),
                   **{"in": {"ctx": cj, "junction_seqs": js, "reads": rd,
                             "strand_hint": _jsonable(strand_hint)}, "out": asm_json(res)})
        _STATE["asm_ids"][id(res)] = (seq, res)
        return res
    AS.Assembler.assemble = assemble

    # ---- structure.classify (+ callbacks)
    class _FinderProxy:
        def __init__(self, inner, log):
            self._inner, self._log = inner, log

        def find(self, seq):
            out = self._inner.find(seq)
            self._log.append([seq, src_json(out)])
            return out

        def __getattr__(self, k):
            return getattr(self._inner, k)

    def classify_factory(orig):
        @functools.wraps(orig)
        def w(res, lib, ctx=None, cfg=None, novel_finder=None, premrna=None, pseudogene=None,
              legacy_class=None):
            novel_log, premrna_log, pgs_log = [], [], []
            nf = _FinderProxy(novel_finder, novel_log) if novel_finder is not None else None
            pm = None
            if premrna is not None:
                def pm(seq, _f=premrna):
                    out = _f(seq)
                    premrna_log.append([seq, out])
                    return out
            pg = pseudogene
            pg_in = None
            if pseudogene:
                pg = list(pseudogene)
                pg_in = {"genes": list(pg[0]), "hits": [list(h) for h in pg[1]]}
                if len(pg) > 2 and pg[2] is not None:
                    def pgs(layouts, _f=pg[2]):
                        out = _f(layouts)
                        pgs_log.append([[layout_json(l) for l in layouts], out])
                        return out
                    pg[2] = pgs
                    pg_in["has_structure_fn"] = True
                pg = tuple(pg)
            ref = _STATE["asm_ids"].get(id(res))
            asm_in = {"assembly_ref": ref[0]} if ref is not None and ref[1] is res else {"assembly": asm_json(res)}
            out = orig(res, lib, ctx, cfg, nf, pm, pg, legacy_class)
            emit("classify", lib=_lib_root(lib), cfg=cfg_json(cfg), **{"in": dict(
                asm_in, ctx=ctx_json(ctx), legacy_class=legacy_class, pseudogene=pg_in,
                has_novel=novel_finder is not None, has_premrna=premrna is not None,
                novel_answers=novel_log, premrna_answers=premrna_log,
                pg_structure_answers=pgs_log), "out": call_json(out)})
            return out
        return w
    _patch(ST, "classify", classify_factory, also=(A,))

    # ---- transduction
    def known_factory(orig):
        @functools.wraps(orig)
        def w(flank_segments, lib):
            out = orig(flank_segments, lib)
            emit("known_source", lib=_lib_root(lib),
                 **{"in": {"segments": [seg_json(s) for s in flank_segments]}, "out": src_json(out)})
            return out
        return w
    _patch(TD, "known_source", known_factory, also=(ST,))

    def ci_factory(orig):
        @functools.wraps(orig)
        def w(seq, cons):
            out = orig(seq, cons)
            emit("cons_identity", **{"in": {"seq": seq, "cons": cons}, "out": out})
            return out
        return w
    _patch(TD, "cons_identity", ci_factory)

    orig_find = TD.NovelSourceFinder.find

    def find(self, seq):
        loc_log = []
        loc = self.locator
        if loc is not None:
            def logged(s, _l=loc):
                out = _l(s)
                loc_log.append([s, _jsonable(out)])
                return out
            self.locator = logged
        try:
            out = orig_find(self, seq)
        finally:
            self.locator = loc
        rm = self.rmsk
        emit("novel_find", lib=_lib_root(self.lib), remap=genome_id(self.genome), cfg=cfg_json(self.cfg),
             **{"in": {"seq": seq, "available": self.available(), "locator_answers": loc_log,
                       "cohort_l1": _jsonable(self.cohort_l1),
                       "rmsk": getattr(rm, "_golden_path", None) if rm is not None else None},
                "out": src_json(out)})
        return out
    TD.NovelSourceFinder.find = find

    orig_rmsk_init = TD.L1Rmsk.__init__

    def rmsk_init(self, path, min_len=5500):
        orig_rmsk_init(self, path, min_len)
        self._golden_path = rel(path)
    TD.L1Rmsk.__init__ = rmsk_init

    # ---- pseudogene
    def _pg_ctx(self, genes):
        return {"exons": {g: [list(x) for x in self.exons.get(g, [])] for g in genes},
                "strands": {g: self.strands[g] for g in genes if g in self.strands}}

    orig_pg_find = PG.ExonJunctionIndex.find

    def pg_find(self, genes, seqs):
        out = orig_pg_find(self, genes, seqs)
        emit("exon_find", remap=genome_id(self.genome), cfg=cfg_json(self.cfg),
             **{"in": dict(_pg_ctx(self, genes), genes=list(genes), seqs=[list(s) for s in seqs]),
                "out": [list(h) for h in out]})
        return out
    PG.ExonJunctionIndex.find = pg_find

    orig_pg_structure = PG.ExonJunctionIndex.structure

    def pg_structure(self, genes, five_prime_seqs, tol=15, probe=25):
        out = orig_pg_structure(self, genes, five_prime_seqs, tol, probe)
        emit("exon_structure", remap=genome_id(self.genome), cfg=cfg_json(self.cfg),
             **{"in": dict(_pg_ctx(self, genes), genes=list(genes), seqs=list(five_prime_seqs),
                           tol=tol, probe=probe), "out": out})
        return out
    PG.ExonJunctionIndex.structure = pg_structure

    # ---- annotator
    orig_premrna_fn = A.RteAnnotator._premrna_fn

    def premrna_fn(self, site):
        fn = orig_premrna_fn(self, site)
        if fn is None:
            return None

        gm = self.gene_model

        class _GmProxy:
            """logs every _candidates / _genic_feature answer (works for duck-typed models)"""
            def __init__(self, inner, log):
                self._inner, self._log = inner, log

            def _candidates(self, contig, p):
                out = self._inner._candidates(contig, p)
                self._log.append(["candidates", [contig, p], _jsonable(out)])
                return out

            def _genic_feature(self, exons, strand, p):
                out = self._inner._genic_feature(exons, strand, p)
                self._log.append(["genic_feature", [_jsonable(exons), strand, p], _jsonable(out)])
                return out

            def __getattr__(self, k):
                return getattr(self._inner, k)

        def w(seq, _self=self):
            log = []
            _self.gene_model = _GmProxy(gm, log) if gm is not None else None
            try:
                out = fn(seq)
            finally:
                _self.gene_model = gm
            model = None
            if gm is not None and hasattr(gm, "genes") and hasattr(gm, "_resolve"):
                key = gm._resolve(site.contig)
                model = {"contig": key, "genes": [[gs, ge, name, strand, [list(x) for x in exons]]
                                                  for gs, ge, name, strand, exons in gm.genes.get(key, [])],
                         "donor_window": gm.donor_window, "acceptor_window": gm.acceptor_window,
                         "ppt_window": gm.ppt_window, "branch_window": gm.branch_window,
                         "prom_up": gm.prom_up}
            emit("premrna", genome=genome_id(self.genome),
                 **{"in": {"seq": seq, "site": [site.contig, site.L, site.R],
                           "window": int(self.cfg.get("rte_premrna_window", 1_000_000)),
                           "gene_model": model, "gene_model_calls": log}, "out": out})
            return out
        return w
    A.RteAnnotator._premrna_fn = premrna_fn

    orig_annotate = A.RteAnnotator.annotate

    def annotate(self, inp, ev=None):
        seq = emit("annotate_begin", **{"in": {"input": inp_json(inp)}})
        _STATE["stack"].append((seq, inp.title))
        try:
            out = orig_annotate(self, inp, ev)
        finally:
            _STATE["stack"].pop()
        evj = None
        if ev is not None:
            evj = {"junctions": [je_json(j) for j in ev.junctions.values()], "n_reads": len(ev.reads)}
            if _STATE["full_reads"]:
                evj["reads"] = [read_json(r) for r in ev.reads]
            capped = self._cap_reads(ev.reads, int(self.cfg.get("rte_max_reads", 400)))
            evj["capped"] = [f"{r.side}|{r.role}|{r.sample}|{r.frag}|{r.r12}" for r in capped]
        emit("annotate", lib=_lib_root(self.lib), genome=genome_id(self.genome),
             remap=genome_id(self.remap), cfg=cfg_json(self.cfg), begin=seq,
             **{"in": {"input": inp_json(inp), "evidence": evj}, "out": rec_json(out)})
        return out
    A.RteAnnotator.annotate = annotate

    orig_key = A.RteAnnotator._annotate_key

    def annotate_key(self, key, inp):
        out = orig_key(self, key, inp)
        emit("annotate_key", **{"in": {"key": key, "n_gt": len(self.gt_reads.get(key) or [])},
                                "out": rec_json(out)})
        return out
    A.RteAnnotator._annotate_key = annotate_key

    orig_all = A.RteAnnotator.annotate_all

    def annotate_all(self, inputs):
        out = orig_all(self, inputs)
        for k, rec in out.items():
            emit("final", has_gt_reads=self.has_gt_reads,
                 **{"in": {"key": k}, "out": rec_json(rec)})
        return out
    A.RteAnnotator.annotate_all = annotate_all


def close():
    if _STATE["fh"] is not None:
        _STATE["fh"].close()
        _STATE["fh"] = None
