# PEAR-TREE - tools/rte TPRT point system + annotate_v2 wiring.
#
#  * score(): every feature is reported as `name:+n`, the total is thresholded, hard artefact
#    signatures cap the call; weights/thresholds are overridable.
#  * on simulated insertions, true TPRT events outscore the artefact classes (chimeric ends,
#    poly-A slippage, 60 bp TSD chimera, fold-back) -- the separation calibrate.py measures.
#  * recurrence pass, and the annotate_v2 integration: one-sided combined records no longer crash
#    read_insertions(); with `rte_library` configured write_table() appends the SPEC columns.
#
# Run:  pytest test/test_rte_score.py
import gzip
import importlib.util
import os
import sys
import tempfile
import types

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rte_sim  # noqa: E402
from tools.rte.annotator import RteAnnotator  # noqa: E402
from tools.rte.record import RteRecord  # noqa: E402
from tools.rte.score import ScoreInput, score, WEIGHTS  # noqa: E402
from tools.rte.sequtil import rc  # noqa: E402

LIB = rte_sim.library()
L1 = LIB.consensus["L1HS"]
CEND = LIB.cons_end["L1HS"]
A = "A" * 25


def run(insert, **kw):
    inp, ev, genome, info = rte_sim.build(insert, **kw)
    return RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome).annotate(inp, ev)


# ------------------------------------------------------------------------------------- unit
def _tp(**kw):
    d = dict(element="L1", structure="TRUNCATED_5P", tsd_len=15, tsd_verified=True, polya_len=25,
             beyond_polya_len=40, beyond_polya_support=3, en_mismatches=0, ends_concordant=True,
             element_identity=0.995, junctions_supported=2, n_samples=2)
    d.update(kw)
    return ScoreInput(**d)


def test_points_are_transparent_and_sum():
    total, points, call = score(_tp())
    parts = dict(p.split(":") for p in points.split(";"))
    assert abs(sum(float(v) for v in parts.values()) - total) < 1e-9
    assert set(parts) >= {"tsd_4_25", "polya_ge10", "beyond_polya", "beyond_polya_2frag", "en_0_1",
                          "ends_concordant", "active_identity_ge98", "junction_supported", "multi_colony"}
    assert call == "TPRT"


def test_en_bins_graded():
    s0 = score(_tp(en_mismatches=0))[0]
    s2 = score(_tp(en_mismatches=2))[0]
    s3 = score(_tp(en_mismatches=3))[0]
    s5 = score(_tp(en_mismatches=5))[0]
    assert s0 > s2 > s3 > s5


def test_beyond_polya_two_fragments_earn_more():
    one = score(_tp(beyond_polya_support=1))[0]
    two = score(_tp(beyond_polya_support=2))[0]
    none = score(_tp(beyond_polya_len=0))[0]
    assert two > one > none


@pytest.mark.parametrize("bad", [dict(tags=["CHIMERIC_ENDS"]), dict(tsd_len=60),
                                 dict(polya_both_sides=True)])
def test_hard_artefact_signatures_cap_the_call(bad):
    total, points, call = score(_tp(**bad))
    assert call in ("UNCERTAIN", "ARTEFACT_LIKE")


def test_slippage_voids_en_points():
    total, points, call = score(_tp(slippage=True, tsd_len=0, beyond_polya_len=0,
                                    ends_concordant=False, element="POLYA_ONLY",
                                    element_identity=0, junctions_supported=1, n_samples=1))
    assert "en_0_1" not in points and "slippage_context" in points
    assert call == "ARTEFACT_LIKE"


def test_twin_priming_590_rule():
    ok = score(_tp(structure="INVERTED_5P", inv_p1=3000))[1]
    bad = score(_tp(structure="INVERTED_5P", inv_p1=300))[1]
    assert "twin_priming_ge590" in ok and "twin_priming_lt590" in bad


def test_weights_and_thresholds_overridable():
    base = score(_tp())[0]
    t2, p2, c2 = score(_tp(), weights={"tsd_4_25": 10.0}, thresholds={"TPRT": 100})
    assert t2 == base + 7 and c2 == "LIKELY_TPRT"
    assert all(isinstance(v, float) for v in WEIGHTS.values())


# ------------------------------------------------------------------------------------- separation
def test_true_insertions_outscore_artefacts():
    tps = [run(L1[:CEND] + A), run(L1[5500:CEND] + A),
           run(rc(L1[CEND - 2014 - 1000:CEND - 2014]) + L1[CEND - 2000:CEND] + A),
           run(LIB.consensus["ALU_YA5"][:LIB.cons_end["ALU_YA5"]] + A, strand=-1, tsd=12),
           run("A" * 30, tsd=12)]
    arts = [run(LIB.consensus["ALU_YA5"][:200] + L1[5700:CEND] + A),       # chimeric ends
            run("A" * 14, tsd=0, ref_homopolymer=15),                     # poly-A slippage
            run(L1[5000:CEND] + A, tsd=60)]                               # TSD 60 bp chimera
    assert all(r.tprt_call == "TPRT" for r in tps), [(r.element, r.tprt_points) for r in tps]
    assert all(r.tprt_call not in ("TPRT", "LIKELY_TPRT") for r in arts), \
        [(r.element, r.tprt_points) for r in arts]
    assert min(r.tprt_score for r in tps) > max(r.tprt_score for r in arts)


def test_recurrence_pass():
    ann = RteAnnotator({"rte_library": rte_sim.FIX, "rte_recurrence_max": 2})
    recs = {}
    for i in range(4):
        r = RteRecord(f"chr1:{i}-{i + 15}", element="L1", structure="TRUNCATED_5P",
                      detail={"j5": 5500 + i % 2}, score_input=_tp())
        ann._score(r)
        recs[r.insertion_id] = r
    ann.recurrence_pass(recs)
    assert all("recurrence" in r.tprt_points for r in recs.values())
    ann2 = RteAnnotator({"rte_library": rte_sim.FIX, "rte_recurrence_max": 5})
    ann2.recurrence_pass(recs := {k: RteRecord(k, element="L1", structure="TRUNCATED_5P",
                                               detail={"j5": 5500}, score_input=_tp())
                                  for k in ("a", "b")})
    assert all(r.score_input.recurrent is False for r in recs.values())


# ------------------------------------------------------------------------------------- annotate_v2
def _load_annotate_v2():
    sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
    try:
        import pysam  # noqa: F401
    except Exception:
        sys.modules["pysam"] = types.ModuleType("pysam")
    if "src.config" not in sys.modules:
        try:
            import src.config  # noqa: F401
        except Exception:
            cfg = types.ModuleType("src.config")
            cfg.CONFIG = {"annotate": {}, "combine_insertions": {}}
            sys.modules["src.config"] = cfg
    spec = importlib.util.spec_from_file_location(
        "annotate_v2_rte", os.path.join(HERE, "..", "tools", "annotate_v2.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_annotate_v2_one_sided_records_and_rte_columns():
    av2 = _load_annotate_v2()
    tmp = tempfile.mkdtemp()
    a = rte_sim.build(L1[5500:CEND] + A, contig="chrA")
    b = rte_sim.build(LIB.consensus["ALU_YA5"][:LIB.cons_end["ALU_YA5"]] + A, strand=-1, tsd=12,
                      contig="chrB")
    combined = os.path.join(tmp, "P1.combined.txt.gz")
    with gzip.open(combined, "wt") as fh:
        for inp, *_ in (a, b):
            fh.write(rte_sim.fastq_record(inp.title, "L", inp.left_seq))
            fh.write(rte_sim.fastq_record(inp.title, "R", inp.right_seq))
        # one-sided record (poly-A anchored left end; only the RIGHT junction has reads)
        fh.write(rte_sim.fastq_record("chrC:polyA_5-900", "R", "ACGTACGTAAGGTTCC" + "ttttttttttttttt"))
    rte_sim.write_sidecars(os.path.join(tmp, "P1"), [a[1], b[1]])
    gfa = os.path.join(tmp, "genome.fa")
    with open(gfa, "w") as fh:
        for x in (a, b):
            fh.write(f">{x[0].title.split(':')[0]}\n{x[3]['ref']}\n")

    vac = object.__new__(av2.VariantAnnotationContainer)
    vac.sample, vac.insertions_file, vac.insertions = "P1", combined, {}
    vac.read_insertions()
    assert len(vac.insertions) == 3
    one = vac.insertions["chrC:polyA_5-900"]
    assert one.left_seq == "" and one.get_fasta()           # no crash on the empty side
    ann_cfg = av2.CONFIG["annotate"]
    saved = dict(ann_cfg)
    try:
        ann_cfg.update({"rte_library": rte_sim.FIX, "genome_2bit": gfa})
        vac.rte_records = vac.run_rte()
    finally:
        ann_cfg.clear()
        ann_cfg.update(saved)
    out = os.path.join(tmp, "P1.annotated.tsv")
    vac.write_table(out)
    rows = [l.rstrip("\n").split("\t") for l in open(out)]
    hdr = rows[0]
    assert hdr[-len(RteRecord.COLUMNS):] == RteRecord.COLUMNS
    by = {r[0]: dict(zip(hdr, r)) for r in rows[1:]}
    ra, rb = by[a[0].title], by[b[0].title]
    assert ra["element"] == "L1" and ra["structure"] == "TRUNCATED_5P" and ra["tsd_len"] == "15"
    assert rb["element"] == "ALU" and rb["tsd_len"] == "12" and rb["tprt_call"] == "TPRT"
    assert by["chrC:polyA_5-900"]["element"] in ("POLYA_ONLY", "UNKNOWN")
    # legacy columns untouched; without rte_library no extra columns
    vac.rte_records = {}
    vac.write_table(out)
    assert open(out).readline().rstrip("\n").split("\t")[-1] == "site_strand"


# ------------------------------------------------------------------------------------- calibrate
def test_calibrate_reports_separation_and_roc():
    from tools.rte.calibrate import calibrate, read_table, format_report
    cases = [  # (truth element, structure, tags, role, insert, build kwargs)
        ("L1", "TRUNCATED_5P", "", "TP", L1[5500:CEND] + A, {}),
        ("L1", "FULL_LENGTH", "", "TP", L1[:CEND] + A, {}),
        ("ALU", "FULL_LENGTH", "", "TP", LIB.consensus["ALU_YA5"][:LIB.cons_end["ALU_YA5"]] + A,
         {"strand": -1, "tsd": 12}),
        ("POLYA_ONLY", "5P_UNRESOLVED", "", "TP", "A" * 30, {"tsd": 12}),
        ("L1", "5P_UNRESOLVED", "CHIMERIC_ENDS", "ARTEFACT",
         LIB.consensus["ALU_YA5"][:200] + L1[5700:CEND] + A, {}),
        ("POLYA_ONLY", "5P_UNRESOLVED", "", "ARTEFACT", "A" * 14, {"tsd": 0, "ref_homopolymer": 15}),
        ("L1", "TRUNCATED_5P", "", "ARTEFACT", L1[5000:CEND] + A, {"tsd": 60}),
    ]
    tmp = tempfile.mkdtemp()
    tp, ap = os.path.join(tmp, "truth.tsv"), os.path.join(tmp, "annot.tsv")
    with open(tp, "w") as t, open(ap, "w") as a:
        t.write("insertion_id\telement\tstructure\ttags\trole\n")
        a.write("\t".join(["locus"] + RteRecord.COLUMNS) + "\n")
        for i, (el, st, tags, role, ins, kw) in enumerate(cases):
            r = run(ins, contig=f"chr{i}", **kw)
            t.write(f"{r.insertion_id}\t{el}\t{st}\t{tags}\t{role}\n")
            a.write("\t".join([r.insertion_id] + r.row()) + "\n")
    rep = calibrate(read_table(tp), read_table(ap))
    assert rep["n_matched"] == len(cases) and rep["n_tp"] == 4 and rep["n_artefact"] == 3
    assert rep["auc"] == 1.0
    assert rep["element_accuracy"] == 1.0
    assert rep["tags"]["CHIMERIC_ENDS"]["recall"] == 1.0
    assert rep["features"]["tsd_gt50"]["frac_ARTEFACT"] > 0 and rep["features"]["tsd_gt50"]["frac_TP"] == 0
    assert "ROC AUC" in format_report(rep)
