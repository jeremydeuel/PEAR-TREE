# PEAR-TREE - tools/rte "annotate round 2" regression tests (plans/tprt_hallmarks/E2E_REPORT.md,
# section "Annotate round 2"). One test per fix found on the E2E's combined TPs:
#   L1_INV_SWITCH (sense 5' piece -> inverted piece), SVA 5' transduction longer than a fragment,
#   local pre-mRNA without a gene model, fold-back of the 5' flank, poly-A noise is not a TD3P,
#   TD3P needs a complex >= 30 bp piece before the TAIL poly-A, pseudogene proof beats an Alu
#   piece inside the mRNA (+ 5' structure from the transcript), EN_INDEPENDENT with a blunt 1 bp
#   "TSD", orphan TD carries TD3P, annotate_v2 counts the genotyper's `insertion` call in TPRT
#   mode, build_gene_model reads the UCSC RefSeq GTF.
#
# Run:  pytest test/test_rte_round2.py
import gzip
import importlib.util
import os
import random
import subprocess
import sys
import types

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import rte_sim  # noqa: E402
from tools.rte.annotator import RteAnnotator  # noqa: E402
from tools.rte.assembly import Assembler, ReadLayout, Segment, SiteContext  # noqa: E402
from tools.rte.genome import FastaGenome  # noqa: E402
from tools.rte.pseudogene import ExonJunctionIndex  # noqa: E402
from tools.rte.sequtil import low_complexity, rc  # noqa: E402
from tools.rte.structure import DEFAULTS, _foldback_5p, _local_templates  # noqa: E402

LIB = rte_sim.library()
L1 = LIB.consensus["L1HS"]
CEND = LIB.cons_end["L1HS"]
A = "A" * 25


def annotate(insert, ann_kw=None, inp_kw=None, cfg=None, **kw):
    inp, ev, genome, info = rte_sim.build(insert, **kw)
    for k, v in (inp_kw or {}).items():
        setattr(inp, k, v)
    c = {"rte_library": rte_sim.FIX}
    c.update(cfg or {})
    ann = RteAnnotator(c, genome=genome, **(ann_kw or {}))
    return ann.annotate(inp, ev)


# ------------------------------------------------------------------------- structure
def test_inverted_switch_sense_piece_first():
    # catalogue 13 as simulated: L1_SWITCH (sense) | L1_INV (anti, further 3') | L1 | polyA
    ins = L1[3000:3250] + rc(L1[3600:4300]) + L1[4300:CEND] + A
    r = annotate(ins)
    assert r.element == "L1"
    assert r.structure == "INVERTED_5P_SWITCH"
    assert "switch" in r.detail


def test_plain_truncation_is_not_a_switch():
    r = annotate(L1[4000:CEND] + A)
    assert r.structure == "TRUNCATED_5P"
    assert "switch" not in r.detail


def test_sva_long_5p_transduction_is_full_length():
    f5 = list(rte_sim.fasta("flanks_5p_sva.fa").values())[0]
    s = LIB.consensus["SVA_F"][:LIB.cons_end["SVA_F"]]
    r = annotate(f5[-900:] + s + A)              # 5' TD longer than a fragment: no read chains it
    assert r.element == "SVA"
    assert r.structure == "FULL_LENGTH"
    assert "TD5P" in r.tags
    assert any(t.startswith("TD5P_SOURCE=") for t in r.tags)


def test_en_independent_blunt_1bp_and_short_3p_end():
    r = annotate(L1[3000:CEND - 40], tsd=1, motif=None)
    assert "EN_INDEPENDENT" in r.tags


# ------------------------------------------------------------------------- tags
def test_polya_noise_is_not_a_transduction():
    # junk between two poly-A runs (SBS jitter / low-quality tail) must not read as a TD3P
    r = annotate(L1[5600:CEND] + "A" * 30 + "GTCAGTTCAG" + "A" * 30)
    assert "TD3P" not in r.tags
    assert "TEMPLATED_LOCAL" not in r.tags


def test_unexplained_tail_before_polya_is_td3p():
    rng = random.Random(3)
    tag = "".join(rng.choice("ACGT") for _ in range(80))
    r = annotate(L1[5600:CEND] + tag + A)
    assert "TD3P" in r.tags
    assert r.detail.get("td_frags", 0) >= 2


def test_orphan_td_carries_td3p():
    fl = rte_sim.fasta("flanks_3p.fa")["TD_UID-5"]
    r = annotate(fl[500:900] + A)
    assert r.element == "ORPHAN_TD"
    assert "TD3P" in r.tags and "TD3P_SOURCE=TD_UID-5" in r.tags
    assert abs(r.detail["td_end"] - 900) <= 5


def test_local_premrna_without_gene_model():
    kw = dict(glen=40000, site=20000)
    _, _, _, info = rte_sim.build(L1[5000:CEND] + A, **kw)
    templ = info["ref"][21500:21800]              # 1.5 kb from the site: beyond templated range
    r = annotate(templ + L1[5000:CEND] + A, **kw)
    assert "PREMRNA_COINSERT" in r.tags
    assert r.detail["premrna"].startswith("local:")
    assert "TEMPLATED_LOCAL" not in r.tags
    assert "L1_MED_DUPLICATION" not in r.tags and "L1_MED_DELETION" not in r.tags


def test_foldback_of_the_5p_flank():
    _, _, _, info = rte_sim.build(L1[5000:CEND] + A)
    R = info["R"]
    fb = rc(info["ref"][R - 3 - 45:R - 3])        # inverted copy of the 45 bp right before the site
    r = annotate(fb + L1[5000:CEND] + A)
    assert "FOLDBACK_INVDUP_5P" in r.tags
    assert "TEMPLATED_LOCAL" not in r.tags


def test_templated_local_still_found():
    _, _, _, info = rte_sim.build(L1[5000:CEND] + A)
    L = info["L"]
    r = annotate(info["ref"][L + 120:L + 160] + L1[5000:CEND] + A)
    assert "TEMPLATED_LOCAL" in r.tags
    assert r.detail["templated_frags"] >= 2


# ------------------------------------------------------------------------- helpers (unit)
def test_low_complexity():
    assert low_complexity("AAAAAAAAAAAGAAAAAAAAAT")
    assert low_complexity("CACACACACACACACACACA")
    assert not low_complexity("ACGTTGCAGTCCATGACGTA")


def _ctx():
    return SiteContext("chrS:1000-1015", "chrS", 1000, 1015, window_start=400, window_seq="N" * 1200)


def test_local_template_needs_two_fragments_and_complexity():
    c = dict(DEFAULTS)
    seq = "ACGTTGCAGTCCATGACGTAGGCATTCAGG"
    one = ReadLayout("r1", "RIGHT", "CLIP", ("S1", "a"), seq,
                     [Segment(0, 30, "LOCAL", "site", 700, 730, 1, 1.0, 30)])
    assert _local_templates([one], _ctx(), c) == {}
    two = ReadLayout("r2", "RIGHT", "CLIP", ("S1", "b"), seq,
                     [Segment(0, 30, "LOCAL", "site", 702, 732, 1, 1.0, 30)])
    out = _local_templates([one, two], _ctx(), c)
    assert out["templated"][1] == 2 and out["templated"][0] <= 250
    polya = ReadLayout("r3", "LEFT", "CLIP", ("S1", "c"), "A" * 28 + "GA",
                       [Segment(0, 30, "LOCAL", "site", 700, 730, 1, 1.0, 30)])
    polyb = ReadLayout("r4", "LEFT", "CLIP", ("S1", "d"), "A" * 28 + "GA",
                       [Segment(0, 30, "LOCAL", "site", 700, 730, 1, 1.0, 30)])
    assert _local_templates([polya, polyb], _ctx(), c) == {}


def test_foldback_helper_geometry():
    c = dict(DEFAULTS)
    lays = []
    for k in "ab":
        lays.append(ReadLayout(k, "RIGHT", "CLIP", ("S1", k), "N" * 150,
                               [Segment(0, 100, "REF", "site", 500, 600, 1, 1.0, 100),
                                Segment(100, 130, "LOCAL", "site", 565, 595, -1, 1.0, 30)]))
    assert len(_foldback_5p(lays, c)) == 2
    # same strand as the flank = not a fold-back
    for lay in lays:
        lay.segments[1].strand = 1
    assert _foldback_5p(lays, c) == []


def test_polya_smoothing_merges_jittered_tail():
    asm = Assembler(LIB)
    lay = ReadLayout("r", "LEFT", "CLIP", ("S1", "a"), "A" * 30 + "GTCAGTTCAG" + "A" * 30 + "C" * 40,
                     [Segment(0, 30, "POLYA", "A", strand=1), Segment(30, 40, "UNKNOWN"),
                      Segment(40, 70, "POLYA", "A", strand=1),
                      Segment(70, 110, "REF", "site", 600, 640, 1, 1.0, 40)])
    asm._smooth_polya(lay)
    assert [s.kind for s in lay.segments] == ["POLYA", "REF"]
    assert lay.segments[0].q_en == 70


def test_split_flank_is_merged_not_templated():
    lay = ReadLayout("r", "RIGHT", "CLIP", ("S1", "a"), "N" * 151,
                     [Segment(0, 22, "REF", "site", 531, 553, 1, 1.0, 22),
                      Segment(22, 30, "POLYA", "A", strand=1),
                      Segment(30, 84, "REF", "site", 561, 615, 1, 1.0, 54),
                      Segment(84, 151, "ELEMENT", "L1HS", 5000, 5067, 1, 1.0, 67)])
    Assembler._mark_local(lay, None)
    assert [s.kind for s in lay.segments] == ["REF", "ELEMENT"]


# ------------------------------------------------------------------------- pseudogene
def _gene():
    rng = random.Random(11)
    g = "".join(rng.choice("ACGT") for _ in range(6000))
    exons = [(500, 700), (1500, 1650), (2800, 3100)]
    remap = FastaGenome(records={"chrG": g})
    return remap, {"PGENE": [("chrG", s, e) for s, e in exons]}, "".join(g[s:e] for s, e in exons)


def test_pseudogene_beats_an_alu_piece_inside_the_mrna():
    remap, exons, mrna = _gene()
    alu = LIB.consensus["ALU_YA5"][:LIB.cons_end["ALU_YA5"]]
    m = mrna[:400] + alu[100:250] + mrna[400:]    # an Alu in the 3' UTR-like part
    r = annotate(m + A, ann_kw={"remap_genome": remap, "exons_by_gene": exons},
                 inp_kw={"pseudogene_genes": ["PGENE"]})
    assert r.element == "PSEUDOGENE"
    assert "EXON_JUNCTION" in r.tags
    assert r.structure == "FULL_LENGTH"
    assert not any(t.startswith("TD3P") for t in r.tags)


def test_pseudogene_5p_truncated_structure():
    remap, exons, mrna = _gene()
    r = annotate(mrna[120:] + A, ann_kw={"remap_genome": remap, "exons_by_gene": exons},
                 inp_kw={"pseudogene_genes": ["PGENE"]})
    assert r.element == "PSEUDOGENE"
    assert r.structure == "TRUNCATED_5P"


def test_exon_index_structure_minus_strand():
    remap, exons, mrna = _gene()
    idx = ExonJunctionIndex(exons, remap, strands={"PGENE": "-"})
    tx = rc(mrna)
    assert idx.mrna("PGENE") == tx
    assert idx.structure(["PGENE"], [tx[:40]]) == "FULL_LENGTH"
    assert idx.structure(["PGENE"], [tx[200:240]]) == "TRUNCATED_5P"
    assert idx.structure(["PGENE"], ["ACGT" * 10]) is None


# ------------------------------------------------------------------------- annotate_v2 genotypes
def _annotate_v2():
    try:
        import pysam  # noqa: F401
    except Exception:
        sys.modules["pysam"] = types.ModuleType("pysam")
    try:
        import src.config  # noqa: F401   (deployment file, untracked)
    except Exception:
        import src
        stub = types.ModuleType("src.config")
        stub.CONFIG = {"annotate": {}}
        sys.modules["src.config"] = stub
        src.config = stub
    spec =importlib.util.spec_from_file_location("annotate_v2_r2", os.path.join(REPO, "tools", "annotate_v2.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_insertion_call_counts_as_carrier_only_in_tprt_mode(tmp_path):
    m = _annotate_v2()
    gt = tmp_path / "P.genotypes.csv.gz"
    with gzip.open(gt, "wt") as fh:
        fh.write("insertion;S1;S2;S3\n")
        fh.write("chr22:100-oneside_100;insertion;wild-type;wild-type\n")
        fh.write("chr22:oneside_500-500;heterozygous;insertion;wild-type\n")
        fh.write("chr22:900-40900;insertion;insertion;wild-type\n")
    A_ = m.CONFIG["annotate"]
    saved = {k: A_.get(k) for k in ("rte_library", "count_insertion_call")}

    def run():
        v = object.__new__(m.VariantAnnotationContainer)
        v.genotyping_file = str(gt)
        v.insertions = {t: m.Insertion(t, "", "acgt") for t in
                        ("chr22:100-oneside_100", "chr22:oneside_500-500", "chr22:900-40900")}
        v.read_genotyping()
        return {k: i.nins for k, i in v.insertions.items()}
    try:
        A_.pop("count_insertion_call", None)
        A_["rte_library"] = None                  # legacy: byte-identical behaviour
        assert run() == {"chr22:oneside_500-500": 1}
        A_["rte_library"] = "resources/rte_library"
        assert run() == {"chr22:100-oneside_100": 1, "chr22:oneside_500-500": 2, "chr22:900-40900": 2}
        A_["count_insertion_call"] = False
        assert run() == {"chr22:oneside_500-500": 1}
    finally:
        for k, v in saved.items():
            if v is None:
                A_.pop(k, None)
            else:
                A_[k] = v


# ------------------------------------------------------------------------- build_gene_model
def test_build_gene_model_reads_ucsc_refseq_gtf(tmp_path):
    gtf = tmp_path / "hs1.ncbiRefSeq.gtf"
    rows = [
        ("chr1", 101, 200, "+", 'gene_id "GA"; transcript_id "NM_1.1"; gene_name "GA";'),
        ("chr1", 301, 400, "+", 'gene_id "GA"; transcript_id "NM_1.1"; gene_name "GA";'),
        ("chr2", 501, 600, "-", 'gene_id "GB"; transcript_id "XM_2.1"; gene_name "GB";'),
        ("chrUn_x", 1, 50, "+", 'gene_id "GC"; transcript_id "NM_3.1"; gene_name "GC";'),
    ]
    with open(gtf, "w") as fh:
        for c, s, e, st, at in rows:
            fh.write(f"{c}\tncbiRefSeq\texon\t{s}\t{e}\t.\t{st}\t.\t{at}\n")
    out = tmp_path / "gm.tsv"
    subprocess.run([sys.executable, os.path.join(REPO, "tools", "build_gene_model.py"), "--curated",
                    str(gtf), str(out)], check=True, capture_output=True)
    got = [l.split("\t") for l in open(out).read().splitlines()]
    assert got == [["chr1", "100", "200", "GA", "+"], ["chr1", "300", "400", "GA", "+"]]
    subprocess.run([sys.executable, os.path.join(REPO, "tools", "build_gene_model.py"),
                    str(gtf), str(out)], check=True, capture_output=True)
    assert any(l.startswith("chr2\t500\t600\tGB\t-") for l in open(out).read().splitlines())
