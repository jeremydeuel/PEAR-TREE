# PEAR-TREE - genotyping reads in annotate_v2's CORE classification (Insertion.gt_query_seqs /
# _gt_supplement, VariantAnnotationContainer.read_gt_core).
#
# The genotype2 extra pass collects reads of joint carriers that did not discover a locus; their
# inserted parts are Dfam-scanned / remapped next to the junction clips. Precedence under test:
# they resolve only an `unknown` / `artefact` junction call, and only when >= 2 GT reads agree
# with the re-run's class; a confident junction call never changes. Logic only (no nhmmscan /
# bowtie2 / pysam runtime), like test_annotate_sva.py.
#
# Run:  pytest test/test_annotate_gt_core.py
import gzip
import importlib.util
import os
import random
import sys
import types

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)
try:
    import pysam  # noqa: F401
except Exception:
    sys.modules["pysam"] = types.ModuleType("pysam")

_spec = importlib.util.spec_from_file_location("annotate_v2", os.path.join(REPO, "tools", "annotate_v2.py"))
annotate_v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(annotate_v2)
Insertion = annotate_v2.Insertion
Dfam_Annotation = annotate_v2.Dfam_Annotation
VAC = annotate_v2.VariantAnnotationContainer

RNG = random.Random(7)


def rnd(n):
    return "".join(RNG.choice("ACGT") for _ in range(n))


def _dfam(model, bits=40.0, strand="+", hmm_start=1):
    return Dfam_Annotation([model, "ACC", "q", str(bits), "1e-10", "0.0", str(hmm_start), str(hmm_start + 99),
                            strand, "1", "100", "1", "100", "6000", "desc"])


LFLANK, RFLANK = rnd(60), rnd(60)          # reference right after L / right before R
ALU = rnd(200)                             # stand-in inserted element sequence


def _ins(left_clip="a" * 20, right_clip=None):
    """junctions in combine's case convention: right = FLANK|ins, left = ins|FLANK (left poly-A)."""
    right_clip = right_clip if right_clip is not None else ALU[:12].lower()
    return Insertion("chr1:1000-1015", left_clip + LFLANK, RFLANK + right_clip)


def test_query_seqs_strip_flank_filter_dedup_and_cap():
    ins = _ins()
    a, b = ALU[:80], ALU[80:170]
    reads = [
        ("RIGHT", "GT_CLIP", RFLANK[-40:] + a),            # insert = a (flank stripped)
        ("RIGHT", "GT_CLIP", RFLANK[-40:] + a),            # duplicate
        ("LEFT", "GT_CLIP", b + LFLANK[:50]),              # insert = b
        ("LEFT", "GT_POLYA", "A" * 40 + LFLANK[:50]),      # poly-A only: dropped
        ("RIGHT", "GT_CLIP", rnd(20) + ALU[:50]),          # no flank anchor: insert unknown, skipped
        ("RIGHT", "GT_MATE", ALU[20:170]),                 # mate: whole read
        ("RIGHT", "GT_MATE", ALU[:25]),                    # too short
        ("RIGHT", "GT_DISC", RFLANK),                      # anchor = flank: never a query
        ("LEFT", "GT_MATE", "ACACACACACACACACACACACACACACACACAC"),   # low entropy
    ]
    q = ins.gt_query_seqs(reads)
    assert q == [("RIGHT", ALU[20:170]), ("LEFT", b), ("RIGHT", a)], q   # mates first, then longest
    annotate_v2.CONFIG["annotate"]["gt_core_max_queries"] = 2
    try:
        assert len(ins.gt_query_seqs(reads)) == 2
    finally:
        del annotate_v2.CONFIG["annotate"]["gt_core_max_queries"]


def test_no_gt_reads_is_unchanged():
    ins = _ins()
    base = ins.conclusion()
    assert ins.gt_core == "" and "[gt:" not in base
    ins.gt_queries = 0
    ins.gt_dfams = [("right", 0, _dfam("AluY"))]      # hits without queries are never used
    assert ins.conclusion() == base


def test_unknown_resolved_by_two_agreeing_gt_reads():
    ins = _ins()
    base = ins.conclusion()
    assert VAC.element_class(base) in ("unknown", "artefact")
    ins.gt_queries = 3
    ins.gt_dfams = [("right", 0, _dfam("AluY", 30)), ("right", 1, _dfam("AluSx", 25)),
                    ("right", 1, _dfam("L1PA2_3end", 10))]          # query 1's best is Alu
    c = ins.conclusion()
    assert c.startswith("polyA <- AluY") and c.endswith("[gt: +ALU from 2 genotyping reads]"), c
    assert VAC.element_class(c) == "ALU"
    assert ins.gt_core.startswith("resolved") and "ALU:2" in ins.gt_core
    # the junction hit lists are restored after the re-run
    assert ins.right_dfams == [] and ins.left_dfams == []


def test_one_gt_read_is_not_enough():
    ins = _ins()
    base = ins.conclusion()
    ins.gt_queries = 4
    ins.gt_dfams = [("right", 0, _dfam("AluY"))]
    assert ins.conclusion() == base
    assert ins.gt_core.startswith("unresolved") and "< support" in ins.gt_core


def test_gt_family_must_match_the_rerun_class():
    # GT hits say L1 on the right, but the right-L1 rule needs a LEFT poly-A this locus lacks
    ins = _ins(left_clip=rnd(20).lower())
    base = ins.conclusion()
    ins.gt_queries = 2
    ins.gt_dfams = [("right", 0, _dfam("L1HS_5end")), ("right", 1, _dfam("L1HS_5end"))]
    assert ins.conclusion() == base
    assert ins.gt_core.startswith("unresolved") and "no junction hallmark" in ins.gt_core


def test_confident_junction_call_never_changes():
    ins = _ins()
    ins.right_dfams = [_dfam("AluY", 50, "+")]
    base_ins = _ins()
    base_ins.right_dfams = [_dfam("AluY", 50, "+")]
    base = base_ins.conclusion()
    assert VAC.element_class(base) == "ALU"
    ins.gt_queries = 5
    ins.gt_dfams = [("right", k, _dfam("L1HS_5end", 80, "+", 3000)) for k in range(5)]
    assert ins.conclusion() == base
    assert ins.gt_core.startswith("discordant ALU: LINE1:5")
    ins.gt_dfams = [("left", k, _dfam("AluY", 30, "+")) for k in range(2)]
    assert ins.conclusion() == base and ins.gt_core.startswith("concordant ALU: ALU:2")


def test_l1_five_prime_extent_and_strands():
    ins = _ins()
    ins.right_dfams = [_dfam("L1HS_3end", 60, "+", 5800)]
    ins.gt_queries = 3
    ins.gt_dfams = [("right", 0, _dfam("L1HS_5end", 50, "+", 3100)),
                    ("right", 1, _dfam("L1HS_5end", 50, "-", 2500))]
    ins.conclusion()
    assert "L1 5' to hmm 2500 (junction 5800)" in ins.gt_core
    assert "L1 both strands on right (inverted segment)" in ins.gt_core


def test_element_class_strips_the_gt_note():
    assert VAC.element_class("polyA <- AluY + (1-100|6000) [gt: +ALU from 3 genotyping reads]") == "ALU"
    # an L1 token inside the note must not reclassify a non-element head
    assert VAC.element_class("unknown [gt: +LINE1 from 3 genotyping reads]") == "unknown"


def test_read_gt_fasta_keeps_gt_roles_only(tmp_path):
    p = tmp_path / "P.insertions.genotype_reads.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">chr1:1000-1015|RIGHT|GT_MATE|S2|f|2\nACGT\n>chr1:1000-1015|RIGHT|CLIP|S2|g|1\nTTTT\n"
                 ">chr9:1-2|LEFT|GT_CLIP|S2|h|1\nGGGG\n")
    got = VAC._read_gt_fasta(str(p), {"chr1:1000-1015"})
    assert got == {"chr1:1000-1015": [("RIGHT", "GT_MATE", "ACGT")]}


def test_alu_vote_blocked_by_any_sva_hit():
    # an SVA's Alu-like segment makes its GT reads vote ALU; any SVA hit blocks ALU adoption
    ins = _ins()
    base = ins.conclusion()
    ins.gt_queries = 3
    ins.gt_dfams = [("right", 0, _dfam("AluY", 30)), ("right", 1, _dfam("AluSx", 25)),
                    ("right", 2, _dfam("SVA_F", 20))]
    assert ins.conclusion() == base
    assert ins.gt_core.startswith("unresolved") and "SVA hits" in ins.gt_core


def test_streamed_bounded_selection_equals_gt_query_seqs():
    """read_gt_core keeps per locus only the gt_keep_bound() smallest distinct candidates while
    streaming the file; the selection must equal gt_query_seqs() over every read (many duplicate
    sequences under different roles / sides, more candidates than the bound)."""
    ins = _ins()
    pool = [ALU[i:i + 40 + (i % 50)] for i in range(0, 150, 3)]
    annotate_v2.CONFIG["annotate"]["gt_core_max_queries"] = 3
    try:
        bound = Insertion.gt_keep_bound()
        for trial in range(30):
            reads = []
            for _ in range(RNG.randint(1, 200)):
                s = RNG.choice(pool)
                role = RNG.choice(["GT_MATE", "GT_CLIP", "GT_POLYA", "GT_DISC"])
                side = RNG.choice(["LEFT", "RIGHT"])
                seq = (RFLANK[-30:] + s) if side == "RIGHT" else (s + LFLANK[:30])
                reads.append((side, role, seq if RNG.random() < 0.8 else s))
            lst = []
            for r in reads:
                c = ins.gt_query_candidate(*r)
                if c is not None:
                    lst.append(c)
                    if len(lst) > 2 * bound:
                        lst[:] = sorted(set(lst))[:bound]
            assert Insertion.gt_select(lst) == ins.gt_query_seqs(reads), trial
    finally:
        del annotate_v2.CONFIG["annotate"]["gt_core_max_queries"]
