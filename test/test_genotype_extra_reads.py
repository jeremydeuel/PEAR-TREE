# PEAR-TREE - genotype2 extra-pass reads (gt_extra_reads): the members table, the carrier-filtered
# merge (tools/genotype_extra_reads.py), their use as classification evidence in tools/rte, and
# that cluster/somatic_table.py's hard evidence rules never count them.
#
# Run:  pytest test/test_genotype_extra_reads.py
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rte_sim  # noqa: E402
from tools import genotype_extra_reads as gx  # noqa: E402
from tools.rte.annotator import RteAnnotator, default_gt_reads, gt_changes  # noqa: E402
from tools.rte.inputs import EvidenceRead  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO, "cluster"))
import somatic_table  # noqa: E402


def _w(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def _r(path):
    with gzip.open(path, "rt") as fh:
        return fh.read()


def test_members_table_from_reads_fa(tmp_path):
    fa = tmp_path / "P.insertions.reads.fa.gz"
    _w(fa, ">chr1:100-115|LEFT|CLIP|S1|a|1\nACGT\n>chr1:100-115|LEFT|MATE|S2|a|2\nACGT\n"
           ">chr1:100-115|RIGHT|CLIP|S1|b|1\nACGT\n>odd|x:5-6|RIGHT|DISC|S3|c|1\nAAAA\n")
    out = tmp_path / "P.members.tsv.gz"
    m = gx.build_members(str(fa), str(out))
    assert m == {"chr1:100-115": {"S1", "S2"}, "odd|x:5-6": {"S3"}}
    assert _r(out) == "locus\tsamples\nchr1:100-115\tS1,S2\nodd|x:5-6\tS3\n"
    assert gx.germline_loci(str(out)) is None, "no colony count: no germline skip"
    # 3 colonies in the run (one without reads anywhere): 2/3 > 0.5 germline, 1/3 not
    gx.build_members(str(fa), str(out), n_colonies=3)
    assert _r(out).startswith("#colonies 3\nlocus\tsamples\n")
    assert gx.germline_loci(str(out)) == {"chr1:100-115"}
    # with 4 colonies 2/4 is NOT > 0.5 (strict)
    gx.build_members(str(fa), str(out), n_colonies=4)
    assert gx.germline_loci(str(out)) == set()
    assert gx.germline_loci(str(out), 0.25) == {"chr1:100-115"}


def test_merge_keeps_joint_carriers_only(tmp_path):
    gdir = tmp_path / "genotypes"
    gdir.mkdir()
    # numeric joint matrix: S1 carries A (0.95), S2 does not (0.5); S2 carries B
    mat = tmp_path / "P.genotypes.csv.gz"
    _w(mat, ";S1;S2\nchr1:1-2;0.9500;0.5000\nchr2:3-4;0.0100;0.9900\n")
    _w(gdir / "S1.txt.gz.extra_reads.fa.gz",
       ">chr1:1-2|RIGHT|GT_CLIP|S1|f1|1\nACGT\n>chr1:1-2|RIGHT|GT_MATE|S1|f2|2\nGGGG\n"
       ">chr2:3-4|LEFT|GT_DISC|S1|f3|1\nTTTT\n")
    _w(gdir / "S2.txt.gz.extra_reads.fa.gz",
       ">chr1:1-2|LEFT|GT_CLIP|S2|f4|1\nCCCC\n>chr2:3-4|LEFT|GT_POLYA|S2|f5|2\nAAAAAAAAAAAA\n")
    _w(gdir / "S3.txt.gz.extra_reads.fa.gz.tmp.1", ">chr1:1-2|LEFT|GT_CLIP|S3|x|1\nA\n")   # ignored
    _w(gdir / "S1.txt.gz", "locus\n")                                                     # not a sidecar
    out = tmp_path / "P.insertions.genotype_reads.fa.gz"
    kept, dropped, dg, n = gx.merge(str(gdir), str(mat), str(out))
    assert (kept, dropped, dg, n) == (3, 2, 0, 2)
    assert _r(out) == (">chr1:1-2|RIGHT|GT_CLIP|S1|f1|1\nACGT\n>chr1:1-2|RIGHT|GT_MATE|S1|f2|2\nGGGG\n"
                       ">chr2:3-4|LEFT|GT_POLYA|S2|f5|2\nAAAAAAAAAAAA\n")
    # members table: chr1:1-2 discovered by 2 of 3 colonies = germline -> dropped (older sidecars)
    mem = tmp_path / "P.members.tsv.gz"
    _w(mem, "#colonies 3\nlocus\tsamples\nchr1:1-2\tS3,S4\nchr2:3-4\tS3\n")
    kept, dropped, dg, n = gx.merge(str(gdir), str(mat), str(out), members=str(mem))
    assert (kept, dropped, dg) == (1, 1, 3)
    assert _r(out) == ">chr2:3-4|LEFT|GT_POLYA|S2|f5|2\nAAAAAAAAAAAA\n"
    # a missing members file is no filter
    assert gx.merge(str(gdir), str(mat), str(out), members=str(tmp_path / "nope.gz"))[:3] == (3, 2, 0)
    assert default_gt_reads("/x/P.combined.txt.gz") == "/x/P.insertions.genotype_reads.fa.gz"


def _gt(reads):
    return [EvidenceRead(r.side, "GT_" + r.role, "S2", r.frag, r.r12, r.seq) for r in reads]


def test_annotate_uses_genotype_reads_and_marks_the_change(tmp_path):
    lib = rte_sim.library()
    l1 = lib.consensus["L1HS"]
    ins = l1[4000:lib.cons_end["L1HS"]] + "A" * 25
    # discovery saw nothing at the 5' junction (as in test_5p_unresolved_without_5p_junction_reads)
    inp, ev, genome, _ = rte_sim.build(ins, exclude=(0,))
    inp.right_seq = ""
    ev.junctions.pop("RIGHT")
    # a carrier that did not discover the locus has the 5' junction reads + mates
    _, full, _, _ = rte_sim.build(ins, samples=("S2",))
    gt = _gt([r for r in full.reads if r.side == "RIGHT"])
    fa = tmp_path / "P.insertions.genotype_reads.fa.gz"
    body = "".join(f">{inp.title}|{r.side}|{r.role}|{r.sample}|{r.frag}|{r.r12}\n{r.seq}\n" for r in gt)
    # a combine role in this file must be ignored (it is not junction evidence of the pipeline)
    body += f">{inp.title}|RIGHT|CLIP|S2|zz|1\n{gt[0].seq}\n"
    _w(fa, body)

    ann = RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome)
    base = ann.annotate(inp, ev)
    assert base.structure == "5P_UNRESOLVED" and base.gt_reads == 0

    ann.evidence = {inp.title: ev}
    ann.load_evidence(gt_reads_path=str(fa), wanted={inp.title})
    assert ann.has_gt_reads and all(r.role.startswith("GT_") for r in ann.gt_reads[inp.title])
    assert len(ann.gt_reads[inp.title]) == len(gt)
    rec = ann.annotate_all({inp.title: inp})[inp.title]
    assert rec.element == "L1" and rec.structure == "TRUNCATED_5P"
    assert abs(rec.detail["j5"] - 4000) <= 5
    assert rec.gt_reads > 0
    assert "structure:5P_UNRESOLVED>TRUNCATED_5P" in rec.gt_changed
    # the combine-only record of the same locus is untouched by the pooled run
    assert ann.evidence[inp.title].reads == ev.reads


def test_annotate_without_genotype_reads_is_unchanged():
    lib = rte_sim.library()
    ins = lib.consensus["L1HS"][5500:lib.cons_end["L1HS"]] + "A" * 25
    inp, ev, genome, _ = rte_sim.build(ins)
    a = RteAnnotator({"rte_library": rte_sim.FIX}, genome=genome)
    r1 = a.annotate(inp, ev)
    a.evidence = {inp.title: ev}
    a.load_evidence(gt_reads_path="/nonexistent.fa.gz")
    r2 = a.annotate_all({inp.title: inp})[inp.title]
    assert not a.has_gt_reads and r1.row() == r2.row() and r2.gt_reads == 0 and r2.gt_changed == ""
    assert gt_changes(r1, r2) == ""


def test_somatic_table_hard_rules_ignore_gt_roles():
    """>= 2 independent junction fragments per end come ONLY from combine roles: GT_* reads (even
    if they were mixed into the reads the report reads) never add support."""
    loc = "chr1:1000-1015"
    clip = "ACGTTGCAAGCTTAGGCATCGATCGGATCCTAGCTAGGCTAACGTAGCTAGCTAGGATCGATTTACG"
    one = [("LEFT", "CLIP", "S1", "f1", "1", clip), ("RIGHT", "CLIP", "S1", "f2", "1", clip)]
    gt = [(side, "GT_" + role, "S1", f"g{i}{side}", "1", clip[i:] + "ACGT" * i)
          for i in range(1, 6) for side in ("LEFT", "RIGHT") for role in ("CLIP", "POLYA", "DISC", "MATE")]
    base = somatic_table.junction_support(one, {}, loc)
    assert somatic_table.junction_support(one + gt, {}, loc) == base
    bad_base, _ = somatic_table.hard_rules(loc, {loc: one}, True)
    bad_gt, _ = somatic_table.hard_rules(loc, {loc: one + gt}, True)
    assert bad_base == bad_gt and any("independent junction fragments" in b for b in bad_gt)
    assert all(not r.startswith("GT_") for r in somatic_table.JUNCTION_ROLES)
