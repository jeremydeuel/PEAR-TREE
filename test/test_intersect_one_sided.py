"""TPRT one-sided / non-TSD loci in combine's intersect_insertions
(plans/tprt_hallmarks/SPEC.md, "Pairing modes and locus names").

Discovery (`one_sided_loci`) names a one-sided locus `contig:L-oneside_L` (real LEFT
junction) or `contig:oneside_R-R` (real RIGHT junction). combine must keep it (pooled over
samples) as a one-real-side record so the genome-aware remap filters still run on it.
Target-site / L1-mediated deletions arrive as plain numeric names with left > right.

Run:  pytest test/test_intersect_one_sided.py
"""
import gzip
import os

from _combine_shim import TEST_CONFIG  # noqa: F401  (registers src/ on sys.path, imports combine)
import combine_insertions_intersect_insertions as ii
from combine_insertions_insertion import (Insertion, TYPE_FULL_INFO, TYPE_LEFT_DISC,
                                          TYPE_RIGHT_DISC)

FLANK = "ACGTTGCATGCAGTCAGTTGACCAGTAGGCATCGATCGGATCCATGCAAGT"  # 51 bp, aligned side
CLIP = "TTTTTTTTTTTTTTTTTTTTGCATCCAGTAAGCT"                       # stored LEFT clip (poly-A tail)


def write_discovery(path, records):
    """records: [(locus, [(field, seq)])] in discovery FASTQ order."""
    with gzip.open(path, "wt") as f:
        for locus, fields in records:
            for field, seq in fields:
                f.write(f"@{locus}:{field}\n{seq}\n+\n{'?' * len(seq)}\n")


def left_one_sided(locus, clip=CLIP, mates=("GGGGCCCCAAAATTTT",)):
    fields = [("LEFT:CLIPPED", clip), ("LEFT:ALIGNED", FLANK)]
    fields += [(f"LEFT:MATE{i}", m) for i, m in enumerate(mates)]
    return (locus, fields)


def right_one_sided(locus):
    return (locus, [("RIGHT:ALIGNED", FLANK), ("RIGHT:CLIPPED", "TTTTTTTTTTTTTTTTGACGTAGG")])


def full(locus):
    return (locus, [("LEFT:CLIPPED", "GACGTAGGCATGCAAT"), ("LEFT:ALIGNED", FLANK),
                    ("RIGHT:ALIGNED", FLANK), ("RIGHT:CLIPPED", "TTTTTTTTTTTTTTTTGACG")])


def parse(tmp_path, name, records):
    p = os.path.join(tmp_path, name)
    write_discovery(p, records)
    return list(Insertion.parseFile(p))


def test_oneside_tokens_parse_onto_one_real_side(tmp_path):
    a = parse(str(tmp_path), "S1.txt.gz", [left_one_sided("chr1:5000-oneside_5000"),
                                          right_one_sided("chr1:oneside_7000-7000")])
    left, right = a
    assert left.type is TYPE_RIGHT_DISC and left.open_side == "RIGHT"
    assert left.left_pos == 5000 and left.right_pos == 5000 and left.right_clipped is None
    assert str(left.left_consensus) == CLIP.lower() + FLANK  # consensus = clip + flank
    assert right.type is TYPE_LEFT_DISC and right.open_side == "LEFT"
    assert right.right_pos == 7000 and right.left_clipped is None


def test_one_sided_loci_are_kept_and_pooled_across_samples(tmp_path):
    t = str(tmp_path)
    s1 = parse(t, "S1.txt.gz", [left_one_sided("chr1:5000-oneside_5000", mates=("AAAACCCCGGGGTTTT",)),
                               full("chr1:9000-9012")])
    s2 = parse(t, "S2.txt.gz", [left_one_sided("chr1:5000-oneside_5000", clip=CLIP + "GGATC",
                                               mates=("CCCCAAAAGGGGTTTT", "TTTTGGGGAAAACCCC")),
                               right_one_sided("chr1:oneside_7000-7000")])
    out = {i.name: i for i in ii.intersect_insertions(s1 + s2, keep_polya_one_sided=False)}
    assert set(out) == {"chr1:5000-oneside_5000", "chr1:oneside_7000-7000", "chr1:9000-9012"}
    one = out["chr1:5000-oneside_5000"]
    assert one.type is TYPE_RIGHT_DISC and one.open_side == "RIGHT"
    assert len(one.left_clipped) == len(CLIP) + 5          # longest real clip is the representative
    assert len(one.left_mates) == 3                         # mates pooled over both samples
    assert sorted(one.files) == ["S1.txt.gz", "S2.txt.gz"]
    assert out["chr1:9000-9012"].type is TYPE_FULL_INFO
    # the one real side is what combine remaps / writes (no crash on the open side)
    assert one.left_consensus.fastq(f"{one.name}:L")


def test_legacy_polya_records_parked_unless_requested(tmp_path):
    rec = ("chr1:3000-polyA_3050", [("LEFT:CLIPPED", "GACGTAGGCATGCAATGG"), ("LEFT:ALIGNED", FLANK),
                                    ("RIGHT:CLIPPED_POLYA", "GCATGCAAGTCCAGT")])
    ins = parse(str(tmp_path), "S1.txt.gz", [rec])
    assert ii.intersect_insertions(list(ins), keep_polya_one_sided=False) == []
    ins = parse(str(tmp_path), "S1.txt.gz", [rec])
    kept = ii.intersect_insertions(list(ins), keep_polya_one_sided=True)
    assert len(kept) == 1
    k = kept[0]
    assert k.name == "chr1:3000-polyA_3050"                # name kept: sidecar link
    assert k.type is TYPE_RIGHT_DISC and k.open_side == "RIGHT" and k.right_pos == 3050
    assert k.right_clipped is None and k.left_pos == 3000


def test_target_site_deletion_name_is_a_full_record(tmp_path):
    # TSD deletion: RIGHT junction 15 bp LEFT of the LEFT junction -> plain numeric, left > right
    ins = parse(str(tmp_path), "S1.txt.gz", [full("chr1:5015-5000")])
    out = ii.intersect_insertions(ins, keep_polya_one_sided=False)
    assert len(out) == 1 and out[0].type is TYPE_FULL_INFO
    assert out[0].left_pos == 5015 and out[0].right_pos == 5000
