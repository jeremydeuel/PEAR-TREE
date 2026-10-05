"""Combine-level TPRT filters from the E2E (plans/tprt_hallmarks/E2E_REPORT.md, problems #1/#3):
reference-tract slippage reject, far-pair (L1DEL/L1DUP) strictness, fuzzy cross-sample merge
with poly-A-aware clip agreement, and re-anchoring of fuzzy-merged evidence.

Run:  pytest test/test_tprt_combine_filters.py
"""
import gzip
import os
import random

import pytest

from _combine_shim import TEST_CONFIG  # noqa: F401  (registers src/ on sys.path)
import combine_insertions_evidence as ev
import combine_insertions_intersect_insertions as ii
import combine_insertions_tprt_filters as tf
from combine_insertions_evidence import EvidenceRow, apply_evidence
from combine_insertions_insertion import Insertion, TYPE_FULL_INFO, TYPE_RIGHT_DISC
from indel_consensus import revcomp

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
_R = random.Random(11)
UNIQUE = "".join(_R.choice("ACGT") for _ in range(200))
ALU_TAIL_RC = "GAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGG"   # rc of an Alu 3' end
CFG = {"polya_min_len": 8}


def line_with_tract(base="T", n=20, inward=True):
    """Outward reference window (junction at 100): a homopolymer tract of n bp on the inward
    side (or outward), unique sequence elsewhere."""
    left = UNIQUE[:100 - n] + base * n if inward else UNIQUE[:100]
    right = UNIQUE[100:200] if inward else base * n + UNIQUE[100:200 - n]
    return left + right, 100


# --------------------------------------------------------------------- slippage

def test_repeat_and_strip():
    line, j = line_with_tract("T", 20)
    a, b, u = tf.repeat_at(line, j)
    assert a <= 80 and (b, u) == (100, "T")
    assert tf.strip_repeat("TTTTTTTGACG", "T") == 7
    assert tf.strip_repeat("TTTATTTTTGACG", "T") == 9          # one sequencing error tolerated
    assert tf.strip_repeat("CACACACAGT", "AC") == 8           # any phase
    assert tf.structured_len("GATTACATTTTTTTTTTT") == 7


def test_slipped_tract_continuation_is_slippage():
    line, j = line_with_tract("T", 20)
    # extra T's, then the reference past the tract (shifted by the slip)
    assert tf.slippage_junction("TTTTT" + line[j:j + 30], line, j, CFG) == "repeat_shifted_reference"
    # nothing but the tract base
    assert tf.slippage_junction("T" * 15, line, j, CFG) == "repeat_only"
    # post-homopolymer SBS junk, still dominated by the tract base
    assert tf.slippage_junction("TTTTTTTTGATCATTTGTTTGTTTTATTCTATTTCTTATTC", line, j, CFG) == "repeat_junk"
    # STR continuation: (CA)n tract, clip = more CA + shifted reference
    s = UNIQUE[:76] + "CA" * 12 + UNIQUE[100:200]
    assert tf.slippage_junction("CACACA" + s[100:130], s, 100, CFG) == "repeat_shifted_reference"


def test_real_insertions_are_not_slippage():
    line, j = line_with_tract("T", 20)
    # poly-A tail (outward poly-T) continuing a reference T tract, then element sequence
    assert tf.slippage_junction("T" * 25 + ALU_TAIL_RC, line, j, CFG) == ""
    # element clip right next to the tract
    assert tf.slippage_junction("GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCA", line, j, CFG) == ""
    # POLYA_ONLY / orphan TD at a site WITHOUT a reference tract: never tested
    plain = UNIQUE[:200]
    assert tf.slippage_junction("T" * 30, plain, 100, CFG) == ""
    assert tf.slippage_junction("T" * 20 + "GATCATTTGTTTGTTTTATTC", plain, 100, CFG) == ""
    # TTTT/AA endonuclease motif (4-6 bp run) is below the 8 bp tract threshold
    short, j2 = line_with_tract("T", 6)
    assert tf.slippage_junction("T" * 30, short, j2, CFG) == ""
    # a clip that opens with its OWN poly-A of another base is not this tract's slippage
    a_line, j3 = line_with_tract("A", 14)
    assert tf.slippage_junction("T" * 12 + "AATTTAAGTTCTAGGGTACATGTGCAC", a_line, j3, CFG) == ""


class FakeMatcher:
    """hit() by substring: 'ELEM' marks element sense, 'MELE' antisense, 'ALUX' another class."""

    def hit(self, seq, with_flanks=False):
        s = seq or ""
        if "GGCCGGGCGCGG" in s:
            return ("L1HS", "L1", "+", 30, 0, 30)
        if revcomp("GGCCGGGCGCGG") in s:
            return ("L1HS", "L1", "-", 30, 0, 30)
        if "GCTCACGCCTGT" in s:
            return ("ALU_Y", "ALU", "+", 30, 0, 30)
        return None


class Rec:
    def __init__(self, side, clip, n_ind=3, seq_with_mates=None):
        from indel_consensus import ConsensusResult
        self.side = side
        self.n_independent = n_ind
        self.combined_consensus = ConsensusResult(seq=clip)
        self.consensus = ConsensusResult(seq=seq_with_mates or clip)
        self.rows = []


class FakeIns:
    def __init__(self, left_pos, right_pos, lc="", rc=""):
        self.reference_name = "chr1"
        self.left_pos, self.right_pos = left_pos, right_pos
        self.left_clipped, self.right_clipped = lc, rc


def test_slippage_check_needs_informative_other_junction():
    seq = UNIQUE[:100] + "T" * 20 + UNIQUE[100:200]        # ref-forward; T tract 100..120

    def fetch(c, s, e):
        return seq[max(0, s):e]
    # RIGHT junction at 120 (end of tract, reads aligned on the left): outward = ref-forward
    slipped = "TTTT" + seq[120:150]
    ins = FakeIns(130, 120)
    other_junk = Rec("LEFT", "ACGATCGTAGCTAGCTAGCATCGA")
    assert ev._slippage_check(ins, [other_junk, Rec("RIGHT", slipped)], CFG, FakeMatcher(), fetch).startswith("slippage:RIGHT")
    other_elem = Rec("LEFT", "GGCCGGGCGCGGTGGCTCACGCC")
    assert ev._slippage_check(ins, [other_elem, Rec("RIGHT", slipped)], CFG, FakeMatcher(), fetch) == ""
    # one-sided locus: no other junction -> rejected
    assert ev._slippage_check(ins, [Rec("RIGHT", slipped)], CFG, FakeMatcher(), fetch) != ""


# --------------------------------------------------------------------- far pairs

FAR_CFG = {"min_independent_fragments": 2}
ELEM = "GGCCGGGCGCGGTGGCTCAC"
TAIL = "T" * 20 + "ACGTAGCTAGG"


def verdict(clips, cols=({"A", "B"}, {"A", "B"}), mates=None, cfg=FAR_CFG):
    return tf.far_pair_verdict({"LEFT": [clips[0]], "RIGHT": [clips[1]]},
                               {"LEFT": cols[0], "RIGHT": cols[1]}, FakeMatcher(), cfg, mates)


def test_far_pair_credible():
    assert verdict((TAIL, ELEM)) == ("", "LEFT")
    assert verdict((ELEM, TAIL)) == ("", "RIGHT")


def test_far_pair_failures():
    assert verdict((ELEM, ELEM))[0] == "no_polarity"
    assert verdict((TAIL, TAIL))[0] == "no_polarity"
    assert verdict((TAIL, ELEM)) == ("", "LEFT")       # no pooled-fragment criterion any more
    assert verdict((TAIL, "ACGATCGATCGTAGCTAGCTAGCATCG")) == ("no_element_on_complex_clip", "LEFT")
    assert verdict((TAIL, revcomp(ELEM)))[0] == "element_antisense"
    assert verdict((TAIL, revcomp(ELEM)), cfg=dict(FAR_CFG, far_pair_allow_antisense=True))[0] == ""
    assert verdict(("T" * 20 + "GCTCACGCCTGTAATC", ELEM))[0] == "element_class_conflict"
    assert verdict((TAIL, ELEM), cols=({"A", "B", "C"}, {"B"}))[0] == "colony_mismatch"


def test_far_pair_short_complex_clip_rescued_by_inside_mates():
    short = "CAATGAGATCACATGGA"                  # 17 bp of element: too short for the library
    assert verdict((TAIL, short))[0] == "no_element_on_complex_clip"
    mates = {"RIGHT": ["AAAA" + ELEM + "CC", ELEM + "TTTT"]}
    assert verdict((TAIL, short), mates=mates)[0] == ""
    assert verdict((TAIL, short), mates={"RIGHT": [ELEM]})[0] == "no_element_on_complex_clip"  # one mate is not enough


def test_colonies_consistent():
    assert tf.colonies_consistent({"A", "B"}, {"A", "B"})
    assert tf.colonies_consistent({"A", "B"}, {"A", "B", "C"})          # one colony missed
    assert not tf.colonies_consistent({"A"}, {"B"})
    assert not tf.colonies_consistent({"B"}, {"A", "B", "C"})
    big = {f"S{i}" for i in range(20)}
    assert tf.colonies_consistent(big, big - {"S1", "S2", "S3"})


def test_far_pair_split_into_one_sided_and_absorbed(monkeypatch):
    """A failing far pair keeps its poly-A junction as a one-sided locus, which joins another
    colony's two-sided call of that junction."""
    cfg = {"min_independent_fragments": 2, "far_pair_strict": True,
           "merge_tolerance_bp": 8}

    def mkrow(locus, side, frag, sample, pos, clip, strand="+"):
        n = 60
        if side == "RIGHT":
            seq, cig, at = "C" * n + clip, f"{n}M{len(clip)}S", n
            p = pos - n
        else:
            seq, cig, at = revcomp(clip) + "G" * n, f"{len(clip)}S{n}M", len(clip)
            p = pos
        return EvidenceRow(sample, dict(locus=locus, side=side, role="CLIP", frag=frag, r12="1", flag="99",
                                        ref="chr1", pos=str(p), strand=strand, outer=str(p + int(frag[-1]) * 25 + (sample == "B")),
                                        mref="*", mpos="-1", mstrand="*", tlen="0", mapq="60", cigar=cig,
                                        clip_at=str(at), seq=seq, qual="I" * len(seq)))
    far, good = "chr1:5000-1000", "chr1:5000-4990"
    rows = {(("A.txt.gz", far), "LEFT"): [mkrow(far, "LEFT", f"a{i}", "A", 5000, TAIL) for i in range(3)],
            (("A.txt.gz", far), "RIGHT"): [mkrow(far, "RIGHT", f"b{i}", "A", 1000, "ACGATCGATCGTAGCTAGCTAGCATCG")
                                           for i in range(3)],
            (("B.txt.gz", good), "LEFT"): [mkrow(good, "LEFT", f"c{i}", "B", 5000, TAIL) for i in range(3)],
            (("B.txt.gz", good), "RIGHT"): [mkrow(good, "RIGHT", f"d{i}", "B", 4990, ELEM) for i in range(3)]}
    monkeypatch.setattr(ev, "load_evidence", lambda files, wanted: (rows, {"A.txt.gz", "B.txt.gz"}))
    import combine_insertions_tprt_filters as tfm
    monkeypatch.setattr(tfm, "LibraryMatcher", lambda *a, **k: FakeMatcher())

    def mk(name, lp, rp, f):
        i = FakeIns(lp, rp, TAIL, "x")
        i.name, i.files, i.type, i.open_side = name, [f], TYPE_FULL_INFO, None
        i.member_loci = [(f, name)]
        i.left_aligned = i.right_aligned = None
        i.left_mates, i.right_mates = [], []
        return i
    a, b = mk(far, 5000, 1000, "A.txt.gz"), mk(good, 5000, 4990, "B.txt.gz")
    kept, records, failed, stats = apply_evidence([a, b], ["A.txt.gz", "B.txt.gz"], cfg)
    assert sorted(k.name for k in kept) == sorted(["chr1:5000-oneside_5000", good]) and not failed
    # folded into the two-sided call only after the remap filters (combine_insertions)
    kept, n = stats["pool"].absorb_one_sided(kept)
    assert [k.name for k in kept] == [good] and n == 1 and "chr1:5000-oneside_5000" not in records
    left = [r for r in records[good] if r.side == "LEFT"][0]
    assert left.n_independent == 6 and left.n_samples == 2        # colony A's poly-A junction pooled
    # without a two-sided call to join, the split one-sided locus is kept on its own
    kept, records, failed, _ = apply_evidence([mk(far, 5000, 1000, "A.txt.gz")], ["A.txt.gz"], cfg)
    assert [k.name for k in kept] == ["chr1:5000-oneside_5000"] and kept[0].type == TYPE_RIGHT_DISC
    assert records["chr1:5000-oneside_5000"][0].fail_reason == "split_from_far_pair:no_element_on_complex_clip"
    # far_pair_split off: dropped
    kept, _, failed, _ = apply_evidence([mk(far, 5000, 1000, "A.txt.gz")], ["A.txt.gz"], dict(cfg, far_pair_split=False))
    assert kept == [] and failed == {far}


def test_library_matcher_real_library():
    pytest.importorskip("mappy")
    m = tf.LibraryMatcher(os.path.join(REPO, "resources", "rte_library"), flanks=False)
    l1 = {}
    with open(os.path.join(REPO, "resources", "rte_library", "consensus.fa")) as fh:
        name = None
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                l1[name] = ""
            else:
                l1[name] += line.strip()
    seg = l1["L1HS"][3000:3030]
    h = m.hit(seg)
    assert h is not None and h[1] == "L1" and h[2] == "+"
    assert m.hit(revcomp(seg))[2] == "-"
    assert m.hit(UNIQUE[:40]) is None
    assert m.hit(seg[:15]) is None                      # too short


# --------------------------------------------------------------------- fuzzy merge

FLANK = "ACGTTGCATGCAGTCAGTTGACCAGTAGGCATCGATCGGATCCATGCAAGT"


def write_discovery(path, records):
    with gzip.open(path, "wt") as f:
        for locus, fields in records:
            for field, seq in fields:
                f.write(f"@{locus}:{field}\n{seq}\n+\n{'?' * len(seq)}\n")


def full(locus, lclip="GACGTAGGCATGCAATCCGT", rclip="TTTTTTTTTTTTTTTTGACG"):
    return (locus, [("LEFT:CLIPPED", lclip), ("LEFT:ALIGNED", FLANK),
                    ("RIGHT:ALIGNED", FLANK), ("RIGHT:CLIPPED", rclip)])


def one_left(locus):
    return (locus, [("LEFT:CLIPPED", "GACGTAGGCATGCAATCCGT"), ("LEFT:ALIGNED", FLANK), ("LEFT:MATE0", "GGGGCCCC")])


def parse(tmp, name, recs):
    p = os.path.join(tmp, name)
    write_discovery(p, recs)
    return list(Insertion.parseFile(p))


def test_fuzzy_merge_nearby_breakpoints(tmp_path):
    t = str(tmp_path)
    recs = (parse(t, "S1.txt.gz", [full("chr1:1000-1015")]) + parse(t, "S2.txt.gz", [full("chr1:1000-1015")])
            + parse(t, "S3.txt.gz", [full("chr1:1002-1018")]))
    out = ii.intersect_insertions(recs, merge_tolerance_bp=0)
    assert sorted(i.name for i in out) == ["chr1:1000-1015", "chr1:1002-1018"]      # legacy: exact names
    recs = (parse(t, "S1.txt.gz", [full("chr1:1000-1015")]) + parse(t, "S2.txt.gz", [full("chr1:1000-1015")])
            + parse(t, "S3.txt.gz", [full("chr1:1002-1018")]))
    out = ii.intersect_insertions(recs, merge_tolerance_bp=5)
    assert [i.name for i in out] == ["chr1:1000-1015"]                              # heaviest name wins
    assert ("S3.txt.gz", "chr1:1002-1018") in out[0].member_loci
    assert set(out[0].files) == {"S1.txt.gz", "S2.txt.gz", "S3.txt.gz"}
    # too far apart / clips disagree: kept apart
    recs = parse(t, "S1.txt.gz", [full("chr1:1000-1015")]) + parse(t, "S2.txt.gz", [full("chr1:1020-1035")])
    assert len(ii.intersect_insertions(recs, merge_tolerance_bp=5)) == 2
    recs = (parse(t, "S1.txt.gz", [full("chr1:1000-1015")])
            + parse(t, "S2.txt.gz", [full("chr1:1001-1015", lclip="CCCCGGGGAAAATTTTCCCCGGGG")]))
    assert len(ii.intersect_insertions(recs, merge_tolerance_bp=5)) == 2


def test_fuzzy_merge_one_sided(tmp_path):
    """intersect merges nearby one-sided loci of the same real side; folding one into a
    two-sided call is left to EvidencePool (after the remap filters)."""
    t = str(tmp_path)
    recs = parse(t, "S1.txt.gz", [full("chr1:1000-1015")]) + parse(t, "S2.txt.gz", [one_left("chr1:1001-oneside_1001")])
    assert len(ii.intersect_insertions(recs, merge_tolerance_bp=5)) == 2
    recs = (parse(t, "S1.txt.gz", [one_left("chr1:1000-oneside_1000")])
            + parse(t, "S2.txt.gz", [one_left("chr1:1003-oneside_1003")]))
    out = ii.intersect_insertions(recs, merge_tolerance_bp=5)
    assert len(out) == 1
    m = [x for x in out[0].member_loci if x[0] != out[0].files[0] or x[1] != out[0].name]
    assert m and out[0].member_sides[m[0]] == ("LEFT",)
    recs = (parse(t, "S1.txt.gz", [one_left("chr1:1000-oneside_1000")])
            + parse(t, "S2.txt.gz", [one_left("chr1:1003-oneside_1003")]))
    assert len(ii.intersect_insertions(recs, merge_tolerance_bp=0)) == 2


def test_polya_aware_clip_agreement(tmp_path):
    """Same junction, poly-A of different length followed by SBS junk: the legacy column check
    drops the locus, the poly-A-aware one keeps it."""
    t = str(tmp_path)
    a = "TTTTTTTGAGACGGA"                      # E2E: chr22:31564814-31564831, two colonies
    b = "TTTTTTTTTGAGA"

    def recs():
        return (parse(t, "S1.txt.gz", [full("chr1:1000-1015", rclip=a)])
                + parse(t, "S2.txt.gz", [full("chr1:1000-1015", rclip=b)]))
    assert ii.intersect_insertions(recs(), polya_aware_clip_agreement=False) == []
    assert len(ii.intersect_insertions(recs(), polya_aware_clip_agreement=True)) == 1
    assert tf.clips_agree([a, b]) and not tf.clips_agree(["GATTACAGATTACAGG", "CCGGTTAACCGGTTAA"])


def test_reanchor_moves_clip_to_insertion_junction():
    r = EvidenceRow("S1", dict(locus="chr1:1000-1018", side="RIGHT", role="CLIP", frag="f", r12="1", flag="99",
                               ref="chr1", pos="958", strand="+", outer="958", mref="*", mpos="-1", mstrand="*",
                               tlen="0", mapq="60", cigar="60M20S", clip_at="60", seq="A" * 80, qual="I" * 80))
    same = ev._reanchor([r], "RIGHT", 1018)
    assert same[0] is r                                          # own locus: untouched
    moved = ev._reanchor([r], "RIGHT", 1015)[0]
    assert moved.clip_at == 57 and r.clip_at == 60               # 3 aligned bases join the clip
    assert ev._locus_junction("chr1:1000-oneside_1000", "LEFT") == 1000


# --------------------------------------------------------------------- shard store

SIDECAR_COLS = ["locus", "side", "role", "frag", "r12", "flag", "ref", "pos", "strand", "outer", "mref",
                "mpos", "mstrand", "tlen", "mapq", "cigar", "clip_at", "seq", "qual"]


def _shard_fixture(tmp_path, n=60):
    """n insertions over two colonies' sidecars (incl. one far pair that splits one-sided and is
    absorbed by the other colony's call of the same junction), every row also written to disk."""
    rng = random.Random(3)
    files = [str(tmp_path / f"{s}.txt.gz") for s in ("A", "B")]
    lines = {f: [] for f in files}
    specs = []

    def add(f, locus, side, frag, pos, clip):
        m = 60
        if side == "RIGHT":
            seq, cig, at, p = "C" * m + clip, f"{m}M{len(clip)}S", m, pos - m
        else:
            seq, cig, at, p = revcomp(clip) + "G" * m, f"{len(clip)}S{m}M", len(clip), pos
        outer = p + 25 * int(frag[1:]) + (7 if f.endswith("B.txt.gz") else 0)
        lines[f].append([locus, side, "CLIP", frag, "1", "99", "chr1", str(p), "+", str(outer), "*", "-1",
                         "*", "0", "60", cig, str(at), seq, "I" * len(seq)])
    for k in range(n):
        f = files[k % 2]
        lp = 1000 + 300 * k
        rp = lp - 10
        name = f"chr1:{lp}-{rp}"
        clip = "".join(rng.choice("ACGT") for _ in range(25))
        for j in range(rng.randint(1, 4)):
            add(f, name, "LEFT", f"l{j}", lp, TAIL)
            add(f, name, "RIGHT", f"r{j}", rp, ELEM + clip)
        specs.append((name, lp, rp, f))
    far, good = "chr1:90000-86000", "chr1:90000-89990"
    for j in range(3):
        add(files[0], far, "LEFT", f"a{j}", 90000, TAIL)
        add(files[0], far, "RIGHT", f"b{j}", 86000, "ACGATCGATCGTAGCTAGCTAGCATCG")
        add(files[1], good, "LEFT", f"c{j}", 90000, TAIL)
        add(files[1], good, "RIGHT", f"d{j}", 89990, ELEM)
    specs += [(far, 90000, 86000, files[0]), (good, 90000, 89990, files[1])]
    for f, ls in lines.items():
        with gzip.open(f + ".evidence.tsv.gz", "wt") as fh:
            fh.write("\t".join(SIDECAR_COLS) + "\n")
            for l in ls:
                fh.write("\t".join(l) + "\n")
    return files, specs


def _shard_insertions(specs):
    out = []
    for name, lp, rp, f in specs:
        i = FakeIns(lp, rp, TAIL, "x")
        i.name, i.files, i.type, i.open_side = name, [os.path.basename(f)], TYPE_FULL_INFO, None
        i.member_loci = [(os.path.basename(f), name)]
        i.left_aligned = i.right_aligned = None
        i.left_mates, i.right_mates = [], []
        out.append(i)
    return out


@pytest.mark.parametrize("threads,gate", [(1, False), (2, False), (2, True)])
def test_shard_store_matches_memory(tmp_path, monkeypatch, threads, gate):
    """Chunked evaluation from on-disk shards (serial and forked workers) gives exactly the
    in-memory result: same kept insertions, same records, same evidence outputs after absorb."""
    import combine_insertions_tprt_filters as tfm
    monkeypatch.setattr(tfm, "LibraryMatcher", lambda *a, **k: FakeMatcher())
    files, specs = _shard_fixture(tmp_path)
    cfg = {"min_independent_fragments": 2, "far_pair_strict": True, "merge_tolerance_bp": 8,
           "indel_aware_consensus": True, "require_independent_fragments": gate}

    def run(shard_dir, thr, tag):
        ins = _shard_insertions(specs)
        kept, records, failed, stats = apply_evidence(ins, files, cfg, ref_fetch=lambda c, s, e: "N" * (e - s),
                                                      threads=thr, shard_dir=shard_dir)
        kept, n_abs = stats["pool"].absorb_one_sided(kept)
        tsv, fa = str(tmp_path / f"{tag}.tsv.gz"), str(tmp_path / f"{tag}.fa.gz")
        ev.write_evidence_outputs(records, [i.name for i in kept] + sorted(failed), tsv, fa,
                                  store=stats.get("store"), preload_from=len(kept))
        stats["store"].cleanup()
        return ([(i.name, str(i.left_clipped), str(i.right_clipped)) for i in kept], n_abs, sorted(failed),
                gzip.open(tsv, "rb").read(), gzip.open(fa, "rb").read())
    mem = run(None, 1, "mem")
    shard = run(str(tmp_path / "shards"), threads, f"shard{threads}")
    assert mem[1] == 1                                   # the split far pair was absorbed
    assert bool(mem[2]) == gate                          # single-fragment loci dropped only by the gate
    assert shard == mem
    assert not os.path.exists(tmp_path / "shards")       # cleaned up
