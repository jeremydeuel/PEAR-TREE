"""Pooled independent-fragment rule + evidence outputs of combine_insertions
(plans/tprt_hallmarks/SPEC.md, "Independence rule (combine)").

Unit tests on src/combine_insertions_evidence.py plus a tiny end-to-end fixture: two fake
samples' discovery `.txt.gz` + `.evidence.tsv.gz` -> combine_insertions() with bowtie2 /
2bit / liftover stubbed (pre-made empty BAMs skip bowtie2; a synthetic genome replaces the
2bit). Checks that without sidecars -- and with sidecars but legacy config -- combined.txt.gz
and genotyping.txt.gz are byte-identical (decompressed), and that the TPRT mode gates and
upgrades the consensus.

Run:  pytest test/test_independence.py
"""
import gzip
import os
import random

import pysam
import pytest

from _combine_shim import GENOME, TEST_CONFIG, combine_insertions as ci_mod, consensus
import combine_insertions_evidence as ev
from combine_insertions_evidence import (EVIDENCE_TSV_COLUMNS, EvidenceRow, apply_evidence,
                                         collapse_fragments, evaluate_junction,
                                         independent_clusters)
from indel_consensus import revcomp
from quality_seq import QualitySeq

HEADER = ["locus", "side", "role", "frag", "r12", "flag", "ref", "pos", "strand", "outer",
          "mref", "mpos", "mstrand", "tlen", "mapq", "cigar", "clip_at", "seq", "qual"]


def row(sample="S1", **kw):
    d = dict(locus="chr1:100-110", side="RIGHT", role="CLIP", frag="f1", r12="1", flag="99",
             ref="chr1", pos="50", strand="+", outer="50", mref="*", mpos="-1", mstrand="*",
             tlen="0", mapq="60", cigar="60M40S", clip_at="60",
             seq="A" * 60 + "GATTACAGATTACAGATTACAGATTACAGATTACAGATTA", qual="")
    d.update({k: str(v) for k, v in kw.items()})
    if not d["qual"]:
        d["qual"] = "I" * len(d["seq"])
    return EvidenceRow(sample, d)


def n_ind(rows):
    clusters, _, _ = independent_clusters(collapse_fragments(rows))
    return len(clusters)


# --------------------------------------------------------------------- unit rules

def test_same_fragment_mates_count_once():
    rows = [row(frag="f1", r12=1, role="CLIP"),
            row(frag="f1", r12=2, role="MATE", ref="*", pos=-1, outer=-1, seq="CCGGTTAA" * 10)]
    frags = collapse_fragments(rows)
    assert len(frags) == 1 and frags[0].mate is not None
    assert n_ind(rows) == 1


def test_mate_only_fragment_is_not_evidence():
    rows = [row(frag="f1"), row(frag="f9", role="MATE", r12=2)]
    assert n_ind(rows) == 1


def test_pcr_duplicate_same_outer_both_ends_collapses():
    rows = [row(frag="f1", mref="chr1", mpos=400, mstrand="-"),
            row(frag="f2", outer=51, mref="chr1", mpos=401, mstrand="-")]
    clusters, n_dup, _ = independent_clusters(collapse_fragments(rows))
    assert len(clusters) == 1 and n_dup == 1


def test_same_outer_different_mate_is_independent():
    rows = [row(frag="f1", mref="chr1", mpos=400, mstrand="-"),
            row(frag="f2", mref="chr1", mpos=470, mstrand="-")]
    assert n_ind(rows) == 2


def test_unmapped_mates_decided_by_mate_sequence():
    m1, m2 = "ACGTTGCA" * 10, "TTGACCAG" * 10
    dup = [row(frag="f1"), row(frag="f1", r12=2, role="MATE", ref="*", pos=-1, outer=-1, seq=m1),
           row(frag="f2"), row(frag="f2", r12=2, role="MATE", ref="*", pos=-1, outer=-1, seq=m1[:-1] + "T")]
    assert n_ind(dup) == 1
    ind = [row(frag="f1"), row(frag="f1", r12=2, role="MATE", ref="*", pos=-1, outer=-1, seq=m1),
           row(frag="f2"), row(frag="f2", r12=2, role="MATE", ref="*", pos=-1, outer=-1, seq=m2)]
    assert n_ind(ind) == 2


def test_different_outer_is_independent_and_dup_flag_ignored():
    rows = [row(frag="f1", flag=1024 + 99), row(frag="f2", outer=60, pos=60, flag=1024 + 99)]
    assert n_ind(rows) == 2
    # the 0x400 flag alone never collapses anything; identical coords+seq without it do
    rows = [row(frag="f1", flag=99), row(frag="f2", flag=99)]
    assert n_ind(rows) == 1


def test_cross_sample_independent_unless_identical():
    a = row(sample="S1", frag="f1", mref="chr1", mpos=400, mstrand="-")
    b = row(sample="S2", frag="f7", outer=51, mref="chr1", mpos=401, mstrand="-")
    clusters, n_dup, n_cross = independent_clusters(collapse_fragments([a, b]))
    assert len(clusters) == 2 and n_cross == 0          # +-2 bp is a dup only within a sample
    c = row(sample="S2", frag="f8", mref="chr1", mpos=400, mstrand="-")
    clusters, _, n_cross = independent_clusters(collapse_fragments([a, c]))
    assert len(clusters) == 1 and n_cross == 1          # exact identity -> flagged + collapsed


# --------------------------------------------------------------------- lenient dedup (2nd pass)

CLIP_TAIL = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTA"


def _shifted(frag, d_outer, d_mate, **kw):
    """A RIGHT CLIP fragment whose read starts d_outer bp later (same junction at 110) and
    whose mapped mate starts d_mate bp later."""
    n_al = 60 - d_outer
    return row(frag=frag, pos=50 + d_outer, outer=50 + d_outer, cigar=f"{n_al}M40S", clip_at=n_al,
               seq="A" * n_al + kw.pop("clip", CLIP_TAIL), mref="chr1", mpos=400 + d_mate,
               mstrand="-", **kw)


@pytest.mark.parametrize("d_outer,d_mate", [(3, 0), (0, 3), (-3, 3), (3, -3)])
def test_lenient_dup_start_end_shift_with_same_mate(d_outer, d_mate):
    rows = [_shifted("f1", 0, 0), _shifted("f2", d_outer, d_mate)]
    stats = {}
    clusters, n_dup, _ = independent_clusters(collapse_fragments(rows), stats=stats)
    assert len(clusters) == 1 and n_dup == 1 and stats == {"n_dup_coord": 1}


def test_lenient_dup_shift_but_different_mate_is_independent():
    rows = [_shifted("f1", 0, 0), _shifted("f2", 3, 60)]
    assert n_ind(rows) == 2


def test_lenient_dup_tolerance_is_configurable():
    rows = [_shifted("f1", 0, 0), _shifted("f2", 3, 3)]
    clusters, _, _ = independent_clusters(collapse_fragments(rows), ev.DedupParams(tol=2))
    assert len(clusters) == 2


def test_lenient_dup_polya_length_and_degenerate_tail():
    """Same molecule read twice: poly-A length differs (SBS jitter) and the sequence 3' of the
    poly-A is low-quality junk in one copy -> still a duplicate (unmapped mates, sequence path)."""
    mate = "ACGTTGCAGGTCCATGATCCAGTTGACCAGTAGGCATCGTTAGCATGCCTAGGATCA"
    c1 = "GGCTCACGCCTGTAATCCCAGC" + "A" * 18 + "GTCAGGATCCTGACCGTTAG"
    c2 = "GGCTCACGCCTGTAATCCCAGC" + "A" * 21 + "GTNAGCATGGTGACTTTAG"
    rows = [row(frag="f1", seq="C" * 60 + c1, cigar=f"60M{len(c1)}S"),
            row(frag="f1", r12=2, role="MATE", ref="*", pos=-1, outer=-1, seq=mate),
            row(frag="f2", outer=52, pos=52, seq="C" * 58 + c2, cigar=f"58M{len(c2)}S", clip_at=58),
            row(frag="f2", r12=2, role="MATE", ref="*", pos=-1, outer=-1, seq=mate[2:] + "TT")]
    stats = {}
    clusters, _, _ = independent_clusters(collapse_fragments(rows), stats=stats)
    assert len(clusters) == 1 and stats == {"n_dup_seq": 1}


def test_lenient_dup_real_clip_difference_is_independent():
    """Same coordinates, but the non-poly-A clip differs by more than the edit budget."""
    other = "GATTCCTGATCAGTTACAGGTTCAGATTACGAATTTGTTA"
    rows = [_shifted("f1", 0, 0), _shifted("f2", 2, 1, clip=other)]
    assert n_ind(rows) == 2


def test_pcr_copies_with_evidence_on_different_reads_collapse():
    """Copy 1: R2 is the junction CLIP read, its mate R1 lies in the flank. Copy 2 (jittered
    ends): R2 fell short of the junction, so R1 is the evidence (DISC) and R2's mate record is
    multi-mapped onto a paralog (MAPQ 0, opposite strand). Same molecule -> one fragment
    (E2E finding: an ART_LONG_TSD chimera survived the gate this way)."""
    r1 = "TGGGGCCCAATGGCCAAGCCTTTTCTTCCCAAATGTCAGGGTCCTGGCACCACAAGG"
    r2 = "GAGGTACAACTAGCATACAGTAAAGTGCATAAATCTTAAGTGCATAGCTTGATGATT"
    copy1 = [row(frag="c1", r12=2, flag=147, strand="-", pos=572, outer=630, cigar="40M17S", clip_at=40,
                 seq=r2, mref="chr1", mpos=505, mstrand="+"),
             row(frag="c1", r12=1, role="MATE", flag=99, strand="+", pos=505, outer=505, cigar="57M",
                 seq=r1, mref="chr1", mpos=572, mstrand="-")]
    copy2 = [row(frag="c2", r12=1, role="DISC", flag=65, strand="+", pos=506, outer=506, cigar="57M",
                 seq=r1, mref="chr9", mpos=900, mstrand="+"),
             row(frag="c2", r12=2, role="MATE", flag=129, ref="chr9", strand="+", pos=900, outer=900,
                 mapq=0, cigar="57M", seq=revcomp(r2[2:]) , mref="chr1", mpos=506, mstrand="+")]
    clusters, n_dup, _ = independent_clusters(collapse_fragments(copy1 + copy2))
    assert len(clusters) == 1 and n_dup == 1


# PD45886b_lo0002 13:46537573-46537585 RIGHT (2026-10-10): one chimeric molecule, two copies.
# Copy 1 (markdup-flagged primary, so only its unflagged SUPPLEMENTARY clip survives): R1
# supplementary at 13:46537545 41M110S, mate R2 at 13:46537303. Copy 2 (kept; its own clip
# failed MAPQ): R2 at 13:46537303 is the DISC read, R1 primary at 5:168032175 is the MATE.
_SUPP_SEQ = ("ACTTATTGGTTAATTTGTTTTTTTTTTTTTGAGACAGGGTCCAACTTAATCAGGAAAGAAAAACTAGAATTCTCAAGGACAAAAATCAC"
             "AAAAGCAGCAAGAACTAGGACACTGGGTCCTGCAGCAAGAGCAGCACCCAGCTGCCCTCACC")
_DISC_SEQ = ("CACAGATTAAGTAATAGCCAGTAAGTGGTGTGCCTACTGAAATCCAGATCATCCGGCCTTACAGACCAAGCTCTTAATCACTTTGTTA"
             "AAGACTCAACCTACACACCTGCATTTGGGTAGATGTCTAGGGGAAAGGGCTGCCAAAATACCC")
_MATE_SEQ = revcomp("AACTTATTGGTTAATTTGTTTTTTTTTTTTGAGACAGGGTCCAACTTAATCAGGAAAGAAAAACTAGAATTCTCAAGGACAAAAATCA"
                    "CAAAAGCAGCAAGAACTAGGACACTGGGTCCTGCAGCAAGAGCAGCACCCAGCTGCCCTCACC")


def _supp_twin(supp_flag=2129, mate_seq=_MATE_SEQ, disc_pos=46537302):
    loc = "13:46537573-46537585"
    copy1 = [row(locus=loc, frag="5df603eabed89565", r12=1, flag=supp_flag, ref="13", pos=46537544,
                 strand="-", outer=46537695, mref="13", mpos=46537302, mstrand="+", cigar="41M110S",
                 clip_at=41, seq=_SUPP_SEQ)]
    copy2 = [row(locus=loc, frag="87593034e6499aee", r12=2, role="DISC", flag=129, ref="13",
                 pos=disc_pos, strand="+", outer=disc_pos, mref="5", mpos=168032174, mstrand="+",
                 cigar="151M", clip_at=-1, seq=_DISC_SEQ),
             row(locus=loc, frag="87593034e6499aee", r12=1, role="MATE", flag=65, ref="5", pos=168032174,
                 strand="+", outer=168032174, mref="13", mpos=disc_pos, mstrand="+", cigar="115M36S",
                 clip_at=-1, seq=mate_seq)]
    return copy1 + copy2


def test_supplementary_clip_and_its_templates_other_read_collapse():
    clusters, n_dup, _ = independent_clusters(collapse_fragments(_supp_twin()))
    assert len(clusters) == 1 and n_dup == 1


def test_supplementary_twin_needs_position_and_read_sequence():
    assert n_ind(_supp_twin(disc_pos=46537320)) == 2                      # DISC 18 bp off the mate
    assert n_ind(_supp_twin(mate_seq="ACGTTGCAGGTCCATGAT" * 8)) == 2      # mate is another read
    assert n_ind(_supp_twin(supp_flag=81)) == 2                           # primary clip: rule off


def test_allele_forward_orientation_of_mates():
    """reads.fa promises allele-forward sequence: a MATE stored on the same strand as its
    partner (placed on a paralog in the opposite orientation) is reverse-complemented."""
    s = "ACGTTGCAGGTCCATGA"
    # partner forward (0x20 unset), mate reverse (0x10 set): FR geometry -> stored as is
    assert ev.allele_forward_seq(row(role="MATE", flag=0x1 | 0x10 | 0x80, seq=s)) == s
    # partner forward, mate also forward (paralog in the other orientation / unmapped as read)
    assert ev.allele_forward_seq(row(role="MATE", flag=0x1 | 0x8 | 0x80, seq=s)) == revcomp(s)
    # partner reverse (0x20), mate forward -> as is
    assert ev.allele_forward_seq(row(role="MATE", flag=0x1 | 0x20 | 0x80, seq=s)) == s
    # junction reads are never touched
    assert ev.allele_forward_seq(row(role="CLIP", flag=0x1 | 0x10, seq=s)) == s


def test_pooling_follows_combine_grouping(monkeypatch):
    """Evidence of every discovery locus merged into one Insertion pools (colony A
    chr1:100-110 + colony B chr1:101-110), not only rows whose locus id equals the name."""
    rows = {(("A.txt.gz", "chr1:100-110"), "RIGHT"): [row(sample="A", frag="a1")],
            (("B.txt.gz", "chr1:101-110"), "RIGHT"): [row(sample="B", frag="b1", outer=55, pos=55)],
            (("A.txt.gz", "chr1:100-110"), "LEFT"): [row(sample="A", side="LEFT", frag="a2", strand="-", outer=200)],
            (("B.txt.gz", "chr1:101-110"), "LEFT"): [row(sample="B", side="LEFT", frag="b2", strand="-", outer=210)]}

    class Ins:
        name = "chr1:100-110"
        files = ["A.txt.gz", "B.txt.gz"]
        member_loci = [("A.txt.gz", "chr1:100-110"), ("B.txt.gz", "chr1:101-110")]
        right_aligned = left_aligned = None

    seen = {}

    def fake_load(files, wanted):
        seen["wanted"] = set(wanted)
        return rows, {"A.txt.gz", "B.txt.gz"}

    monkeypatch.setattr(ev, "load_evidence", fake_load)
    cfg = {"min_independent_fragments": 2}
    kept, records, failed, _ = apply_evidence([Ins()], ["A.txt.gz", "B.txt.gz"], cfg)
    assert set(Ins.member_loci) <= seen["wanted"]
    assert len(kept) == 1 and not failed
    r = [x for x in records["chr1:100-110"] if x.side == "RIGHT"][0]
    assert r.n_independent == 2 and r.n_samples == 2
    assert r.member_loci == "chr1:100-110,chr1:101-110"


def test_insertion_iadd_records_member_loci():
    from combine_insertions_insertion import Insertion
    q = lambda s: QualitySeq(s, [30] * len(s))
    data = lambda: {"LEFT:MATE": [], "RIGHT:MATE": [], "LEFT:CLIPPED": q("acgtacgtacgt"),
                    "LEFT:ALIGNED": q("ACGTACGTAAAA"), "RIGHT:CLIPPED": q("ttttacgtacgt"),
                    "RIGHT:ALIGNED": q("GGGGACGTACGT")}
    a = Insertion("chr1", "100", "110", data(), "/x/A.txt.gz")
    b = Insertion("chr1", "101", "110", data(), "/x/B.txt.gz")
    a += b
    assert a.member_loci == [("A.txt.gz", "chr1:100-110"), ("B.txt.gz", "chr1:101-110")]


# --------------------------------------------------------------------- SHORT overhang reads

_SR = random.Random(77)
SREF = "".join(_SR.choice("ACGT") for _ in range(2000))
SREF = SREF[:994] + "AAAAAA" + SREF[1000:]        # reference A-tract ending at 1000 (inward)
SJ = 1300                                          # RIGHT junction used by most cases
ELEM_CLIP = "GGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGC"
SHORT_CFG = {"count_short_overhang": True, "min_independent_fragments": 2}


def _sfetch(contig, start, end):
    return SREF[max(0, start):end]


def _clip_frag(frag="c1", j=SJ, clip=ELEM_CLIP, n_al=60, sample="S1"):
    return row(sample=sample, frag=frag, ref="chr1", pos=j - n_al, outer=j - n_al,
               cigar=f"{n_al}M{len(clip)}S", clip_at=n_al, seq=SREF[j - n_al:j] + clip)


def _short_frag(over, frag="s1", j=SJ, n_al=90, clipped=True, sample="S1"):
    cig = f"{n_al}M{len(over)}S" if clipped else f"{n_al + len(over)}M"
    return row(sample=sample, frag=frag, role="SHORT", ref="chr1", pos=j - n_al, outer=j - n_al,
               cigar=cig, clip_at=n_al, seq=SREF[j - n_al:j] + over)


def _short_junction(rows, cfg=SHORT_CFG):
    return evaluate_junction("chr1:1290-1300", "RIGHT", rows, cfg, _sfetch)


def test_short_overhang_matching_consensus_counts():
    assert sum(a != b for a, b in zip(ELEM_CLIP[:8], SREF[SJ:SJ + 8])) >= 2
    rec = _short_junction([_clip_frag(), _short_frag(ELEM_CLIP[:8])])
    assert rec.n_short_used == 1 and rec.n_independent == 2 and rec.n_independent_no_short == 1
    assert rec.n_short_mate_inside == 1          # mate unmapped (no MATE row, mref '*')
    # off by default: SHORT rows are invisible
    rec = _short_junction([_clip_frag(), _short_frag(ELEM_CLIP[:8])], {"min_independent_fragments": 2})
    assert rec.n_short_used == 0 and rec.n_independent == 1 and all(r.role != "SHORT" for r in rec.rows)


def test_short_overhang_too_short_not_counted():
    rec = _short_junction([_clip_frag(), _short_frag(ELEM_CLIP[:3])])
    assert rec.n_short_used == 0 and rec.n_short_rejected == 1 and rec.n_independent == 1
    assert rec.short_reasons == {"overhang_too_short": 1}


def test_short_overhang_matching_reference_not_counted():
    rec = _short_junction([_clip_frag(), _short_frag(SREF[SJ:SJ + 8], clipped=False)])
    assert rec.n_short_used == 0 and rec.n_short_rejected == 1
    # even when the consensus itself starts like the reference (microhomology), an overhang
    # that equals the reference proves nothing
    rec = _short_junction([_clip_frag(clip=SREF[SJ:SJ + 8] + ELEM_CLIP),
                           _short_frag(SREF[SJ:SJ + 8], clipped=False)])
    assert rec.n_short_used == 0 and rec.short_reasons == {"matches_reference": 1}


def test_short_polya_overhang_next_to_reference_a_tract_not_counted():
    clip = "A" * 20 + ELEM_CLIP
    rec = evaluate_junction("chr1:990-1000", "RIGHT",
                            [_clip_frag(j=1000, clip=clip), _short_frag("A" * 8, j=1000)], SHORT_CFG, _sfetch)
    assert rec.n_short_used == 0 and rec.short_reasons == {"ref_homopolymer": 1}


def test_short_read_through_slipped_reference_not_counted():
    """Slippage: the CLIP reads' clip is shifted reference (a slipped homopolymer); an unclipped
    SHORT read carries the slip as a deletion next to the junction. Its overhang IS reference
    (only shifted by the deletion) and must not rescue the junction (E2E finding)."""
    slipped = SREF[SJ + 3:SJ + 43]
    short = row(frag="s1", role="SHORT", ref="chr1", pos=SJ - 90, outer=SJ - 90, cigar="90M3D10M",
                clip_at=90, seq=SREF[SJ - 90:SJ] + SREF[SJ + 3:SJ + 13])
    rec = _short_junction([_clip_frag(clip=slipped), short])
    assert rec.n_short_used == 0 and rec.short_reasons == {"matches_reference": 1}


def test_short_only_junction_not_supported():
    rec = _short_junction([_short_frag(ELEM_CLIP[:8]), _short_frag(ELEM_CLIP[:9], frag="s2", n_al=80)])
    assert rec.n_independent == 0 and rec.short_reasons == {"no_clip_fragment": 2}


def test_short_read_duplicate_of_clip_fragment_collapses():
    """A SHORT read of the same molecule as the CLIP read (start shifted 2 bp, same bases)."""
    short = row(frag="s1", role="SHORT", ref="chr1", pos=SJ - 58, outer=SJ - 58, cigar="58M8S",
                clip_at=58, seq=SREF[SJ - 58:SJ] + ELEM_CLIP[:8])
    rec = _short_junction([_clip_frag(), short])
    assert rec.n_short_used == 1 and rec.n_independent == 1 and rec.n_duplicates == 1


def test_one_sided_locus_evaluates_only_its_real_side(monkeypatch):
    """`oneside_` / Feature-A loci: the open end has no reads by construction; only the real
    side is evaluated. Below 2 pooled fragments it is reported (supported=0), never dropped."""
    rows = {(("A.txt.gz", "chr1:100-oneside_100"), "LEFT"): [
        row(sample="A", side="LEFT", frag="l1", strand="-", outer=200),
        row(sample="A", side="LEFT", frag="l2", strand="-", outer=230)]}
    monkeypatch.setattr(ev, "load_evidence", lambda files, wanted: (rows, {"A.txt.gz"}))
    cfg = {"min_independent_fragments": 2}

    def mk(open_side, typ):
        class Ins:
            name = "chr1:100-oneside_100"
            files = ["A.txt.gz"]
            right_aligned = left_aligned = None
        Ins.open_side, Ins.type = open_side, typ
        return Ins()
    kept, records, failed, _ = apply_evidence([mk("RIGHT", 4)], ["A.txt.gz"], cfg)
    assert len(kept) == 1 and [r.side for r in records["chr1:100-oneside_100"]] == ["LEFT"]
    kept, _, _, _ = apply_evidence([mk(None, 4)], ["A.txt.gz"], cfg)      # Feature-A disc end
    assert len(kept) == 1
    rows[(("A.txt.gz", "chr1:100-oneside_100"), "LEFT")].pop()
    kept, records, failed, _ = apply_evidence([mk("RIGHT", 4)], ["A.txt.gz"], cfg)
    assert len(kept) == 1 and not failed and records["chr1:100-oneside_100"][0].supported == 0


def test_far_flank_trim_for_short_insertions():
    """A short insertion spanned by the junction reads: the clip = insert + far flank; the
    remap copy is cut where the far flank starts (else it maps next to the breakpoint)."""
    from combine_insertions_insertion import Insertion
    q = lambda s: QualitySeq(s, [30] * len(s))
    g = GENOME["chr1"]
    L, R = 5000, 5012                      # 12 bp TSD
    ins_seq = "GGCTCACGCCTGTAATCC" + "A" * 15
    right_clip = ins_seq + g[L:L + 30]     # outward from R ... runs into ref[L:]
    left_clip_fwd = g[R - 30:R] + ins_seq  # ref-forward: ref[..R] + insert (as discovery writes)
    data = {"LEFT:MATE": [], "RIGHT:MATE": [],
            "LEFT:CLIPPED": q(left_clip_fwd.lower()), "LEFT:ALIGNED": q(g[L:L + 40]),
            "RIGHT:CLIPPED": q(right_clip.lower()), "RIGHT:ALIGNED": q(g[R - 40:R])}
    i = Insertion("chr1", str(L), str(R), data, "/x/A.txt.gz")
    assert str(ci_mod._far_flank_trimmed(i, "R")).upper() == ins_seq
    assert str(ci_mod._far_flank_trimmed(i, "L")).upper() == revcomp(ins_seq)
    # long insertion (no far flank inside the clip) is untouched
    data["RIGHT:CLIPPED"] = q((ins_seq * 3).lower())
    i = Insertion("chr1", str(L), str(R), data, "/x/A.txt.gz")
    assert str(ci_mod._far_flank_trimmed(i, "R")).upper() == ins_seq * 3


def test_polya_end_is_reported_not_gated(monkeypatch):
    """RIGHT has 2 independent fragments, the poly-A (LEFT, outward poly-T) end only one: the
    junction is reported supported=0, the insertion kept (no pooled fragment gate)."""
    elem = "GGCTCACGCCTGTAATCCCGGATCCAGT"
    left_clip_fwd = elem[::-1][:0] + revcomp("T" * 15 + revcomp(elem))  # ref-forward [clip]
    rows = {
        ("chr1:100-110", "RIGHT"): [row(frag="r1", seq="C" * 60 + elem + "TTGCA" * 3),
                                    row(frag="r2", outer=40, pos=40, seq="C" * 60 + elem + "TTGCA" * 3)],
        ("chr1:100-110", "LEFT"): [row(side="LEFT", frag="l1", strand="-", outer=200, cigar="43S60M",
                                       clip_at=len(left_clip_fwd), seq=left_clip_fwd + "G" * 60)],
    }

    class Ins:
        name = "chr1:100-110"
        files = ["S1.txt.gz"]
        right_aligned = None
        left_aligned = None

    monkeypatch.setattr(ev, "load_evidence", lambda files, wanted: (rows, {"S1.txt.gz"}))
    cfg = {"min_independent_fragments": 2}
    kept, records, failed, stats = apply_evidence([Ins()], ["S1.txt.gz"], cfg)
    assert len(kept) == 1 and not failed and not stats["reasons"]
    left = [r for r in records["chr1:100-110"] if r.side == "LEFT"][0]
    assert left.polya_end == 1 and left.supported == 0
    # the optional pooled gate drops it
    kept, records, failed, stats = apply_evidence([Ins()], ["S1.txt.gz"], dict(cfg, require_independent_fragments=True))
    assert kept == [] and failed == {"chr1:100-110"} and stats["reasons"] == {"LEFT(polyA)": 1}


def test_left_clip_orientation():
    r = row(side="LEFT", seq="ACGTT" + "G" * 20, clip_at=5, cigar="5S20M")
    s, q = r.outward_clip()
    assert s == revcomp("ACGTT")
    r = row(side="LEFT", seq="ACGTT" + "G" * 20, clip_at=-1, cigar="5S20M")
    assert r.outward_clip()[0] == revcomp("ACGTT")


# --------------------------------------------------------------------- end-to-end

ELEM = "GGCTCACGCCTGTAATCCCG"
BEYOND = "GTCAGGATCCTGACCGTTAGCCAGTCTTGC"
READ_LEN = 100


def _ins_seq(pa):
    return ELEM + "A" * pa + BEYOND


def _q(n):
    return "".join(chr(33 + 35) for _ in range(n))


class Sample:
    def __init__(self, name):
        self.name = name
        self.rows = []
        self.clips = {}  # (locus, side) -> [outward clip]
        self.n = 0

    def frag_id(self):
        self.n += 1
        return f"{hash((self.name, self.n)) & 0xFFFFFFFFFFFFFFFF:016x}"

    def add_right(self, locus, jpos, pa, n_al, mate_off, frag=None):
        """RIGHT junction read: [aligned ref ending at jpos][clip into the insertion]."""
        frag = frag or self.frag_id()
        ins = _ins_seq(pa)
        aligned = GENOME["chr1"][jpos - n_al:jpos]
        clip = ins[:READ_LEN - n_al]
        seq = aligned + clip
        self.rows.append(dict(locus=locus, side="RIGHT", role="CLIP", frag=frag, r12=1, flag=73,
                              ref="chr1", pos=jpos - n_al, strand="+", outer=jpos - n_al, mref="*",
                              mpos=-1, mstrand="*", tlen=0, mapq=60, cigar=f"{n_al}M{len(clip)}S",
                              clip_at=n_al, seq=seq, qual=_q(len(seq))))
        downstream = ins + GENOME["chr1"][jpos - 15:jpos + 400]
        mate = downstream[mate_off:mate_off + READ_LEN]
        self.rows.append(dict(locus=locus, side="RIGHT", role="MATE", frag=frag, r12=2, flag=133,
                              ref="*", pos=-1, strand="*", outer=-1, mref="chr1", mpos=jpos - n_al,
                              mstrand="+", tlen=0, mapq=0, cigar="*", clip_at=-1, seq=mate,
                              qual=_q(len(mate))))
        self.clips.setdefault((locus, "RIGHT"), []).append(clip)
        return frag

    def add_left(self, locus, jpos, pa, n_al, frag=None):
        """LEFT junction read: [clip = insertion end][aligned ref starting at jpos]."""
        frag = frag or self.frag_id()
        ins = _ins_seq(pa)
        clip_fwd = ins[-(READ_LEN - n_al):]
        aligned = GENOME["chr1"][jpos:jpos + n_al]
        seq = clip_fwd + aligned
        self.rows.append(dict(locus=locus, side="LEFT", role="CLIP", frag=frag, r12=1, flag=89,
                              ref="chr1", pos=jpos, strand="-", outer=jpos + n_al, mref="*", mpos=-1,
                              mstrand="*", tlen=0, mapq=60, cigar=f"{len(clip_fwd)}S{n_al}M",
                              clip_at=len(clip_fwd), seq=seq, qual=_q(len(seq))))
        self.clips.setdefault((locus, "LEFT"), []).append(revcomp(clip_fwd))
        return frag

    def write(self, d, with_sidecar=True):
        """discovery-like `.txt.gz`: per locus the sample's own (legacy column-wise) clip
        consensus, as discovery emits it -- truncated by the poly-A jitter."""
        txt = os.path.join(d, f"{self.name}.txt.gz")
        loci = sorted({k[0] for k in self.clips})
        with gzip.open(txt, "wt") as f:
            for locus in loci:
                l, r = (int(x) for x in locus.split(":")[1].split("-"))
                lc = consensus.find_consensus([QualitySeq(c, [30] * len(c)) for c in self.clips[(locus, "LEFT")]])
                rc = consensus.find_consensus([QualitySeq(c, [30] * len(c)) for c in self.clips[(locus, "RIGHT")]])
                la = GENOME["chr1"][l:l + 40]
                ra = GENOME["chr1"][r - 40:r]
                f.write(lc.revcomp().fastq(f"{locus}:LEFT:CLIPPED"))
                f.write(QualitySeq(la, [30] * 40).fastq(f"{locus}:LEFT:ALIGNED"))
                f.write(QualitySeq(ra, [30] * 40).fastq(f"{locus}:RIGHT:ALIGNED"))
                f.write(rc.fastq(f"{locus}:RIGHT:CLIPPED"))
        if with_sidecar:
            with gzip.open(os.path.join(d, f"{self.name}.evidence.tsv.gz"), "wt") as f:
                f.write("\t".join(HEADER) + "\n")
                for rr in self.rows:
                    f.write("\t".join(str(rr[h]) for h in HEADER) + "\n")
        return txt


X = "chr1:10000-10015"   # in both samples -> supported
Y = "chr1:20000-20012"   # sample A only; RIGHT = two PCR duplicates -> 1 independent


def build_fixture(d, with_sidecar=True):
    rng = random.Random(5)
    a, b = Sample("sampleA"), Sample("sampleB")
    # X RIGHT: A has 4 fragments, two of them PCR duplicates (same read + same mate);
    # poly-A lengths differ within each sample so discovery's column-wise consensus
    # (emulated in Sample.write) truncates at the poly-A, as on real data. Unmapped mates start
    # past the poly-A at distinct offsets: the lenient dedup compares sequence only up to a long
    # poly-A, so mates inside ELEM + poly-A would be indistinguishable (= duplicates).
    a.add_right(X, 10015, 18, 12, 45)
    f = a.add_right(X, 10015, 20, 14, 60)
    a.add_right(X, 10015, 20, 14, 60)       # PCR dup of f (new qname, same molecule)
    a.add_right(X, 10015, 16, 15, 75)
    b.add_right(X, 10015, 16, 13, 48)
    b.add_right(X, 10015, 19, 12, 85)
    # X LEFT: one fragment per sample
    a.add_left(X, 10000, 17, 30)
    b.add_left(X, 10000, 21, 34)
    # Y: sample A only
    a.add_right(Y, 20012, 18, 20, 4)
    a.add_right(Y, 20012, 18, 20, 4)        # duplicate
    a.add_left(Y, 20000, 18, 30)
    a.add_left(Y, 20000, 19, 40)
    files = [a.write(d, with_sidecar), b.write(d, with_sidecar)]
    return files


def run_combine(d, files, monkeypatch, cfg_extra=None):
    monkeypatch.setattr(ci_mod.pyliftover, "LiftOver", lambda path: None)
    cfg = TEST_CONFIG["combine_insertions"]
    saved = dict(cfg)
    cfg.update(cfg_extra or {})
    stem = os.path.join(d, "PT")
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": len(GENOME["chr1"])}]}
    for p in (f"{stem}.bam", f"{stem}.insertionsonly.bam"):
        with pysam.AlignmentFile(p, "wb", header=hdr):
            pass
    try:
        ci_mod.combine_insertions(files, f"{stem}.genotyping.txt.gz", f"{stem}.combined.txt.gz",
                                  f"{stem}.fq.gz", f"{stem}.bam", 1)
    finally:
        cfg.clear()
        cfg.update(saved)
    out = {}
    for k in ("combined.txt.gz", "genotyping.txt.gz", "insertions.evidence.tsv.gz", "insertions.reads.fa.gz"):
        p = f"{stem}.{k}"
        out[k] = gzip.open(p, "rb").read() if os.path.exists(p) else None
    return out


def _fastq_records(blob):
    lines = blob.decode().splitlines()
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines), 4)}


def _tsv(blob):
    lines = blob.decode().splitlines()
    h = lines[0].split("\t")
    return [dict(zip(h, l.split("\t"))) for l in lines[1:]]


def test_end_to_end_byte_identity_and_tprt_mode(tmp_path, monkeypatch):
    d0, d1, d2 = (tmp_path / x for x in ("nosidecar", "legacycfg", "tprt"))
    for x in (d0, d1, d2):
        x.mkdir()
    legacy = run_combine(str(d0), build_fixture(str(d0), with_sidecar=False), monkeypatch)
    assert legacy["insertions.evidence.tsv.gz"] is None and legacy["insertions.reads.fa.gz"] is None
    recs = _fastq_records(legacy["combined.txt.gz"])
    assert {k.rsplit(":", 1)[0] for k in recs} == {X, Y}
    # legacy longest-clip consensus loses the beyond-poly-A sequence
    assert BEYOND.lower() not in recs[f"{X}:R"]

    # sidecars present, legacy config -> legacy outputs byte-identical (decompressed;
    # the gzip header carries an mtime), new evidence files written
    side = run_combine(str(d1), build_fixture(str(d1)), monkeypatch)
    assert side["combined.txt.gz"] == legacy["combined.txt.gz"]
    assert side["genotyping.txt.gz"] == legacy["genotyping.txt.gz"]
    rows = _tsv(side["insertions.evidence.tsv.gz"])
    assert list(rows[0].keys()) == EVIDENCE_TSV_COLUMNS
    by = {(r["insertion_id"], r["side"]): r for r in rows}
    assert by[(X, "RIGHT")]["n_independent"] == "5"       # A: 4 frags - 1 dup; B: 2
    assert by[(X, "RIGHT")]["n_duplicates"] == "1"
    assert by[(X, "RIGHT")]["n_samples"] == "2"
    assert by[(X, "RIGHT")]["n_mates"] == "6"
    assert by[(X, "LEFT")]["n_independent"] == "2"        # pooled 1 + 1
    assert by[(Y, "RIGHT")]["n_independent"] == "1"
    assert by[(X, "RIGHT")]["supported"] == "1" and by[(Y, "RIGHT")]["supported"] == "0"
    # beyond-poly-A recovered exactly; the overlapping (unmapped) mates even extend it past
    # the insertion into the TSD/flank that follows it
    beyond = by[(X, "RIGHT")]["beyond_polya"]
    assert beyond.startswith(BEYOND.lower())
    assert GENOME["chr1"][10000:10400].lower().startswith(beyond[len(BEYOND):])
    assert len(beyond) > len(BEYOND)
    assert int(by[(X, "RIGHT")]["beyond_polya_support"]) >= 2
    assert BEYOND.lower() in by[(X, "RIGHT")]["clip_consensus"]
    assert by[(X, "RIGHT")]["polya_len_range"] in ("16-20", "16-19", "18-20")
    assert len(by[(X, "RIGHT")]["consensus_depth"].split(",")) == sum(
        c.islower() for c in by[(X, "RIGHT")]["clip_consensus"])
    fa = side["insertions.reads.fa.gz"].decode()
    assert f">{X}|RIGHT|CLIP|sampleA|" in fa and f">{X}|RIGHT|MATE|sampleB|" in fa

    # TPRT mode: no pooled fragment gate -- Y stays (supported=0 in the evidence TSV); X's
    # combined clip upgraded to the indel-aware consensus carrying the beyond-poly-A sequence.
    tprt = run_combine(str(d2), build_fixture(str(d2)), monkeypatch, {"indel_aware_consensus": True})
    recs = _fastq_records(tprt["combined.txt.gz"])
    assert {k.rsplit(":", 1)[0] for k in recs} == {X, Y}
    assert BEYOND.lower() in recs[f"{X}:R"]
    assert f">{Y}" in tprt["genotyping.txt.gz"].decode()
    assert f">{X}" in tprt["genotyping.txt.gz"].decode()
    rows = _tsv(tprt["insertions.evidence.tsv.gz"])
    assert {(r["insertion_id"], r["supported"]) for r in rows if r["side"] == "RIGHT"} == {(X, "1"), (Y, "0")}
