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


def test_polya_end_is_gated(monkeypatch):
    """RIGHT has 2 independent fragments, the poly-A (LEFT, outward poly-T) end only one."""
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
    cfg = {"require_independent_fragments": True, "min_independent_fragments": 2}
    kept, records, failed, stats = apply_evidence([Ins()], ["S1.txt.gz"], cfg)
    assert kept == [] and failed == {"chr1:100-110"}
    left = [r for r in records["chr1:100-110"] if r.side == "LEFT"][0]
    assert left.polya_end == 1 and left.supported == 0
    assert stats["reasons"] == {"LEFT(polyA)": 1}
    # same data, gate off -> kept, still reported
    kept, records, failed, _ = apply_evidence([Ins()], ["S1.txt.gz"], {"min_independent_fragments": 2})
    assert len(kept) == 1 and not failed


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
    # (emulated in Sample.write) truncates at the poly-A, as on real data
    a.add_right(X, 10015, 18, 12, 5)
    f = a.add_right(X, 10015, 20, 14, 9)
    a.add_right(X, 10015, 20, 14, 9)        # PCR dup of f (new qname, same molecule)
    a.add_right(X, 10015, 16, 15, 7)
    b.add_right(X, 10015, 16, 13, 3)
    b.add_right(X, 10015, 19, 12, 11)
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

    # TPRT mode: Y gated out (still in the evidence TSV with supported=0); X's combined
    # clip upgraded to the indel-aware consensus carrying the beyond-poly-A sequence.
    tprt = run_combine(str(d2), build_fixture(str(d2)), monkeypatch,
                       {"require_independent_fragments": True, "indel_aware_consensus": True})
    recs = _fastq_records(tprt["combined.txt.gz"])
    assert {k.rsplit(":", 1)[0] for k in recs} == {X}
    assert BEYOND.lower() in recs[f"{X}:R"]
    assert f">{Y}" not in tprt["genotyping.txt.gz"].decode()
    assert f">{X}" in tprt["genotyping.txt.gz"].decode()
    rows = _tsv(tprt["insertions.evidence.tsv.gz"])
    assert {(r["insertion_id"], r["supported"]) for r in rows if r["side"] == "RIGHT"} == {(X, "1"), (Y, "0")}
