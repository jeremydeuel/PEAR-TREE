"""Structure checks for the insertion-type simulators (test/simlib, val1 --types, fullstack).

Runs on the synthetic stand-in library (no genome needed). Set PEARTREE_HS1_2BIT +
PEARTREE_HS1_RMSK to additionally exercise the real hs1 extraction.

    venv/bin/python -m pytest test/test_simulate_types.py -q
"""
import os
import random
import statistics
import subprocess
import sys

import pysam
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from simlib import library as L  # noqa: E402
from simlib import models as M  # noqa: E402
from simlib.reads import FragmentSampler, Hap  # noqa: E402
from simlib.seqs import en_mismatches, revcomp, rnd_seq  # noqa: E402
from simlib.val1render import PairWriter, prepare_event, render_sample  # noqa: E402


@pytest.fixture(scope="module")
def lib():
    hs1, rmsk = os.environ.get("PEARTREE_HS1_2BIT"), os.environ.get("PEARTREE_HS1_RMSK")
    if hs1 and rmsk and os.path.exists(hs1) and os.path.exists(rmsk):
        return L.load_library(None, hs1, rmsk)
    return L.RteLibrary.synthetic()


@pytest.fixture(scope="module")
def ctx(lib):
    rng = random.Random(7)
    genes = [g for g in (L.gene_from_sequence(f"g{i}", "syn", rnd_seq(rng, 30000), 0, rng)
                         for i in range(6)) if g]
    return M.Ctx(lib, genes)


def apply(key, ctx, strand, seed=1, N=12000, nick=6000):
    rng = random.Random(seed)
    ev = M.build_event(key, rng, ctx)
    ref = rnd_seq(rng, N)
    if ev.plant_en:
        from simlib.seqs import degenerate_en_motif
        ref = M.plant_en_motif(ref, nick, strand, degenerate_en_motif(rng, ev.en_mm))
    return ev, ref, M.apply_event(ref, nick, strand, ev)


@pytest.mark.parametrize("strand", "+-")
@pytest.mark.parametrize("key", ["L1_FULL", "L1_TRUNC", "L1_INV", "ALU_YA5", "SVA_F", "POLYA_ONLY",
                                 "L1_TD3P", "ORPHAN_TD3P", "PSEUDOGENE"])
def test_tsd_present_twice(key, strand, ctx):
    ev, ref, alt = apply(key, ctx, strand)
    t = ev.target_len
    assert ev.target == "TSD" and 4 <= t <= 25
    assert alt.right - alt.left == t
    tsd = ref[alt.left:alt.right]
    assert tsd == alt.tsd_seq and len(tsd) == t
    x0, x1 = alt.x_span
    assert alt.alt[x0 - t:x0] == tsd and alt.alt[x1:x1 + t] == tsd       # duplicated
    assert alt.alt[:x0] == ref[:alt.right] and alt.alt[x1:] == ref[alt.left:]
    # inserted sequence = element-sense insert, reverse-complemented for '-' insertions
    expect = ev.ins if strand == "+" else revcomp(ev.ins)
    assert alt.alt[x0:x1] == expect


@pytest.mark.parametrize("strand", "+-")
def test_polya_side_and_en_motif(strand, ctx):
    for seed in range(20):
        ev, ref, alt = apply("L1_TRUNC", ctx, strand, seed=seed)
        x0, x1 = alt.x_span
        if strand == "+":
            assert alt.polya_side == "LEFT" and alt.alt[x1 - 10:x1].count("A") >= 9
        else:
            assert alt.polya_side == "RIGHT" and alt.alt[x0:x0 + 10].count("T") >= 9
        assert alt.en_mm == ev.en_mm                       # the planted degenerate motif
        assert en_mismatches(alt.en_motif[1:]) == ev.en_mm


def test_tsd_length_distribution(ctx):
    rng = random.Random(3)
    ts = [M.build_event("L1_TRUNC", rng, ctx).target_len for _ in range(400)]
    assert min(ts) >= 4 and max(ts) <= 25
    assert 13 <= statistics.median(ts) <= 17


def test_twin_priming_constraints(ctx, lib):
    rng = random.Random(11)
    juncs, ratios = {}, []
    for _ in range(300):
        ev = M.build_event("L1_INV", rng, ctx)
        el = next(e for e in lib.elements["L1"] if e.id == ev.element_id)
        i = ev.info
        assert i["inv_breakpoint"] >= 590                   # never in the first 590 bp
        a, b, c = i["inv_a"], i["inv_b"], i["inv_breakpoint"]
        inv_part = ev.ins[:b - a]
        assert inv_part == revcomp(el.seq[a:b])
        assert ev.ins[b - a:b - a + (len(el.seq) - c)] == el.seq[c:]
        if i["inv_junction"] == "deletion":
            assert b < c
        elif i["inv_junction"] == "duplication":
            assert b > c
        else:
            assert b == c
        juncs[i["inv_junction"]] = juncs.get(i["inv_junction"], 0) + 1
        ratios.append(i["fwd_len"] / i["inv_len"])
    assert 0.55 <= juncs["deletion"] / 300 <= 0.77
    assert 1.6 <= statistics.median(ratios) <= 3.2


def test_twin_priming_switch(ctx, lib):
    rng = random.Random(5)
    for _ in range(50):
        ev = M.build_event("L1_INV_SWITCH", rng, ctx)
        assert ev.structure == "INVERTED_5P_SWITCH" and ev.type_id == 13
        labels = [p[0] for p in ev.parts]
        assert labels[:3] == ["L1_SWITCH", "L1_INV", "L1"]
        el = next(e for e in lib.elements["L1"] if e.id == ev.element_id)
        a2, b2 = ev.info["switch_a"], ev.info["switch_b"]
        assert b2 <= ev.info["inv_a"] + 40 and ev.ins.startswith(el.seq[a2:b2])


def test_transduction_tag_is_real_flank(ctx, lib):
    rng = random.Random(9)
    for key in ("L1_TD3P", "SVA_TD3P", "ORPHAN_TD3P"):
        for _ in range(30):
            ev = M.build_event(key, rng, ctx)
            klass = "SVA" if key.startswith("SVA") else "L1"
            src = next(s for s in lib.sources[klass] if s.id == ev.source_id)
            td = next(p for p in ev.parts if p[0] == "TD3P")
            tag = ev.ins[td[1]:td[2]]
            s0, e0 = map(int, td[3].split(":")[1].split("-"))
            assert tag == src.flank3[s0:e0]                  # exactly the real downstream flank
            assert any(abs(e0 - x) <= 3 for x in src.endpoints) or e0 == len(src.flank3)
            assert f"TD3P_SOURCE={src.id}" in ev.tags and "TD3P" in ev.tags
            assert ev.ins[td[2]:].count("A") >= 0.9 * (len(ev.ins) - td[2])   # poly-A after tag
            if key != "ORPHAN_TD3P":       # partnered: the 5' part is the SAME source element
                assert ev.element_id == src.id
            else:
                assert ev.element == "ORPHAN_TD"


def test_transduction_endpoints_fixed_per_source(lib):
    for s in lib.sources["L1"]:
        assert 1 <= len(s.endpoints) <= 3
        assert s.endpoints == L._poly_a_signal_endpoints(s.id, s.flank3)   # deterministic


def test_sva_td5p(ctx, lib):
    rng = random.Random(2)
    ev = M.build_event("SVA_TD5P", rng, ctx)
    src = next(s for s in lib.sources["SVA"] if s.id == ev.source_id)
    n5 = ev.info["td5_len"]
    assert ev.ins[:n5] == src.flank5[-n5:] and ev.ins[n5:].startswith(src.element_seq[:50])
    assert "TD5P" in ev.tags


def test_pseudogene_has_exon_junction_decoy_not(ctx):
    rng = random.Random(4)
    for _ in range(20):
        ev = M.build_event("PSEUDOGENE", rng, ctx)
        g = next(x for x in ctx.genes if x.id == ev.info["gene"])
        exons = [p for p in ev.parts if p[0].startswith("EXON")]
        assert len(exons) >= 2 and "EXON_JUNCTION" in ev.tags and ev.element == "PSEUDOGENE"
        for (l1, a1, b1, _), (l2, a2, b2, _) in zip(exons, exons[1:]):
            i1, i2 = int(l1[4:]) - 1, int(l2[4:]) - 1
            assert i2 == i1 + 1 and b1 == a2
            junction = g.exon_seqs[i1][-15:] + g.exon_seqs[i2][:15]
            assert junction in ev.ins and junction not in g.premrna   # spliced, not genomic
        d = M.build_event("PSEUDOGENE_DECOY", rng, ctx)
        gd = next(x for x in ctx.genes if x.id == d.info["gene"])
        body = d.ins[:d.parts[-1][1]]
        assert body in gd.premrna                       # contiguous genomic: NO exon-exon junction
        assert "EXON_JUNCTION" not in d.tags and d.element != "PSEUDOGENE"


def test_tsd_deletion(ctx):
    for strand in "+-":
        ev, ref, alt = apply("L1_TSD_DELETION", ctx, strand, seed=8)
        d = ev.target_len
        assert ev.target == "DEL" and "TSD_DELETION" in ev.tags
        assert alt.left - alt.right == d and alt.tsd_seq == ""
        x0, x1 = alt.x_span
        assert alt.alt[:x0] == ref[:alt.right] and alt.alt[x1:] == ref[alt.left:]


def test_l1_mediated_deletion_microhomology(ctx):
    for strand in "+-":
        ev, ref, alt = apply("L1_MED_DELETION", ctx, strand, seed=6, N=60000, nick=40000)
        assert ev.target == "L1DEL" and 1 <= ev.mh <= 5 and "L1_MED_DELETION" in ev.tags
        assert alt.left - alt.right == ev.target_len           # deleted segment, no TSD
        assert len(alt.mh_seq) == ev.mh
        x0, x1 = alt.x_span
        if strand == "+":
            # the element's first mh bases are supplied by the reference (microhomology)
            assert alt.alt[x0 - ev.mh:x0] == alt.mh_seq
            assert alt.alt[x0:x1] == ev.ins[ev.mh:]


def test_l1_mediated_duplication(ctx):
    ev, ref, alt = apply("L1_MED_DUPLICATION", ctx, "+", seed=6, N=30000, nick=12000)
    assert ev.target == "TSD" and ev.target_len >= 50 and "L1_MED_DUPLICATION" in ev.tags
    x0, x1 = alt.x_span
    t = ev.target_len
    assert alt.alt[x0 - t:x0] == alt.alt[x1:x1 + t]


def test_en_independent(ctx):
    rng = random.Random(1)
    for _ in range(20):
        ev = M.build_event("EN_INDEPENDENT", rng, ctx)
        assert ev.polya_len == 0 and not ev.plant_en and ev.target in ("NONE", "DEL")
        assert ev.ins[-12:].count("A") < 10
        assert "EN_INDEPENDENT" in ev.tags and ev.info["l1_b"] < len(
            next(e for e in ctx.lib.elements["L1"] if e.id == ev.element_id).seq) - 20


@pytest.mark.parametrize("strand", "+-")
def test_templated_local_and_foldback(ctx, strand):
    for seed in range(10):
        ev, ref, alt = apply("TEMPLATED_LOCAL", ctx, strand, seed=seed)
        lab, a, b, info = alt.parts[0]
        assert lab == "TEMPLATED" and b - a == ev.info["templ_len"] < 250
        ref_p = ref if strand == "+" else revcomp(ref)
        r0, r1 = map(int, info[4:-1].split("-"))
        x5 = (alt.right if strand == "+" else len(ref) - alt.left)   # 5' junction, '+' frame
        assert min(abs(r0 - x5), abs(r1 - x5)) <= 15                 # copied from <= 15 bp of the site
        seg = ref_p[r0:r1] if info.endswith("+") else revcomp(ref_p[r0:r1])
        assert alt.x_seq.startswith(seg)
        ev, ref, alt = apply("FOLDBACK_INVDUP_5P", ctx, strand, seed=seed)
        ref_p = ref if strand == "+" else revcomp(ref)
        x5 = alt.right if strand == "+" else len(ref) - alt.left
        f, g = ev.info["foldback_len"], ev.info["foldback_gap"]
        assert alt.x_seq.startswith(revcomp(ref_p[x5 - g - f:x5 - g]))   # inverted dup of upstream


def test_premrna_coinsert(ctx):
    ev, ref, alt = apply("PREMRNA_COINSERT", ctx, "+", seed=3)
    lab, a, b, info = alt.parts[0]
    assert lab == "PREMRNA" and "PREMRNA_COINSERT" in ev.tags
    r0, r1 = map(int, info[4:-1].split("-"))
    assert alt.x_seq.startswith(ref[r0:r1]) and abs(r0 - alt.right) >= 300


def _render(ctx, key, strand, pcr=0.0, depth=30, vaf=0.5, seed=1, jitter=1.0):
    rng = random.Random(seed)
    ev = M.build_event(key, rng, ctx)
    prep = prepare_event(random.Random(seed + 1), ev, strand, ctx.lib)
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "1", "LN": 10 ** 7}]})
    w = PairWriter(hdr, 0, "1", 100000, 151)
    s = FragmentSampler(rng, jitter_scale=jitter, pcr_dup_frac=pcr)
    counts = render_sample(rng, prep, w, s, depth, vaf, True, "t")
    return ev, prep, w.records, counts


def test_pcr_duplicates_unflagged_identical_outer_coords(ctx):
    ev, prep, recs, counts = _render(ctx, "L1_TRUNC", "+", pcr=1.0)
    assert not any(r.flag & 0x400 for r in recs)
    groups = {}
    for r in recs:
        if r.is_supplementary or r.is_unmapped:
            continue
        base = r.query_name.split("_d")[0]
        groups.setdefault((base, r.is_read1), []).append(r)
    dup_groups = [g for g in groups.values() if len(g) > 1]
    assert dup_groups
    for g in dup_groups:
        # identical outer (unclipped 5') coordinates and mate placement
        outers = {(r.reference_end + (r.query_length - r.query_alignment_end) if r.is_reverse
                   else r.reference_start - r.query_alignment_start, r.next_reference_start) for r in g}
        # same template ends; homopolymer jitter near a read end can move the unclipped
        # coordinate by a few bp (hence the SPEC's +-2 tolerance), never more
        xs = [o[0] for o in outers]
        assert max(xs) - min(xs) <= 4
    assert counts["L_reads"] > counts["L_frags"]       # reads incl. duplicates > fragments


def test_polya_jitter_per_read(ctx):
    seq = rnd_seq(random.Random(1), 300) + "A" * 60 + rnd_seq(random.Random(2), 300)
    hap = Hap(seq)
    n0 = next(r1 - r0 for r0, r1, b in hap.runs if r1 - r0 >= 60)      # true run incl. flank A's
    s = FragmentSampler(random.Random(3), jitter_scale=1.0, burst_scale=0.0, error_rate=0.0)
    lens = []
    for _ in range(300):
        (e, sq, q, f), _ = s.pair(hap, 250, 600)
        run = max((len(x) for x in sq.split("C") + sq.split("G") + sq.split("T") if set(x) == {"A"}),
                  default=0)
        full = sq.find("A" * 30) > 0 and sq[sq.find("A" * 30):].lstrip("A")
        if full:
            lens.append(run)
    assert lens and statistics.pstdev(lens) > 1.0       # +-a few bases for a 60-mer
    assert n0 - 2 <= statistics.mean(lens) <= n0 + 2                    # unbiased


def test_mates_and_clips(ctx):
    ev, prep, recs, counts = _render(ctx, "L1_TRUNC", "+", depth=40)
    clipped = [r for r in recs if not r.is_unmapped and "S" in r.cigarstring and not r.is_supplementary]
    assert clipped
    names = {r.query_name for r in recs}
    for n in list(names)[:50]:
        rr = [r for r in recs if r.query_name == n and not r.is_supplementary]
        assert len(rr) == 2 and {r.is_read1 for r in rr} == {True, False}
    # some mates fall inside the element: unmapped, placed at the mapped mate
    assert any(r.is_unmapped and r.reference_start >= 0 for r in recs)
    assert counts["L_frags"] >= 2 and counts["R_frags"] >= 2


def test_val1_cli_multisample(tmp_path, lib):
    out = tmp_path / "s.bam"
    truth = tmp_path / "t.tsv"
    cmd = [sys.executable, os.path.join(HERE, "val1", "simulate.py"), "--out-bam", str(out),
           "--out-truth", str(truth), "--types", "L1_TRUNC,ORPHAN_TD3P,ART_LIGATION_PCR",
           "--n-per-type", "2", "--samples", "2", "--n-l1", "1", "--n-erv", "0",
           "--n-l1-variant", "0", "--n-alu", "0", "--n-sva", "0", "--n-erv-variant", "0",
           "--n-hervk113", "0", "--n-hervk117", "0"]
    subprocess.run(cmd, check=True, capture_output=True)
    for i in (1, 2):
        p = tmp_path / f"s.S{i}.bam"
        assert p.exists() and pysam.AlignmentFile(str(p)).mapped > 0
    rows = [l.rstrip("\n").split("\t") for l in open(truth)]
    hdr = rows[0]
    assert hdr[:8] == ["contig", "left", "right", "class", "tsd", "alt_reads", "ref_reads", "vaf"]
    for c in ("role", "type_id", "element", "structure", "tags", "samples"):
        assert c in hdr
    roles = {r[hdr.index("variant")]: r[hdr.index("role")] for r in rows[1:]}
    assert roles["ART_LIGATION_PCR"] == "ARTEFACT" and roles["L1_TRUNC"] == "TP"
