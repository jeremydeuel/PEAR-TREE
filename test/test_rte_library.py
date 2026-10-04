# PEAR-TREE - RTE reference library sanity tests (resources/rte_library/, tools/rte_library/)
#
# Locks in the orientation / content contract of the committed library that annotate relies on:
#   * the L1HS consensus is sense (starts with the 5'UTR GGGGGAGGAGCCAAGATGGCCGAATAGGAACAGCTCCGG,
#     ends at the poly-A signal) and carries ORF1 1,017 nt / ORF2 3,828 nt landmarks
#   * the Alu consensus carries the Pol III A box and B box
#   * every intact L1 is sense (5' end matches the consensus 5'UTR)
#   * 3' flanks are downstream of the source on the correct strand (+ and -), including
#     insertion-point (non-reference) sources, and the flank-window arithmetic is right on a toy
#     genome
#   * the Ta/pre-Ta typer and the novel-source tiering behave as documented
# No genome is needed: everything is checked against the committed files.
#
# Run with:  python test/test_rte_library.py     (from the repo root)
#        or:  pytest test/test_rte_library.py

import csv
import gzip
import hashlib
import os
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LIB = os.path.join(REPO, "resources", "rte_library")
sys.path.insert(0, os.path.join(REPO, "tools", "rte_library"))

from common import read_fasta, revcomp  # noqa: E402

try:
    import edlib  # noqa: F401
    HAVE_EDLIB = True
except ImportError:
    HAVE_EDLIB = False

L1_5UTR = "GGGGGAGGAGCCAAGATGGCCGAATAGGAACAGCTCCGG"


def _tsv(name):
    with open(os.path.join(LIB, name)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _cons():
    return {n: s for n, _, s in read_fasta(os.path.join(LIB, "consensus.fa"))}


def _landmarks():
    out = {}
    for r in _tsv("consensus_landmarks.tsv"):
        out.setdefault(r["consensus"], {}).setdefault(r["feature"], []).append(
            (int(r["start"]), int(r["end"])))
    return out


def _flank_headers():
    hdr = {}
    with gzip.open(os.path.join(LIB, "flanks_3p.fa.gz"), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name, desc = line[1:].rstrip("\n").split(" ", 1)
                hdr[name] = desc
    return hdr


def _parse_flank_desc(desc):
    m = re.match(r"(hs1|hg38):(\S+):(\d+)-(\d+)\(([+-])\)", desc)
    return m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)


def _end_matches(s, c, n=20, shifts=3, max_mm=5):
    """3' end of s equals the consensus 3' end within `max_mm` mismatches at a shift of <= 3 bp
    (sense check that tolerates point mutations and a ragged poly-A boundary)."""
    for d in range(-shifts, shifts + 1):
        a = s[len(s) - n + d:len(s) + d] if d <= 0 else s[len(s) - n:]
        b = c[-n:] if d <= 0 else c[len(c) - n - d:len(c) - d]
        if len(a) == len(b) == n and sum(x != y for x, y in zip(a, b)) <= max_mm:
            return True
    return False


# ------------------------------------------------------------------------------------ files
def test_manifest_matches_files():
    rows = _tsv("manifest.tsv")
    assert rows, "empty manifest"
    total = 0
    for r in rows:
        p = os.path.join(LIB, r["file"])
        assert os.path.exists(p), r["file"]
        assert hashlib.md5(open(p, "rb").read()).hexdigest() == r["md5"], "md5 drift: " + r["file"]
        total += int(r["bytes"])
    assert total < 10 * 1024 * 1024, "library must stay < 10 MB (is %d)" % total
    for f in ("l1_intact.fa", "alu_y_intact.fa", "sva_intact.fa", "consensus.fa", "active.tsv",
              "transduction_sources.tsv", "flanks_3p.fa.gz", "flanks_5p_sva.fa.gz",
              "consensus_landmarks.tsv", "README.md"):
        assert os.path.exists(os.path.join(LIB, f)), f


# ------------------------------------------------------------------------------------ consensus
def test_l1hs_consensus_is_sense_and_has_landmarks():
    c = _cons()["L1HS"]
    assert 5900 <= len(c) <= 6100
    # starts with the 5'UTR (allow <= 2 mismatches)
    mm = sum(1 for a, b in zip(c[:len(L1_5UTR)], L1_5UTR) if a != b)
    assert mm <= 2, c[:60]
    # ends at the canonical pA signal: AATAAA whose last bases are the poly-A itself
    assert (c[-20:] + "AAA").endswith("AATAAA"), c[-30:]
    lm = _landmarks()["L1HS"]
    o1, o2 = lm["ORF1"][0], lm["ORF2"][0]
    assert o1[1] - o1[0] + 1 == 1017 and o2[1] - o2[0] + 1 == 3828
    assert c[o1[0] - 1:o1[0] + 2] == "ATG" and c[o2[0] - 1:o2[0] + 2] == "ATG"
    assert c[o2[1] - 3:o2[1]] in ("TAA", "TAG", "TGA")
    assert lm["5UTR"][0] == (1, o1[0] - 1)
    assert lm["POLYA_SIGNAL"][0][1] == len(c)
    assert o2[1] < lm["3UTR"][0][1] == len(c)


def test_alu_consensus_has_a_and_b_box():
    cons = _cons()
    lm = _landmarks()
    for name in ("ALU_Y", "ALU_YA5", "ALU_YB8"):
        c = cons[name]
        assert 270 <= len(c) <= 300, (name, len(c))
        assert re.search("GGCTCACGCC", c[:40]), name                    # A box
        assert re.search("GAG[AT]TCGAGAC", c[50:110]), name             # B box
        assert not c.endswith("AAAAA"), "consensus must be poly-A stripped"
        assert "A_BOX" in lm[name] and "B_BOX" in lm[name]
        assert lm[name]["A_BOX"][0][0] < lm[name]["B_BOX"][0][0] < lm[name]["A_RICH_LINKER"][0][0]


def test_sva_consensus_starts_with_hexamer():
    cons = _cons()
    for name in ("SVA_E", "SVA_F"):
        assert "CCCTCT" in cons[name][:60], name
        assert "AATAAA" in cons[name][-80:] + "AAA", name     # may span the poly-A boundary


def test_ta_typer_on_consensus():
    if not HAVE_EDLIB:
        return
    from build import ta_status
    cons = _cons()
    assert ta_status(cons["L1HS"]) == "Ta"
    assert ta_status(cons["L1PA2"]) == "nonTa"
    assert ta_status(cons["L1PA3"]) == "nonTa"


# ------------------------------------------------------------------------------------ intact sets
def test_l1_intact_sense_and_annotated():
    rows = {r["id"]: r for r in _tsv("l1_intact.tsv")}
    seqs = read_fasta(os.path.join(LIB, "l1_intact.fa"))
    assert len(rows) == len(seqs) == 146
    c = _cons()["L1HS"]
    for name, _, s in seqs:
        r = rows[name]
        # L1Base FLI-L1 = both ORFs intact by L1Base's criteria; for the L1HS class we require
        # full-length ORFs (3 L1PA2 copies carry small in-frame ORF2 deletions)
        if r["subfamily_call"] == "L1HS":
            assert int(r["orf1_nt"]) >= 1014 and int(r["orf2_nt"]) >= 3822, name
        else:
            assert int(r["orf1_nt"]) >= 1000 and int(r["orf2_nt"]) >= 3400, name
        assert r["ta_status"] in ("Ta", "preTa", "nonTa", "ambiguous", "unknown")
        assert len(s) == int(r["length"])
        # sense: the 3' end carries the L1 3'UTR end ("...CTTAGAGTATAAT"); trimming can be
        # off by a base or two at the poly-A boundary
        assert _end_matches(s, c), (name, s[-20:])
    if HAVE_EDLIB:
        import edlib
        # 5' end = L1HS 5'UTR (L1PA2 5'UTRs diverge; their sense is proven by the forward-
        # strand ORFs above)
        for name, _, s in seqs:
            if rows[name]["subfamily_call"] != "L1HS" or int(rows[name]["l1hs_cons_start"]) > 50:
                continue                     # L1Base FLI-L1 = intact ORFs; some lack 5'UTR
            fwd = edlib.align(s[:300], c[:600], mode="HW")["editDistance"]
            rev = edlib.align(revcomp(s[:300]), c, mode="HW")["editDistance"]
            assert fwd < 40 and fwd < rev, (name, fwd, rev)
    young = [r for r in rows.values() if r["subfamily_call"] == "L1HS"]
    assert len(young) >= 100


def test_alu_sva_intact_sense():
    for fa, motif in (("alu_y_intact.fa", "GGCTCAC"), ("sva_intact.fa", "CCCTCT")):
        seqs = read_fasta(os.path.join(LIB, fa))
        hits = sum(1 for _, _, s in seqs if motif in s[:80].upper())
        assert hits >= 0.8 * len(seqs), (fa, hits, len(seqs))


def test_active_has_identity():
    rows = _tsv("active.tsv")
    assert len(rows) >= 100
    for r in rows:
        if r["identity_to_consensus"] != ".":
            assert 0.9 <= float(r["identity_to_consensus"]) <= 1.0, r
    assert any(r["tier"] == "hot_source" for r in rows)
    assert any(r["tier"] == "L1HS_Ta_intact" for r in rows)


# ------------------------------------------------------------------------------------ sources
def test_sources_ids_flanks_and_hot_loci():
    rows = _tsv("transduction_sources.tsv")
    hdr = _flank_headers()
    ids = [r["id"] for r in rows]
    assert len(ids) == len(set(ids)), "duplicate source ids"
    for r in rows:
        for nm in r["flank_3p"].split(","):
            assert nm in hdr, nm
        assert r["hotness"] in ("hot", "strong", "active", "none_reported", "candidate")
    # the hottest known source (22q12.1, TTC28 intron) is there, reference, + strand
    hot = [r for r in rows if r["band"] == "22q12.1" and int(r["n_daughters"]) >= 500]
    assert len(hot) == 1 and hot[0]["reference"] == "yes" and hot[0]["strand"] == "+"
    for band in ("Xp22.2", "6p24.1", "1p12", "2q24.1"):
        assert any(r["band"] == band and r["hotness"] == "hot" for r in rows), band


def test_flanks_are_downstream_on_the_correct_strand():
    rows = _tsv("transduction_sources.tsv")
    hdr = _flank_headers()
    seen = {"+": 0, "-": 0, "point+": 0, "point-": 0}
    for r in rows:
        if r["hs1_chrom"] == "." or r["hs1_strand"] not in "+-" or "," in r["flank_3p"]:
            continue
        g, chrom, f0, f1, fst = _parse_flank_desc(hdr[r["flank_3p"]])
        assert g == "hs1" and chrom == r["hs1_chrom"] and fst == r["hs1_strand"]
        s, e = int(r["hs1_start"]), int(r["hs1_end"])
        if r["hs1_status"] == "insertion_point":
            # s == e == 0-based junction offset
            if fst == "+":
                assert f0 == s + 1, r["id"]
                seen["point+"] += 1
            else:
                assert f1 == s, r["id"]
                seen["point-"] += 1
        else:
            if fst == "+":
                assert f0 == e + 1, r["id"]           # starts right after the element's 3' end
                seen["+"] += 1
            else:
                assert f1 == s - 1, r["id"]           # ends right before the element's 5'-most base
                seen["-"] += 1
        assert f1 - f0 + 1 == int(r["flank_3p_len"]) or f0 == 1
    assert min(seen.values()) > 10, seen


def test_flank_window_arithmetic_on_toy_genome():
    from sources import flank_window
    #          0123456789012345678901234567890123456789
    genome = "AAAACCCCGGGGTTTTACGTACGTGGGGCCCCAAAATTTT"
    # element on + at 1-based 9..16 (GGGGTTTT): downstream 3' flank = genome[16:24]
    f0, f1 = flank_window(9, 16, "+", 8, "3p")
    assert genome[f0:f1] == "ACGTACGT"
    # element on - at 1-based 9..16: its 3' end is the left edge -> flank = rc(genome[0:8])
    f0, f1 = flank_window(9, 16, "-", 8, "3p")
    assert (f0, f1) == (0, 8) and revcomp(genome[f0:f1]) == "GGGGTTTT"
    # 5' flank of a + element is upstream
    f0, f1 = flank_window(9, 16, "+", 8, "5p")
    assert genome[f0:f1] == "AAAACCCC"
    # insertion point (junction offset 16): + flank starts there, - flank ends there
    assert flank_window(16, 16, "+", 8, "3p", point=True) == (16, 24)
    assert flank_window(16, 16, "-", 8, "3p", point=True) == (8, 16)


def test_novel_source_tiers():
    from add_source import source_tier
    assert source_tier(6000, 0.992, 1) == ("A", [])
    t, f = source_tier(6000, 0.97, 1)
    assert t == "B" and f
    assert source_tier(6000, 0.97, 2) == ("B", [])
    t, f = source_tier(6000, 0.92, 5)
    assert t == "reject" and f
    t, f = source_tier(4000, 0.995, 5)
    assert f and "length" in f[0]


if __name__ == "__main__":
    fails = 0
    for k, v in sorted(globals().items()):
        if k.startswith("test_") and callable(v):
            try:
                v()
                print("ok   ", k)
            except Exception as e:  # noqa: BLE001
                fails += 1
                print("FAIL ", k, repr(e))
    sys.exit(1 if fails else 0)
