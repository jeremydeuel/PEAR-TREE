# PEAR-TREE - genotyping of the TPRT locus kinds (one-sided, target-site deletion, blunt,
# L1-mediated deletion / duplication): the one-sided contract extension
# (src/genotyping_contract_oneside.py) and combine_genotypes' handling of the new names.
#
# Run with:  pytest test/test_genotype_tprt_loci.py

import gzip
import os
import sys


import pytest

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import genotyping_contract_oneside as gco  # noqa: E402  (imports no config at module level)

NAMES = {
    "chr22:1000-1015": "TSD",
    "chr22:1000-980": "TSD_DELETION",
    "chr22:1000-1000": "BLUNT",
    "chr22:1000-1001": "BLUNT",
    "chr22:5000-1000": "L1_MED_DELETION",
    "chr22:1000-1300": "L1_MED_DUPLICATION",
    "chr22:1000-oneside_1000": "ONE_SIDED",
    "chr22:oneside_2000-2000": "ONE_SIDED",
    "HLA-A*01:01:01:100-112": "TSD",
    "chr22:polyA_100-200": "OTHER",
}


TEST_CONFIG = ('CONFIG = {"combine_genotypes": {"min_wild-types": 1, "min_insertions": 1, '
               '"max_artefact": 2, "max_na": 2, "min_best_score": 0, "min_dispersion": 0}}\n')


@pytest.fixture(scope="module")
def cg(tmp_path_factory):
    """combine_genotypes with a test `config` (src/config.py is a gitignored deploy copy). The
    config is a real file on sys.path, not a sys.modules entry: collect_genotype's Pool workers
    (spawned on macOS) re-import combine_genotypes and inherit only sys.path."""
    pytest.importorskip("pandas")
    d = tmp_path_factory.mktemp("cgconfig")
    (d / "config.py").write_text(TEST_CONFIG)
    saved = sys.modules.pop("config", None)
    sys.path.insert(0, str(d))
    sys.modules.pop("combine_genotypes", None)
    import combine_genotypes
    yield combine_genotypes
    sys.path.remove(str(d))
    sys.modules.pop("config", None)
    sys.modules.pop("combine_genotypes", None)
    if saved is not None:
        sys.modules["config"] = saved


def test_locus_kind(cg):
    for name, kind in NAMES.items():
        assert cg.locus_kind(name) == kind, name


def _write_genotypes(path, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("insertion\tgenotype\tscore_genotype\tscore_alternative\tcoverage\tn_alt\tn_ref\tn_art\n")
        for name, gt, n_alt, n_ref in rows:
            fh.write(f"{name}\t{gt}\t2000\t100\t{n_alt + n_ref}\t{n_alt}\t{n_ref}\t0\n")


def test_collect_genotype_keeps_new_names_verbatim(cg, tmp_path, capsys):
    """One-sided / far-pair / target-site-deletion names pass through combine_genotypes as
    matrix row keys unchanged (annotate joins <patient>.genotypes.csv.gz on them), and the
    clade gates work on them like on any TSD locus."""
    names = ["chr22:1000-oneside_1000", "chr22:oneside_2000-2000", "chr22:5000-1000",
             "chr22:1000-1300", "chr22:3000-2990", "chr22:7000-7015"]
    files = []
    for k in range(4):
        carrier = k < 2                      # colonies 0,1 carry every locus, 2,3 are wild-type
        rows = [(n, "heterozygous" if carrier else "wild-type", 6 if carrier else 0, 6 if carrier else 12)
                for n in names]
        # one locus that is het everywhere -> no wild-type colony -> removed (min_wild-types 1)
        rows.append(("chr22:9000-oneside_9000", "heterozygous", 6, 6))
        p = tmp_path / f"S{k + 1}.txt.gz"
        _write_genotypes(p, rows)
        files.append(str(p))
    out = tmp_path / "P1.genotypes.csv.gz"
    cg.collect_genotype(files, str(out), 1)
    log = capsys.readouterr().out
    assert "per locus kind" in log and "ONE_SIDED" in log
    with gzip.open(out, "rt") as fh:
        lines = [line.rstrip("\n").split(";") for line in fh]
    assert lines[0][1:] == ["S1", "S2", "S3", "S4"]
    kept = {row[0]: row[1:] for row in lines[1:]}
    assert set(kept) == set(names)
    assert kept["chr22:1000-oneside_1000"] == ["heterozygous", "heterozygous", "wild-type", "wild-type"]


def test_kind_summary_counts(cg):
    idx = ["chr22:1000-oneside_1000", "chr22:oneside_5-5", "chr22:5000-1000", "chr22:100-112"]
    s = cg.kind_summary(idx, [True, False, False, True])
    assert s.loc["ONE_SIDED"].tolist() == [2, 1, 1]
    assert s.loc["L1_MED_DELETION"].tolist() == [1, 0, 1]
    assert s.loc["TSD"].tolist() == [1, 1, 0]


# ---------------------------------------------------------------- one-sided contract extension
GENOME = "".join("ACGTTGCAAGCTTCGA"[(i * 7) % 16] for i in range(4000))


def _get_sequence(contig, start, end):
    return GENOME[start:end] if contig == "chr22" and 0 <= start < end else ""


def _score(seqs):
    from collections import Counter
    n = min(len(s) for s in seqs)
    same = sum(1 for i in range(n) if Counter(s[i] for s in seqs).most_common(1)[0][1] > len(seqs) / 2)
    return (same - 2 * (n - same)) / n


def _fq(name, seq):
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n"


def test_one_sided_side_parse():
    assert gco.one_sided_side("chr22:1000-oneside_1000") == ("LEFT", "chr22", 1000)
    assert gco.one_sided_side("chr22:oneside_2000-2000") == ("RIGHT", "chr22", 2000)
    assert gco.one_sided_side("chr22:1000-1015") is None


def test_build_entry_left_and_right():
    clip_l = "ttttttttgcgcgcgcatat"            # lower-case clip, junction-adjacent at its 3' end
    sides = {"L": clip_l + GENOME[1000:1060]}
    text, why = gco.build_entry("chr22:1000-oneside_1000", sides, _get_sequence, 12, _score)
    assert why is None
    assert text == (">chr22:1000-oneside_1000\n@LEFT_INSERTION\n" + clip_l[-12:].upper()
                    + "\n@LEFT_REFERENCE\n" + GENOME[988:1000] + "\n")
    clip_r = "aaaaaaaaaaaaaaaaacgt"
    sides = {"R": GENOME[1940:2000] + clip_r}
    text, why = gco.build_entry("chr22:oneside_2000-2000", sides, _get_sequence, 12, _score)
    assert why is None
    assert "@RIGHT_INSERTION\n" + clip_r[:12].upper() + "\n@RIGHT_REFERENCE\n" + GENOME[2000:2012] in text
    assert "LEFT" not in text
    # a clip identical to the reference is excluded, like combine excludes it
    sides = {"L": GENOME[988:1000].lower() + GENOME[1000:1060]}
    assert gco.build_entry("chr22:1000-oneside_1000", sides, _get_sequence, 12, _score) == \
        (None, "clip similar to reference")
    # unknown contig / long contig names
    assert gco.build_entry("chrUn_KI270302v1:5-oneside_5", {"L": "acgt" + "ACGT"}, _get_sequence, 12, _score)[1] == "contig"


def test_extend_contract_appends_after_unchanged_base(tmp_path):
    base = ">chr22:100-112\n@RIGHT_INSERTION\nAAAA\n@RIGHT_REFERENCE\nCCCC\n"
    contract = tmp_path / "P1.genotyping.txt.gz"
    with gzip.open(contract, "wt") as fh:
        fh.write(base)
    combined = tmp_path / "P1.combined.txt.gz"
    with gzip.open(combined, "wt") as fh:
        fh.write(_fq("chr22:100-112:L", "acgtacgtacgt" + GENOME[100:140]))
        fh.write(_fq("chr22:100-112:R", GENOME[60:112] + "tttt"))
        fh.write(_fq("chr22:1000-oneside_1000:L", "ttttttttgcgcgcgcatat" + GENOME[1000:1060]))
        fh.write(_fq("chr22:oneside_2000-2000:R", GENOME[1940:2000] + "aaaaaaaaaaaaaaaaacgt"))
    out = tmp_path / "P1.genotyping.tprt.txt.gz"
    n, excluded = gco.extend_contract(str(contract), str(combined), str(out), _get_sequence, 12, _score)
    assert (n, excluded) == (2, {})
    with gzip.open(out, "rt") as fh:
        text = fh.read()
    assert text.startswith(base)
    assert [line for line in text.splitlines() if line.startswith(">")] == \
        [">chr22:100-112", ">chr22:1000-oneside_1000", ">chr22:oneside_2000-2000"]
