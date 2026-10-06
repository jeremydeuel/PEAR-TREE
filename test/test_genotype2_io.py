"""tools/genotype2_io: format detection (per-colony files, patient matrix), buckets and the
joint table, for legacy call strings and peartree-genotype2 numeric output.

Run:  pytest test/test_genotype2_io.py
"""
import gzip
import os
import sys

import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, REPO)
from tools import genotype2_io as GIO  # noqa: E402


def _gz(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return str(path)


V2_ROWS = ("\t".join(GIO.V2_HEADER) + "\n"
           "chr1:100-120\tTSD\tok\t30\t12\t10\t3\t0\t0\t6\t6\t0.52\t0.0000\t0.9990\t0.0010\t300\t0\t200\t99\t500\t400\n"
           "chr1:200-220\tTSD\tok\t30\t0\t25\t1\t0\t0\t0\t0\t0.00\t0.9500\t0.0500\t0.0000\t0\t13\t300\t13\t0\t900\n"
           "chr1:300-320\tTSD\tok\t4\t1\t2\t1\t0\t0\t1\t0\t0.30\t0.4000\t0.5000\t0.1000\t3\t0\t9\t3\t20\t40\n"
           "chr1:400-420\tTSD\thigh_coverage\t801" + "\t" * 17 + "\n")
LEGACY_ROWS = ("insertion\tgenotype\tscore_genotype\tscore_alternative\tcoverage\tn_alt\tn_ref\tn_art\n"
               "chr1:100-120\theterozygous\t10\t20\t22\t12\t10\t0\n"
               "chr1:200-220\twild-type\t10\t20\t25\t0\t25\t0\n"
               "chr1:300-320\twild-type?\t1\t2\t3\t1\t2\t0\n")


def test_colony_format_and_rows(tmp_path):
    v2 = _gz(tmp_path / "S1.txt.gz", V2_ROWS)
    lg = _gz(tmp_path / "S2.txt.gz", LEGACY_ROWS)
    assert GIO.colony_format(v2) == GIO.FMT_V2 and GIO.colony_format(lg) == GIO.FMT_LEGACY
    fmt, rows = GIO.read_colony_rows(v2)
    assert [GIO.v2_row_bucket(r) for r in rows.values()] == ["present", "absent", "ambiguous", "high_coverage"]
    fmt, rows = GIO.read_colony_rows(lg)
    assert [GIO.colony_row_bucket(r, fmt) for r in rows.values()] == ["present", "absent", "ambiguous"]
    df = GIO.read_colony_df(v2)
    assert df.attrs["format"] == "v2" and df.loc["chr1:400-420", "n_alt"] != df.loc["chr1:400-420", "n_alt"]  # NaN
    assert GIO.read_colony_df(lg).index.name == "locus"
    bad = _gz(tmp_path / "x.txt.gz", "a\tb\n")
    with pytest.raises(ValueError):
        GIO.colony_format(bad)
    assert set(GIO.list_colony_files(str(tmp_path))) == {"S1", "S2", "x"}


def test_matrix_format_buckets_and_joint(tmp_path):
    num = _gz(tmp_path / "P.genotypes.csv.gz", ";S1;S2;S3\nchr1:100-120;0.9500;;0.0500\nchr1:200-220;0.5000;0.1000;1.0000\n")
    calls = _gz(tmp_path / "Q.genotypes.csv.gz", "insertion;S1;S2\nchr1:100-120;heterozygous;wild-type\n")
    empty = _gz(tmp_path / "R.joint_matrix.csv.gz", ";S1;S2\n")
    assert GIO.matrix_format(num) == "numeric" and GIO.matrix_format(calls) == "calls"
    assert GIO.matrix_format(empty) == "numeric"
    fmt, cols, rows = GIO.read_matrix(num)
    assert cols == ["S1", "S2", "S3"]
    assert [GIO.matrix_cell_bucket(c, fmt) for c in rows["chr1:100-120"]] == ["present", "no_data", "absent"]
    assert [GIO.p_carrier_bucket(c) for c in rows["chr1:200-220"]] == ["ambiguous", "absent", "present"]
    m = GIO.read_matrix_df(num)
    assert m.attrs["format"] == "numeric" and m.isna().sum().sum() == 1
    assert GIO.read_matrix(str(tmp_path / "missing.csv.gz")) == (None, [], None)
    assert GIO.joint_tsv_for(num) is None
    (tmp_path / "P.joint.tsv").write_text("locus\tbest\tcarriers\tn_carriers\tlog10_bf_tree\n"
                                          "chr1:100-120\tN2\tS1,S3\t2\t1.2\nchr1:200-220\tS1\tS1\t1\t0.5\n"
                                          "chr1:300-320\tNOISE\t\t0\t-2\n")
    assert GIO.joint_tsv_for(num) == str(tmp_path / "P.joint.tsv")
    j = GIO.read_joint(GIO.joint_tsv_for(num))
    assert [GIO.joint_class(r) for r in j.values()] == ["clade", "private", "NOISE"]
    assert GIO.read_joint_df(GIO.joint_tsv_for(num)).loc["chr1:100-120", "best"] == "N2"


def test_legacy_bucket_vocabularies():
    assert GIO.legacy_bucket("insertion") == "present"
    assert GIO.legacy_bucket("wild-type?") == "ambiguous"
    assert GIO.legacy_bucket("wild-type?", wildtype={"wild-type", "wild-type?"}) == "absent"
    assert GIO.legacy_bucket("no-coverage") == "no_data" and GIO.legacy_bucket("artefact") == "ambiguous"
    assert GIO.P_ABSENT_MATRIX == pytest.approx(1 - GIO.P_CARRIER)
