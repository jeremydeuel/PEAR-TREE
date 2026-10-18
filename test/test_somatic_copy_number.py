"""cluster/somatic_table.py copy-number rule: a germline insertion on a homolog lost in a clade."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cluster"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from somatic_table import cn_masked, cn_region, read_centromeres, read_copy_number  # noqa: E402

CEN = read_centromeres(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "resources", "centromeres.tsv"),
                       "GRCh38")
# ((a,b),(c,d)) and ((e,f),(g,h)): A = (a,b,c,d), B = (e,f,g,h)
TREE = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"
TIPS = list("abcdefgh")


def events(rows):
    d = tempfile.mkdtemp()
    t, c = os.path.join(d, "t.tree"), os.path.join(d, "copy_number.tsv")
    with open(t, "w") as fh:
        fh.write(TREE)
    with open(c, "w") as fh:
        fh.write("# comment\nregion\tclade\tevent\tsource\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return read_copy_number(c, t, TIPS, CEN)


class Region(unittest.TestCase):
    def test_forms(self):
        self.assertEqual(cn_region("7:100-200", CEN), ("7", 100, 200))
        cs, ce, ln = CEN["7"]
        self.assertEqual(cn_region("chr7", CEN), ("7", 0, ln))
        self.assertEqual(cn_region("7p", CEN), ("7", 0, cs))
        self.assertEqual(cn_region("7q", CEN), ("7", ce, ln))
        with self.assertRaises(ValueError):
            cn_region("7r", CEN)


class Masked(unittest.TestCase):
    def test_complement_of_lost_clade_is_masked(self):
        ev = events([("7", "c,d", "loss", "paper")])            # clade (c,d)
        carriers = list("abefgh")
        self.assertIn("loss 7", cn_masked("7:1000-1015", carriers, TIPS, ev))
        self.assertIn("loss 7", cn_masked("chr7:1000-1015", carriers, TIPS, ev))

    def test_other_chromosome_or_arm_not_masked(self):
        ev = events([("7q", "c,d", "loss", "")])
        self.assertEqual(cn_masked("8:1000-1015", list("abefgh"), TIPS, ev), "")
        self.assertEqual(cn_masked("7:1000-1015", list("abefgh"), TIPS, ev), "")   # 7p is intact

    def test_mrca_expands_clade(self):
        ev = events([("7", "a,d", "LOH", "")])                    # MRCA(a,d) = (a,b,c,d)
        self.assertEqual(ev[0][3], frozenset("abcd"))
        self.assertTrue(cn_masked("7:1-2", list("efgh"), TIPS, ev))

    def test_real_somatic_clade_not_masked(self):
        ev = events([("7", "c,d", "loss", "")])
        self.assertEqual(cn_masked("7:1-2", list("ef"), TIPS, ev), "")       # somatic in (e,f)
        self.assertEqual(cn_masked("7:1-2", list("abcdefgh"), TIPS, ev), "")  # the other homolog was lost

    def test_carrier_inside_lost_clade_breaks_pattern(self):
        ev = events([("7", "a,d", "loss", "")])
        self.assertEqual(cn_masked("7:1-2", list("aefgh"), TIPS, ev), "")    # 1 of 4 lost tips carries

    def test_two_clades_lose_different_homologs(self):
        ev = events([("7", "a,b", "loss", ""), ("7", "g,h", "loss", "")])
        self.assertTrue(cn_masked("7:1-2", list("cdefgh"), TIPS, ev))       # (g,h) lost the other homolog
        self.assertTrue(cn_masked("7:1-2", list("cdef"), TIPS, ev))         # both lost this one

    def test_tolerance_is_proportional(self):
        ev = events([("7", "a,f", "loss", "")])                   # MRCA(a,f) = the root
        self.assertEqual(ev[0][3], frozenset(TIPS))
        ev = events([("7", "a,c", "loss", "")])                   # (a,b,c,d) lost -> expect (e,f,g,h)
        self.assertEqual(cn_masked("7:1-2", list("efg"), TIPS, ev), "")     # 1 of 4 missing: < 20 tips = exact
        many = [f"t{i}" for i in range(40)]
        ev2 = [("7", 0, 10 ** 9, frozenset(many[:20]), "loss 7")]
        self.assertTrue(cn_masked("7:1-2", many[20:39], many, ev2))          # 1 of 20 missing = 5 %
        self.assertEqual(cn_masked("7:1-2", many[20:38], many, ev2), "")

    def test_gain_ignored(self):
        self.assertEqual(events([("7", "c,d", "gain", "")]), [])

    def test_unknown_tip_rejected(self):
        with self.assertRaises(ValueError):
            events([("7", "c,zz", "loss", "")])


if __name__ == "__main__":
    unittest.main()
