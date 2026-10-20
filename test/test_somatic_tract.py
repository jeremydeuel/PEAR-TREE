"""cluster/somatic_table.py poly-A tract rule: a 'TSD' that is a reference homopolymer (read-through)."""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cluster"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from somatic_table import rte_clip, rte_index, tract_tsd  # noqa: E402
from tools.rte.genome import FastaGenome  # noqa: E402

LIB = rte_index(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "resources", "rte_library"))
ALU = dict((n, s) for n, s, _ in LIB)["ALU_Y"]
rng = random.Random(7)
UP = "".join(rng.choice("ACGT") for _ in range(300))     # unique flanks (no long runs by chance)
DOWN = "".join(rng.choice("ACGT") for _ in range(300))
UP, DOWN = UP.replace("AAAA", "ACGA").replace("TTTT", "TGCT"), DOWN.replace("AAAA", "ACGA").replace("TTTT", "TGCT")
UP, DOWN = UP[:-1] + "C", "G" + DOWN[1:]          # the tract is exactly A18
OFF = 1000
# reference: flank | A18 | flank; the tract is 0-based [1300, 1318)
GENOME = FastaGenome(records={f"chr1:{OFF}-{OFF + 618}": UP + "A" * 18 + DOWN})


def case(lo, hi, lclip, rclip):
    loc = f"chr1:{lo}-{hi}"
    cons = {(loc, "L"): lclip.lower() + "ACGT", (loc, "R"): "ACGT" + rclip.lower()}
    return tract_tsd(loc, cons, {}, GENOME, LIB)


class Tract(unittest.TestCase):
    def test_tsd_is_tract_capped(self):
        v, note = case(1300, 1318, "TGTTCTTCCCCCCC", "AAAAAAAAGGGGGG")
        self.assertEqual(v, "cap")
        self.assertIn("A18", note)

    def test_tolerance_and_one_end_flag(self):
        self.assertEqual(case(1302, 1316, "CCCCCC", "GGGGG")[0], "cap")
        self.assertEqual(case(1300, 1330, "CCCCCC", "GGGGG")[0], "flag")
        self.assertEqual(case(1250, 1270, "CCCCCC", "GGGGG")[0], "")

    def test_rte_clip_exempts(self):
        # a genuine Alu: one clip is the young-Alu 5' head, the other its poly-A
        self.assertEqual(case(1300, 1318, ALU[:120], "AAAAAAAAAAAA")[0], "")
        # an Alu tail that only starts after extra tract bases still counts (the tract run is stripped)
        self.assertEqual(case(1300, 1318, "CCCCCC", "A" * 6 + ALU[150:280])[0], "")

    def test_short_tract_ignored(self):
        g = FastaGenome(records={f"chr1:{OFF}-{OFF + 608}": UP + "A" * 8 + DOWN})
        loc = "chr1:1300-1308"
        self.assertEqual(tract_tsd(loc, {}, {}, g, LIB)[0], "")

    def test_rte_clip(self):
        self.assertTrue(rte_clip(ALU[40:140], LIB).startswith("ALU"))
        mut = "".join(("G" if c == "A" else c) if i % 9 == 0 else c for i, c in enumerate(ALU[40:140]))
        self.assertTrue(rte_clip(mut, LIB))                       # ~90 % identity still RTE
        self.assertEqual(rte_clip(UP[:100], LIB), "")
        self.assertEqual(rte_clip(ALU[40:60], LIB), "")             # < 25 bp


if __name__ == "__main__":
    unittest.main()
