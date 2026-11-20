import sys, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from run_pipeline import find_orfs, gc_percent, reverse_complement


class TestPipelineFunctions(unittest.TestCase):
    def test_gc_percent(self):
        self.assertEqual(gc_percent("GGCC"), 100)
        self.assertEqual(gc_percent("ATAT"), 0)
        self.assertEqual(gc_percent("ACGTNN"), 50)  # N ignored
        self.assertIsNone(gc_percent("NNN"))

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")

    def test_forward_orf_coordinates(self):
        orf = "ATG" + "GCT" * 98 + "TAA"      # 3 + 294 + 3 = 300 nt
        seq = "CCCC" + orf + "CCCC"
        found = find_orfs(seq, min_len=300)
        self.assertEqual(len(found), 1)
        o = found[0]
        self.assertEqual((o["strand"], o["start"], o["end"], o["length_nt"]), ("+", 5, 304, 300))

    def test_reverse_orf_coordinates(self):
        orf = "ATG" + "GCT" * 98 + "TAA"
        seq = "CCCC" + reverse_complement(orf) + "CCCC"
        found = find_orfs(seq, min_len=300)
        self.assertEqual(len(found), 1)
        o = found[0]
        self.assertEqual((o["strand"], o["start"], o["end"]), ("-", 5, 304))

    def test_short_orf_ignored(self):
        self.assertEqual(find_orfs("ATG" + "GCT" * 10 + "TAA", min_len=300), [])


if __name__ == "__main__":
    unittest.main()
