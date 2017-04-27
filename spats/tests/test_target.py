
import unittest

from spats.util import reverse_complement
from spats.target import Target

TARGET_SRP = "ATCGGGGGCTCTGTTGGTTCTCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCCGTCCTTGGTGCCCGAGTCAG"
TARGET_5S = open("test/5s/5s.fa", 'rb').read().split('\n')[1]

class TestTargetIndexing(unittest.TestCase):
    def test_index_srp(self):
        target = Target("SRP", TARGET_SRP)
        target.index()
        self.assertEqual(8, target.longest_self_match())
    def test_index_5S(self):
        target = Target("5S", TARGET_5S)
        target.index()
        self.assertEqual(10, target.longest_self_match())
    def test_5s_cases(self):
        tgt = Target("5S", TARGET_5S)
        tgt.index()
        self.assertEqual((0, 8, 0), tgt.find_partial("GGATGCCTTTTTTTTTTTTTTTTTTTTTTTTTTTT"))
        self.assertEqual((0, 8, 135), tgt.find_partial("CCAAGGACTGGAAGATCGGAAGAGCGTCGTGTAGG"))
        self.assertEqual((11, 20, 123), tgt.find_partial("CGGGCACCAAGCTGACTCGGGCACCAAGGAC"))


class SRPTargetTest(unittest.TestCase):
    def setUp(self):
        self.target = Target("SRP", TARGET_SRP)
        self.target.index()
    def tearDown(self):
        self.target = None

class TestTarget(SRPTargetTest):
    def test_exact(self):
        self.assertEqual(self.target.find_exact('ATCGGGGGCT'), 0)
        self.assertEqual(self.target.find_exact('TCTGTTGGTTCTC'), 9)
        self.assertEqual(self.target.find_exact('CCCCCCCCCCC'), -1)
        self.assertEqual(self.target.find_exact('CCCCCCCCCCC'), -1)
        self.assertEqual(self.target.find_exact('TCTGTTGGTTCCC'), -1)
    def test_partial(self):
        self.assertEqual(self.target.find_partial('ATCGGGGGCTCTGTTGGTT'), (0, 19, 0))
        old_min = self.target.minimum_match_length
        self.target.minimum_match_length = 12
        self.assertEqual(self.target.find_partial('ATGGGGGGCTCTGTTGGTT'), (3, 16, 3))
        self.assertEqual(self.target.find_partial('CCCCC' + 'CAGCCAAGGCAGATGA' + 'GGGGG'), (5, 16, 64))
        self.target.minimum_match_length = old_min
    def test_SRPs(self):
        self.assertEqual(self.target.find_partial(reverse_complement("GGGCCTGACTCGGGCACCAAGGACGGGTGGGGGCC")), (0, 31, 106)) #R1 0_0
        self.assertEqual(self.target.find_partial(reverse_complement("CCCGCTGACTCGGGCACCAAGGACGGGTGGGGGCC")), (0, 31, 106)) #R1 0_1
        self.assertEqual(self.target.find_partial(reverse_complement("CCCGCTGACTCGGGCACCAAGGACGGGTGGGGGCC")), (0, 31, 106)) #R1 106
        self.assertEqual(self.target.find_partial(reverse_complement("GGGCCTGACTCGGGCACCAAGGACGGGTGGGGGCA")), (1, 30, 107)) #R1 107
        self.assertEqual(self.target.find_partial(reverse_complement("GGGCCTGACTCGGGCACCAAGGACAGATCGGAAGA")), (11, 20, 117)) #R1 117
        self.assertEqual(self.target.find_partial("ATCGGGGGCTCTGTTGGTTCTCCCGCAACGCTACT"), (0, 35, 0)) #R2 0_0
        self.assertEqual(self.target.find_partial("ATCGGGGGCTCTGTTGGTTCTCCCGCAACGCTACT"), (0, 35, 0)) #R2 0_1
        self.assertEqual(self.target.find_partial("GCAGGGCCCCCACCCGTCCTTGGTGCCCGAGTCAG"), (0, 35, 102)) #R2 102
        self.assertEqual(self.target.find_partial("CAGGGCCCCCACCCGTCCTTGGTGCCCGAGTCAGG"), (0, 34, 103)) #R2 103
        self.assertEqual(self.target.find_partial("GGCCCCCACCCGTCCTTGGTGCCCGAGTCAGGCCC"), (0, 31, 106)) #R2 106
        self.assertEqual(self.target.find_partial("GCCCCCACCCGTCCTTGGTGCCCGAGTCAGGCCCA"), (0, 30, 107)) #R2 107
        self.assertEqual(self.target.find_partial("GTCCTTGGTGCCCGAGTCAGGCCCAGATCGGAAGA"), (0, 20, 117)) #R2 117
