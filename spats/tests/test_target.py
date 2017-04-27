
import unittest

from spats.parse import fasta_parse
from spats.target import Targets
from spats.util import reverse_complement

TARGET_SRP = "ATCGGGGGCTCTGTTGGTTCTCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCCGTCCTTGGTGCCCGAGTCAG"
TARGET_5S = open("test/5s/5s.fa", 'rb').read().split('\n')[1]

class TestMultiple(unittest.TestCase):

    def test_multiple(self):
        target = Targets()
        target.addTarget("SRP", TARGET_SRP)
        target.addTarget("5S", TARGET_5S)
        target.index()
        self.assertEqual(10, target.longest_self_match())
        tgt, s, l, i = target.find_partial("GGATGCCTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
        self.assertEqual(("5S", 0, 8, 0), (tgt.name, s, l, i))
        tgt, s, l, i = target.find_partial('ATGGGGGGCTCTGTTGGTT')
        self.assertEqual(("SRP", 3, 16, 3), (tgt.name, s, l, i))

    def test_multiple_parse(self):
        target = Targets()
        for name, seq in fasta_parse("test/panel_RNAs/panel_RNAs_complete.fa"):
            target.addTarget(name, seq)
        self.assertEqual(14, len(target.targets))
        target.index()
        self.assertEqual(169, target.longest_self_match())


class TestTargetIndexing(unittest.TestCase):
    def test_index_srp(self):
        target = Targets()
        target.addTarget("SRP", TARGET_SRP)
        target.index()
        self.assertEqual(8, target.longest_self_match())
    def test_index_5S(self):
        target = Targets()
        target.addTarget("5S", TARGET_5S)
        target.index()
        self.assertEqual(10, target.longest_self_match())
    def test_5s_cases(self):
        tgt = Targets()
        tgt.addTarget("5S", TARGET_5S)
        tgt.index()
        self.assertEqual((0, 8, 0), tgt.find_partial("GGATGCCTTTTTTTTTTTTTTTTTTTTTTTTTTTT")[1:])
        self.assertEqual((0, 8, 135), tgt.find_partial("CCAAGGACTGGAAGATCGGAAGAGCGTCGTGTAGG")[1:])
        self.assertEqual((11, 20, 123), tgt.find_partial("CGGGCACCAAGCTGACTCGGGCACCAAGGAC")[1:])

class SRPTargetTest(unittest.TestCase):
    def setUp(self):
        self.target = Targets()
        self.target.addTarget("SRP", TARGET_SRP)
        self.target.index()
    def tearDown(self):
        self.target = None

class TestTarget(SRPTargetTest):
    def test_exact(self):
        self.assertEqual(self.target.find_exact('ATCGGGGGCT')[1], 0)
        self.assertEqual(self.target.find_exact('TCTGTTGGTTCTC')[1], 9)
        self.assertEqual(self.target.find_exact('CCCCCCCCCCC')[0], None)
        self.assertEqual(self.target.find_exact('CCCCCCCCCCC')[0], None)
        self.assertEqual(self.target.find_exact('TCTGTTGGTTCCC')[0], None)
    def test_partial(self):
        self.assertEqual(self.target.find_partial('ATCGGGGGCTCTGTTGGTT')[1:], (0, 19, 0))
        old_min = self.target.minimum_match_length
        self.target.minimum_match_length = 12
        self.assertEqual(self.target.find_partial('ATGGGGGGCTCTGTTGGTT')[1:], (3, 16, 3))
        self.assertEqual(self.target.find_partial('CCCCC' + 'CAGCCAAGGCAGATGA' + 'GGGGG')[1:], (5, 16, 64))
        self.target.minimum_match_length = old_min
    def test_SRPs(self):
        self.assertEqual(self.target.find_partial(reverse_complement("GGGCCTGACTCGGGCACCAAGGACGGGTGGGGGCC"))[1:], (0, 31, 106)) #R1 0_0
        self.assertEqual(self.target.find_partial(reverse_complement("CCCGCTGACTCGGGCACCAAGGACGGGTGGGGGCC"))[1:], (0, 31, 106)) #R1 0_1
        self.assertEqual(self.target.find_partial(reverse_complement("CCCGCTGACTCGGGCACCAAGGACGGGTGGGGGCC"))[1:], (0, 31, 106)) #R1 106
        self.assertEqual(self.target.find_partial(reverse_complement("GGGCCTGACTCGGGCACCAAGGACGGGTGGGGGCA"))[1:], (1, 30, 107)) #R1 107
        self.assertEqual(self.target.find_partial(reverse_complement("GGGCCTGACTCGGGCACCAAGGACAGATCGGAAGA"))[1:], (11, 20, 117)) #R1 117
        self.assertEqual(self.target.find_partial("ATCGGGGGCTCTGTTGGTTCTCCCGCAACGCTACT")[1:], (0, 35, 0)) #R2 0_0
        self.assertEqual(self.target.find_partial("ATCGGGGGCTCTGTTGGTTCTCCCGCAACGCTACT")[1:], (0, 35, 0)) #R2 0_1
        self.assertEqual(self.target.find_partial("GCAGGGCCCCCACCCGTCCTTGGTGCCCGAGTCAG")[1:], (0, 35, 102)) #R2 102
        self.assertEqual(self.target.find_partial("CAGGGCCCCCACCCGTCCTTGGTGCCCGAGTCAGG")[1:], (0, 34, 103)) #R2 103
        self.assertEqual(self.target.find_partial("GGCCCCCACCCGTCCTTGGTGCCCGAGTCAGGCCC")[1:], (0, 31, 106)) #R2 106
        self.assertEqual(self.target.find_partial("GCCCCCACCCGTCCTTGGTGCCCGAGTCAGGCCCA")[1:], (0, 30, 107)) #R2 107
        self.assertEqual(self.target.find_partial("GTCCTTGGTGCCCGAGTCAGGCCCAGATCGGAAGA")[1:], (0, 20, 117)) #R2 117
