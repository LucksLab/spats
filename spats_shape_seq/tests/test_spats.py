
import unittest


from spats_shape_seq.util import reverse_complement, string_find_errors, string_match_errors
from spats_shape_seq.mask import longest_match

class TestUtils(unittest.TestCase):
    def test_longest_match(self):
        self.assertEqual((0, 1), longest_match("ATC", (0, 1), "TTATGA", (2, 1)))
        self.assertEqual((0, 0), longest_match("ATC", (0, 1), "TTAGGA", (2, 1)))
        self.assertEqual((0, 2), longest_match("ATC", (0, 1), "TTATCA", (2, 1)))
        self.assertEqual((0, 2), longest_match("GATC", (1, 1), "TTATCA", (2, 1)))
        self.assertEqual((1, 2), longest_match("GATC", (1, 1), "TGATCA", (2, 1)))
    def test_string_match(self):
        self.assertEqual([], string_match_errors("GATC", "GATC"))
        self.assertEqual([2], string_match_errors("GATC", "GACC"))
        self.assertEqual([0], string_match_errors("GATC", "AATC"))
        self.assertEqual([0, 3], string_match_errors("GATC", "CATG"))
        self.assertEqual(range(4), string_match_errors("GATC", "CTAG"))
    def test_reverse_complement(self):
        self.assertEqual("GATC", reverse_complement("GATC"))
        self.assertEqual("CGTCCAA", reverse_complement("TTGGACG"))
        self.assertEqual("CAACAGAGCCCCCGAT", reverse_complement("ATCGGGGGCTCTGTTG"))
        self.assertEqual("GATNC", reverse_complement("GNATC"))
    def test_string_find(self):
        m = string_find_errors("TTTT", "TTATAGGCGATGGAGTTCGCCATAAACGCTGCTTAGCTAATGACTCCTACCAGTATCACTACTGGTAGGAGTCTATTTTTTTAGGAGGAAGGATCTATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTC", 0, 1)
        self.assertTrue(len(m) == 1  and  len(set(m) & set([75, 76, 77, 78, 115])) == 1)
        m = string_find_errors("TTTT", "TTATAGGCGATGGAGTTCGCCATAAACGCTGCTTAGCTAATGACTCCTACCAGTATCACTACTGGTAGGAGTCTATTTTTTTAGGAGGAAGGATCTATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTC", 0, 4)
        self.assertTrue(len(m) == 4  and  len(set(m) & set([75, 76, 77, 78, 115])) == 4)
        m = string_find_errors("TTTT", "TTATAGGCGATGGAGTTCGCCATAAACGCTGCTTAGCTAATGACTCCTACCAGTATCACTACTGGTAGGAGTCTATTTTTTTAGGAGGAAGGATCTATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTC", 0, 6)
        self.assertTrue(len(m) == 5  and  len(set(m) & set([75, 76, 77, 78, 115])) == 5)

        m = string_find_errors("ACAT", "ATCGGGGGCTCTGTTGGTTCCCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCCGGGTCGGCATGGCATCTCCACCTCCTCGCGGT", 1, 1)
        self.assertTrue(len(m) == 1  and  len(set(m) & set([74, 123, 128, 136])) == 1)
        m = string_find_errors("ACAT", "ATCGGGGGCTCTGTTGGTTCCCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCCGGGTCGGCATGGCATCTCCACCTCCTCGCGGT", 1, 3)
        self.assertTrue(len(m) == 3  and  len(set(m) & set([74, 123, 128, 136])) == 3)
        m = string_find_errors("ACAT", "ATCGGGGGCTCTGTTGGTTCCCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCCGGGTCGGCATGGCATCTCCACCTCCTCGCGGT", 1, 4)
        self.assertTrue(len(m) == 4  and  len(set(m) & set([74, 123, 128, 136])) == 4)

        m = string_find_errors("GCAT", "AAAACCCCGGGGTTTTATATACGTCAGCCC", 2, 1)
        self.assertTrue(len(m) == 1  and  len(set(m) & set([9, 10, 11, 14, 16, 20, 23, 26])) == 1)
        m = string_find_errors("GCAT", "AAAACCCCGGGGTTTTATATACGTCAGCCC", 2, 4)
        self.assertTrue(len(m) == 4  and  len(set(m) & set([9, 10, 11, 14, 16, 20, 23, 26])) == 4)
        m = string_find_errors("GCAT", "AAAACCCCGGGGTTTTATATACGTCAGCCC", 2, 8)
        self.assertTrue(len(m) == 8  and  len(set(m) & set([9, 10, 11, 14, 16, 20, 23, 26])) == 8)
