
import unittest

from spats_shape_seq.util import reverse_complement, string_find_errors, string_match_errors, align_strings, char_sim
from spats_shape_seq.mask import longest_match, base_similarity_ind


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

    def test_align_strings(self):
        ## Note more testing of this is done in the json test suite.
        ## Values used by spats:
        MV = 3
        MMC = 2
        GOC = 5
        GEC = 1
        simfn = lambda a,b: base_similarity_ind(a, b, MV, MMC, .5*MV)

        R = "GGMCSCGATGCCGNACGATKTAAGTCCGAGCATCAACTATGCCCTACCTGCTTCGRCCGATAAAGCTTTCAAWAGACGAYAAT"
        T = "GGACCCGATGCCGGACGAAAGTCCGCGCATCAACTATGCCTCTACCTGCTTCGGCCGATAAAGCCGACGATAATACTCCCAAAGCCC"
        a = align_strings(R, T, simfn, GOC, GEC)
        self.assertEqual(a.score, 179.5)
        self.assertEqual(a.target_match_start, 0)
        self.assertEqual(a.target_match_end, 73)
        self.assertEqual(a.src_match_start, 0)
        self.assertEqual(a.src_match_end, 82)
        self.assertEqual(a.mismatched, [25])
        self.assertEqual(a._indels_as_dict(), { 18: { "insert_type": True, "seq": "TKT" },
                                                40: { "insert_type": False, "seq": "T" },
                                                64: { "insert_type": True, "seq": "TTT" },
                                                65: { "insert_type": True, "seq": "AAWA" } })

        ## The following tests exercise the 'penalize_ends' part of align_strings()

        ## First we want to try both sides of:
        ##     2*GOC + GEC*(len(prefix1)+len(prefix2)) < MMC * min(len(prefix1), len(prefix2))
        ##     -->  MMC = 5 and 4
        R = "CCCCCAAAAAAAAAAAAAAAAAAAAAACCCC"
        T =  "TTTTAAAAAAAAAAAAAAAAAAAAAATTTTT"

        MMC = 5    # all indels
        simfn = lambda a,b: char_sim(a, b, MV, MMC)
        a = align_strings(R, T, simfn, GOC, GEC)
        self.assertEqual(a.score, 32.0)
        self.assertEqual(a.target_match_start, 0)
        self.assertEqual(a.target_match_end, 29)
        self.assertEqual(a.src_match_start, 1)
        self.assertEqual(a.src_match_end, 30)
        self.assertEqual(len(a.mismatched), 0)
        self.assertEqual(a._indels_as_dict(), { 0: { 'insert_type': True, 'seq': "CCCCC" },
                                                3: { 'insert_type': False, 'seq': "TTTT" },
                                               26: { 'insert_type': True, 'seq': "CCCC" },
                                               30: { 'insert_type': False, 'seq': "TTTTT" } })

        MMC = 4    # all mismatches
        simfn = lambda a,b: char_sim(a, b, MV, MMC)
        a = align_strings(R, T, simfn, GOC, GEC)
        self.assertEqual(a.score, 34.0)
        self.assertEqual(a.target_match_start, 0)
        self.assertEqual(a.target_match_end, 29)
        self.assertEqual(a.src_match_start, 1)
        self.assertEqual(a.src_match_end, 30)
        self.assertEqual(len(a.indels), 0)
        self.assertTrue(set(a.mismatched), set([0, 1, 2, 3, 26, 27, 28, 29]))

        R = "TTTCCCCCAAAAGACGATAAT"
        T = "CCCCCGACGATAATACTCCCAAAGCCCACCCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
        MMC = 2
        simfn = lambda a,b: char_sim(a, b, MV, MMC)
        a = align_strings(R, T, simfn, GOC, GEC)
        self.assertEqual(a.score, 34.0)
        self.assertEqual(a.target_match_start, 0)
        self.assertEqual(a.target_match_end, 13)
        self.assertEqual(a.src_match_start, 3)
        self.assertEqual(a.src_match_end, 20)
        self.assertEqual(a._indels_as_dict(), { 5: { 'insert_type': True, 'seq': "AAAA" } })
        self.assertEqual(len(a.mismatched), 0)
