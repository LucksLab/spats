
import unittest

from spats_common import longest_match, reverse_complement, Target



class TestLongestMatch(unittest.TestCase):
    def test_longest_match(self):
        self.assertEqual(longest_match("abc", (0, 1), "frabda", (2, 1)), (0, 1))
        self.assertEqual(longest_match("abc", (0, 1), "fradda", (2, 1)), (0, 0))
        self.assertEqual(longest_match("abc", (0, 1), "frabca", (2, 1)), (0, 2))
        self.assertEqual(longest_match("xabc", (1, 1), "frabca", (2, 1)), (0, 2))
        self.assertEqual(longest_match("xabc", (1, 1), "fxabca", (2, 1)), (1, 2))


class SRPTargetTest(unittest.TestCase):
    def setUp(self):
        self.target = Target("SRP",
                             "ATCGGGGGCTCTGTTG" + 
                             "GTTCTCCCGCAACGCT" + 
                             "ACTCTGTTTACCAGGT" + 
                             "CAGGTCCGGAAGGAAG" + 
                             "CAGCCAAGGCAGATGA" + 
                             "CGCGTGTGCCGGGATG" + 
                             "TAGCTGGCAGGGCCCC" + 
                             "CACCCGTCCTTGGTGC" + 
                             "CCGAGTCAG")
        self.target.index()
    def tearDown(self):
        self.target = None


class TestTarget(SRPTargetTest):
    def test_exact(self):
        self.assertEqual(self.target.find_exact('ATCGGGGGCT'), 0)
        self.assertEqual(self.target.find_exact('TCTGTTGGTTCTC'), 9)
        self.assertEqual(self.target.find_exact('CCCCCCCCCCC'), None)
        self.assertEqual(self.target.find_exact('CCCCCCCCCCC'), None)
        self.assertEqual(self.target.find_exact('TCTGTTGGTTCCC'), None)
    def test_partial(self):
        self.assertEqual(self.target.find_partial('ATCGGGGGCTCTGTTGGTT'), (0, 19, 0))
        self.assertEqual(self.target.find_partial('ATGGGGGGCTCTGTTGGTT', 12), (3, 16, 3))
        self.assertEqual(self.target.find_partial('CCCCC' + 'CAGCCAAGGCAGATGA' + 'GGGGG', 12), (5, 16, 64))
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

class TestMisc(SRPTargetTest):
    def test_id_to_site(self):
        from spats_common import id_to_site
        bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/t7/"
        id_to_site(bp + "NOMASK_1.fq", bp + "RRRY.sam", bp + "YYYR.sam", 143)

    def test_make_subset(self):
        from spats_common import make_subset
        bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
        if False:
            make_subset(bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                        bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
                        bp + "t7/5s_dev_diff_ids.out",
                        bp + "5s_errors")
        else:
            make_subset(bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                        bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
                        bp + "t7/5s_dev_diff_sample.ids",
                        bp + "5s_sample")

    def test_misc(self):
        from spats_common import spats
        if False:
            bp = "/Users/jbrink/mos/tasks/1RwIBa/refactor/spats/test/Read_Mapping/"
            spats(bp + "SRP_All_Stops.fa", bp + "SRP_All_Stops_R1.fq", bp + "SRP_All_Stops_R2.fq", bp + "t5")
        elif False:
            bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
            spats(bp + "5s/5S.fa",
                  bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                  bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
                  bp + "t7")
        elif True:
            bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
            spats(bp + "5s/5S.fa",
                  bp + "5s_sample/filtered_R1.fq",
                  bp + "5s_sample/filtered_R2.fq",
                  bp + "t9")
        else:
            bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/5sq_dev/"
            spats(bp + "5S.fa",
                  bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq", 
                  bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq", 
                  bp + "t2")
