
import unittest


from spats_clean import longest_match, hamming_distance, reverse_complement
class TestUtils(unittest.TestCase):
    def test_longest_match(self):
        self.assertEqual(longest_match("ATC", (0, 1), "TTATGA", (2, 1)), (0, 1))
        self.assertEqual(longest_match("ATC", (0, 1), "TTAGGA", (2, 1)), (0, 0))
        self.assertEqual(longest_match("ATC", (0, 1), "TTATCA", (2, 1)), (0, 2))
        self.assertEqual(longest_match("GATC", (1, 1), "TTATCA", (2, 1)), (0, 2))
        self.assertEqual(longest_match("GATC", (1, 1), "TGATCA", (2, 1)), (1, 2))
    def test_hamming_distance(self):
        self.assertEqual(hamming_distance("GATC", "GATC"), 0)
        self.assertEqual(hamming_distance("GATC", "GACC"), 1)
        self.assertEqual(hamming_distance("GATC", "AATC"), 1)
        self.assertEqual(hamming_distance("GATC", "CATG"), 2)
        self.assertEqual(hamming_distance("GATC", "CTAG"), 4)
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("GATC"), "GATC")
        self.assertEqual(reverse_complement("TTGGACG"), "CGTCCAA")
        self.assertEqual(reverse_complement("ATCGGGGGCTCTGTTG"), "CAACAGAGCCCCCGAT")


from spats_clean import Target
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


from spats_clean import reverse_complement
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


class Target5STest(unittest.TestCase):
    def setUp(self):
        from spats_clean import Spats
        self.spats = Spats("test/5s/5s.fa", "test/5s")
        self.spats.setup()
    def tearDown(self):
        self.spats = None

no_error_cases = [
    [ "1101:11562:1050", "AAACGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCACCTGACCCCATGCCGAACTCAGAAGTGAAACG", 'RRRY', 29 ],
    [ "1101:20069:1063", "TTTAGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "TCCCACCTGACCCCATGCCGAACTCAGAAGTGAAA", 'YYYR', 27 ],
    [ "21189", "TTTGGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACCAAAAGATCGGAAGA", 'YYYR', 123 ],
    [ "18333", "GAGTGTCCTTGGTGCCCGAGTCAGTGGTAGATCGG", "ACCACTGACTCGGGCACCAAGGACACTCAGATCGG", 'RRRY', None ],
    [ "1101:10021:3261", "AAGCGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCTGACCCCATGCCGAACTCAGAAGTGAAACCCCG", 'RRRY', None ],
    [ "1101:10505:2593", "TCTGGTCCTTGGTGCCCGAGTAGATCGGAAGAGAC", "ACTCGGGCACCAAGGACCAGAAGATCGGAAGAGCG", 'YYYR', None ],
    [ "1109:9248:13419", "AGATGTCCTTGGTGCCCGAGTCAGAAGATCGGGAA", "TCTGACTCGGGCACCAAGGACATCTAGATCGGAAG", 'RRRY', None ],
    [ "1101:10051:23846", "CTTAGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACTAAGAGATCGGAAAA", 'YYYR', None ],
    [ "1101:13433:5831", "TTCAGTCCTTGGTGCCCGAGTCAGATAGATCGGAA", "ATCTGACTCGGGCACCAAGGACTGAAAGATCGAAA", 'YYYR', None ],
    [ "1102:6599:2593", "AAGTGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACACTTAGATCGGAGAC", 'RRRY', None ],
    [ "1101:12888:8140", "GGATGTCCTTGGTGCCCGAGTCAGATGCCAGATCG", "GGCATCTGACTCGGGCACCAAGGACATACAGATCG", 'RRRY', 118 ],
    # next 2 are good tests for matching the wrong substring if your min_len is too small and you don't keep looking in find_partial
    [ "1101:10652:13566", "GAATGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATG", 'RRRY', 64 ],
    [ "1101:13864:21135", "GGGTGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "GCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCAT", 'RRRY', 63 ],
]

spats_v102_match_cases = [

    # technically, this should be None, but is 32 to match spats v1.0.2, xref spats_config.allow_errors_in_last_four_of_R2
    [ "1101:10021:3261", "AAGCGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCTGACCCCATGCCGAACTCAGAAGTGAAACCCCG", 'RRRY', 32 ],

    # technically, this should be None, but is 123 to match spats v1.0.2, xref spats_config.allow_errors_in_last_four_of_R2 (this time mimsatch is in adapter trim)
    [ "1101:10051:23846", "CTTAGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACTAAGAGATCGGAAAA", 'YYYR', 123 ],

    # only misses one at the last 4 of adapter, but is (presumably) discarded by 1.0.2 b/c there's not enough adapter to match after trimming last 4
    [ "1101:13433:5831", "TTCAGTCCTTGGTGCCCGAGTCAGATAGATCGGAA", "ATCTGACTCGGGCACCAAGGACTGAAAGATCGAAA", 'YYYR', None ],

    # only misses in the last 4 of adapter and has sufficient length, but misses 3 bp, so is (presumably) discarded by 1.0.2 for that reason?
    [ "1102:6599:2593", "AAGTGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACACTTAGATCGGAGAC", 'RRRY', None ],

    # this next one is impossible to do without recreating v1.0.2 bugs
    [ "1101:10582:1913", "AAACGTCCTTGGTGCCCGAGTCAGAGATCGAAGAG", "CTGACTCGGGCACCAAGGGCGTGTATATCGGAAGA", 'RRRY', 123 ],

    # next one may need more work, but ultimately fails b/c of a one-off toggle
    [ "1101:12888:8140", "GGATGTCCTTGGTGCCCGAGTCAGATGCCAGATCG", "GGCATCTGACTCGGGCACCAAGGACATACAGATCG", 'RRRY', 118 ],
]

# cases that require the ability to handle one or more nt's toggled/mistranscribed
toggle_cases = [
    [ "1101:10051:23846", "CTTAGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACTAAGAGATCGGAAAA", 'YYYR', 123 ],  # requires adapter-trimming to be able to handle toggle
#    [ "1101:10337:2165", "AAGTGTCCTTGGTGCCCGAGTCAGATGCCTGGCCG", "GCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATG", 'RRRY', ??? ],
]    

diagram_cases = [
    no_error_cases[0],
    [ "1101:11562:1050 mask tweaked", "CAACGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCACCTGACCCCATGCCGAACTCAGAAGTGAAACG", ],
    [ "1101:11562:1050 R2 tweaked", "AAACGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCATCTTACCCTTTTCCGTACTCTTAAGTGTAATG" ],
]

from spats_clean import Pair, spats_config
from diagram import diagram

class TestPairs(Target5STest):

    def pair_for_case(self, case):
        pair = Pair()
        pair.set_from_data(case[0], case[1], case[2])
        return pair

    def run_case(self, case, show_diagram = False):
        pair = self.pair_for_case(case)
        self.spats.process_pair(pair)
        if show_diagram:
            print diagram(self.spats._target, pair)
        self.assertEqual(pair.mask.chars, case[3], msg = case[0])
        self.assertEqual(pair.site, case[4], msg = case[0])
        return pair
        
    def test_case(self):
        spats_config.debug = True
        for case in no_error_cases + spats_v102_match_cases:
            if case[0].startswith("*"):
                self.run_case(case, show_diagram = True)
        
    def test_adapter_trim(self):
        self.run_case(pair_cases[3])
        #print diagram(self.spats._target, pair)

    def test_pairs(self):
        for case in no_error_cases:
            self.run_case(case)
        print "Ran {} pair->site cases.".format(len(no_error_cases))

    def test_diagram(self):
        pair = Pair()
        for case in diagram_cases:
            pair.set_from_data(case[0], case[1], case[2])
            self.spats.process_pair(pair)
            print diagram(self.spats._target, pair)
            print "\n\n"


#TODO: DELME
# just keeping for some usage examples of file-grepping of prev.gen. tools
class TestMisc(unittest.TestCase): # inherit unittest.TestCase to use
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
        if False:
            bp = "/Users/jbrink/mos/tasks/1RwIBa/refactor/spats/test/Read_Mapping/"
            self.spats(bp + "SRP_All_Stops.fa", bp + "SRP_All_Stops_R1.fq", bp + "SRP_All_Stops_R2.fq", bp + "t5")
        elif True:
            bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
            self.spats(bp + "5s/5S.fa",
                       bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq",
                       bp + "5s/data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq",
                       bp + "t11")
        elif False:
            bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/"
            self.spats(bp + "5s/5S.fa",
                       bp + "5s_sample/filtered_R1.fq",
                       bp + "5s_sample/filtered_R2.fq",
                       bp + "t11")
        else:
            bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/5sq_dev/"
            self.spats(bp + "5S.fa",
                       bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq", 
                       bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq", 
                       bp + "t4")

    def spats(self, target, r1, r2, out, show_sites = True):
        from spats_clean import Spats, spats_config
        s = Spats(target, out)
        s.setup()
        if show_sites:
            spats_config.show_id_to_site = True
        s.process_pair_data(r1, r2)
        if not show_sites:
            s.compute_profiles()
            s.write_reactivities()
              
    def test_refactor(self):
        from spats_clean import Spats
        bp = "/Users/jbrink/mos/tasks/1RwIBa/tmp/5sq_dev/"
        out = bp + "t3/"
        s = Spats(bp + "5S.fa", out)
        s.setup()
        s.process_pair_data(bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R1_001.fastq", 
                            bp + "data/17571-AD1AW-KEW11-5S-2p1-18x-23FEB15-GGCTAC_S10_L001_R2_001.fastq")
        s.compute_profiles()
        s.write_reactivities()
        import subprocess
        subprocess.check_call(["diff", bp + "t2/rx.out", out + "/rx.out"])
        print "Diff OK"
