
import unittest


from spats_clean import longest_match, reverse_complement, string_match_errors
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
    [ "1101:10652:13566", "GAATGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATG", 'RRRY', 64 ], # tests for matching the wrong substring if your min_len
    [ "1101:13864:21135", "GGGTGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "GCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCAT", 'RRRY', 63 ], # is too small and you don't keep looking in find_partial
    [ "1101:11920:1274", "CTTAGTCCTTGGTGCCCGAGTCAGCTTGGTGCCCG", "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCT", 'YYYR', None ], # similarly if you don't keep checking all sites
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
        self.assertEqual(case[3], pair.mask.chars, msg = case[0])
        self.assertEqual(case[4], pair.site, msg = case[0])
        return pair
        
    def test_case(self):
        spats_config.debug = True
        for case in no_error_cases + spats_v102_match_cases:
            if case[0].startswith("*"):
                self.run_case(case, show_diagram = True)
        spats_config.debug = False

    def test_pairs(self):
        for case in no_error_cases:
            self.run_case(case)
        print "Ran {} pair->site cases.".format(len(no_error_cases))

    def test_check_all_sites(self):
        pair = Pair()
        pair.set_from_data("x", "CTTAGTCCTTGGTGCCCGAGTCAGCTTGGTGCCCG", "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCT")
        self.spats.process_pair(pair)
        self.assertEqual(143, pair.r1.right)
        self.assertEqual(6, len(pair.r1.match_errors))
