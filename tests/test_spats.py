
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
    def test_5s_cases(self):
        tgt = Target("5S", "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCTGACTCGGGCACCAAGGAC")
        tgt.index()
        self.assertEqual((0, 8, 0), tgt.find_partial("GGATGCCTTTTTTTTTTTTTTTTTTTTTTTTTTTT", 8))
        self.assertEqual((0, 8, 135), tgt.find_partial("CCAAGGACTGGAAGATCGGAAGAGCGTCGTGTAGG", 8))



class Target5STest(unittest.TestCase):
    def setUp(self):
        from spats_clean import Spats
        self.spats = Spats("test/5s/5s.fa", "test/5s")
        self.spats.setup()
    def tearDown(self):
        self.spats = None

no_error_cases = [
    [ "1101:11562:1050", "AAACGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCACCTGACCCCATGCCGAACTCAGAAGTGAAACG", 29 ],
    [ "1101:20069:1063", "TTTAGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "TCCCACCTGACCCCATGCCGAACTCAGAAGTGAAA", 27 ],
    [ "21189", "TTTGGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACCAAAAGATCGGAAGA", 123 ],
    [ "18333", "GAGTGTCCTTGGTGCCCGAGTCAGTGGTAGATCGG", "ACCACTGACTCGGGCACCAAGGACACTCAGATCGG", None ],
    [ "1101:10021:3261", "AAGCGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCTGACCCCATGCCGAACTCAGAAGTGAAACCCCG", None ],
    [ "1101:10505:2593", "TCTGGTCCTTGGTGCCCGAGTAGATCGGAAGAGAC", "ACTCGGGCACCAAGGACCAGAAGATCGGAAGAGCG", None ],
    [ "1109:9248:13419", "AGATGTCCTTGGTGCCCGAGTCAGAAGATCGGGAA", "TCTGACTCGGGCACCAAGGACATCTAGATCGGAAG", None ],
    [ "1101:10051:23846", "CTTAGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACTAAGAGATCGGAAAA", None ],
    [ "1101:13433:5831", "TTCAGTCCTTGGTGCCCGAGTCAGATAGATCGGAA", "ATCTGACTCGGGCACCAAGGACTGAAAGATCGAAA", None ],
    [ "1102:6599:2593", "AAGTGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACACTTAGATCGGAGAC", None ],
    [ "1101:12888:8140", "GGATGTCCTTGGTGCCCGAGTCAGATGCCAGATCG", "GGCATCTGACTCGGGCACCAAGGACATACAGATCG", 118 ],
    [ "1101:10652:13566", "GAATGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATG", 64 ], # tests for matching the wrong substring if your min_len
    [ "1101:13864:21135", "GGGTGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "GCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCAT", 63 ], # is too small and you don't keep looking in find_partial
    [ "1101:11920:1274", "CTTAGTCCTTGGTGCCCGAGTCAGCTTGGTGCCCG", "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCT", None ], # similarly if you don't keep checking all sites
]

spats_v102_match_cases = [

    # technically, this should be None, but is 32 to match spats v1.0.2, xref spats_config.allow_errors_in_last_four_of_R2
    [ "1101:10021:3261", "AAGCGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "CCTGACCCCATGCCGAACTCAGAAGTGAAACCCCG", 32 ],

    # technically, this should be None, but is 123 to match spats v1.0.2, xref spats_config.allow_errors_in_last_four_of_R2 (this time mimsatch is in adapter trim)
    [ "1101:10051:23846", "CTTAGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACTAAGAGATCGGAAAA", 123 ],

    # this next one is impossible to do without recreating v1.0.2 bugs
    [ "1101:10582:1913", "AAACGTCCTTGGTGCCCGAGTCAGAGATCGAAGAG", "CTGACTCGGGCACCAAGGGCGTGTATATCGGAAGA", 123 ],

    # next one may need more work, but ultimately fails b/c of a one-off toggle
    [ "1101:12888:8140", "GGATGTCCTTGGTGCCCGAGTCAGATGCCAGATCG", "GGCATCTGACTCGGGCACCAAGGACATACAGATCG", 118 ],

    # only misses one at the last 4 of adapter, but is (presumably) discarded by 1.0.2 b/c there's not enough adapter to match after trimming last 4
    [ "1101:13433:5831", "TTCAGTCCTTGGTGCCCGAGTCAGATAGATCGGAA", "ATCTGACTCGGGCACCAAGGACTGAAAGATCGAAA", None ],

    # 
    [ '1101:15138:1004', 'NTTAGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG', 'NCCGAACTCAGAAGTGAANCGCCGTAGCGCNGANG', None ],

    # one bp error in last 4, ignored by v102
    [ '2119:9713:17009', 'TTCGGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG', 'GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCTCCT', 0 ],

    # only misses in the last 4 of adapter and has sufficient length, but misses 3 bp, so is (presumably) discarded by 1.0.2 for that reason?
    [ "1102:6599:2593", "AAGTGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACACTTAGATCGGAGAC", None ],

    # v102 gets these two one off-by-one due to two bp toggles where adapter_trimmer doesn't trim the right amount
    [ '1114:24625:21410', 'AAACGTCCTTGGTGCCCGAGTCAATCGGAAGAGCA', 'ACTCGGGCACCAAGGACGCTTAGATCGGAAGAGCG', 126 ],
    [ '1109:25722:16247', 'CTCAGTCCTTGGTGCCCGAGTCAATCGGAAGAGCA', 'ACTCGGGCACCAAGGACTGAGAGATCGGCAGAGCG', 126 ],

    # adapter
    [ '1110:22635:4995', 'TTTAGTCCTTGGTGCCCGAGATCGGAAGAGCACAC', 'CGGGCACCAAGGACTAAAAGATCGGAAGAGCGTCG', 129 ],


    # one bp toggle error in len-9 adapter trim -- adapter_trimmer rejects it
    [ '1101:13433:5831', 'TTCAGTCCTTGGTGCCCGAGTCAGATAGATCGGAA', 'ATCTGACTCGGGCACCAAGGACTGAAAGATCGAAA', None ],

    # good match len-9 adapter_trim -- adapter_trimmer accepts it
    [ '1113:10835:22556', 'AAGTGTCCTTGGTGCCCGAGTCAGATAGATCGGAA', 'ATCTGACTCGGGCACCAAGGACACTTAGATCGGAA', 121 ],

    # one insertion error in R1 adapter_b causing many missed bp -- adapter_trimmer accepts it?
    [ '1103:19743:16573', 'GAGTGTCCTTGGTGCCCGAGTATCGGAAGAGCACA', 'CGGGCACCAAGGACACTCAGATCGGAAGAGCGTCG', 129 ],

    # in dev, missing from v102
    [ '1101:15979:21832', 'GAATGTCCTTGGTGCCCGAGTCAGATGCAGAACGG', 'GCATCTGACTCGGGCACCAAGGACATTCAGATCGG', None ],

    # this one is accepted by v102 b/c there's only 4 adapter, and it's chomped off by the time it gets to NOMASK
    [ '1101:10344:11542', 'TCTAGTCCTTGGTGCCCGAGTCAGATGCCTGAGAT', 'CAGGCATCTGACTCGGGCACCAAGGACTAGAAGAT', 116 ],

    # these are accepted by v102, good-match len-7 adapter
    [ '1101:11816:8298', 'CCCGGTCCTTGGTGCCCGAGTCAGATGCAGATCGG', 'GCATCTGACTCGGGCACCAAGGACCGGGAGATCGG', 119 ],
    [ '1101:11998:14960', 'AGACGTCCTTGGTGCCCGAGTCAGATGCAGATCGG', 'GCATCTGACTCGGGCACCAAGGACGTCTAGATCGG', 119 ],

    # rejected by v102, len-11 adapter with two missing bp
    [ '1101:25841:19393', 'CTTAGTCCTTGGTGCCCGAGTCAGAGACCGGAAGA', 'CTGACTCGGGCACCAAGGACTAAGGGAGCGGAAGA', None ],

    # accepted by v102, presumably due to trimming off the N of R2
    [ '1102:14595:1033', 'NTCAGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG', 'GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACNT', 0 ],

    # accepted by v102, due to trimming off bp toggle in last 4 of R2 (after adapter trimming)
    [ '2116:14830:22969', 'CTCAGTCCTTGGTGCCCGAGTCAGGATCGGAAGAG', 'TGACTCGGGCACCAAAGACTGAGAGATCGGAAGAG', 124 ],

    # i think this one caused a bug in dev...should be fixed now
    [ '1101:10344:11542', 'TCTAGTCCTTGGTGCCCGAGTCAGATGCCTGAGAT', 'CAGGCATCTGACTCGGGCACCAAGGACTAGAAGAT', 116 ],

    # accepted by v102, despite match error in last 4 of R2 (after the adapter trim)
    [ '1115:24186:4558', 'TTTGGTCCTTGGTGCCCGAGTAGATCGGAAGAGCA', 'ACTCGGGCACCAAAGACCAAAAGATCGGAAGAGCG', 126 ],

    # accepted by v102, matched at site 135, match length is only 8
    # note in reactivities.out, there are almost no matches after site 128
    # we can probably do better than this...
    [ '1109:22737:14675', 'TCCAGTCCTTGGAGATCGGAAGAGCACACGTCTGA', 'CCAAGGACTGGAAGATCGGAAGAGCGTCGTGTAGG', 135 ],

    # accepted by v102, despite insertion error on R1 adapter
    # in dev, we reject it due to many mismatches in adapter_trimmer (we're not trying to find insertions...)
    [ '1103:19743:16573', 'GAGTGTCCTTGGTGCCCGAGTATCGGAAGAGCACA', 'CGGGCACCAAGGACACTCAGATCGGAAGAGCGTCG', 129 ],

    # similar to above but a deletion error on R1 adapter
    [ '1108:11212:15952', 'TTTGGTCCTTGGTGCCCGAGTCAGAGACGGAAGAG', 'CTGACTCGGGCACCAAGGACCAAAAGATCGGAAGA', 123 ],

    # similar to above but a deletion error on R2 adapter
    [ '1107:6287:7763', 'TCTGGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA', 'CTGACTCGGGCACCAAGGACCAGAAGATCGAAGAG', 123 ],


    # below is a spot-check of cases that dev (with v102 compat options) gets but 5s does not
    [ '*1102:16621:23746', 'TCCAGTCCTTGGTGCCCGAGTCAGGATCGGAAGAG', 'TGACTCGGGCACCAAGGCCTGGAAGAACGGAAGAA', None ],
    [ '*1102:27276:10366', 'TTTAGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA', 'CTGACTCGGGCACCAAGGACTAAAAGATTGGAAAA', None ],
    [ '*1105:12183:24798', 'AGATGTCCTTGGTGCCCGAGTCAGATCGGAAGAGC', 'GACTCGGGCACCAAGGACATCTAGATCGGAAAACC', None ],
    [ '*1105:28564:18308', 'CTTAGTCCTTGGTGCCCGAGTAGATCGGAAGAGCA', 'ACTCGGGCACCAAGGACTAAGGGATAGGAAGAGCG', None ],
    [ '*1109:14013:17918', 'GGGTGTCCTTGGTGCCCGAGTCAGATGCCTAGATC', 'AGGCATCTGACTCGGGCACCAAGGCCACCCAGATC', None ],
    [ '*1116:21540:21211', 'AGATGTCCTTGGTGCCCGAGTCAGATGCCTGAGAT', 'CAGGCATCTGACTCGGGCACCAAGCACATCTAGAT', None ],
    [ '*2104:18173:12075', 'AAGCGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG', 'ACTGCCAGGCATCTGACTCGGGCACCAAGGGCGCT', None ],
    [ '*2113:28288:17700', 'TCCGGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG', 'AACTGCCAGGCATCTGACTCGGGCACCAACGACCG', None ],
    [ '*2118:26021:18628', 'TTTAGTCCTTGGTGCCCGAGTCAGATGCCTGGCAA', 'TGCCAGGCATCTGACTCGGGCACCAAGGCCTAAAA', None ],
    [ '*1114:21343:8367', 'TTCAGTCCTTGGTGCCCGAGTGATCGGAAGAGCAC', 'CTCGGGCACCAAGGACTGAAAGCTCGGAAGAGCGA', None ],    

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
        self.assertEqual(case[3], pair.site, msg = case[0])
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

    def test_find_partial_weird_case(self):
        pair = Pair()
        pair.set_from_data("x", 'CTCAGTCCTTGGTGCCCGAGTCAGGATCGGAAGAG', 'TGACTCGGGCACCAAAGACTGAGAGATCGGAAGAG')
        self.spats.process_pair(pair)
        print "{} / {}".format(pair.site, pair.failure)
