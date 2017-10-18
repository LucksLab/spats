import unittest

from spats_shape_seq import Spats
from spats_shape_seq.pair import Pair

# [ id, r1, r2, end, site, muts ]
cases = [

    # indexing convention, xref https://trello.com/c/2qIGo9ZR/201-stop-map-mutation-indexing-convention
    [ 'x1', 'TTCAGTCCTTGGTGCCCGAGTCAGCCAGCTACATCCCGGCACACGCGTCATCTGCCTTGGCTGCTTCCTTCCGGA', 'AGGTCAGATCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCTGACTCGGGCACCAA', 103, 44, [52] ],

    # R1/R2 disagree on overlap; similar to x3, may need to factor in quality at some point
    [ 'x2', 'AGGCGTCCTTGGTGCCCGAGTCAGCCTTGGCTGCTTCCTTCCGGACCTGACCTGGTAAACAGAGTAGCGTTGCGG', 'ATCGGGGGCTCTGTTGGTTCCCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGTCT', None, None, [], ],

    # note: on this one, the quality score on R1 vs. R2 indicates that mut 73 really is a misread on R2, so this case may change in the future (to be 88, 0, [88])...
    # for now it should be tossed b/c R1 and R2 disagree where they overlap
    # xref https://trello.com/c/35mBHvPA/197-stop-map-r1-r2-disagree-case
    [ 'x3', 'GAATGTCCTTGGTGCCCGAGTCAGGACACGCGTCATCTGCCTTGGCTGCTTCCTTCCGGACCTGACCTGGTAAAC', 'ATCGGGGGCTCTGTTGGTTCCCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGTCA', None, None, [] ],

    [ 'x4', 'GAGCGTCCTTGGTGCCCGAGTCAGATGCCGACCCGGGTGGGGGCCCTGCCAGCTACATCCCGGCACACGCGTCAT', 'TAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCCG', 127, 43, [44] ],

    [ 'x5', 'GAATGTCCTTGGTGCCCGAGTCAGTCCTTGGTGCCCGAGTCAGTCCTTGGTTCCCGAGTCACTCCTTTGTTCCCC', 'AGGACTGACTCGGGCACCAAGGACTTTCTCGTTCACCTATTTCTTTCTCTTCCCCCTTTTTCTTTCTCTTTCTCC', None, None, [] ],

    # this one is slightly ambiguous due to mutation at the site. the code reports it the same as the rest for consistency, 
    # but it's counted separately in counters for further analysis
    # xref https://trello.com/c/FulYfVjT/200-stop-map-mutation-on-edge-case
    ['18112', 'AAATGTCCTTGGTGCCCGAGTCAGATCTGCCTTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGA', 'TAAGGCAGATCTGACTCGGGCACCAAGGACATTTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCG', 78, 68, [69] ],

    # relatively short target match with one bp toggle
    [ '11555', 'CTCAGTCCTTGGTGCCCGAGTCAGTGAGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTC', 'AGCTCACTGACTCGGGCACCAAGGACTGAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGG', 50, 44, [47] ],

    # matches more than one site-end-toggle combination in the target
    # https://trello.com/c/fxS3L7aM/202-stop-map-ambiguous-match-minimum-target-match-length
    [ '16610', 'AAGCGTCCTTGGTGCCCGAGTCAGTGGAGGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCT', 'ACCTCCACTGACTCGGGCACCAAGGACGCTTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTG', None, None, [] ],

]


class TestMutPairs(unittest.TestCase):
    
    def setUp(self):
        self.spats = Spats()
        self.spats.run.cotrans = True
        self.spats.run.cotrans_linker = 'CTGACTCGGGCACCAAGGAC'
        self.spats.run.count_mutations = True
        self.spats.run.allowed_target_errors = 1
        self.spats.run.adapter_b = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
        self.setup_processor()
        self.spats.addTargets("test/mut/mut_single.fa")

    def setup_processor(self):
        self.spats.run.algorithm = "find_partial"

    def tearDown(self):
        self.spats = None

    def pair_for_case(self, case):
        pair = Pair()
        pair.set_from_data(case[0], case[1], case[2])
        return pair

    def run_case(self, case):
        pair = self.pair_for_case(case)
        self.spats.process_pair(pair)
        self.assertEqual(case[4], pair.site, "res={} != {} ({}, {})".format(pair.site, case[4], self.__class__.__name__, case[0]))
        if pair.site:
            self.assertEqual(case[3], pair.end)
            self.assertEqual(case[5], pair.mutations)
        return pair

    def test_pairs(self):
        self.spats.run.pair_length = len(cases[0][1])
        for case in cases:
            self.run_case(case)
        print("Ran {} pair->site cases.".format(len(cases)))


class TestMutPairsLookup(TestMutPairs):
    
    def setup_processor(self):
        self.spats.run.algorithm = "lookup"
