import unittest

from spats_shape_seq import Spats
from spats_shape_seq.pair import Pair

# [ id, r1, r2, end, site, muts ]
cases = [

    [ 'x1', 'TTCAGTCCTTGGTGCCCGAGTCAGCCAGCTACATCCCGGCACACGCGTCATCTGCCTTGGCTGCTTCCTTCCGGA', 'AGGTCAGATCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCTGACTCGGGCACCAA', 103, 44, [51] ],

    [ 'x2', 'AGGCGTCCTTGGTGCCCGAGTCAGCCTTGGCTGCTTCCTTCCGGACCTGACCTGGTAAACAGAGTAGCGTTGCGG', 'ATCGGGGGCTCTGTTGGTTCCCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGTCT', 73, 0, [72], ],

    # note: on this one, the quality score on R1 vs. R2 indicates that mut 72 really is a misread on R2, so this case may change in the future...
    [ 'x3', 'GAATGTCCTTGGTGCCCGAGTCAGGACACGCGTCATCTGCCTTGGCTGCTTCCTTCCGGACCTGACCTGGTAAAC', 'ATCGGGGGCTCTGTTGGTTCCCCCGCAACGCTACTCTGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGTCA', 88, 0, [72, 87] ],

    [ 'x4', 'GAGCGTCCTTGGTGCCCGAGTCAGATGCCGACCCGGGTGGGGGCCCTGCCAGCTACATCCCGGCACACGCGTCAT', 'TAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCAGATGACGCGTGTGCCGGGATGTAGCTGGCAGGGCCCCCACCCG', 127, 43, [43] ],

    [ 'x5', 'GAATGTCCTTGGTGCCCGAGTCAGTCCTTGGTGCCCGAGTCAGTCCTTGGTTCCCGAGTCACTCCTTTGTTCCCC', 'AGGACTGACTCGGGCACCAAGGACTTTCTCGTTCACCTATTTCTTTCTCTTCCCCCTTTTTCTTTCTCTTTCTCC', None, None, [] ],

    # this one is slightly ambiguous due to mutation at the site. we count it for the site and a mutation, for now.
    ['18112', 'AAATGTCCTTGGTGCCCGAGTCAGATCTGCCTTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGA', 'TAAGGCAGATCTGACTCGGGCACCAAGGACATTTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCG', 78, 68, [68] ],

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
