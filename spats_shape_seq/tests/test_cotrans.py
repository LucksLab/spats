import unittest

from spats_shape_seq import Spats
from spats_shape_seq.pair import Pair


# [ id, r1, r2, end, site ]
cases = [

    [ "1116:19486:8968", "TCCGGTCCTTGGTGCCCGAGTCAGTCCTTCCTCCTA", "GAGTCTATTTTTTTAGGAGGAAGGACTGACTCGGGC", 93, 68 ],
    # R2 partial linker case

    [ "1116:16151:46609", "GGGTGTCCTTGGTGCCCGAGTCAGAAAAGTTCTTCT", "TCTATGAGCAAAGGAGAAGAACTTTTCTGACTCGGG", 119, 93 ],
    # R2 partial linker case

    [ "1116:2824:48570", "GGGTGTCCTTGGTGCCCGAGTCAGGTTCTTCTCCTT", "TACTGGTAGGAGTCTATTTTTTTAGGAGGAAGGATA", None, None ],
    # 1 bp toggle in last 4 of R2 (xref v102 compat cases)

    [ "301028", "AAGTGTCCTTGGTGCCCGAGTCAGAGATAGATCGGA", "ATCTCTGACTCGGGCACCAAGGACACTTAGATCGGA", 96, 92 ],
    # adapter trim case

    [ "360389", "GAATGTCCTTGGTGCCCGAGTCAGAAAAAAATTTTT", "ATGGAGTTCGCCATAAACGCTGCTTAGCTAATGACT", None, None ],
    # R1 no-match case

    [ "683779", "TCCGGTCCTTGGTGCCCGAGTCAGAAAAAAATAGAA", "TCTATTTTTTTCTGACTCGGGCACCAAGGACCGGAA", 82, 71 ],
    # TTATAGGC....CACTACTGGTAGGAGTCTATTTTTTTAGGAGGAAGGATCTATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTC
    # r2:                        TCTATTTTTTT.CTGACTCGGGCACCAAGGAC.CGGA.A
    # LNK:                                  .CTGACTCGGGCACCAAGGAC
    # rc                         TCTATTTTTTT.CTGACTCGGGCACCAAGGAC
    # r1:                   TCCG.GTCCTTGGTGCCCGAGTCAG.AAAAAAATAGA.A
    # one bp adapter trim case

    [ "31631284", "TTCAGTCCTTGGTGCCCGAGTCAGAGATAGATCGGA", "ATCTCTGACTCGGGCACCAATGACCGGAAGATCGGA", None, None ],
    # r2: ATCTCTGACTCGGGCACCAATGAC.CGGA.AGATCGGA
    # LNK:   .CTGACTCGGGCACCAAGGAC.
    # one bp toggle in linker !

    [ "7232", "AGGTGTCCTTGGTGCCCGAGTCAGTAGCTAAGAAAT", "TTATAGGCGATGGAGTTCGCCATAAACGCTGCTTAG", None, None ],
    #     TTATAGGCGATGGAGTTCGCCATAAACGCTGCTTAGCTA.ATGACTCCTACCAGTATCACTACTGGTAGGAGTCTATTTTTTTAGGAGGAAGGATCTATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTC
    # rc                             ATTTCTTAGCTA.CTGACTCGGGCACCAAGGAC ==> match errors
    # LNK:                                       .CTGACTCGGGCACCAAGGAC
    # r1:                       AGGT.GTCCTTGGTGCCCGAGTCAGTAGCTAAGAAAT

    [ "16845404", "AAATGTCCTTGGTGCCCGAGTCAGACTGGTAGGAGT", "TCTTATAGGCGATGGAGTTCGCCATAAACGCTGCTT", None, None ],
    # R2 left of zero case

    [ "24102328", "AAGCGTCCTTGGTGCCCGAGTCAGGAGTCATAGATC", "ATGACTCCTGACTCGGGCACCAAGGACGCTTAGATC", 46, 39 ],
    # f_wt: TTATAGGCGATGGAGTTCGCCATAAACGCTGCTTAGCTAATGACTCCTACCAGTATCACTACTGGTAGGAGTCTATTTTTTTAGGAGGAAG.>GATCTATGA<.GCAAAGGAG...
    # r2:                                          ATGACTC.CTGACTCGGGCACCAAGGAC.GCTT.AGATC
    #                                                     .CTGACTCGGGCACCAAGGAC.
    # rc                                     GATCT.ATGACTC.CTGACTCGGGCACCAAGGAC
    # r1:                                     AAGC.GTCCTTGGTGCCCGAGTCAG.GAGTCAT.AGATC
    # corner case b/c longest_match on rc(r1) could match the target spot, or the spot indicated by ".>xx<." in target,
    # both are 9bp long. make sure the code matches the correct one.

    [ "51216106", "GGGTGTCCTTGGTGCCCGAGTCAGATTAGCTAAGCA", "AGCTAATCTGACTCGGGCACCAAGGACGCTGCTTAG", None, None ],
    #f_wt:  TTATAGGCGATGGAGTTCGCCATAAAC->GC.TGCTTAGCTAAT.GACTCCTACCAGTATCACTACTGGTAGGAGTCTATTTTTTTAGGAGGAAGGATCTATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTC
    #r2:                                         AGCTAAT.CTGACTCGGGCACCAAGGAC.GCTGCTTAG  ==> R2 is misaligned with R1, when it includes the linker
    #LNK:                                                CTGACTCGGGCACCAAGGAC                it should left-aligned with R1. also, note that the
    #rc                                     TGCTTAGCTAAT.CTGACTCGGGCACCAAGGAC                right side of R2 matches the target, to the right of
    #r1:                               GGGT.GTCCTTGGTGCCCGAGTCAGATTAGCTAAGCA                 the "->". clearly messed up, should be rejected.


    [ "41823514", "GAATGTCCTTGGTGCCCGAGTCAGAACTCCAAGATC", "TGGAGTTCTGACTCGGGCACCAAGGACATTCAGATC", None, None ],
    #    TTATAGGCGA.TGGAGTT.CGCC....TTCAC.TGGAGTT.GTC  ==> multiple R1 match
    #r2:            TGGAGTT.CTGACTCGGGCACCAAGGAC.ATTC.AGATC
    #LNK:                  .CTGACTCGGGCACCAAGGAC.

    [ "180", "AAGCTGTCCTTGGTGCCCGAGTCAGGAAAAGTTCTT", "TTTTTTTAGGAGGAAGGATCTATGAGCAAAGGAGAA", None, None ],
    # f_wt:  TTAT....GAGTCTATTTTTTTAGGAGGAAGGATCTATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTC
    # LNK:                                                               .CTGACTCGGGCACCAAGGAC.
    # rc                                                      AAGAACTTTTC.CTGACTCGGGCACCAAGGAC.A ==> R1 over right edge (of linker)
    # r1:                                                 AAGC.TGTCCTTGGTGCCCGAGTCAGGAAAAGTTCTT

    [ "67219", "GAGTGTCCTTGGTGCCCGAGTCAGTCGACAACTCCA", "TTATAGGCGATGGAGTTCGCCATAAACGCTGCTTAG", None, None ],
    # f_wt: TTAT...TCACTGGAGTTGTC
    # LNK:                         .CTGACTCGGGCACCAAGGAC
    # rc:              TGGAGTTGTCGA.CTGACTCGGGCACCAAGGAC  ==> R1 over right edge (of target)
    # r1: GAGT.GTCCTTGGTGCCCGAGTCAGTCGACAACTCCA                       

    [ "58726", "GGATGTCCTTGGTGCCCGAGTCAGCCTTAGATCGGA", "AAGGCTGACTCGGGCACCAAGGACATCCAGATCGGA", None, None ],
    # LNK:                   .CTGACTCGGGCACCAAGGAC
    #r2:                 AAGG.CTGACTCGGGCACCAAGGAC.ATCC.AGATCGGA
    #f_wt:   TTAT...GAGG.AAGG.ATCTATGAGCA.AAGG.AGAAGAACTTTTCACTGGAGTTGTC => multiple R1 match

    [ "188425", "GGACGTCCTTGGTGCCCGAGTCAGTATAGATCGGAA", "ATACTGACTCGGGCACCAAGGACTTCCAGATCGGAA", None, None ],
    #f_wt:    TT.ATA.GGCGATGGAGTTCGCC.ATA.AACGCTGCTTA...
    #r2:  ATA.CTGACTCGGGCACCAAGGACTTCC.AGATCGGAA
    # should fail b/c of short match and multiple R1

    [ "jjb_683779'", "TCCGGTCCTTGGTGCCCGAGTCAGAAAAAAATAGAA", "TCTATTTTTTTCTGACTCGGGCACCAAGGACCGTAA", 82, 71 ],
    # TTATAGGC....CACTACTGGTAGGAGTCTATTTTTTTAGGAGGAAGGATCTATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTC
    # r2:                        TCTATTTTTTT.CTGACTCGGGCACCAAGGAC.CGTA.A
    # r1:                   TCCG.GTCCTTGGTGCCCGAGTCAG.AAAAAAATAGA.A
    # 1-bp toggle in R2 of handle -- for now, confirms that we ignore this

    [ "jjb_L21", "GGACGTCCTTGGTGCCCGAGTCAGGGCGAACTAGAT", "AGTTCGCCCTGACTCGGGCACCAAGGACGTCCAGAT", 21, 13 ],
    # testing cotrans minimum length -- L=21 passes

    [ "jjb_L20", "GGACGTCCTTGGTGCCCGAGTCAGGCGAACTCAGAT", "GAGTTCGCCTGACTCGGGCACCAAGGACGTCCAGAT", 20, 12 ],
    # testing cotrans minimum length -- L=20 passes

    [ "jjb_L19", "GGACGTCCTTGGTGCCCGAGTCAGCGAACTCCAGAT", "GGAGTTCGCTGACTCGGGCACCAAGGACGTCCAGAT", None, None ],
    # testing cotrans minimum length -- L=19 fails

    [ "406149", "AGGTGTCCTTGGTGCCCGAGTCAGGACAACTCCAGT", "TTATAGGCGATGGAGTTCGCCATAAACGCTGCTTAG", 132, 0 ],
    # R1 on end of target sequence

    [ "70394", "AAGCGTCCTTGGTGCCCGAGTCAGTTGAGATCGGAA", "CAACTGACTCGGGCACCAAGGACCCTTAGATCGGAA", 104, 101 ],
    # very short (3bp) target match

    [ "89185", "TCCAGTCCTTGGTGCCCGAGTCAGCTAAGCAGCGTT", "AATGACTCCTACCAGTATCACTACTGGTAGGAGTCT", None, None ],
    # R2 is to the right of R1 -- invalid pair

    [ "3185000", "GAACGTCCTTGGTGCCCGAGTCAGGTTTATGGCGAA", "TCGCCATAAACCTGACTCGGGCACCAAGGACGTTCC", None, None ],
    # adapter trim corner case (R2 requires a single bp trim, but it's a 1bp toggle)

    [ "jjb_3185000'", "GAACGTCCTTGGTGCCCGAGTCAGGTTTATGGCGAA", "TCGCCATAAACCTGACTCGGGCACCAAGGACGTTCA", 27, 16 ],
    # adapter trim corner case (R2 requires a trim -- which now matches R1)
]

v102_compat_cases = [
    [ "1116:2824:48570", "GGGTGTCCTTGGTGCCCGAGTCAGGTTCTTCTCCTT", "TACTGGTAGGAGTCTATTTTTTTAGGAGGAAGGATA", 115, 59 ],
    [ "1011640", "GGACGTCCTTGGTGCCCGAGTCAGTAGCTAAGCAGC", "AACGCTGCTTAGCTACTGACTCGGGCACCAAGTACG", 39, 24 ],
]

class TestPairsPartial(unittest.TestCase):
    
    def setUp(self):
        self.spats = Spats()
        self.spats.run.cotrans = True
        self.spats.run.cotrans_linker = 'CTGACTCGGGCACCAAGGAC'
        self.setup_processor()
        self.spats.addTargets("test/cotrans/cotrans_single.fa")

    def setup_processor(self):
        self.run_compat = True
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
        return pair

    def test_pairs(self):
        for case in cases:
            self.run_case(case)
        self.spats.run._p_v102_compat = True
        if self.run_compat:
            for case in v102_compat_cases:
                self.run_case(case)
        print("Ran {} pair->site cases.".format(len(cases)))


class TestPairsLookup(TestPairsPartial):
    
    def setup_processor(self):
        self.run_compat = False
        self.spats.run.algorithm = "lookup"
