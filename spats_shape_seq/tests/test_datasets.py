import nose
import os
import shutil
import unittest

from spats_shape_seq import Spats
from spats_shape_seq.db import PairDB
from spats_shape_seq.pair import Pair
from spats_shape_seq.processor import Failures

cases = [ "5s", "cotrans" ] #, "panel_RNAs" ]
algorithms = [ "find_partial", "lookup", "native" ]

class TestDatasets(unittest.TestCase):
    
    def test_datasets(self):
        for case in cases:
            for alg in algorithms:
                if os.environ.get('SKIP_SLOW_TESTS') and alg == 'native':
                    raise nose.SkipTest('skipping slow tests')
                self.run_dataset(case, alg)
        print("Ran {} datasets.".format(len(cases)))

    def run_dataset(self, case, algorithm):
        bp = "test/{}/".format(case)
        test_file = bp + "test.spats.tmp"
        try:
            shutil.copyfile(bp + "ds.spats", test_file)
            db = PairDB(test_file)
            s = Spats()
            db.load_run(s.run)
            if not s.run.cotrans and algorithm == "native":
                return
            s.run.writeback_results = True
            s.run.result_set_name = "test"
            s.run.algorithm = algorithm
            s.run.quiet = True
            s.loadTargets(db)
            if not s._processor.exists():
                # just ignore the native test if it's not available
                self.assertEqual("native", algorithm)
                return
            s.process_pair_db(db, batch_size = 1024) # small batch_size just to exercise multiprocessing code
            msg = None
            count = 0
            for res in db.differing_results("test", "test_validation"):
                msg = str([str(x) for x in res])
                count += 1
            self.assertEqual(0, count, "{} differing results: {} / {} \n{}".format(count, case, algorithm, msg))
        finally:
            if os.path.exists(test_file):
                os.remove(test_file)
        


overlap_cases = [
    [ "o1", "GGGCGTCCTTGGTGCCCGAGTCAGAAGATCGGAAGA", "ACTGACTCGGGCACCAAGGACGCCCAGATCGGAAGA", None, Failures.r1_r2_overlap ],
    [ "o2", "GGGCGTCCTTGGTGCCCGAGTCAGAAGATCGGAAGA", "TCTGACTCGGGCACCAAGGACGCCCAGATCGGAAGA", 116, None ],
]

class TestOverlap(unittest.TestCase):

    def tearDown(self):
        self.spats = None

    def pair_for_case(self, case):
        pair = Pair()
        pair.set_from_data(case[0], case[1], case[2])
        return pair

    def run_case(self, case):
        pair = self.pair_for_case(case)
        print('running: {} / {}'.format(case[0], self.spats.run.algorithm))
        self.spats.process_pair(pair)
        self.assertEqual(case[3], pair.site, "site res={} != {} ({}, {}, {}, {})".format(pair.site, case[3], self.__class__.__name__, case[0], self.spats.run.algorithm, pair.failure))
        self.assertEqual(case[4], pair.failure, "failure res={} != {} ({}, {}, {})".format(pair.failure, case[4], self.__class__.__name__, case[0], self.spats.run.algorithm))

    def test_pairs(self):
        for alg in algorithms:
            if alg == 'native':
                continue
            self.run_algorithm(alg)

    def run_algorithm(self, alg):
        from spats_shape_seq import Spats
        self.spats = Spats()
        self.spats.run.algorithm = alg
        self.spats.run.count_mutations = True
        self.spats.run.allowed_target_errors = True
        self.spats.addTargets("test/SRP/SRP.fa")
        self.run_pairs()

    def run_pairs(self):
        for case in overlap_cases:
            self.run_case(case)
        print("Ran {} prefix test cases.".format(len(cases)))
