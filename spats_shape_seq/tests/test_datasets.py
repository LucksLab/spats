import os
import shutil
import unittest

from spats_shape_seq import Spats
from spats_shape_seq.db import PairDB

cases = [ "5s" , "cotrans" ] #, "panel_RNAs" ]
algorithms = [ "find_partial", "lookup"  ]

class TestDatasets(unittest.TestCase):
    
    def test_datasets(self):
        for case in cases:
            for alg in algorithms:
                self.run_dataset(case, alg)
        print "Ran {} datasets.".format(len(cases))

    def run_dataset(self, case, algorithm):
        bp = "test/{}/".format(case)
        test_file = bp + "test.spats.tmp"
        try:
            shutil.copyfile(bp + "ds.spats", test_file)
            db = PairDB(test_file)
            s = Spats()
            db.load_run(s.run)
            s.run.writeback_results = True
            s.run.result_set_name = "test"
            s.run.algorithm = algorithm
            s.run.quiet = True
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
        


