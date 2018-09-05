import sys
import unittest
import simplejson     # Use so don't have to bother with unicode->str 
from spats_shape_seq import Spats
from spats_shape_seq.pair import Pair



class TestHarness:
    class SpatsTestResults(unittest.TextTestResult):
        # Taken from legacy python unittest
        class WritelnDecorator:
            """Used to decorate file-like objects with a handy 'writeln' method"""
            def __init__(self, stream):
                self.stream = stream
            def __getattr__(self, attr):
                return getattr(self.stream,attr)
            def writeln(self, arg = None):
                if arg: self.write(arg)
                self.write('\n') # text-mode streams translate to \r\n if needed
        def __init__(self, outer):
            super(TestHarness.SpatsTestResults, self).__init__(TestHarness.SpatsTestResults.WritelnDecorator(sys.stdout), True, 1)
        # TAI:  can override addFailure(self, test, err) and # addError(self, test, err) here if we want 

    def __init__(self, testfile = "tests.json"):
        self.test_results = self.SpatsTestResults(self)
        with open(testfile) as TF:
            self.tests = simplejson.load(TF)

    def run_testsets(self):
        for testset_dict in self.tests["tests"]:
            self.TestSet(self, testset_dict).run(self.test_results)   # TODO:  use a unittest.TestResult object

    def print_stats(self):
        print("\n\nSUMMARY")
        print("Tests Sets Run:  {}".format(self.test_results.testsRun))
        print("All Test Sets Passed?  {}".format(self.test_results.wasSuccessful()))
        if not self.test_results.wasSuccessful():
            self.test_results.printErrors()


    # TAI:  Use a unittest.TestSuite object
    class TestSet(unittest.TestCase):
        def __init__(self, outer, testset_dict):
            super(TestHarness.TestSet, self).__init__()
            self.outer = outer
            self.testset_dict = testset_dict
            self.set_name = self.testset_dict['set_name']
            self.algorithms = [ "find_partial" ]

        def setUp(self):
            #print("setting up options for test set '{}'...".format(self.set_name))
            self.spats = Spats()
            for key, value in self.outer.tests['default_opts'].iteritems():
                if (key == 'algorithms'):
                    self.algorithms = value
                else:
                    setattr(self.spats.run, key, value)
            for key, value in self.testset_dict.get('run_opts', {}).iteritems():
                if (key == 'algorithms'):
                    self.algorithms = value
                else:
                    setattr(self.spats.run, key, value)
            all_targets = self.outer.tests['targets']
            for target in self.testset_dict['targets']:
                self.spats.addTarget(target, all_targets[target])

        def tearDown(self):
            self.spats = None

        def runTest(self):
            for algorithm in self.algorithms:
                #print("running tests for test set '{}' with algorithm '{}'...".format(self.set_name, algorithm))
                self.spats.run.algorithm = algorithm
                for case in self.testset_dict['tests']:
                    try:
                        self._run_case(case)
                    except Exception as e:
                        # TODO:  This isn't right.  If we catch here, it wont update the test results.  (but if we don't, it aborts in the middle of the test cases for this set.  --> need to use a suite.
                        print("exception caught on test '{}': {}".format(case['id'], e))

        def _run_case(self, case):
            #print("running case {} ({})...".format(case['id'], case.get('comment', '')))
            pair = self._pair_for_case(case)
            self.spats.counters.reset()
            self.spats.process_pair(pair)
            self._check_expects(case['expect'], pair, case['id'])

        def _pair_for_case(self, case):
            pair = Pair()
            pair.set_from_data(case['id'], case['r1'], case['r2'])
            if 'r1_quality' in case:
                pair.r1.quality = case['r1_quality']
            else:
                pair.r1.quality = 'K' * len(case['r1'])
            if 'r2_quality' in case:
                pair.r2.quality = case['r2_quality']
            else:
                pair.r2.quality = 'K' * len(case['r2'])
            return pair

        def _check_expects(self, expects, pair, caseid):
            msg = "testset='{}', test id='{}' failed:  ".format(self.set_name, caseid)
            if expects['site'] is None:
                self.assertIs(pair.site, None, msg + "pair.site={} when expecting none.".format(pair.site))
            else:
                self.assertIsNot(pair.site, None, msg + "pair.site is none when expecting {}.".format(expects['site']))
                self.assertEqual(expects['site'], pair.site, msg + "pair.site={} != expect.site={}".format(pair.site, expects['site']))
                if 'end' in expects:
                    self.assertEqual(expects['end'], pair.end, msg + "pair.end={} != expect.end={}".format(pair.end, expects['end']))
            if 'muts' in expects:
                if expects['muts'] is not None  and  len(expects['muts']) > 0:
                    self.assertEqual(sorted(expects['muts']), sorted(pair.mutations) if pair.mutations else pair.mutations, msg + "mismatching mutations:  expected={}, pair.mutations={}".format(expects['muts'], pair.mutations))
                else:
                    self.assertTrue(pair.mutations is None  or len(pair.mutations) == 0, msg + "unexpected mutations: {}".format(pair.mutations))
            if 'counters' in expects:
                for counter, value in expects['counters'].iteritems():
                    self.assertEqual(getattr(self.spats.counters, counter), value, msg + "counter '{}' value off: expected={} != got={}".format(counter, value, getattr(self.spats.counters, counter)))
            # TODO:  not sure about this one...
            if 'pair.target' in expects:
                self.assertEqual(pair.target, expects['pair.target'], msg + "pair.target={} != expect.pair.target={}".format(pair.target, expects['pair.target']))



if __name__ == "__main__":
    th = TestHarness()
    th.run_testsets()
    th.print_stats()
