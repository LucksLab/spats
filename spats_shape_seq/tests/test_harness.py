import json
import os
import sys
import unittest

from collections import defaultdict
from spats_shape_seq import Spats
from spats_shape_seq.pair import Pair



class TestHarness:
    def __init__(self, testfile = None):
        testfile = testfile or os.path.join(os.path.dirname(__file__), "tests.json")
        self.test_results = self.SpatsTestResults(self)
        self.tests = json.loads(open(testfile, 'r').read())
        self.all_targets = { str(key) : str(val) for key, val in self.tests['targets'].iteritems() }

    def run_testsets(self):
        for testset_dict in self.tests["tests"]:
            self.test_results.current_testset = str(testset_dict['set_name'])
            TestHarness.SpatsTestSet(self, testset_dict).run(self.test_results)
            self.test_results.test_sets_run += 1
        self.print_stats()
        if not self.test_results.wasSuccessful():
            raise Exception('Test harness failure')

    def print_stats(self):
        self.test_results.print_stats()


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

        def __init__(self, outer, verbose = False):
            super(TestHarness.SpatsTestResults, self).__init__(TestHarness.SpatsTestResults.WritelnDecorator(sys.stdout), True, 2 if verbose else 1)
            self.test_sets_run = 0
            self.tests_passed = 0
            self.current_testset = ""
            self.testset_success = defaultdict(int)
            self.testset_failures = defaultdict(int)
            self.testset_errors = defaultdict(int)
            self.failure_cases = []

        def addFailure(self, test, err):
            super(TestHarness.SpatsTestResults, self).addFailure(test, err)
            print(' ==> FAIL: {}'.format(err[1]))
            self.testset_failures[self.current_testset] += 1
            self.failure_cases.append('{}/{}'.format(self.current_testset, str(test.case_dict['id'])))
            
        def addError(self, test, err):
            import traceback
            super(TestHarness.SpatsTestResults, self).addError(test, err)
            print(' ==> ERROR: {}/{} {}\n{}'.format(self.current_testset, str(test.case_dict['id']), err[1], traceback.format_exc(err[2])))
            self.testset_errors[self.current_testset] += 1
            self.failure_cases.append('{}/{}'.format(self.current_testset, str(test.case_dict['id'])))

        def addSuccess(self, test):
            self.testset_success[self.current_testset] += 1
            self.tests_passed += 1

        def print_stats(self):
            print("\n\nSUMMARY")
            print("Tests Sets Run:  {}".format(self.test_sets_run))
            print("Total Tests Run:  {}".format(self.testsRun))
            print("Total Tests Passed:  {}".format(self.tests_passed))
            print("All Test Sets Passed?  {}".format(self.wasSuccessful()))
            if not self.wasSuccessful():
                if len(self.testset_failures) > 0:
                    print("Test Sets with failures ({}):".format(sum(self.testset_failures.values())))
                    print(str(dict(self.testset_failures)))
                if len(self.testset_errors) > 0:
                    print("Test Sets with errors ({}):".format(sum(self.testset_errors.values())))
                    print(str(dict(self.testset_errors)))


    class SpatsTestSet(unittest.TestSuite):
        def __init__(self, outer, testset_dict):
            super(TestHarness.SpatsTestSet, self).__init__()
            self.outer = outer
            self.testset_dict = testset_dict
            self.set_name = self.testset_dict['set_name']
            self.algorithms = [ "find_partial" ]
            self._add_all_testcases()

        def _add_all_testcases(self):
            for case in self.testset_dict['tests']:
                self.addTest(TestHarness.SpatsTestCase(self, case))

        def spats_setUp(self, spatso):
            self.algorithms = [ "find_partial", "lookup" ]
            for key, value in self.testset_dict.get('run_opts', {}).iteritems():
                if (key == 'algorithms'):
                    self.algorithms = value
                else:
                    setattr(spatso.run, key, value)
            for target in self.testset_dict['targets']:
                spatso.addTarget(target, self.outer.all_targets[target])


    class SpatsTestCase(unittest.TestCase):
        def __init__(self, test_set, case_dict):
            super(TestHarness.SpatsTestCase, self).__init__()
            self.test_set = test_set
            self.case_dict = case_dict

        def setUp(self):
            try:
                self.spats = Spats()
                self.test_set.spats_setUp(self.spats)
            except Exception as e:
                print("exception caught on testset '{}' setup : {}".format(self.test_set.set_name, e))

        def tearDown(self):
            self.spats = None

        def runTest(self):
            for algorithm in self.test_set.algorithms:
                self.spats.run.algorithm = algorithm
                self._run_case(self.case_dict)

        def _run_case(self, case):
            pair = self._pair_for_case(case)
            self.spats.counters.reset()
            self.spats.process_pair(pair)
            self._check_expects(case['expect'], pair, case['id'])

        def _pair_for_case(self, case):
            pair = Pair()
            pair.set_from_data(str(case['id']), str(case['r1']), str(case['r2']))
            if 'r1_quality' in case:
                pair.r1.quality = str(case['r1_quality'])
            if 'r2_quality' in case:
                pair.r2.quality = str(case['r2_quality'])
            return pair

        def _check_expects(self, expects, pair, caseid):
            msg = "testset='{}', test id='{}', algorithm='{}' failed: ".format(self.test_set.set_name, caseid, self.spats.run.algorithm)
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


def run_harness():
    th = TestHarness()
    th.run_testsets()
    th.print_stats()
