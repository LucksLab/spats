import json
import os
import sys
import unittest

from collections import defaultdict
from spats_shape_seq import Spats
from spats_shape_seq.pair import Pair


class SpatsCase(object):

    def __init__(self, jsonDict):
        self.id = str(jsonDict['id'])
        self.r1 = str(jsonDict['r1'])
        self.r2 = str(jsonDict['r2'])
        self.r1_quality = str(jsonDict.get('r1_quality', '')) or None
        self.r2_quality = str(jsonDict.get('r2_quality', '')) or None
        self.comment = jsonDict.get('comment')
        self.expect = jsonDict.get('expect', {})
        self.run_opts = jsonDict.get('run_opts', {})
        self.targets = jsonDict.get('targets', {})
        self.set_name = jsonDict.get('set_name')

    def clone(self):
        return SpatsCase(self.jsonDict())

    def jsonDict(self):
        info = self.caseDict()
        if self.run_opts:
            info["run_opts"] = { k : v for k, v in self.run_opts.iteritems() }
        if self.targets:
            info["targets"] = { k : v for k, v in self.targets.iteritems() }
        if self.set_name:
            info["set_name"] = self.set_name
        if self.comment:
            info["comment"] = self.comment
        return info

    def caseDict(self):
        info = { "id" : self.id, "r1" : self.r1, "r2" : self.r2 }
        if self.r1_quality:
            info["r1_quality"] = self.r1_quality
            info["r2_quality"] = self.r2_quality
        if self.expect:
            info["expect"] = { k : v for k, v in self.expect.iteritems() }
        return info

    def pair(self):
        pair = Pair()
        pair.set_from_data(self.id, self.r1, self.r2)
        if self.r1_quality:
            pair.r1.quality = self.r1_quality
            pair.r2.quality = self.r2_quality
        return pair


class SpatsCaseSet(object):

    def __init__(self, reg, jsonDict):
        self.reg = reg
        self.name = str(jsonDict['set_name'])
        self.tests = [ SpatsCase(d) for d in jsonDict['tests'] ]
        self.targets = jsonDict['targets']
        self.run_opts = jsonDict.get('run_opts', {})

    def extract_case(self, case_id):
        for case in self.tests:
            if case.id == case_id:
                clone = case.clone()
                for k, v in self.run_opts.iteritems():
                    if k not in clone.run_opts:
                        clone.run_opts[k] = v
                for t in self.targets:
                    if t not in clone.targets:
                        clone.targets[t] = self.reg.targets[t]
                clone.set_name = self.name
                return clone
        return None


class SpatsCaseRegistry(object):

    def __init__(self, path = None):
        self.testfile = path or os.path.join(os.path.dirname(__file__), "tests.json")
        self.tests = json.loads(open(self.testfile, 'r').read())
        self.testsets = [ SpatsCaseSet(self, d) for d in self.tests["tests"] ]
        self.targets = { str(key) : str(val) for key, val in self.tests["targets"].iteritems() }

    def extract_case(self, case_id):
        for testset in self.testsets:
            case = testset.extract_case(case_id)
            if case:
                return case
        return None

    def add_case(self, case_json):
        for key in [ 'id', 'r1', 'r2', 'run_opts', 'set_name', 'targets' ]:
            if key not in case_json:
                raise Exception("Test case missing required key '{}'".format(key))
        case = SpatsCase(case_json)
        testset = None
        for ts in self.tests["tests"]:
            if case.set_name == ts["set_name"]:
                testset = ts
                break
        if testset:
            if case.run_opts != testset['run_opts']:
                raise Exception('Case options do not match set options {}'.format(testset('run_opts')))
            if case.targets.keys() != testset['targets']:
                raise Exception('Case targets do not match set targets {}'.format(testset('targets')))
            testset['tests'].append(case.caseDict())
        else:
            print("Creating new test_set '{}'...".format(case.set_name))
            self.tests["tests"].append({ "set_name" : case.set_name,
                                         "run_opts" : case.run_opts,
                                         "targets" : case.targets.keys(),
                                         "tests" : [ case.caseDict() ] })

        open(self.testfile, 'w').write(json.dumps(self.tests, sort_keys = True, indent = 4, separators = (',', ': ')))


def registry():
    return SpatsCaseRegistry()


class TestHarness:

    def __init__(self, testfile = None):
        self.registry = SpatsCaseRegistry(testfile)
        self.test_results = self.SpatsTestResults(self)

    def run_testsets(self):
        for testset in self.registry.testsets:
            self.test_results.current_testset = testset.name
            TestHarness.SpatsTestSet(self, testset).run(self.test_results)
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
            self.failure_cases.append(test.case.id)
            
        def addError(self, test, err):
            import traceback
            super(TestHarness.SpatsTestResults, self).addError(test, err)
            print(' ==> ERROR: {}/{} {}\n{}'.format(self.current_testset, test.case.id, err[1], traceback.format_exc(err[2])))
            self.testset_errors[self.current_testset] += 1
            self.failure_cases.append(test.case.id)

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
                    print(dict(self.testset_failures))
                if len(self.testset_errors) > 0:
                    print("Test Sets with errors ({}):".format(sum(self.testset_errors.values())))
                    print(dict(self.testset_errors))
                print("Failed tests:")
                for f in self.failure_cases:
                    print("  make tc.{}".format(f))

    class SpatsTestSet(unittest.TestSuite):

        def __init__(self, outer, testset):
            super(TestHarness.SpatsTestSet, self).__init__()
            self.outer = outer
            self.testset = testset
            self.algorithms = [ "find_partial", "lookup" ]
            self._add_all_testcases()

        @property
        def name(self):
            return self.testset.name

        def _add_all_testcases(self):
            for case in self.testset.tests:
                self.addTest(TestHarness.SpatsTestCase(self, case))

        def spats_setUp(self, spatso):
            self.algorithms = [ "find_partial", "lookup" ]
            for key, value in self.testset.run_opts.iteritems():
                if key == 'algorithms':
                    self.algorithms = value
                else:
                    setattr(spatso.run, key, value)
            for target in self.testset.targets:
                spatso.addTarget(target, self.outer.registry.targets[target])


    class SpatsTestCase(unittest.TestCase):
        def __init__(self, test_set, case):
            super(TestHarness.SpatsTestCase, self).__init__()
            self.test_set = test_set
            self.case = case

        def setUp(self):
            try:
                self.spats = Spats()
                self.test_set.spats_setUp(self.spats)
            except Exception as e:
                print("exception caught on testset '{}' setup : {}".format(self.test_set.name, e))
                raise e

        def tearDown(self):
            self.spats = None

        def runTest(self):
            for algorithm in self.test_set.algorithms:
                self.spats.run.algorithm = algorithm
                self._run_case(self.case)

        def _run_case(self, case):
            pair = case.pair()
            self.spats.counters.reset()
            self.spats.process_pair(pair)
            self._check_expects(case, pair)

        def _check_expects(self, case, pair):
            expects = case.expect
            msg = "testset='{}', test id='{}', algorithm='{}' failed: ".format(self.test_set.name, case.id, self.spats.run.algorithm)
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
                    self.assertEqual(getattr(self.spats.counters, str(counter)), value, msg + "counter '{}' value off: expected={} != got={}".format(counter, value, getattr(self.spats.counters, counter)))
            if 'pair.target' in expects:
                tname = pair.target.name if pair.target else None
                self.assertEqual(tname, expects['pair.target'], msg + "pair.target={} != expect.pair.target={}".format(tname, expects['pair.target']))


def run_harness():
    th = TestHarness()
    th.run_testsets()
    th.print_stats()
