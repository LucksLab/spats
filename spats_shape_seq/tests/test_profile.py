import nose
import os
import subprocess
import unittest

from spats_shape_seq.tool import SpatsTool

class TestProfile(unittest.TestCase):

    def test_profile(self):
        testdir = os.getcwd()
        try:
            with open(os.devnull, 'w') as FNULL:
                os.chdir(os.path.join(testdir, "test", "profile"))
                st = SpatsTool(os.getcwd())
                st._run(["run"])
                st._run(["dump", "run"])
                ref = open("ref_output.csv", 'r').read()
                cur = open("handmade_test.csv", 'r').read()
                self.assertEqual(ref, cur)
        finally:
            os.chdir(testdir)

