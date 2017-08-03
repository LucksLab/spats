import ast
import os
import subprocess
from processor import PairProcessor, Failures

# wrapper for unit-testing the native code. extremely slow for anything but the smallest datasets.
class CotransNativeProcessor(PairProcessor):

    def prepare(self):
        pass

    def exists(self):
        import spats_shape_seq
        bin_path = os.path.join(os.path.dirname(spats_shape_seq.__file__), "../native/bin/unittest")
        return os.path.exists(bin_path)

    def process_pair(self, pair):
        try:
            import spats_shape_seq
            basepath = os.path.dirname(spats_shape_seq.__file__)
            bin_path = os.path.join(basepath, "../native/bin/unittest")
            tgt_path = os.path.join(basepath, "../test/cotrans/cotrans_single.fa")
            res = subprocess.check_output([bin_path, tgt_path, pair.r1.original_seq, pair.r2.original_seq])
            val = ast.literal_eval(res)
        except:
            import traceback
            print(traceback.format_exc())
            pair.failure = Failures.nomatch
            return

        if not self._match_mask(pair):
            return

        pair.end = val[1]
        pair.target = self._targets.targets[0]
        pair.site = val[2]
        self.counters.register_count(pair)

