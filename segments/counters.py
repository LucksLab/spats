import re


def _dictIncr(d, key, m = 1):
    v = d.get(key, 0) + m
    d[key] = v
    return v

class Counters(object):

    def __init__(self, run = None):
        self._run = run
        self.reset()

    def reset(self):
        self._counts = {}

    def __getattr__(self, key):
        if key.startswith('_'):
            return super.__getattr__(self, key)
        else:
            return self._counts.get(key, 0)

    def __setattr__(self, key, value):
        if key.startswith('_'):
            super.__setattr__(self, key, value)
        else:
            self._counts[key] = value

    def countsDict(self):
        return { key : value for key, value in self._counts.items() if not key.startswith('_') }

    def incrementKey(self, counterKey, multiplicity = 1):
        return _dictIncr(self._counts, counterKey, multiplicity)
    
    def getKey(self, counterKey):
        return self._counts.get(counterKey, 0)

    def setKey(self, counterKey, value):
        self._counts[counterKey] = value

    def report(self, skip = None):
        counters = self
        total = counters.total_pairs
        warn_keys = [ "multiple_R1_match", ]
        skip_keys = skip or [ "_barcode:", "ins:", "del:", "mut:" ]
        skipped_some = False
        countinfo = self.countsDict()
        for key in sorted(countinfo.keys(), key = lambda k : "{:04d}:{}".format(1000000 - countinfo[k], k)):
            skipped = False
            for s in skip_keys:
                if s in key:
                    skipped = True
                    break
            if skipped:
                skipped_some = True
                continue
            print("  {}{} : {} ({:.1f}%)".format("*** " if key in warn_keys else "", key, countinfo[key], 100.0 * (float(countinfo[key])/float(total)) if total else 0))
        if 0 and skipped_some:
            print("[INFO] Some counters not printed above.")
