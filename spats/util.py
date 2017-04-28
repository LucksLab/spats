
import string

from config import spats_config


def _warn(stuff):
    print stuff

def _debug(stuff):
    if spats_config.debug:
        print stuff

rev_comp_complementor = string.maketrans("ATCGatcg", "TAGCtagc")

def reverse_complement(seq):
    return seq.translate(rev_comp_complementor)[::-1]

# hamming distance with tracking and shortcut-out
def string_match_errors(substr, target_str, max_errors = None):
    errors = []
    # TODO: do we need this "min" here?
    for i in range(min(len(substr), len(target_str))):
        if substr[i] != target_str[i]:
            errors.append(i)
            if max_errors and len(errors) >= max_errors:
                return errors
    return errors

class Counters(object):

    def __init__(self):
        self.reset()

    def reset(self):
        self._counts = {}

    def __getattr__(self, key):
        if key == "_counts":
            return super.__getattr__(self, key)
        return self._counts.get(key, 0)

    def __setattr__(self, key, value):
        if key == "_counts":
            super.__setattr__(self, key, value)
        self._counts[key] = value

    def counts_dict(self):
        return { key : value for key, value in self._counts.iteritems() if key != "_counts" }
