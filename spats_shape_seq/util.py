
import string

_debug_run = None

def _set_debug(run):
    global _debug_run
    _debug_run = run

def _warn(stuff):
    if _debug_run:
        _debug_run.log.write(str(stuff) + "\n")

def _debug(stuff):
    if _debug_run and _debug_run.debug:
        _debug_run.log.write(str(stuff) + "\n")


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

def _slow_string_find_errors(substr, target_str, max_errors = 2, max_indices = 2):
    assert(max_errors > 1)
    assert(max_indices > 0)
    if max_errors >= len(substr):
        raise Exception("more errors allowed than len of substr being sought!")
    result = []
    ls = len(substr)
    for t in range(len(target_str) - ls + 1):
        num_errors = 0
        s = 0
        while s < ls:
            if target_str[t + s] != substr[s]:
                num_errors += 1
                if num_errors > max_errors:
                    s = ls
            s += 1
        if s == ls:
            result.append(t)
            if len(result) >= max_indices:
                break
    return result

def string_find_errors(substr, target_str, max_errors = 0, max_indices = 2):
    assert(max_errors >= 0)
    assert(max_indices > 0)    # TAI:  consider using max_indices == -1 to return all
    me = min(max_errors, 1)

    half_sublen = (len(substr) >> 1)
    first_half = substr[:half_sublen]
    second_half = substr[half_sublen:]

    candidates = []
    index = 0
    while len(candidates) < max_indices:
        index = target_str.find(first_half, index)
        if -1 == index or index + len(substr) > len(target_str):
            break
        errors = string_match_errors(second_half, target_str[index + len(first_half):], me + 1)
        if len(errors) <= me:
            candidates.append(index)
        index = index + 1

    index = 0
    while len(candidates) < max_indices:
        index = target_str.find(second_half, index)
        if -1 == index:
            break
        if index - len(first_half) >= 0:
            errors = string_match_errors(first_half, target_str[index - len(first_half):], me + 1)
            if len(errors) <= me:
                candidates.append(index - len(first_half))
        index = index + 1

    result = sorted(list(set(candidates)))

    if len(result) < max_indices  and  max_errors > 1:
        result += _slow_string_find_errors(substr, target_str, max_errors, max_indices - len(result))

    return result

class Colors(object):

    def __init__(self):
        self.green = [0,1,0]
        self.black = [0,0,0]
        self.cyan = [0,1,1]
        self.red = [1,0,0]
        self.blue = [0,0,1]

# to allow setting of arbitrary dot properties
class SimpleObject(object):
    pass
