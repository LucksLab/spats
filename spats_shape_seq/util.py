
import string

_debug_run = None

def _set_debug(run):
    global _debug_run
    _debug_run = run

def _warn(stuff):
    _debug_run.log.write(str(stuff) + "\n")

def _debug(stuff):
    if _debug_run.debug:
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
