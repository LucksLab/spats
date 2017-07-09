
import string

_is_debug = False

def _set_debug(dbg):
    global _is_debug
    _is_debug = dbg

def _warn(stuff):
    spats_config.log.write(str(stuff) + "\n")

def _debug(stuff):
    if _is_debug:
        spats_config.log.write(str(stuff) + "\n")


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
