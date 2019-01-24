
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

def max_element(l, zero_based = True):
    ''' @return pair containing the index and value respectively of the maximum element of list l '''
    return max(enumerate(l, 0 if zero_based else 1), key=lambda x:x[1])

def min_element(l, zero_based = True):
    ''' @return pair containing the index and value respectively of the minimum element of list l '''
    return min(enumerate(l, 0 if zero_based else 1), key=lambda x:x[1])

rev_comp_complementor = string.maketrans("ATCGatcg", "TAGCtagc")

def reverse_complement(seq):
    return seq.translate(rev_comp_complementor)[::-1]

# hamming distance with tracking and shortcut-out
def string_match_errors(substr, target_str, max_errors = None):
    errors = []
    # TODO: do we need this "min" here?
    for i in xrange(min(len(substr), len(target_str))):
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
    for t in xrange(len(target_str) - ls + 1):
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

def string_edit_distance(s1, s2, substitution_cost = 2, insert_delete_cost = 1):
    '''
     Standard Levenshtein string edit distance if substitution cost is either 1x or 2x times insertion/deletion cost.
     The algorithm runs in O(M*N) time (where M is the legnth of s1 and N is the length of S2) and O(max(M,N)) space.
      @param  s1                     First string in the comparison
      @param  s2                     Second string in the comparison
      @param  substitution_cost      The amount to penalize a character substitution (usually 1 or 2)
      @param  insert_delete_cost     The amount to penalize a character insertion or deletion (usually 1)
      @return tuple of:
          - the cost of the character operations (substitution, insertion, deletion) to convert s1 into s2
          - maximum edit distance possible between s1 and s2
    '''
    m = len(s1) + 1
    n = len(s2) + 1
    max_possible = substitution_cost * min(m, n) + insert_delete_cost * (max(m, n) - min(m, n))
    if 1 == max_possible:
        assert(False)
        return 0, max_possible
    prev_row = [j for j in xrange(n)]
    cur_row = [0] * n
    for i in xrange(1, m):
        cur_row[0] = i
        for j in xrange(1, n):
            cur_row[j] = prev_row[j-1] if s1[i-1] == s2[j-1] else min(insert_delete_cost + min(prev_row[j], cur_row[j-1]), substitution_cost + prev_row[j-1])
        cur_row, prev_row = prev_row, cur_row
    return prev_row[n-1], max_possible


def string_edit_distance2(s1, s2, substitution_cost = 2, insert_delete_cost = 1):
    '''
     Standard Levenshtein string edit distance if substitution cost is either 1x or 2x times insertion/deletion cost.
     This version keeps track of the operations required such that it can return the number of consecutive edits result.
     The algorithm runs in O(M*N) time (where M is the legnth of s1 and N is the length of S2) and O(M*N) space.
      @param  s1                     First string in the comparison
      @param  s2                     Second string in the comparison
      @param  substitution_cost      The amount to penalize a character substitution (usually 1 or 2)
      @param  insert_delete_cost     The amount to penalize a character insertion or deletion (usually 1)
      @return tuple of:
          - the cost of the character operations (substitution, insertion, deletion) to convert s1 into s2
          - maximum edit distance possible between s1 and s2
          - the number of *regions* in the two strings that are different
    '''
    m = len(s1) + 1
    n = len(s2) + 1
    max_possible = substitution_cost * min(m, n) + insert_delete_cost * (max(m, n) - min(m, n))
    if 1 == max_possible:
        assert(False)
        return 0, max_possible, 0

    prev_row = [j for j in xrange(n)]
    cur_row = [0] * n
    B = [[0 for j in xrange(n)] for i in xrange(m)]
    for j in xrange(n):
        B[0][j] = 2

    for i in xrange(1, m):
        cur_row[0] = i
        B[i][0] = 1
        for j in xrange(1, n):
            if s1[i-1] == s2[j-1]:
                cur_row[j] = prev_row[j-1]
            else:
                B[i][j], cur_row[j] = min_element([ insert_delete_cost + prev_row[j], insert_delete_cost + cur_row[j-1], substitution_cost + prev_row[j-1] ], False)
        cur_row, prev_row = prev_row, cur_row

    num_consecutive_edits = 0
    i = m - 1
    j = n - 1
    inOp = False
    while i > 0  or  j > 0:
        if B[i][j] == 0:
            i -= 1
            j -= 1
            inOp = False
        else:
            if not inOp:
                num_consecutive_edits += 1
                inOp = True
            if B[i][j] == 1:
                i -= 1
            elif B[i][j] == 2:
                j -= 1
            elif B[i][j] == 3:
                i -= 1
                j -= 1
            else:
                assert(False)
                i -= 1

    return prev_row[n-1], max_possible, num_consecutive_edits


class Indel:
    def __init__(self, insert_type, length = 1):
        self.insert_type = insert_type    # True for inserts, False for deletes
        self.length = length

def find_indels(source, target, match_value = 2, mismatch_cost = 2, gap_open_cost = 6, gap_extend_cost = 1):
    """
     Find the indels (insertions and deletions) in a source string relative to a target string.
     First, heuristically tries to find the best alignment of two strings of length M and N, 
     allowing for insertions and deletions using the affine-gap Smith-Waterman algorithm.
     Once an alignment is found, returns the indel mutations with respect to the target.
     Deletions are indexed at the last item deleted in target.
     Insertions are indexed at the first spot in the target to be moved.
     Warning:  This function has not been optimized; use with care on long strings.

     @param  source              Source string
     @param  target              Target string
     @param  match_value         Score bump/reward for each matching spot
     @param  mismatch_cost       Penalty for each mismatching spot (should be non-negative)
     @param  gap_open_cost       Penalty for opening a gap (should be non-negative)
     @param  gap_extend_cost     Penalty for extending an already-open gap (should be non-negative)
     @return a dictionary mapping indices in the target string to Indel objects
    """
    m = len(source) + 1
    n = len(target) + 1

    H = [[0.0 for c in xrange(n)] for r in xrange(m)]
    P = [[(r-1, c-1) for c in xrange(n)] for r in xrange(m)]

    maxH = 0.0
    maxs = (0, 0)
    choices = [0.0] * 4
    for i in xrange(1, m):
        for j in xrange(1, n):
            choices[1] = H[i-1][j-1] + (match_value if source[i-1] == target[j-1] else -mismatch_cost)
            ## TAI:  The next two lines should be considered for # optimization...
            ki, choices[2] = max_element([(H[k][j] - (gap_open_cost + gap_extend_cost * (i - k - 1))) for k in xrange(i)])
            kj, choices[3] = max_element([(H[i][k] - (gap_open_cost + gap_extend_cost * (j - k - 1))) for k in xrange(j)])
            c, H[i][j] = max_element(choices)
            if c == 1:
                P[i][j] = (i-1, j-1)
            elif c == 2:
                P[i][j] = (ki, j)
            elif c == 3:
                P[i][j] = (i, kj)
            if H[i][j] >= maxH:
                maxs = (i, j)
                maxH = H[i][j]

    i, j = maxs
    indels = {}
    cur_idj = -1
    while i > 0  and  j > 0  and  H[i][j] > 0.0:
        lasti, lastj = i, j
        i, j = P[i][j]
        deli, delj = lasti - i, lastj - j
        if deli and  delj:
            cur_idj = -1
        elif deli:
            if cur_idj >= 0  and  indels[cur_idj].insert_type:
                indels[cur_idj].length += deli
            else:
                cur_idj = lastj                         # count insertion at first index in target to be moved
                indels[cur_idj] = Indel(True, deli)
        elif delj:
            if cur_idj >= 0  and  not indels[cur_idj].insert_type:
                indels[cur_idj].length += delj
            else:
                cur_idj = lastj - 1                     # count deletion at index of last item deleted in target
                indels[cur_idj] = Indel(False, delj)
        else:
            raise Exception("alignment failed")
    return indels


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

