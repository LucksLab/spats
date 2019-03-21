
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

def min_element(l, base = 0):
    ''' @return pair containing the index and value respectively of the minimum element of list l '''
    return min(enumerate(l, base), key=lambda x:x[1])

rev_comp_complementor = string.maketrans("ATCGatcg", "TAGCtagc")

def reverse_complement(seq):
    return str(seq).translate(rev_comp_complementor)[::-1]

def string_find_with_overlap(needle, haystack):
    hlen = len(haystack)
    nlen = len(needle)
    if hlen >= nlen:
        h = haystack.find(needle)
        if -1 != h:
            return h
    h = min(hlen - nlen + 1, 0)
    n = hlen - h
    while h < hlen and n >= 0:
        if needle[:n] == haystack[h:]:
            return h
        h += 1
        n -= 1
    return -1

# hamming distance with tracking and shortcut-out
def string_match_errors(substr, target_str, max_errors = None):
    errors = []
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
        index += 1

    index = 0
    while len(candidates) < max_indices:
        index = target_str.find(second_half, index)
        if -1 == index:
            break
        if index - len(first_half) >= 0:
            errors = string_match_errors(first_half, target_str[index - len(first_half):], me + 1)
            if len(errors) <= me:
                candidates.append(index - len(first_half))
        index += 1

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
                # TAI:  can optimize away min_element here...
                B[i][j], cur_row[j] = min_element([ insert_delete_cost + prev_row[j], insert_delete_cost + cur_row[j-1], substitution_cost + prev_row[j-1] ], 1)
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



def objdict_to_dict(objdict):
    """ convenience for debugging/comparing dictionaries of objects """
    return dict(zip(map(str, objdict.keys()), map(vars, objdict.values())))

def objdict_as_str(objdict, indent = 3):
    """ convenience for debugging/comparing dictionaries of objects """
    import json
    return json.dumps(objdict_to_dict(objdict), indent=indent, sort_keys=True)


class Indel:
    ''' deltas from target to source string '''
    def __init__(self, insert_type, seq, src_index):
        self.insert_type = insert_type    # True for inserts, False for deletes
        self.seq = seq
        self.src_index = src_index        # for easy reference into source string (index inserts to *left* side)


class AlignmentParams:
    '''
     Parameters used for align_strings() algorithm.
     @param  simfn                Similarity function for string elements
     @param  gap_open_cost        Penalty for opening a gap (should be non-negative)
     @param  gap_extend_cost      Penalty for extending an already-open gap (should be non-negative)
     @param  front_biased         Set to True to bias towards aligning fronts
     @param  penalize_ends        If True, will penalize regions before/after alignment for mismatches
     @param  penalize_front_clip  If True and penalize_ends, will penalize front end if a clip is required
     @param  penalize_back_clip   If True and penalize_ends, will penalize back end if a clip is required
    '''
    def __init__(self, simfn = lambda n1, n2: AlignmentParams.char_sim(n1, n2),
                       gap_open_cost = 6, gap_extend_cost = 1,
                       front_biased = True, penalize_ends = True,
                       penalize_front_clip = True, penalize_back_clip = False):
        self.simfn = simfn
        self.gap_open_cost = gap_open_cost
        self.gap_extend_cost = gap_extend_cost
        self.front_biased = front_biased
        self.penalize_ends = penalize_ends
        self.penalize_front_clip = penalize_front_clip
        self.penalize_back_clip = penalize_back_clip

    @staticmethod
    def char_sim(sc, tc, match_value = 2, mismatch_cost = 2):
        """
         Basic character similarity.
         @param  sc             source character
         @param  tc             target character
         @param  match_value    Score bump/reward if they match
         @param  mismatch_cost  Penalty if they don't match (should be non-negative)
         @return similarity between sc and tc
        """
        return match_value if sc == tc else -mismatch_cost

    def __str__(self):
        """ convenience for debugging """
        return str(vars(self))


class Alignment:
    def __init__(self, params, score, target_len, target_start, target_end, src_len, src_start, src_end, indels, mismatched, max_run):
        self.params = params                    # AlignmentParams object used
        self.score = score                      # Smith-Waterman alignment score
        self.orig_target_len = target_len       # the full length of the original target string
        self.target_match_start = target_start  # the first indice in target to which source maps ("site")
        self.target_match_end = target_end      # the last indice in target to which source maps ("end")
        self.target_match_len = target_end - target_start + 1
        self.orig_source_len = src_len          # the full length of the original source string 
        self.src_match_start = src_start        # the first indice in source to match the target
        self.src_match_end = src_end            # the last indice in source to match the target
        self.indels = indels                    # dict mapping indices in the target string to Indel objects
        self.mismatched = mismatched            # list of indices in target that don't match source 
        self.max_run = max_run                  # longest sequence where target matched source in aligned region

    @property
    def indels_delta(self):
        return (self.src_match_end - self.src_match_start) - (self.target_match_end - self.target_match_start)

    def __str__(self):
        """ convenience for debugging """
        return "score={}\ntarget_match_start={}, target_match_end={} (len={})\nsrc_match_start={}, src_match_end={}\nmax_run={}\nmismatched={}\nindels_delta={}\nindels={}\nparams={}".format(self.score, self.target_match_start, self.target_match_end, self.target_match_len, self.src_match_start, self.src_match_end, self.max_run, self.mismatched, self.indels_delta, objdict_as_str(self.indels), self.params) 

    def indels_as_dict(self):
        """ convenience for testing """
        return objdict_to_dict(self.indels)

    def flip(self):
        """
         Flip alignment results for strings that were aligned in reverse.
         @param  slen   Length of the original source string
         @param  tlen   Length of the original target string
        """
        tlen = self.orig_target_len
        slen = self.orig_source_len
        self.target_match_start, self.target_match_end = tlen - self.target_match_end - 1, tlen - self.target_match_start - 1
        self.src_match_start, self.src_match_end = slen - self.src_match_end - 1, slen - self.src_match_start - 1
        self.mismatched = [ (tlen - m - 1) for m in self.mismatched ]
        newindels = {}
        for i, indel in self.indels.iteritems():
            indel.seq = indel.seq[::-1]
            if indel.insert_type:
                indel.src_index = slen - indel.src_index - len(indel.seq)
                newindels[tlen - i] = indel
            else:
                indel.src_index = slen - indel.src_index
                newindels[len(indel.seq) + tlen - i - 2] = indel
        self.indels = newindels


def align_strings(source, target, params = AlignmentParams()):
    """
     Find the indels (insertions and deletions) in a source string relative to a target string.
     First, heuristically tries to find the best alignment of two strings of length M and N, 
     allowing for insertions and deletions using the affine-gap Smith-Waterman (SW) algorithm.
     Once an alignment is found, returns the indel mutations with respect to the target.
     Deletions are indexed at the last item deleted in target.
     Insertions are indexed at the first spot in the target to be moved.
     NOTE: all indices returned are 0-based.
     Warning:  This function has not been optimized; use with care on long strings.
     TAI:  Consider option to produce CIGAR format.

     @param  source  Source string
     @param  target  Target string
     @param  params  AlignmentParams object controlling alignment
     @return an Alignment object
    """
    m = len(source) + 1
    n = len(target) + 1

    H = [[0.0]*n for r in xrange(m)]
    P = [[(0, 0)]*n for r in xrange(m)]

    simfn = params.simfn
    gap_open_cost = params.gap_open_cost
    gap_extend_cost = params.gap_extend_cost
    front_biased = params.front_biased

    ## Dynamic-programming - build up costs from left
    maxH = 0.0
    maxs = [0, 0]
    colmax = [0.0]*n
    colmaxi = [0]*n
    for i in xrange(1, m):
        imo = i -1
        Hi = H[i]
        rowmax, rowmaxj = 0.0, 0
        for j in xrange(1, n):
            jmo = j - 1
            h = H[imo][jmo] + simfn(source[imo], target[jmo])
            h2 = colmax[j] - gap_open_cost - gap_extend_cost * (imo - colmaxi[j])
            h3 = rowmax - gap_open_cost - gap_extend_cost * (jmo - rowmaxj)
            P[i][j] = (imo, jmo)
            if h2 >= h:
                h = h2
                P[i][j] = (colmaxi[j], j)
            if h3 >= h:
                h = h3
                P[i][j] = (i, rowmaxj)
            Hi[j] = max(h, 0.0)    # omit the max0 (and backtrack from corner) for Needleman-Wunsch instead
            if h >= maxH  and  (h > maxH  or  not front_biased  or  abs(i - j) < abs(maxs[0] - maxs[1])):
                maxH = h
                maxs = [i, j]
            if Hi[j] > h2 + gap_open_cost - gap_extend_cost:
                # above same as:  colmax[j] - gap_extend_cost * (i - colmaxi[j])
                colmax[j], colmaxi[j] = Hi[j], i
            if Hi[j] > h3 + gap_open_cost - gap_extend_cost:
                # above same as:  rowmax - gap_extend_cost * (j - rowmaxj)
                rowmax, rowmaxj = Hi[j], j

    ## Now backtrack from max Hij...
    i, j = maxs
    indels = {}
    mismatches = []
    cur_indel = None
    cur_run = 0
    max_run = 0
    score = 0.0
    penalize_ends = params.penalize_ends
    while i > 0  and  j > 0  and  (H[i][j] > 0.0  or  penalize_ends):
        lasti, lastj = i, j
        i, j = P[i][j]
        deli, delj = lasti - i, lastj - j
        if deli and delj:
            #assert(deli == 1  and delj == 1)
            cur_indel = None
            delscore = H[lasti][lastj] - H[i][j]
            if delscore <= 0.0:
                mismatches.append(j)
                if cur_run > max_run:
                    max_run = cur_run
                cur_run = 0
                score += (simfn(source[i], target[j]) if delscore == 0.0 else delscore)
            else:
                cur_run += 1
                score += delscore
        elif deli:
            if cur_indel  and  cur_indel.insert_type:
                cur_indel.seq = source[i:lasti] + cur_indel.seq
                cur_indel.src_index = i
                score -= gap_extend_cost * deli
            else:
                # count insertion at first index in target to be moved
                cur_indel = indels[lastj] = Indel(True, source[i:lasti], i)
                score -= (gap_open_cost + gap_extend_cost * (deli - 1))
            if cur_run > max_run:
                max_run = cur_run
            cur_run = 0
        elif delj:
            if cur_indel  and  not cur_indel.insert_type:
                cur_indel.seq = target[j:lastj] + cur_indel.seq
                score -= gap_extend_cost * delj
            else:
                # count deletion at index of last item deleted in target
                cur_indel = indels[lastj - 1] = Indel(False, target[j:lastj], i)
                score -= (gap_open_cost + gap_extend_cost * (delj - 1))
            if cur_run > max_run:
                max_run = cur_run
            cur_run = 0
        else:
            raise Exception("alignment failed")

    if cur_run > max_run:
        max_run = cur_run

    i = min(i, m - 2)
    j = min(j, n - 2)

    if penalize_ends:
        if i > 0 and j > 0:
            prei, prej = i, j
            prefix_match_score = 0
            prefix_mismatches = []
            while i > 0 and j > 0:
                i -= 1
                j -= 1
                s = simfn(source[i], target[j])
                prefix_match_score += s
                if s <= 0.0:
                    prefix_mismatches.append(j)
            if params.penalize_front_clip and (i + j) > 0:
                prefix_match_score -= (gap_open_cost + (i + j - 1) * gap_extend_cost)
            in_del_cost = -2 * gap_open_cost - gap_extend_cost * max(prei + prej - 2, 0)
            #in_del_cost = -gap_open_cost - gap_extend_cost * max(prei - 1, prej - 1, 0)
            if prefix_match_score > in_del_cost:
                score += prefix_match_score
                mismatches += prefix_mismatches
            else:
                score += in_del_cost
                indels[0] = Indel(True, source[:prei], 0)
                indels[prej - 1] = Indel(False, target[:prej], prei)
                i = j = 0
        elif params.penalize_front_clip:
            if i > 0:
                indels[0] = Indel(True, source[:i], 0)
                score -= (gap_open_cost + (i - 1) * gap_extend_cost)
                i = 0
            elif j > 0:
                indels[j - 1] = Indel(False, target[:j], 0)
                score -= (gap_open_cost + (j - 1) * gap_extend_cost)
                j = 0

        if maxs[0] < m - 1 and maxs[1] < n - 1:
            suffi, suffj = maxs
            suffix_match_score = 0
            suffix_mismatches = []
            while maxs[0] < m - 1 and maxs[1] < n - 1:
                s = simfn(source[maxs[0]], target[maxs[1]])
                suffix_match_score += s
                if s <= 0.0:
                    suffix_mismatches.append(maxs[1])
                maxs[0] += 1
                maxs[1] += 1
            if params.penalize_back_clip and (m + n - 2 - maxs[0] - maxs[1]) > 0:
                suffix_match_score -= (gap_open_cost + (m + n - 3 - maxs[0] - maxs[1]) * gap_extend_cost)
            il = m - 1 - suffi
            jl = n - 1 - suffj
            in_del_cost = -2 * gap_open_cost - gap_extend_cost * max(il + jl - 2, 0)
            #in_del_cost = -gap_open_cost - gap_extend_cost * max(il - 1, jl - 1, 0)
            if suffix_match_score > in_del_cost:
                score += suffix_match_score
                mismatches += suffix_mismatches
            else:
                score += in_del_cost
                indels[suffj] = Indel(True, source[suffi:], suffi)
                indels[n - 2] = Indel(False, target[suffj:], m - 1)
                maxs[0] = m - 1
                maxs[1] = n - 1
        elif params.penalize_back_clip:
            if maxs[0] < m - 1:
                indels[maxs[1]] = Indel(True, source[maxs[0]:], maxs[0])
                score -= (gap_open_cost + (m - maxs[0] - 2) * gap_extend_cost)
                maxs[0] = m - 1
            elif maxs[1] < n - 1:
                indels[n - 2] = Indel(False, target[maxs[1]:], m - 1)
                score -= (gap_open_cost + (n - maxs[1] - 2) * gap_extend_cost)
                maxs[1] = n - 1

    return Alignment(params, score, n - 1, j, maxs[1] - 1, m - 1, i, maxs[0] - 1, indels, mismatches, max_run)


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

