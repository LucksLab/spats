
import json
import os
import shutil

from .logging import LoggingClass
from . import model as M, fs

_speedup = 1

nuc_A = (1)
nuc_C = (1<<1)
nuc_G = (1<<2)
nuc_T = (1<<3)

char_to_mask = { 'A' : nuc_A,
                 'C' : nuc_C,
                 'G' : nuc_G,
                 'T' : nuc_T,
                 'U' : nuc_T,
                 'Y' : (nuc_C | nuc_T),
                 'R' : (nuc_G | nuc_A),
                 'S' : (nuc_C | nuc_G),
                 'W' : (nuc_A | nuc_T),
                 'K' : (nuc_G | nuc_T),
                 'M' : (nuc_A | nuc_C),
                 'B' : (nuc_C | nuc_G | nuc_T),
                 'D' : (nuc_A | nuc_G | nuc_T),
                 'H' : (nuc_A | nuc_C | nuc_T),
                 'V' : (nuc_A | nuc_C | nuc_G),
                 'N' : (nuc_A | nuc_C | nuc_G | nuc_T),
}

rev_comp_complementor = str.maketrans("ATCGatcg", "TAGCtagc")

# xref https://en.wikipedia.org/wiki/FASTQ_format#Encoding
def phredToQuality(ch):
    assert(ord(ch) >= 33 and ord(ch) <= 33 + 42)
    return ord(ch) - 33

def qualityToPhred(q):
    assert(q >= 0 and q <= 42)
    return chr(q + 33)

def reverse_complement(seq):
    return str(seq).translate(rev_comp_complementor)[::-1]
rc = reverse_complement

def crossProduct(lefts, rights):
    if not rights:
        return lefts
    cross = []
    for right in rights:
        cross += [ l + [ right ] for l in lefts ]
    assert(len(cross) == len(lefts) * len(rights))
    return cross

# returns (left, right), where 'left' is the max number of chars extending to the left,
# and 'right' is the max number of chars extending to the right, s.t. s1 matches s2
# when the passed-in ranges (pos, len) are extended to the left and right the
# indicated amounts.
# code is relatively wordy for maximal optimization as this is on the hot path.
#@profile
def longest_match(s1, range1, s2, range2):

    left1 = range1[0]
    left2 = range2[0]
    lmax = min(left1, left2)
    right1 = left1 + range1[1]
    right2 = left2 + range2[1]
    s1len = len(s1)
    s2len = len(s2)
    rmax = min(s1len - right1, s2len - right2)
    if s1[left1:right1] != s2[left2:right2]:
        raise Exception("longest_match must already start with a match")

    left = 1
    right = 0
    c1 = 0
    c2 = 0
    m1 = 0
    m2 = 0
    while True:
        if left > lmax:
            break
        c1 = s1[left1 - left]
        c2 = s2[left2 - left]
        if c1 != c2:
            m1 = char_to_mask[c1]
            m2 = char_to_mask[c2]
            if 0 == (m1 & m2):
                break
        left += 1

    while True:
        if right >= rmax:
            break
        c1 = s1[right1 + right]
        c2 = s2[right2 + right]
        if c1 != c2:
            m1 = char_to_mask[c1]
            m2 = char_to_mask[c2]
            if 0 == (m1 & m2):
                break
        right += 1

    return left - 1, right

def longestExactMatch(s1, range1, s2, range2):

    left1 = range1[0]
    left2 = range2[0]
    lmax = min(left1, left2)
    right1 = left1 + range1[1]
    right2 = left2 + range2[1]
    s1len = len(s1)
    s2len = len(s2)
    rmax = min(s1len - right1, s2len - right2)
    if s1[left1:right1] != s2[left2:right2]:
        raise Exception("longest_match must already start with a match")

    left = 1
    right = 0
    c1 = 0
    c2 = 0
    m1 = 0
    m2 = 0
    while True:
        if left > lmax:
            break
        c1 = s1[left1 - left]
        c2 = s2[left2 - left]
        if c1 != c2:
            break
        left += 1

    while True:
        if right >= rmax:
            break
        c1 = s1[right1 + right]
        c2 = s2[right2 + right]
        if c1 != c2:
            break
        right += 1

    return left - 1, right

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
    prev_row = [j for j in range(n)]
    cur_row = [0] * n
    for i in range(1, m):
        cur_row[0] = i
        for j in range(1, n):
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

    prev_row = [j for j in range(n)]
    cur_row = [0] * n
    B = [[0 for j in range(n)] for i in range(m)]
    for j in range(n):
        B[0][j] = 2

    for i in range(1, m):
        cur_row[0] = i
        B[i][0] = 1
        for j in range(1, n):
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



class Indel(M.Indel):
    ''' deltas from target to source string '''

    def __init__(self, insert_type, seq, src_index, targetIndex, errorIndex = None):
        M.Indel.__init__(self,
                         insertType = insert_type,
                         seq = seq,
                         sourceIndex = src_index,
                         targetIndex = targetIndex)
        self.errorIndex = errorIndex or (targetIndex + (0 if insert_type else (1 - len(seq))))
        self.substSize = 0

    def fromModelObject(mobj):
        res = Indel(mobj.insertType, mobj.seq, mobj.sourceIndex, mobj.targetIndex, mobj.errorIndex)
        res.substSize = mobj.substSize
        return res

    def __repr__(self):
        return "{}:{}@{}".format(f"S-{self.substSize}" if self.substSize else ("I" if self.insert_type else "D"), self.seq, self.errorIndex)

    def description(self):
        return "{}[{},{}]".format(self, self.errorIndex, self.substSize)

    @property
    def mutType(self):
        return fs.substitution if self.substSize else (fs.insertion if self.insertType else fs.deletion)

    @property
    def insert_type(self):
        return self.insertType

    @property
    def src_index(self):
        return self.sourceIndex

    @property
    def size(self):
        return len(self.seq)

    @property
    def delta(self):
        return (self.size - self.substSize) if self.insertType else (0 - self.size)

    def apply(self, source):
        assert(not self.substSize)
        seq, srci = self.seq, self.src_index
        if self.insert_type:
            return source[:srci] + seq + source[srci:]
        else:
            endStart = srci + len(seq)
            assert(seq == source[srci:endStart])
            return source[:srci] + source[endStart:]

    def shift(self, srcIndexDelta, shiftAll = False):
        i = Indel(self.insert_type, self.seq, self.src_index + srcIndexDelta, self.targetIndex, self.errorIndex)
        if shiftAll:
            i.targetIndex += srcIndexDelta
            i.errorIndex += srcIndexDelta
        i.substSize = self.substSize
        return i

    def shiftAll(self, delta):
        i = Indel(self.insert_type, self.seq, self.src_index + delta, self.targetIndex + delta, self.errorIndex + delta)
        i.substSize = self.substSize
        return i

    def shiftToGap(self, gap):
        i = Indel(self.insert_type, self.seq, self.src_index + gap.queryStart, self.targetIndex + gap.segmentStart, self.errorIndex + gap.segmentStart)
        i.substSize = self.substSize
        return i

    def equivalent(self, other):
        return (self.insertType == other.insertType) and (self.seq == other.seq) and (self.sourceIndex == other.sourceIndex) and \
            (self.errorIndex == other.errorIndex) and (self.substSize == other.substSize)


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
        self.src_match_len = src_end - src_start + 1
        self.indels = indels                    # dict mapping indices in the target string to Indel objects
        self.mismatched = mismatched            # list of indices in target that don't match source 
        self.max_run = max_run                  # longest sequence where target matched source in aligned region

    @property
    def indels_delta(self):
        return (self.src_match_end - self.src_match_start) - (self.target_match_end - self.target_match_start)

    def __str__(self):
        """ convenience for debugging """
        return "score={}\ntarget_match_start={}, target_match_end={} (len={})\nsrc_match_start={}, src_match_end={}\nmax_run={}\nmismatched={}\nindels_delta={}\nindels={}\nparams={}".format(self.score, self.target_match_start, self.target_match_end, self.target_match_len, self.src_match_start, self.src_match_end, self.max_run, self.mismatched, self.indels_delta, self.indels, self.params) 

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
        for i, indel in self.indels.items():
            indel.seq = indel.seq[::-1]
            if indel.insert_type:
                indel.src_index = slen - indel.src_index - len(indel.seq)
                newindels[tlen - i] = indel
            else:
                indel.src_index = slen - indel.src_index
                newindels[len(indel.seq) + tlen - i - 2] = indel
        self.indels = newindels


def align_strings(source, target, params = None, hack = False):
    if hack:
        pre = "LMNOPQRSUVWXYZ"
        suf = pre[::-1]
        a = align_strings(pre + source + suf, pre + target + suf, params = params)
        plen = len(pre)
        a.mismatched = [ e - plen for e in a.mismatched ]
        a.indels = { k - plen : v.shiftAll(-plen) for k, v in a.indels.items() }
        a.src_match_end -= 2 * plen
        a.target_match_end -= 2 * plen
        a.score -= 2 * (2 * len(pre)) # match reward of 2 for each char in prefix and suffix
        return a
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
    params = params or AlignmentParams()
    m = len(source) + 1
    n = len(target) + 1

    H = [[0.0]*n for r in range(m)]
    P = [[(0, 0)]*n for r in range(m)]

    simfn = params.simfn
    gap_open_cost = params.gap_open_cost
    gap_extend_cost = params.gap_extend_cost
    front_biased = params.front_biased

    ## Dynamic-programming - build up costs from left
    maxH = 0.0
    maxs = [0, 0]
    colmax = [0.0]*n
    colmaxi = [0]*n
    for i in range(1, m):
        imo = i -1
        Hi = H[i]
        rowmax, rowmaxj = 0.0, 0
        for j in range(1, n):
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
                cur_indel = indels[lastj] = Indel(True, source[i:lasti], i, lastj)
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
                cur_indel = indels[lastj - 1] = Indel(False, target[j:lastj], i, lastj - 1)
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
        # TAI:  should this be "or" instead of "and"?  or should we compare an indel to muts in the fall-through else?
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
                indels[0] = Indel(True, source[:prei], 0, 0)
                indels[prej - 1] = Indel(False, target[:prej], prei, prej - 1)
                i = j = 0
        elif params.penalize_front_clip:
            if i > 0:
                indels[0] = Indel(True, source[:i], 0, 0)
                score -= (gap_open_cost + (i - 1) * gap_extend_cost)
                i = 0
            elif j > 0:
                indels[j - 1] = Indel(False, target[:j], 0, j - 1)
                score -= (gap_open_cost + (j - 1) * gap_extend_cost)
                j = 0

        # TAI:  should this be "or" instead of "and"?  or should we compare an indel to muts in the fall-through else?
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
                indels[suffj] = Indel(True, source[suffi:], suffi, suffj)
                indels[n - 2] = Indel(False, target[suffj:], m - 1, n - 2)
                maxs[0] = m - 1
                maxs[1] = n - 1
        elif params.penalize_back_clip:
            if maxs[0] < m - 1:
                indels[maxs[1]] = Indel(True, source[maxs[0]:], maxs[0], maxs[1])
                score -= (gap_open_cost + (m - maxs[0] - 2) * gap_extend_cost)
                maxs[0] = m - 1
            elif maxs[1] < n - 1:
                indels[n - 2] = Indel(False, target[maxs[1]:], m - 1, n - 2)
                score -= (gap_open_cost + (n - maxs[1] - 2) * gap_extend_cost)
                maxs[1] = n - 1

    return Alignment(params, score, n - 1, j, maxs[1] - 1, m - 1, i, maxs[0] - 1, indels, mismatches, max_run)


def testIndelApply():
    i = Indel(True, "abcd", 4, 0)
    assert(i.apply("GHIJKLMN") == "GHIJabcdKLMN")
    d = Indel(False, "def", 3, 0)
    assert(d.apply("abcdefghi") == "abcghi")


def eaTest():

    #case = [ R1/R2,
    #         target seq,
    #         matchStart, matchLen, matchIndex,
    #         expectedMatchStart, expectedMatchLen, expectedMatchIndex,
    #         expectedMutSpots, expectedIndelSpots, gap_open_cost ]
    cases = [ [ "TGAACAGCGACTAGGCTTCCCAAAGTACCAGTTTGCCAC",
    #0           |               |
                "TGAACAGCGACTAGGCTCTTCCGATCT",
                0, 17, 0,
                0, 17, 0,
                [], [] ],

              [ "GTTCGGAGAAGCATGACGGACAAGTACAAGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCC",
    #1                                       !|                                            |
                                            "GGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCCTTCTGGGCAGTTGATGCC",
                29, 46, 1,
                29, 46, 1,
                [], [] ],

              [ "GTTCGGAGAAGCATGACGGACAAGTACAAGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCC",
    #2           |                      |!                                                 |
                "GTTCGGAGAAGCATGACGGACAAGAACAAGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCC",
    #                                        !|                                            |
                0, 24, 0,
                0, 75, 0,
                [24], [] ],

              [ "GTTCGGAGAAGCATGACGGACAAGTACAAGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCC",
    #3           |                       !|                                                |
                "GTTCGGAGAAGCATGACGGACAAGAACAAGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCC",
    #                                        !|                                            |
                25, 50, 25,
                0, 75, 0,
                [24], [] ],

              [ "TTTTTGATAAGCATGACGGACAAGTACAGGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCC",
    #4                | !|                  |
                     "GAGAAGCATGACGGACAAGTACA",
                8, 20, 3,
                5, 23, 0,
                [2], [] ],

              [ "GGACGGGAAGCATGACGGACAAGTACAGGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCCT",
    #5               |!|                   |
                    "GAGAAGCATGACGGACAAGTACA",
                6, 21, 2,
                6, 21, 2,
                [], [] ],

              [ "GGACGGGAAGCATGACGGACAAGTACAGGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCCT",
    #6              | !|                   |
                   "CGAGAAGCATGACGGACAAGTACA",
                6, 21, 3,
                3, 24, 0,
                [2], [] ],

              [ "ACTATGAGAGCATGACGGACAAGTACAGGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCCT",
    #7                  -|                 |
                    "GAGAAGCATGACGGACAAGTACA",
                8, 19, 4,
                8, 19, 4,
                [], [], 8 ], # set the gap open cost high so this indel doesn't take

              [ "ACTATGAGAGCATGACGGACAAGTACAGGCTGCACCTGTCAGTGGCCGACCTCCTCTTTGTCATCACGCTTCCCT",
    #8                  -|                 |
                    "GAGAAGCATGACGGACAAGTACA",
                8, 19, 4,
                5, 23, 0,
                [], [4], 5 ], # set gap open cost to 5 so that this indel takes

              ["CAGGGGAAAAAAGGAGAAGGGGGAAAAGGAAGAAAGGGAAAAGAAGGGAAGAGGGAGGACAAAGAGGAGGGAGGC",
    #9                                       ! !!  !! !!!|        |! !  !            |!!  |
                                          "GGCATCAACTGCCCAGAAGGGAAGCGTGATGACAAAGAGGAGGTCGGCCACTGACAGGTGCAGC",
    #                                                  14^        ^23                ^42  ^47
               41, 10, 14,
               41, 34, 14,
               [10, 12, 15, 29, 30], [], 5 ],


              # same seqs as last case, but somehow this one extends too far left (8 muts and 6 matches to left of site 14)
              # (it's because match_start=57 now, so 14-30 will bump up the alignment score.).
              ["CAGGGGAAAAAAGGAGAAGGGGGAAAAGGAAGAAAGGGAAAAGAAGGGAAGAGGGAGGACAAAGAGGAGGGAGGC",
    #10                                      ! !!  !! !!!|         ! !  !|           |!!  |
                                          "GGCATCAACTGCCCAGAAGGGAAGCGTGATGACAAAGAGGAGGTCGGCCACTGACAGGTGCAGC",
    #                                                  14^               ^30         ^42  ^47
               57, 13, 30,
               41, 34, 14,
               [10, 12, 15, 29, 30], [], 5 ],

              # want this to extend past the deletion b/c plenty of match, but then stop
              # before all of the muts
              # TAI:  consider with a very high mut (mismatch) cost
              [ "TGAACAGCGACT GGCTCTTCCCAAAGTACCAGTTTGCCACGGCATCAACTGCCCAGAAGGGAAGCGTGATGACAA",
    #11          |          |-        |! !!!
                "TGAACAGCGACTAGGCTCTTCCGATCT",
                0, 12, 0,
                0, 22, 0,
                [], [12], 5 ],

              [ "TATGGGAGAAGCATGACGGACAAGTACAGGCTGCACCTGTCAGGGGACGACAGAAGCGGGGTAAGAACGATGACC",
    #12                                      |             |!  !    !!!!! !!                                                                 
                                            "GGCTGCACCTGTCAGTGGCCGACCTCCTCTT",
                28, 15, 0,
                28, 23, 0,
                [15, 18], [], 5],

    ]

    ap = AlignmentParams(penalize_ends = False, penalize_front_clip = False, penalize_back_clip = False)
    assertCases = True

    for idx in range(len(cases)):
        case = cases[idx]
        rN, target = [ x.replace(' ', '') for x in case[:2] ]
        matchStart, matchLen, matchIndex = case[2:5]
        rEnd, targetEnd = matchStart + matchLen, matchIndex + matchLen
        assert rN[matchStart:rEnd] == target[matchIndex:targetEnd], "{} != {}".format(rN[matchStart:rEnd], target[matchIndex:targetEnd])
        # should always start with a maximal match
        if matchStart > 0:
            assert(rN[matchStart - 1] != target[matchIndex - 1])
        if rEnd < len(rN) and targetEnd < len(target):
            assert(rN[rEnd] != target[targetEnd])
        if len(case) > 10:
            ap.gap_open_cost = case[10]
        else:
            ap.gap_open_cost = 6
        a = alignment = align_strings(rN, target, ap)
        matchErrors = sorted([ x - a.target_match_start for x in a.mismatched ])
        indels = { k - a.target_match_start : v for k, v in a.indels.items() }
        print("Case {}:".format(idx), a.src_match_start, a.target_match_len, a.target_match_start, matchErrors, indels)
        #print(a)
        if assertCases:
            assert a.src_match_start == case[5], "{} != {}".format(a.src_match_start, case[5])
            assert a.target_match_len == case[6], "{} != {}".format(a.target_match_len, case[6])
            assert a.target_match_start == case[7], "{} != {}".format(a.target_match_index, case[7])
            assert matchErrors == case[8], "{} != {}".format(matchErrors, case[8])
            assert sorted(list(indels.keys())) == case[9], "{} != {}".format(sorted(list(indels.keys())), case[9])


def dataWithFileAtPath(path):
    with open(path, 'rb') as f:
        return f.read()

def strWithFileAtPath(path):
    with open(path, 'r') as f:
        return f.read()

def writeDataToFileAtPath(inData, path, skipIfNoDiff = False, overwriteReadonly = False):
    if isinstance(inData, str):
        data = inData.encode()
    else:
        data = inData
    if os.path.exists(path):
        if skipIfNoDiff:
            try:
                existing = dataWithFileAtPath(path)
                if existing == data:
                    return False
            except:
                pass
        if overwriteReadonly:
            os.chmod(path, 0o644) # for windoze
            os.remove(path)
    with open(path, 'wb') as f:
        f.write(data)
    return True

def displaySize(num, cost = False):
    if num < .0001:
        return "{:.08f}".format(num)
    elif num < .01:
        return "{:.04f}".format(num)
    elif num < 1:
        return "{:.02f}".format(num)
    elif num < 1000:
        return "{:.1f}".format(num)
    elif num < 1000 * 1000:
        return "{}k".format(int(num / 1000))
    elif num < 1000 * 1000 * 1000:
        return "{}m".format(int(num / (1000 * 1000)))
    elif num < 1000 * 1000 * 1000 * 1000:
        return "{}{}".format(int(num / (1000 * 1000 * 1000)), "b" if cost else "g")
    else:
        return "{}t".format(int(num / (1000 * 1000 * 1000 * 1000)))

def flatten(listOfLists):
    return [item for sublist in listOfLists for item in sublist]

def jsonAtPath(path, emptyIfMissing = False):
    if emptyIfMissing  and  not os.path.exists(path):
        return {}
    return json.loads(strWithFileAtPath(path))

def writeJsonToPath(dictionary, path, pretty = False, skipIfNoDiff = False, overwriteReadonly = False):
    unserializable = None
    try:
        data = json.dumps(dictionary, sort_keys = pretty, indent = 4 if pretty else None, separators = (',', ': ') if pretty else None)
    except Exception as e:
        if str(e).endswith('JSON serializable'):
            unserializable = e
        else:
            raise
    if unserializable:
        raise Exception("unserializable object in dictionary: {}".format(str(unserializable)))
    return writeDataToFileAtPath(data,
                                 path,
                                 skipIfNoDiff = skipIfNoDiff,
                                 overwriteReadonly = overwriteReadonly)
def prettyJson(dict):
    return json.dumps(dict, sort_keys = True, indent = 4, separators = (',', ': '))

def ensureFolderExists(inPath, wipe = False):
    path = os.path.expanduser(inPath)
    exists = os.path.exists(path)
    if wipe and exists:
        wipePath(path)
        exists = False
    if not exists:
        os.makedirs(path)
    return path

def csvForRows(rows, skipQuotes = False):
    formatter = '{}' if skipQuotes else '"{}"'
    return '\n'.join([','.join([formatter.format('' if x == None else x) for x in r]) for r in rows])

def csvAtPath(path, expectSameRowCounts = True, csvStr = None):
    rows = []
    for line in (csvStr or strWithFileAtPath(path)).split('\n'):
        if not line.strip():
            continue
        row = []
        left = line
        while left:
            if left.startswith('"'):
                next = left.find('"', 1)
                if -1 == next:
                    raise Exception("Bad entry: row {}, col {}".format(len(rows), 1 + len(row)))
                row.append(left[1:next])
                left = left[next + 1:]
            else:
                next = left.find(',')
                if -1 == next:
                    row.append(left)
                    left = None
                else:
                    row.append(left[:next])
                    left = left[next:]
            if left and left[0] == ',':
                left = left[1:]
        rows.append(row)
    if rows and expectSameRowCounts:
        checkLen = len(rows[0])
        for i in range(len(rows)):
            if len(rows[i]) != checkLen:
                raise Exception("Unexpected row count: row {} has {} items, row 1 has {}".format(i + 1, len(rows[i]), checkLen))
    return rows

def writeCsvToPath(rows, path, skipQuotes = False):
    writeDataToFileAtPath(csvForRows(rows, skipQuotes = skipQuotes), path)

def movePath(srcPath, dstPath, ignoreExisting = False):
    assert(os.path.exists(srcPath))
    if not ignoreExisting:
        assert(not os.path.exists(dstPath)) # pass ignoreExisting = True
    shutil.move(srcPath, dstPath)
