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

class Mask(object):

    def __init__(self, chars):
        self.chars = chars.upper()
        self.values = [ char_to_mask[ch] for ch in chars ]
        self.total = 0
        self.kept = 0
        self._counts = {}

    def matches(self, seq):
        # raises if len(seq) < len(self.values)
        maskvals = self.values
        for i in range(len(maskvals)):
            seqval = char_to_mask[seq[i]]
            if 0 == (seqval & maskvals[i]):
                return False
        return True

    def register_count(self, target, site, multiplicity = 1):
        self.counts(target)[site] += multiplicity
        self.kept += multiplicity

    def counts(self, target):
        counts = self._counts.get(target.name)
        if not counts:
            counts = [ 0 for x in range(target.n + 1) ] # TODO: numpy.empty(n, dtype=int) ??
            self._counts[target.name] = counts
        return counts

    def count_data(self):
        return self._counts

    def update_with_count_data(self, count_data, target_map):
        for key, values in count_data.iteritems():
            our_counts = self.counts(target_map[key])
            for j in range(len(values)):
                our_counts[j] += values[j]


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
