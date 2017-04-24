import math
import os
import string
import sys
from itertools import izip


class SpatsConfig(object):

    def __init__(self):

        self.adapter_t = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"  # JJB: shows revcomped on the end of R2, from the back
        self.adapter_b = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"                          # JJB: shows as-is on the end of R1, from the front

        self.debug = False
        self.minimum_target_match_length = 8 # want to make this as big as possible for speed.
        self.minimum_adapter_len = 0
        self.allow_indeterminate = False
        self.show_progress = True
        self.show_id_to_site = False

        # TODO: would it be better to base these off of some kind of acceptable %/freq, to control for error:length ratios?
        self.allowed_adapter_errors = 0
        self.allowed_target_errors = 0

        # DELME: options to attempt v102 compliance:
        #self.allowed_adapter_errors = 1
        #self.minimum_adapter_len = 10
        self.minimum_adapter_matches = 8
        self.allow_errors_in_last_4_of_R2 = True
        self.ignore_minimal_adapter = True


spats_config = SpatsConfig()


def _warn(stuff):
    print stuff

def _debug(stuff):
    if spats_config.debug:
        print stuff



rev_comp_complementor = string.maketrans("ATCGatcg", "TAGCtagc")

def reverse_complement(seq):
    return seq.translate(rev_comp_complementor)[::-1]


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

def matches_mask(seq, maskvals):
    # assume len(seq) >= len(maskvals) -- otherwise this will raise
    for i in range(len(maskvals)):
        seqval = char_to_mask[seq[i]]
        if 0 == (seqval & maskvals[i]):
            return False
    return True

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


class Mask(object):

    def __init__(self, chars):
        self.chars = chars
        self.values = [ char_to_mask[ch] for ch in chars ]


class FastqRecord(object):

    def __init__(self):
        self.recordNumber = 0
        self.reset()

    def reset(self):
        self.identifier = None
        self.sequence = None
        self.identifier2 = None
        self.quality = None

    def read(self, infile):
        first = infile.readline()
        if not first or len(first) == 0:
            self.reset()
            return
        self.identifier = first.lstrip('@').rstrip('\r\n').split(' ')[0]
        self.sequence = infile.readline().rstrip('\r\n')
        self.identifier2 = infile.readline().rstrip('\r\n')
        self.quality = infile.readline().rstrip('\r\n')

    def write(self, outfile):
        for line in [ self.identifier, self.sequence, self.identifier2, self.quality ]:
            outfile.write(line)
            outfile.write('\n')

    def reverse_complement(self):
        self.sequence = reverse_complement(self.sequence)
        self.quality = self.quality[::-1]


def fasta_parse(target_path):
    pairs = []
    with open(target_path, 'rb') as infile:
        while True:
            line = infile.readline()
            if 0 == len(line):
                break
            name = line.strip('>\n')
            line = infile.readline()
            seq = line.strip()
            if name and seq:
                pairs.append((name, seq))
    return pairs


class Sequence(object):

    def __init__(self):
        self._reset(None)

    def _reset(self, seq):
        self._seq = seq
        self._length = len(seq) if seq else None
        self.match_start = None
        self.match_len = None
        self.match_index = None
        self._ltrim = 0
        self._rtrim = 0
        self.adapter_errors = []
        self.match_errors = []

    def set_seq(self, seq):
        self._reset(seq)

    @property
    def original_seq(self):
        return self._seq

    @property
    def original_len(self):
        return self._length

    @property
    def subsequence(self):
        if self._rtrim:
            return self._seq[self._ltrim:-self._rtrim]
        else:
            return self._seq[self._ltrim:]

    @property
    def seq_len(self):
        return self._length - self._ltrim - self._rtrim

    @property
    def reverse_complement(self):
        return reverse_complement(self.subsequence)

    @property
    def matched(self):
        return bool(self.match_len)

    @property
    def left(self):
        return self.match_index if self.match_len else None

    @property
    def right(self):
        return self.match_index + self.match_len if self.match_len else None

    @property
    def trimmed(self):
        return (0 != self._rtrim)

    def find_in_target(self, target, reverse_complement = False):
        seq = self.reverse_complement if reverse_complement else self.subsequence
        self.match_start, self.match_len, self.match_index = target.find_partial(seq, spats_config.minimum_target_match_length)

    def trim(self, length, reverse_complement = False):
        self._rtrim = length
        if reverse_complement:
            delta = self.match_len - self.seq_len
            if delta > 0:
                _debug("trim reducing original match_len {} -> {}".format(self.match_len, self.seq_len))
                self.match_len = self.seq_len
                self.match_start += delta
                self.match_index += delta

    def match_to_seq(self, reverse_complement = False):
        if reverse_complement:
            self.match_index -= (self._rtrim - self.match_start)
            self.match_start = self._rtrim
            self.match_len = self.seq_len
        else:
            self.match_index -= self.match_start
            self.match_start = 0
            self.match_len = self.seq_len

class Pair(object):

    def __init__(self):
        self.reset()

    def reset(self):
        self.identifier = None
        self.r1 = Sequence()
        self.r2 = Sequence()
        self.site = None
        self.mask = None
        self.failure = None

    def set_from_data(self, identifier, r1_seq, r2_seq):
        self.reset()
        self.identifier = identifier
        self.r1.set_seq(r1_seq)
        self.r2.set_seq(r2_seq)

    def set_from_records(self, r1_record, r2_record):
        if not r1_record.identifier or r1_record.identifier != r2_record.identifier:
            raise Exception("Invalid record IDs for pair: {}, {}".format(r1_record.identifier, r2_record.identifier))
        self.reset()
        self.identifier = r1_record.identifier
        self.r1.set_seq(r1_record.sequence)
        self.r2.set_seq(r2_record.sequence)

    def is_determinate(self):
        return set(self.r1.original_seq + self.r2.original_seq) <= set('ACGT')

    def set_mask(self, mask):
        self.mask = mask
        mask.total += 1
        self.r1._ltrim = 4

    @property
    def matched(self):
        return (self.r1.match_len and self.r2.match_len)

    @property
    def has_site(self):
        return bool(self.site is not None)

    @property
    def left(self):
        return self.r2.left

    @property
    def right(self):
        return self.r1.right

    def register_count(self):
        self.mask.kept += 1
        self.site = self.left
        self.mask.counts[self.site] += 1



# returns (left, right), where 'left' is the max number of chars extending to the left,
# and 'right' is the max number of chars extending to the right, s.t. s1 matches s2
# when the passed-in ranges (pos, len) are extended to the left and right the
# indicated amounts.
#@profile
def faster_longest_match(s1, range1, s2, range2):
    left1 = range1[0]
    left2 = range2[0]
    lmax = max(left1, left2)
    right1 = left1 + range1[1]
    right2 = left2 + range2[1]
    s1len = len(s1)
    s2len = len(s2)
    rmax = min(s1len - right1, s2len - right2) - 1
    if s1[left1:right1] != s2[left2:right2]:
        raise Exception("longest_match must already start with a match")
    left = 1
    right = 0
    c1 = 0
    c2 = 0
    m1 = 0
    m2 = 0
    #while left1 >= 0 and left2 >= 0 and 0 != (char_to_mask[s1[left1]] & char_to_mask[s2[left2]]):
    #    left1 -= 1
    #    left2 -= 1
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

    #while right1 <= s1len and right2 <= s2len and 0 != (char_to_mask[s1[right1 - 1]] & char_to_mask[s2[right2 - 1]]):
    #    right1 += 1
    #    right2 += 1
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

    #return range1[0] - left1 - 1, right1 - range1[1] - range1[0] - 1
    return left - 1, right

# returns (left, right), where 'left' is the max number of chars extending to the left,
# and 'right' is the max number of chars extending to the right, s.t. s1 matches s2
# when the passed-in ranges (pos, len) are extended to the left and right the
# indicated amounts.
#@profile
def longest_match(s1, range1, s2, range2):
    left1 = range1[0]
    left2 = range2[0]
    right1 = left1 + range1[1]
    right2 = left2 + range2[1]
    s1len = len(s1)
    s2len = len(s2)
    if s1[left1:right1] != s2[left2:right2]:
        raise Exception("longest_match must already start with a match")
    while left1 >= 0 and left2 >= 0 and 0 != (char_to_mask[s1[left1]] & char_to_mask[s2[left2]]):
        left1 -= 1
        left2 -= 1
    while right1 <= s1len and right2 <= s2len and 0 != (char_to_mask[s1[right1 - 1]] & char_to_mask[s2[right2 - 1]]):
        right1 += 1
        right2 += 1
    return range1[0] - left1 - 1, right1 - range1[1] - range1[0] - 1



class Target(object):

    def __init__(self, name, seq, index_word_length = 8):
        self.name = name
        self.seq = seq
        self.n = len(seq)
        self._index = None
        self.index_word_length = index_word_length
        self._warned = False

    def index(self):
        seq = self.seq
        index = {}
        word_len = self.index_word_length
        for i in range(len(seq) - word_len + 1):
            key = seq[i:(i + word_len)]
            sites = index.get(key)
            if not sites:
                sites = []
                index[key] = sites
            sites.append(i)
        self._index = index

    def find_exact(self, query):
        assert(self._index)
        word_len = self.index_word_length
        query_len = len(query)
        if query_len < word_len:
            raise Exception("Query too short: len({}) < {}".format(query, word_len))
        query_key = query[0:word_len]
        for index in self._index.get(query_key, []):
            if self.seq[index:index+query_len] == query:
                return index
        return None

    # returns (query_start_index, match_len, sequence_index), where:
    #  query_start_index: the index into the query where the match starts
    #  match_len: the length of the match
    #  sequence_index: the index into the target sequence where the match starts
    def find_partial(self, query, minimum_length = None):
        assert(self._index)
        min_len = minimum_length or (self.index_word_length << 1)
        word_len = self.index_word_length
        if min_len < word_len:
            raise Exception("minimum_length too short: {} < {}".format(min_len, word_len))
        check_every = min_len - word_len # norah has proved that this guarantees finding a match if it exists
        if check_every < 4 and not self._warned:
            print "Warning: minimum_length {} is not much longer than index length {}".format(min_len, word_len)
            self._warned = True
        query_len = len(query)
        check_sites = range(0, query_len - check_every, max(check_every, 1))
        check_sites.append(query_len - check_every)
        candidate = (None, None, None)
        # NOTE: it's important to check all sites, and all hits -- to find the longest match.
        for site in check_sites:
            site_key = query[site:site+word_len]
            #print "CS: {}, {}".format(site, site_key)
            for index in self._index.get(site_key, []):
                #print "GOT: " + str(index)
                left, right = longest_match(query, (site, word_len), self.seq, (index, word_len))
                total_len = left + right + word_len
                #print "extends: <--{}, -->{} / {} ({})".format(left, right, total_len, min_len)
                if total_len >= min_len:
                    if total_len == len(query):
                        # we can return immediately if we've got a full match...
                        return site - left, total_len, index - left
                    elif not candidate[1] or total_len > candidate[1]:
                        # ...otherwise, keep it if it's the best match so far
                        candidate = (site - left, total_len, index - left)
                        #print "C: {}".format(candidate)
        return candidate


class Spats(object):

    def __init__(self, target_path, output_folder):
        self.target_path = target_path
        self.output_folder = output_folder

        # user-configurable parameters
        self.masks = [ 'RRRY', 'YYYR' ]
        self.quiet = False
        self.adapter_t = spats_config.adapter_t
        self.adapter_b = spats_config.adapter_b

        # private vars
        self._targets = None
        self._masks = None

        self.total_pairs = self.processed_pairs = self.chucked_pairs = 0

    def setup(self):
        self._indexTargets()
        self._setupMasks()
        self.adapter_t_rc = reverse_complement(self.adapter_t)

    def _indexTargets(self):
        if self._targets:
            return
        targets = []
        for name, seq in fasta_parse(self.target_path):
            target = Target(name, seq)
            target.index()
            targets.append(target)
        self._targets = targets

        # TODO: handle multiple targets. for now just use the one
        self._target = targets[0]

    def _setupMasks(self):
        if self._masks:
            return
        masks = map(Mask, self.masks)
        for mask in masks:
            mask.counts = [ 0 for x in range(self._target.n + 1) ] # TODO: numpy.empty(num_sites, dtype=int) is better
            mask.total = 0
            mask.kept = 0
        self._masks = masks
        self._maskSize = max([ len(m.chars) for m in masks ])

    def _match_mask(self, pair):
        seq = pair.r1.original_seq
        if len(seq) <= self._maskSize:
            _warn("sequence {} is too short for masking".format(seq))
            return
        for mask in self._masks:
            if matches_mask(seq, mask.values):
                pair.set_mask(mask)
                return

    def _find_matches(self, pair):
        pair.r1.find_in_target(self._target, reverse_complement = True)
        pair.r2.find_in_target(self._target)
        _debug([pair.r1.match_start, pair.r1.match_len, pair.r1.match_index, "--", pair.r2.match_start, pair.r2.match_len, pair.r2.match_index])

    def _trim_adapters(self, pair):

        # if we're here, we know that R2 hangs over the right edge. first, figure out by how much

        r2_seq = pair.r2.subsequence
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        r2_len_should_be = self._target.n - r2_start_in_target
        r2_length_to_trim = pair.r2.seq_len - r2_len_should_be
        pair.r2.trim(r2_length_to_trim)
        _debug("R2 trim: {}, {}, {}".format(r2_start_in_target, r2_len_should_be, r2_length_to_trim))

        if spats_config.minimum_adapter_len and r2_length_to_trim - 4 < spats_config.minimum_adapter_len:
            _debug("  !! v102 minimum adapter len {}".format(r2_length_to_trim - 4))
            return False

        if r2_length_to_trim < 4:
            # TODO: should we verify that this matches RC of R1 handle, and register errors for bp that don't?
            # for now, just ignore this part
            # also, this means that there's nothing to trim from R1, so we're done
            return True

        if spats_config.ignore_minimal_adapter and r2_length_to_trim <= 8:
            _debug("  !! v102 auto-trimming short adapter")
            pair.r1.trim(r2_length_to_trim - 4, reverse_complement = True)
            return True

        # find out how good of a match the end of R2 is for adapter_t_rc
        r2_adapter_match = r2_seq[4-r2_length_to_trim:]
        pair.r2.adapter_errors = string_match_errors(r2_adapter_match, self.adapter_t_rc)
        _debug("  check = {}, errors = {}".format(r2_adapter_match, pair.r2.adapter_errors))
        if spats_config.minimum_adapter_matches and len(pair.r2.adapter_errors) > spats_config.allowed_adapter_errors:
            if r2_length_to_trim - 4 - len(pair.r2.adapter_errors) >= spats_config.minimum_adapter_matches:
                _debug(" !! v102 accepting slightly mismatched R2 adapter with minimum len")
            else:
                _debug(" !! v102 rejecting short/mismatched R2 adapter")
                return False
        elif len(pair.r2.adapter_errors) > spats_config.allowed_adapter_errors:
            return False

        # now, same thing on r1 (smaller trim b/c of no handle, hence -4)
        r1_seq = pair.r1.subsequence
        r1_length_to_trim = r2_length_to_trim - 4
        r1_adapter_match = r1_seq[-r1_length_to_trim:]
        pair.r1.trim(r1_length_to_trim, reverse_complement = True)
        pair.r1.adapter_errors = string_match_errors(r1_adapter_match, self.adapter_b)
        _debug("  R1 check = {}, errors = {}".format(r1_adapter_match, pair.r1.adapter_errors))
        if spats_config.minimum_adapter_matches and len(pair.r1.adapter_errors) > spats_config.allowed_adapter_errors:
            if r1_length_to_trim - len(pair.r1.adapter_errors) >= spats_config.minimum_adapter_matches:
                _debug(" !! v102 accepting slightly mismatched R1 adapter with minimum len")
            else:
                _debug(" !! v102 rejecting short/mismatched R1 adapter")
                return False
        elif len(pair.r1.adapter_errors) > spats_config.allowed_adapter_errors:
            return False

        _debug("successful adapter trim of {}/{} bp from R1/R2".format(pair.r1._rtrim, pair.r2._rtrim))

        return True

    def process_pair(self, pair):
        _debug("> processing " + pair.identifier + "\n  --> " + pair.r1.original_seq + " , " + pair.r2.original_seq)
        _debug("  rc(R1): {}".format(pair.r1.reverse_complement))
        self._process_pair(pair)
        if pair.failure:
            _debug(pair.failure)
        else:
            assert(pair.has_site)
            _debug("  ===> KEPT {}-{}".format(pair.left, pair.right))
            if spats_config.show_id_to_site:
                print "{} --> {}".format(pair.identifier, pair.site) #, pair.mask.chars)


    def _process_pair(self, pair):

        if not spats_config.allow_indeterminate  and  not pair.is_determinate():
            pair.failure = "indeterminate sequence failure"
            return

        self._match_mask(pair)
        if not pair.mask:
            pair.failure = "mask failure"
            return

        self._find_matches(pair)
        if not pair.matched:
            pair.failure = "no match"
            return

        # this is where R2 should start (if not a complete match, then r2.match_start will be > 0)
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        if r2_start_in_target < 0:
            pair.failure = "R2 to left of site 0 failure on R2 for: {}".format(pair.identifier)
            return
        elif r2_start_in_target + pair.r2.original_len <= self._target.n:
            # we're in the middle -- no problem
            pass
        elif not self._trim_adapters(pair):
            # we're over the right edge, and adapter trimming failed
            pair.failure = pair.failure or "adapter trim failure"
            return
        else:
            # we're at the right and trimming went ok, cool
            pass

        if pair.r2.match_len == pair.r2.seq_len  and  pair.r1.match_len == pair.r1.seq_len:
            # everything that wasn't adapter trimmed was matched -- nothing to do
            pass
        else:
            # set the match to be the rest of the (possibly trimmed) sequence, and count errors
            pair.r1.match_to_seq(reverse_complement = True)
            pair.r2.match_to_seq()
            target = self._target.seq
            pair.r1.match_errors = string_match_errors(pair.r1.reverse_complement, target[pair.r1.match_index:])
            pair.r2.match_errors = string_match_errors(pair.r2.subsequence, target[pair.r2.match_index:])

            if spats_config.allow_errors_in_last_4_of_R2:
                errors = [ e for e in pair.r2.match_errors if e < pair.r2.seq_len - 4]
                _debug([pair.r2.seq_len, pair.r2.match_errors, errors ])
                if len(errors) < len(pair.r2.match_errors):
                    _debug("   !! v102 allowing errors in last 4 of R2")
                pair.r2.match_errors = errors

            if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > spats_config.allowed_target_errors:
                if pair.r1.match_errors:
                    _debug("R1 errors: {}".format(pair.r1.match_errors))
                if pair.r2.match_errors:
                    _debug("R2 errors: {}".format(pair.r2.match_errors))
                pair.failure = "match errors failure"
                return

        n = self._target.n
        assert(pair.matched and pair.left >= 0 and pair.left <= n)

        # NOTE: this might change later due to "starts"
        if pair.right != n:
            pair.failure = "R1 right edge failure: {} - {}, n={}".format(pair.left, pair.right, n)
            return

        pair.register_count()


    def process_pair_data(self, data_r1_path, data_r2_path):
        total_pairs = 0
        chucked_pairs = 0
        processed_pairs = 0

        with open(data_r1_path, 'rb') as r1_in:
            with open(data_r2_path, 'rb') as r2_in:
                r1_record = FastqRecord()
                r2_record = FastqRecord()
                pair = Pair()
                while True:
                    r1_record.read(r1_in)
                    if not r1_record.identifier:
                        break
                    r2_record.read(r2_in)
                    pair.set_from_records(r1_record, r2_record)
                    total_pairs += 1
                    pair.numeric_id = total_pairs

                    self.process_pair(pair)
                    if not pair.mask:
                        chucked_pairs += 1
                    elif pair.has_site:
                        processed_pairs +=1

                    if spats_config.show_progress and 0 == total_pairs % 20000:
                        sys.stderr.write('.')
                        sys.stderr.flush()

        self.total_pairs += total_pairs
        self.chucked_pairs += chucked_pairs
        self.processed_pairs += processed_pairs

        if not self.quiet:
            self.report_counts()

    def report_counts(self):
        m0 = self._masks[0]
        m1 = self._masks[1]
        format_str = "Processed {tot} properly paired fragments, " + \
                     "kept {kept0}/{tot0} ({pct0:.1f}%) treated, " + \
                     "{kept1}/{tot1} ({pct1:1f}%) untreated"
        print format_str.format(tot = m0.total + m1.total,
                                kept0 = m0.kept,
                                tot0 = m0.total,
                                pct0 = (100.0 * float(m0.kept)) / float(m0.total),
                                kept1 = m1.kept,
                                tot1 = m1.total,
                                pct1 = (100.0 * float(m1.kept)) / float(m1.total))

    def compute_profiles(self):
        # TODO: use numpy here ?
        n = self._target.n
        treated_counts = self._masks[0].counts
        untreated_counts = self._masks[1].counts
        betas = [ 0 for x in range(n+1) ]
        thetas = [ 0 for x in range(n+1) ]
        treated_sum = 0.0    # keep a running sum
        untreated_sum = 0.0  # for both channels
        running_c_sum = 0.0  # can also do it for c

        for k in range(n):
            X_k = float(treated_counts[k])
            Y_k = float(untreated_counts[k])
            treated_sum += X_k    #treated_sum = float(sum(treated_counts[:(k + 1)]))
            untreated_sum += Y_k  #untreated_sum = float(sum(untreated_counts[:(k + 1)]))
            if 0 == treated_sum  or  0 == untreated_sum:
                betas[k] = 0
                thetas[k] = 0
            else:
                Xbit = (X_k / treated_sum)
                Ybit = (Y_k / untreated_sum)
                if Ybit >= 1:
                    betas[k] = 0
                    thetas[k] = 0
                else:
                    betas[k] = max(0, (Xbit - Ybit) / (1 - Ybit))
                    thetas[k] = math.log(1.0 - Ybit) - math.log(1.0 - Xbit)
                    running_c_sum -= math.log(1.0 - betas[k])

        c = running_c_sum
        c_factor = 1.0 / c
        for k in range(n+1):
            thetas[k] = max(c_factor * thetas[k], 0)
        self.betas = betas
        self.thetas = thetas
        self.c = c

    def write_reactivities(self):
        out_path = os.path.join(self.output_folder, "rx.out")
        n = self._target.n
        treated_counts = self._masks[0].counts
        untreated_counts = self._masks[1].counts
        with open(out_path, 'wb') as outfile:
            outfile.write('sequence\trt_start\tfive_prime_offset\tnucleotide\ttreated_mods\tuntreated_mods\tbeta\ttheta\tc\n')
            format_str = "{name}\t{rt}\t".format(name = self._target.name, rt = n - 1) + "{i}\t{nuc}\t{tm}\t{um}\t{b}\t{th}" + "\t{c:.5f}\n".format(c = self.c)
            # TODO: xref https://trello.com/c/OtbxyiYt/23-3-nt-missing-from-reactivities-out
            # looks like we may want this to be range(n), chopping off was unintentional bug of previous version
            for i in range(n - 1):
                outfile.write(format_str.format(i = i,
                                                nuc = self._target.seq[i - 1] if i > 0 else '*',
                                                tm = treated_counts[i],
                                                um = untreated_counts[i],
                                                b = self.betas[i] if i > 0 else '-',
                                                th = self.thetas[i] if i > 0 else '-'))
