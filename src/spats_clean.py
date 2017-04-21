import math
import os
import string
import sys


class SpatsConfig(object):

    def __init__(self):
        self.minimum_target_match_length = 8
        self.debug = False
        self.allow_indeterminate = False
        self.show_progress = True
        self.show_id_to_site = False

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
        self.seq = seq
        self.length = len(seq) if seq else None
        self._rc = None
        self.match_start = None
        self.match_len = None
        self.match_index = None
        self.left = None
        self.right = None

    def set_seq(self, seq):
        self._reset(seq)

    @property
    def reverse_complement(self):
        if not self._rc:
            self._rc = reverse_complement(self.seq)
        return self._rc

    def find_in_target(self, target, reverse_complement = False):
        seq = self.reverse_complement if reverse_complement else self.seq
        self.match_start, self.match_len, self.match_index = target.find_partial(seq, spats_config.minimum_target_match_length)
        if self.match_len:
            self.left = self.match_index
            self.right = self.match_index + self.match_len


class Pair(object):

    def __init__(self):
        self.reset()

    def reset(self):
        self.identifier = None
        self.r1 = Sequence()
        self.r2 = Sequence()
        self.site = None

    def set_from_data(self, identifier, r1_seq, r2_seq):
        self.identifier = identifier
        self.r1.set_seq(r1_seq)
        self.r2.set_seq(r2__seq)

    def set_from_records(self, r1_record, r2_record):
        if not r1_record.identifier or r1_record.identifier != r2_record.identifier:
            raise Exception("Invalid record IDs for pair: {}, {}".format(r1_record.identifier, r2_record.identifier))
        self.reset()
        self.identifier = r1_record.identifier
        self.r1.set_seq(r1_record.sequence)
        self.r2.set_seq(r2_record.sequence)

    def is_determinate(self):
        return set(self.r1.seq + self.r2.seq) <= set('ACGT')

    @property
    def matched(self):
        return (self.r1.match_len and self.r2.match_len)

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
def longest_match(s1, range1, s2, range2):
    left1 = range1[0]
    left2 = range2[0]
    right1 = left1 + range1[1]
    right2 = left2 + range2[1]
    if s1[left1:right1] != s2[left2:right2]:
        raise Exception("longest_match must already start with a match")
    while left1 >= 0 and left2 >= 0 and 0 != (char_to_mask[s1[left1]] & char_to_mask[s2[left2]]):
        left1 -= 1
        left2 -= 1
    while right1 <= len(s1) and right2 <= len(s2) and 0 != (char_to_mask[s1[right1 - 1]] & char_to_mask[s2[right2 - 1]]):
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
        for i in range(len(seq) - word_len):
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
        for site in check_sites:
            site_key = query[site:site+word_len]
            #print "CS: {}, {}".format(site, site_key)
            for index in self._index.get(site_key, []):
                #print "GOT: " + str(index)
                left, right = longest_match(query, (site, word_len), self.seq, (index, word_len))
                #print "extends: <--{}, -->{}".format(left, right)
                total_len = left + right + word_len
                if total_len >= min_len:
                    return site - left, total_len, index - left
        return None, None, None


class Spats(object):

    def __init__(self, target_path, data_r1_path, data_r2_path, output_folder):
        self.target_path = target_path
        self.data_r1_path = data_r1_path
        self.data_r2_path = data_r2_path
        self.output_folder = output_folder

        # user-configurable parameters
        self.masks = [ 'RRRY', 'YYYR' ]
        self.quiet = False
        self.adapter_t = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"  # JJB: shows revcomped on the end of R2, from the back
        self.adapter_b = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"                          # JJB: shows as-is on the end of R1, from the front

        # private vars
        self._targets = None
        self._masks = None
        self._adapter_t_index = None

        self._comb_ids = open("/Users/jbrink/mos/tasks/1RwIBa/tmp/5sq_dev/tmp/comb_R1.ids", 'rb').read().split('\n')

    def setup(self):
        self._indexTargets()
        self._indexAdapters()
        self._setupMasks()

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

    def _indexAdapters(self):
        if self._adapter_t_index:
            return
        self._adapter_t_index = Target("adapter_t", self.adapter_t, index_word_length = 4)
        self._adapter_t_index.index()
        self._adapter_b_index = Target("adapter_b", self.adapter_b, index_word_length = 4)
        self._adapter_b_index.index()

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

    def _cleanup(self):
        pass

    def _match_mask(self, pair):
        if len(pair.r1.seq) <= self._maskSize:
            _warn("sequence {} is too short for masking".format(seq))
            return
        for mask in self._masks:
            if matches_mask(pair.r1.seq, mask.values):
                mask.total += 1
                pair.mask = mask
                return

    def _find_matches(self, pair):
        rev1 = pair.r1.reverse_complement
        _debug("rc(R1): {}".format(rev1))
        pair.r1.find_in_target(self._target, reverse_complement = True)
        pair.r2.find_in_target(self._target)
        _debug([pair.r1.match_start, pair.r1.match_len, pair.r1.match_index, "--", pair.r2.match_start, pair.r2.match_len, pair.r2.match_index])

    def _trim_adapters(self, pair):
        # TODO: needs a lot of work
        start1 = pair.r1.match_start
        if start1 > 4:
            # don't try trimming with less than 4, assume it's just a handle
            tail_start = (start1 - 4) if (start1 > 8) else start1
            rest = pair.r1.seq[-tail_start:]
            _debug("R1 trim {}: {}".format(pair.identifier, rest))
            rstart, rlen, aindex = self._adapter_b_index.find_partial(rest, 4)
            _debug((rstart, rlen, aindex))
            if not rlen:
                _debug("No adapter found to trim")
                return False
            shifted = pair.r1.seq[-(tail_start+aindex):]
            if self.adapter_b.startswith(shifted):
                # ok, this is a good trim
                _debug("Found adapter_b trim L{} on: {}".format(tail_start + aindex, pair.identifier))
                pair.r1.match_start = tail_start + aindex + 4 # add 4 for the handle
                delta = pair.r1.match_start - start1
                pair.r1.match_len -= delta
                pair.r1.match_index += delta
                _debug("new start1={}, len1={}, index1={}".format(pair.r1.match_start, pair.r1.match_len, pair.r1.match_index))
            else:
                _debug("Can't trim R1 on: {}".format(pair.identifier))
                return False

        if pair.r2.match_start + pair.r2.match_len != pair.r2.length:
            rest = pair.r2.seq[pair.r2.match_start + pair.r2.match_len:]
            _debug("R2 off: " + rest)
            # if the the leftover is <= 4, just assume it's the handle and move on; so only need to check for adapter if len > 4
            if len(rest) > 4:
                possible_adapter = reverse_complement(rest[4:])
                _debug(possible_adapter)
                if not self.adapter_t.endswith(possible_adapter):
                    _debug("untrimmable R2 failure")
                    return False
                else:
                    _debug("accepting trimmed R2")

        return True

    def process_pair(self, pair):

        #spats_config.debug = ("10157" in r1_record.identifier)# or "14617" in r1_record.identifier)# or "19869" in r1_record.identifier or "22589" in r1_record.identifier)

        _debug("> processing " + pair.identifier + "\n  --> " + pair.r1.seq + " , " + pair.r2.seq)

        self._match_mask(pair)
        if not pair.mask:
            _debug("mask failure")
            return

        self._find_matches(pair)
        if not pair.matched:
            _debug("no match")
            return

        if not self._trim_adapters(pair):
            _debug("adapter trim failure")
            return

        if not spats_config.allow_indeterminate  and  not pair.is_determinate():
            _debug("indeterminate sequence failure")
            return

        if pair.r2.match_start > 0 and pair.r2.match_index > 0:
            _debug("inside 5S failure on R2 for: {}".format(pair.identifier))
            return

        if pair.r2.match_start > 0:
            _debug("start too soon failure on R2 for: {}".format(pair.identifier))
            return

        # JBL - only register this fragment if the left read is within the sequence
        # and the right read aligns with the end of the RNA (or RNA subsequence)
        n = self._target.n
        _debug([pair.left, pair.right, n])
        if pair.left < 0 or pair.left > n or pair.right != n:
            _debug("LR failure")
            return

        _debug("  ===> KEPT {}-{}".format(pair.left, pair.right))
        pair.register_count()

        if spats_config.show_id_to_site:
            print "{} --> {} ({})".format(pair.identifier, pair.site, pair.mask.chars)


    def get_counts(self):
        total_pairs = 0
        chucked_pairs = 0
        process_pairs = 0

        with open(self.data_r1_path, 'rb') as r1_in:
            with open(self.data_r2_path, 'rb') as r2_in:
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
                    elif if pair.site != None:
                        processed_pairs +=1

                    if spats_config.show_progress and 0 == total_pairs % 20000:
                        sys.stdout.write('.')
                        sys.stdout.flush()

        self.total_pairs = total_pairs
        self.chucked_pairs = chucked_pairs
        self.processed_pairs = processed_pairs
        self._cleanup()
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

    def run(self):
        self.setup()
        self.get_counts()
        self.compute_profiles()
        self.write_reactivities()
