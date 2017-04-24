
import math
import os
import string
import sys


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
        self.identifier = first.rstrip('\r\n')
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

    def relabel_reads(self):
        self.recordNumber += 1
        self.identifier2 += self.identifier[1:]
        self.identifier = "@{}".format(self.recordNumber)

    def strip_and_relabel(self, size, total):
        self.identifier = "@{}".format(total)
        self.sequence = self.sequence[size:] if size > 0 else self.sequence[:size]
        self.quality = self.quality[size:] if size > 0 else self.quality[:size]
        self.identifier2 = "+"


def fastq_transform(input_path, output_path, record_transformer):
    with open(input_path, 'rb') as infile:
        with open(output_path, 'wb') as outfile:
            record = FastqRecord()
            while True:
                record.read(infile)
                if not record.identifier:
                    break
                record_transformer(record)
                if record.identifier:
                    record.write(outfile)

def rev_comp(inputfile, output_revcomp = None):
    output_path = output_revcomp or os.path.splitext(inputfile)[0] + '_rc.fq'
    fastq_transform(inputfile, output_path, lambda x : x.reverse_complement())

def relabel_reads(handle_reads_path, nonhandle_reads_path, output_path):
    transformer = lambda x : x.relabel_reads()
    fastq_transform(handle_reads_path, os.path.join(output_path, "NOMASK_1.fq"), transformer)
    fastq_transform(nonhandle_reads_path, os.path.join(output_path, "NOMASK_2.fq"), transformer)

def filter_reads(handle_reads_path, nonhandle_reads_path, mask_strings, output_folder, relabel = True):
    total = 0
    num_chucked = 0
    masks = [ Mask(m) for m in mask_strings ]
    maskSize = max([ len(m.chars) for m in masks ])
    with open(handle_reads_path, 'rb') as handle_in:
        with open(nonhandle_reads_path, 'rb') as nonhandle_in:
            for mask in masks:
                mask.handle_out = open(os.path.join(output_folder, "{}_1.fq".format(mask.chars)), 'wb')
                mask.nonhandle_out = open(os.path.join(output_folder, "{}_2.fq".format(mask.chars)), 'wb')
            try:
                handle_record = FastqRecord()
                nonhandle_record = FastqRecord()

                while True:
                    handle_record.read(handle_in)
                    nonhandle_record.read(nonhandle_in)
                    if not handle_record.identifier or not nonhandle_record.identifier:
                        break
                    total += 1
                    matched_mask = None
                    if min(len(handle_record.sequence), len(nonhandle_record.sequence)) > maskSize + 1:
                        for mask in masks:
                            if matches_mask(handle_record.sequence, mask.values):
                                matched_mask = mask
                                break
                    if not matched_mask:
                        num_chucked += 1
                    else:
                        sz = len(mask.chars)
                        handle_record.strip_and_relabel(sz, total)
                        handle_record.write(matched_mask.handle_out)
                        nonhandle_record.strip_and_relabel(-sz, total)
                        nonhandle_record.write(matched_mask.nonhandle_out)
            finally:
                for mask in masks:
                    mask.handle_out.close()
                    mask.nonhandle_out.close()

    print "Kept {} of {} reads".format(total - num_chucked, total)

class SamRecord(object):

    def parse(self, line):
        bits = line.split("\t")
        if len(bits) < 6:
            self.identifier = None
            return
        self.identifier = bits[0]
        self.flags = int(bits[1])
        # bits[2] is the target name
        self.left = int(bits[3]) - 1 # TODO: subtract 1 is required, explain
        self.quality = int(bits[4])
        lengthPart = bits[5][:-1]
        self.length = int(lengthPart if len(lengthPart) > 0 else 0)
        if self.length > 0 and "M" != bits[5][-1]:
            raise Exception("Expected M on SAM length field, got: {}".format(bits[5]))
        # rest of bits are not used afaict
        self.right = self.left + self.length

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

def get_counts(sam_path, target_length):
    num_sites = target_length + 1
    counts = [ 0 for x in range(num_sites) ] # TODO: numpy.empty(num_sites, dtype=int) is better
    num_fragments = 0
    num_chucked = 0
    r1 = SamRecord()
    r2 = SamRecord()
    with open(sam_path, 'rb') as infile:
        while True:
            line = infile.readline()
            if 0 == len(line):
                break
            if line.startswith('@'):
                continue
            r1.parse(line)
            r2.parse(infile.readline())
            if not r1.identifier or r1.identifier != r2.identifier:
                raise Exception("Parse error?") # might just want continue?
            num_fragments += 1
            left = min(r1.left, r2.left)
            right = max(r1.right, r2.right)
            if left >= 0 and left < num_sites and right == target_length:
                # JBL - only register this fragment if the left read is within the sequence and the right read aligns with the end of the RNA (or RNA subsequence)
                counts[left] += 1
            else:
                num_chucked += 1
    return counts, num_fragments, num_fragments - num_chucked

def compute_from_counts(num_sites, treated_counts, untreated_counts):
    n = int(num_sites)
    betas = [ 0 for x in range(n+1) ]
    thetas = [ 0 for x in range(n+1) ]
    treated_sum = 0.0    # keeping a running sum for both
    untreated_sum = 0.0  # channels is much faster
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
        thetas[k] = c_factor * thetas[k]

    return betas, thetas, c

def compute_profiles(target_path, treated_path, untreated_path, output_dir):
    pairs = fasta_parse(target_path)

    # TODO: handle multiple sequences?

    name = pairs[0][0]
    seq = pairs[0][1]
    n = len(seq)
    treated_counts, treated_total, treated_kept = get_counts(treated_path, n)
    untreated_counts, untreated_total, untreated_kept = get_counts(untreated_path, n)
    print "Processed {} properly paired fragments, kept {}/{} ({:.2f}%) treated, {}/{} ({:2f}%) untreated".format(treated_total + untreated_total,
                                                                                                                  treated_kept,
                                                                                                                  treated_total,
                                                                                                                  (100.0 * float(treated_kept)) / float(treated_total),
                                                                                                                  untreated_kept,
                                                                                                                  untreated_total,
                                                                                                                  (100.0 * float(untreated_kept)) / float(untreated_total))

    betas, thetas, c = compute_from_counts(n, treated_counts, untreated_counts)
    write_reactivities(os.path.join(output_dir, "reactivities.out"), name, seq, treated_counts, untreated_counts, betas, thetas, c)

def write_reactivities(out_path, name, seq, treated_counts, untreated_counts, betas, thetas, c):
    n = len(seq)
    with open(out_path, 'wb') as outfile:
        outfile.write('sequence\trt_start\tfive_prime_offset\tnucleotide\ttreated_mods\tuntreated_mods\tbeta\ttheta\tc\n')
        format_str = "{name}\t{rt}\t".format(name = name, rt = n - 1) + "{i}\t{nuc}\t{tm}\t{um}\t{b}\t{th}" + "\t{c:.5f}\n".format(c = c)
        # TODO: xref https://trello.com/c/OtbxyiYt/23-3-nt-missing-from-reactivities-out
        # looks like we may want this to be range(n), chopping off was unintentional bug of previous version
        for i in range(n - 1):
            outfile.write(format_str.format(i = i,
                                            nuc = seq[i - 1] if i > 0 else '*',
                                            tm = treated_counts[i],
                                            um = untreated_counts[i],
                                            b = betas[i] if i > 0 else '-',
                                            th = thetas[i] if i > 0 else '-'))



###############################################################

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
        self.show_progress = True
        self.write_intermediate_outputs = False
        self.show_id_to_site = False
        self.quiet = False
        self.debug = False
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
        if self.write_intermediate_outputs:
            self.combined_r1_out = open(os.path.join(self.output_folder, "comb_R1.fq"), 'wb')
            self.combined_r2_out = open(os.path.join(self.output_folder, "comb_R2.fq"), 'wb')

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
        self._n = len(self._target.seq)

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
            mask.counts = [ 0 for x in range(self._n + 1) ] # TODO: numpy.empty(num_sites, dtype=int) is better
            mask.total = 0
            mask.kept = 0
            if self.write_intermediate_outputs:
                mask.handle_out = open(os.path.join(self.output_folder, "{}_1.fq".format(mask.chars)), 'wb')
                mask.nonhandle_out = open(os.path.join(self.output_folder, "{}_2.fq".format(mask.chars)), 'wb')
        self._masks = masks
        self._maskSize = max([ len(m.chars) for m in masks ])

    def _cleanup(self):
        if not self.write_intermediate_outputs:
            return
        for mask in self._masks:
            mask.handle_out.close()
            mask.nonhandle_out.close()
        self.combined_r1_out.close()
        self.combined_r2_out.close()

    def _match_mask(self, seq):
        if len(seq) > self._maskSize + 1:
            for mask in self._masks:
                if matches_mask(seq, mask.values):
                    return mask
        return None

    def _DBG(self, stuff):
        if self.debug:
            print stuff

    def _process_pair(self, r1_record, r2_record, numeric_pair_id = 0):
        mask = self._match_mask(r1_record.sequence)
        if not mask:
            return False

        n = self._n
        mask.total += 1
        self.debug = True#("10157" in r1_record.identifier)# or "14617" in r1_record.identifier)# or "19869" in r1_record.identifier or "22589" in r1_record.identifier)
        self._DBG("> processing " + r1_record.identifier + "\n  --> " + r1_record.sequence + " , " + r2_record.sequence)
        rev1 = reverse_complement(r1_record.sequence)
        self._DBG(rev1)
        r1len, r2len = map(len, [r1_record.sequence, r2_record.sequence])
        start1, len1, index1 = self._target.find_partial(rev1, 8)
        start2, len2, index2 = self._target.find_partial(r2_record.sequence, 8)
        self._DBG([start1, len1, index1, "--", start2, len2, index2, r2len])
        #print ["**  " if r1_record.identifier in self._comb_ids else "  ", start1, len1, index1, "--", start2, len2, index2, "-" if index2 is None else index2 + len2, r1_record.identifier]
        if not len1 or not len2:
            self._DBG("len failure")
            return True

        if start1 > 4:
            # don't try trimming with less than 4, assume it's just a handle
            tail_start = (start1 - 4) if (start1 > 8) else start1
            rest = r1_record.sequence[-tail_start:]
            self._DBG("R1 trim {}: {}".format(r1_record.identifier, rest))
            rstart, rlen, aindex = self._adapter_b_index.find_partial(rest, 4)
            self._DBG((rstart, rlen, aindex))
            if not rlen:
                return True
            shifted = r1_record.sequence[-(tail_start+aindex):]
            if self.adapter_b.startswith(shifted):
                # ok, this is a good trim
                self._DBG("Found adapter_b trim L{} on: {}".format(tail_start + aindex, r1_record.identifier))
                new_start = tail_start + aindex + 4 # add 4 for the handle
                delta = new_start - start1
                start1 = new_start
                len1 -= delta
                index1 += delta
                self._DBG("new start1={}, len1={}, index1={}".format(start1, len1, index1))
            else:
                self._DBG("Can't trim R1 on: {}".format(r1_record.identifier))
                return True

        if start2 + len2 != r2len:
            rest = r2_record.sequence[start2 + len2:]
            self._DBG("R2 off: " + rest)
            # if the the leftover is <= 4, just assume it's the handle and move on; so only need to check for adapter if len > 4
            if len(rest) > 4:
                possible_adapter = reverse_complement(rest[4:])
                self._DBG(possible_adapter)
                if not self.adapter_t.endswith(possible_adapter):
                    self._DBG("untrimmable R2 failure")
                    return True
                else:
                    self._DBG("accepting trimmed R2")

        if not set(r1_record.sequence + r2_record.sequence) <= set('ACGT'):
            # TODO: we should keep these -- but current adapter_trimmer does not
            self._DBG("indeterminate sequence failure")
            return True

        if start2 > 0 and index2 > 0:
            self._DBG("inside 5S failure on R2 for: {}".format(r1_record.identifier))
            return True

        if start2 > 0:
            self._DBG("start too soon failure on R2 for: {}".format(r1_record.identifier))
            return True

        if self.write_intermediate_outputs:
            r1_record.write(self.combined_r1_out)
            r2_record.write(self.combined_r2_out)

        left = min(index1, index2)
        right = max(index1 + len1, index2 + len2)
        self._DBG([left, right, n])

        # JBL - only register this fragment if the left read is within the sequence
        # and the right read aligns with the end of the RNA (or RNA subsequence)
        if left < 0 or left > n or right != n:
            self._DBG("LR failure")
            return True

        self._DBG("  ===> KEPT {}-{}".format(left, right))
        mask.kept += 1
        mask.counts[left] += 1

        if self.show_id_to_site:
            print "{} --> {} ({})".format(r1_record.identifier, left, mask.chars)

        if not self.write_intermediate_outputs:
            return True

        self._DBG(r1_record.identifier)
        if len1 < r1len:
            self._DBG(rev1[start1 : start1 + len1])
            self._DBG(reverse_complement(rev1[start1 : start1 + len1]))
            r1_record.sequence = reverse_complement(rev1[start1 : start1 + len1])
            r1_record.quality = r1_record.quality[r1len - start1 - len1 : r1len - start1]
        if len2 < len(r2_record.sequence):
            r2_record.sequence = r2_record.sequence[start2 : start2 + len2]
            r2_record.quality = r2_record.quality[start2 : start2 + len2]
        r1_record.identifier = r2_record.identifier = "@{}".format(numeric_pair_id)
        #if self.debug: exit(1)
        if self.write_intermediate_outputs:
            r1_record.write(mask.handle_out)
            r2_record.write(mask.nonhandle_out)
        return True

    def get_counts(self):
        total_pairs = 0
        chucked_pairs = 0

        with open(self.data_r1_path, 'rb') as r1_in:
            with open(self.data_r2_path, 'rb') as r2_in:
                r1_record = FastqRecord()
                r2_record = FastqRecord()
                while True:
                    r1_record.read(r1_in)
                    r2_record.read(r2_in)
                    if not r1_record.identifier or not r2_record.identifier:
                        break

                    total_pairs += 1
                    if self.show_progress and 0 == total_pairs % 20000:
                        sys.stdout.write('.')
                        sys.stdout.flush()

                    if not self._process_pair(r1_record, r2_record, total_pairs):
                        chucked_pairs += 1

        self.total_pairs = total_pairs
        self.chucked_pairs = chucked_pairs
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
        n = self._n
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
        n = self._n
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


def spats(target_path, r1_path, r2_path, output_folder):
    s = Spats(target_path, r1_path, r2_path, output_folder)
    s.write_intermediate_outputs = True
    s.show_id_to_site = True
    s.run()

def make_subset(r1_path, r2_path, id_list_path, output_folder):
    ids = set(open(id_list_path, 'rb').read().split('\n'))
    n = 0
    with open(os.path.join(output_folder, "filtered_R1.fq"), 'wb') as r1_out:
        with open(os.path.join(output_folder, "filtered_R2.fq"), 'wb') as r2_out:
            with open(r1_path, 'rb') as r1_in:
                with open(r2_path, 'rb') as r2_in:
                    r1_record = FastqRecord()
                    r2_record = FastqRecord()
                    while True:
                        r1_record.read(r1_in)
                        r2_record.read(r2_in)
                        if not r1_record.identifier or not r2_record.identifier:
                            break
                        rid = r1_record.identifier.split(' ')[0][1:]
                        if rid in ids:
                            n += 1
                            r1_record.write(r1_out)
                            r2_record.write(r2_out)
    print "Filtered {} records.".format(n)

def id_to_site(nomask_path, treated_sam, untreated_sam, target_length):
    nid_to_id = {}
    with open(nomask_path, 'rb') as infile:
        r = FastqRecord()
        while True:
            r.read(infile)
            if not r.identifier:
                break
            nid = int(r.identifier.lstrip('@'))
            if nid < 525 and nid > 420:
                print "{} -> {}".format(nid, r.identifier2.lstrip("+"))
            nid_to_id[nid] = r.identifier2.lstrip('+')
    num_sites = target_length + 1
    r1 = SamRecord()
    r2 = SamRecord()
    for sam_path in [ treated_sam, untreated_sam ]:
        print sam_path
        with open(sam_path, 'rb') as infile:
            while True:
                line = infile.readline()
                if 0 == len(line):
                    break
                if line.startswith('@'):
                    continue
                r1.parse(line)
                r2.parse(infile.readline())
                if not r1.identifier or r1.identifier != r2.identifier:
                    raise Exception("Parse error?") # might just want continue?
                left = min(r1.left, r2.left)
                right = max(r1.right, r2.right)
                if left >= 0 and left < num_sites and right == target_length:
                    print nid_to_id[int(r1.identifier)] + " --> " + str(left)

    
        
def old_id_to_site(sam_path, target_length):
    num_sites = target_length + 1
    r1 = SamRecord()
    r2 = SamRecord()
    with open(sam_path, 'rb') as infile:
        while True:
            line = infile.readline()
            if 0 == len(line):
                break
            if line.startswith('@'):
                continue
            r1.parse(line)
            r2.parse(infile.readline())
            if not r1.identifier or r1.identifier != r2.identifier:
                print r1.identifier
                print r2.identifier
                raise Exception("Parse error?") # might just want continue?
            left = min(r1.left, r2.left)
            right = max(r1.right, r2.right)
            if left >= 0 and left < num_sites and right == target_length:
                print r1.identifier + " --> " + str(left)

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "rc":
        print reverse_complement(sys.argv[2])
    else:
        print "Unknown cmd: {}".format(sys.argv[1])
