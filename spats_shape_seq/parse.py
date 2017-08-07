
import os


# not currently used in spats, but potentially useful for tools
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
        if not first:
            self.reset()
            return
        self.parse([ first, infile.readline(), infile.readline(), infile.readline() ])

    def parse(self, lines):
        self.identifier = lines[0].lstrip('@').split(' ')[0]
        self.sequence = lines[1].rstrip('\r\n')
        self.identifier2 = lines[2].rstrip('\r\n')
        self.quality = lines[3].rstrip('\r\n')

    def write(self, outfile):
        for line in [ self.identifier, self.sequence, self.identifier2, self.quality ]:
            outfile.write(line)
            outfile.write('\n')

    def reverse_complement(self):
        self.sequence = reverse_complement(self.sequence)
        self.quality = self.quality[::-1]


class FastFastqParser(object):

    def __init__(self, r1_path, r2_path):
        self.r1_path = r1_path
        self.r2_path = r2_path

    def pair_length(self):
        with open(self.r1_path, 'rb') as r1_in:
            with open(self.r2_path, 'rb') as r2_in:
                r1_in.readline()
                r1_first = r1_in.readline().strip('\r\n')
                r2_in.readline()
                r2_first = r2_in.readline().strip('\r\n')
                pair_length = len(r1_first)
                if pair_length != len(r2_first):
                    raise Exception("Unexpected pair length mismatch in R1 vs R2: {} / {}".format(pair_length, len(r2_first)))
                return pair_length

    def appx_number_of_pairs(self):
        with open(self.r1_path, 'rb') as r1_in:
            # the +1 is since first records tend to be short, and we'd rather underestimate than overestimate
            frag_len = 1 + len(r1_in.readline()) + len(r1_in.readline()) + len(r1_in.readline()) + len(r1_in.readline())
        return int(float(os.path.getsize(self.r1_path)) / float(frag_len))

    def __enter__(self):
        self.r1_in = open(self.r1_path, 'rb')
        self.r2_in = open(self.r2_path, 'rb')
        self.r1_iter = iter(self.r1_in)
        self.r2_iter = iter(self.r2_in)
        return self

    def __exit__(self, type, value, traceback):
        self.r1_in.close()
        self.r2_in.close()
        self.r1_in = None
        self.r2_in = None
        self.r1_iter = None
        self.r2_iter = None

    def iterator(self, batch_size):
        while True:
            batch = self.iterator_read(batch_size)
            if batch:
                yield batch
            else:
                return

    # kept separate from other read fns for speed
    def iterator_read(self, batch_size):
        pairs = []
        r1_iter = self.r1_iter
        r2_iter = self.r2_iter
        count = 0
        try:
            while count < batch_size:
                R1_id = r1_iter.next() #.split(' ')[0]
                R1_seq = r1_iter.next().rstrip('\n\r')
                r1_iter.next()
                r1_iter.next()
                R2_id = r2_iter.next() #.split(' ')[0]
                R2_seq = r2_iter.next().rstrip('\n\r')
                r2_iter.next()
                r2_iter.next()
                if 0 == count:
                    # good enough to just spot-check this, and improve parsing speed by skipping most of the time
                    R1_id = R1_id.split(' ')[0]
                    R2_id = R2_id.split(' ')[0]
                    if R1_id != R2_id:
                        raise Exception("Malformed input files, id mismatch: {} != {}".format(R1_id, R2_id))
                pairs.append((1, R1_seq, R2_seq, 0))
                count += 1
        except StopIteration:
            pass
        return pairs

    # returns a list of (id, r1, r2), of length <= max_num_pairs, len<max_num_pairs iff eof
    def read(self, max_num_pairs):
        pairs = []
        count = 0
        r1_iter = self.r1_iter
        r2_iter = self.r2_iter
        try:
            while count < max_num_pairs:
                R1_id = r1_iter.next().split(' ')[0]
                R1_seq = r1_iter.next().rstrip('\n\r')
                r1_iter.next()
                r1_iter.next()
                R2_id = r2_iter.next().split(' ')[0]
                R2_seq = r2_iter.next().rstrip('\n\r')
                r2_iter.next()
                r2_iter.next()
                if R1_id != R2_id:
                    raise Exception("Malformed input files, id mismatch: {} != {}".format(R1_id, R2_id))
                pairs.append((R1_id.lstrip('@'), R1_seq, R2_seq))
                count += 1
        except StopIteration:
            pass
        return pairs, count

    # returns a list of (numeric_id, r1, r2, original_id), of length <= max_num_pairs, len<max_num_pairs iff eof
    # separate function in order to keep read() optimized for standard case
    def read_nomask(self, max_num_pairs):
        pairs = []
        count = 0
        r1_iter = self.r1_iter
        r2_iter = self.r2_iter
        try:
            while count < max_num_pairs:
                R1_numeric_id = int(r1_iter.next().strip('@\n\r'))
                R1_seq = r1_iter.next().rstrip('\n\r')
                R1_original_id = r1_iter.next().strip('+\n\r')
                r1_iter.next()
                R2_numeric_id = int(r2_iter.next().strip('@\n\r'))
                R2_seq = r2_iter.next().rstrip('\n\r')
                R2_original_id = r2_iter.next().strip('+\n\r')
                r2_iter.next()
                if R1_numeric_id != R2_numeric_id or R1_original_id != R2_original_id:
                    raise Exception("Malformed NOMASK files, id mismatch: ({},{}) != ({},{})".format(R1_numeric_id, R1_original_id, R2_numeric_id, R2_original_id))
                pairs.append((R1_numeric_id, R1_seq, R2_seq, R1_original_id))
                count += 1
        except StopIteration:
            pass
        return pairs, count


def fasta_parse(target_path):
    pairs = []
    with open(target_path, 'rb') as infile:
        def nextline():
            while True:
                l = infile.readline()
                if len(l) == 0:
                    return l
                l = l.strip('>\n')
                if 0 < len(l):
                    return l
        while True:
            name = nextline()
            if not name:
                break
            seq = nextline()
            if seq:
                seq = seq.upper().replace('U', 'T')
            if name and seq:
                pairs.append((name, seq))
    return pairs


class SamRecord(object):

    def parse(self, line):
        bits = line.split("\t")
        if len(bits) < 6:
            self.identifier = None
            return
        self.identifier = bits[0]
        self.flags = int(bits[1])
        self.target_name = bits[2]
        self.left = int(bits[3]) - 1 # TODO: subtract 1 is required, explain
        self.quality = int(bits[4])
        if self.target_name == '*':
            self.left = -1
            self.target_name = None
            self.right = -1
            return
        lengthPart = bits[5][:-1]
        self.length = int(lengthPart if len(lengthPart) > 0 else 0)
        if self.length > 0 and "M" != bits[5][-1]:
            raise Exception("Expected M on SAM length field, got: {}".format(bits[5]))
        # rest of bits are not used afaict
        self.right = self.left + self.length


class SamParser(object):

    def __init__(self, path, target_map):
        self.sam_path = path
        self.target_map = target_map

    def __enter__(self):
        self.sam_in = open(self.sam_path, 'rb')
        self.sam_iter = iter(self.sam_in)
        return self

    def __exit__(self, type, value, traceback):
        self.sam_in.close()
        self.sam_in = None
        self.sam_iter = None

    # returns a list of (target, site, end, mask, numeric_id) -- the 'mask' is convenience to not have to recreate lists when inserting in DB
    def read(self, max_num_pairs, mask, cotrans = False, single_target = None):
        pairs = []
        count = 0
        sam_iter = self.sam_iter
        r1 = SamRecord()
        r2 = SamRecord()

        def nextline():
            while True:
                l = sam_iter.next()
                if not l.startswith('@'):
                    return l

        try:
            while count < max_num_pairs:
                r1.parse(nextline())
                r2.parse(nextline())
                if not r1.identifier or r1.identifier != r2.identifier:
                    raise Exception("Parse error?") # might just want continue?
                if not r1.target_name or single_target and single_target != r1.target_name:
                    continue
                left = min(r1.left, r2.left)
                right = max(r1.right, r2.right)
                if cotrans:
                    target_length = 20 + int(r1.target_name[r1.target_name.rfind('_')+1:-2])
                    if left >= 0 and left <= target_length and right == target_length:
                        pairs.append((left, right - 20, mask, r1.identifier))
                else:
                    target_length = len(self.target_map[r1.target_name])
                    if left >= 0 and left <= target_length and right == target_length:
                        pairs.append((r1.target_name, left, right, mask, r1.identifier))
                count += 1
        except StopIteration:
            pass
        return pairs, count


# return list of (target, rt_start, site, nuc, treated_count, untreated_count, beta, theta, c)
# all strings
def reactivities_parse(path):
    sites = []
    with open(path, 'rb') as infile:
        infile.readline() # ignore header line
        while True:
            line = infile.readline()
            if not line:
                break
            sites.append(line.split('\t'))
    return sites
