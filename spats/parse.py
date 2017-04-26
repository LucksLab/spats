
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

    # returns a list of (id, r1, r2), of length <= max_num_pairs, len<max_num_pairs iff eof
    def read(self, max_num_pairs):
        pairs = []
        count = 0
        r1_iter = self.r1_iter
        r2_iter = self.r2_iter
        try:
            while count < max_num_pairs:
                R1_id = r1_iter.next()
                R1_seq = r1_iter.next()
                r1_iter.next()
                r1_iter.next()
                r2_iter.next()
                R2_seq = r2_iter.next()
                r2_iter.next()
                r2_iter.next()
                pairs.append((R1_id, R1_seq, R2_seq))
                count += 1
        except StopIteration:
            pass
        return pairs, count


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

