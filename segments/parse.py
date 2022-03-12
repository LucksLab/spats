import os


class FastFastqParser(object):

    def __init__(self, r1_path, r2_path, parse_quality = False):
        self.r1_path = r1_path
        self.r2_path = r2_path
        self.parse_quality = parse_quality

    def pair_length(self):
        with open(self.r1_path, 'rb') as r1_in:
            with open(self.r2_path, 'rb') as r2_in:
                r1_in.readline()
                r1_first = r1_in.readline().strip('\r\n')
                r2_in.readline()
                r2_first = r2_in.readline().strip('\r\n')
                pair_length = len(r1_first)
                if pair_length != len(r2_first):
                    print("Warning: pair length mismatch in R1 vs R2: {} / {}".format(pair_length, len(r2_first)))
                    return -1
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
        include_quality = self.parse_quality
        try:
            while count < batch_size:
                R1_id = next(r1_iter).decode() #.split(' ')[0]
                R1_seq = next(r1_iter).decode().rstrip('\n\r')
                next(r1_iter)
                R1_q = next(r1_iter).decode()
                R2_id = next(r2_iter).decode() #.split(' ')[0]
                R2_seq = next(r2_iter).decode().rstrip('\n\r')
                next(r2_iter)
                R2_q = next(r2_iter).decode()
                if 0 == count:
                    # good enough to just spot-check this, and improve parsing speed by skipping most of the time
                    R1_id = R1_id.split(' ')[0]
                    R2_id = R2_id.split(' ')[0]
                    if R1_id != R2_id:
                        raise Exception("Malformed input files, id mismatch: {} != {}".format(R1_id, R2_id))
                if include_quality:
                    pairs.append((1, R1_seq, R2_seq, R1_id.split(' ')[0], R1_q.rstrip('\n\r'), R2_q.rstrip('\n\r')))
                else:
                    pairs.append((1, R1_seq, R2_seq, str(count)))
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
                R1_id = next(r1_iter).decode().split(' ')[0]
                R1_seq = next(r1_iter).decode().rstrip('\n\r')
                next(r1_iter)
                next(r1_iter)
                R2_id = next(r2_iter).decode().split(' ')[0]
                R2_seq = next(r2_iter).decode().rstrip('\n\r')
                next(r2_iter)
                next(r2_iter)
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
                R1_original_id = next(r1_iter).decode().strip('+\n\r')
                R1_numeric_id = int(R1_original_id.split(':')[0].strip("@M"))
                R1_seq = next(r1_iter).decode().rstrip('\n\r')
                next(r1_iter)
                R2_original_id = next(r2_iter).decode().strip('+\n\r')
                R2_numeric_id = int(R2_original_id.split(':')[0].strip("@M"))
                R2_seq = next(r2_iter).decode().rstrip('\n\r')
                next(r2_iter)
                if R1_numeric_id != R2_numeric_id or R1_original_id != R2_original_id:
                    raise Exception("Malformed NOMASK files, id mismatch: ({},{}) != ({},{})".format(R1_numeric_id, R1_original_id, R2_numeric_id, R2_original_id))
                pairs.append((R1_numeric_id, R1_seq, R2_seq, R1_original_id))
                count += 1
        except StopIteration:
            pass
        return pairs, count
