
from sequence import Sequence

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
