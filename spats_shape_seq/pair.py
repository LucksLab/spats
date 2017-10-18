
from sequence import Sequence

class Pair(object):

    def __init__(self):
        self.reset()

    def reset(self):
        self.identifier = None
        self.r1 = Sequence()
        self.r2 = Sequence()
        self.mask = None
        self.target = None
        self.site = None
        self.failure = None
        self.multiplicity = 1
        self.tags = None
        self._end = -1
        self.mutations = None

    def set_from_data(self, identifier, r1_seq, r2_seq, multiplicity = 1):
        self.reset()
        self.identifier = identifier
        self.r1.set_seq(r1_seq)
        self.r2.set_seq(r2_seq)
        self.multiplicity = multiplicity

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
        self.r1._ltrim = 4

    def check_mutations(self):
        mutations = set()
        for seq in [ self.r1, self.r2 ]:
            for err_index in seq.match_errors:
                # +1 for M_j indexing convention, xref https://trello.com/c/2qIGo9ZR/201-stop-map-mutation-indexing-convention
                mutations.add(seq.match_index + err_index + 1)
        if mutations:
            self.mutations = list(mutations)

    @property
    def matched(self):
        return (self.target and self.r1.match_len and self.r2.match_len)

    @property
    def has_site(self):
        return (self.target and self.site is not None)

    @property
    def left(self):
        return self.r2.left

    @property
    def right(self):
        return self.r1.right

    @property
    def end(self):
        return self.right if self._end == -1 else self._end

    @end.setter
    def end(self, val):
        self._end = val
