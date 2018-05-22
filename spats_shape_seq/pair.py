
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

    def check_mutation_quality(self, minimum_quality_score):
        if not minimum_quality_score or not self.mutations:
            return 0
        removed = []
        r1_start = self.r1.match_index
        r2_start = self.r2.match_index
        seq_len = self.r2.original_len
        for mut in self.mutations:
            q1 = q2 = None
            if mut < r2_start + seq_len + 1:
                q1 = ord(self.r2.quality[mut - r2_start - 1])
            if mut > r1_start:
                q2 = ord(self.r1.quality[::-1][mut - r1_start - 1])
            # note: this assumes agreement. TODO: handle disagreement, xref check_overlap below
            q = max(q1, q2) if (q1 and q2) else (q1 if q1 else q2)
            if q < minimum_quality_score:
                removed.append(mut)
        for mut in removed:
            self.mutations.remove(mut)
        return len(removed)

    def check_overlap(self):
        # note: in the case that overlaps disagree, we may decide it one way or the other via quality
        # xref https://trello.com/c/35mBHvPA/197-stop-map-r1-r2-disagree-case
        r2_match_len = self.r2.match_len
        overlap_index = self.r1.match_index
        overlap_len = self.r2.match_index + r2_match_len - overlap_index
        if overlap_len > 0:
            #print('O: {} / {}'.format(self.r1.reverse_complement[:overlap_len], self.r2.subsequence[r2_match_len-overlap_len:r2_match_len]))
            if self.r1.reverse_complement[:overlap_len] != self.r2.subsequence[r2_match_len-overlap_len:r2_match_len]:
                return False
        return True

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
    def length(self):
        return self.r1.right - self.r2.left

    @property
    def end(self):
        return self.right if self._end == -1 else self._end

    @end.setter
    def end(self, val):
        self._end = val
