
from sequence import Sequence

class Pair(object):

    def __init__(self):
        self.reset()

    def __repr__(self):
        return 'T={}/E={}/S={}/m={}'.format(self.target.name if self.target else None, self.end, self.site, self.mutations)

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
        self.removed_mutations = None
        self.linker = None
        self.edge_mut = None
        self.dumbbell = None

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
        self.r1._ltrim = mask.length()

    def check_mutations(self):
        mutations = set()
        for seq in [ self.r1, self.r2 ]:
            for err_index in seq.match_errors:
                # +1 for M_j indexing convention, xref https://trello.com/c/2qIGo9ZR/201-stop-map-mutation-indexing-convention
                mutations.add(seq.match_index + err_index + 1)
        if mutations:
            self.mutations = list(mutations)

    def check_prefix_quality(self, prefix_length, minimum_quality_score):
        if minimum_quality_score is None:
            return True
        assert(prefix_length > 0)
        # prefix is always on the left end of R2
        prefix_quality = self.r2.quality[0:prefix_length]
        for q in prefix_quality:
            if ord(q) < minimum_quality_score + ord('!'):
                return False
        return True

    def check_mutation_quality(self, minimum_quality_score):
        if minimum_quality_score is None  or  not self.mutations or (not self.r1.quality and not self.r2.quality):
            return 0
        removed = []
        r1_start = self.r1.match_index
        r2_start = self.r2.match_index
        seq_len = self.r2.original_len
        min_quality = minimum_quality_score + ord('!')
        for mut in self.mutations:
            q1 = q2 = None
            nt1 = nt2 = None
            if mut < r2_start + seq_len + 1 - self.r2._ltrim:
                idx = mut - r2_start + self.r2._ltrim - 1
                if idx >= 0:
                    q2 = ord(self.r2.quality[idx])
                    nt2 = self.r2.original_seq[idx]
            if mut > r1_start:
                idx = mut - r1_start - 1
                if idx >= 0:
                    q1 = ord(self.r1.quality[::-1][idx + self.r1._rtrim])
                    nt1 = self.r1.reverse_complement[idx]
            if q1 and q2:
                # check for agreement
                if nt1 == nt2:
                    q = max(q1, q2)
                else:
                    # xref https://trello.com/c/usT0vTiG/308-discordant-mutations-case-handling-and-unit-tests
                    #    and associated test cases
                    ref = self.target.seq[mut]     # TODO STEVE:  why isn't this mut-1?  (unit tests still pass)
                    if q1 < min_quality and q2 >= min_quality:
                        if nt2 == ref:
                            q = q1 # signals to ignore the mutation from low-quality q1
                        else:
                            q = q2 # signals to count the mutation from high-quality q2
                    elif q2 < min_quality and q1 >= min_quality:
                        if nt1 == ref:
                            q = q2 # signals to ignore the mutation from low-quality q2
                        else:
                            q = q1 # signals to count the mutation from high-quality q1
                    elif q1 < min_quality and q2 < min_quality:
                        q = max(q1, q2) # both quality are low, we're going to throw it away
                    else:
                        # both high quality and disagree, signal low quality which has effect of ignoring the mut
                        q = 0
            elif q1:
                q = q1
            elif q2:
                q = q2
            if q < min_quality:
                removed.append(mut)
        for mut in removed:
            self.mutations.remove(mut)
        if removed:
            self.removed_mutations = removed
        return len(removed)

    def check_overlap(self):
        # note: in the case that overlaps disagree, we may decide it one way or the other via quality
        # xref https://trello.com/c/35mBHvPA/197-stop-map-r1-r2-disagree-case
        r2_match_len = self.r2.match_len
        overlap_index = max(self.r1.match_index, self.r2.match_index)
        overlap_len = self.r2.match_index + r2_match_len - overlap_index
        if overlap_len > 0:
            r1start = max(self.r2.match_index - self.r1.match_index, 0)
            r1_part = self.r1.reverse_complement[r1start:r1start+overlap_len]
            r2_part = self.r2.subsequence[r2_match_len-overlap_len:r2_match_len]
            #print('O: {} / {}'.format(r1_part, r2_part))
            if r1_part != r2_part:
                return False
        return True

    @property
    def matched(self):
        return (self.target and self.r1.match_len and self.r2.match_len)

    @property
    def has_site(self):
        return bool(self.target and self.site is not None)

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
