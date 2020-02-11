
from sequence import Sequence
from util import reverse_complement
from mask import PLUS_PLACEHOLDER, MINUS_PLACEHOLDER


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
        self.mask_label = ""
        self._target = None
        self.site = None
        self.interesting = False
        self.failure = None
        self.multiplicity = 1
        self.tags = None
        self._end = -1
        self.mutations = None
        self.removed_mutations = None
        self.linker = None
        self.edge_mut = None
        self.edge_indel = False
        self.ambig_indel = False
        self.dumbbell = None
        self._fully_matched = False
        self._indels_match = None

    def set_from_data(self, identifier, r1_seq, r2_seq, multiplicity = 1):
        self.reset()
        self.identifier = identifier
        self.r1.set_seq(r1_seq, True)
        self.r2.set_seq(r2_seq)
        self.multiplicity = multiplicity

    def set_from_records(self, r1_record, r2_record):
        if not r1_record.identifier or r1_record.identifier != r2_record.identifier:
            raise Exception("Invalid record IDs for pair: {}, {}".format(r1_record.identifier, r2_record.identifier))
        self.reset()
        self.identifier = r1_record.identifier
        self.r1.set_seq(r1_record.sequence, True)
        self.r2.set_seq(r2_record.sequence)

    def is_determinate(self):
        return set(self.r1.original_seq + self.r2.original_seq) <= set('ACGT')

    def set_mask(self, mask):
        self.mask = mask
        self.r1.ltrim = mask.length()
        self.mask_label = mask.empty_place_holder if mask.empty_place_holder else mask.chars

    def check_mutations(self):
        mutations = set()
        for seq in [ self.r1, self.r2 ]:
            for err_index in seq.match_errors:
                # +1 for M_j indexing convention, xref https://trello.com/c/2qIGo9ZR/201-stop-map-mutation-indexing-convention
                mutations.add(seq.match_index + err_index + 1)
        if mutations:
            self.mutations = list(mutations)

    def check_prefix_quality(self, prefix_length, minimum_quality_score, offset = 0):
        if minimum_quality_score is None:
            return True
        assert(prefix_length > 0)
        # prefix is always on the left end of R2
        prefix_quality = self.r2.quality[offset:offset + prefix_length]
        for q in prefix_quality:
            if ord(q) < minimum_quality_score + ord('!'):
                return False
        return True

    def check_mutation_quality(self, minimum_quality_score):
        if minimum_quality_score is None  or  not self.mutations or (not self.r1.quality and not self.r2.quality):
            return 0
        removed = []
        r1_start = self.r1.match_index
        if self.r1.indels:
            r1rcseq, r1rqual = self.r1.apply_indels()
        else:
            r1rcseq, r1rqual = self.r1.reverse_complement, self.r1.reverse_quality
        if self.r2.indels:
            r2_start = self.r2.match_index + 1
            r2_end = r2_start + self.r2.match_len
            # apply_indels() uses subsequence, so no need to correct for ltrim above
            r2seq, r2qual = self.r2.apply_indels() 
        else:
            r2_start = self.r2.match_index - self.r2.ltrim + 1
            r2_end = r2_start + self.r2.original_len
            r2seq, r2qual = self.r2.original_seq, self.r2.quality
        min_quality = minimum_quality_score + ord('!')
        for mut in self.mutations:
            q = q1 = q2 = None
            nt1 = nt2 = None
            if mut < r2_end  and  mut >= r2_start:
                idx = mut - r2_start
                q2 = ord(r2qual[idx])
                nt2 = r2seq[idx]
            if mut > r1_start:
                idx = mut - r1_start - 1
                if idx < len(r1rqual):
                    q1 = ord(r1rqual[idx])
                    nt1 = r1rcseq[idx]
            if q1 and q2:
                # check for agreement
                if nt1 == nt2:
                    q = max(q1, q2)
                else:
                    # xref https://trello.com/c/usT0vTiG/308-discordant-mutations-case-handling-and-unit-tests
                    #    and associated test cases
                    ref = self.target.seq[mut - 1]
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
            if not q or q < min_quality:
                removed.append(mut)
        for mut in removed:
            self.mutations.remove(mut)
        if removed:
            self.removed_mutations = removed
        return len(removed)


    def check_overlap(self, ignore_indels = False):
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
                if ignore_indels or (not self.r1.indels and not self.r2.indels):
                    return False
                newr1,_ = self.r1.apply_indels()
                newr2,_ = self.r2.apply_indels()
                return newr1[r1start:r1start+overlap_len] == newr2[r2_match_len-overlap_len:r2_match_len]
        return True

    @property
    def target(self):
        return self._target

    @target.setter
    def target(self, val):
        self._target = val
        if val:
            self.r1.target_len = val.n
            self.r2.target_len = val.n
        else:
            self.r1.target_len = None
            self.r2.target_len = None

    @property
    def matched(self):
        return (self.target and self.r1.match_len and self.r2.match_len)

    @property
    def fully_matched(self):
        if not self._fully_matched and self.r1.fully_matched and self.r2.fully_matched:
            self._fully_matched = True
        return self._fully_matched

    @fully_matched.setter
    def fully_matched(self, val):
        self._fully_matched = val
        self.r1.trim_to_match()
        self.r2.trim_to_match()

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

    @property
    def indels_match(self):
        if self._indels_match is None:
            if not self.r1.indels and not self.r2.indels:
                self._indels_match = True
            elif len(self.r1.indels) != len(self.r2.indels):
                self._indels_match = False
            else:
                r1indels = dict(zip(self.r1.indels.keys(), map(vars, self.r1.indels.values())))
                r2indels = dict(zip(self.r2.indels.keys(), map(vars, self.r2.indels.values())))
                self._indels_match = (r1indels == r2indels)
        return self._indels_match

