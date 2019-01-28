
from util import _debug, reverse_complement

class Sequence(object):

    def __init__(self):
        self._reset(None)

    def _reset(self, seq):
        self._seq = seq.upper() if seq else None
        self._length = len(seq) if seq else None
        self.match_start = None
        self.match_len = None
        self.match_index = None
        self._ltrim = 0
        self._rtrim = 0
        self.adapter_errors = []
        self.match_errors = []
        self.indels = {}
        self.indels_delta = 0
        self.quality = None

    def set_seq(self, seq):
        self._reset(seq)

    @property
    def original_seq(self):
        return self._seq

    @property
    def original_len(self):
        return self._length

    @property
    def subsequence(self):
        if self._rtrim:
            return self._seq[self._ltrim:-self._rtrim]
        else:
            return self._seq[self._ltrim:]

    @property
    def subquality(self):
        if self._rtrim:
            return self.quality[self._ltrim:-self._rtrim]
        else:
            return self.quality[self._ltrim:]

    @property
    def seq_len(self):
        return self._length - self._ltrim - self._rtrim

    @property
    def reverse_complement(self):
        return reverse_complement(self.subsequence)

    @property
    def reverse_complement_quality(self):
        return reverse_complement(self.subquality)

    @property
    def matched(self):
        return bool(self.match_len)

    @property
    def left(self):
        return self.match_index if self.match_len else None

    @property
    def right(self):
        return self.match_index + self.match_len if self.match_len else None

    @property
    def trimmed(self):
        return (0 != self._rtrim)

    def find_in_targets(self, targets, reverse_complement = False, force_target = None, min_length_override = 0):
        seq = self.reverse_complement if reverse_complement else self.subsequence
        target, self.match_start, self.match_len, self.match_index = targets.find_partial(seq, force_target = force_target, min_length_override = min_length_override)
        return target

    def trim(self, length, reverse_complement = False):
        self._rtrim = length
        if reverse_complement:
            delta = self._rtrim - self.match_start
            _debug("trim reducing original match_len {} -> {}".format(self.match_len, self.match_len - delta))
            self.match_len -= delta
            self.match_start += delta
            self.match_index += delta
            return True
        return False

    def match_to_seq(self, reverse_complement = False):
        if reverse_complement:
            if 0 == self._rtrim:
                # we've already done this if we've trimmed
                self.match_index -= self.match_start
                self.match_start = self._rtrim
            self.match_len = self.seq_len
        else:
            self.match_index -= (self.match_start - self._ltrim)
            self.match_start = self._ltrim
            self.match_len = self.seq_len
        _debug(["M2S:", self.match_index, self.match_len, self.match_start, "-- ", self._rtrim])

    def matched_alignment(self, alignment, target, prefix = 0, suffix = 0):
        self.match_index = alignment.target_match_start
        self.match_len = alignment.target_match_len
        self.match_start = alignment.src_match_start
        self.indels = alignment.indels
        self.match_errors = alignment.mismatched
        self.indels_delta = alignment.indels_delta
        if prefix:
            ## Hack for adapter_b
            if self.match_index < prefix:
                delta = prefix - self.match_index
                self.match_len -= delta
                self.match_start += delta
                self.match_index += delta
                self._rtrim = self.match_start
            self.match_index -= prefix
            new_indels = {}
            for indind in self.indels.keys():
                newind = indind - prefix
                if newind > 0:
                    new_indels[newind] = self.indels[indind]
                    indlen = len(new_indels[newind].seq)
                    if not new_indels[newind].insert_type and newind - indlen < 0:
                        self.match_index += indlen
                        self.match_len -= indlen
            self.indels = new_indels
            for i in xrange(len(self.match_errors)):
                self.match_errors[i] -= prefix
        if suffix  and  self.match_index + self.match_len > target.n:
            ## Hack for adapter_t and mask
            self.trim(self.match_index + self.match_len - target.n)
            self.match_len = target.n - self.match_index
            self.trim_indels(target.n)

    def shift_errors(self, target_len):
        nme = []
        for err in self.match_errors:
            if err < target_len:
                nme.append(err - self.match_index) 
            else:
                self.adapter_errors.append(err - target_len)
        self.match_errors = nme

    def trim_indels(self, maxind):
        ## TAI: trim match_errors too?
        for indind in self.indels.keys():
            if indind >= maxind:
                del self.indels[indind]
