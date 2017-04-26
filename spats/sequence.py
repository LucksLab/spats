
from config import spats_config
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
    def seq_len(self):
        return self._length - self._ltrim - self._rtrim

    @property
    def reverse_complement(self):
        return reverse_complement(self.subsequence)

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

    def find_in_target(self, target, reverse_complement = False):
        seq = self.reverse_complement if reverse_complement else self.subsequence
        self.match_start, self.match_len, self.match_index = target.find_partial(seq, 2 * spats_config.minimum_target_match_length)
        if not self.match_len:
            # it's much faster to search for longer partial matches, then fall back on the minimum
            self.match_start, self.match_len, self.match_index = target.find_partial(seq, spats_config.minimum_target_match_length)

    def trim(self, length, reverse_complement = False):
        self._rtrim = length
        if reverse_complement:
            delta = self.match_len - self.seq_len
            if delta > 0:
                _debug("trim reducing original match_len {} -> {}".format(self.match_len, self.seq_len))
                self.match_len = self.seq_len
                self.match_start += delta
                self.match_index += delta

    def match_to_seq(self, reverse_complement = False):
        if reverse_complement:
            self.match_index -= (self._rtrim - self.match_start)
            self.match_start = self._rtrim
            self.match_len = self.seq_len
        else:
            self.match_index -= self.match_start
            self.match_start = 0
            self.match_len = self.seq_len
