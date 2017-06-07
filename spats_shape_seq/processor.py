
from util import _warn, _debug, Counters, reverse_complement, string_match_errors


class PairProcessor(object):

    def __init__(self, run, targets, masks):
        self._run = run
        self._targets = targets
        self._masks = masks
        self.__adapter_t_rc = 0
        self.uses_tags = False
        self.counters = Counters()
        self.prepare()

    @property
    def _adapter_t_rc(self):
        if not self.__adapter_t_rc:
            self.__adapter_t_rc = reverse_complement(self._run.adapter_t)
        return self.__adapter_t_rc

    def _match_mask(self, pair):
        seq = pair.r1.original_seq
        for mask in self._masks:
            if mask.matches(seq):
                pair.set_mask(mask)
                return True
        self.counters.mask_failure += pair.multiplicity
        pair.failure = "mask failure"
        return False

    def _check_indeterminate(self, pair):
        if not self._run.allow_indeterminate  and  not pair.is_determinate():
            pair.failure = "indeterminate sequence failure"
            self.counters.indeterminate += pair.multiplicity
            return False
        return True

    def prepare(self):
        raise Exception("subclasses must override")

    def process_pair(self, pair):
        raise Exception("subclasses must override")
