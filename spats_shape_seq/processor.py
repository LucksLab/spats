
from counters import Counters
from mask import match_mask_optimized
from util import _warn, _debug, reverse_complement, string_match_errors

class Failures(object):
    mask = "mask failure"
    indeterminate = "indeterminate sequence failure"
    nomatch = "no match"
    adapter_trim = "adapter trim failure"
    mismatch = "R1/R2 mismatch"
    multiple_R1 = "multiple R1 match"
    linker = "linker match failure"
    left_of_zero = "R2 to left of site 0 failure"
    match_errors = "match errors failure"
    right_edge = "R1 right edge failure"
    cotrans_min = "cotrans minimum"
    r1_r2_overlap = "R1/R2 disagree on overlap"
    extra_r1 = "extra bp in R1"

    @staticmethod
    def all_failures():
        failures = []
        for attr in dir(Failures):
            if attr == "all_failures" or attr.startswith('_'):
                continue
            failures.append(getattr(Failures, attr))
        return failures


class PairProcessor(object):

    def __init__(self, run, targets, masks):
        self._run = run
        self._targets = targets
        self._masks = masks
        self.__adapter_t_rc = 0
        self.uses_tags = False
        self.counters = Counters(self._run)
        self._match_mask = self._match_mask_optimized if (run.masks[0] == 'RRRY' and run.masks[1] == 'YYYR') else self._match_mask_general
        if self._match_mask != self._match_mask_optimized:
            print("Warning: not using optimized mask match.")
        run.apply_config_restrictions()
        self.prepare()

    def exists(self):
        return True

    @property
    def _adapter_t_rc(self):
        if not self.__adapter_t_rc:
            self.__adapter_t_rc = reverse_complement(self._run.adapter_t)
        return self.__adapter_t_rc

    def _match_mask_general(self, pair):
        seq = pair.r1.original_seq
        for mask in self._masks:
            if mask.matches(seq):
                pair.set_mask(mask)
                self.counters.increment_mask(mask, pair.multiplicity)
                return True
        self.counters.mask_failure += pair.multiplicity
        pair.failure = Failures.mask
        return False

    # optimized version for RRRY/YYYR
    def _match_mask_optimized(self, pair):
        mask = match_mask_optimized(pair.r1.original_seq[:4], self._masks)
        if mask:
            pair.set_mask(mask)
            self.counters.increment_mask(mask, pair.multiplicity)
            return True
        else:
            self.counters.mask_failure += pair.multiplicity
            pair.failure = Failures.mask
            return False

    def _check_indeterminate(self, pair):
        if not self._run.allow_indeterminate  and  not pair.is_determinate():
            pair.failure = Failures.indeterminate
            self.counters.indeterminate += pair.multiplicity
            return False
        return True

    def prepare(self):
        raise Exception("subclasses must override")

    def process_pair(self, pair):
        raise Exception("subclasses must override")
