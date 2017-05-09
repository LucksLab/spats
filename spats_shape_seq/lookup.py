
from processor import PairProcessor
from util import _warn, _debug, reverse_complement


class LookupProcessor(PairProcessor):

    def prepare(self):
        self._targets.index()
        self._targets.build_lookups(self._run.adapter_b)

    #@profile
    def process_pair(self, pair):

        if not self._check_indeterminate(pair) or not self._match_mask(pair):
            return

        if not self._run.allow_indeterminate  and  not pair.is_determinate():
            pair.failure = "indeterminate sequence failure"
            self.counters.indeterminate += pair.multiplicity
            return

        targets = self._targets
        r1_res = targets.r1_lookup.get(pair.r1.original_seq[4:])
        if not r1_res or not r1_res[0]:
            pair.failure = "no match"
            return

        site = -1
        target = r1_res[0]
        r2len = pair.r2.original_len
        if r1_res[1] == None:
            lookup = targets.r2_lookup[target.name]
            r2_match_len = len(lookup.keys()[0])
            r2_key = pair.r2.original_seq[:r2_match_len]
            r2_res = lookup.get(r2_key)
            if r2_res is not None:
                match_site = r2_res
            else:
                pair.failure = "no match"
                return
        else:
            match_site = r1_res[1]

        # need to check R2 against expectation
        match_len = min(r2len, target.n - match_site)
        if pair.r2.original_seq[:match_len] == target.seq[match_site:match_site + match_len]:
            adapter_len = r2len - match_len - 4
            if adapter_len <= 0 or self._adapter_t_rc[:adapter_len] == pair.r2.original_seq[-adapter_len:]:
                site = match_site
            else:
                pair.failure = "adapter trim failure"
                return
        else:
            pair.failure = "no match"
            return

        if site is None or site == -1:
            raise Exception("No site?")

        # in rare cases, need to double-check what R1 should be based on R2
        if site != r1_res[1]:
            match_len = min(pair.r1.original_len - 4, target.n - match_site)
            adapter_len = pair.r1.original_len - match_len - 4
            r1_match = reverse_complement(target.seq[-match_len:]) + self._run.adapter_b[:adapter_len]
            if pair.r1.original_seq[4:] != r1_match:
                pair.failure = "R1/R2 mismatch"
                return

        pair.target = target
        pair.site = site
        n = target.n

        assert(pair.site is not None)
        pair.register_count()
        self.counters.processed_pairs += pair.multiplicity
