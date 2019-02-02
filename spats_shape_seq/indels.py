
from processor import PairProcessor, Failures
from mask import base_similarity_ind
from util import _warn, _debug, reverse_complement, char_sim


class IndelsProcessor(PairProcessor):
    ''' This will be slow. '''

    def prepare(self):
        self._targets.index()

    def _find_matches(self, pair):
        # use R1 to determine which target
        target = pair.r1.find_in_targets(self._targets)
        if isinstance(target, list):
            _debug("dropping pair due to multiple R1 match")
            pair.failure = Failures.multiple_R1
            self.counters.multiple_R1_match += pair.multiplicity
            return
        if target:
            pair.target = pair.r2.find_in_targets(self._targets, force_target = target)
        _debug([pair.r1.match_start, pair.r1.match_len, pair.r1.match_index, "--", pair.r2.match_start, pair.r2.match_len, pair.r2.match_index])


    def process_pair(self, pair):
        if not self._check_indeterminate(pair) or not self._match_mask(pair):
            return

        self._find_matches(pair)
        if not pair.matched:
            self.counters.unmatched += pair.multiplicity
            pair.failure = Failures.nomatch
            return

        run = self._run

        dumblen = 0
        if run.dumbbell:
            if not pair.r2.original_seq.startswith(run.dumbbell):
                pair.failure = Failures.dumbbell
                return
            dumblen = len(run.dumbbell)
            pair.r2.ltrim += dumblen
            dbi = pair.r1.original_seq.find(reverse_complement(run.dumbbell))
            if -1 == dbi:
                pair.failure = Failures.dumbbell
                return
            pair.r1.rtrim += pair.r1.original_len - dbi

        if run.cotrans:
            if run.dumbbell:
                raise Exception("can't use Indel algorithm with both cotrans and a dumbbell")
            pair.r1.ltrim += len(run.cotrans_linker)
            pair.r2.rtrim += len(run.cotrans_linker)

        pair.r1.match_start -= pair.r1.rtrim
        pair.r2.match_start -= pair.r2.ltrim

        masklen = pair.mask.length()
        r2suffix = reverse_complement(pair.r1.original_seq[:masklen]) + reverse_complement(run.adapter_t) 
        if run.allow_indeterminate:
            simfn = lambda nt1, nt2: base_similarity_ind(nt1, nt2, run.indel_match_value, run.indel_mismatch_cost, .5 * run.indel_match_value)
        else:
            simfn = lambda nt1, nt2: char_sim(nt1, nt2, run.indel_match_value, run.indel_mismatch_cost) 

        pair.r2.align_with_target(pair.target, simfn, run.indel_gap_open_cost, run.indel_gap_extend_cost, r2suffix)
        r2_adapter_trim = max(0, pair.r2.match_index + pair.r2.match_len - pair.target.n)
        r1_adapter_trim = pair.r1.seq_len - (pair.target.n - pair.r2.match_index)
        if r1_adapter_trim > 0 and not dumblen:
            pair.r1.rtrim += r1_adapter_trim
            pair.r1.match_start -= r1_adapter_trim
        pair.r1.align_with_target(pair.target, simfn, run.indel_gap_open_cost, run.indel_gap_extend_cost)

        if run.minimum_adapter_len and r1_adapter_trim < run.minimum_adapter_len:
            _debug("  !! v102 minimum adapter len {}".format(r1_adapter_trim))
            pair.failure = Failures.adapter_trim
            self.counters.adapter_trim_failure += pair.multiplicity
            return 

        pair.r1.shift_errors(pair.target.n)   # make errors relative to match_index, not 0

        pair.r2.rtrim += r2_adapter_trim
        pair.r2.match_len -= r2_adapter_trim
        adapter_indels = pair.r2.trim_indels(pair.target.n)
        pair.r2.shift_errors(pair.target.n, masklen)
        if run.dumbbell:
            pair.dumbbell = pair.r2.match_index - dumblen

        if adapter_indels + len(pair.r2.adapter_errors) > run.allowed_adapter_errors:
            pair.failure = Failures.adapter_trim
            self.counters.adapter_trim_failure += pair.multiplicity
            return

        if pair.r2.match_start > pair.r2.match_index:
            assert(pair.r2.match_index == 0)    # align_strings() will ensure this if penalize_ends is True
            r2_start_in_target = pair.r2.match_index - pair.r2.match_start    # will be negative
            self.counters.left_of_target += pair.multiplicity
            prefix = pair.r2.original_seq[dumblen:dumblen - r2_start_in_target]
            if run.count_left_prefixes:
                if (run.mutations_require_quality_score is None) or pair.check_prefix_quality(0 - r2_start_in_target, run.mutations_require_quality_score):
                    self.counters.register_prefix(prefix, pair)
                else:
                    self.counters.low_quality_prefixes += pair.multiplicity
            if run.collapse_left_prefixes and (not run.collapse_only_prefixes or prefix in run._p_collapse_only_prefix_list):
                pair.r2.ltrim -= r2_start_in_target
                pair.r2.match_start += r2_start_in_target
            else:
                pair.failure = Failures.left_of_zero
                return

        if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > run.allowed_target_errors:
            if pair.r1.match_errors:
                _debug("R1 errors: {}".format(pair.r1.match_errors))
            if pair.r2.match_errors:
                _debug("R2 errors: {}".format(pair.r2.match_errors))
            if run._p_v102_compat and not (pair.r1.match_errors + [e for e in pair.r2.match_errors if e < pair.r2.original_len - masklen]):
                _debug("** v102 compat, allowing")
            else:
                pair.failure = Failures.match_errors
                self.counters.match_errors += pair.multiplicity
                return

        if pair.right != pair.target.n and not run.allow_multiple_rt_starts:
            pair.failure = Failures.right_edge
            self.counters.r1_not_on_right_edge += pair.multiplicity
            return

        if run.ignore_stops_with_mismatched_overlap  and  not pair.check_overlap(True):
            pair.failure = Failures.r1_r2_overlap
            self.counters.r1_r2_overlap += pair.multiplicity
            return

        if run.count_mutations:
            pair.check_mutations()
            self.counters.register_mut_count(pair)
            if pair.mutations and len(pair.mutations) > run.allowed_target_errors:
                pair.failure = Failures.match_errors
                self.counters.match_errors += pair.multiplicity
                return
            self.counters.low_quality_muts += pair.check_mutation_quality(run.mutations_require_quality_score)

        if run.count_only_full_reads and pair.left != 0:
            pair.failure = Failures.not_full_read
            return

        if run.count_indels:
            self.counters.r1_indels += (pair.multiplicity * len(pair.r1.indels))
            self.counters.r2_indels += (pair.multiplicity * len(pair.r2.indels))
            if not pair.indels_match:
                self.counters.mismatching_indel_pairs += pair.multiplicity

        pair.site = pair.left
        self.counters.register_count(pair)

