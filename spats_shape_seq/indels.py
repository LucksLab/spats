
from processor import PairProcessor, Failures
from mask import base_similarity_ind
from util import _warn, _debug, reverse_complement, char_sim, align_strings


class IndelsProcessor(PairProcessor):
    ''' This will be slow. '''

    def prepare(self):
        pass

    def process_pair(self, pair):
        if not self._check_indeterminate(pair) or not self._match_mask(pair):
            return

        run = self._run

        dumblen = 0
        if run.dumbbell:
            if not pair.r2.original_seq.startswith(run.dumbbell):
                pair.failure = Failures.dumbbell
                return
            dumblen = len(run.dumbbell)
            pair.r2._ltrim += dumblen
            dbi = pair.r1.original_seq.find(reverse_complement(run.dumbbell))
            if -1 == dbi:
                pair.failure = Failures.dumbbell
                return
            pair.r1._rtrim = pair.r1.original_len - dbi

        if run.cotrans:
            if run.dumbbell:
                raise Exception("can't use Indel algorithm with both cotrans and a dumbbell")
            pair.r1._ltrim += len(run.cotrans_linker)
            pair.r2._rtrim += len(run.cotrans_linker)

        masklen = pair.mask.length()
        simfn = char_sim if not run.allow_indeterminate else lambda nt1, nt2: base_similarity_ind(nt1, nt2, run.indel_match_value, run.indel_mismatch_cost, .5 * run.indel_match_value)
        r2suffix =  reverse_complement(pair.r1.original_seq[:masklen]) + reverse_complement(run.adapter_t) 
        max_score = -1
        for target in self._targets.targets:
            align_r2 = align_strings(pair.r2.subsequence, target.seq + r2suffix, simfn, run.indel_gap_open_cost, run.indel_gap_extend_cost)
            r2_adapter_trim = max(0, align_r2.target_match_end - target.n + 1)
            r1_adapter_trim = max(0, pair.r1.seq_len - (target.n - align_r2.target_match_start)) if not dumblen else 0
            align_r1 = align_strings(pair.r1.reverse_complement[r1_adapter_trim:], target.seq, simfn, run.indel_gap_open_cost, run.indel_gap_extend_cost)
            score = align_r1.score + align_r2.score
            if score > max_score:
                max_score = score
                best_align_r1, best_align_r2 = align_r1, align_r2
                best_r1_adapter_trim, best_r2_adapter_trim = r1_adapter_trim, r2_adapter_trim
                pair.target = target

        if max_score <= 0  or  best_align_r1.target_match_len == 0  or  best_align_r2.target_match_len <= best_r2_adapter_trim  or  min(best_align_r1.max_run, best_align_r2.max_run) < run.minimum_target_match_length:
            pair.failure = Failures.nomatch
            self.counters.unmatched += pair.multiplicity
            return
        if run.minimum_adapter_len and best_r1_adapter_trim < run.minimum_adapter_len:
            _debug("  !! v102 minimum adapter len {}".format(best_r1_adapter_trim))
            pair.failure = Failures.adapter_trim
            self.counters.adapter_trim_failure += pair.multiplicity
            return 

        pair.r1.matched_alignment(best_align_r1, pair.target)
        pair.r1._rtrim += best_r1_adapter_trim
        pair.r1.shift_errors(pair.target.n)   # make errors relative to match_index, not 0

        pair.r2.matched_alignment(best_align_r2, pair.target)
        pair.r2._rtrim += best_r2_adapter_trim
        pair.r2.match_len -= best_r2_adapter_trim
        adapter_indels = pair.r2.trim_indels(pair.target.n)
        pair.r2.shift_errors(pair.target.n + masklen)   # note: assuming no errors in the mask on R2
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
                pair.r2._ltrim -= r2_start_in_target
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
            if not pair.indels_match():
                self.counters.mismatching_indel_pairs += pair.multiplicity

        pair.site = pair.left
        self.counters.register_count(pair)

