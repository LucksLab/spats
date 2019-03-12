
from processor import PairProcessor, Failures
from mask import base_similarity_ind
from util import _warn, _debug, reverse_complement, string_match_errors, string_find_errors, string_find_with_overlap, AlignmentParams


class PartialFindProcessor(PairProcessor):

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


    def _recheck_targets(self, pair):
        if pair.r1.rtrim and len(self._targets.targets) > 1:
            # since we trimmed down our R1 due to an adapter (or dumbbell) at some point,
            # need to see if that means the leftover matches multiple targets;
            # if so, need to reject this pair.
            target, _, _, _ = self._targets.find_partial(pair.r1.reverse_complement, min_length_override = pair.r1.match_len)
            if not target or isinstance(target, list):
                # note that target may be None if pair.r1.match_len is less than the index word length
                _debug("dropping pair due to multiple R1 match after adapter trim")
                pair.target = None
                pair.failure = Failures.multiple_R1
                self.counters.multiple_R1_match += pair.multiplicity
                return False
        return True


    def _cotrans_find_short_matches(self, pair):
        run = self._run
        r2li = pair.r2.subsequence.find(run.cotrans_linker)
        if r2li <= 0:
            self.counters.unmatched += pair.multiplicity
            pair.failure = Failures.nomatch
            return False
        target = self._targets.targets[0]
        r2spots = string_find_errors(pair.r2.subsequence[:r2li], target.seq, run.allowed_target_errors)
        r1spots = string_find_errors(reverse_complement(pair.r1.subsequence[:r2li]), target.seq, run.allowed_target_errors)
        if len(r1spots) == 0 or len(r2spots) == 0 or len(r1spots) != len(r2spots):
            self.counters.unmatched += pair.multiplicity
            pair.failure = Failures.nomatch
            return False
        elif len(r1spots) > 1 or len(r2spots) > 1:
            pair.failure = Failures.multiple_R1
            return False
        pair.target = target
        pair.r1.match_index = r1spots[0]
        pair.r1.match_start = pair.r1.seq_len - r2li
        pair.r1.match_len = r2li
        pair.r2.match_index = r2spots[0]
        pair.r2.match_start = 0
        pair.r2.match_len = r2li
        pair.fully_matched = True    # also sets the trim
        pair.linker = pair.r1.match_index + pair.r1.match_len
        return True


    ## Indels complicate trimming (of R2) greatly!
    ## For exmaple, prior to trimming linker, we might have:
    ##    RR*RRRL*LLLll
    ##    TTTTTTTTTttttttt   -> match_len = 9, indels_delta=-2
    ## After trimming 6 linker bps, we should get:
    ##    RR*RRR
    ##    TTTTTTtttttttttt   -> match_len = 6, indels_delta=-1
    ## Note that indels_delta for all of r2 is insufficient to compute the delta to match_len!
    ## So we have to chop and then trim indels.
    def _trim_adapters(self, pair, r1too):
        run = self._run
        r2_adapter_trimmed = None
        masklen = pair.mask.length()
        if run.cotrans:
            r2_linker_start = string_find_with_overlap(run.cotrans_linker, pair.r2.subsequence)
            ## Since we don't know r1.indels_delta yet, use r2's as an estimate,
            ## but be conservative and only take into account deletes.
            ## (It's ok if we trim a bit too much of R2.)
            if -1 != r2_linker_start  and  pair.r2.match_index + r2_linker_start + min(0, pair.r2.indels_delta) >= pair.r1.right_est:
                full_trim = pair.r2.seq_len - r2_linker_start
                r2_adapter_len = full_trim - len(run.cotrans_linker) - masklen
                if r2_adapter_len > 0:
                    r2_adapter_trimmed = pair.r2.subsequence[-r2_adapter_len:]
                pair.r2.rtrim += full_trim    # also updates indels and match_len in r2
            # TODO:  don't let r2 go beyond (trimmed) r1 end in cotrans case
        else:
            # trim everything beyond the end of the target (including mask if there)
            full_trim = pair.r2.right_est - pair.target.n
            if full_trim > 0:
                if full_trim > masklen:
                    r2_adapter_trimmed = pair.r2.subsequence[masklen - full_trim:]
                pair.r2.rtrim += full_trim    # also updates indels and match_len in r2
        if r2_adapter_trimmed:
            if run.minimum_adapter_len  and  len(r2_adapter_trimmed) < run.minimum_adapter_len:
                return False
            pair.r2.adapter_errors = string_match_errors(r2_adapter_trimmed, reverse_complement(run.adapter_t))
            if len(pair.r2.adapter_errors) > run.allowed_adapter_errors:
                return False
            self.counters.adapter_trimmed += pair.multiplicity
        if not r1too:
            return True
        ## Note: we really shouldn't trim r1 prior to aligning b/c it could have big inserts, but
        ## trimming greatly saves alignment time, so compromise by assuming that indels (mostly) match
        ## between reads so use R2's indels_delta to buffer trim amount.
        if run.cotrans:
            possible_matchlen = r2_linker_start if -1 != r2_linker_start else pair.r2.match_len
        else:
            possible_matchlen = min(pair.target.n - max(0, pair.r2.match_index - pair.r2.match_start), pair.r1.seq_len + pair.r1.match_index - pair.r1.match_start)
        r1_adapter_trim = pair.r1.seq_len - possible_matchlen - max(pair.r2.indels_delta, 0)
        dumblen = len(run.dumbbell) if run.dumbbell else 0
        if run.minimum_adapter_len  and  r1_adapter_trim - dumblen < run.minimum_adapter_len:
            return False
        if r1_adapter_trim > 0:
            r1_adapter = pair.r1.subsequence[dumblen - r1_adapter_trim:]
            pair.r1.adapter_errors = string_match_errors(r1_adapter, run.adapter_b)
            if len(pair.r1.adapter_errors) > run.allowed_adapter_errors:
                return False
            pair.r1.rtrim += r1_adapter_trim
            pair.r1.match_start -= r1_adapter_trim
            if pair.r1.match_start < 0:
                pair.r1.match_index -= pair.r1.match_start
                pair.r1.match_len += pair.r1.match_start
                pair.r1.match_start = 0
            #self.counters.adapter_trimmed += pair.multiplicity   # Note:  only counting R2 trimming 
        return self._recheck_targets(pair)


    def _extend_match(self, pair, r1_fully_rtrimmed):
        run = self._run
        if run.allow_indeterminate:
            simfn = lambda nt1, nt2: base_similarity_ind(nt1, nt2, run.indel_match_value, run.indel_mismatch_cost, .5 * run.indel_match_value)
        else:
            simfn = lambda nt1, nt2: AlignmentParams.char_sim(nt1, nt2, run.indel_match_value, run.indel_mismatch_cost) 
        ap = AlignmentParams(simfn, run.indel_gap_open_cost, run.indel_gap_extend_cost)

        masklen = pair.mask.length()
        dumblen = len(run.dumbbell) if run.dumbbell else 0

        ## Align remaining R2 first to find adapter overhang.
        ## TAI: probably better way of doing adapter alignment  (see indels_7 testcase)
        r2suffix = "" if run.cotrans else reverse_complement(pair.r1.original_seq[:masklen]) + reverse_complement(run.adapter_t)
        pair.r2.align_with_target(pair.target, ap, r2suffix)
        if run.dumbbell:
            pair.dumbbell = pair.r2.match_index - dumblen

        ## Trim the adapters off both R1 and R2, and remaining linker off R2
        if not self._trim_adapters(pair, not r1_fully_rtrimmed):
            if not pair.failure:
                pair.failure = Failures.adapter_trim
                self.counters.adapter_trim_failure += pair.multiplicity
            return False

        ## Now align remaining R1 if necessary
        if not pair.r1.fully_matched:
            pair.r1.align_with_target(pair.target, ap)
            # we may have not trimmed enough
            r1_overhang = pair.r1.right_est - pair.target.n
            if r1_overhang > 0:
                pair.r1.ltrim += r1_overhang    # also updates match_len
            if pair.linker:
                pair.r1.match_len = pair.linker - pair.r1.match_start - pair.r1.indels_delta
                pair.linker += pair.r1.match_index - pair.r1.rtrim - pair.r1.indels_delta  # linker now equates to 1 beyond the end of the match in the *target* (where the linker starts)  TAI:  keep in R1 coordinates?

        return True


    def _verify_full_match(self, pair):
        run = self._run
        maxmatch =  min(pair.target.n, pair.r1.match_index + pair.r1.match_len)
        r2_adapter_trimmed = ""
        if pair.r2.match_index + pair.r2.seq_len > maxmatch:
            masklen = pair.mask.length()
            trimmed = pair.r2.right_est - maxmatch
            if run.cotrans:
                linker_len = len(run.cotrans_linker)
                suffix = pair.r2.subsequence[-trimmed:]
                if ((trimmed <  linker_len and not run.cotrans_linker.startswith(suffix))
                    or (trimmed > linker_len and not suffix.startswith(run.cotrans_linker))):
                    pair.failure = Failures.linker
                    return False
                if trimmed > linker_len + masklen:
                    r2_adapter_trimmed = pair.r2.subsequence[masklen + linker_len - trimmed:]
            elif trimmed > masklen:
                r2_adapter_trimmed = pair.r2.subsequence[masklen - trimmed:]
            if trimmed > 0:
                pair.r2.rtrim = trimmed
            if r2_adapter_trimmed:
                if run.minimum_adapter_len  and  len(r2_adapter_trimmed) < run.minimum_adapter_len:
                    pair.failure = Failures.adapter_trim
                    self.counters.adapter_trim_failure += pair.multiplicity
                    return False
                pair.r2.adapter_errors = string_match_errors(r2_adapter_trimmed, reverse_complement(run.adapter_t))
                if len(pair.r2.adapter_errors) > run.allowed_adapter_errors:
                    pair.failure = Failures.adapter_trim
                    self.counters.adapter_trim_failure += pair.multiplicity
                    return False
                self.counters.adapter_trimmed += pair.multiplicity
        if pair.r1.match_index == 0  and  pair.r1.match_start > 0  and  not run.dumbbell:
            # Note: if we had a full dumbbell, it would have been entirely trimmed already
            # So we either have no dumbbell and an adapter or a partial dumbbell only.
            r1_adapter_trimmed = pair.r1.reverse_complement[-pair.r1.match_start:]
            pair.r1.rtrim += pair.r1.match_start
            pair.r1.match_start = 0
            if r1_adapter_trimmed:
                if run.minimum_adapter_len  and  len(r1_adapter_trim) < run.minimum_adapter_len:
                    pair.failure = Failures.adapter_trim
                    self.counters.adapter_trim_failure += pair.multiplicity
                    return False
                pair.r1.adapter_errors = string_match_errors(r1_adapter_trimmed, run.adapter_b)
                if len(pair.r1.adapter_errors) > run.allowed_adapter_errors:
                    pair.failure = Failures.adapter_trim
                    self.counters.adapter_trim_failure += pair.multiplicity
                    return False
            #self.counters.adapter_trimmed += pair.multiplicity    # Note: only counting R2 trimming
        return self._recheck_targets(pair)


    ## Generally, experiments look like:
    ##    R1':    AAAAAAAAAAAAAA.DDDDDDDD.TTTTTTTTTTTTTTTTT.LLLLLLL.MMMM...................
    ##    R2:     ...............DDDDDDDD.TTTTTTTTTTTTTTTTT.LLLLLLL.MMMM.AAAAAAAAAAAAAAAAAA
    ## Where:
    ##    A+ = adapters (adapter_b on left and adapter_t on right)
    ##    D+ = optional dumbbell (forward/R2 reads will always start here if present)
    ##    T+ = region to be matched/aligned with a reference target;
    ##         if no cotrans/linker, then always goes all the way to the right end
    ##         of the reference target, but does not necessarily start at the beginning.
    ##    L+ = optionl linker in the case of cotrans experiments (R2 read may end before or during this)
    ##    M+ = mask/handle designating treatment (reverse/R1 reads will always start here)

    def process_pair(self, pair):
        if not self._check_indeterminate(pair) or not self._match_mask(pair):
            return

        run = self._run
        masklen = pair.mask.length()
        pair.r2.auto_adjust_match = True
        #pair.r1.auto_adjust_match = True

        ## Pre-trim dumbbell from R2 immediately since it won't help with alignment
        ## and simplifies dealing with match_index, match_start, etc.
        ## (It's also a bit faster for SW.)
        r1_fully_rtrimmed = False
        dumblen = 0
        if run.dumbbell:
            if not pair.r2.original_seq.startswith(run.dumbbell):
                pair.failure = Failures.dumbbell
                self.counters.dumbbell_failures += pair.multiplicity
                return
            dumblen = len(run.dumbbell)
            pair.r2.ltrim = dumblen
            pair.dumbbell = -dumblen
            # Also may as well trim R1 dumbbell if it's there in full...
            # TAI:  this may get the adapter too. if so, do we need to count that?  (no. only counting R2s.)
            dbi = pair.r1.original_seq.rfind(reverse_complement(run.dumbbell), masklen + run.minimum_target_match_length)
            if -1 != dbi:
                pair.r1.rtrim = pair.r1.original_len - dbi
                r1_fully_rtrimmed = True

        ## Also pre-trim any (full) adapter on R1 for the same reason
        if not r1_fully_rtrimmed  and  run.adapter_b:
            adi = pair.r1.original_seq.rfind(run.adapter_b, masklen + run.minimum_target_match_length)
            if -1 != adi:
                pair.r1.rtrim = pair.r1.original_len - adi
                r1_fully_rtrimmed = True
                self.counters.adapter_trimmed += pair.multiplicity

        ## Finally, pre-trim any linker from R1
        if run.cotrans:
            if pair.r1.subsequence.startswith(reverse_complement(run.cotrans_linker)):
                pair.r1.ltrim += len(run.cotrans_linker)
                pair.linker = pair.r1.seq_len   # linker is currently in R1 coordinates
            else:
                pair.failure = Failures.linker
                return

        ## Now find best exactly matching substring with target
        self._find_matches(pair)
        if not pair.matched:
            if not run.cotrans or pair.r1.ltrim == 0:
                self.counters.unmatched += pair.multiplicity
                pair.failure = Failures.nomatch
                return
            if not self._cotrans_find_short_matches(pair):
                return

        ## And extend the match if necessary using string alignment on the rest
        aligned = False
        if pair.fully_matched:
            # Don't just trim adapter here b/c in cotrans case it has already been done
            if not self._verify_full_match(pair):
                return
        elif run.handle_indels:
            if not self._extend_match(pair, r1_fully_rtrimmed):
                return
            aligned = True
        elif not self._trim_adapters(pair, not r1_fully_rtrimmed):
            if not pair.failure:
                pair.failure = Failures.adapter_trim
                self.counters.adapter_trim_failure += pair.multiplicity
            return

        counted_prefix = None
        if pair.r2.match_start > pair.r2.match_index:
            assert(not aligned or pair.r2.match_index == 0)  # align_strings() will ensure this if penalize_ends is True
            r2_start_in_target = pair.r2.match_index - pair.r2.match_start    # will be negative
            self.counters.left_of_target += pair.multiplicity
            prefix = pair.r2.original_seq[dumblen:dumblen - r2_start_in_target]
            if run.count_left_prefixes:
                if (run.mutations_require_quality_score is None) or pair.check_prefix_quality(dumblen - r2_start_in_target, run.mutations_require_quality_score, dumblen):
                    self.counters.register_prefix(prefix, pair)
                    counted_prefix = prefix
                else:
                    self.counters.low_quality_prefixes += pair.multiplicity
            if run.collapse_left_prefixes and (not run.collapse_only_prefixes or prefix in run._p_collapse_only_prefix_list):
                pair.r2.ltrim -= r2_start_in_target
                pair.r2.match_start += r2_start_in_target
                pair.r2.shift_indels(r2_start_in_target)
                if pair.r1.rtrim and pair.r1.match_index == pair.r2.match_index:
                    pair.r1.rtrim -= r2_start_in_target
                    pair.r1.match_start += r2_start_in_target
                    pair.r1.shift_indels(r2_start_in_target)
            else:
                pair.failure = Failures.left_of_zero
                self.counters.left_of_zero += pair.multiplicity
                return
        elif run.count_left_prefixes:
            counted_prefix = "NONE"
            self.counters.register_prefix(counted_prefix, pair)

        if run.cotrans and pair.right < run.cotrans_minimum_length:
            pair.failure = Failures.cotrans_min
            return

        if not aligned:
            pair.r1.match_to_seq()
            pair.r2.match_to_seq()
            pair.r1.match_errors = string_match_errors(pair.r1.reverse_complement, pair.target.seq[pair.r1.match_index:])
            pair.r2.match_errors = string_match_errors(pair.r2.subsequence, pair.target.seq[pair.r2.match_index:])

        if run._p_rois  and  (pair.r1.error_in_region(run._p_rois) or pair.r2.error_in_region(run._p_rois)):
            pair.interesting = True
            self.counters.interesting_pairs += pair.multiplicity

        if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > run.allowed_target_errors:
            if pair.r1.match_errors:
                _debug("R1 errors: {}".format(pair.r1.match_errors))
            if pair.r2.match_errors:
                _debug("R2 errors: {}".format(pair.r2.match_errors))
            pair.failure = Failures.match_errors
            self.counters.match_errors += pair.multiplicity
            return

        if ((not run.cotrans and pair.right != pair.target.n and not run.allow_multiple_rt_starts)  or
            (run.cotrans  and  (pair.right > pair.target.n or pair.right < pair.left or
             pair.r1.rtrim + pair.r1.match_start + pair.r1.match_len + pair.r1.indels_delta + pair.r1.ltrim < pair.r1.original_len))):
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

        if run.handle_indels and not pair.indels_match:
            self.counters.mismatching_indel_pairs += pair.multiplicity
            # TAI:  might bail on these

        pair.site = pair.left
        if run._p_rois and not pair.interesting:
            for roi in run._p_rois:
                if roi[0] <= pair.site  and  pair.site <= roi[1]:
                    pair.interesting = True
                    self.counters.interesting_pairs += pair.multiplicity
                    break
        self.counters.register_count(pair)
        if counted_prefix:
            self.counters.register_mapped_prefix(counted_prefix, pair)
        if run.count_mutations:
            self.counters.register_mapped_mut_count(pair)

