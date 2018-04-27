
from processor import PairProcessor, Failures
from util import _warn, _debug, string_match_errors, string_find_errors


class PartialFindProcessor(PairProcessor):

    def prepare(self):
        self._targets.index()

    def _find_matches(self, pair):
        # use R1 to determine which target
        target = pair.r1.find_in_targets(self._targets, reverse_complement = True)
        if isinstance(target, list):
            _debug("dropping pair due to multiple R1 match")
            pair.failure = Failures.multiple_R1
            self.counters.multiple_R1_match += pair.multiplicity
            return
        if target:
            pair.target = pair.r2.find_in_targets(self._targets, force_target = target)
        _debug([pair.r1.match_start, pair.r1.match_len, pair.r1.match_index, "--", pair.r2.match_start, pair.r2.match_len, pair.r2.match_index])

    def _trim_adapters(self, pair):

        # if we're here, we know that R2 hangs over the right edge. first, figure out by how much

        r2_seq = pair.r2.subsequence
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        r2_len_should_be = pair.target.n - r2_start_in_target
        r2_length_to_trim = pair.r2.seq_len - r2_len_should_be
        pair.r2.trim(r2_length_to_trim)
        _debug("R2 trim: {}, {}, {}".format(r2_start_in_target, r2_len_should_be, r2_length_to_trim))

        if self._run.minimum_adapter_len and r2_length_to_trim - 4 < self._run.minimum_adapter_len:
            _debug("  !! v102 minimum adapter len {}".format(r2_length_to_trim - 4))
            return False

        if r2_length_to_trim <= 4:
            # TODO: should we verify that this matches RC of R1 handle, and register errors for bp that don't?
            # for now, just ignore this part
            # also, this means that there's nothing to trim from R1, so we're done
            return True

        # find out how good of a match the end of R2 is for adapter_t_rc
        r2_adapter_match = r2_seq[4-r2_length_to_trim:]
        pair.r2.adapter_errors = string_match_errors(r2_adapter_match, self._adapter_t_rc)
        _debug("  check = {}, errors = {}".format(r2_adapter_match, pair.r2.adapter_errors))
        if len(pair.r2.adapter_errors) > self._run.allowed_adapter_errors:
            return False

        # now, same thing on r1 (smaller trim b/c of no handle, hence -4)
        r1_seq = pair.r1.subsequence
        r1_length_to_trim = r2_length_to_trim - 4
        r1_adapter_match = r1_seq[-r1_length_to_trim:]
        r1_match_trimmed = pair.r1.trim(r1_length_to_trim, reverse_complement = True)
        pair.r1.adapter_errors = string_match_errors(r1_adapter_match, self._run.adapter_b)
        _debug("  R1 check = {}, errors = {}".format(r1_adapter_match, pair.r1.adapter_errors))
        if len(pair.r1.adapter_errors) > self._run.allowed_adapter_errors:
             return False

        if r1_match_trimmed and len(self._targets.targets) > 1:
            # ok, we trimmed down our R1 due to adapters. need to see if that means the leftover matches
            # multiple targets; if so, need to reject this pair.
            target = pair.r1.find_in_targets(self._targets, reverse_complement = True, min_length_override = pair.r1.match_len)
            if not target or isinstance(target, list):
                # note that target may be None if pair.r1.match_len is less than the index word length
                _debug("dropping pair due to multiple R1 match after adapter trim")
                print("multiple R1 trim: " + pair.r1.original_seq)
                print("{}".format([r1_match_trimmed, r1_adapter_match, r1_length_to_trim]))
                pair.target = None
                pair.failure = Failures.multiple_R1
                self.counters.multiple_R1_match += pair.multiplicity
                return False

        _debug("successful adapter trim of {}/{} bp from R1/R2".format(pair.r1._rtrim, pair.r2._rtrim))

        return True

    def process_pair(self, pair):

        if not self._check_indeterminate(pair) or not self._match_mask(pair):
            return

        self._find_matches(pair)

        if not pair.matched:
            self.counters.unmatched += pair.multiplicity
            pair.failure = Failures.nomatch
            return

        run = self._run

        # this is where R2 should start (if not a complete match, then r2.match_start will be > 0)
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        if r2_start_in_target < 0:
            self.counters.left_of_target += pair.multiplicity
            if run.count_left_prefixes:
                prefix = pair.r2.original_seq[0:0 - r2_start_in_target]
                self.counters.register_prefix(prefix, pair)
            if run.collapse_left_prefixes and (not run.collapse_only_prefixes or prefix in run._p_collapse_only_prefix_list):
                pair.r2._ltrim = 0 - r2_start_in_target
            else:
                pair.failure = Failures.left_of_zero
                return
        elif r2_start_in_target + pair.r2.original_len <= pair.target.n:
            # we're in the middle -- no problem
            pass
        elif not self._trim_adapters(pair):
            # we're over the right edge, and adapter trimming failed
            pair.failure = pair.failure or Failures.adapter_trim
            self.counters.adapter_trim_failure += pair.multiplicity
            return
        else:
            # we're at the right and trimming went ok, cool
            self.counters.adapter_trimmed += pair.multiplicity

        # set the match to be the rest of the (possibly trimmed) sequence, and count errors
        pair.r1.match_to_seq(reverse_complement = True)
        pair.r2.match_to_seq()
        target_seq = pair.target.seq
        r1_matcher = pair.r1.reverse_complement
        pair.r1.match_errors = string_match_errors(r1_matcher, target_seq[pair.r1.match_index:])
        pair.r2.match_errors = string_match_errors(pair.r2.subsequence, target_seq[pair.r2.match_index:])

        if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > run.allowed_target_errors:
            if pair.r1.match_errors:
                _debug("R1 errors: {}".format(pair.r1.match_errors))
            if pair.r2.match_errors:
                _debug("R2 errors: {}".format(pair.r2.match_errors))
            if run._p_v102_compat and not (pair.r1.match_errors + [e for e in pair.r2.match_errors if e < pair.r2.original_len - 4]):
                _debug("** v102 compat, allowing")
            else:
                pair.failure = Failures.match_errors
                self.counters.match_errors += pair.multiplicity
                return

        n = pair.target.n
        assert(pair.matched and pair.left >= 0 and pair.left <= n)

        if pair.right != n and not run.allow_multiple_rt_starts:
            pair.failure = Failures.right_edge
            self.counters.r1_not_on_right_edge += pair.multiplicity
            return

        if not pair.check_overlap():
            pair.failure = Failures.r1_r2_overlap
            return

        if run.count_mutations:
            pair.check_mutations()
            self.counters.low_quality_muts += pair.check_mutation_quality(self._run.mutations_require_quality_score)
            if pair.mutations and len(pair.mutations) > run.allowed_target_errors:
                pair.failure = Failures.match_errors
                self.counters.match_errors += pair.multiplicity
                return

        pair.site = pair.site or pair.left
        self.counters.register_count(pair)


#
# _process_pair_cotrans algorithm
#
# overview:
# (1) find linker in R1
# (2) find_partial R2 in target
#   (2a) if no match, find linker in R2
#     (2a1) if no linker, no match, return
#     (2a2) if linker, then determine target index, and goto 3
#   (2b) if match, the find_partial of R1 (w/o linker) in target
#     (2b1) if no match, return
#     (2b2) if match, then goto 3
# (3) R1/R2 are now fully determined, check for match/errors and return
#
# details:
#
# (1) this is just a string search of rc(R1) for the linker. it has to
#     be there to be a valid fragment.
#
# (2) this is the same algorithm (find_partial) as was used before for
#     non-cotrans experiments: it finds the longest subsequence of R2
#     that matches the target (not including any linker).
#
#   (2a) if we can't find a long enough match (>=10bp), it could be
#        because R1 and R2 are not far apart, and therefore R2
#        contains the linker. so we search for it (simple string
#        search)
#
#     (2a1) if R2 has no linker and doesn't have a long enough
#           subsequence matching the target, then it's not a match
#
#     (2a2) if the linker is found in R2, then everything to the left
#           of the linker is the part of the target to search. simple
#           string search will find that. if we don't find it -- or
#           find it multiple times (can be the case if it's very
#           short, say 3bp) -- then it's not a match,
#           return. otherwise, we now know where in the target R2
#           aligns, and therefore the linker, and therefore R1. go to
#           step [3]
#
#   (2b) in this case, we've found a subsequence match to the target
#        in R2, so we now know where it aligns. do the same
#        find_partial on R1 to find out where it aligns.
#
#     (2b1) if R1 doesn't have a sufficient-length matching
#           subsequence in the target, then it's not a match. (it
#           can't be too short, b/c we've found R2, and the only way
#           R1 can have a match that's too short is if R1 and R2
#           overlap, which corresponds to case [2a2] above.)
#
#     (2b2) if we do find a match, then we've aligned both R1 and R2
#           to the target. go to step [3]
#
# (3) we now know where R1 and R2 align with the target. if they are
#     near each other, R2 will include some or all of the linker;
#     verify that matches properly, if that's the case. if they are
#     even closer together, then adapter trimming will be required: in
#     which case, we know exactly how much of the adapters should be
#     on the end, we can check that it matches and trim it
#     off. finally, we make sure that R1 and R2 match the target
#     everywhere they're supposed to based on alignment. if all of
#     that passes, the a count is registered at the appropriate length
#     (right side where R1 matches the target) and site (left side of
#     where R2 matches the target).
#


class CotransPartialFindProcessor(PairProcessor):

    def prepare(self):
        self._targets.index()
        if 1 != len(self._targets.targets):
            raise Exception("multiple cotrans targets?")
        if 0 != self._run.allowed_adapter_errors:
            print("Warning: cotrans match w/adapter errors NYI")

    def process_pair(self, pair):

        # algorithm:
        # (1) find linker in R1
        # (2) find_partial R2 in target
        #   (2a) if no match, find linker in R2
        #     (2a1) if no linker, no match, return
        #     (2a2) if linker, then determine target index, and goto 3
        #   (2b) if match, the find_partial of R1 (w/o linker) in target
        #     (2b1) if no match, return
        #     (2b2) if match, then goto 3
        # (3) R1/R2 are now fully determined, check for match/errors and return
        #
        # xref full algorithm description at the end of the file for more detail
        #

        if not self._match_mask(pair)  or  not self._check_indeterminate(pair):
            return

        run = self._run

        linker = run.cotrans_linker
        linker_len = len(linker)
        target = self._targets.targets[0]
        tseq = target.seq

        r1_len = pair.r1.seq_len
        r1_rc = pair.r1.reverse_complement

        # (1)
        lindex = r1_rc.rfind(linker)
        if -1 == lindex:
            pair.failure = Failures.linker
            return

        r2_seq = pair.r2.subsequence
        r2_len = pair.r2.original_len
        rtrim = 0
        # (2)
        pair.r2.find_in_targets(self._targets, force_target = target)
        if not pair.r2.match_len:
            # (2a) if no match, then only chance is short match b/c of linker/trim
            l2index = r2_seq.rfind(linker)
            # TODO match errors in linker
            if -1 == l2index:
                pair.failure = Failures.nomatch
                return

            target_match_len = l2index
            target_match = r2_seq[:l2index]

            indices = string_find_errors(target_match, tseq, run.allowed_target_errors)
            if 0 == len(indices):
                pair.failure = Failures.nomatch
                return
            elif 1 == len(indices):
                index = indices[0]
            else:
                pair.failure = Failures.multiple_R1
                return

            # (2a2)
            rtrim = r2_len - l2index - linker_len - 4

            pair.r1.match_start = rtrim
            pair.r1.match_len = lindex - rtrim
            pair.r1.match_index = index

            pair.r2.match_start = 0
            pair.r2.match_len = l2index
            pair.r2.match_index = index

        else:

            # (2b)
            pair.r1.find_in_targets(self._targets, reverse_complement = True, force_target = target)
            if not pair.r1.match_len:
                pair.failure = Failures.nomatch
                return

            linker_start = pair.r1.match_index + (lindex - pair.r1.match_start)
            check_len = 0
            if pair.r2.match_index + (r2_len - pair.r2.match_start) > linker_start:
                linker_in_r2_idx = linker_start - pair.r2.match_index + pair.r2.match_start
                linker_check = r2_seq[linker_in_r2_idx:]
                pair.r2.linker_errors = string_match_errors(linker_check, linker)
                if len(pair.r2.linker_errors) > run.allowed_adapter_errors:
                    if run._p_v102_compat and not [e for e in pair.r2.linker_errors if e + linker_in_r2_idx < pair.r2.original_len - 4]:
                        _debug("** v102 compat, allowing")
                    else:
                        pair.failure = Failures.linker
                        return
                check_len = len(linker_check)
                if check_len > linker_len + 4:
                    rtrim = check_len - linker_len - 4

            delta = pair.r1.match_start - rtrim
            pair.r1.match_start -= delta
            pair.r1.match_len = lindex - pair.r1.match_start
            pair.r1.match_index -= delta

            delta = pair.r2.match_start
            pair.r2.match_start = 0
            pair.r2.match_index -= delta
            pair.r2.match_len = r2_len - check_len


        # (3) R1/R2 are now fully determined, check for match/errors and return
        pair.target = target

        # trim case
        if rtrim > 0:
            pair.r1.adapter_errors = string_match_errors(pair.r1.original_seq[-rtrim:], run.adapter_b)
            pair.r2.adapter_errors = string_match_errors(pair.r2.original_seq[-rtrim:], self._adapter_t_rc)
            if max(len(pair.r2.adapter_errors), len(pair.r1.adapter_errors)) > self._run.allowed_adapter_errors:
                if run._p_v102_compat and not [e for e in (pair.r1.adapter_errors + pair.r2.adapter_errors) if e < rtrim - 4]:
                    _debug("** v102 compat, allowing adapter errors in last 4")
                else:
                    pair.failure = Failures.adapter_trim
                    return
            pair.r1.trim(rtrim)
            pair.r2.trim(rtrim + 4) # also trim rc of handle

        if pair.right > target.n  or  pair.right < pair.left:
            pair.failure = Failures.right_edge
            return

        if pair.left < 0:
            self.counters.left_of_target += pair.multiplicity
            if run.count_left_prefixes:
                prefix = pair.r2.original_seq[0:0 - pair.left]
                self.counters.register_prefix(prefix, pair)
            if run.collapse_left_prefixes and (not run.collapse_only_prefixes or prefix in run._p_collapse_only_prefix_list):
                pair.r2._ltrim = 0 - pair.left
                pair.r2.match_to_seq()
            else:
                pair.failure = Failures.left_of_zero
                return

        if pair.right < run.cotrans_minimum_length:
            pair.failure = Failures.cotrans_min
            return

        if pair.r1.match_start + pair.r1.match_len + linker_len + rtrim < r1_len:
            pair.failure = Failures.right_edge
            return

        if pair.r1.match_len + linker_len != pair.r1.seq_len:
            pair.failure = Failures.extra_r1
            return

        pair.r1.match_errors = string_match_errors(pair.r1.reverse_complement[:pair.r1.match_len], tseq[pair.r1.match_index:])
        pair.r2.match_errors = string_match_errors(pair.r2.subsequence[:pair.r2.match_len], tseq[pair.r2.match_index:])
        if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > run.allowed_target_errors:
            if len(pair.r1.match_errors) > run.allowed_target_errors:
                _debug("R1 errors: {}".format(pair.r1.match_errors))
            if len(pair.r2.match_errors) > run.allowed_target_errors:
                _debug("R2 errors: {}".format(pair.r2.match_errors))
            if run._p_v102_compat and not pair.r1.match_errors and not [e for e in pair.r2.match_errors if e < pair.r2.original_len - 4]:
                _debug("** v102 compat, allowing")
            else:
                pair.failure = Failures.match_errors
                self.counters.match_errors += pair.multiplicity
                return

        if not pair.check_overlap():
            pair.failure = Failures.r1_r2_overlap
            return

        if run.count_mutations:
            pair.check_mutations()
            self.counters.low_quality_muts += pair.check_mutation_quality(run.mutations_require_quality_score)
            if pair.mutations and len(pair.mutations) > run.allowed_target_errors:
                pair.failure = Failures.match_errors
                self.counters.match_errors += pair.multiplicity
                return

        pair.site = pair.left
        self.counters.register_count(pair)
