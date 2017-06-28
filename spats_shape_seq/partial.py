
from processor import PairProcessor
from util import _warn, _debug, string_match_errors


class PartialFindProcessor(PairProcessor):

    def prepare(self):
        self._targets.index()

    def _find_matches(self, pair):
        # use R1 to determine which target
        target = pair.r1.find_in_targets(self._targets, reverse_complement = True)
        if isinstance(target, list):
            _debug("dropping pair due to multiple R1 match")
            pair.failure = "multiple R1 match"
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

        if r1_match_trimmed:
            # ok, we trimmed down our R1 due to adapters. need to see if that means the leftover matches
            # multiple targets; if so, need to reject this pair.
            target = pair.r1.find_in_targets(self._targets, reverse_complement = True, min_length_override = pair.r1.match_len)
            if not target or isinstance(target, list):
                # note that target may be None if pair.r1.match_len is less than the index word length
                _debug("dropping pair due to multiple R1 match after adapter trim")
                pair.target = None
                pair.failure = "multiple R1 match"
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
            pair.failure = "no match"
            return


        if self._run.cotrans:
            # ok, R1 and R2 match
            # R1 should have linker on the end, adjust for that
            r1_rc = pair.r1.reverse_complement
            r1_end = pair.r1.match_start + pair.r1.match_len 
            linker = self._run.cotrans_linker
            linker_len = len(linker)
            if len(r1_rc) - r1_end < linker_len:
                pair.failure = "R1 linker length failure"
                return
            pair.r1.linker_errors = string_match_errors(linker, r1_rc[r1_end:r1_end+linker_len])
            if len(pair.r1.linker_errors) > self._run.allowed_adapter_errors:
                pair.failure = "R1 linker match failure"
                return
            pair.r1.match_len += linker_len


        # this is where R2 should start (if not a complete match, then r2.match_start will be > 0)
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        if r2_start_in_target < 0:
            pair.failure = "R2 to left of site 0 failure on R2 for: {}".format(pair.identifier)
            self.counters.left_of_target += pair.multiplicity
            return
        elif r2_start_in_target + pair.r2.original_len <= pair.target.n:
            # we're in the middle -- no problem
            pass
        elif not self._trim_adapters(pair):
            # we're over the right edge, and adapter trimming failed
            pair.failure = pair.failure or "adapter trim failure"
            self.counters.adapter_trim_failure += pair.multiplicity
            return
        else:
            # we're at the right and trimming went ok, cool
            self.counters.adapter_trimmed += pair.multiplicity

        if pair.r2.match_len == pair.r2.seq_len  and  pair.r1.match_len == pair.r1.seq_len:
            # everything that wasn't adapter trimmed was matched -- nothing to do
            pass
        else:
            # set the match to be the rest of the (possibly trimmed) sequence, and count errors
            pair.r1.match_to_seq(reverse_complement = True)
            pair.r2.match_to_seq()
            target_seq = pair.target.seq
            pair.r1.match_errors = string_match_errors(pair.r1.reverse_complement, target_seq[pair.r1.match_index:])
            pair.r2.match_errors = string_match_errors(pair.r2.subsequence, target_seq[pair.r2.match_index:])

            if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > self._run.allowed_target_errors:
                if pair.r1.match_errors:
                    _debug("R1 errors: {}".format(pair.r1.match_errors))
                if pair.r2.match_errors:
                    _debug("R2 errors: {}".format(pair.r2.match_errors))
                pair.failure = "match errors failure"
                self.counters.match_errors += pair.multiplicity
                return

        n = pair.target.n
        assert(pair.matched and pair.left >= 0 and pair.left <= n)

        if self._run.cotrans:
            # we already verified it matches up with end (linker) above.
            pass
        else:
            # NOTE: this might change later due to "starts"
            if pair.right != n:
                pair.failure = "R1 right edge failure: {} - {}, n={}".format(pair.left, pair.right, n)
                self.counters.r1_not_on_right_edge += pair.multiplicity
                return

        pair.register_count(pair.right)
        self.counters.processed_pairs += pair.multiplicity
