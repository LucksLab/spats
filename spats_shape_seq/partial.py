
from processor import PairProcessor, Failures
from util import _warn, _debug, string_match_errors


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

        if r1_match_trimmed:
            # ok, we trimmed down our R1 due to adapters. need to see if that means the leftover matches
            # multiple targets; if so, need to reject this pair.
            target = pair.r1.find_in_targets(self._targets, reverse_complement = True, min_length_override = pair.r1.match_len)
            if not target or isinstance(target, list):
                # note that target may be None if pair.r1.match_len is less than the index word length
                _debug("dropping pair due to multiple R1 match after adapter trim")
                print "multiple R1 trim: " + pair.r1.original_seq
                print [r1_match_trimmed, r1_adapter_match, r1_length_to_trim]
                pair.target = None
                pair.failure = Failures.multiple_R1
                self.counters.multiple_R1_match += pair.multiplicity
                return False

        _debug("successful adapter trim of {}/{} bp from R1/R2".format(pair.r1._rtrim, pair.r2._rtrim))

        return True

    def _try_cotrans_adapter(self, pair):

        run = self._run
        linker = run.cotrans_linker
        linker_len = len(linker)

        r1_rc = pair.r1.reverse_complement
        lindex = r1_rc.rfind(linker)
        if -1 == lindex:
            pair.failure = Failures.linker
            return

        r2_seq = pair.r2.original_seq
        lindex = r2_seq.find(linker)
        if -1 == lindex:
            pair.failure = Failures.nomatch
            return

        # if R2 has the full linker, then this may be an adapter-trim case
        adapter_start = lindex + 20 + 4
        r2_len = pair.r2.original_len
        trim_length = 0
        if adapter_start < r2_len:
            if self._adapter_t_rc.startswith(r2_seq[adapter_start:]):
                trim_length = r2_len - adapter_start
                pair.r2.trim(trim_length + 4) # +4 b/c of handle
                if run.adapter_b.startswith(pair.r1.original_seq[-trim_length:]):
                    pair.r1.trim(trim_length)
                else:
                    pair.failure = Failures.adapter_trim
                    return
            else:
                pair.failure = Failures.adapter_trim
                return
        target_match_len = r2_len - linker_len - trim_length - 4
        if target_match_len < 4:
            pair.failure = Failures.nomatch
            return
        target_match = r2_seq[:target_match_len]
        if pair.r1.reverse_complement[:target_match_len] != target_match:
            pair.failure = Failures.mismatch
            return

        target = self._targets.targets[0]
        tseq = target.seq
        index = tseq.find(target_match)
        if -1 == index:
            pair.failure = Failures.nomatch
            return
        if -1 != tseq.find(target_match, index + 1):
            pair.failure = Failures.multiple_R1
            return

        pair.r2.match_index = index
        pair.r2.match_len = target_match_len
        pair.r2.match_start = 0
        pair.r1.match_index = index
        pair.r1.match_len = target_match_len
        pair.r1.match_start = trim_length
        pair.target = target
        pair.site = pair.site or pair.left
        self.counters.register_count(pair)


    def _process_pair_cotrans(self, pair):
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

        run = self._run
        # TODO: xref TODO below
        if 0 != run.allowed_target_errors:
            raise Exception("match w/errors NYI")
        if 1 != len(self._targets.targets):
            raise Exception("multiple cotrans targets?")

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

            # TODO: w/match errors
            index = tseq.find(target_match)
            if -1 == index:
                pair.failure = Failures.nomatch
                return
            if -1 != tseq.find(target_match, index + 1):
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
                linker_in_r2_idx = linker_start - pair.r2.match_index
                linker_check = r2_seq[linker_in_r2_idx:]
                pair.r2.linker_errors = string_match_errors(linker_check, linker)
                if len(pair.r2.linker_errors) > run.allowed_adapter_errors:
                    if run._v102_compat and not [e for e in pair.r2.linker_errors if e + linker_in_r2_idx < pair.r2.original_len - 4]:
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
                if run._v102_compat and not [e for e in (pair.r1.adapter_errors + pair.r2.adapter_errors) if e < rtrim - 4]:
                    _debug("** v102 compat, allowing adapter errors in last 4")
                else:
                    pair.failure = Failures.adapter_trim
                    return
            pair.r1.trim(rtrim)
            pair.r2.trim(rtrim + 4) # also trim rc of handle

        pair.r1.match_errors = string_match_errors(pair.r1.reverse_complement[:pair.r1.match_len], tseq[pair.r1.match_index:])
        pair.r2.match_errors = string_match_errors(pair.r2.subsequence[:pair.r2.match_len], tseq[pair.r2.match_index:])
        if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > run.allowed_target_errors:
            if len(pair.r1.match_errors) > run.allowed_target_errors:
                _debug("R1 errors: {}".format(pair.r1.match_errors))
            if len(pair.r2.match_errors) > run.allowed_target_errors:
                _debug("R2 errors: {}".format(pair.r2.match_errors))
            if run._v102_compat and not pair.r1.match_errors and not [e for e in pair.r2.match_errors if e < pair.r2.original_len - 4]:
                _debug("** v102 compat, allowing")
            else:
                pair.failure = Failures.match_errors
                self.counters.match_errors += pair.multiplicity
                return

        if pair.left < 0:
            pair.failure = Failures.left_of_zero
            return

        if pair.right < run.cotrans_minimum_length:
            pair.failure = Failures.cotrans_min
            return

        pair.site = pair.left
        self.counters.register_count(pair)


    def process_pair(self, pair):

        if not self._check_indeterminate(pair) or not self._match_mask(pair):
            return

        run = self._run
        cotrans = run.cotrans
        if cotrans:
            self._process_pair_cotrans(pair)
            return

        self._find_matches(pair)
        if not pair.matched:
            if cotrans:
                self._try_cotrans_adapter(pair)
                return
            self.counters.unmatched += pair.multiplicity
            pair.failure = Failures.nomatch
            return

        linker = run.cotrans_linker or ""
        linker_len = len(linker)

        if cotrans:
            # ok, R1 and R2 match
            # R1 should have linker on the end, adjust for that
            r1_rc = pair.r1.reverse_complement
            r1_end = pair.r1.match_start + pair.r1.match_len 
            delta = 0
            if len(r1_rc) - r1_end < linker_len:
                delta = linker_len - (len(r1_rc) - r1_end)
                pair.r1.match_len -= delta
                r1_end -= delta
            pair.r1.linker_errors = string_match_errors(linker, r1_rc[r1_end:r1_end+linker_len])
            if len(pair.r1.linker_errors) > run.allowed_adapter_errors:
                pair.failure = Failures.linker
                return
            pair.r1.match_len += linker_len

            if pair.r2.match_len < pair.r2.seq_len and pair.r2.match_index + pair.r2.seq_len > pair.r1.match_index:
                # r2 overlaps the linker. check and extend as necessary
                pair.r2.match_len -= delta
                pair.r2.linker_errors = string_match_errors(pair.r2.original_seq[pair.r2.match_start + pair.r2.match_len:], linker)
                if len(pair.r2.linker_errors) > run.allowed_adapter_errors:
                    if run._v102_compat and not [e for e in pair.r2.linker_errors if e + pair.r2.match_start + pair.r2.match_len < pair.r2.original_len - 4]:
                        _debug("** v102 compat, allowing")
                    else:
                        pair.failure = Failures.linker
                        return
                pair.r2.match_len = min(pair.r2.original_len - pair.r2.match_start, pair.r2.match_len + linker_len + 4)


        # this is where R2 should start (if not a complete match, then r2.match_start will be > 0)
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        if r2_start_in_target < 0:
            pair.failure = Failures.left_of_zero
            self.counters.left_of_target += pair.multiplicity
            return
        elif r2_start_in_target + pair.r2.original_len <= pair.target.n:
            # we're in the middle -- no problem
            pass
        elif cotrans and r2_start_in_target + pair.r2.original_len <= pair.target.n + linker_len:
            # we're still in the middle due to cotrans linker
            pass
        elif not self._trim_adapters(pair):
            # we're over the right edge, and adapter trimming failed
            pair.failure = pair.failure or Failures.adapter_trim
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
            r1_matcher = pair.r1.reverse_complement
            if cotrans:
                r1_matcher = r1_matcher[:-linker_len]
            pair.r1.match_errors = string_match_errors(r1_matcher, target_seq[pair.r1.match_index:])
            pair.r2.match_errors = string_match_errors(pair.r2.subsequence, target_seq[pair.r2.match_index:])

            if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > run.allowed_target_errors:
                if pair.r1.match_errors:
                    _debug("R1 errors: {}".format(pair.r1.match_errors))
                if pair.r2.match_errors:
                    _debug("R2 errors: {}".format(pair.r2.match_errors))
                if run._v102_compat and not (pair.r1.match_errors + [e for e in pair.r2.match_errors if e < pair.r2.original_len - 4]):
                    _debug("** v102 compat, allowing")
                else:
                    pair.failure = Failures.match_errors
                    self.counters.match_errors += pair.multiplicity
                    return

        n = pair.target.n
        assert(pair.matched and pair.left >= 0 and pair.left <= n)

        if cotrans:
            # we already verified it matches up with end (linker) above. just make sure about minimum len now.
            if pair.right <= run.cotrans_minimum_length:
                pair.failure = Failures.cotrans_min
                return

            # ok this is a match! subtract the linker length (L should not include linker)
            pair.r1.match_len -= linker_len

        else:
            # NOTE: this might change later due to "starts"
            if pair.right != n:
                pair.failure = Failures.right_edge
                self.counters.r1_not_on_right_edge += pair.multiplicity
                return

        pair.site = pair.site or pair.left
        self.counters.register_count(pair)
