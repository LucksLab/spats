
from util import _debug, reverse_complement, align_strings, Indel

class Sequence(object):

    def __init__(self):
        self._reset(None, False)

    def _reset(self, seq, needs_rc):
        self._seq = seq.upper() if seq else None
        self._needs_rc = needs_rc
        self._seq_rc = None
        self._length = len(seq) if seq else None
        self._ltrim = 0
        self._rtrim = 0
        self.quality = None
        self.match_start = None
        self.match_len = None
        self.match_index = None
        self.target_len = None
        self.adapter_errors = []
        self.match_errors = []
        self.indels = {}              # note: .seq in Indel objects will be reverse-complemented for R1
        self.indels_delta = 0
        self._seq_with_indels = None
        self._quality_with_indels = None
        self.auto_adjust_match = False

    def set_seq(self, seq, needs_reverse_complement = False):
        self._reset(seq, needs_reverse_complement)

    def _debug_print(self):
        from util import objdict_as_str
        print("\nOrig Seq = {}".format(self.original_seq))
        print("   len={}, needs_RC={}".format(len(self.original_seq), self._needs_rc))
        print("   ltrim={}, rtrim={}".format(self.ltrim, self.rtrim))
        print("   sub={}".format(self.reverse_complement if self._needs_rc else self.subsequence))
        print("   sublen={}".format(self.seq_len))
        print("Source match_start={}".format(self.match_start))
        print("Target length={}".format(self.target_len))
        print("Target match_index={}, match_len={}".format(self.match_index, self.match_len))
        print("Match Errors at: {}".format(self.match_errors))
        print("Indels (indels_delta={}):  {}".format(self.indels_delta, objdict_as_str(self.indels)))
        print(" ")

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
    def subquality(self):
        if self._rtrim:
            return self.quality[self._ltrim:-self._rtrim]
        else:
            return self.quality[self._ltrim:]

    @property
    def seq_len(self):
        return self._length - self._ltrim - self._rtrim

    @property
    def reverse_complement(self):
        if not self._seq_rc:
            self._seq_rc = reverse_complement(self.subsequence)
        return self._seq_rc

    @property
    def reverse_quality(self):
        return self.subquality[::-1]

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
    def right_est(self):
        ''' Use this if match_len is not correct (cf. after trimming w/ indels) '''
        return self.match_index - self.match_start + self.seq_len - self.indels_delta

    @property
    def trimmed(self):
        return (0 != self._rtrim)

    @property
    def ltrim(self):
        return self._ltrim

    @ltrim.setter
    def ltrim(self, val):
        self._ltrim = val
        self._seq_rc = None

    @property
    def rtrim(self):
        return self._rtrim

    @rtrim.setter
    def rtrim(self, val):
        changed = (val != self._rtrim)
        self._rtrim = val
        if changed:
            self._seq_rc = None
            if self.auto_adjust_match:
                self.trim_indels()          # should only matter for R2 since R1 won't have indels
                self.match_errors = [ err for err in self.match_errors if err < self.seq_len - self.indels_delta ]
                if self.match_len and self.target_len:
                    self.match_len = min(self.match_len, self.target_len - self.match_index, self.seq_len - self.indels_delta)

    def find_in_targets(self, targets, force_target = None, min_length_override = 0):
        seq = self.reverse_complement if self._needs_rc else self.subsequence
        target, self.match_start, self.match_len, self.match_index = targets.find_partial(seq, force_target = force_target, min_length_override = min_length_override)
        return target

    @property
    def fully_matched(self):
        return self.match_start == 0 and self.match_len == min(self.seq_len + self.indels_delta, self.target_len) 

    def trim(self, length):
        self.rtrim = length
        if self._needs_rc:
            delta = self._rtrim - self.match_start
            _debug("trim reducing original match_len {} -> {}".format(self.match_len, self.match_len - delta))
            self.match_len -= delta
            self.match_start += delta
            self.match_index += delta
            return True
        return False

    def trim_to_match(self):
        if self._needs_rc:
            if self.match_start > 0:
                self.rtrim = self.match_start
                self.match_start = 0
            self.ltrim += (self.seq_len - self.match_len - self.indels_delta )
        else:
            if self.match_start > 0:
                self.ltrim = self.match_start
                self.match_start = 0
            self.rtrim += (self.seq_len - self.match_len - self.indels_delta )

    def match_to_seq(self):
        if self._needs_rc:
            if 0 == self._rtrim:
                # we've already done this if we've trimmed
                self.match_index -= self.match_start
                self.match_start = self._rtrim
            self.match_len = self.seq_len
        else:
            # TODO:  The next 2 lines are wrong if match_start is relative to subsequence.
            self.match_index -= (self.match_start - self._ltrim)
            self.match_start = self._ltrim
            self.match_len = self.seq_len
        _debug(["M2S:", self.match_index, self.match_len, self.match_start, "-- ", self._rtrim])

    def align_with_target(self, target, simfn, gap_open_cost, gap_extend_cost, suffix = ""):
        read_end = self.match_start + self.match_len
        if self.match_start > 0  and  self.match_index > 0:
            front_read = self.reverse_complement[:self.match_start] if self._needs_rc else self.subsequence[:self.match_start]
            front_target = target.seq[:self.match_index]
            # Reverse strings for front-bias since their ends are more likely to align...
            align_front = align_strings(front_read[::-1], front_target[::-1], simfn, gap_open_cost, gap_extend_cost)
            align_front.flip()

            good_alignment = True
            tlen = len(front_target)
            if tlen - align_front.target_match_end - 1 > 0:
                delseq = front_target[align_front.target_match_end + 1:]
                existing_del = align_front.indels.get(align_front.target_match_end)
                oc = gap_extend_cost if existing_del else gap_open_cost
                if 0.0 < align_front.score - oc - (len(delseq) - 1) * gap_extend_cost:
                    if existing_del:
                        assert(not existing_del.insert_type)
                        delseq = existing_del.seq + delseq
                        del align_front.indels[align_front.target_match_end]
                    self.indels[tlen - 1] = Indel(False, delseq, self.match_start)
                    self.indels_delta -= len(delseq)
                else:
                    good_alignment = False
            elif len(front_read) - align_front.src_match_end - 1 > 0:
                insseq = front_read[align_front.src_match_end + 1:]
                if 0.0 < align_front.score - gap_open_cost - (len(insseq) - 1) * gap_extend_cost:
                    self.indels[tlen] = Indel(True, insseq, align_front.src_match_end + 1)
                    self.indels_delta += len(insseq)
                else:
                    good_alignment = False
            if good_alignment:
                self.match_index = align_front.target_match_start
                self.match_start = align_front.src_match_start
                self.match_len += tlen - align_front.target_match_start
                self.indels_delta += align_front.indels_delta
                self.indels.update(align_front.indels)
                self.match_errors += [ err - self.match_index for err in align_front.mismatched ]
            else:
                m = min(self.match_start, self.match_index)
                self.match_index -= m
                self.match_start -= m
                self.match_len += m
                self.match_errors += [ e for e in xrange(m) if front_read[e + self.match_start] != front_target[e + self.match_index] ]

        target_end = self.match_index + self.match_len
        if read_end < self.seq_len and (target_end < target.n  or  len(suffix) > 0):
            back_read = self.reverse_complement[read_end:] if self._needs_rc else self.subsequence[read_end:]
            back_target = target.seq[target_end:] + suffix
            align_back = align_strings(back_read, back_target, simfn, gap_open_cost, gap_extend_cost)

            good_alignment = True
            if align_back.target_match_start > 0:
                delseq = back_target[:align_back.target_match_start]
                if 0.0 < align_back.score - gap_open_cost - (len(delseq) - 1) * gap_extend_cost:
                    align_back.indels[align_back.target_match_start - 1] = Indel(False, delseq, 0)
                    align_back.indels_delta -= len(delseq)
                else:
                    good_alignment = False
            elif align_back.src_match_start > 0:
                insseq = back_read[:align_back.src_match_start]
                existing_ins = align_back.indels.get(0)
                oc = gap_extend_cost if existing_ins else gap_open_cost
                if 0.0 < align_back.score - oc - (len(insseq) - 1) * gap_extend_cost:
                    if existing_ins:
                        existing_ins.seq += insseq
                    else:
                        align_back.indels[0] = Indel(True, insseq, 0)
                    align_back.indels_delta += len(insseq)
                else:
                    good_alignment = False
            if good_alignment:
                self.match_len += align_back.target_match_end + 1
                self.indels_delta += align_back.indels_delta
                self.update_shifted_indels(align_back.indels, target_end, read_end)
                self.match_errors += [ (err + target_end - self.match_index) for err in align_back.mismatched if (err + target_end) < self.target_len ]
            else:
                m = min(len(back_read), len(target.seq[target_end:]))
                self.match_len += m
                self.match_errors += [ e + target_end - self.match_index for e in xrange(m) if back_read[e] != back_target[e] ]

    def update_shifted_indels(self, newindels, key_delta, val_delta):
        for indind,indel in newindels.iteritems():
            indel.src_index += val_delta
            self.indels[indind + key_delta] = indel

    def trim_indels(self):
        trim_start = self.seq_len
        new_indels = {}
        for indind,indel in self.indels.iteritems():
            if indel.insert_type:
                if indel.src_index > trim_start:
                    self._seq_with_indels = None
                    self.indels_delta -= len(indel.seq)
                elif indel.src_index + len(indel.seq) > trim_start:
                    self._seq_with_indels = None
                    indel.seq = indel.seq[:trim_start - indel.src_index]
                    self.indels_delta -= len(indel.seq)
                    new_indels[indind] = indel    # NOTE: This means we can now have inserts at the end of the target match!
                else:
                    new_indels[indind] = indel
            elif indel.src_index >= trim_start:
                self._seq_with_indels = None
                self.indels_delta += len(indel.seq)
            else:
                new_indels[indind] = indel
        self.indels = new_indels

    def apply_indels(self):
        if self._seq_with_indels:
            return self._seq_with_indels, self._quality_with_indels
        if self._needs_rc:
            seq = self.reverse_complement
            qual = self.reverse_quality if self.quality else None
        else:
            seq = self.subsequence
            qual = self.subquality if self.quality else None
        if not self.indels:
            self._seq_with_indels = seq
            self._quality_with_indels = qual
            return seq, qual
        nsl = []
        nql = []
        sind = self.match_start
        lind = self.match_start + self.match_len + self.indels_delta
        delta = 0
        self._debug_print()
        for i in xrange(self.match_index, self.match_index + self.match_len):
            indel = self.indels.get(i, None)
            if indel:
                if indel.insert_type:
                    sind += len(indel.seq)
                    if sind < lind:
                        nsl.append(seq[sind])
                        if qual:
                            nql.append(qual[sind])
                        sind += 1
                else:
                    dind = delta + len(nsl) - len(indel.seq) + 1
                    nsl = nsl[:dind] + list(indel.seq)       # already reverse-complemented for R1
                    nql = nql[:dind] + ['!']*len(indel.seq)  # quality shouldn't matter b/c this can never be a mut
                    sind -= (len(indel.seq) - 1)
            elif sind < lind:
                nsl.append(seq[sind])
                if qual:
                    nql.append(qual[sind])
                sind += 1
            else:
                delta += 1
        self._seq_with_indels = "".join(nsl)
        self._quality_with_indels = "".join(nql) if qual else None
        return self._seq_with_indels, self._quality_with_indels

