from .util import align_strings, Indel, string_match_errors

class Sequence(object):

    def __init__(self, seq = None):
        self._reset(seq)

    def _reset(self, seq):
        self._seq = seq.upper() if seq else None
        self._length = len(seq) if seq else None
        self.quality = None
        self.match_start = None
        self.match_len = None
        self.match_index = None
        self.target = None
        self.match_errors = []
        self.indels = {}
        self.indels_delta = 0
        self._seq_with_indels = None
        self._quality_with_indels = None

    @property
    def target_len(self):
        return len(self.target) if self.target else 0

    def set_seq(self, seq):
        assert(seq)
        self._reset(seq)

    def _debug_print(self):
        print("\nOrig Seq = {}".format(self.original_seq))
        print("   len={}, needs_RC={}".format(len(self.original_seq), self._needs_rc))
        print("   ltrim={}, rtrim={}".format(self.ltrim, self.rtrim))
        print("   sublen={}".format(self.seq_len))
        print("Source match_start={}".format(self.match_start))
        print("Target length={}".format(self.target_len))
        print("Target match_index={}, match_len={}".format(self.match_index, self.match_len))
        print("Match Errors at: {}".format(self.match_errors))
        print("Indels (indels_delta={}):  {}".format(self.indels_delta, self.indels))
        print(" ")

    @property
    def original_seq(self):
        return self._seq

    @property
    def original_len(self):
        return self._length

    @property
    def seq_len(self):
        return self._length

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
    def fully_matched(self):
        return ((self.match_start == 0 and self.match_len == min(self.seq_len + self.indels_delta, self.target_len))  or
               (self.match_index == 0 and ((self.match_start + self.match_len == self.seq_len) or self.match_len == self.target_len)))

    def matchedLen(self):
        return self.match_len + self.indels_delta

    def match_to_seq(self):
        self.match_index -= self.match_start
        self.match_start = 0
        self.match_len = self.seq_len

    def resolve_ambig_indels(self, target):
        if not self.indels:
            return
        newindels = {}
        prevind = -1
        seq = self._seq
        for indind,indel in sorted(self.indels.items()):
            indel.ambiguous = False
            if indel.insert_type:
                ## TAI: Busan & Weeks, 2018 is about deletes, but go ahead and do same for inserts.
                while indind > 0  and  (indel.seq[-1] == target[indind - 1]  or  (indind - 1 - self.match_index) in self.match_errors):
                    indel.seq = seq[indel.src_index - 1] + indel.seq[:-1]
                    indind -= 1
                    indel.src_index -= 1
                    indel.ambiguous = True
                if not indel.ambiguous  and  ((indind < len(target)  and  target[indind] == indel.seq[0])  or  (indind - self.match_index) in self.match_errors):
                    # in case it started to the left already (from stitching)
                    indel.ambiguous = True
                newindels[indind] = indel
            else:
                ## Use heuristic from Busan & Weeks, 2018:  ambiguous deletes aligned to 5' side.
                ilen = len(indel.seq)
                indst = indind - ilen
                while indind >= ilen + prevind:
                    if target[indst:indind] == indel.seq:
                        pass
                    elif seq[indel.src_index - 1] == indel.seq[-1]:
                        ## For case like CT** --> C**T --> **CT  (target=CTCT)
                        indel.seq = target[indst:indind]
                    elif (indst - self.match_index) in self.match_errors:
                        indel.seq = target[indst:indind]
                        self.match_errors[self.match_errors.index(indst - self.match_index)] += len(indel.seq)
                    else:
                        break
                    indel.src_index -= 1
                    indel.ambiguous = True
                    indind -= 1
                    indst -= 1
                if not indel.ambiguous  and  ((indind < len(target) - 1  and  target[(indst + 2):(indind + 2)] == indel.seq)  or  (indind + 1 - self.match_index) in self.match_errors):
                    # in case it started to the left already (from stitching)
                    indel.ambiguous = True
                newindels[indind] = indel
            prevind = indind
        self.indels = newindels

    def _check_front_indel(self, alignment, read, target):
        tlen = len(target)
        ap = alignment.params
        if tlen - alignment.target_match_end - 1 > 0:
            delseq = target[alignment.target_match_end + 1:]
            existing_del = alignment.indels.get(alignment.target_match_end)
            oc = ap.gap_extend_cost if existing_del else ap.gap_open_cost
            if alignment.score <= oc + ap.gap_extend_cost * (len(delseq) - 1):
                return False
            assert(alignment.src_match_end == len(read) - 1)
            if existing_del:
                assert(not existing_del.insert_type)
                delseq = existing_del.seq + delseq
                del alignment.indels[alignment.target_match_end]
            self.indels[tlen - 1] = Indel(False, delseq, self.match_start)
            self.indels_delta -= len(delseq)
        elif len(read) - alignment.src_match_end - 1 > 0:
            insseq = read[alignment.src_match_end + 1:]
            if alignment.score <= ap.gap_open_cost + ap.gap_extend_cost * (len(insseq) - 1):
                return False
            assert(alignment.target_match_end == tlen - 1)
            self.indels[tlen] = Indel(True, insseq, alignment.src_match_end + 1)
            self.indels_delta += len(insseq)
        return True

    def _check_back_indel(self, alignment, read, target):
        ap = alignment.params
        if alignment.target_match_start > 0:
            delseq = target[:alignment.target_match_start]
            if alignment.score <= ap.gap_open_cost + ap.gap_extend_cost * (len(delseq) - 1):
                return False
            assert(alignment.src_match_start == 0)
            alignment.indels[alignment.target_match_start - 1] = Indel(False, delseq, 0)
            self.indels_delta -= len(delseq)
        elif alignment.src_match_start > 0:
            insseq = read[:alignment.src_match_start]
            if alignment.score <= ap.gap_open_cost + ap.gap_extend_cost * (len(insseq) - 1):
                return False
            assert(alignment.target_match_start == 0)
            alignment.indels.setdefault(0, Indel(True, "", 0)).seq += insseq
            self.indels_delta += len(insseq)
        return True

    def extend_alignment(self, target, ap, suffix = ""):
        self.target = target
        read_end = self.match_start + self.match_len + self.indels_delta
        if self.match_start > 0  and  self.match_index > 0:
            ## Extend towards left/front/5' end...
            front_read = self._seq[:self.match_start]
            front_target = target[:self.match_index]
            # Reverse strings for front-bias since their ends are more likely to align...
            align_front = align_strings(front_read[::-1], front_target[::-1], ap)
            align_front.flip()

            if align_front.score > 0.0  and  self._check_front_indel(align_front, front_read, front_target):
                self.match_index = align_front.target_match_start
                self.match_start = align_front.src_match_start
                self.match_len += len(front_target) - align_front.target_match_start
                self.indels_delta += align_front.indels_delta
                self.indels.update(align_front.indels)
                self.match_errors += [ err - self.match_index for err in align_front.mismatched ]
            elif ap.penalize_ends and ap.penalize_front_clip:  # TAI: add a different flag for this
                m = min(self.match_start, self.match_index)
                self.match_index -= m
                self.match_start -= m
                self.match_len += m
                self.match_errors += [ e for e in range(m) if front_read[e + self.match_start] != front_target[e + self.match_index] ]

        target_end = self.match_index + self.match_len
        if read_end < self.seq_len and (target_end < len(target)  or  len(suffix) > 0):
            ## Extend towards right/rear/3' end...
            back_read = self._seq[read_end:]
            back_target = target[target_end:] + suffix
            align_back = align_strings(back_read, back_target, ap)

            if align_back.score > 0.0  and  self._check_back_indel(align_back, back_read, back_target):
                self.match_len += align_back.target_match_end + 1
                self.indels_delta += align_back.indels_delta
                self.update_shifted_indels(align_back.indels, target_end, read_end)
                self.match_errors += [ (err + target_end - self.match_index) for err in align_back.mismatched if (err + target_end) < self.target_len ]
            elif ap.penalize_ends and ap.penalize_back_clip:  # TAI: add a different flag for this
                m = min(len(back_read), len(target[target_end:]))
                self.match_len += m
                self.match_errors += [ e + target_end - self.match_index for e in range(m) if back_read[e] != back_target[e] ]

        self.resolve_ambig_indels(target + suffix)

    def extendAlignmentToBestScore(self, target, ap):
        matchStart, matchLen, matchIndex = self.match_start, self.match_len, self.match_index
        originalLeft = self.match_start
        originalRight = self.match_start + self.match_len
        self.extend_alignment(target, ap)
        left = self.match_start
        right = self.match_start + self.match_len
        mutIndexes = self.match_errors

        # not efficient, but demonstrates
        mutCost, matchScore = -2, 2
        bestLeft = (0, left)
        scoreDelta = 0
        for idx in range(originalLeft - left - 1):
            if idx in mutIndexes:
                scoreDelta -= mutCost
            elif idx in self.indels:
                NYI # could just break/ignore? or deal with them -- wait for a case to come up
                break
            else:
                scoreDelta -= matchScore
            #print("  L @ {}/{} ~~> {}".format(idx, left + idx + 1, scoreDelta))
            if scoreDelta >= bestLeft[0]:
                bestLeft = (scoreDelta, left + idx + 1)
        #print("bestLeft:", bestLeft, left)

        bestRight = (0, right)
        scoreDelta = 0
        for idx in range(right - originalRight - 1):
            r = self.match_len - idx - 1
            if r in mutIndexes:
                scoreDelta -= mutCost
            elif r in self.indels:
                NYI # could just break/ignore? or deal with them -- wait for a case to come up
                break
            else:
                scoreDelta -= matchScore
            #print("  R @ {}/{} ~~> {}".format(idx, self.match_start + r, scoreDelta))
            if scoreDelta >= bestRight[0]:
                bestRight = (scoreDelta, self.match_start + r)
        #print("bestRight:", bestRight, right)

        if bestLeft[1] != left:
            delta = bestLeft[1] - left
            assert(delta > 0)
            self.match_start += delta
            self.match_len -= delta
            self.match_index += delta
            self.match_errors = [ e - delta for e in self.match_errors if e >= delta ]
            self.indels = { k - delta : v for k, v in self.indels.items() if k >= delta }

        if bestRight[1] != right:
            assert(right > bestRight[1])
            self.match_len -= (right - bestRight[1])
            self.match_errors = [ e for e in self.match_errors if e < self.match_len ]
            self.indels = { k : v for k, v in self.indels.items() if k < self.match_len }


    def update_shifted_indels(self, newindels, key_delta, val_delta):
        for indind,indel in newindels.items():
            indel.src_index += val_delta
            self.indels[indind + key_delta] = indel

    def shift_indels(self, val_delta):
        for _,indel in self.indels.items():
            indel.src_index += val_delta

    def trim_indels(self):
        trim_start = self.seq_len
        new_indels = {}
        for indind,indel in self.indels.items():
            if indel.insert_type:
                if indel.src_index > trim_start:
                    self._seq_with_indels = None
                    self.indels_delta -= len(indel.seq)
                elif indel.src_index + len(indel.seq) > trim_start:
                    self._seq_with_indels = None
                    self.indels_delta -= (len(indel.seq) - trim_start + indel.src_index)
                    indel.seq = indel.seq[:trim_start - indel.src_index]
                    new_indels[indind] = indel    # NOTE: This means we can now have inserts at the end of the target match!
                else:
                    new_indels[indind] = indel
            elif indel.src_index >= trim_start:
                self._seq_with_indels = None
                self.indels_delta += len(indel.seq)
            else:
                new_indels[indind] = indel
        self.indels = new_indels

    def indels_delta_before(self, target_index):
        result = 0
        for indind,indel in self.indels.items():
            if indind <= target_index:
                if indel.insert_type:
                    result += len(indel.seq)
                else:
                    result -= len(indel.seq)
        return result

    def apply_indels(self):
        if self._seq_with_indels:
            return self._seq_with_indels, self._quality_with_indels
        seq = self._seq
        qual = self.quality
        if not self.indels:
            self._seq_with_indels = seq
            self._quality_with_indels = qual
            return seq, qual
        nsl = []
        nql = []
        sind = self.match_start
        lind = self.match_start + self.match_len + self.indels_delta
        delta = 0
        for i in range(self.match_index, self.match_index + self.match_len):
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
                    if qual:
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
