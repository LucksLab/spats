
from . import model as M
from .diagram import AlignmentDiagram, SegmentDiagram
from .index import AlignmentFinder
from .failures import Failures
from .logging import LoggingClass
from .sequence import Sequence
from .util import AlignmentParams, Indel, align_strings, string_match_errors, longest_match, crossProduct, _speedup


class SegmentResult(M.SegmentResult):

    def __repr__(self):
        if self.unmatched:
            return "|??? @ {}|".format(self.segment.key)
        base = "|{}-{} @ {}:{}".format(self.queryStart, self.queryEnd, self.segment.key, self.segmentStart)
        if self.matchTruncatedCount:
            base += " T{}".format(self.matchTruncatedCount)
        if self.errors:
            base += " E{}".format(self.errors)
        if self.indels:
            base += " I{}".format([ e.errorIndex for e in self.indels ])
        base += "|"
        return base

    @property
    def unmatched(self):
        return (self.queryStart == -1)

    @property
    def length(self):
        return self.queryEnd - self.queryStart - self.indelsDelta

    def isFullLength(self):
        return self.length == self.segment.minLength

    @property
    def segmentEnd(self):
        return self.segmentStart + self.length

    def errorCount(self):
        return len(self.errors or []) + len(self.indels or [])

    def indelsCount(self):
        return len(self.indels or [])

    def rightExtendWithAlignment(self, rng, alignment):
        assert(rng.segment.key == self.segment.key)
        assert(rng.start >= self.segmentEnd) # f"{rng.start}/{self.segmentEnd} -- {self}/{rng.json()}"
        assert(self.segment.wildcard or (self.segment.characters[rng.start:rng.end] == rng.text))
        amount = rng.end - rng.start
        myShift = self.length
        newErrors = [ e + myShift for e in alignment.mismatched if e < amount ]
        newIndels = [ i.shiftAll(myShift) for i in alignment.indels.values() if i.errorIndex < amount ]
        newDelta = sum([ i.delta for i in newIndels ])
        #print("*** REWL in", self, amount, myShift, newErrors, newIndels, newDelta, self.indelsDelta)
        self.queryEnd += (amount + newDelta)
        if self.segment.wildcard:
            assert(not newIndels)
        else:
            self.errors += newErrors
            self.indels += newIndels
            self.indelsDelta += newDelta
        #print("*** REWL out", self, self.indelsDelta, self.length, self.segmentEnd, self.segment.minEnd)
        alignment.mismatched = [ e - amount for e in alignment.mismatched if e >= amount ]
        alignment.indels = { l - amount : i.shiftAll(-amount) for l, i in alignment.indels.items() if i.errorIndex >= amount }
        return newDelta

    def leftExtendWithAlignment(self, rng, alignment):
        assert(rng.segment.key == self.segment.key)
        assert(rng.end <= self.segmentStart)
        assert(self.segment.wildcard or (self.segment.characters[rng.start:rng.end] == rng.text))
        amount = rng.end - rng.start
        newErrors = [ e for e in alignment.mismatched if e < amount ]
        newIndels = [ i for i in alignment.indels.values() if i.errorIndex <= amount ]
        newDelta = sum([ i.delta for i in newIndels ])
        #print("*** LEWL in", self, amount, newErrors, newIndels, newDelta)
        self.queryStart -= (amount + newDelta)
        self.segmentStart -= amount
        if self.segment.wildcard:
            assert(not newIndels)
        else:
            self.errors = newErrors + self.errors
            self.indels = newIndels + self.indels
            self.indelsDelta += newDelta
        #print("*** LEWL out", self)
        alignment.mismatched = [ e - amount for e in alignment.mismatched if e >= amount ]
        alignment.indels = { l - amount : i.shiftAll(-amount) for l, i in alignment.indels.items() if i.errorIndex > amount }
        return newDelta

    def leftEdgePossibleRanges(self):
        segd = self.segment
        assert(not segd.wildcard)
        chars = segd.characters
        possibles = []
        for start in range(0, min(self.segmentStart, segd.maxStart + 1)):
            end = self.segmentStart
            slen = end - start
            if slen and (slen + self.length <= segd.maxLength):
                possibles.append(M.Range(segment = segd, start = start, end = end, text = chars[start:end]))
        return possibles

    def rightEdgePossibleRanges(self):
        possibles = []
        segd = self.segment
        if segd.wildcard:
            NYI # needs checking against a real case
            possibles.append(M.Range(segment = segd, start = 0, end = segd.minLength))
            return possibles
        chars = segd.characters
        for end in range(max(self.segmentEnd, segd.minEnd), len(chars) + 1):
            start = self.segmentEnd
            slen = end - start
            if slen and (slen + self.length <= segd.maxLength):
                possibles.append(M.Range(segment = segd, start = start, end = end, text = chars[start:end]))
        return possibles

    def _possibleRangesForSegment(self, segd):
        if segd.wildcard:
            return [ M.Range(segment = segd, end = sl, text = sl * "?") for sl in range(segd.minLength, segd.maxLength + 1) ]
        possibles = []
        chars = segd.characters
        for start in range(0, segd.maxStart + 1):
            for end in range(segd.minEnd, len(segd.characters) + 1):
                slen = end - start
                if (slen >= segd.minLength) and (slen <= segd.maxLength):
                    possibles.append(M.Range(segment = segd, start = start, end = end, text = chars[start:end]))
        return possibles

    def canShiftRight(self, query, segd):
        return (query[self.queryStart + 1:self.queryEnd + 1] == segd.characters)

    def shiftRight(self, amount):
        self.queryStart += amount
        self.queryEnd += amount


class SequenceResult(M.SequenceResult):

    def __init__(self, sequence, seqd):
        M.SequenceResult.__init__(self, sequence = sequence)
        self.keyedSegments = { s.key : s for s in seqd.segments }

    def displayString(self):
        lines = []
        if self.success:
            lines.append("SUCCESS:")
        else:
            lines.append("FAIL: {}".format(self.failure))
        seq = self.sequence.characters
        for sr in self.segments:
            lines.append(" {}[{}] @ {}-{}: {}".format(sr.segment.key, sr.segmentStart, sr.queryStart, sr.queryEnd, seq[sr.queryStart:sr.queryEnd]))
        return "\n".join(lines)

    def segmentForResult(self, sr):
        return self.keyedSegments[sr.segmentKey]

    def segmentWithKey(self, key):
        for sr in self.segments:
            if sr.segment.key == key:
                return sr
        return None

    def subsequenceForSegment(self, srOrSrKey):
        sr = srOrSrKey
        if isinstance(srOrSrKey, str):
            sr = self.segmentWithKey(srOrSrKey)
        return self.sequence.characters[sr.queryStart:sr.queryEnd]

    def errorCount(self):
        return sum([ sr.errorCount() for sr in self.segments])

    def unmatchedPrefix(self):
        seq = self.sequence.characters
        if not self.segments:
            return seq
        return seq[:self.segments[0].queryStart]

    def unmatchedSuffix(self):
        seq = self.sequence.characters
        if not self.segments:
            return seq
        return seq[self.segments[-1].queryEnd:]

    def leftGap(self, segment):
        assert(segment in self.segments)
        i = self.segments.index(segment)
        if i == 0:
            return segment.queryStart
        else:
            return segment.queryStart - self.segments[i - 1].queryEnd

    def rightGap(self, segment):
        assert(segment in self.segments)
        i = self.segments.index(segment)
        if i == len(self.segments) - 1:
            return len(self.sequence.characters) - segment.queryEnd
        else:
            return self.segments[i + 1].queryEnd - segment.queryStart

    def hasErrors(self):
        for sr in self.segments:
            if sr.errors or sr.indels:
                return True
        return False

    def hasIndels(self):
        for sr in self.segments:
            if sr.indels:
                return True
        return False

    def matchTruncatedCount(self):
        return sum([sr.matchTruncatedCount for sr in self.segments])

    def indelsCount(self):
        return sum([sr.indelsCount() for sr in self.segments ])

    def validate(self):
        seq = self.sequence.characters
        for sr in self.segments:
            if sr.unmatched:
                raise Exception("segment with no match?", sr.json())
            queryPart = seq[sr.queryStart:sr.queryEnd]
            seg = sr.segment
            segPart = seg.characters[sr.segmentStart:sr.segmentStart + sr.length]
            if sr.errors:
                # TODO: apply indels
                if 1:
                    continue
                # not quite right, but can do better later...
                lastError = max(sr.errors)
                firstError = min(sr.errors)
                segPart = segPart[:firstError] + segPart[lastError + 1:]
                queryPart = queryPart[:firstError] + queryPart[lastError + 1:]
            elif sr.indels:
                # TODO
                #print("indels:", sr.indels)
                continue
            if segPart == queryPart:
                continue
            if ("N" not in segPart) and ("N" not in queryPart):
                raise Exception("Validation failure: {} \n{}  != \n{} {} {}".format(self.sequence.name, queryPart, segPart, sr.errors, sr.json()))

    def diagram(self):
        return SegmentDiagram(self).make()


class SegmentProcessor(LoggingClass):

    def __init__(self, seqd):
        LoggingClass.__init__(self)
        #self.setLevelDebug()
        #self.setLevelVerbose()
        self._validate(seqd)
        self.seqd = seqd
        self._metadata = [ M.SegmentMetadata(segd = s) for s in seqd.segments ]
        self._setup()

    def _validate(self, seqd):
        assert(isinstance(seqd, M.SequenceDescriptor))
        segs = seqd.segments
        assert(segs)
        assert(len(set([s.key for s in segs])) == len(segs))
        for s in segs:
            self._validateSeg(s)

    def _validateSeg(self, s):
        try:
            assert(s.key)
            assert(s.minLength <= s.maxLength)
            if s.wildcard:
                assert((not s.characters) and (not s.index))
                assert(0 == s.maxStart)
                assert(0 == s.minEnd)
                assert(0 == s.maxAllowedErrors)
                assert(not s.handleIndels)
                #assert(s.minLength == s.maxLength) # variable-length wildcard NYI
            else:
                assert(s.characters)
                assert(s.minLength <= len(s.characters))
                assert(s.maxLength <= len(s.characters))
                if s.minLength == s.maxLength:
                    assert(0 == s.maxStart)
                    assert(s.minLength == s.minEnd)
        except:
            self.warn("Invalid segment descriptor:", s.json())
            raise

    def _validateResult(self, res):
        cur = 0
        for sr in res.segments:
            assert(sr and isinstance(sr, M.SegmentResult))
            segd = sr.segment
            if not _speedup:
                self.verbose("validate:", sr.json())
                self.verbose("  errors: ", sr.errors)
                self.verbose("  indels: ", sr.indels)
            assert(cur == sr.queryStart)
            assert(sr.queryEnd >= sr.queryStart)
            assert(sr.length >= segd.minLength)
            assert(sr.length <= segd.maxLength)
            assert(sr.segmentStart <= segd.maxStart)
            assert(sr.segmentEnd >= segd.minEnd)
            for e in (sr.errors or []):
                assert(e >= 0 and e < sr.length)
            assert(len(sr.errors or []) + len(sr.indels or []) <= segd.maxAllowedErrors)
            cur = sr.queryEnd
            if segd.wildcard:
                continue
            source = segd.characters[sr.segmentStart:sr.segmentEnd]
            dest = res.sequence.characters[sr.queryStart:sr.queryEnd]
            if sr.errors or sr.indels:
                # TODO: assert the applied alignment matches
                pass
            else:
                assert(source == dest), (segd.key, source, "/", dest)
        assert(cur == len(res.sequence.characters))

    def _setup(self):
        for i in range(len(self._metadata)):
            md = self._metadata[i]
            if md.segd.index and (not md.aligner):
                md.aligner = AlignmentFinder(md.segd.characters)
                if len(set(md.segd.characters)) == 1:
                    # special case handling for long constant
                    # sequences, e.g. the loop in SpatsUMI; the index
                    # will just return the first occurrence, so we may
                    # be able to shift it, if necessary, to fill gaps.
                    # note only shifts to the *right* will be
                    # possible, since the index will find the first to
                    # the left
                    md.shiftableConstant = True
                    assert(0 == md.segd.maxAllowedErrors) # otherwise NYI

        self.ap = AlignmentParams(simfn = lambda n1, n2: 2 if n1 == n2 else (0 if (n1 == '?' or n2 == '?') else -2),
                                  gap_open_cost = 5, gap_extend_cost = 1,
                                  penalize_ends = True, penalize_front_clip = True, penalize_back_clip = True)

    def _findAlignment(self, md, query, result):
        if not md.aligner:
            return None
        candidate = md.aligner.align(query)
        if not candidate:
            self.debug(" _findAlignment no match: {}".format(md.segd.key))
            return None

        if not _speedup:
            self.debug("alignment:\n" + AlignmentDiagram(query, md.aligner, candidate).make())
        if (md.segd.maxAllowedErrors <= 0) and (len(candidate.segments) > 1):
            # if we're not allowed errors, use the longest match
            ssort = sorted(candidate.segments, key = lambda s : s.length, reverse = True)
            candidate.segments = ssort[:1]
            if not _speedup:
                self.debug(" --> revised to:\n" + AlignmentDiagram(query, md.aligner, candidate).make())

        first = candidate.segments[0]
        last = candidate.segments[-1]
        ss, qs = first.segmentIndex, first.queryIndex
        se, qe = last.segmentIndex + last.length, last.queryIndex + last.length
        return SegmentResult(segment = md.segd,
                             alignment = candidate,
                             segmentStart = ss,
                             queryStart = qs,
                             queryEnd = qe,
                             errors = [],
                           indels = [],
                             indelsDelta = (qe - qs) - (se - ss))


    def _findAlignedSegments(self, query, result):
        # find all indexed matches
        result.segments = []
        lastHit = None
        leftQuerySpot, rightQuerySpot = 0, 0
        for md in self._metadata:
            segd = md.segd
            if not md.aligner:
                result.segments.append(None)
                continue
            sr = self._findAlignment(md, query, result)
            result.segments.append(sr)
            if not sr:
                if segd.minLength > 0:
                    result.failure = Failures.nomatch
                    result.failedSegment = segd.key
                    return
                continue
            if (sr.queryStart < leftQuerySpot) or (sr.queryEnd < rightQuerySpot):
                result.failure = Failures.outOfOrder
                result.failedSegment = segd.key
                return
            if (segd.minLength > 0) and (sr.length < (segd.minLength >> 3)):
                # this detects very short matches, that are almost certainly not extendable
                result.failure = Failures.shortMatch
                result.failedSegment = segd.key
                return

            assert(not sr.errorCount())
            self._verboseState(result, "_findAlignment")
            # if we have inside gaps, fix them now
            if len(sr.alignment.segments) > 1:
                self._alignSegmentGaps(query, sr)
                self._verboseState(result, "_alignSegmentGaps")

            leftQuerySpot = sr.queryStart
            rightQuerySpot = sr.queryEnd

        self.debug("_findAlignedSegments", result.segments)
        if not [ x for x in result.segments if x ]:
            # no matches, drop it
            result.failure = Failures.nomatch

    def _shiftAsNecessary(self, query, result):
        shifted = False
        for idx in range(len(result.segments)):
            md = self._metadata[idx]
            if not md.shiftableConstant:
                continue
            sr = result.segments[idx]
            if not sr:
                continue
            if 0 == idx:
                self.debug("Not shifting shiftable as there's no room to the left")
                continue

            # ok, we have a shiftable constant. let's see if neighbors
            # want us to shift.
            prevMd = self._metadata[idx - 1]
            prevSr = result.segments[idx - 1]
            nextMd = self._metadata[idx + 1] if idx + 1 < len(self._metadata) else None
            nextSr = result.segments[idx + 1] if idx + 1 < len(self._metadata) else None
            prevFixed = result.segments[idx - 2].queryEnd if idx >= 2 else 0
            nextFixed = result.segments[idx + 2].queryStart if idx + 2 < len(self._metadata) else len(query)

            while True:
                # shift as long as it seems like it's better if we do
                prevGap = sr.queryStart - (prevSr.queryEnd if prevSr else prevFixed)
                nextGap = (nextSr.queryStart if nextSr else nextFixed) - sr.queryEnd
                self.debug("_san", md.segd.key, prevGap, nextGap)
                if (prevGap < prevMd.segd.minLength) and (nextGap > nextMd.segd.maxLength) and sr.canShiftRight(query, md.segd):
                    sr.shiftRight(1)
                    shifted = True
                else:
                    break

            self.debug("_san out, gaps are now", prevGap, nextGap)
            if shifted:
                self._verboseState(result, "_shift")

        return shifted

    def _alignGaps(self, query, result):
        curQuery = 0
        leftSegmentIndex = -1
        for idx in range(len(result.segments)):
            sr = result.segments[idx]
            if not sr:
                continue
            if curQuery > sr.queryStart:
                # TODO: this is not necessarily fatal, but wait to see a case where it happens
                result.failure = Failures.interestingQueryOverlap
                result.failedSegment = sr.segment.key
            elif curQuery < sr.queryStart:
                self._alignInterSegmentGap(M.Gap(queryStart = curQuery, queryEnd = sr.queryStart,
                                                 segmentStart = leftSegmentIndex, segmentEnd = idx),
                                           query,
                                           result)

            if result.failure:
                return
            curQuery = sr.queryEnd
            leftSegmentIndex = idx
        if curQuery < len(query):
            self._alignInterSegmentGap(M.Gap(queryStart = curQuery, queryEnd = len(query),
                                             segmentStart = leftSegmentIndex, segmentEnd = len(self._metadata)),
                                       query,
                                       result)

    def _possibleRangesForSegment(self, segd):
        if segd.wildcard:
            return [ M.Range(segment = segd, end = sl, text = sl * "?") for sl in range(segd.minLength, segd.maxLength + 1) ]
        possibles = []
        chars = segd.characters
        for start in range(0, segd.maxStart + 1):
            for end in range(segd.minEnd, len(segd.characters) + 1):
                slen = end - start
                if (slen >= segd.minLength) and (slen <= segd.maxLength):
                    possibles.append(M.Range(segment = segd, start = start, end = end, text = chars[start:end]))
        return possibles

    def _alignInterSegmentGap(self, gap, query, result):
        # | ..left.. | |m0| |m1| | ..right.. |
        if not _speedup:
            self.debug("_aisg in", gap.json())
        srSegds = [ (None, self._metadata[i].segd, i) for i in range(gap.segmentStart + 1, gap.segmentEnd) ]
        if gap.segmentStart != -1:
            sr = result.segments[gap.segmentStart]
            srSegds = [ (sr, sr.segment, gap.segmentStart) ] + srSegds
        if gap.segmentEnd != len(self._metadata):
            sr = result.segments[gap.segmentEnd]
            srSegds = srSegds + [ (sr, sr.segment, gap.segmentEnd) ]

        # determine all possible segment texts
        allWildcard = True
        segmentRanges = [ [] ]
        for sr, segd, idx in srSegds:
            self.debug("_aisg seg", sr, idx, segd.key)
            allWildcard = allWildcard and segd.wildcard
            if sr:
                if idx == gap.segmentEnd:
                    possibles = sr.leftEdgePossibleRanges()
                else:
                    possibles = sr.rightEdgePossibleRanges()
            else:
                possibles = self._possibleRangesForSegment(segd)
            self.debug("  -> possibles:", [r.json() for r in possibles])
            segmentRanges = crossProduct(segmentRanges, possibles)
        assert(segmentRanges)
        self.debug("_aisg ranges:", len(segmentRanges), len(segmentRanges[0]), segmentRanges[:4])

        # find the best scoring
        queryText = query[gap.queryStart:gap.queryEnd]
        best = None
        assert(len(segmentRanges) < 10) # just want to see a case when we hit this, can probably do a lot of optimizing
        for ranges in segmentRanges:
            segText = "".join([ r.text for r in ranges ])
            a = align_strings(queryText, segText, self.ap, hack = True)
            if (not best) or a.score > best[1].score:
                best = (ranges, a)

        # apply it
        assert(best)
        ranges, a = best
        self.debug("_aisg alignment:\n", queryText, "\n", segText, "\n", a)
        if a.src_match_start != 0 or a.target_match_start != 0 or \
           (a.src_match_end + 1) != len(queryText) or (a.target_match_end + 1 != len(segText)):
            result.failure = Failures.gapAlignment
            result.failedSegment = ":".join([ sr[1].key for sr in srSegds ])
            return

        curQuery = gap.queryStart
        for rng in ranges:
            if not rng.text:
                continue
            for ss in srSegds:
                if ss[1] == rng.segment:
                    sr, segd, idx = ss
                    break
            assert(segd)
            if not sr:
                rngSize = rng.end - rng.start
                if segd.wildcard and [ i for i in a.indels.values() if i.errorIndex < rngSize ]:
                    if (idx + 1 < len(self._metadata)) and self._metadata[idx + 1].shiftableConstant and self._shiftAsNecessary(query, result):
                        # ok we shifted, start over; this recursion will end when it works, or we can't shift anymore
                        self.debug("recursively re-aligning gaps after a shift...")
                        self._alignGaps(query, result)
                        return
                    if not _speedup:
                        self.debug("wildcard has indels, failure", a.indels, rng.json())
                    result.failure = Failures.gapTooSmall if a.indels_delta < 0 else Failures.gapTooBig
                    result.failedSegment = segd.key
                    return
                sr = SegmentResult(segment = segd, queryStart = curQuery, queryEnd = curQuery)
                result.segments[idx] = sr
            assert(segd == sr.segment)
            if idx == gap.segmentEnd:
                delta = sr.leftExtendWithAlignment(rng, a)
            else:
                delta = sr.rightExtendWithAlignment(rng, a)
            curQuery = curQuery + (rng.end - rng.start) + delta
            if not _speedup:
                self.debug("_aisg", len(a.indels), curQuery, sr, gap.json())
                self._verboseState(result, "_aisg")
        if a.indels:
            assert(1 == len(a.indels))
            i = a.indels[0]
            last = srSegds[-1][0]
            if not last:
                result.failure = Failures.gapTooSmall if a.indels_delta < 0 else Failures.gapTooBig
                result.failedSegment = sr.segment.key
                return
            assert((0 == i.errorIndex) and last)
            delta = last.leftExtendWithAlignment(M.Range(segment = last.segment, text = ""), a)
            curQuery = curQuery + delta
            self._verboseState(result, "_aisg edge indel")
        assert(curQuery == gap.queryEnd), f"{curQuery} / {gap.queryEnd}"

    def _alignSegmentGaps(self, query, sr):
        for gap in sr.alignment.gaps():
            queryText = query[gap.queryStart:gap.queryEnd]
            segText = sr.segment.characters[gap.segmentStart:gap.segmentEnd]
            if not segText:
                sr.indels.append(Indel(True, queryText, gap.queryStart, gap.segmentStart, errorIndex = gap.segmentStart))
                continue
            if not queryText:
                sr.indels.append(Indel(False, segText, gap.queryStart, gap.segmentStart, errorIndex = gap.segmentStart))
                continue
            if (1 == len(segText)) and (1 == len(queryText)):
                # subst of one is a mut
                sr.errors.append(gap.segmentStart)
                continue
            a = align_strings(queryText, segText, self.ap, hack = True)
            if not _speedup:
                self.debug(a.score, (self.ap.gap_open_cost + len(segText) + len(queryText)))
            # TODO: is this the right metric?
            #if a.score <= 0 - (self.ap.gap_open_cost + len(segText) + len(queryText)):
            if len(a.indels) + len(a.mismatched) >= 3:
                # the idea being, one subst is more plausible than 3+ other events
                # 2 other events (e.g., insert and mut) is about as plausible
                self.debug("_asg using a subst for poor alignment", a.score)
                subst = Indel(True, queryText, gap.queryStart, gap.segmentStart, errorIndex = gap.segmentStart)
                subst.substSize = len(segText)
                sr.indels.append(subst)
                continue
            if not _speedup:
                self.debug("_asg", sr, gap.json(), list(a.indels.values()), "\n", queryText, "\n", segText, "\n", a)
            sr.errors += [ e + gap.segmentStart for e in a.mismatched ]
            sr.indels += [ i.shiftToGap(gap) for i in a.indels.values() ]
        assert(sr.indelsDelta == sum([i.delta for i in sr.indels]))

    def _processInternal(self, querySeq):
        assert(isinstance(querySeq, M.Sequence))
        result = SequenceResult(querySeq, self.seqd)
        query = querySeq.characters
        self._findAlignedSegments(query, result)
        if not result.failure:
            self._shiftAsNecessary(query, result)
        if not result.failure:
            self._alignGaps(query, result)
        if not result.failure:
            self._checkSegmentRequirements(result)
        if not result.failure:
            self._checkForEdgeIndels(result)
        return result

    def _checkForEdgeIndels(self, res):
        assert(not res.failure)
        for sr in res.segments:
            if (not sr) or (not sr.indels):
                continue
            if sr.segment.wildcard:
                res.failure = Failures.wildcardIndel
                res.failedSegment = sr.segment.key
                return
            for i in sr.indels:
                if (sr.segmentStart == i.errorIndex) or \
                   (sr.segmentEnd <= i.errorIndex + (i.substSize if i.insertType else i.size)):
                    res.failure = Failures.edgeIndel
                    res.failedSegment = sr.segment.key
                    return

    def process(self, querySeq):
        result = self._processInternal(querySeq)
        result.segments = [ sr for sr in result.segments if (sr and sr.length > 0) ]
        for sr in result.segments:
            sr.errors = sorted(sr.errors)
            sr.indels = sorted(sr.indels, key = lambda i : i.errorIndex)
        if not result.failure:
            try:
                self._validateResult(result)
                result.success = True
            except Exception as e:
                result.failure = Failures.validation
                self.warn("validation failure, please fix", str(e))
                raise
        self.debug("process out: f=", result.failure)
        return result

    def _verboseState(self, result, tag):
        if (not _speedup) and self.isLevelVerbose():
            self.verbose(tag + "\n" + SegmentDiagram(result).make())

    def _checkSegmentRequirements(self, res):
        assert(not res.failure)
        for md, sr in zip(self._metadata, res.segments):
            if not sr:
                if md.segd.minLength > 0:
                    res.failure = Failures.minSegmentLength
                    res.failedSegment = md.segd.key
                    return
                continue
            self.debug("CSR", md.segd.key, "/",
                       md.segd.maxStart, sr.segmentStart, "/",
                       md.segd.minEnd, sr.segmentEnd)
            if sr.segmentStart > md.segd.maxStart:
                res.failure = Failures.maxStart
            elif sr.segmentEnd < md.segd.minEnd:
                res.failure = Failures.minEnd
            elif sr.errorCount() > md.segd.maxAllowedErrors:
                res.failure = Failures.matchErrors
            elif sr.length < md.segd.minLength:
                res.failure = Failures.minSegmentLength
            elif sr.length > md.segd.minLength:
                res.failure = Failures.maxSegmentLength
            else:
                continue
            res.failedSegment = md.segd.key
            return
