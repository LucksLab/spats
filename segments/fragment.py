

from . import model as M
from .diagram import AlignmentDiagram, FragmentDiagram
from .failures import Failures
from .index import AlignmentFinder
from .logging import LoggingClass
from .sequence import Sequence
from .util import AlignmentParams, align_strings, string_match_errors, longest_match, reverse_complement, phredToQuality, _speedup


class FragmentMaker(LoggingClass):

    def __init__(self, r1IsRC = False):
        LoggingClass.__init__(self)
        #self.setLevelDebug()

        # configurable options
        self.r1IsRC = r1IsRC
        self.minimumExpectedOverlap = 8
        self.maxOverlap = 10 * self.minimumExpectedOverlap # TBD
        self.maxAllowedErrors = 3
        self.maxResolvedErrors = 10
        self.handleIndels = False
        self.allowIndeterminate = False # LL EXP 6: True
        self.allowContained = False # Barcode "First Batch - Bad-20220629"
        self.allowNonOverlapping = True
        self._qualityWarned = False
        self._fastHits = self._slowHits = self._noHits = 0
        self.qualityDeltaRequiredToResolveError = 10 # indicates 1/100 chance of being incorrect
        self.ap = AlignmentParams(gap_open_cost = 5, penalize_ends = False, penalize_front_clip = False, penalize_back_clip = False)

    def make(self, r1s, r2s, **kwargs):
        frag = M.Fragment(seq = M.Sequence(**kwargs), r1 = r1s, r2 = r2s, r1IsRC = self.r1IsRC)
        frag.seq.name = frag.seq.name or r1s.name.replace("/R1", "/F")
        frag.seq.sequenceId = frag.seq.sequenceId or r1s.sequenceId
        left, right = r1s.characters, r2s.characters
        leftq, rightq = r1s.quality, r2s.quality
        if self.r1IsRC:
            left = reverse_complement(left)
            leftq = leftq[::-1]
        else:
            right = reverse_complement(right)
            rightq = rightq[::-1]

        self._makeOverlapped(frag, left, right, leftq, rightq)

        if self.allowNonOverlapping and (not frag.success):
            # ok, if we have this much overlap (should indicate a run
            # of 8 or more, or multiple slightly shorter runs), then
            # we should consider this a genuine failure, rather than
            # consider this a nonoverlapping fragment.
            #
            # TODO: this may still need some work, but seems close for
            # now...maybe the thing to do is get info from the
            # fragment match and then redo the fragment alignment
            # based on that data?
            FRAGMENT_OVERLAP_FAILURE_SCORE_THRESHOLD = 64

            self.debugf("nonoverlapping fragment [{}], score={}", frag.failure, frag.alignment.score if frag.alignment else 0)
            if frag.alignment and (frag.alignment.score > FRAGMENT_OVERLAP_FAILURE_SCORE_THRESHOLD) and (frag.failure != Failures.indelsInFragment):
                self.debugf("considering fragment @{} with sufficient alignment score {} a true failure", frag.seq.sequenceId, frag.alignment.score)
                return frag

            # then we assume there was no overlap
            frag.errors = None
            frag.qualityResolved = None

            realIndel = False
            if (frag.failure == Failures.indelsInFragment) and (frag.overlap):
                # TODO: temp heuristic, could do better: from observed
                # cases, poor indel matches seem to be in the ~.02
                # range, and a poor but positive match
                # (s5_65k[7430]/F) is ~.1 -- so let's try .08 and see
                # how it goes. [ perfect is 1.0, matching a run of
                # half is 0.5^2 = 0.25 ]
                # some examples:
                #  negative:
                #    - s5_65k[3435]/F = 0.0541
                #    - s5_65k[1120]/F = 0.0248
                #    - s5_65k[1654]/F = 0.0197
                #    - s5_65k[220]/F  = 0.0185
                #    - s5_65k[6919]/F = 0.0159
                #  positive:
                #    - s5_65k[7430]/F = 0.1089
                #    - s5_65k[4305]/F = 0.1023
                REAL_INDEL_SCORE_RATIO = .08
                self.debug("indel check", frag.alignment.matchSize(), frag.alignment.score,
                           frag.overlap, "{:.4f}".format(frag.alignment.score / (frag.overlap * frag.overlap)))
                if frag.alignment.score / (frag.overlap * frag.overlap) > REAL_INDEL_SCORE_RATIO:
                    realIndel = True
            if realIndel:
                frag.success = False
            else:
                frag.overlap = 0
                frag.success = True
                frag.failure = None
            frag.seq.characters = (right + left[frag.overlap:]) if frag.reverseOrder else (left + right[frag.overlap:])

        if (not self.allowIndeterminate) and ('N' in frag.seq.characters): #(len(set(left)) > 4 or len(set(right)) > 4):
            frag.failure = Failures.indeterminate
            frag.success = False

        self.debug("frag out", frag.success, frag.overlap, len(r1s.characters), len(r2s.characters), " / ", len(frag.seq.characters), 
                   len(r1s.characters) + len(r2s.characters) - frag.overlap)

        assert(len(frag.seq.characters) == len(r1s.characters) + len(r2s.characters) - frag.overlap)
        return frag

    def alignmentHitReport(self):
        f, s, n = self._fastHits, self._slowHits, self._noHits
        div = (f + s + n) / 100
        if div:
            self.info("Alignment hits: {} ({:.2f}%) fast, {} ({:.2f}%) slow, {} ({:.2f}%) none".format(f, f/div, s, s/div, n, n/div))

    def _makeOverlapped(self, frag, left, right, leftq, rightq):
        # ok, we assume we've got the following picture:
        #  |       left            |
        #                | overlap |
        #                |            right      |
        #  |              fragment               |

        # we're trying to find the overlap. start with
        # `minimumExpectedOverlap`-sized chunks at the beginning of
        # `right`, and see if they are in left.

        # apparently, though, sometimes left and right are reversed.
        # so we look over the full right sequence

        # TODO: need to handle short overlap, which probably can't be
        # determined without looking at the experiment
        # descriptor. although, in that case, no overlap is probably
        # ok. TBD

        # try a fast aligner first...
        af = AlignmentFinder(left, minimumSegmentLength = 10)
        al = af.align(right)
        if not al:
            # ..and if it fails, use a conservative aligner
            af = AlignmentFinder(left, indexWordSize = 5, minimumSegmentLength = 5)
            al = af.align(right)
            if al:
                self._slowHits += 1
            else:
                self._noHits += 1
        else:
            self._fastHits += 1

        if not al:
            # we couldn't find a chunk, abort
            if not _speedup:
                self.debug("  no overlap")
                self.debug("    detail:\n{}\n{}\n{}".format(left, "" ,right))
            # TODO: in some applications this is ok
            frag.failure = Failures.noOverlap
            return frag

        if not _speedup:
            self.debug("alignment:\n" + AlignmentDiagram(right, af, al).make())

        frag.alignment = al
        overlap = al.segmentLength()
        frag.overlap = overlap
        leftStart, rightStart = al.segmentStart(), al.queryStart()
        segChars, queryChars = left, right

        if overlap != al.queryLength():
            frag.failure = Failures.indelsInFragment
            if leftStart > rightStart:
                od = rightStart + (len(left) - leftStart - overlap) + max(0, al.queryLength() - al.segmentLength())
                self.debug("indels", al.matchSize(), "/", frag.overlap, od, "->", frag.overlap + od, al.queryLength(), leftStart, rightStart)
                if frag.overlap + od > len(left):
                    frag.overlap = 0
                else:
                    frag.overlap += od
            else:
               frag.overlap = 0
            return frag

        if leftStart > rightStart:
            if not _speedup:
                self.debug("    detail:\n{}\n{}\n{}{}".format(left, (" " * leftStart) + "|" + (" " * (overlap - 2)) + "|", " " * + (leftStart - rightStart), right))
            frag.overlap = rightStart + (len(left) - leftStart)
            frag.seq.characters = left + right[frag.overlap:]
        else:
            if not _speedup:
                self.debug("    detail:\n{}{}\n{}\n{}".format(" " * (rightStart - leftStart), left, (" " * rightStart) + "|" + (" " * (overlap - 2)) + "|", right))
            if overlap > leftStart + (len(right) - rightStart - overlap):
                # if the overlap is significant, mark this as seeming like reverse order, maybe dig into it later...
                self.debug(" *** R1/R2 seem to overlap in reverse order...")
                frag.reverseOrder = True
                frag.overlap = leftStart + (len(right) - rightStart)
                frag.seq.characters = right + left[frag.overlap:]
            else:
                frag.failure = Failures.noOverlap
                return frag

        contained = False
        if (frag.overlap >= len(left)) or (frag.overlap >= len(right)):
            # case where R1/R2 is contained in the other
            contained = True
            if not self.allowContained:
                frag.failure = Failures.contained
                return frag
            if frag.overlap >= len(left):
                frag.overlap = len(left)
                frag.seq.characters = right
            else:
                frag.overlap = len(right)
                frag.seq.characters = left

        assert(len(frag.seq.characters) == len(left) + len(right) - frag.overlap)
        if not _speedup:
            self.debug("after left/right:\n" + FragmentDiagram(frag).make())

        gaps = al.gaps() or []
        if leftStart > rightStart:
            if 0 < rightStart:
                gaps = [ M.Gap(queryStart = al.queryStart() - rightStart,
                               queryEnd = al.queryStart(),
                               segmentStart = al.segmentStart() - rightStart,
                               segmentEnd = al.segmentStart()) ] + gaps
            if len(left) > leftStart + overlap:
                gaps.append(M.Gap(queryStart = al.queryEnd(),
                                  queryEnd = al.queryEnd() + len(left) - leftStart - overlap,
                                  segmentStart = al.segmentEnd(),
                                  segmentEnd = al.segmentEnd() + len(left) - leftStart - overlap))
            if (not gaps) and (left[-overlap:] != right[:overlap]):
                self.debug("ML={}\n{}\n{}".format(overlap, left[-overlap:], right[:overlap]))
                assert(("N" in left[-overlap:]) or ("N" in right[:overlap]))
                if not self._qualityWarned:
                    self.warn("TODO: should see if the N are corrected by quality...")
                    self._qualityWarned = True
                frag.failure = Failures.indeterminate
                return frag
        else:
            if 0 < leftStart:
                gaps = [ M.Gap(queryStart = al.queryStart() - leftStart,
                               queryEnd = al.queryStart(),
                               segmentStart = al.segmentStart() - leftStart,
                               segmentEnd = al.segmentStart()) ] + gaps
            if len(right) > rightStart + overlap:
                gaps.append(M.Gap(queryStart = al.queryEnd(),
                                  queryEnd = al.queryEnd() + len(right) - rightStart - overlap,
                                  segmentStart = al.segmentEnd(),
                                  segmentEnd = al.segmentEnd() + len(right) - rightStart - overlap))
            if (not gaps) and (right[-overlap:] != left[:overlap]):
                self.debug("ML={}\n{}\n{}".format(overlap, left[overlap:], right[:-overlap]))
                assert(("N" in right[-overlap:]) or ("N" in left[:overlap]))
                if not self._qualityWarned:
                    self.warn("TODO: should see if the N are corrected by quality...")
                    self._qualityWarned = True
                frag.failure = Failures.indeterminate
                return frag

        if not gaps:
            if not _speedup:
                self.debug("  perfect  match:", overlap)
            frag.success = True
            return frag

        overlap = frag.overlap
        if leftStart > rightStart:
            rightStart = 0
            leftStart = len(left) - overlap
        else:
            leftStart = 0
            rightStart = len(right) - overlap
        if not _speedup:
            self.debug(overlap, leftStart, rightStart, "gaps:", [ g.json() for g in gaps ])

        errors = []
        for gap in gaps:
            qs, qe, ss, se = gap.queryStart, gap.queryEnd, gap.segmentStart, gap.segmentEnd
            gapSize = qe - qs
            if gapSize != se - ss:
                frag.failure = Failures.indelsInFragment
                return frag
            for i in range(gapSize):
                if (qs + i >= len(queryChars)) or (ss + i >= len(segChars)):
                    assert(contained) # should only ever happen here
                    break
                if queryChars[qs + i] != segChars[ss + i]:
                    errors.append(min(qs, ss) + i)
                    if not _speedup:
                        self.debug(" added error: {} / {}".format(errors, gap.json()), qs, ss, queryChars[qs:qs + 4], segChars[ss:ss+4])
                    if len(errors) > (self.maxAllowedErrors + self.maxResolvedErrors):
                        frag.errors = errors
                        frag.failure = Failures.overlapErrors
                        return frag

        if errors:
            frag.errors = errors
            if not _speedup:
                self.debug("after errors: {}\n".format(errors) + FragmentDiagram(frag).make())

        qerrors = []
        qresolved = []
        fchars = frag.seq.characters
        for e in errors:
            lqp, rqp = leftq[leftStart + e], rightq[rightStart + e]
            lq, rq = phredToQuality(lqp), phredToQuality(rqp)
            delta = 0
            res = None
            if lq > rq:
                res = left[leftStart + e]
                delta = lq - rq
            else:
                delta = rq - lq
                res = right[rightStart + e]
            if frag.reverseOrder:
                residx = rightStart + e
            else:
                residx = leftStart + e
            if not _speedup:
                self.debug("qr error", leftStart, rightStart, e, residx, fchars[residx], "->", res)
            fchars = fchars[:residx] + res + fchars[residx + 1:]
            if delta >= self.qualityDeltaRequiredToResolveError:
                qresolved.append(e)
            else:
                qerrors.append(e)
            if not _speedup:
                self.debug("fragment Q resolve: ", lqp, rqp, lq, rq, "~~> ", res, qresolved, qerrors)
        if not _speedup:
            self.debug("quality resolved {}/{} errors".format(len(qresolved), len(errors)))

        frag.seq.characters = fchars
        frag.errorCount = len(qerrors)
        frag.errors = qerrors
        frag.qualityResolved = qresolved

        if not _speedup:
            self.debug("after qresolve:\n" + FragmentDiagram(frag).make())

        if (frag.errorCount > self.maxAllowedErrors) or (len(qresolved) > self.maxResolvedErrors):
            frag.failure = Failures.matchErrors
            return frag

        # TODO: on disagreement, we should take the one that matches the segd?

        frag.success = True
        return frag

