
from .failures import Failures
from .util import reverse_complement, M, LoggingClass, phredToQuality

CUT_SITE = 26

class _DiagramSupport(LoggingClass):

    def __init__(self):
        LoggingClass.__init__(self)
        self.setLevelDebug()

    def _setup(self, fullWidth):
        self.lines = []
        self.fullWidth = fullWidth

    def _add(self, txt):
        self.lines.append(txt + (" " * max(0, self.fullWidth - len(txt))))

    def _spacer(self):
        self._add("")

    def _placeOn(self, txt, base, loc):
        if loc > len(base):
            res = base + (' ' * (loc - len(base))) + txt
        else:
            def chmerge(t, b):
                if t not in ' _':
                    return t
                elif b not in ' _':
                    return b
                elif t == '_' or b == '_':
                    return '_'
                return ' '
            merged = ''.join([ chmerge(t, b) for t, b in zip(txt, base[loc:loc + len(txt)] + ' ' * len(txt)) ])
            res = base[:loc] + merged + base[loc + len(txt):]
        return res

    def _centerOn(self, txt, base):
        return self._placeOn(txt, base, max(0, len(base) - len(txt)) >> 1)

    def _center(self, txt):
        return (" " * max(0, self.center - (len(txt) >> 1))) + txt

    def _addCentered(self, txt):
        self._add(self._center(txt))

    def _trimmable(self, lines, idx):
        for l in lines:
            if (idx < len(l)) and (l[idx] not in [ ' ', '-' ]):
                return False
        return True

    def _trim(self):
        lines = self.lines
        cur = 0
        margin = 2
        lastKeep = 0
        trimmed = 0
        maxLen = max([len(l) for l in lines])
        trimmables = [ i for i in range(maxLen) if self._trimmable(lines, i) ]
        tset = set(trimmables)
        margined = []
        for i in trimmables:
            ok = True
            for j in range(margin):
                if ((i + 1 + j) not in tset) or ((i - 1 - j) not in tset):
                    ok = False
            if ok:
                margined.append(i)
        trimmed = []
        for l in lines:
            delta = 0
            for i in margined:
                l = l[:i-delta] + l[i-delta+1:]
                delta += 1
            trimmed.append(l)
        self.lines = trimmed

    def _mergeLines(self, allWordLines, linesToMerge, newLines):
        if not linesToMerge:
            assert(not allWordLines)
            return newLines, newLines
        merged = []
        for cl, nl in zip(linesToMerge, newLines):
            if self._overlaps(cl, nl):
                return allWordLines + newLines, newLines
            merged.append(self._mergeLine(cl, nl))
        lc = len(linesToMerge)
        assert(allWordLines[0 - lc:] == linesToMerge)
        return allWordLines[:0 - lc] + merged, merged

    def _mergeLine(self, cl, nl):
        mlen = max(len(cl), len(nl))
        if len(cl) < mlen:
            cl = cl + (' ' * (mlen - len(cl)))
        if len(nl) < mlen:
            nl = nl + (' ' * (mlen - len(nl)))
        return ''.join([ cl[i] if ' ' == nl[i] else nl[i] for i in range(mlen) ])

    def _overlaps(self, cl, nl):
        for i in range(min(len(cl), len(nl))):
            if cl[i] != ' ' and nl[i] != ' ' and cl[i] != nl[i]:
                return True
        return False

    def _gap(self, lines, ranges, minToGap = 8, minToKeep = 4, scoreSpacing = 5):
        cur = minToKeep
        elidedRegions = []

        ranges = sorted(ranges, key = lambda r : r[0])
        for rloc, rlen in ranges + [ (len(lines[0]), 0) ]:
            #print(loc, cur, loc - cur - minToKeep)
            rend = rloc + rlen + scoreSpacing + minToKeep
            if rloc - cur - minToKeep > minToGap:
                elidedRegions.append((cur, rloc - minToKeep))
            cur = rend
        #print(elidedRegions)

        gapLine = ''
        numLine = ''
        for gapStart, gapEnd in elidedRegions:
            #print("gapping", gapStart - len(gapLine), gapEnd-gapStart)
            gapLine += '-' * (gapStart - len(gapLine))
            gapLine += 'X' * (gapEnd - gapStart)
            assert(len(gapLine) == gapEnd)
            numLine += str(len(numLine))
            numLine += ' ' * (gapEnd - len(numLine))
        gapLine += '-' * (len(lines[0]) - len(gapLine))
        numLine += str(len(numLine))
        numLine += ' ' * (len(lines[0]) - len(numLine))
        #print(len(gapLine), len(lines[0]))
        #print(elidedRegions)
        #assert(len(gapLine) == len(lines[0]))
        #lines = [ numLine, gapLine ] + lines + [ gapLine ]
        lines = [ numLine ] + lines

        #print(elidedRegions)
        delta = 0
        for gapStart, gapEnd in elidedRegions:
            degapped = []
            #print("gap: ", gapStart, gapEnd)
            for line in lines:
                #print("  ", delta, gapStart, gapEnd)
                start, end = line[:gapStart + delta], line[gapEnd + delta:]
                #print("    ", len(start), len(end))
                insertChar = '.' if end and end[0] != ' ' else ' '
                insert = (insertChar * 2)
                degapped.append(start + insert + end)
            lines = degapped
            delta += len(insert) - (gapEnd - gapStart)
        return lines


class SegmentDiagram(_DiagramSupport):

    def __init__(self, result, fragment = None, classifyResult = None):
        _DiagramSupport.__init__(self)
        self.result = result
        self.fragment = fragment
        self.classifyResult = classifyResult
        self.prefix = ' ' * 20
        self._reads = False
        self._alignment = None

    def make(self):
        self._setup(128)
        plen = len(self.prefix)

        res = self.result
        seqObj = res.sequence
        f = self.fragment
        l = self.prefix + seqObj.characters
        hdr = "{}:".format((f.seq.name or f.seq.identifier) if f else (seqObj.name or seqObj.identifier))
        if len(hdr) >= len(self.prefix):
            hdr = hdr[2-len(self.prefix):]
        l = self._placeOn(hdr, l, 0)
        self._add(l)

        for sr in res.segments:
            if not sr:
                continue
            sr.shiftedErrors = sr.errors or []
            if not sr.indels:
                continue
            for i in sr.indels:
                i.shiftedIndex = i.errorIndex
            for i in sr.indels:
                sr.shiftedErrors = [ e + (i.delta if e > i.errorIndex else 0) for e in sr.shiftedErrors ]
                for other in sr.indels:
                    if other.errorIndex > i.errorIndex:
                        other.shiftedIndex += i.delta

        for sridx in range(len(res.segments)):
            sr = res.segments[sridx]
            if not sr:
                continue
            seg = sr.segment
            if seg.wildcard:
                segLen = sr.length
            else:
                segLen = len(seg.characters)
            matchLen = sr.queryEnd - sr.queryStart
            l = "| {}:".format(seg.key)
            if sr.matchTruncatedCount:
                l += "*{}".format(sr.matchTruncatedCount)
            start = sr.segmentStart
            end = start + matchLen - sr.indelsDelta
            segPart = ("N" * sr.length) if seg.wildcard else seg.characters[start:end]
            for e in (sr.indels or []):
                loc = e.shiftedIndex
                if e.substSize:
                    segPart = segPart[:loc] + (" " * (len(e.seq))) + segPart[loc + e.substSize:]
                elif e.insertType:
                    segPart = segPart[:loc] + (" " * len(e.seq)) + segPart[loc:]
                else:
                    segPart = segPart[:loc] + segPart[loc + len(e.seq):]
            segIndex = plen + sr.queryStart
            spacer = "|"
            bars = []
            for idx2 in range(sridx + 1, len(res.segments)):
                sr2 = res.segments[idx2]
                if not sr2:
                    continue
                bars += [ plen + sr2.queryStart, plen + sr2.queryEnd - 1 ]
                for e in (sr2.shiftedErrors or []):
                    bars.append(plen + sr2.queryStart + e)
                for i in (sr2.indels or []):
                    bars.append(plen + sr2.queryStart + i.shiftedIndex)
            for bar in bars:
                spacer = self._placeOn("|", spacer, bar)
                l = self._placeOn("|", l, bar)
            for bar in bars + [ plen + sr.queryStart, plen + sr.queryEnd - 1 ]:
                spacer = self._placeOn("|", spacer, bar)
            for e in (sr.shiftedErrors or []):
                spacer = self._placeOn("!", spacer, plen + sr.queryStart + e)
            for e in (sr.indels or []):
                if e.substSize:
                    spacer = self._placeOn("+{}-{}".format(len(e.seq), e.substSize), spacer, plen + sr.queryStart + e.shiftedIndex)
                else:
                    sign = "+" if e.insertType else "-"
                    spacer = self._placeOn("{}{}".format(sign, len(e.seq)), spacer, plen + sr.queryStart + e.shiftedIndex)
            if start > 0:
                assert(not seg.wildcard)
                extra = min(start, 4)
                left = seg.characters[start - extra:start].lower()
                if extra < start:
                    left = "..." + left
                segPart = left + segPart
                segIndex -= len(left)
            if end < segLen:
                assert(not seg.wildcard)
                extra = min(segLen - end, 16)
                right = seg.characters[end:end+extra].lower()
                if end + extra < segLen:
                    right += "..."
                segPart = segPart + right
            l = self._placeOn(segPart, l, segIndex)
            self._add(spacer)
            self._add(l)
            if sr.segment.minLength == sr.segment.maxLength:
                continue
            start, end = sr.segmentStart, sr.segmentStart + matchLen - sr.indelsDelta
            l = "|"
            for bar in bars:
                l = self._placeOn("|", l, bar)
            if "target" == seg.key:
                l = self._placeOn("X", l, plen + sr.queryStart + CUT_SITE - sr.segmentStart)
            starter = str(start) + "^"
            l = self._placeOn(starter, l, plen + sr.queryStart - len(starter) + 1)
            ender = "^" + str(end)
            l = self._placeOn(ender, l, plen + sr.queryEnd - 1)
            if sr.errors:
                errors = sr.errors if len(sr.errors) <= 3 else (sr.errors[0], sr.errors[-1])
                for error in errors:
                    erridx = sr.segmentStart + error
                    l = self._placeOn("^" + str(erridx), l, plen + sr.queryStart + error)
            if sr.indels:
                errors = [ x.errorIndex for x in sr.indels ][:3]
                for error in errors:
                    erridx = sr.segmentStart + error
                    l = self._placeOn("^" + str(erridx), l, plen + sr.queryStart + error)
            self._add(l)

        self._add("|")
        if self._reads:
            self._add("\---> {} reads matched".format(len(res.segments)))
        elif self._alignment:
            self._add("\---> score={}".format(self._alignment.score))
        elif res.success:
            l = "\---> SUCCESS"
            if self.classifyResult:
                l += " [{}{}]".format(self.classifyResult.category, " / {}".format(self.classifyResult.description) if self.classifyResult.description else "")
            if hasattr(res, 'errorCount') and res.errorCount():
                l += " ({} errors)".format(res.errorCount())
            self._add(l)
                
        elif res.failure:
            self._add("\---> FAIL: {} [ {}{} ]".format(Failures.description(res.failure),
                                                       res.failure,
                                                       " / {}".format(res.failedSegment) if res.failedSegment else ""))
        else:
            self._add("\---> ???")

        return '\n'.join(self.lines)


class ReadsDiagram(SegmentDiagram):
    def __init__(self, sequence, readsResults):
        _DiagramSupport.__init__(self)
        # TODO: cleanup; for now just translate, same code should work
        sr = M.SequenceResult(sequence = sequence, success = True)
        for rr in readsResults:
            sr.segments.append(M.SegmentResult(segment = M.SegmentDescriptor(key = rr.item.key,
                                                                             characters = rr.item.characters),
                                               queryStart = rr.queryStart,
                                               queryEnd = rr.queryEnd,
                                               segmentStart = rr.itemStart,
                                               errors = rr.errors,
                                               indels = rr.indels,
                                               indelsDelta = rr.indelsDelta))
        SegmentDiagram.__init__(self, sr)
        self._reads = True

class AlignmentDiagram(SegmentDiagram):
    def __init__(self, query, index, a):
        _DiagramSupport.__init__(self)
        # TODO: cleanup; for now just translate, same code should work
        sr = M.SequenceResult(sequence = M.Sequence(characters = query))
        for aseg in a.segments:
            sr.segments.append(M.SegmentResult(segment = M.SegmentDescriptor(key = "q{}-{}, s{}-{}".format(aseg.queryIndex,
                                                                                                           aseg.queryIndex + aseg.length,
                                                                                                           aseg.segmentIndex,
                                                                                                           aseg.segmentIndex + aseg.length),
                                                                             characters = index.characters),
                                               queryStart = aseg.queryIndex,
                                               queryEnd = aseg.queryIndex + aseg.length,
                                               segmentStart = aseg.segmentIndex))
        SegmentDiagram.__init__(self, sr)
        self._alignment = a


class FragmentDiagram(_DiagramSupport):

    def __init__(self, fragment):
        _DiagramSupport.__init__(self)
        self.fragment = fragment
        self.prefix = ' ' * 20

    def make(self):
        self._setup(128)
        plen = len(self.prefix)
        f = self.fragment
        lefts, rights = f.r1, f.r2
        left, right = lefts.characters, rights.characters
        leftq, rightq = lefts.quality, rights.quality
        if f.r1IsRC:
            left = reverse_complement(left)
            leftq = leftq[::-1]
        else:
            right = reverse_complement(right)
            rightq = rightq[::-1]
 
        l = self._placeOn("{}:".format(f.seq.name or f.seq.identifier), self.prefix, 0)
        self._add(l)
        r1Start = 0
        r2Start = 0
        overlapStart = 0
        if f.overlap:
            if f.reverseOrder:
                r1Start = len(right) - f.overlap
                overlapStart = r1Start
            else:
                r2Start = len(left) - f.overlap
                overlapStart = r2Start
        l = self._placeOn("  rc(R1):" if f.r1IsRC else "  R1:", "|", 0)
        l = self._placeOn(left, l, plen + r1Start)
        self._add(l)

        l = "|"
        if f.overlap:
            l = self._placeOn("|", l, plen + r1Start)
            if f.reverseOrder:
                l = self._placeOn("|", l, plen + r1Start + f.overlap - 1)
            else:
                l = self._placeOn("|", l, plen + r2Start)
            l = self._placeOn("|", l, plen + r1Start + len(f.r1.characters) - 1)

        errsAndResolved = sorted((f.errors or []) + (f.qualityResolved or []))
        if not errsAndResolved:
            self._add(l)
        else:
            leftqs = { err : leftq[(0 if f.reverseOrder else len(left) - f.overlap) + err] for err in errsAndResolved }
            rightqs = { err : rightq[(len(right) - f.overlap if f.reverseOrder else 0) + err] for err in errsAndResolved }
            rql = l
            lql = l
            last = 0
            for err in errsAndResolved:
                lql = self._placeOn(leftqs[err], lql, plen + overlapStart + err)
                rql = self._placeOn(rightqs[err], rql, plen + overlapStart + err)
                l = self._placeOn("!" if (err in f.errors) else ("v" if rightqs[err] > leftqs[err] else "^"), l, plen + overlapStart + err)
                last = max(last, plen + overlapStart + err)
            self._add(self._placeOn("(qual):", lql, 4))
            self._add(self._placeOn(str([ phredToQuality(leftqs[err]) - phredToQuality(rightqs[err]) for err in errsAndResolved]), l, last + 3))
            self._add(self._placeOn("(qual):", rql, 4))

        l = self._placeOn("  R2:" if f.r1IsRC else "  rc(R2):", "|", 0)
        l = self._placeOn(right, l, plen + r2Start)
        if f.overlap:
            if f.reverseOrder:
                l = self._placeOn("|", l, plen + r1Start + len(f.r1.characters) - 1)
            else:
                l = self._placeOn("|", l, plen)
        self._add(l)

        l = "|"
        if f.overlap:
            l = self._placeOn("|", l, plen + r1Start)
            l = self._placeOn("|", l, plen + r1Start + len(f.r1.characters) - 1)
            l = self._placeOn("|", l, plen + r2Start)
            l = self._placeOn("|", l, plen + r2Start + len(f.r2.characters) - 1)
            if f.reverseOrder:
                l = self._placeOn("|", l, plen + r1Start + f.overlap - 1)
            else:
                l = self._placeOn("|", l, plen + r2Start + f.overlap - 1)
        for qr in (f.qualityResolved or []):
            l = self._placeOn("|", l, plen + qr + (r1Start if f.reverseOrder else r2Start))
        self._add(l)

        if f.overlap:
            l = self._placeOn("  F:", "|", 0)
            l = self._placeOn(f.seq.characters, l, plen)
            self._add(l)
            self._add("|")

        if f.success:
            if f.overlap:
                l = "\---> SUCCESS: {} overlap, {} total length".format(f.overlap, len(f.seq.characters))
            else:
                l = "\---> SUCCESS (no overlap)"
            if f.errorCount:
                l += " ({} errors)".format(f.errorCount)
            if f.reverseOrder:
                l += "  (** reverse order)"
            self._add(l)
        else:
            self._add("\---> FAIL: fragment / {}".format(f.failure))
        return '\n'.join(self.lines)

    def foo():
#|        |R1 |
#|   start| ! |     veidx
#|        | target  |
#|    sidx^     | ! |
#\-> site, end ({} errors)

        if 0:
            self._add("\---> SUCCESS: {}-{}{}{}".format(site, end, "  ({} errors)".format(errorCount) if errorCount else ""))
        else:
            self._add("\---> FAIL: {}".format(self.r1.failure or self.r2.failure))
