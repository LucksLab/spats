
from .diagram import AlignmentDiagram
from .util import longest_match, longestExactMatch, LoggingClass, M, _speedup


class Index(LoggingClass):

    def __init__(self, characters, index_word_length = 8):
        LoggingClass.__init__(self)
        #self.setLevelDebug()
        assert(characters)
        index_word_length = min(len(characters), index_word_length)
        self.characters = characters
        self._index = {}
        self._index_word_length = index_word_length
        self._minimum_length = index_word_length
        assert(self._index_word_length > 0)
        assert(self._minimum_length > 0)

    def build(self):
        index = {}
        word_len = self._index_word_length
        seq = self.characters
        n = len(seq)
        for i in range(n - word_len + 1):
            key = seq[i:(i + word_len)]
            sites = index.get(key)
            if not sites:
                sites = []
                index[key] = sites
            sites.append(i)
        assert(index)
        self._index = index

    def find_partial(self, query, min_length_override = 0):
        # it's much faster to search for longer partial matches, then fall back on the minimum
        candidate = self._find_partial(query, min_length_override, multiple = 2)
        candidate = candidate if candidate[1] else self._find_partial(query, min_length_override, multiple = 1)
        return candidate if candidate[1] else None

    #@profile
    def _find_partial(self, query, min_length_override, multiple):
        min_len = min_length_override or self._minimum_length * multiple
        word_len = self._index_word_length
        check_every = max(min_len - word_len, 1) # norah has proved that this guarantees finding a match if it exists
        query_len = len(query)
        last = query_len - max(check_every, word_len)
        check_sites = list(range(0, last, check_every))
        check_sites.append(last)
        candidate = [None, None, None]
        # NOTE: it's important to check all sites, and all hits -- to find the longest match.
        for site in check_sites:
            site_key = query[site:site+word_len]
            #self.verbose("CS: {}, {}".format(site, site_key))
            for index in self._index.get(site_key, []):
                self.debug("GOT: " + str(index))
                left, right = longest_match(query, (site, word_len), self.characters, (index, word_len))
                total_len = left + right + word_len
                self.debug("extends: <--{}, -->{} / {} ({})".format(left, right, total_len, min_len))
                if total_len >= min_len:
                    if not candidate[1] or total_len > candidate[1]:
                        # keep it if it's the best match so far
                        # candidate = [ start, matchLen, segmentStart ]
                        candidate = [site - left, total_len, index - left]
                        self.debug("C: {}".format(candidate))
        return candidate


    def find_all_partials(self, query, min_length_override = None):
        assert(query)
        min_len = min_length_override or self._minimum_length
        word_len = self._index_word_length
        check_every = max(min_len - word_len, 1) # norah has proved that this guarantees finding a match if it exists
        query_len = len(query)
        last = query_len - max(check_every, word_len)
        check_sites = list(range(0, last, check_every))
        check_sites.append(last)
        self.debug("FAP", word_len, min_len, check_every, query_len)
        assert((word_len > 0) and (min_len > 0) and (query_len > 0))
        candidates = []
        max_hit_right = 0
        # NOTE: it's important to check all sites, and all hits -- to find the longest match.
        for site in check_sites:
            if site < max_hit_right:
                continue
            site_key = query[site:site+word_len]
            #self.verbose("CS: {}, {}".format(site, site_key))
            for index in self._index.get(site_key, []):
                self.debug("GOT: " + site_key + "@" + str(index))
                left, right = longest_match(query, (site, word_len), self.characters, (index, word_len))
                total_len = left + right + word_len
                cleft = site - left
                cright = cleft + total_len
                self.debug("extends: <--{}, -->{} / {} ({})".format(cleft, cright, total_len, min_len))
                if total_len >= min_len:
                    candidates.append([cleft, total_len, cright, index - left])
                    self.debug("  -> candidates", candidates)
                    assert(len(candidates) < 10)
                    max_hit_right = max(cright, max_hit_right)
        return candidates


# use this on a constant string (segment)
# pass in the query string to make an alignment
class AlignmentFinder(LoggingClass):

    def __init__(self, characters, indexWordSize = 6, minimumSegmentLength = None):
        LoggingClass.__init__(self)
        #self.setLevelDebug()
        #self.setLevelVerbose()
        assert(characters)
        self.characters = characters
        self._index = {}
        self._indexWordSize = min(len(characters), indexWordSize)
        assert(self._indexWordSize > 0)
        self._minimumSegmentLength = minimumSegmentLength or (self._indexWordSize + 0)
        assert(self._minimumSegmentLength > 0)
        self._build()

    def _build(self):
        index = {}
        wordSize = self._indexWordSize
        seq = self.characters
        n = len(seq)
        for i in range(n - wordSize + 1):
            key = seq[i:(i + wordSize)]
            index.setdefault(key, []).append(i)
        assert(index)
        self._index = index

    def align(self, query):
        # 1) find the longest segments
        # 2) if they don't overlap, use them all
        # 3) if they overlap, score all and use the max
        # 4) iterate on subsegments
        candidates = [ None ]
        best = None
        while candidates:
            candidate = candidates[0]
            if candidate:
                if (not best) or (candidate.score > best.score):
                    best = candidate
                if not _speedup:
                    self.debug("candidate:\n"+ AlignmentDiagram(query, self, candidate).make())
            candidates = self._findBestSegments(query, candidate) + candidates[1:]
        return best

    def _makeAlignment(self, segments, soFar = None, score = None):
        if soFar:
            segs = sorted(soFar.segments + segments, key = lambda s : s.queryIndex)
        else:
            segs = segments
        assert(segs)
        return Alignment(segments = segs, score = score or sum([s.length * s.length for s in segs]))

    def _findBestSegments(self, query, alignmentSoFar):
        gaps = self._findGaps(query, alignmentSoFar)
        if not _speedup:
            self.debug("_fBS", len(alignmentSoFar.segments) if alignmentSoFar else 0, " / gaps =", len(gaps))
        segs = []
        for gap in gaps:
            segs += self._findBestSegmentsInGap(query, alignmentSoFar, gap)
        segs = sorted(segs, key = lambda s: s.length, reverse = True)
        if (not _speedup) and segs:
            self.debug("FBS:\n" + AlignmentDiagram(query, self, self._makeAlignment(segs, score = -1)).make())
        alignments = []
        for seg in segs:
            alignments.append(self._makeAlignment([seg], alignmentSoFar))
        return alignments

    def _findGaps(self, query, a):
        curQuery = 0
        curSegment = 0
        gaps = []
        for aseg in (a.segments if a else []):
            qi, si = aseg.queryIndex, aseg.segmentIndex
            gaps.append(M.Gap(queryStart = curQuery, queryEnd = qi,
                              segmentStart = curSegment, segmentEnd = si))
            curQuery = qi + aseg.length
            curSegment = si + aseg.length
        gaps.append(M.Gap(queryStart = curQuery, queryEnd = len(query),
                          segmentStart = curSegment, segmentEnd = len(self.characters)))
        if not _speedup:
            self.debug("_findGaps", len(gaps), [ "{} / {}".format(g.queryEnd - g.queryStart, g.segmentEnd - g.segmentStart) for g in gaps ])
        minLen = self._minimumSegmentLength
        return [ g for g in gaps if (g.queryEnd - g.queryStart >= minLen) and (g.segmentEnd - g.segmentStart >= minLen) ]

    def _findBestSegmentsInGap(self, query, alignmentSoFar, gap):
        wordSize = self._indexWordSize
        minLen = self._minimumSegmentLength
        checkEvery = max(minLen - wordSize, 1)
        checkSites = list(range(gap.queryStart, gap.queryEnd, checkEvery))
        last = gap.queryEnd - checkEvery
        if last != checkSites[-1]:
            checkSites.append(last)
        candidates = []
        if not _speedup:
            self.debug("_FBSIG", gap.json(skipTypes = True))
        for site in checkSites:
            siteKey = query[site:site + wordSize]
            #self.verbose("CS: {}, {}".format(site, siteKey))
            for segIndex in self._index.get(siteKey, []):
                if (segIndex < gap.segmentStart) or (segIndex + wordSize > gap.segmentEnd):
                    #self.verbose("segment outside of gap, skip")
                    continue
                alreadyHit = False
                for c in candidates:
                    if (site >= c.queryIndex) and (site + wordSize <= c.queryEnd) and \
                       (segIndex >= c.segmentIndex) and (segIndex + wordSize <= c.segmentEnd):
                        alreadyHit = True
                        break
                if alreadyHit:
                    #self.verbose("skipping due to already hit", segIndex - left, segIndex - left + matchLen)
                    continue
                left, right = longestExactMatch(query, (site, wordSize), self.characters, (segIndex, wordSize))
                if not _speedup:
                    self.debug("GOT: {} @ s={}, q={}".format(siteKey, segIndex, site))
                if site - left < gap.queryStart:
                    left = site - gap.queryStart
                if segIndex - left < gap.segmentStart:
                    left = gap.segmentStart - segIndex
                if site + wordSize + right > gap.queryEnd:
                    right = gap.queryEnd - (site + wordSize)
                if segIndex + wordSize + right > gap.segmentEnd:
                    right = gap.segmentEnd - (segIndex + wordSize)
                matchLen = left + right + wordSize
                if not _speedup:
                    self.debug("extends: <--{}, -->{} / {} ({})".format(left, right, matchLen, minLen))
                if matchLen < minLen:
                    continue
                for c in candidates:
                    if (site - left >= c.queryIndex) and (site - left + matchLen <= c.queryEnd) and \
                       (segIndex - left >= c.segmentIndex) and (segIndex - left + matchLen <= c.segmentEnd):
                        alreadyHit = True
                        break
                if alreadyHit:
                    if not _speedup:
                        self.debug("skipping due to already hit", segIndex - left, segIndex - left + matchLen)
                    continue
                newSeg = AlignedSegment(queryIndex = site - left, segmentIndex = segIndex - left, length = matchLen)
                if not _speedup:
                    self.debug(newSeg.json())

                assert(newSeg.queryIndex >= gap.queryStart and newSeg.queryEnd <= gap.queryEnd)
                assert(newSeg.segmentIndex >= gap.segmentStart and newSeg.segmentEnd <= gap.segmentEnd)
                if 1 or _speedup: # TODO: when well-validated...
                    qcheck = query[newSeg.queryIndex:newSeg.queryIndex+matchLen]
                    scheck = self.characters[newSeg.segmentIndex:newSeg.segmentIndex+matchLen]
                    if qcheck != scheck:
                        assert(not [ 1 for qc, sc in zip(qcheck, scheck) if qc != sc and qc != "N" and sc != "N"])
                candidates.append(newSeg)
                if not _speedup:
                    self.debug("  -> candidates", [c.json(skipTypes = True) for c in candidates])
                if len(candidates) > 24:
                    self.debug("high candidate count, aborting", len(candidates), query)
                    return []

        if not _speedup:
            self.debug("_fBSIG returning {} candidates".format(len(candidates)))
        return candidates



class AlignedSegment(M.AlignedSegment):

    @property
    def queryStart(self):
        return self.queryIndex

    @property
    def queryEnd(self):
        return self.queryIndex + self.length

    @property
    def segmentStart(self):
        return self.segmentIndex

    @property
    def segmentEnd(self):
        return self.segmentIndex + self.length


class Alignment(M.Alignment):

    def gaps(self):
        gaps = []
        for idx in range(len(self.segments) - 1):
            left = self.segments[idx]
            right = self.segments[idx + 1]
            gaps.append(M.Gap(queryStart = left.queryEnd,
                              queryEnd = right.queryStart,
                              segmentStart = left.segmentEnd,
                              segmentEnd = right.segmentStart))
        return gaps

    def segmentLength(self):
        return self.segmentEnd() - self.segmentStart()

    def queryLength(self):
        return self.queryEnd() - self.queryStart()

    def queryStart(self):
        return self.segments[0].queryStart

    def queryEnd(self):
        return self.segments[-1].queryEnd

    def segmentStart(self):
        return self.segments[0].segmentStart

    def segmentEnd(self):
        return self.segments[-1].segmentEnd

    def matchSize(self):
        return sum([s.length for s in self.segments])
