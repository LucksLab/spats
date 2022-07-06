
from .diagram import FragmentDiagram, SegmentDiagram, ReadsDiagram
from .failures import Failures
from .fragment import FragmentMaker
from .reads import ReadsAnalyzer
from .results import ResultsDB
from .segment import SegmentProcessor
from .util import jsonAtPath, writeJsonToPath, LoggingClass, reverse_complement, M, fs, displaySize, _speedup


COMPARISON_SKIP_SEQUENCE_IDS = {
    's1' : [],
    's2' : [],
    's3' : [],
    's4' : [
        102981, # this is one is a nasty perfect-storm; could get it right, not worth it atm
        146936, # this is fine, too hard to say why...
        182554, 185710, # correct now
        233520, # a subst instead of 3 muts, fine
    ],
    's5' : [ 1038, 1707, 2310 ],
}

# indeterminate diagnostics (xref `indcheck`)
g_countIndeterminate = True
g_indeterminatesOnly = False

class FragmentSetProcessor(LoggingClass):

    def __init__(self, experiment, fragSet):
        LoggingClass.__init__(self, initLogging = True)
        #self.setLevelDebug()
        self.experiment = experiment
        self.fragSet = fragSet
        self._compareSkipSequenceIds = COMPARISON_SKIP_SEQUENCE_IDS.get(self.fragSet.key[:2], [])
        self.maker = FragmentMaker()
        self.maker.r1IsRC = self.experiment.r1IsRC
        self.descriptor = experiment.fragmentDescriptor()
        self.classifier = experiment.classifier(self.descriptor)
        self.processor = SegmentProcessor(self.descriptor)
        self.reportEvery = (1 << 10)
        self.limit = 0
        self.skipCacheHits = False
        self.useFragmentCache = True
        self.barcodePrune = True
        self.barcodePruneSelectFirst = False # this is a major speedup, but seems to produce different #s, haven't investigated why
        self._barcodeInfo = {}
        self._classifyLookup = {}
        self._reads = None
        self.results = None
        self._resultsDb = None
        self._writeResultsEvery = (1 << 18)
        self._diagram = False
        self._allDiagrams = False
        self._diagramMatching = None
        self._diagramBarcode = None
        self._skipCommon = False
        self._keepResults = True
        self._fragmentCache = { }
        if _speedup:
            self.info("*** Speedup is turned on, debugging info will be limited.")

    def autosave(self):
        self._keepResults = True
        self._resultsDb = ResultsDB(self.key)
        self._resultsDb.create()

    @property
    def key(self):
        return self.fragSet.key

    def debugMode(self):
        assert(not _speedup)
        self.setLevelDebug()
        self.maker.setLevelDebug()
        self.processor.setLevelVerbose()

    def _setupReads(self):
        if self._reads:
            return
        self.debug("creating reads analyzer for", self.experiment.name)
        self._reads = ReadsAnalyzer()
        self._reads.addExperiment(self.experiment)

    def _makeResult(self, r1s, r2s, segs = True, reads = False):
        #if r1s.sequenceId == 11722:
        #    self.debugMode()
        #elif r1s.sequenceId > 11722:
        #    STOP
        cls = self.classifier
        ctrs = cls.counters
        ctrs.totalPairs += 1

        if g_countIndeterminate:
            nc1 = r1s.characters.count('N')
            nc2 = r2s.characters.count('N')
            ctrs.incrementKey('R1_N:{}'.format(nc1))
            ctrs.incrementKey('R2_N:{}'.format(nc2))
            ctrs.incrementKey('N:{}'.format(nc1 + nc2))
            if g_indeterminatesOnly:
                return None

        if self.barcodePruneSelectFirst:
            barcode = self.experiment.barcode(r1s, r2s)
            if barcode in self._barcodeInfo:
                self._barcodeInfo[barcode].multiplicity += 1
                ctrs.barcodeSkipped += 1
                return None

        fr = M.FragmentResult(r1 = r1s, r2 = r2s)
        fr.barcode = None # temp hack

        fkey = None
        rkey = self.experiment.r1r2CacheKey(r1s, r2s) if self.useFragmentCache else None

        if rkey and (rkey in self._fragmentCache):
            if self.skipCacheHits:
                return None
            self.debug("r1r2 cache hit")
            ctrs.r1r2CacheHit += 1
            hit = self._fragmentCache[rkey]
            fr.f, fr.s, fr.c, fr.r = hit.f, hit.s, hit.c, hit.r
            if hasattr(fr.c, 'read'):
                fr.c.read.sequenceId = r1s.sequenceId
        elif rkey:
            self.debug("r1r2 cache miss")
            ctrs.uniqueR1R2 += 1

        if not fr.f:
            fr.f = self.maker.make(r1s, r2s)
            fkey = self.experiment.fragmentCacheKey(fr.f) if self.useFragmentCache else None
            if fkey and (fkey in self._fragmentCache):
                if self.skipCacheHits:
                    return None
                ctrs.fragmentCacheHit += 1
                hit = self._fragmentCache[fkey]
                fr.f, fr.s, fr.c, fr.r = hit.f, hit.s, hit.c, hit.r
                if hasattr(fr.c, 'read'):
                    fr.c.read.sequenceId = r1s.sequenceId
                self.debug("fragcache hit", hit.r1.sequenceId)
                fr.s.sequence = fr.f.seq
            elif self.useFragmentCache:
                self.debug("fragcache miss")
                ctrs.uniqueFragment += 1

        if segs:
            if not fr.s:
                fr.s = self.processor.process(fr.f.seq) if fr.f.success else M.SequenceResult(success = False, failure = Failures.fragment)
            if not fr.c:
                fr.c = cls.classify(fr.f, fr.s)

            fr.barcode = cls.getBarcode(fr)
            # pretty sure we don't need/use this counter, and causes counters to be huge
            #ctrs.incrementKey("barcode:{}".format(fr.barcode))
            if self.barcodePrune:
                if self.barcodePruneSelectFirst:
                    assert(fr.barcode not in self._barcodeInfo)
                    self._barcodeInfo[fr.barcode] = fr
                    fr.multiplicity = 1
                else:
                    self._barcodeInfo.setdefault(fr.barcode, []).append(fr)
            else:
                cls.count(fr)

        if reads and (not fr.r) and (fr.f.seq.characters):
            fr.r = self._reads.reads(fr.f.seq.characters)

        if rkey:
            self._fragmentCache[rkey] = fr
        elif fkey:
            self._fragmentCache[fkey] = fr

        if fr.f.reverseOrder:
            ctrs.fragmentReverseOrder += 1

        return fr

    def _makeReads(self, r1s, r2s):
        return self._makeResult(r1s, r2s, segs = False, reads = True)

    def _makeReadsAndResult(self, r1s, r2s):
        return self._makeResult(r1s, r2s, segs = True, reads = True)

    def _compare(self, r1s, r2s):
        cls = self.classifier
        ctrs = cls.counters
        fr = self._makeResult(r1s, r2s)
        read = fr.c.read
        if read.sequenceId < self._comparisonStart:
            ctrs.incrementKey("compare:skipped")
            return fr
        if read.sequenceId in self._compareSkipSequenceIds:
            self.debugf("skipping sequenceId @{}..", read.sequenceId)
            ctrs.incrementKey("compare:skipped")
            return fr
        comp = self._comparisonReads[read.sequenceId]
        self._comparisonReads[read.sequenceId] = None
        assert(comp.sequenceId == read.sequenceId)
        if comp.success != read.success:
            if read.success:
                if read.overlap > 16:
                    ctrs.incrementKey("compare:overlapOverrule")
                    ctrs.incrementKey("compare:overlapOverrule:{}".format(read.overlap))
                    self.debugf("ignoring fail -> success @{} due to large overlap ({})", read.sequenceId, read.overlap)
                    return fr
                lgInserts = sum([ m.size for m in read.mutations if m.mutType == fs.insertion])
                if lgInserts >= 30:
                    ctrs.incrementKey("compare:lgInsert:{}".format(lgInserts))
                    self.debugf("ignoring fail -> success @{} due to large insertion ({})", read.sequenceId, lgInserts)
                    return fr
                lgSubs = [ m for m in read.mutations if m.mutType == fs.substitution and m.size > 30 ]
                if lgSubs:
                    ctrs.incrementKey("compare:lgSub:{}".format(lgSubs[0].size))
                    self.debugf("ignoring fail -> success @{} due to large insertion ({})", read.sequenceId, lgSubs[0].json())
                    return fr
                edgeMuts = [ m for m in read.mutations if m.site <= 1 ]
                if edgeMuts:
                    ctrs.incrementKey("compare:edgeMuts")
                    self.debugf("ignoring fail -> success @{} due to edge muts ({})", read.sequenceId, [ m.json() for m in edgeMuts ])
                    return fr
                if read.mutations and [ m for m in read.mutations if (m.mutType == fs.mutation and (m.site in [0, 1, 17, 18, 45, 46, 63, 64])) ]:
                    ctrs.incrementKey("compare:r1r2EdgeMuts")
                    self.debugf("ignoring fail -> success @{} due to could not complete R1/R2 edge issue ({})", read.sequenceId, read.json())
                    return fr
                if fr.f.qualityResolved and [ q for q in fr.f.qualityResolved if q in [ 0, 1, fr.f.overlap - 2, fr.f.overlap - 1 ] ]:
                    ctrs.incrementKey("compare:r1r2qualityResolvedEdgeMuts")
                    self.debugf("ignoring fail -> success @{} due to quality resolved could not complete R1/R2 edge issue ({})", read.sequenceId, fr.f.qualityResolved)
                    return fr
                fp, rp = fr.s.segmentWithKey("fPrimer"), fr.s.segmentWithKey("rPrimer")
                if fp.errorCount() or rp.errorCount():
                    ctrs.incrementKey("compare:primer")
                    self.debugf("ignoring fail -> success @{} due to primer errors ({} / {})", read.sequenceId, fp.errors, rp.errors)
                    return fr
            else:
                if fr.f.failure == Failures.indeterminate:
                    # we could address these; look overlap ignoring Ns
                    ctrs.incrementKey("compare:indeterminateFragment")
                    self.debugf("ignoring success -> fail @{} due to indeterminate fragment", read.sequenceId)
                    return fr
                if (fr.s.failure == Failures.minEnd) or (fr.s.failure == Failures.maxStart):
                    ctrs.incrementKey("compare:targetEdge")
                    self.warnf("ignoring success -> fail @{} due to previous lack of checking for target edges", read.sequenceId)
                    return fr
            raise Exception("Primary read comparison error @{}: {} != {}".format(read.sequenceId, read.json(), comp.json()))
        if read.success:
            def dsz(m):
                if m.mutType == fs.mutation:     return 0
                if m.mutType == fs.insertion:    return m.size
                if m.mutType == fs.deletion:     return 0 - m.size
                if m.mutType == fs.substitution: return m.size - m.substSize
                ASSERT_NOT_REACHED
            overlapDelta = sum([ dsz(m) for m in read.mutations] )
            if comp.overlap != read.overlap:
                if comp.overlap + overlapDelta != read.overlap:
                    ctrs.incrementKey("compare:overlapMismatch")
                    if (not comp.mutations):
                        raise Exception("Overlap mismatch @{}: {} != {} / {}".format(read.sequenceId, read.overlap, comp.overlap, overlapDelta))
                    else:
                        self.warnf("Overlap mismatch @{}: {} != {} / {}", read.sequenceId, read.overlap, comp.overlap, overlapDelta)

            if len(comp.mutations or []) != len(read.mutations or []):
                cmuts = []
                delSizes = {}
                insSizes = {}
                for cm in comp.mutations:
                    if cm.mutType == 1376: # mut
                        if (fr.f.qualityResolved or fr.f.errors) and (len(r1s.characters) - (5 + 23 + cm.site)) <= read.overlap and ((len(r1s.characters) - (5 + 23 + cm.site) - overlapDelta) >= 0):
                            self.debugf("ignoring compared mut @{} that seems to be on resolved overlap  {}", read.sequenceId, cm.json())
                            ctrs.incrementKey("compare:resolvedOverlap")
                            continue
                        if [ x for x in comp.mutations if (x.site == cm.site - 1) and (x.mutType == 1376) ]:
                            self.debugf("ignoring compared mut @{} that seems to be a mut dup {}", read.sequenceId, cm.json())
                            ctrs.incrementKey("compare:dupMut")
                            continue
                    if cm.mutType == 1374: # del
                        if cm.size in delSizes:
                            ctrs.incrementKey("compare:dupDel")
                            self.debugf("ignoring compared mut @{} that seems to be a deletion dup  {}", read.sequenceId, cm.json())
                            continue
                        delSizes[cm.size] = cm
                    if cm.mutType == 1375: # ins
                        if cm.size in insSizes:
                            ctrs.incrementKey("compare:dupIns")
                            self.debugf("ignoring compared mut @{} that seems to be an insertion dup  {}", read.sequenceId, cm.json())
                            continue
                        insSizes[cm.size] = cm
                    cmuts.append(cm)
                if len(cmuts) != len(read.mutations):
                    rmuts = []
                    for rm in read.mutations:
                        if rm.mutType == fs.mutation and (len(r1s.characters) - (5 + 23 + rm.site)) <= read.overlap and ((len(r1s.characters) - (5 + 23 + rm.site)) >= 0):
                            self.debugf("ignoring read mut @{} that seems to be on resolved overlap  {}", read.sequenceId, rm.json())
                            ctrs.incrementKey("compare:resolvedNewOverlap")
                            continue
                        rmuts.append(rm)
                    if len(rmuts) != len(cmuts):
                        raise Exception("Mut count failure @{}: {} != {}".format(read.sequenceId, read.json(), comp.json()))
                comp.mutations = cmuts
            if comp.mutations:
                for rm, cm in zip(sorted(read.mutations, key = lambda m : m.site), sorted(comp.mutations, key = lambda m : m.site)):
                    cm.mutType = self._mutTranslator[cm.mutType]
                    if rm.mutType != cm.mutType or rm.size != cm.size:
                        raise Exception("Mut basic failure @{}: {} != {}".format(read.sequenceId, rm.json(), cm.json()))
                    if abs(rm.site - cm.site) > 3 or (rm.mutType == fs.mutation and abs(rm.site - cm.site) > 1):
                        raise Exception("Mut site failure @{}: {} != {}".format(read.sequenceId, rm.json(), cm.json()))
        return fr

    def _doRun(self, handler):
        appxTotal = self.fragSet.approximatePairCount()
        self.info("Running {}, ~{} fragments".format(self.key, displaySize(appxTotal)))
        reportEvery = self.reportEvery
        if appxTotal / 50 > reportEvery:
            reportEvery = int(appxTotal / 50)
        count = 0
        results = []
        shownKeys = set()
        showBarcode = self._diagramBarcode
        for r1s, r2s in self.fragSet.iterator():
            count += 1
            if 0 == (count % reportEvery):
                self.info("  {}@{}  [{:.2f}%, {}]".format(self.key, count, 100 * count / max(count, appxTotal), len(self._fragmentCache)))
            if showBarcode:
                barcode = r1s.characters[:5] + "_" + reverse_complement(r2s.characters[:5])
                if barcode != showBarcode:
                    continue
            try:
                if self._diagramMatching:
                    self.classifier.counters.reset()
                res = handler(r1s, r2s)
            except:
                self.warn("Error {}/@{}".format(self.key, r1s.sequenceId))
                raise
            if not res:
                assert(self.skipCacheHits or self.barcodePruneSelectFirst or g_indeterminatesOnly)
                continue

            #if r1s.sequenceId == 11722:
            #    self._allDiagrams = self._diagram = True

            if self._diagram or self._diagramMatching:
                show = self._allDiagrams or showBarcode
                if (not show) and (res.s):
                    # TODO: res.description?
                    resKey = "{}:{}:{}:{}:{}::{}".format(res.f.success,
                                                         res.f.failure,
                                                         res.s.success,
                                                         len(res.s.segments),
                                                         res.s.failure,
                                                         res.s.failedSegment)
                    show = resKey not in shownKeys
                    shownKeys.add(resKey)
                if show and self._skipCommon and ((res.c.category == "match" and res.c.description == "full") or (res.c.category == "barcodeSkipped")):
                    show = False
                if self._diagramMatching:
                    if self.barcodePrune and res.s and res.c:
                        # otherwise we have to wait for the counting at the end
                        self.classifier.count(res)
                    show = (self.classifier.counters.getKey(self._diagramMatching) > 0)
                    if (not show) and self._reads:
                        show = self._reads.tagMatch(self._diagramMatching, res.r)
                if show:
                    print("\n" + self.diagram(res))
                    print(" ".join([ f"{k}={v}" for k, v in self.classifier.counters._counts.items() if ":" not in k ]))
                    input("...")
            if self._keepResults:
                results.append(res)
            if self._resultsDb and len(results) >= self._writeResultsEvery:
                self._resultsDb.write(results)
                results = []
            if self.limit and (len(results) >= self.limit):
                self.info("Limited to {} results".format(self.limit))
                break

        self._count = count
        if self.barcodePrune:
            self._countBarcodeResults()
        self.maker.alignmentHitReport()
        self.results = results

    def _countBarcodeResults(self):
        cls = self.classifier
        ctrs = cls.counters
        topBarcode = 0
        for barcode, fresults in self._barcodeInfo.items():
            ctrs.uniqueBarcode += 1
            if self.barcodePruneSelectFirst:
                cls.count(fresults)
                topBarcode = max(topBarcode, fresults.multiplicity)
                continue
            passing = [ fr for fr in fresults if fr.c and fr.c.success ]
            if not passing:
                cls.count(fresults[0])
                continue
            topBarcode = max(topBarcode, len(passing))
            if len(passing) > 1:
                ctrs.incrementKey("barcodeSkipped", multiplicity = len(passing) - 1)
            cls.count(passing[0])
        if topBarcode:
            ctrs.setKey("topBarcode", topBarcode)

    def reads(self):
        self._setupReads()
        self._keepResults = False
        self._doRun(self._makeReads)
        total = self._count
        imap = self._reads.hits()
        self.info("Reads tags:")
        for k in sorted(imap.keys(), key = lambda k : 0 - imap[k]):
            self.info("  {}/{}: {}  ({:.2f}%)".format(self.key, k, imap[k], 100 * imap[k] / total))

    def run(self):
        self._doRun(self._makeResult)

    def readsAndRun(self):
        self._setupReads()
        self._doRun(self._makeReadsAndResult)

    def compare(self, compareKey, start = 0):
        self._mutTranslator = {
            1374: 1386,  # deletion
            1375: 1385,  # insertion
            1376: 1384,  # mutation
        }
        self._comparisonStart = start
        self._comparisonReads = [ M.NGSRead().loadFromJson(r) for r in jsonAtPath("{}_reads.json".format(compareKey))['reads'] ]
        try:
            self._doRun(self._compare)
        finally:
            cd = self.classifier.counters.countsDict()
            compares = { k[8:] : cd[k] for k in sorted(cd.keys(), key = lambda k : 0 - cd[k]) if k.startswith('compare:') }
            self.info("Comparison counts:", compares)

    def saveResults(self):
        self.classifier.saveCounts(self.key)
        if self.results:
            if self._resultsDb:
                self._resultsDb.write(self.results)
                self._resultsDb.finish()
            else:
                ResultsDB(self.key).save(self.results, self._reads)
        if 0:
            fn = "{}_results.json".format(self.key)
            self.info("Saving results to: {}".format(fn))
            trimmed = []
            for r in self.results:
                info = r.json(skipTypes = True)
                if "s" in info:
                    for s in info["s"]["segments"]:
                        del s["segment"]
                trimmed.append(info)
            writeJsonToPath({ "results" : trimmed }, fn)
        
    def resultWithId(self, sequenceId, resultOnly = False):
        if self.results:
            return self.results[sequenceId]
        if not resultOnly:
            self._setupReads()
        r1s, r2s = self.fragSet.pairWithId(sequenceId)
        return self._makeResult(r1s, r2s) if resultOnly else self._makeReadsAndResult(r1s, r2s)

    def resultForPair(self, r1s, r2s):
        self._setupReads()
        return self._makeReadsAndResult(r1s, r2s)

    def diagram(self, res):
        assert(isinstance(res, M.FragmentResult))
        d = ""
        if res.f:
            d = "{}\n".format(FragmentDiagram(res.f).make())
        if res.r:
            d += "\n{}\n".format(ReadsDiagram(res.f.seq, res.r).make())
        if res.s and ((not res.f) or res.f.success):
            d += "\n{}\n".format(SegmentDiagram(res.s, fragment = res.f, classifyResult = res.c).make())
        return d
