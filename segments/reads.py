
from .index import Index
from .util import LoggingClass, M, reverse_complement, flatten


class ReadsAnalyzer(LoggingClass):

    def __init__(self):
        LoggingClass.__init__(self)
        #self.setLevelDebug()
        self._items = []

    def reset(self):
        for i in self._items:
            i.hits = 0

    def addExperiment(self, exp, includeReverseComplements = True, indexedOnly = False):
        self.debug("ae", exp, exp.fragmentDescriptor())
        rcs = []
        usedAlready = set()
        for sd in exp.fragmentDescriptor().segments:
            if sd.wildcard or (indexedOnly and (not sd.index)):
                continue
            if sd.characters in usedAlready:
                continue
            self.add(sd.key, sd.characters)
            rcs.append(sd)
            usedAlready.add(sd.characters)
        if not includeReverseComplements:
            return
        for sd in rcs:
            self.add("{}_rc".format(sd.key), reverse_complement(sd.characters))
        exp.addExtraReads(self)

    def add(self, key, characters, minimumWordLength = 10):
        item = M.ReadsItem(key = key, characters = characters)
        self.debug("building index for {}, len={} [{}]".format(key, len(characters), characters[:64]))
        item.index = Index(characters)
        if minimumWordLength:
            item.index._minimum_length = minimumWordLength
        item.index.build()
        self._items.append(item)

    def keys(self):
        return [ i.key for i in self._items ]

    def hits(self):
        return { i.key : i.hits for i in self._items }

    def reads(self, query):
        self.debug("reads query:", query)
        items = flatten([ self._findItem(i, query) for i in self._items ])
        assert(len(items) < 20)
        return [ i for i in items if i ]

    def _findItem(self, item, query):
        candidates = item.index.find_all_partials(query)
        res = []
        for candidate in candidates:
            queryStart, matchLen, queryEnd, itemStart = candidate
            item.hits += 1
            res.append(M.ReadsResult(item = item,
                                     itemStart = itemStart,
                                     queryStart = queryStart,
                                     queryEnd = queryEnd))
        return res

    def tagMatch(self, tagQuery, reads):
        rdTags = set([r.item.key for  r in reads ])
        for t in tagQuery.split(":"):
            if t.startswith("-"):
                if t[1:] in rdTags:
                    return False
            else:
                if t not in rdTags:
                    return False
        return True
