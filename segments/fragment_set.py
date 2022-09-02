import os

from . import model as M
from .logging import LoggingClass
from .parse import FastFastqParser
from .sequence import Sequence



class FragmentSet(LoggingClass):

    def __init__(self, key, path = None, r1 = None, r2 = None, singleFragment = None):
        LoggingClass.__init__(self)
        #self.setLevelDebug()
        self.key = key
        if singleFragment:
            assert(isinstance(M.Fragment, singleFragment))
            self._singleFragment = singleFragment
        else:
            self._singleFragment = None
            self.r1Path = r1 or os.path.join(path or ".", "{}_r1.fastq".format(key))
            self.r2Path = r2 or os.path.join(path or ".", "{}_r2.fastq".format(key))
            if not os.path.exists(self.r1Path):
                self.r1Path = self.r1Path.replace('_r1.', '_R1.')
            if not os.path.exists(self.r2Path):
                self.r2Path = self.r2Path.replace('_r2.', '_R2.')
            if (not os.path.exists(self.r1Path)) or (not os.path.exists(self.r2Path)):
                raise Exception("Could not locate R1/R2 at {}/{}".format(self.r1Path, self.r2Path))
            self.batchSize = 131072
            self._curBatch = None

    def _makePair(self, item, sequenceId):
        _, r1, r2, ident, r1q, r2q = item
        r1s = M.Sequence(sequenceId = sequenceId,
                         name = "{}[{}]/R1".format(self.key, sequenceId),
                         identifier = ident,
                         characters = r1,
                         quality = r1q)
        r2s = M.Sequence(sequenceId = sequenceId,
                         name = "{}[{}]/R2".format(self.key, sequenceId),
                         identifier = ident,
                         characters = r2,
                         quality = r2q)
        return (r1s, r2s,)

    def _parser(self):
        assert(not self._singleFragment)
        return FastFastqParser(self.r1Path, self.r2Path, parse_quality = True)

    def iterator(self):
        if self._singleFragment:
            yield self.singlePair()
            return
        sequenceId = 0
        with self._parser() as p:
            for batch in p.iterator(self.batchSize):
                for item in batch:
                    yield self._makePair(item, sequenceId)
                    sequenceId += 1

    def pairWithId(self, sequenceId):
        assert(sequenceId >= 0)
        if self._singleFragment:
            assert(0 == sequenceId)
            return self.singlePair()

        if self._curBatch:
            # caching to help regressions
            if (sequenceId >= self._batchStart) and (sequenceId < self._batchStart + len(self._curBatch)):
                return self._makePair(self._curBatch[sequenceId - self._batchStart], sequenceId)

        rem = sequenceId
        batchStart = 0
        with self._parser() as p:
            for batch in p.iterator(self.batchSize):
                if rem <= len(batch):
                    self._curBatch = batch
                    self._batchStart = batchStart
                    return self._makePair(batch[rem], sequenceId)
                rem -= len(batch)
                batchStart += len(batch)
            raise Exception("SequenceId not found: {}".format(sequenceId))

    def singlePair(self):
        return (self._singleFragment.r1, self._singleFragment.r2,)

    def approximatePairCount(self):
        return 1 if self._singleFragment else self._parser().appx_number_of_pairs()
