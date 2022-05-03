
from . import model as M
from .counters import Counters
from .failures import Failures
from .util import reverse_complement, LoggingClass, writeJsonToPath


class Classifier(LoggingClass):

    def __init__(self, descriptor):
        LoggingClass.__init__(self)
        self.counters = Counters()
        self.descriptor = descriptor

    def categorize(self, r1, r2):
        pass

    def classifyResult(self, fragment, sequenceResult):
        pass

    def saveCounts(self, key):
        fn = "{}_counts.json".format(key)
        self.info("Saving counts to: {}".format(fn))
        writeJsonToPath(self.counters.countsDict(), fn)

    def saveReads(self, key):
        pass

    def getBarcode(self, fr):
        return None
