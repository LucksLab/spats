

class Experiment(object):

    def __init__(self, config = None):
        self._config = config or {}

    def configure(self, config):
        # we want to return new instances if there is different configuration
        return self.__class__(config = config)

    @property
    def name(self):
        return self.__class__.__name__


    # required overrides

    def fragmentDescriptor(self):
        pass

    def classifier(self, descriptor):
        pass


    # optional overrides

    def aliases(self):
        return []

    @property
    def r1IsRC(self):
        return False

    def barcode(self, r1s, r2s):
        # return the barcode; only used if barcode pruning. return
        # what it would be on a match; it'll only be a prune skip if
        # there's a matching one already in the cache.
        return None

    def r1r2CacheKey(self, r1s, r2s):
        # R1/R2 cache keys are a bit risky, because the same R1/R2
        # might produce different fragments with different
        # quality. you can turn it on, but you might get slightly
        # different results.
        return None

    def fragmentCacheKey(self, frag):
        return frag.seq.characters

    def addExtraReads(self, reads):
        pass
