
from . import model as M
from .util import reverse_complement


class Descriptor(object):

    def __init__(self):
        pass

    def makeSegd(self, key, nts = None, wildcard = False, fixedLength = 0, **kwargs):
        segd = M.SegmentDescriptor(key = key, characters = nts, wildcard = wildcard, **kwargs)
        if wildcard:
            if fixedLength:
                assert("minLength" not in kwargs)
                assert("maxLength" not in kwargs)
                segd.minLength = segd.maxLength = fixedLength
        else:
            assert(nts)
            segd.maxLength = len(nts)
            if fixedLength:
                assert("minEnd" not in kwargs)
                assert("minLength" not in kwargs)
                segd.minEnd = segd.minLength = len(nts)
            else:
                segd.maxStart = kwargs.get("maxStart", segd.maxLength)
        return segd
