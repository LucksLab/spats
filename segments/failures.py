
class Failures(object):
    
    _failures = {

        # fragment
        "fragment" : "failed to create fragment",
        "indelsInFragment" : "fragment overlap has indels",
        "noOverlap" : "no overlap",
        "contained" : "one of R1/R2 is contained in the other",
        "overlapErrors" : "overlap alignment failed",
        "indeterminate" : "indeterminate sequence",

        # segment
        "nomatch" : "no match",
        "matchErrors" : "match errors failure",
        "minSegmentLength" : "segment below minimum length",
        "maxSegmentLength" : "segment above maximum length",
        "gapTooBig" : "gap too big to fill",
        "gapTooSmall" : "gap too small to fill",
        "validation" : "validation failure",
        "outOfOrder" : "segments found out-of-order",
        "gapAlignment" : "could not align a gap between two segments",
        "shortMatch" : "match was too short to extend",
        "minEnd" : "segment could not be aligned to right edge",
        "maxStart" : "segment could not be aligned to left edge",
        "edgeIndel" : "ambiguous indel on the edge of a segment",
        "wildcardIndel" : "invalid indel on a wildcard segment",
        "interestingQueryOverlap" : "aligned segments overlap",

        # spats
        "mask" : "mask failure",
        "prefix" : "prefix failure",

    }

    @staticmethod
    def _setup():
        for key, value in Failures._failures.items():
            setattr(Failures, key, key)

    @staticmethod
    def allFailureKeys():
        return Failures._failures.keys()

    @staticmethod
    def description(failure):
        return Failures._failures[failure]


Failures._setup()
