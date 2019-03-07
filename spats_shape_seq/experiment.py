# large number of objects is intentional: conceptual clarity comes
# from having each object do its own task. some of these may
# eventually collapse and/or go away, but useful at least now for
# communication.

# a sequence represents a fixed string of ACGT used for experimental
# processing. includes targets, adapters, linkers, dumbbells, etc..
class Sequence(self):

    TARGET    = "target"
    ADAPTERT  = "adapter_t"
    ADAPTERB  = "adapter_b"
    LINKER    = "linker"
    DUMBBELL  = "dumbbell"
    
    def __init__(self, seq, type):
        assert(type in [ Sequence.TARGET,
                         Sequence.ADAPTERT,
                         Sequence.ADAPTERB,
                         Sequence.LINKER,
                         Sequence.DUMBBELL,
                         ])
        self.sequence = seq
        self.type = type

    @property
    def string(self):
        return self.sequence


class Target(Sequence): pass


# core idea: things should be based on fragments, not pairs. a
# fragment is a straightforward string of chars, which can be
# converted to (and ideally from) an R1/R2 pair. it's much easier to
# describe how fragments can be built up than it is to describe valid
# values of R1/R2.
#
# fragments are broken up into sections, see below. for now, we are
# dealing only with "perfect fragments": fragments that are valid in
# the context of the experiment, with no mutations or indels.
#
# for example:
#    TGAACAGCGACTAGGCTCTTCA GCGAGAGTAGGGAACTGCCAGGCATC CTGACTCGGGCACCAAGGAC AGATCGG
#    [  dumbbell          ] [       target           ] [    linker        ] [adapter]
# this is a fragment with 4 sections. 1st and 3rd are exact matches, 
# the 2nd is a substring match, and the last is a prefix match.
#
# note that it is easy to read stop data directly off of a fragment
# (look for the target section and find its sequence_start).

class Fragment(object):

    def __init__(self):
        self._sections = []
        self._keyedSections = {}

    @property
    def isTransformed(self):
        return False

    def appendSection(self, sec):
        self._sections.append(sec)
        self._keyedSections[sec.type] = sec

    @property
    def string(self):
        return "".join([s.string for s in self._sections ])

    def clone(self):
        f = Fragment()
        for s in self._sections:
            f.appendSection(f)
        return f

    def pair(self, experiment):
        return Pair.PairFromFragment(self, experiment)


# a fragment section defines a known part of a fragment, in reference
# to one or more of the base sequences.
class FragmentSection(object):

    def __init__(self, seq, start, length):
        self.sequence = seq
        self.sequence_start = start
        self.length = length

    @property
    def sequence_end(self):
        return self.sequence_start + length

    @property
    def string(self):
        return self.sequence.string[self.sequence_start:self.sequence_end]


# like usual, showing how to go from a fragment to a pair
class Pair(object):

    def __init__(self, r1, r2):
        self.r1 = r1
        self.r2 = r2

    def PairFromFragment(fragment, experiment):
        fragmentString = fragment.string
        r2 = fragmentString[:self.experiment.r2Length]
        r1 = reverse_complement(fragmentString[-self.experiment.r1Length:])
        return Pair(r1, r2)


# core idea: we can use data-driven descriptors of the expected
# subsections of a fragment, instead of handling each possible case in
# code.
class SectionMatcher(object):

    MATCH_EXACT      = "exact"
    MATCH_SUBSTRING  = "substring"
    MATCH_PREFIX     = "prefix"
    MATCH_SUFFIX     = "suffix"

    def __init__(self, seq, required = True, matchType = SectionMatcher.MATCH_EXACT):
        assert(matchType in [ SectionMatcher.MATCH_EXACT, SectionMatcher.MATCH_SUBSTRING, SectionMatcher.MATCH_PREFIX, SectionMatcher.MATCH_SUFFIX ])
        self.sequence = seq
        self.required = required
        self.matchType = matchType

    @property
    def sequenceType(self):
        return self.sequence.type


# similar to run.py, describes the (user-configurable) metadata for an experiment
class Experiment(object):

    def __init__(self):
        self.cotrans = False
        self.cotrans_linker = Sequence("...")
        self.use_dumbbell = False
        self.dumbbell = Sequence("...")
        self.adapter_t = Sequence("...")
        self.adapter_b = Sequence("...")
        self.target = Target("...")

        self.r1Length = 32
        self.r2Length = 36
        ...

    def descriptor(self):
        # core idea: here's where we build our data-driven description
        # of what fragments in this experiment look like
        maxLen = len(target)
        if self.cotrans: etc...
        fd = FragmentDescriptor(self.pairLength, maxLen)
        if self.use_dumbbell:
            fd.addMatcher(SectionMatcher(self.dumbbell, True, SectionMatcher.MATCH_EXACT))
        if self.cotrans:
            fd.addMatcher(SectionMatcher(self.target, True, SectionMatcher.MATCH_SUBSTRING))
            fd.addMatcher(SectionMatcher(self.cotrans_linker, True, SectionMatcher.MATCH_EXACT))
        else:
            fd.addMatcher(SectionMatcher(self.target, True, SectionMatcher.MATCH_SUFFIX))
        fd.addMatcher(SectionMatcher(self.adapter_t, False, SectionMatcher.MATCH_PREFIX))
        return fd

    # sanity check on an experiment to make sure that perfect matches are not ambiguous
    def validateExperiment(self):
        fd = self.descriptor()
        frags = fd.perfectFragments()
        fragMap = {}
        for f in frags:
            p = f.pair()
            key = p.r1 + p.r2
            if key in fragMap:
                raise Exception("two fragments produce same r1/r2 pair: {}, {}, {}, {}", f, fragMap[key], p.r1, p.r2)
            fragMap[key] = f
        return # ok, no dups, passes validation


# core idea: this contains the metadata that describe the possible
# fragments in an experiment. can be used to generate all possible
# perfect fragments and, therefore, R1/R2 pairs -- in other words, the
# basis for the lookup algorithm.
class FragmentDescriptor(object):

    def __init__(self, minLen, maxLen):
        self.minLength = minLen
        self.maxLength = maxLen
        self._matchers = []
        self._keyedMatchers = {}

    @property
    def matchers(self):
        return self._matchers

    def addMatcher(self, m):
        self._matchers.append(m)
        self._keyedMatchers[m.sequenceType] = m

    def perfectFragments(self, experiment):
        requiredLength = 0
        for m in self._matchers:
            if m.required:
                if m.type == SectionMatcher.MATCH_EXACT:
                    requiredLength += len(m.sequence.string)
        frags = []
        self._recursivePerfectFrags(self._matchers, Fragment(), requiredLength, frags)
        return frags

    def _recursivePerfectFrags(self, matchersLeft, workingFrag, requiredLength, fragsSoFar):
        TODO
        if 0 == len(matchersLeft):
            return
        m = matchersLeft[0]
        if not m.required:
            _recursivePerfectFrags(matchersLeft[1:], workingFrag, requiredLength, fragsSoFar)
        if m.type == SectionMatcher.MATCH_EXACT:
            workingFrag.addsection #...
        else:
            # get the possible ranges based on whether it's substring/prefix/suffix


# a tranformation represents a change from a perfect fragment. we
# handle mutations, inserts, and deletes.
class Transformation(object):

    XFRM_MUTATION  = "mut"
    XFRM_INSERT    = "ins"
    XFRM_DELETE    = "del"

    def __init__(self, type, location, parameter):
        self.type = type
        self.location = location
        self.parameter = parameter

    def applyTo(self, string):
        return getattr(self, "_apply_{}".format(self.type))(string)

    def _apply_mut(self, string):
        param = self.parameter
        loc = self.location
        assert(1 == len(param))
        if string[loc] == param:
            raise Exception("requested mutation does not change result")
        return string[:loc] + param + string[loc + 1:]

    def _apply_ins(self, string):
        param = self.parameter
        loc = self.location
        assert(0 < len(param))
        return string[:loc] + param + string[loc:]

    def _apply_del(self, string):
        param = self.parameter
        loc = self.location
        assert(0 < param)
        if loc + param >= len(string):
            raise Exception("requested delete too large")
        return string[:loc] + string[loc + param:]


# a TransformedFragment is a wrapper around a fragment, which includes
# information about the mutations applied to get to a target
# string. note that this refers to the original (perfect) fragment,
# which means it's easy to read stop data off of this; and mut data is
# clear from the transformations.
class TransformedFragment(object):

    def __init__(self, frag):
        self.fragment = frag
        self._transformations = []

    @property
    def isTransformed(self):
        return bool(self._transformations)

    def addTransformation(self, xfrm):
        self._transformations.append(xfrm)

    @property
    def string(self):
        string = self.fragment.string
        for xfrm in self._transformations:
            string = xfrm.applyTo(string)
        return string

    def pair(self, experiment):
        return Pair.PairFromFragment(self, experiment)

    def mutations(self):
        # TODO: this isn't quite correct, since an insert can affect downstream locations
        # muts need to be relative to the target
        # but it gives the idea...
        return [ xfrm.location for xfrm in self._transformations ]


# the idea of this is that we can generate (many of) the possible
# transformed fragments in an experiment. for some parameters we can
# use this as the basis of an index to look up results with mutations
# (lookup algorithm). 
class TransformedFragmentGenerator(object):

    def __init__(self, experiment):
        self.experiment = experiment
        self.maxInsertLength = 2
        self.maxDeleteLength = 4
        self.maxTransformations = 1
        self._cachedPossibles = None


    # what kind of transformations are biologically likely/interesting?
    #  - for example: can a site be mutated twice? is it possible/likely experimentally?
    #  - probably we're interested in: "some shit happened at site N", not exactly what that was
    #    - ie, whether it's an insert of 8NT at site 7, or a mut at site 7, doesn't affect stats

    def generate(self, includePerfect = True):
        assert(1 == self.maxTransformations) #TBD

        xfrms = []

        for frag in baseFrags:

            if includePerfect:
                xfrms.append(TransformedFragment(frag))  # original frag

            # biologically, do we care about mutations in:
            #  - adapter, dumbbell, linker, etc
            #  - my sense is "no", and those should be thrown out; but needs discussion
            #    - and, we should be able to answer questions like: "how pairs had muts in the adapter?"
            # for now, this assumes we only care about muts in the target
            fragString = frag.string
            targetInFragRange = frag.getTargetRange()

            for idx in targetInFragRange:
                cur = fragString[idx]
                for ch in [ "A", "C", "G", "T" ]:
                    if ch == cur:
                        continue
                xf = TransformedFragment(frag)
                xf.addTransformation(Transformation(Transformation.MUTATION, idx, ch))
                xfrms.append(xf)

                for insert in self._possibleInsertions():
                    xf = TransformedFragment(frag)
                    xf.addTransformation(Transformation(Transformation.INSERT, idx, insert))
                    xfrms.append(xf)

                for delLen in range(1, self.maxDeleteLength + 1):
                    if idx + delLen >= targetInFragRange.stop:
                        continue
                    xf = Transformation(frag)
                    xf.addTransformation(Transformation(Transformation.DELETE, idx, delLen))
                    xfrms.append(xf)

        return xfrms

    def _possibleInsertions(self):
        if not self._cachedPossibles:
            self._cachedPossibles = self._possibleStringsOfLength(self.maxInsertLength)
        return self._cachedPossibles

    def _possibleStringsOfLength(self, length):
        if 1 == length:
            return [ "A", "C", "G", "T" ]
        strs = []
        for s in self._possibleStringsOfLength(length - 1):
            for ch in self._possibleStringsOfLength(1):
                strs.append(s + ch)
        return strs


# core idea: we can also use this to detect where ambiguities are
# inherent in the experiment with mutations. for example, we know an
# edge toggle at stop n will be indistinguishable from an insert with
# stop n-1...but there will be many other things like this. using this
# tool can point out all the cases where amiguity are present. then
# perhaps we could use biology to determine which of the alternatives
# is likely, or assign them probabilities?
def experimentAmbiguity(exp):
    fd = exp.descriptor()
    tfg = TransformedFragmentGenerator(exp)
    xfrms = tfg.generate()

    keyMap = {}
    for xf in xfrms:
        pair = xf.pair(exp)
        key = pair.r1 + pair.r2
        if key in keyMap:
            existing = keyMap[key]
            if not existing.isTransformed:
                print("Transformed {} is ambiguous with perfect {}".format(xf, existing))
            else:
                print("Transformations {} and {} are ambiguous".format(xf, existing))
        else:
            keyMap[key] = xf


class PairResult(object):

    def __init__(self, experiment, pair, fragment):
        self.experiment = experiment
        self.pair = pair
        self.fragment = fragment

    @property
    def matched(self):
        return bool(self.fragment)

    @property
    def stop(self):
        NYI # read it off the fragment

    @property
    def mutations(self):
        if not self.fragment.isTransformed:
            return None
        

# core idea: create an index of all valid (perfect) fragments, and for
# a given r1/r2 pair, determine the one which is the best hit. after
# that, use S/W to figure out a likely set of transformations.
class IndexedPairProcessor(object):

    def __init__(self, experiment):
        self.experiment = experiment
        self.index()

    def _keyForPair(self, pair):
        # can probably do a lot better than this, based on what subset of this
        # is unique in the experiment
        return pair.r1 + pair.r2

    def index():
        index = {}
        fd = self.experiment.descriptor()
        for frag in fd.perfectFragments():
            index[self._keyForPair(frag.pair(self.experiment))] = frag
        self._index = index

    def processPair(self, pair):
        key = self._keyForPair(pair)
        exact = self._index[key]
        if exact:
            return PairResult(self.experiment, pair, exact)

        # no perfect match, so find MEM for r1/r2 in the index

        # then use S/W to create a TransformedFragment off of the base Fragment with best match

        # then make a result with it
        return PairResult(self.experiment, pair, xfrmed)


# core idea: use the metadata of the SectionMatchers to determine the
# most likely match and transformations for a given r1/r2 pair.
class DescriptorPairProcessor(object):

    def __init__(self, experiment):
        self.experiment = experiment
        self.descriptor = experiment.descriptor()
        self.index()

    def index(self):
        indices = {}
        for matcher in self.descriptor.matchers:
            index = None # TODO: build an index for the sequence referred to by that matcher, similar to find_partial
            indices[matcher.sequenceType] = index
        self._indices = indices

    def processPair(self, pair):
        remaining = pair.r2
        usedR1 = False
        # TODO: if r1 and r2 have (statistically) significant overlap, then assume they are joined and set remaining to the full fragment
        # otherwise, assume there is a gap (within the target) and leave usedR1=False

        # first, let each section match independently, searching for
        # maximal matches (at most two, e.g., if there's a mut in the
        # middle). taking into account match type
        # (exact/substring/prefix/suffix).
        matches = []
        for matcher in self.descriptor.matchers:
            index = self._indices(matcher.sequenceType)
            res = find_partial(remaining, index)
            if not res and matcher.required:
                return PairResult(self.experiment, pair, None)
            matches.append(res)

        # maybe need to reconcile overlapping matches...

        # now build a fragment out of the matches...

        # and use S/W to find the transformations necessary...

        # this one feels pretty complicated? but also doable, hard to say



# NOTES:
# - handles should be treated separately, using a preprocessing filter that separates into channel files
#   - so there is no handle processing in the core code
# - we should check with luckslab: if your experimental can prefer
#   stops over muts, it makes for *much* clearer data analysis...for example see the next bullet:
# - luckslab discussion/analysis point:
#   - an insert at stop n-1 has a 1/4 chance of being
#     indistinguishable from an unmutated stop n -- this has to affect
#     statistics if inserts are at all likely. proper analysis would
#     be something like: if you see a certain prevalence of 1-nt
#     inserts at site n, take 1/3 of that number and assume that many
#     of the unmutated stops and site n+1 are in fact 1-nt
#     inserts...
#   - this kind of shit must be happening everywhere
#   - prefer stops-only experiments!!
# - another lab qn: what is significance of r1/r2 mismatch?
#   - we probably want to count them, yes -- at least as a measure of data quality (?)
#   - but doesn't this mean experiment error, rather than biologically interesting?
# - lab qn: are multiple transformations interesting or no?
#   - is the experiment likely to have multiples, or just one per fragment?



# SAH notes:
# - can we just build an RE and execute it?
#   - no, b/c of indels
