from . import fs


class SegmentMetadata(object):

    def __init__(self, segd = None, aligner = None, requiredLeft = None, requiredRight = None, availableLeft = None, availableRight = None, shiftableConstant = None):
        self.segd = segd or None  # type SegmentDescriptor
        self.aligner = aligner or None  # type Object
        self.requiredLeft = requiredLeft or 0  # type Integer
        self.requiredRight = requiredRight or 0  # type Integer
        self.availableLeft = availableLeft or 0  # type Integer
        self.availableRight = availableRight or 0  # type Integer
        self.shiftableConstant = shiftableConstant or False  # type Boolean

    @property
    def typeName(self):
        return "SegmentMetadata"

    @property
    def fsType(self):
        return fs.SegmentMetadata

    def defaultDict(self):
        return {
            'segd' : self.segd or None,
            'aligner' : self.aligner or None,
            'requiredLeft' : self.requiredLeft or 0,
            'requiredRight' : self.requiredRight or 0,
            'availableLeft' : self.availableLeft or 0,
            'availableRight' : self.availableRight or 0,
            'shiftableConstant' : self.shiftableConstant or False,
        }

    def _description(self):
        return "SegmentMetadata: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return SegmentMetadata()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.segd = SegmentDescriptor().loadFromJson(json.get('segd'), skipNull = skipNull) if ((not skipNull) or json.get('segd')) else None
        self.aligner = Object().loadFromJson(json.get('aligner'), skipNull = skipNull) if ((not skipNull) or json.get('aligner')) else None
        self.requiredLeft = json.get('requiredLeft', 0)
        self.requiredRight = json.get('requiredRight', 0)
        self.availableLeft = json.get('availableLeft', 0)
        self.availableRight = json.get('availableRight', 0)
        self.shiftableConstant = json.get('shiftableConstant', False)
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.segd != None) if minimal else (self.segd)): d['segd'] = self.segd.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.segd, 'json') else id(self.segd)
        if ((self.aligner != None) if minimal else (self.aligner)): d['aligner'] = self.aligner.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.aligner, 'json') else id(self.aligner)
        if ((self.requiredLeft != None) if minimal else (self.requiredLeft)): d['requiredLeft'] = self.requiredLeft
        if ((self.requiredRight != None) if minimal else (self.requiredRight)): d['requiredRight'] = self.requiredRight
        if ((self.availableLeft != None) if minimal else (self.availableLeft)): d['availableLeft'] = self.availableLeft
        if ((self.availableRight != None) if minimal else (self.availableRight)): d['availableRight'] = self.availableRight
        if ((self.shiftableConstant != None) if minimal else (self.shiftableConstant)): d['shiftableConstant'] = self.shiftableConstant
        return d

class SegmentDescriptor(object):

    def __init__(self, key = None, characters = None, index = None, wildcard = None, maxStart = None, minEnd = None, minLength = None, maxLength = None, maxAllowedErrors = None, handleIndels = None, partner = None):
        self.key = key or ''  # type String
        self.characters = characters or ''  # type String
        self.index = index or False  # type Boolean
        self.wildcard = wildcard or False  # type Boolean
        self.maxStart = maxStart or 0  # type Integer
        self.minEnd = minEnd or 0  # type Integer
        self.minLength = minLength or 0  # type Integer
        self.maxLength = maxLength or 0  # type Integer
        self.maxAllowedErrors = maxAllowedErrors or 0  # type Integer
        self.handleIndels = handleIndels or False  # type Boolean
        self.partner = partner or ''  # type String

    @property
    def typeName(self):
        return "SegmentDescriptor"

    @property
    def fsType(self):
        return fs.SegmentDescriptor

    def defaultDict(self):
        return {
            'key' : self.key or '',
            'characters' : self.characters or '',
            'index' : self.index or False,
            'wildcard' : self.wildcard or False,
            'maxStart' : self.maxStart or 0,
            'minEnd' : self.minEnd or 0,
            'minLength' : self.minLength or 0,
            'maxLength' : self.maxLength or 0,
            'maxAllowedErrors' : self.maxAllowedErrors or 0,
            'handleIndels' : self.handleIndels or False,
            'partner' : self.partner or '',
        }

    def _description(self):
        return "SegmentDescriptor: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return SegmentDescriptor()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.key = json.get('key', '')
        self.characters = json.get('characters', '')
        self.index = json.get('index', False)
        self.wildcard = json.get('wildcard', False)
        self.maxStart = json.get('maxStart', 0)
        self.minEnd = json.get('minEnd', 0)
        self.minLength = json.get('minLength', 0)
        self.maxLength = json.get('maxLength', 0)
        self.maxAllowedErrors = json.get('maxAllowedErrors', 0)
        self.handleIndels = json.get('handleIndels', False)
        self.partner = json.get('partner', '')
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.key != None) if minimal else (self.key)): d['key'] = self.key
        if ((self.characters != None) if minimal else (self.characters)): d['characters'] = self.characters
        if ((self.index != None) if minimal else (self.index)): d['index'] = self.index
        if ((self.wildcard != None) if minimal else (self.wildcard)): d['wildcard'] = self.wildcard
        if ((self.maxStart != None) if minimal else (self.maxStart)): d['maxStart'] = self.maxStart
        if ((self.minEnd != None) if minimal else (self.minEnd)): d['minEnd'] = self.minEnd
        if ((self.minLength != None) if minimal else (self.minLength)): d['minLength'] = self.minLength
        if ((self.maxLength != None) if minimal else (self.maxLength)): d['maxLength'] = self.maxLength
        if ((self.maxAllowedErrors != None) if minimal else (self.maxAllowedErrors)): d['maxAllowedErrors'] = self.maxAllowedErrors
        if ((self.handleIndels != None) if minimal else (self.handleIndels)): d['handleIndels'] = self.handleIndels
        if ((self.partner != None) if minimal else (self.partner)): d['partner'] = self.partner
        return d

class SequenceDescriptor(object):

    def __init__(self, segments = None):
        self.segments = segments or []  # type [SegmentDescriptor]

    @property
    def typeName(self):
        return "SequenceDescriptor"

    @property
    def fsType(self):
        return fs.SequenceDescriptor

    def defaultDict(self):
        return {
            'segments' : self.segments or [],
        }

    def _description(self):
        return "SequenceDescriptor: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return SequenceDescriptor()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.segments = [ SegmentDescriptor().loadFromJson(x, skipNull = skipNull) for x in json.get('segments') or [] ]
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.segments != None) if minimal else (self.segments)): d['segments'] = [ x.json(skipTypes = skipTypes, minimal = minimal) for x in self.segments ]
        return d

class ReadsItem(object):

    def __init__(self, key = None, characters = None, index = None, hits = None):
        self.key = key or ''  # type String
        self.characters = characters or ''  # type String
        self.index = index or None  # type Object
        self.hits = hits or 0  # type Integer

    @property
    def typeName(self):
        return "ReadsItem"

    @property
    def fsType(self):
        return fs.ReadsItem

    def defaultDict(self):
        return {
            'key' : self.key or '',
            'characters' : self.characters or '',
            'index' : self.index or None,
            'hits' : self.hits or 0,
        }

    def _description(self):
        return "ReadsItem: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return ReadsItem()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.key = json.get('key', '')
        self.characters = json.get('characters', '')
        self.index = Object().loadFromJson(json.get('index'), skipNull = skipNull) if ((not skipNull) or json.get('index')) else None
        self.hits = json.get('hits', 0)
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.key != None) if minimal else (self.key)): d['key'] = self.key
        if ((self.characters != None) if minimal else (self.characters)): d['characters'] = self.characters
        if ((self.index != None) if minimal else (self.index)): d['index'] = self.index.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.index, 'json') else id(self.index)
        if ((self.hits != None) if minimal else (self.hits)): d['hits'] = self.hits
        return d

class ReadsResult(object):

    def __init__(self, item = None, queryStart = None, queryEnd = None, itemStart = None, errors = None, indels = None, indelsDelta = None):
        self.item = item or None  # type ReadsItem
        self.queryStart = queryStart or 0  # type Integer
        self.queryEnd = queryEnd or 0  # type Integer
        self.itemStart = itemStart or 0  # type Integer
        self.errors = errors or []  # type [Integer]
        self.indels = indels or []  # type [Indel]
        self.indelsDelta = indelsDelta or 0  # type Integer

    @property
    def typeName(self):
        return "ReadsResult"

    @property
    def fsType(self):
        return fs.ReadsResult

    def defaultDict(self):
        return {
            'item' : self.item or None,
            'queryStart' : self.queryStart or 0,
            'queryEnd' : self.queryEnd or 0,
            'itemStart' : self.itemStart or 0,
            'errors' : self.errors or [],
            'indels' : self.indels or [],
            'indelsDelta' : self.indelsDelta or 0,
        }

    def _description(self):
        return "ReadsResult: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return ReadsResult()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.item = ReadsItem().loadFromJson(json.get('item'), skipNull = skipNull) if ((not skipNull) or json.get('item')) else None
        self.queryStart = json.get('queryStart', 0)
        self.queryEnd = json.get('queryEnd', 0)
        self.itemStart = json.get('itemStart', 0)
        self.errors = json.get('errors')
        self.indels = [ Indel().loadFromJson(x, skipNull = skipNull) for x in json.get('indels') or [] ]
        self.indelsDelta = json.get('indelsDelta', 0)
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.item != None) if minimal else (self.item)): d['item'] = self.item.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.item, 'json') else id(self.item)
        if ((self.queryStart != None) if minimal else (self.queryStart)): d['queryStart'] = self.queryStart
        if ((self.queryEnd != None) if minimal else (self.queryEnd)): d['queryEnd'] = self.queryEnd
        if ((self.itemStart != None) if minimal else (self.itemStart)): d['itemStart'] = self.itemStart
        if ((self.errors != None) if minimal else (self.errors)): d['errors'] = self.errors
        if ((self.indels != None) if minimal else (self.indels)): d['indels'] = [ x.json(skipTypes = skipTypes, minimal = minimal) for x in self.indels ]
        if ((self.indelsDelta != None) if minimal else (self.indelsDelta)): d['indelsDelta'] = self.indelsDelta
        return d

class Indel(object):

    def __init__(self, insertType = None, seq = None, sourceIndex = None, targetIndex = None, errorIndex = None, substSize = None):
        self.insertType = insertType or False  # type Boolean
        self.seq = seq or ''  # type String
        self.sourceIndex = sourceIndex or 0  # type Integer
        self.targetIndex = targetIndex or 0  # type Integer
        self.errorIndex = errorIndex or 0  # type Integer
        self.substSize = substSize or 0  # type Integer

    @property
    def typeName(self):
        return "Indel"

    @property
    def fsType(self):
        return fs.Indel

    def defaultDict(self):
        return {
            'insertType' : self.insertType or False,
            'seq' : self.seq or '',
            'sourceIndex' : self.sourceIndex or 0,
            'targetIndex' : self.targetIndex or 0,
            'errorIndex' : self.errorIndex or 0,
            'substSize' : self.substSize or 0,
        }

    def _description(self):
        return "Indel: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return Indel()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.insertType = json.get('insertType', False)
        self.seq = json.get('seq', '')
        self.sourceIndex = json.get('sourceIndex', 0)
        self.targetIndex = json.get('targetIndex', 0)
        self.errorIndex = json.get('errorIndex', 0)
        self.substSize = json.get('substSize', 0)
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.insertType != None) if minimal else (self.insertType)): d['insertType'] = self.insertType
        if ((self.seq != None) if minimal else (self.seq)): d['seq'] = self.seq
        if ((self.sourceIndex != None) if minimal else (self.sourceIndex)): d['sourceIndex'] = self.sourceIndex
        if ((self.targetIndex != None) if minimal else (self.targetIndex)): d['targetIndex'] = self.targetIndex
        if ((self.errorIndex != None) if minimal else (self.errorIndex)): d['errorIndex'] = self.errorIndex
        if ((self.substSize != None) if minimal else (self.substSize)): d['substSize'] = self.substSize
        return d

class SegmentResult(object):

    def __init__(self, segment = None, queryStart = None, queryEnd = None, segmentStart = None, matchTruncatedCount = None, errors = None, indels = None, indelsDelta = None, alignment = None):
        self.segment = segment or None  # type SegmentDescriptor
        self.queryStart = queryStart or 0  # type Integer
        self.queryEnd = queryEnd or 0  # type Integer
        self.segmentStart = segmentStart or 0  # type Integer
        self.matchTruncatedCount = matchTruncatedCount or 0  # type Integer
        self.errors = errors or []  # type [Integer]
        self.indels = indels or []  # type [Indel]
        self.indelsDelta = indelsDelta or 0  # type Integer
        self.alignment = alignment or None  # type Alignment

    @property
    def typeName(self):
        return "SegmentResult"

    @property
    def fsType(self):
        return fs.SegmentResult

    def defaultDict(self):
        return {
            'segment' : self.segment or None,
            'queryStart' : self.queryStart or 0,
            'queryEnd' : self.queryEnd or 0,
            'segmentStart' : self.segmentStart or 0,
            'matchTruncatedCount' : self.matchTruncatedCount or 0,
            'errors' : self.errors or [],
            'indels' : self.indels or [],
            'indelsDelta' : self.indelsDelta or 0,
            'alignment' : self.alignment or None,
        }

    def _description(self):
        return "SegmentResult: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return SegmentResult()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.segment = SegmentDescriptor().loadFromJson(json.get('segment'), skipNull = skipNull) if ((not skipNull) or json.get('segment')) else None
        self.queryStart = json.get('queryStart', 0)
        self.queryEnd = json.get('queryEnd', 0)
        self.segmentStart = json.get('segmentStart', 0)
        self.matchTruncatedCount = json.get('matchTruncatedCount', 0)
        self.errors = json.get('errors')
        self.indels = [ Indel().loadFromJson(x, skipNull = skipNull) for x in json.get('indels') or [] ]
        self.indelsDelta = json.get('indelsDelta', 0)
        self.alignment = Alignment().loadFromJson(json.get('alignment'), skipNull = skipNull) if ((not skipNull) or json.get('alignment')) else None
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.segment != None) if minimal else (self.segment)): d['segment'] = self.segment.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.segment, 'json') else id(self.segment)
        if ((self.queryStart != None) if minimal else (self.queryStart)): d['queryStart'] = self.queryStart
        if ((self.queryEnd != None) if minimal else (self.queryEnd)): d['queryEnd'] = self.queryEnd
        if ((self.segmentStart != None) if minimal else (self.segmentStart)): d['segmentStart'] = self.segmentStart
        if ((self.matchTruncatedCount != None) if minimal else (self.matchTruncatedCount)): d['matchTruncatedCount'] = self.matchTruncatedCount
        if ((self.errors != None) if minimal else (self.errors)): d['errors'] = self.errors
        if ((self.indels != None) if minimal else (self.indels)): d['indels'] = [ x.json(skipTypes = skipTypes, minimal = minimal) for x in self.indels ]
        if ((self.indelsDelta != None) if minimal else (self.indelsDelta)): d['indelsDelta'] = self.indelsDelta
        if ((self.alignment != None) if minimal else (self.alignment)): d['alignment'] = self.alignment.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.alignment, 'json') else id(self.alignment)
        return d

class SequenceResult(object):

    def __init__(self, sequence = None, success = None, failure = None, failedSegment = None, segments = None):
        self.sequence = sequence or None  # type Sequence
        self.success = success or False  # type Boolean
        self.failure = failure or ''  # type String
        self.failedSegment = failedSegment or ''  # type String
        self.segments = segments or []  # type [SegmentResult]

    @property
    def typeName(self):
        return "SequenceResult"

    @property
    def fsType(self):
        return fs.SequenceResult

    def defaultDict(self):
        return {
            'sequence' : self.sequence or None,
            'success' : self.success or False,
            'failure' : self.failure or '',
            'failedSegment' : self.failedSegment or '',
            'segments' : self.segments or [],
        }

    def _description(self):
        return "SequenceResult: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return SequenceResult()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.sequence = Sequence().loadFromJson(json.get('sequence'), skipNull = skipNull) if ((not skipNull) or json.get('sequence')) else None
        self.success = json.get('success', False)
        self.failure = json.get('failure', '')
        self.failedSegment = json.get('failedSegment', '')
        self.segments = [ SegmentResult().loadFromJson(x, skipNull = skipNull) for x in json.get('segments') or [] ]
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.sequence != None) if minimal else (self.sequence)): d['sequence'] = self.sequence.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.sequence, 'json') else id(self.sequence)
        if ((self.success != None) if minimal else (self.success)): d['success'] = self.success
        if ((self.failure != None) if minimal else (self.failure)): d['failure'] = self.failure
        if ((self.failedSegment != None) if minimal else (self.failedSegment)): d['failedSegment'] = self.failedSegment
        if ((self.segments != None) if minimal else (self.segments)): d['segments'] = [ x.json(skipTypes = skipTypes, minimal = minimal) for x in self.segments ]
        return d

class Fragment(object):

    def __init__(self, seq = None, r1 = None, r2 = None, r1IsRC = None, success = None, failure = None, overlap = None, reverseOrder = None, alignment = None, errorCount = None, errors = None, indels = None, indelsDelta = None, qualityResolved = None):
        self.seq = seq or None  # type Sequence
        self.r1 = r1 or None  # type Sequence
        self.r2 = r2 or None  # type Sequence
        self.r1IsRC = r1IsRC or False  # type Boolean
        self.success = success or False  # type Boolean
        self.failure = failure or ''  # type String
        self.overlap = overlap or 0  # type Integer
        self.reverseOrder = reverseOrder or False  # type Boolean
        self.alignment = alignment or None  # type Alignment
        self.errorCount = errorCount or 0  # type Integer
        self.errors = errors or []  # type [Integer]
        self.indels = indels or []  # type [Indel]
        self.indelsDelta = indelsDelta or 0  # type Integer
        self.qualityResolved = qualityResolved or []  # type [Integer]

    @property
    def typeName(self):
        return "Fragment"

    @property
    def fsType(self):
        return fs.Fragment

    def defaultDict(self):
        return {
            'seq' : self.seq or None,
            'r1' : self.r1 or None,
            'r2' : self.r2 or None,
            'r1IsRC' : self.r1IsRC or False,
            'success' : self.success or False,
            'failure' : self.failure or '',
            'overlap' : self.overlap or 0,
            'reverseOrder' : self.reverseOrder or False,
            'alignment' : self.alignment or None,
            'errorCount' : self.errorCount or 0,
            'errors' : self.errors or [],
            'indels' : self.indels or [],
            'indelsDelta' : self.indelsDelta or 0,
            'qualityResolved' : self.qualityResolved or [],
        }

    def _description(self):
        return "Fragment: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return Fragment()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.seq = Sequence().loadFromJson(json.get('seq'), skipNull = skipNull) if ((not skipNull) or json.get('seq')) else None
        self.r1 = Sequence().loadFromJson(json.get('r1'), skipNull = skipNull) if ((not skipNull) or json.get('r1')) else None
        self.r2 = Sequence().loadFromJson(json.get('r2'), skipNull = skipNull) if ((not skipNull) or json.get('r2')) else None
        self.r1IsRC = json.get('r1IsRC', False)
        self.success = json.get('success', False)
        self.failure = json.get('failure', '')
        self.overlap = json.get('overlap', 0)
        self.reverseOrder = json.get('reverseOrder', False)
        self.alignment = Alignment().loadFromJson(json.get('alignment'), skipNull = skipNull) if ((not skipNull) or json.get('alignment')) else None
        self.errorCount = json.get('errorCount', 0)
        self.errors = json.get('errors')
        self.indels = [ Indel().loadFromJson(x, skipNull = skipNull) for x in json.get('indels') or [] ]
        self.indelsDelta = json.get('indelsDelta', 0)
        self.qualityResolved = json.get('qualityResolved')
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.seq != None) if minimal else (self.seq)): d['seq'] = self.seq.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.seq, 'json') else id(self.seq)
        if ((self.r1 != None) if minimal else (self.r1)): d['r1'] = self.r1.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.r1, 'json') else id(self.r1)
        if ((self.r2 != None) if minimal else (self.r2)): d['r2'] = self.r2.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.r2, 'json') else id(self.r2)
        if ((self.r1IsRC != None) if minimal else (self.r1IsRC)): d['r1IsRC'] = self.r1IsRC
        if ((self.success != None) if minimal else (self.success)): d['success'] = self.success
        if ((self.failure != None) if minimal else (self.failure)): d['failure'] = self.failure
        if ((self.overlap != None) if minimal else (self.overlap)): d['overlap'] = self.overlap
        if ((self.reverseOrder != None) if minimal else (self.reverseOrder)): d['reverseOrder'] = self.reverseOrder
        if ((self.alignment != None) if minimal else (self.alignment)): d['alignment'] = self.alignment.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.alignment, 'json') else id(self.alignment)
        if ((self.errorCount != None) if minimal else (self.errorCount)): d['errorCount'] = self.errorCount
        if ((self.errors != None) if minimal else (self.errors)): d['errors'] = self.errors
        if ((self.indels != None) if minimal else (self.indels)): d['indels'] = [ x.json(skipTypes = skipTypes, minimal = minimal) for x in self.indels ]
        if ((self.indelsDelta != None) if minimal else (self.indelsDelta)): d['indelsDelta'] = self.indelsDelta
        if ((self.qualityResolved != None) if minimal else (self.qualityResolved)): d['qualityResolved'] = self.qualityResolved
        return d

class FragmentResult(object):

    def __init__(self, r1 = None, r2 = None, f = None, barcode = None, s = None, c = None, r = None):
        self.r1 = r1 or None  # type Sequence
        self.r2 = r2 or None  # type Sequence
        self.f = f or None  # type Fragment
        self.barcode = barcode or ''  # type String
        self.s = s or None  # type SequenceResult
        self.c = c or None  # type ClassifyResult
        self.r = r or []  # type [ReadsResult]

    @property
    def typeName(self):
        return "FragmentResult"

    @property
    def fsType(self):
        return fs.FragmentResult

    def defaultDict(self):
        return {
            'r1' : self.r1 or None,
            'r2' : self.r2 or None,
            'f' : self.f or None,
            'barcode' : self.barcode or '',
            's' : self.s or None,
            'c' : self.c or None,
            'r' : self.r or [],
        }

    def _description(self):
        return "FragmentResult: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return FragmentResult()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.r1 = Sequence().loadFromJson(json.get('r1'), skipNull = skipNull) if ((not skipNull) or json.get('r1')) else None
        self.r2 = Sequence().loadFromJson(json.get('r2'), skipNull = skipNull) if ((not skipNull) or json.get('r2')) else None
        self.f = Fragment().loadFromJson(json.get('f'), skipNull = skipNull) if ((not skipNull) or json.get('f')) else None
        self.barcode = json.get('barcode', '')
        self.s = SequenceResult().loadFromJson(json.get('s'), skipNull = skipNull) if ((not skipNull) or json.get('s')) else None
        self.c = ClassifyResult().loadFromJson(json.get('c'), skipNull = skipNull) if ((not skipNull) or json.get('c')) else None
        self.r = [ ReadsResult().loadFromJson(x, skipNull = skipNull) for x in json.get('r') or [] ]
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.r1 != None) if minimal else (self.r1)): d['r1'] = self.r1.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.r1, 'json') else id(self.r1)
        if ((self.r2 != None) if minimal else (self.r2)): d['r2'] = self.r2.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.r2, 'json') else id(self.r2)
        if ((self.f != None) if minimal else (self.f)): d['f'] = self.f.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.f, 'json') else id(self.f)
        if ((self.barcode != None) if minimal else (self.barcode)): d['barcode'] = self.barcode
        if ((self.s != None) if minimal else (self.s)): d['s'] = self.s.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.s, 'json') else id(self.s)
        if ((self.c != None) if minimal else (self.c)): d['c'] = self.c.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.c, 'json') else id(self.c)
        if ((self.r != None) if minimal else (self.r)): d['r'] = [ x.json(skipTypes = skipTypes, minimal = minimal) for x in self.r ]
        return d

class Range(object):

    def __init__(self, start = None, end = None, segment = None, text = None):
        self.start = start or 0  # type Integer
        self.end = end or 0  # type Integer
        self.segment = segment or None  # type SegmentDescriptor
        self.text = text or ''  # type String

    @property
    def typeName(self):
        return "Range"

    @property
    def fsType(self):
        return fs.Range

    def defaultDict(self):
        return {
            'start' : self.start or 0,
            'end' : self.end or 0,
            'segment' : self.segment or None,
            'text' : self.text or '',
        }

    def _description(self):
        return "Range: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return Range()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.start = json.get('start', 0)
        self.end = json.get('end', 0)
        self.segment = SegmentDescriptor().loadFromJson(json.get('segment'), skipNull = skipNull) if ((not skipNull) or json.get('segment')) else None
        self.text = json.get('text', '')
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.start != None) if minimal else (self.start)): d['start'] = self.start
        if ((self.end != None) if minimal else (self.end)): d['end'] = self.end
        if ((self.segment != None) if minimal else (self.segment)): d['segment'] = self.segment.json(skipTypes = skipTypes, minimal = minimal) if hasattr(self.segment, 'json') else id(self.segment)
        if ((self.text != None) if minimal else (self.text)): d['text'] = self.text
        return d

class Gap(object):

    def __init__(self, queryStart = None, queryEnd = None, segmentStart = None, segmentEnd = None):
        self.queryStart = queryStart or 0  # type Integer
        self.queryEnd = queryEnd or 0  # type Integer
        self.segmentStart = segmentStart or 0  # type Integer
        self.segmentEnd = segmentEnd or 0  # type Integer

    @property
    def typeName(self):
        return "Gap"

    @property
    def fsType(self):
        return fs.Gap

    def defaultDict(self):
        return {
            'queryStart' : self.queryStart or 0,
            'queryEnd' : self.queryEnd or 0,
            'segmentStart' : self.segmentStart or 0,
            'segmentEnd' : self.segmentEnd or 0,
        }

    def _description(self):
        return "Gap: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return Gap()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.queryStart = json.get('queryStart', 0)
        self.queryEnd = json.get('queryEnd', 0)
        self.segmentStart = json.get('segmentStart', 0)
        self.segmentEnd = json.get('segmentEnd', 0)
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.queryStart != None) if minimal else (self.queryStart)): d['queryStart'] = self.queryStart
        if ((self.queryEnd != None) if minimal else (self.queryEnd)): d['queryEnd'] = self.queryEnd
        if ((self.segmentStart != None) if minimal else (self.segmentStart)): d['segmentStart'] = self.segmentStart
        if ((self.segmentEnd != None) if minimal else (self.segmentEnd)): d['segmentEnd'] = self.segmentEnd
        return d

class AlignedSegment(object):

    def __init__(self, queryIndex = None, segmentIndex = None, length = None):
        self.queryIndex = queryIndex or 0  # type Integer
        self.segmentIndex = segmentIndex or 0  # type Integer
        self.length = length or 0  # type Integer

    @property
    def typeName(self):
        return "AlignedSegment"

    @property
    def fsType(self):
        return fs.AlignedSegment

    def defaultDict(self):
        return {
            'queryIndex' : self.queryIndex or 0,
            'segmentIndex' : self.segmentIndex or 0,
            'length' : self.length or 0,
        }

    def _description(self):
        return "AlignedSegment: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return AlignedSegment()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.queryIndex = json.get('queryIndex', 0)
        self.segmentIndex = json.get('segmentIndex', 0)
        self.length = json.get('length', 0)
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.queryIndex != None) if minimal else (self.queryIndex)): d['queryIndex'] = self.queryIndex
        if ((self.segmentIndex != None) if minimal else (self.segmentIndex)): d['segmentIndex'] = self.segmentIndex
        if ((self.length != None) if minimal else (self.length)): d['length'] = self.length
        return d

class Alignment(object):

    def __init__(self, segments = None, score = None):
        self.segments = segments or []  # type [AlignedSegment]
        self.score = score or 0  # type Integer

    @property
    def typeName(self):
        return "Alignment"

    @property
    def fsType(self):
        return fs.Alignment

    def defaultDict(self):
        return {
            'segments' : self.segments or [],
            'score' : self.score or 0,
        }

    def _description(self):
        return "Alignment: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return Alignment()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.segments = [ AlignedSegment().loadFromJson(x, skipNull = skipNull) for x in json.get('segments') or [] ]
        self.score = json.get('score', 0)
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.segments != None) if minimal else (self.segments)): d['segments'] = [ x.json(skipTypes = skipTypes, minimal = minimal) for x in self.segments ]
        if ((self.score != None) if minimal else (self.score)): d['score'] = self.score
        return d

class Sequence(object):

    def __init__(self, sequenceId = None, name = None, identifier = None, characters = None, quality = None):
        self.sequenceId = sequenceId or 0  # type Integer
        self.name = name or ''  # type String
        self.identifier = identifier or ''  # type String
        self.characters = characters or ''  # type String
        self.quality = quality or ''  # type String

    @property
    def typeName(self):
        return "Sequence"

    @property
    def fsType(self):
        return fs.Sequence

    def defaultDict(self):
        return {
            'sequenceId' : self.sequenceId or 0,
            'name' : self.name or '',
            'identifier' : self.identifier or '',
            'characters' : self.characters or '',
            'quality' : self.quality or '',
        }

    def _description(self):
        return "Sequence: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return Sequence()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.sequenceId = json.get('sequenceId', 0)
        self.name = json.get('name', '')
        self.identifier = json.get('identifier', '')
        self.characters = json.get('characters', '')
        self.quality = json.get('quality', '')
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.sequenceId != None) if minimal else (self.sequenceId)): d['sequenceId'] = self.sequenceId
        if ((self.name != None) if minimal else (self.name)): d['name'] = self.name
        if ((self.identifier != None) if minimal else (self.identifier)): d['identifier'] = self.identifier
        if ((self.characters != None) if minimal else (self.characters)): d['characters'] = self.characters
        if ((self.quality != None) if minimal else (self.quality)): d['quality'] = self.quality
        return d
