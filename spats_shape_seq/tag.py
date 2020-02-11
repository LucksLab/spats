
from lookup import LookupProcessor
from processor import PairProcessor, Failures
from target import Targets
from util import _warn, _debug, string_match_errors
from mask import PLUS_PLACEHOLDER, MINUS_PLACEHOLDER

TAG_MATCH = "match"
TAG_ADAPTER = "adapter"
TAG_INDETERMINATE = "indeterminate"
TAG_TREATED = "treated"
TAG_UNTREATED = "untreated"
TAG_MASK_FAILURE = "mask_failure"
TAG_ADAPTER_ERRORS = "adapter_errors"
TAG_MATCH_ERRORS = "match_errors"
TAG_LINKER = "linker_cotrans"
TAG_INDELS = "indels"
TAG_INDELS_MISMATCHED = "indels_mismatched"
TAG_UNKNOWN = "unknown"
TAG_INTERESTING = "interesting"
ALL_TAGS = [ TAG_MATCH,
             TAG_ADAPTER,
             TAG_INDETERMINATE,
             TAG_TREATED,
             TAG_UNTREATED,
             TAG_MASK_FAILURE,
             TAG_ADAPTER_ERRORS,
             TAG_MATCH_ERRORS,
             TAG_LINKER,
             TAG_INDELS,
             TAG_INDELS_MISMATCHED,
             TAG_UNKNOWN,
             TAG_INTERESTING,
]
MASK_TO_TAG = { "RRRY": TAG_TREATED, "YYYR": TAG_UNTREATED, PLUS_PLACEHOLDER: TAG_TREATED, MINUS_PLACEHOLDER: TAG_UNTREATED }


class TagProcessor(PairProcessor):

    def prepare(self):
        self.uses_tags = True
        self._base_processor = self._run._get_base_processor_class()(self._run, self._targets, self._masks)
        self.counters = self._base_processor.counters
        self._base_processor.prepare()
        self._tag_targets = Targets()
        self._tag_targets_indexed = False
        self._extra_tags = []
        self._plugins = {}

    def reset_counts(self):
        self._base_processor.reset_counts()
        self.counters = self._base_processor.counters

    def setup_tags(self, pair_db):
        pair_db.setup_tags()
        pair_db.add_tags(ALL_TAGS)
        pair_db.add_tags([t.name for t in self._targets.targets])
        pair_db.add_tags(Failures.all_failures())
        pair_db.add_tags(self._extra_tags)
        pair_db.add_tags(self._plugins.keys())
        self._tagmap = pair_db.tagmap()

    def addTagTarget(self, name, seq):
        self._extra_tags.append(name)
        self._tag_targets.addTarget(name, seq)

    def addTagPlugin(self, tag, handler):
        self._plugins[tag] = handler

    def _longest_unmatched(self, seq, tags):
        lstart, llen = 0, 0
        idx = 0
        while idx <= len(tags):
            if 0 == idx:
                start = 0
                tlen = tags[idx][1] if len(tags) > 0 else len(seq)
            elif idx == len(tags):
                tag = tags[idx - 1]
                start = tag[1] + tag[2]
                tlen = len(seq) - start
            else:
                t1 = tags[idx - 1]
                t2 = tags[idx]
                start = t1[1] + t1[2]
                tlen = t2[1] - start
            if tlen > llen:
                lstart = start
                llen = tlen
            idx += 1
        return lstart, llen

    def _find_tags(self, pair):
        # this pair wasn't a match, we want to all long-enough subsequences found (even if overlaps)
        for seq, start, tags in [ (pair.r1.original_seq, pair.mask.length(), pair.r1.tags), (pair.r2.original_seq, 0, pair.r2.tags) ]:
            tagmap = {}
            for target, match_start, match_len, match_index in self._tag_targets.find_partial_all(seq):
                if (0 == start and target.name == "adapter_b") or \
                   (pair.mask.length() == start and target.name == "adapter_t_rc"):
                    continue
                tag = (target.name, match_start, match_len, match_index)
                cur = tagmap.get(target.name)
                if not cur or tag[2] > cur[2]:
                    tagmap[target.name] = tag
            tags.extend(sorted(tagmap.values(), key = lambda x : x[1]))

    def _find_tag_names(self, pair):
        tags = []
        for seq in [ pair.r1.original_seq, pair.r2.original_seq ]:
            for target, match_start, match_len, match_index in self._tag_targets.find_partial_all(seq):
                tags.append(target.name)
        return tags

    # note targets should include all RC's as they may be found in r1/r2
    def _find_tags_old(self, pair):
        targets = self._targets
        unknown_tag = "???"
        pair.r1.tags = [ (pair.mask_label, 0, pair.mask.length(), 0) ]
        pair.r2.tags = []
        for seq, start, tags in [ (pair.r1.original_seq, pair.mask.length(), pair.r1.tags), (pair.r2.original_seq, 0, pair.r2.tags) ]:
            index = start
            length = len(seq)
            while True:
                print("  LU: {}".format(tags))
                lu_start, lu_len = self._longest_unmatched(seq, tags)
                print("  {}, {}".format(lu_start, lu_len))
                if lu_len < targets._minimum_length:
                    break
                target, match_start, match_len, match_index = targets.find_partial(seq[lu_start:lu_start + lu_len])
                if target and match_len:
                    tname = target[0].name + "*" if isinstance(target, list) else target.name
                    tag = (tname, match_start + lu_start, match_len, match_index)
                    idx = 0
                    while idx <= len(tags):
                        if idx == len(tags) or tag[1] < tags[idx][1]:
                            tags.insert(idx, tag)
                            break
                        idx += 1
                else:
                    break

        #if 0
        if True:
            return
        for seq, start, tags in [ (pair.r1.original_seq, pair.mask.length(), pair.r1.tags), (pair.r2.original_seq, 0, pair.r2.tags) ]:
            index = start
            length = len(seq)
            while index < length:
                target, match_start, match_len, match_index = targets.find_partial_prefix(seq[index:])
                if target and not isinstance(target, list) and match_len:
                    tags.append( (target.name, match_start + index, match_len, match_index) )
                    index += match_len
                elif len(tags) > 0 and tags[-1][0] == unknown_tag:
                    tags[-1] = ( unknown_tag, tags[-1][1], tags[-1][2] + 1, 0)
                    index += 1
                else:
                    tags.append( ( unknown_tag, index, 1, 0) )
                    index += 1
        #endif

    def _tag_handles(self, pair):
        pair.r1.tags.append( (pair.mask_label, 0, pair.mask.length(), 0) )
        if pair.r2.rtrim > 0:
            start = pair.r2.original_len - pair.r2.rtrim
            pair.r2.tags.append( (pair.mask_label, start, min(pair.mask.length(), pair.r2.rtrim), 0) )

    def _tag_match(self, pair):
        if self._run.cotrans:
            linker_len = len(self._run.cotrans_linker)
            pair.r1.tags.append( (pair.target.name + "_rc", pair.r1.ltrim + linker_len, pair.r1.seq_len - linker_len, pair.target.n - pair.r1.match_index - pair.r1.match_len + linker_len) )
            pair.r1.tags.append( (TAG_LINKER, pair.r1.ltrim, linker_len, 0) )
        else:
            pair.r1.tags.append( (pair.target.name + "_rc", pair.r1.ltrim, pair.r1.seq_len, pair.target.n - pair.r1.match_index - pair.r1.match_len) )
        if pair.r1.rtrim > 0:
            pair.r1.tags.append( ("adapter_b", pair.r1.original_len - pair.r1.rtrim, pair.r1.rtrim, 0) )
        pair.r2.tags.insert(0, (pair.target.name, pair.r2.ltrim, pair.r2.seq_len, pair.r2.match_index) )
        if pair.r2.rtrim > pair.mask.length():
            pair.r2.tags.append( ("adapter_t_rc", pair.r2.original_len - pair.r2.rtrim + pair.mask.length(), pair.r2.rtrim - pair.mask.length(), 0) )

    def process_pair(self, pair):
        self._base_processor.process_pair(pair)
        tags = []
        if pair.has_site:
            tags.append(TAG_MATCH)
            tags.append(pair.target.name)
            if pair.r1.match_errors or pair.r2.match_errors:
                tags.append(TAG_MATCH_ERRORS)
            if pair.r1.rtrim or pair.r2.rtrim:
                tags.append(TAG_ADAPTER)
            if pair.r1.adapter_errors or pair.r2.adapter_errors:
                tags.append(TAG_ADAPTER_ERRORS)
            if self._run.cotrans:
                tags.append(TAG_LINKER)
        else:
            if not self._tag_targets_indexed:
                self._tag_targets._minimum_length = self._run.minimum_tag_match_length
                self._tag_targets._index_word_length = self._run.minimum_tag_match_length
                self._tag_targets.index()
                self._tag_targets_indexed = True
            for t in self._find_tag_names(pair):
                tname = t[:-3] if t.endswith("_rc") else t
                if tname.startswith("adapter"):
                    tags.append(TAG_ADAPTER)
                else:
                    tags.append(tname)
            if 0 == len(tags):
                tags.append(TAG_UNKNOWN)
            if pair.failure:
                tags.append(pair.failure)
        if not pair.is_determinate():
            tags.append(TAG_INDETERMINATE)
        if pair.mask:
            masktag = MASK_TO_TAG.get(pair.mask_label)
            if masktag:
                tags.append(masktag)
        else:
            tags.append(TAG_MASK_FAILURE)
        if pair.interesting:
            tags.append(TAG_INTERESTING)
        if pair.r1.indels or pair.r2.indels:
            tags.append(TAG_INDELS)
            if not pair.indels_match:
                tags.append(TAG_INDELS_MISMATCHED)

        pair.tags = [self._tagmap[t] for t in set(tags)]

        for tag, handler in self._plugins.iteritems():
            if handler(pair, tags):
                pair.tags.append(self._tagmap[tag])


    def process_pair_detail(self, pair):

        self._base_processor.process_pair(pair)

        pair.r1.tags = []
        pair.r2.tags = []

        if pair.mask:
            self._tag_handles(pair)
        if pair.has_site:
            self._tag_match(pair)
            return

        if not self._tag_targets_indexed:
            self._tag_targets._minimum_length = self._run.minimum_tag_match_length
            self._tag_targets._index_word_length = self._run.minimum_tag_match_length
            self._tag_targets.index()
            self._tag_targets_indexed = True

        self._find_tags(pair)

