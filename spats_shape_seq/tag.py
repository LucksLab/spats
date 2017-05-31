
from partial import PartialFindProcessor
from processor import PairProcessor
from target import Targets
from util import _warn, _debug, string_match_errors


class TagProcessor(PairProcessor):

    def prepare(self):
        self._partial = PartialFindProcessor(self._run, self._targets, self._masks)
        self._partial.prepare()
        self._tag_targets = Targets()
        self._tag_targets_indexed = False

    def addTagTarget(self, name, seq):
        self._tag_targets.addTarget(name, seq)

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
        for seq, start, tags in [ (pair.r1.original_seq, 4, pair.r1.tags), (pair.r2.original_seq, 0, pair.r2.tags) ]:
            for target, match_start, match_len, match_index in self._tag_targets.find_partial_all(seq):
                if (0 == start and target.name == "adapter_b") or \
                   (4 == start and target.name == "adapter_t_rc"):
                    continue
                tag = (target.name, match_start, match_len, match_index)
                idx = 0
                while idx <= len(tags):
                    if idx == len(tags) or tag[1] < tags[idx][1]:
                        tags.insert(idx, tag)
                        break
                    idx += 1


    # note targets should include all RC's as they may be found in r1/r2
    def _find_tags_old(self, pair):
        targets = self._targets
        unknown_tag = "???"
        pair.r1.tags = [ (pair.mask.chars, 0, 4, 0) ]
        pair.r2.tags = []
        for seq, start, tags in [ (pair.r1.original_seq, 4, pair.r1.tags), (pair.r2.original_seq, 0, pair.r2.tags) ]:
            index = start
            length = len(seq)
            while True:
                print "  LU: {}".format(tags)
                lu_start, lu_len = self._longest_unmatched(seq, tags)
                print "  {}, {}".format(lu_start, lu_len)
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
        for seq, start, tags in [ (pair.r1.original_seq, 4, pair.r1.tags), (pair.r2.original_seq, 0, pair.r2.tags) ]:
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
        pair.r1.tags.append( (pair.mask.chars, 0, 4, 0) )
        if pair.r2._rtrim > 0:
            start = pair.r2.original_len - pair.r2._rtrim
            pair.r2.tags.append( (pair.mask.chars, start, min(4, pair.r2._rtrim), 0) )

    def _tag_match(self, pair):
        pair.r1.tags.append( (pair.target.name, pair.r1._ltrim, pair.r1.seq_len, pair.r1.match_index) )
        if pair.r1._rtrim > 0:
            pair.r1.tags.append( ("adapter_b", pair.r1.original_len - pair.r1._rtrim, pair.r1._rtrim, 0) )
        pair.r2.tags.insert(0, (pair.target.name, pair.r2._ltrim, pair.r2.seq_len, pair.r2.match_index) )
        if pair.r2._rtrim > 4:
            pair.r2.tags.append( ("adapter_t_rc", pair.r2.original_len - pair.r2._rtrim + 4, pair.r2._rtrim - 4, 0) )

    def process_pair(self, pair):

        self._partial.process_pair(pair)

        pair.r1.tags = []
        pair.r2.tags = []

        if pair.mask:
            self._tag_handles(pair)
        if pair.has_site:
            self._tag_match(pair)
            return

        if not self._tag_targets_indexed:
            self._tag_targets._minimum_length = 6
            self._tag_targets._index_word_length = 6
            self._tag_targets.index()

        self._find_tags(pair)

