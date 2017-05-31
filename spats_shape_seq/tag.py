
from processor import PairProcessor
from util import _warn, _debug, string_match_errors


class TagProcessor(PairProcessor):

    def prepare(self):
        self._targets.index()

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

    # note targets should include all RC's as they may be found in r1/r2
    def _find_tags(self, pair):
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

    def process_pair(self, pair):

        if not self._match_mask(pair):
            return

        self._find_tags(pair)
        pair.failure = "tags"

        print pair.r1.original_seq
        print pair.r1.tags
        print pair.r2.original_seq
        print pair.r2.tags
        print "-----------------------------"
