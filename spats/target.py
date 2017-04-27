
from mask import longest_match

class _Target(object):

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.n = len(seq)


class Targets(object):

    def __init__(self, index_word_length = 8):
        self._index = None
        self._index_word_length = index_word_length
        self._minimum_length = index_word_length
        self.targets = []

    def addTarget(self, name, seq):
        self.targets.append(_Target(name, seq))

    @property
    def minimum_match_length(self):
        return self._minimum_length

    @minimum_match_length.setter
    def minimum_match_length(self, min_length):
        if min_length < self._index_word_length:
            raise Exception("minimum_length too short: {} < {}".format(min_length, self._index_word_length))
        self._minimum_length = min_length
        # disable this warning for now
        if False and min_length - self._index_word_length < 4:
            print "Warning: minimum_length {} is not much longer than index length {}".format(min_len, word_len)

    def index(self):
        index = {}
        word_len = self._index_word_length
        for target in self.targets:
            seq = target.seq
            n = target.n
            for i in range(n - word_len + 1):
                key = seq[i:(i + word_len)]
                sites = index.get(key)
                if not sites:
                    sites = []
                    index[key] = sites
                sites.append((target, i))
        self._index = index

    def longest_self_match(self, minimum_length = None):
        min_len = self._minimum_length
        candidate = (None, None, None)
        for target in self.targets:
            seq = target.seq
            seq_len = target.n
            index = 0
            while True:
                query = seq[index:index+min_len]
                match_target, match_index = self.find_exact(query, exclude = (target, index))
                if match_target:
                    assert(match_target != target or match_index != index)
                    left, right = longest_match(seq, (index, min_len), match_target.seq, (match_index, min_len))
                    total_len = min_len + left + right
                    if not candidate[1] or total_len > candidate[1]:
                        candidate = (index, total_len, match_index)
                        #print "C {} {} matches {}".format(target.name, candidate, match_target.name)
                index += 1
                if index >= seq_len - max(min_len, candidate[1]):
                    break
        # note that we can't say there's no self-match below min_len,
        # so if we didn't find one, just return that
        return candidate[1] or min_len

    def find_exact(self, query, exclude = (None, -1)):
        word_len = self._index_word_length
        query_len = len(query)
        if query_len < word_len:
            raise Exception("Query too short: len({}) < {}".format(query, word_len))
        query_key = query[0:word_len]
        for target, index in self._index.get(query_key, []):
            if index == exclude[1] and target == exclude[0]:
                continue
            if target.seq[index:index+query_len] == query:
                return target, index
        return None, -1

    # returns (query_start_index, match_len, sequence_index), where:
    #  query_start_index: the index into the query where the match starts
    #  match_len: the length of the match
    #  sequence_index: the index into the target sequence where the match starts
    def find_partial(self, query, force_target = None):
        # it's much faster to search for longer partial matches, then fall back on the minimum
        candidate = self._find_partial(query, force_target, multiple = 2)
        return candidate if candidate[1] else self._find_partial(query, force_target, multiple = 1)

    #@profile
    def _find_partial(self, query, force_target, multiple):
        min_len = self._minimum_length * multiple
        word_len = self._index_word_length
        check_every = max(min_len - word_len, 1) # norah has proved that this guarantees finding a match if it exists
        query_len = len(query)
        last = query_len - max(check_every, word_len)
        check_sites = range(0, last, check_every)
        check_sites.append(last)
        candidate = (None, None, None, None)
        # NOTE: it's important to check all sites, and all hits -- to find the longest match.
        for site in check_sites:
            site_key = query[site:site+word_len]
            #print "CS: {}, {}".format(site, site_key)
            for target, index in self._index.get(site_key, []):
                if force_target and target != force_target:
                    continue
                #print "GOT: " + str(index)
                left, right = longest_match(query, (site, word_len), target.seq, (index, word_len))
                total_len = left + right + word_len
                #print "extends: <--{}, -->{} / {} ({})".format(left, right, total_len, min_len)
                if total_len >= min_len:
                    if total_len >= query_len:
                        # we can return immediately if we've got a full match
                        # BUT: xref notes -- need to make sure there can't be multiple in target
                        # some target-analysis preprocessing will be useful here
                        return target, site - left, total_len, index - left
                    elif not candidate[2] or total_len > candidate[2]:
                        # ...otherwise, keep it if it's the best match so far
                        candidate = (target, site - left, total_len, index - left)
                        #print "C: {}".format(candidate)
        return candidate
