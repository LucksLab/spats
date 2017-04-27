
from mask import longest_match

class Target(object):

    def __init__(self, name, seq, index_word_length = 8):
        self.name = name
        self.seq = seq
        self.n = len(seq)
        self._index = None
        self._index_word_length = index_word_length
        self._minimum_length = index_word_length

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
        seq = self.seq
        index = {}
        word_len = self._index_word_length
        for i in range(len(seq) - word_len + 1):
            key = seq[i:(i + word_len)]
            sites = index.get(key)
            if not sites:
                sites = []
                index[key] = sites
            sites.append(i)
        self._index = index

    def longest_self_match(self, minimum_length = None):
        seq = self.seq
        seq_len = len(self.seq)
        min_len = self._minimum_length
        candidate = (None, None, None)
        index = 0
        while True:
            query = seq[index:index+min_len]
            match_index = self.find_exact(query, exclude = index)
            if match_index >= 0:
                assert(match_index != index)
                left, right = longest_match(seq, (index, min_len), seq, (match_index, min_len))
                total_len = min_len + left + right
                if not candidate[1] or total_len > candidate[1]:
                    candidate = (index, total_len, match_index)
            index += 1
            if index >= seq_len - max(min_len, candidate[1]):
                break
        # note that we can't say there's no self-match below min_len,
        # so if we didn't find one, just return that
        print "LM: {}".format(candidate)
        return candidate[1] or min_len

    def find_exact(self, query, exclude = None):
        word_len = self._index_word_length
        query_len = len(query)
        if query_len < word_len:
            raise Exception("Query too short: len({}) < {}".format(query, word_len))
        query_key = query[0:word_len]
        for index in self._index.get(query_key, []):
            if index == exclude and exclude is not None:
                continue
            if self.seq[index:index+query_len] == query:
                return index
        return -1

    # returns (query_start_index, match_len, sequence_index), where:
    #  query_start_index: the index into the query where the match starts
    #  match_len: the length of the match
    #  sequence_index: the index into the target sequence where the match starts
    def find_partial(self, query):
        # it's much faster to search for longer partial matches, then fall back on the minimum
        candidate = self._find_partial(query, multiple = 2)
        return candidate if candidate[1] else self._find_partial(query, multiple = 1)

    #@profile
    def _find_partial(self, query, multiple):
        min_len = self._minimum_length * multiple
        word_len = self._index_word_length
        check_every = max(min_len - word_len, 1) # norah has proved that this guarantees finding a match if it exists
        query_len = len(query)
        last = query_len - max(check_every, word_len)
        check_sites = range(0, last, check_every)
        check_sites.append(last)
        candidate = (None, None, None)
        # NOTE: it's important to check all sites, and all hits -- to find the longest match.
        for site in check_sites:
            site_key = query[site:site+word_len]
            #print "CS: {}, {}".format(site, site_key)
            for index in self._index.get(site_key, []):
                #print "GOT: " + str(index)
                left, right = longest_match(query, (site, word_len), self.seq, (index, word_len))
                total_len = left + right + word_len
                #print "extends: <--{}, -->{} / {} ({})".format(left, right, total_len, min_len)
                if total_len >= min_len:
                    if total_len >= query_len:
                        # we can return immediately if we've got a full match
                        # BUT: xref notes -- need to make sure there can't be multiple in target
                        # some target-analysis preprocessing will be useful here
                        return site - left, total_len, index - left
                    elif not candidate[1] or total_len > candidate[1]:
                        # ...otherwise, keep it if it's the best match so far
                        candidate = (site - left, total_len, index - left)
                        #print "C: {}".format(candidate)
        return candidate
