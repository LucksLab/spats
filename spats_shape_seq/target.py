
from mask import longest_match
from util import reverse_complement, _debug, _warn

class _Target(object):

    def __init__(self, name, seq, rowid):
        self.rowid = rowid
        self.name = name
        self.seq = str(seq).upper()
        self.n = len(self.seq)


class Targets(object):

    def __init__(self, index_word_length = 8):
        self._index = None
        self._index_word_length = index_word_length
        self._minimum_length = index_word_length
        self.targets = []

    def addTarget(self, name, seq, rowid = -1):
        self.targets.append(_Target(name, seq, rowid))

    def merge_target(self, name, seq, rowid = -1):
        # Don't add "duplicate" targets.
        for t in self.targets:
            if t.name == name and t.seq == seq:
                break
        else:
            self.targets.append(_Target(name, seq, rowid))

    @property
    def minimum_match_length(self):
        return self._minimum_length

    @minimum_match_length.setter
    def minimum_match_length(self, min_length):
        if min_length < self._index_word_length:
            raise Exception("minimum_length too short: {} < {}".format(min_length, self._index_word_length))
        if min_length > 35:
            raise Exception("minimum length {} is too long, do you need to manually set?".format(min_length))
        self._minimum_length = min_length
        # disable this warning for now
        if False and min_length - self._index_word_length < 4:
            print("Warning: minimum_length {} is not much longer than index length {}".format(min_len, word_len))

    def index(self):
        index = {}
        word_len = self._index_word_length
        for target in self.targets:
            seq = target.seq
            n = target.n
            for i in xrange(n - word_len + 1):
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
                        #print("C {} {} matches {}".format(target.name, candidate, match_target.name))
                index += 1
                if index >= seq_len - max(min_len, candidate[1]):
                    break
        # note that we can't say there's no self-match below min_len,
        # so if we didn't find one, just return that
        return candidate[1] or min_len

    def longest_target_self_matches(self, minimum_length = None):
        min_len = self._minimum_length
        matches = {}
        for target in self.targets:
            seq = target.seq
            seq_len = target.n
            index = 0
            candidate = (None, None, None)
            while True:
                query = seq[index:index+min_len]
                match_target, match_index = self.find_exact(query, exclude = (target, index))
                if match_target == target:
                    assert(match_index != index)
                    left, right = longest_match(seq, (index, min_len), match_target.seq, (match_index, min_len))
                    total_len = min_len + left + right
                    if not candidate[1] or total_len > candidate[1]:
                        candidate = (index, total_len, match_index)
                        #print("C {} {} matches {}".format(target.name, candidate, match_target.name))
                index += 1
                if index >= seq_len - max(min_len, candidate[1]):
                    break
            # note that we can't say there's no self-match below min_len,
            # so if we didn't find one, just return that
            matches[target.name] = candidate[1] or min_len
        return matches

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

    def find_partial_prefix(self, query):
        min_len = self._minimum_length
        word_len = self._index_word_length
        query_len = len(query)
        candidate = [None, None, None, None]
        if query_len < word_len:
            return candidate
        site_key = query[:word_len]
        #print("  fpp query: " + query + " --> " + site_key)
        for target, index in self._index.get(site_key, []):
            left, right = longest_match(query, (0, word_len), target.seq, (index, word_len))
            total_len = left + right + word_len
            #print("    extends: <--{}, -->{} / {} ({})".format(left, right, total_len, min_len))
            if total_len >= min_len:
                if not candidate[2] or total_len > candidate[2]:
                    # keep it if it's the best match so far
                    candidate = [target, 0 - left, total_len, index - left]
                    #print("C: {}".format(candidate))
                elif total_len == candidate[2] and target != candidate[0]:
                    # need to keep track if multiple candidates have this same max length
                    if isinstance(candidate[0], list):
                        candidate[0] = candidate[0] if target in candidate[0] else [ target ] + candidate[0]
                    else:
                        candidate[0] = [ target, candidate[0] ]
        return candidate

    def find_partial_all(self, query, min_length = 0):
        min_len = min_length or self._minimum_length
        word_len = self._index_word_length
        check_every = max(min_len - word_len, 1) # norah has proved that this guarantees finding a match if it exists
        query_len = len(query)
        last = query_len - max(check_every, word_len)
        check_sites = range(0, last, check_every)
        check_sites.append(last)
        candidates = []
        # NOTE: it's important to check all sites, and all hits -- to find the longest match.
        for site in check_sites:
            site_key = query[site:site+word_len]
            #print("CS: {}, {}".format(site, site_key))
            for target, index in self._index.get(site_key, []):
                #print("GOT: " + str(index) + " -- " + target.name)
                left, right = longest_match(query, (site, word_len), target.seq, (index, word_len))
                total_len = left + right + word_len
                #print("extends: <--{}, -->{} / {} ({})".format(left, right, total_len, min_len))
                if total_len >= min_len:
                    candidate = [target, site - left, total_len, index - left]
                    if candidate not in candidates:
                        candidates.append(candidate)
        return candidates
        
    # returns ([target or targets], query_start_index, match_len, sequence_index), where:
    #  target or targets: if a single target match, then the target; otherwise a list of all matched targets, rest of params corresponding to first one
    #  query_start_index: the index into the query where the match starts
    #  match_len: the length of the match
    #  sequence_index: the index into the target sequence where the match starts
    def find_partial(self, query, force_target = None, min_length_override = 0):
        # it's much faster to search for longer partial matches, then fall back on the minimum
        candidate = self._find_partial(query, force_target, min_length_override, multiple = 2)
        return candidate if candidate[0] else self._find_partial(query, force_target, min_length_override, multiple = 1)

    #@profile
    def _find_partial(self, query, force_target, min_length_override, multiple):
        min_len = min_length_override or self._minimum_length * multiple
        word_len = self._index_word_length
        check_every = max(min_len - word_len, 1) # norah has proved that this guarantees finding a match if it exists
        query_len = len(query)
        last = query_len - max(check_every, word_len)
        check_sites = range(0, last, check_every)
        check_sites.append(last)
        candidate = [None, None, None, None]
        # NOTE: it's important to check all sites, and all hits -- to find the longest match.
        for site in check_sites:
            site_key = query[site:site+word_len]
            #print("CS: {}, {}".format(site, site_key))
            for target, index in self._index.get(site_key, []):
                if force_target and target != force_target:
                    continue
                #print("GOT: " + str(index) + " -- " + target.name)
                left, right = longest_match(query, (site, word_len), target.seq, (index, word_len))
                total_len = left + right + word_len
                #print("extends: <--{}, -->{} / {} ({})".format(left, right, total_len, min_len))
                if total_len >= min_len:
                    if not candidate[2] or total_len > candidate[2]:
                        # keep it if it's the best match so far
                        candidate = [target, site - left, total_len, index - left]
                        #print("C: {}".format(candidate))
                    elif total_len == candidate[2] and target != candidate[0]:
                        # need to keep track if multiple candidates have this same max length
                        if isinstance(candidate[0], list):
                            candidate[0] = candidate[0] if target in candidate[0] else [ target ] + candidate[0]
                        else:
                            candidate[0] = [ target, candidate[0] ]
        return candidate

    def build_cotrans_lookups(self, run):
        # for cotrans experiments, R1 includes a linker and is relatively restricted
        # store RC in table so that we can directly lookup R1[4:]

        # TODO: this could be set to run.minimum_target_match_length, but given that there's linker
        # involved, it makes sense to assume that a match including linker will hit that minimum.
        # so for this we can include very small bits of the actual target.
        minimum_target_length = 3

        linker = run.cotrans_linker
        linker_len = len(linker)
        r1_table = {}

        if 1 != len(self.targets):
            raise Exception("cotrans requires one target")
        target = self.targets[0]
        tseq = target.seq
        tlen = len(tseq)
        adapter_b = run.adapter_b
        pair_len = run.pair_length
        assert(0 < pair_len)
        masklen = 4    # TODO
        r1_match_len = pair_len - masklen

        for end in xrange(minimum_target_length, tlen + 1):
            target_subseq = tseq[:end]
            for i in xrange(0, r1_match_len - linker_len - minimum_target_length + 1):
                tstart = i - (r1_match_len - linker_len)
                if tstart + end < 0:
                    continue
                target_bit = target_subseq[tstart:]
                r1_rc_match = target_bit + linker
                r1_match = reverse_complement(r1_rc_match) + adapter_b[:i]
                entries = r1_table.get(r1_match, [])
                entries.append( (target, end, i, []) ) # target, end, amount of adapter to trim, mutations
                r1_table[r1_match] = entries

                if run.count_mutations:
                    bit_len = len(target_bit)
                    for toggle_idx in xrange(bit_len):
                        for nt in [ 'A', 'C', 'G', 'T' ]:
                            if target_bit[toggle_idx] == nt:
                                continue
                            mutated_bit = target_bit[:toggle_idx] + nt + target_bit[toggle_idx + 1:]
                            mutated_rc_match = mutated_bit + linker
                            mutated_match = reverse_complement(mutated_rc_match) + adapter_b[:i]
                            entries = r1_table.get(mutated_match, [])
                            entries.append( (target, end, i, [end - (bit_len - toggle_idx) + 1]) )
                            r1_table[mutated_match] = entries

        self.r1_lookup = r1_table
        self._build_R1_aliases(adapter_b, r1_match_len)

        # we only need to build R2 lookups for full sequences (excepting linker)
        # trim cases are just tested against R1
        self._build_R2_lookup(pair_len - linker_len - masklen, run.count_mutations)


    def build_lookups(self, run, length = None, end_only = True):
        use_length = length or 35
        if use_length < 0:
            raise Exception('Cannot build lookups on variable-length inputs. Use find_partial processor.')
        masklen = 4  # TODO
        self._build_R1_lookup(run.adapter_b, use_length - masklen, end_only, run.count_mutations, run.dumbbell)
        self._build_R2_lookup(use_length - masklen, run.count_mutations, run.dumbbell)

    def _build_R1_lookup(self, adapter_b, length = 31, end_only = True, mutations = False, dumbbell = None):
        # we can pre-build the set of all possible (error-free) R1, b/c:
        #  - R1 has to include the right-most nt
        #  - R1 can include some adapter-b (or dumbbell) off the end
        #  - this is done for each target
        #  - note that in cases where R1 includes some (enough) adapter, then position and content of R2 is determined
        # note that this does *not* include the handle.
        r1_table = {}
        use_aliases = False
        for target in self.targets:
            tlen = target.n
            rc_tgt = reverse_complement(target.seq)
            rc_dumbbell = reverse_complement(dumbbell) if dumbbell else None
            tcandidates = 0
            for i in xrange(1, length + 1):
                if rc_dumbbell:
                    if length - i <= len(rc_dumbbell):
                        r1_candidate = rc_tgt[:i] + rc_dumbbell[:length - i]
                    else:
                        r1_candidate = rc_tgt[:i] + rc_dumbbell + adapter_b[:length - len(rc_dumbbell) - i]
                else:
                    r1_candidate = rc_tgt[:i] + adapter_b[:length - i]
                res = (target, None if i == length else tlen - i, length - i, []) # target, end, amount of adapter to trim, mutations
                existing = r1_table.get(r1_candidate)
                if existing:
                    existing.append(res)
                else:
                    r1_table[r1_candidate] = [ res ]
                tcandidates += 1
                if mutations:
                    for toggle_idx in xrange(i):
                        for nt in [ 'A', 'C', 'G', 'T' ]:
                            if r1_candidate[toggle_idx] == nt:
                                continue
                            mutated_bit = r1_candidate[:toggle_idx] + nt + r1_candidate[toggle_idx + 1:]
                            mres = (res[0], res[1], res[2], [ tlen - toggle_idx ])
                            existing = r1_table.get(mutated_bit)
                            if existing:
                                existing.append(mres)
                            else:
                                r1_table[mutated_bit] = [ mres ]

            if 0 == tcandidates:
                _warn("!! No R1 match candidates for {}".format(target.name))

        self.r1_lookup = r1_table
        self._build_R1_aliases(adapter_b, length)


    # if we don't have enough adapter_b, then make aliases for short lookup keys
    def _build_R1_aliases(self, adapter_b, length):

        self.r1_aliases = None
        self.r1_lookup_length = length

        use_aliases = False
        for key in self.r1_lookup.keys():
            if len(key) < length:
                use_aliases = True
                break
        if not use_aliases:
            return

        minimum_length = 1 + len(adapter_b)
        if use_aliases:
            r1_aliases = {}
            for key in self.r1_lookup.keys():
                if len(key) == length:
                    continue
                alias = key[:minimum_length]
                if alias in r1_aliases:
                    r1_aliases[alias].append(key)
                else:
                    r1_aliases[alias] = [ key ]
            self.r1_aliases = r1_aliases
            self.r1_lookup_length = minimum_length


    def _build_R2_lookup(self, length = 35, mutations = False, dumbbell = None):
        # for the R2 table, we only care about R2's that are in the sequence
        # when R2 needs adapter trimming, R1 will determine that
        self_matches = self.longest_target_self_matches()
        #for tname in self_matches.keys():
        #    print("{} : {}".format(tname, self_matches[tname]))
        r2_full_table = {}
        r2_match_lengths = {}
        for target in self.targets:
            mlen = length if mutations else self_matches[target.name] + 1
            if dumbbell:
                mlen += len(dumbbell)
            if length < mlen:
                raise Exception("R2 length not long enough for target self-match ({} / {})".format(length, mlen))
            r2_table = {}
            r2_full_table[target.name] = r2_table
            r2_match_lengths[target.name] = mlen
            tlen = target.n
            tgt_seq = target.seq
            for i in xrange(tlen - mlen + 1):
                if dumbbell:
                    r2_candidate = dumbbell + tgt_seq[i:i+mlen-len(dumbbell)]
                else:
                    r2_candidate = tgt_seq[i:i+mlen]
                if r2_table.get(r2_candidate):
                    raise Exception("indeterminate R2 candidate {} in target?".format(r2_candidate))
                r2_table[r2_candidate] = (i, [])

                if mutations:
                    bit_len = len(r2_candidate)
                    for toggle_idx in xrange(bit_len):
                        for nt in [ 'A', 'C', 'G', 'T' ]:
                            if r2_candidate[toggle_idx] == nt:
                                continue
                            mutated_bit = r2_candidate[:toggle_idx] + nt + r2_candidate[toggle_idx + 1:]
                            if r2_table.get(mutated_bit):
                                raise Exception("indeterminate R2 candidate {} in target?".format(mutated_bit))
                            r2_table[mutated_bit] = (i, [ i + toggle_idx + 1 ])

        self.r2_lookup = r2_full_table
        self.r2_match_lengths = r2_match_lengths

    def lookup_r1(self, seq):
        res = self.r1_lookup.get(seq)
        if not res and self.r1_aliases:
            keylist = self.r1_aliases.get(seq[:self.r1_lookup_length])
            if keylist:
                for key in keylist:
                    if key and seq.startswith(key):
                        res = self.r1_lookup.get(key)
                        break
        return res

    def lookup_r2(self, target_name, seq):
        lookup = self.r2_lookup[target_name]
        r2_match_len = self.r2_match_lengths[target_name]
        return lookup.get(seq[:r2_match_len])
