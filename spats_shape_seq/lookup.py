#
# lookup algorithm implementation
#
# general idea: there are only so many valid possibilities for R1 (w/o
# handle): so, just create a dictionary whose keys are the possible
# R1[4:], and whose values indicate the associated information (which
# target, length L of subsequence for cotrans, how much adapter to
# trim, etc). to process a pair, first lookup R1 in the dictionary,
# and then use a similar lookup on R2 (or else R2 is determined by the
# presence of adapter).
#
# cotrans detail:
#
# - possible R1 keys (i.e., R1[4:]) are of the form:
#     rc([target][linker]) + adapter_b
#   the total length is 32, so len(adapter) + len(target) must add up to len(R1) - 4 - len(linker)
#
# - xref target.py:build_cotrans_lookups, where the lookup table is
#   built.  basically, just iterate over all possible L, and then over
#   the valid amounts of target/adapter, to form the table
#
# - note that some keys are ambiguous: for example, if the target-part
#   starts with "T", then when there's no adapter, the key will end in
#   "A". but the len=1 adapter key will also end in "A", and otherwise
#   be identical (since adapter_b starts with "A"). so,
#   _process_pair_cotrans iterates through all of these "hits" and see
#   which one works.
#
# - in some cases the keys are ambiguous: the same key and R2 could be
#   at two different spots. that fails with the "multiple_R1" error
#   code.
#
# - once R1 is looked up, there are two possibilities for R2:
#
#   (a) it's fully in the sequence: then we look it up using the
#       r2_lookup table, and we're done
#
#   (b) it overlaps linker/adapter: in which case, the site it has to
#       start at is determined by R1. so then it only remains to check
#       that R2 has no errors.
#
# - that's pretty much it. in the common case, this is just two
#   dictionary lookups and a few string compares, so it's very
#   fast. see _try_lookup_hit() below for details.
#

from processor import PairProcessor, Failures
from util import _warn, _debug, reverse_complement, string_match_errors

class LookupProcessor(PairProcessor):

    def prepare(self):
        targets = self._targets
        targets.index()
        targets.build_lookups(self._run, length = self._run.pair_length)
        if not self._run.quiet:
            print("Lookup table: {} R1 entries, {} R2 entries for {} targets.".format(len(targets.r1_lookup),
                                                                                      sum(map(len, targets.r2_lookup.values())),
                                                                                      len(targets.r2_lookup)))


    #@profile
    def process_pair(self, pair):

        if not self._match_mask(pair):
            return

        r1_res = self._targets.r1_lookup.get(pair.r1.original_seq[4:])
        if not r1_res:
            pair.failure = Failures.nomatch
            return

        if len(r1_res) > 1:
            for trim in [ h[2] for h in r1_res ]:
                if len([ h[1] for h in r1_res if h[2] == trim]) > 1:
                    pair.failure = Failures.multiple_R1
                    return

        for hit in r1_res:
            self._try_lookup_hit(pair, hit)
            if pair.has_site:
                return

    #@profile
    def _try_lookup_hit(self, pair, r1_res):

        targets = self._targets
        run = self._run

        site = -1
        target = r1_res[0]
        r2len = pair.r2.original_len
        r2_mutations = []
        if r1_res[1] == None:
            lookup = targets.r2_lookup[target.name]
            r2_match_len = targets.r2_match_lengths[target.name]
            r2_key = pair.r2.original_seq[:r2_match_len]
            r2_res = lookup.get(r2_key)
            if r2_res is not None:
                match_site = r2_res[0]
            else:
                pair.failure = Failures.nomatch
                return
        else:
            match_site = r1_res[1]

        # need to check R2 against expectation
        match_len = min(r2len, target.n - match_site)
        if match_len > 0:
            pair.r2.match_errors = string_match_errors(pair.r2.original_seq[:match_len], target.seq[match_site:match_site + match_len])
            #+1 for M_j indexing convention, xref https://trello.com/c/2qIGo9ZR/201-stop-map-mutation-indexing-convention
            r2_mutations = map(lambda x : x + match_site + 1, pair.r2.match_errors)
        if match_len <= 0  or  len(pair.r2.match_errors) > run.allowed_target_errors:
            pair.failure = Failures.match_errors
            return

        adapter_len = r2len - match_len - 4
        if adapter_len > 0:
            pair.r2.adapter_errors = string_match_errors(self._adapter_t_rc[:adapter_len], pair.r2.original_seq[-adapter_len:])
            if len(pair.r2.adapter_errors) > run.allowed_adapter_errors:
                pair.failure = Failures.adapter_trim
                return
        site = match_site

        # in rare cases, need to double-check what R1 should be based on R2
        if site != r1_res[1]:
            match_len = min(pair.r1.original_len - 4, target.n - match_site)
            adapter_len = pair.r1.original_len - match_len - 4
            pair.r1.match_errors = string_match_errors(reverse_complement(target.seq[-match_len:]), pair.r1.original_seq[4:4+match_len])
            if len(pair.r1.match_errors) > run.allowed_target_errors:
                pair.failure = Failures.match_errors
                return
            if adapter_len > 0:
                pair.r1.adapter_errors = string_match_errors(pair.r1.original_seq[-adapter_len:], self._run.adapter_b[:adapter_len])
                if len(pair.r1.adapter_errors) > run.allowed_adapter_errors:
                    pair.failure = Failures.adapter_trim
                    return

        if not self._check_indeterminate(pair):
            return

        pair.r1.match_len = min(pair.r1.original_len - 4, target.n - match_site)
        pair.r1.match_index = r1_res[1] or (target.n - pair.r1.match_len)
        pair.r1._rtrim = r1_res[2]
        pair.r2.match_index = site
        pair.r2.match_len = match_len

        if not pair.check_overlap():
            pair.failure = Failures.r1_r2_overlap
            return

        if r2_mutations or r1_res[3]:
            pair.mutations = list(set(r2_mutations + r1_res[3]))
            self.counters.low_quality_muts += pair.check_mutation_quality(run.mutations_require_quality_score)
            if pair.mutations and len(pair.mutations) > run.allowed_target_errors:
                pair.failure = Failures.match_errors
                return

        pair.target = target
        pair.site = site
        pair.end = target.n
        pair.failure = None
        self.counters.register_count(pair)


class CotransLookupProcessor(PairProcessor):

    def prepare(self):
        targets = self._targets
        targets.index()
        targets.build_cotrans_lookups(self._run)
        target = targets.targets[0]
        self.r2_lookup = targets.r2_lookup[target.name]
        self.r2_match_len = targets.r2_match_lengths[target.name]
        if not self._run.quiet:
            print("Lookup table: {} R1 entries, {} R2 entries.".format(len(targets.r1_lookup),
                                                                       sum(map(len, targets.r2_lookup.values()))))

    #@profile
    def process_pair(self, pair):

        if not self._match_mask(pair):
            return

        r1_res = self._targets.r1_lookup.get(pair.r1.original_seq[4:])
        if not r1_res:
            pair.failure = Failures.nomatch
            return
        if r1_res[0][1] < self._run.cotrans_minimum_length:
            pair.failure = Failures.cotrans_min
            return
        if len(r1_res) > 1:
            for trim in [ h[2] for h in r1_res ]:
                if len([ h[1] for h in r1_res if h[2] == trim]) > 1:
                    pair.failure = Failures.multiple_R1
                    return
        for hit in r1_res:
            self._try_lookup_hit(pair, hit)
            if pair.has_site:
                return

    #@profile
    def _try_lookup_hit(self, pair, r1_res):

        run = self._run
        linker = run.cotrans_linker
        linker_len = len(linker)
        r2_seq = pair.r2.original_seq
        pair_len = pair.r2.original_len
        target = r1_res[0]
        tseq = target.seq
        L = r1_res[1]
        trim = r1_res[2]
        site = -1
        r2_mutations = []

        if 0 == trim:
            r2_match_len = self.r2_match_len
            r2_res = self.r2_lookup.get(r2_seq[:r2_match_len])
            if r2_res is not None:
                site = r2_res[0]
                r2_mutations = r2_res[1]
            else:
                pair.failure = Failures.nomatch
                return
        else:
            site = L - (pair_len - linker_len - 4) + trim

        # now need to verify R2 contents
        # R2: [target][linker][handle][adapter]

        target_match_len = min(pair_len, L - site)
        if target_match_len > 0:
            pair.r2.match_errors = string_match_errors(r2_seq[:target_match_len], tseq[site:site + target_match_len])
            #+1 for M_j indexing convention, xref https://trello.com/c/2qIGo9ZR/201-stop-map-mutation-indexing-convention
            r2_mutations = map(lambda x : x + site + 1, pair.r2.match_errors)
        if target_match_len <= 0  or  len(pair.r2.match_errors) > run.allowed_target_errors:
            pair.failure = Failures.match_errors
            return

        if target_match_len < pair_len:
            linker_match_len = min(linker_len, pair_len - target_match_len)
            if r2_seq[target_match_len:target_match_len + linker_match_len] != linker[:linker_match_len]:
                pair.failure = Failures.linker
                return

            if trim > 0 and r2_seq[-trim:] != self._adapter_t_rc[:trim]:
                pair.failure = Failures.adapter_trim
                return
            if pair_len - target_match_len - linker_len - 4 > trim:
                pair.failure = Failures.adapter_trim
                return

        if not self._check_indeterminate(pair):
            return

        pair.r1.match_index = L - (pair_len - linker_len - 4) + trim
        pair.r1._rtrim = trim
        pair.r2.match_index = site
        pair.r2.match_len = target_match_len
        if not pair.check_overlap():
            pair.failure = Failures.r1_r2_overlap
            return

        if r2_mutations or r1_res[3]:
            pair.mutations = list(set(r2_mutations + r1_res[3]))
            self.counters.low_quality_muts += pair.check_mutation_quality(run.mutations_require_quality_score)
            if pair.mutations and len(pair.mutations) > run.allowed_target_errors:
                pair.failure = Failures.match_errors
                return

        pair.end = L
        pair.target = target
        pair.site = site
        pair.failure = None
        self.counters.register_count(pair)
