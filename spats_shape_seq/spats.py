
# spats

# primary sequencing problem:
#
#  given R1-R2 read pair, determine the fragment location in the target sequence


# algorithm overview:
#
#  excepting pairs that correspond to very short fragments (which
#  could be handled as a special case, TBD), we can determine the
#  location of the, e.g., R2 subfragment by finding the longest
#  substring of R2 that matches some subsequence of the target. unless
#  there are many evenly-spaced transcription errors (in which case we
#  probably want to discard the pair anyway), we're guaranteed to find
#  some suitably-long (~8+ bp) matching subsequence. (this can be done
#  efficiently.)
#
#  once we have the matching R1/R2 subsequences, everything else about
#  the pair is determined (details/examples below). then we can
#  compare how many errors the input data have with the expectation,
#  and choose to keep/discard the fragment accordingly.
#
#  note that this will work regardless of how much adapter is included
#  in the fragment. (i.e., will work fine on case 2-I and 2-II of
#  previous adapter trimming algorithm.)


# conceptual algorithm:
#   (xref diagram examples below)

#  (a) find the longest subsequence of R2 that matches target T -- let
#      the results be s_2, l_2, i_2:
#         - s_2: the index in R2 where the match starts
#         - l_2: the length of the matched subsequence
#         - i_2: the index in T where the match starts
#      in other words, R2[s_2:s_2+l_2] == T[i_2:i_2+l_2]

#  (b) let a_2 be the index in T where the left end of R2 would end up
#      when the matching subsequences are aligned. in code,
#        a_2 = i_2 - s_2
#      note that if s_2 = 0 (the start of R2 is part of the matching 
#      subsequence), then a_2 = i_2. similarly, let b_2 correspond to 
#      the right end of R2, 
#        b_2 = a_2 + len(R2)
#      note that if all of R2 matches, then b_2 = i_2 + l_2
      
#  (c) let H = R1[:4], the first four bp of R1 (used for determining
#      the handle).  let R1' be the reverse complement of the rest of
#      R1, i.e., 
#        R1' = reverse_complement(R1[4:]). 
#      find the longest subsequence of R1' that matches T, and denote
#      s_1, l_1, i_1, a_1, and b_1 accordingly.

#  (d) if a_i > 0 and b_i <= len(T) for i=1,2: then both R1 and R2 lie
#      within T. let TR_i be the corresponding subsequences of T:
#        TR_2 = T[a_2:b_2]
#        TR_1 = reverse_complement(T[a_1:b_1])
#      then go to step [f]

#   (e) otherwise, some parts of R1 and R2 lie outside of T. if a_i <
#       0, then it's outside of T (off the left end), so should be
#       discarded. otherwise, at least one of the b_i > len(T), and
#       this corresponds to the adapter trimming case. set:
#          T_2 = T + reverse_complement(H) + reverse_complement(adapter_t)
#       (recall H = R1[:4]). let 
#         TR_2 = T_2[a_2:b_2]
#       for the R1 part, we take as much of T as we can, and reverse complement:
#          tmp = reverse_complement(T[a_1:min(b_1, len(T))])
#       and then append however much of adapter_b is necessary to get TR_1
#         TR_1 = tmp + adapter_b[0:len(R1')-len(tmp)]

#   (f) the number of errors can now be observed by comparing with the inputs, i.e.:
#        E_2 = |TR_2 - R2|, where |x-y| is hamming distance (number of mismatched bp)
#        E_1 = |TR_1 - R1|

#   (g) if the E_i are within configured range, increment the count
#       for handle H for a_2 (left end of R2), with "start" position
#       b_1 (right end of R1 -- for now, always len(T), but easy to
#       extend to multiple start sites).


# algorithm diagrammatic examples:


#  (Ex1) R2 = CCACCTGACCCCATGCCGAACTCAGAAGTGAAACG, R1 = AAACGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG
#    easy case: all of R1 and R2 match, and lie within T:
#
#                                                                                                                RRRY               R1
#                                                                                                                AAAC.GTCCTTGGTGCCCGAGTCAGATGCCTGGCAG
#                                                   R2                                                                |          (revcomp)          |
#                                  CCACCTGACCCCATGCCGAACTCAGAAGTGAAACG                                                CTGCCAGGCATCTGACTCGGGCACCAAGGAC 
# T   ggatgcctggcggccgtagcgcggtggtcCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGccgtagcgccgatggtagtgtggggtctccccatgcgagagtagggaaCTGCCAGGCATCTGACTCGGGCACCAAGGAC       
#                                  ^--   L = 29                                                                                           R = 143 --^
#
#  so we count site L = a_2 = 29 for handle RRRY.
#

#  (Ex2) R2 = CCAGCTGTCCCCATGCCGAAGTCAGAAATGAAACG, R1 = AAACCTGTCAGGCATCTGACACGGGCACAAAGGAC
#    similar to (1) but now some bp's are toggled, so we only match parts of R2 and R1':
#  
#
#                                          s        R2                                                                    s      R1'
#                                  --------CCCCATGCCGAA---------------                                                ----CAGGCATCTGAC--------------- 
# T   ggatgcctggcggccgtagcgcggtggtcccacctgaCCCCATGCCGAActcagaagtgaaacgccgtagcgccgatggtagtgtggggtctccccatgcgagagtagggaactgcCAGGCATCTGACtcgggcaccaaggac       
#                                  a       i                         b                                                a   i                         b
#
# we've found s_i, l_i, i_i, a_i, b_i, and this is case (d) above, so we can create TR2, TR1 based on T:
#                           TR_2:  CCACCTGACCCCATGCCGAACTCAGAAGTGAAACG                                         TR_1:  CTGCCAGGCATCTGACTCGGGCACCAAGGAC
# and compare:                        !   !            !      !                                                          !            !       !   
#                           R2:    CCAGCTGTCCCCATGCCGAAGTCAGAAATGAAACG                                          R1':  CTGTCAGGCATCTGACACGGGCACAAAGGAC
#
# E_1 = 3, E_2 = 4, and we can count this match (still at L = a_2 = 29, R = b_1 = 143) based on configuration. 
#


#  (Ex3) R2 = CTCGAGCACCAAGGACTGAAAGCTCGGAAGAGCGA, R1 = TTCAGTTCTTGGTGCCCGAGTGATCGGAAGAGCAC       
#    on adapter trimming, need to build up TR_1, TR_2 with adapter (and handle complement). first, the matching substrings:
#            
#                                                   R2  :              -----GCACCAAGGAC-------------------
#    T      ggatgcctggcggc..tggggtctccccatgcgagagtagggaactgccaggcatctgaCTCGGGCACCAAGGAC
#                                                   R1' :              CTCGGGCACCAAG-----------------
#
# now, to build T_2, we append RC(H) and RC(adapter_t):
#   H = TTCA, RC(H) = TGAA, adapter_t = AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
#   RC(adapter_t) = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#  T_2 = GGATGCCTGGCGGC..TGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCTGACTCGGGCACCAAGGAC + TGAA + AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#   so, to get TR_2:
#                                                        R2  :      -----GCACCAAGGAC-------------------
#   T_2: ggatgcctggcggc..tggggtctccccatgcgagagtagggaactgccaggcatctgactcggGCACCAAGGACtgaaagatcggaagagcgtcgtgtagggaaagagtgtagatctcggtggtcgccgtatcatt
#                                                      TR_2  :      CTCGGGCACCAAGGACTGAAAGATCGGAAGAGCGT
# and compare to R2:                                                    !                 !           !
#                                                                   CTCGAGCACCAAGGACTGAAAGCTCGGAAGAGCGA
# that gives E_2 = 3 (1 error in the target, 2 in the adapter).
# now, for TR_1, first we need to compare R1' to T_1 = T + adapter_b = T + AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
#                                                        R1' :      CTCGGGCACCAAG------------------
#   T_1: ggatgcctggcggc..tggggtctccccatgcgagagtagggaactgccaggcatctgaCTCGGGCACCAAGgac 
#                                                                   CTCGGGCACCAAGGAC
# and then take the reverse complement                              GTCCTTGGTGCCCGAG
# and then append adapter_b:                                        GTCCTTGGTGCCCGAGagatcggaagagcacacgtctgaactccagtcac
# and trim to get                                       TR_1 :      GTCCTTGGTGCCCGAGAGATCGGAAGAGCAC
# and compare to R1:                                                  !             !              
#                                                              TTCA.GTTCTTGGTGCCCGAGTGATCGGAAGAGCAC             
#                                                                   GTGCTCTTCCGATCACTCGGGCACCAAGAAC
# which gives E_2 = 2 (1 error in the target, 1 in the adapter).
# the match site is at L = a_2 = 127, R = b_1 = 143



# notes:
#
# the implementation below does not follow the algorithm literally,
# but effectively comes to the same thing. it allows for separate
# configuration of allowable errors within the target and within the
# adapter. also, right now errors on the part of T_2 that corresponds
# to reverse_complement(H) in the trimming case are ignored (easily
# configurable).
#
# the algorithm is implemented in the _process_pair method of the
# Spats class.
#
# the key to the algorithm is step (a), which is implemented by the
# Targets class, in particular Targets.find_partial
#
# s_i, l_i, i_i are referred to as match_start, match_len, and
# match_index in the code (in Sequence class). a_i and b_i are
# referred to as left and right.
#
# right now, the code is a bit more cluttered than it should be, in an
# attempt to conform as much as possible to v1.0.2 behavior.
#


import multiprocessing
import Queue
import time

from config import spats_config
from db import PairDB
from mask import Mask, longest_match
from pair import Pair
from parse import FastFastqParser, fasta_parse
from profiles import Profiles
from target import Targets
from util import _warn, _debug, Counters, reverse_complement, string_match_errors


class Spats(object):

    def __init__(self):
        self._targets = Targets()
        self._masks = []
        self._maskSize = 0
        self._adapter_t_rc = 0
        self._profiles = None
        self.counters = Counters()
        self.writeback_results = False

    def addMasks(self, *args):
        for arg in args:
            mask = Mask(arg)
            self._masks.append(mask)
            self._maskSize = max(self._maskSize, len(arg))

    def addTargets(self, *target_paths):
        target = self._targets
        for path in target_paths:
            for name, seq in fasta_parse(path):
                target.addTarget(name, seq)
        if not target.targets:
            raise Exception("didn't get any targets!")
        if spats_config.minimum_target_match_length:
            target.minimum_match_length = spats_config.minimum_target_match_length
            target.index()
        else:
            target.index()
            target.minimum_match_length = 1 + target.longest_self_match()

    @property
    def adapter_t_rc(self):
        if not self._adapter_t_rc:
            self._adapter_t_rc = reverse_complement(spats_config.adapter_t)
        return self._adapter_t_rc

    def _match_mask(self, pair):
        seq = pair.r1.original_seq
        if len(seq) <= self._maskSize:
            _warn("sequence {} is too short for masking".format(seq))
            return
        for mask in self._masks:
            if mask.matches(seq):
                pair.set_mask(mask)
                return

    def _find_matches(self, pair):
        # use R1 to determine which target
        target = pair.r1.find_in_targets(self._targets, reverse_complement = True)
        if isinstance(target, list):
            _debug("dropping pair due to multiple R1 match")
            pair.failure = "multiple R1 match"
            self.counters.multiple_R1_match += pair.multiplicity
            return
        if target:
            pair.target = pair.r2.find_in_targets(self._targets, force_target = target)
        _debug([pair.r1.match_start, pair.r1.match_len, pair.r1.match_index, "--", pair.r2.match_start, pair.r2.match_len, pair.r2.match_index])

    def _trim_adapters(self, pair):

        # if we're here, we know that R2 hangs over the right edge. first, figure out by how much

        r2_seq = pair.r2.subsequence
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        r2_len_should_be = pair.target.n - r2_start_in_target
        r2_length_to_trim = pair.r2.seq_len - r2_len_should_be
        pair.r2.trim(r2_length_to_trim)
        _debug("R2 trim: {}, {}, {}".format(r2_start_in_target, r2_len_should_be, r2_length_to_trim))

        if spats_config.minimum_adapter_len and r2_length_to_trim - 4 < spats_config.minimum_adapter_len:
            _debug("  !! v102 minimum adapter len {}".format(r2_length_to_trim - 4))
            return False

        if r2_length_to_trim <= 4:
            # TODO: should we verify that this matches RC of R1 handle, and register errors for bp that don't?
            # for now, just ignore this part
            # also, this means that there's nothing to trim from R1, so we're done
            return True

        # find out how good of a match the end of R2 is for adapter_t_rc
        r2_adapter_match = r2_seq[4-r2_length_to_trim:]
        pair.r2.adapter_errors = string_match_errors(r2_adapter_match, self.adapter_t_rc)
        _debug("  check = {}, errors = {}".format(r2_adapter_match, pair.r2.adapter_errors))
        if len(pair.r2.adapter_errors) > spats_config.allowed_adapter_errors:
            return False

        # now, same thing on r1 (smaller trim b/c of no handle, hence -4)
        r1_seq = pair.r1.subsequence
        r1_length_to_trim = r2_length_to_trim - 4
        r1_adapter_match = r1_seq[-r1_length_to_trim:]
        r1_match_trimmed = pair.r1.trim(r1_length_to_trim, reverse_complement = True)
        pair.r1.adapter_errors = string_match_errors(r1_adapter_match, spats_config.adapter_b)
        _debug("  R1 check = {}, errors = {}".format(r1_adapter_match, pair.r1.adapter_errors))
        if len(pair.r1.adapter_errors) > spats_config.allowed_adapter_errors:
            return False

        if r1_match_trimmed:
            # ok, we trimmed down our R1 due to adapters. need to see if that means the leftover matches
            # multiple targets; if so, need to reject this pair.
            target = pair.r1.find_in_targets(self._targets, reverse_complement = True, min_length_override = pair.r1.seq_len)
            if isinstance(target, list):
                _debug("dropping pair due to multiple R1 match after adapter trim")
                pair.target = None
                pair.failure = "multiple R1 match"
                self.counters.multiple_R1_match += pair.multiplicity
                return False

        _debug("successful adapter trim of {}/{} bp from R1/R2".format(pair.r1._rtrim, pair.r2._rtrim))

        return True

    def process_pair(self, pair):
        _debug("> processing " + pair.identifier + "\n  --> " + pair.r1.original_seq + " , " + pair.r2.original_seq)
        _debug("  rc(R1): {}".format(pair.r1.reverse_complement))
        self._process_pair(pair)
        if pair.failure:
            _debug(pair.failure)
        else:
            assert(pair.has_site)
            _debug("  ===> KEPT {}-{}".format(pair.left, pair.right))
            if spats_config.show_id_to_site:
                print "{} --> {}".format(pair.identifier, pair.site) #, pair.mask.chars)

    #@profile
    def _process_pair(self, pair):

        if not spats_config.allow_indeterminate  and  not pair.is_determinate():
            pair.failure = "indeterminate sequence failure"
            self.counters.indeterminate += pair.multiplicity
            return

        self._match_mask(pair)
        if not pair.mask:
            self.counters.mask_failure += pair.multiplicity
            pair.failure = "mask failure"
            return

        self._find_matches(pair)
        if not pair.matched:
            self.counters.unmatched += pair.multiplicity
            pair.failure = "no match"
            return

        # this is where R2 should start (if not a complete match, then r2.match_start will be > 0)
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        if r2_start_in_target < 0:
            pair.failure = "R2 to left of site 0 failure on R2 for: {}".format(pair.identifier)
            self.counters.left_of_target += pair.multiplicity
            return
        elif r2_start_in_target + pair.r2.original_len <= pair.target.n:
            # we're in the middle -- no problem
            pass
        elif not self._trim_adapters(pair):
            # we're over the right edge, and adapter trimming failed
            pair.failure = pair.failure or "adapter trim failure"
            self.counters.adapter_trim_failure += pair.multiplicity
            return
        else:
            # we're at the right and trimming went ok, cool
            self.counters.adapter_trimmed += pair.multiplicity

        if pair.r2.match_len == pair.r2.seq_len  and  pair.r1.match_len == pair.r1.seq_len:
            # everything that wasn't adapter trimmed was matched -- nothing to do
            pass
        else:
            # set the match to be the rest of the (possibly trimmed) sequence, and count errors
            pair.r1.match_to_seq(reverse_complement = True)
            pair.r2.match_to_seq()
            target_seq = pair.target.seq
            pair.r1.match_errors = string_match_errors(pair.r1.reverse_complement, target_seq[pair.r1.match_index:])
            pair.r2.match_errors = string_match_errors(pair.r2.subsequence, target_seq[pair.r2.match_index:])

            if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > spats_config.allowed_target_errors:
                if pair.r1.match_errors:
                    _debug("R1 errors: {}".format(pair.r1.match_errors))
                if pair.r2.match_errors:
                    _debug("R2 errors: {}".format(pair.r2.match_errors))
                pair.failure = "match errors failure"
                self.counters.match_errors += pair.multiplicity
                return

        n = pair.target.n
        assert(pair.matched and pair.left >= 0 and pair.left <= n)

        # NOTE: this might change later due to "starts"
        if pair.right != n:
            pair.failure = "R1 right edge failure: {} - {}, n={}".format(pair.left, pair.right, n)
            self.counters.r1_not_on_right_edge += pair.multiplicity
            return

        pair.register_count()
        self.counters.processed_pairs += pair.multiplicity

    def memory_db_from_pairs(self, data_r1_path, data_r2_path):
        start = time.time()
        db = PairDB()
        total_pairs = db.parse(data_r1_path, data_r2_path)
        report = "Parsed {} records in {:.1f}s".format(total_pairs, time.time() - start)

        # unclear if this helps, but potentially useful for further analysis later, and doesn't cost *too* much
        # but if it's holding this up, nuke it
        db.index()
        report += ", indexed in {:.1f}s".format(time.time() - start)

        if spats_config.quiet:
            _debug(report)
        else:
            print report
        return db

    def process_pair_data(self, data_r1_path, data_r2_path):
        self.process_pair_db(self.memory_db_from_pairs(data_r1_path, data_r2_path))
    
    #@profile
    def process_pair_db(self, pair_db):

        start = time.time()
        writeback = self.writeback_results
        if writeback:
            pair_db.prepare_results()

        pairs_to_do = multiprocessing.Queue()

        #@profile
        def worker(x):
            pair = Pair()
            writeback = self.writeback_results
            while True:
                pairs = pairs_to_do.get()
                if not pairs:
                    break
                results = []
                for lines in pairs:
                    pair.set_from_data('', lines[1], lines[2], lines[0])
                    self.process_pair(pair)
                    if writeback:
                        results.append((lines[3],
                                        pair.target.name if pair.target else None,
                                        pair.site if pair.has_site else -1,
                                        pair.mask.chars if pair.mask else None,
                                        pair.multiplicity,
                                        pair.failure))
                                        
                if writeback:
                    pair_db.add_results(results)

            pairs_to_do.put((self.counters._counts, [(m.total, m.kept, m.count_data()) for m in self._masks]))

        threads = []
        num_workers = max(1, spats_config.num_workers)
        for i in range(num_workers):
            thd = multiprocessing.Process(target = worker, args = (i,))
            threads.append(thd)
            thd.start()
        _debug("created {} workers".format(num_workers))

        db_fn = pair_db.all_pairs if writeback else pair_db.unique_pairs_with_counts
        for pair_info in db_fn(batch_size = 16384):
            if not writeback:
                self.counters.unique_pairs += len(pair_info)
            pairs_to_do.put(pair_info)

        for thd in threads:
            pairs_to_do.put(None) # just a dummy object to signal we're done
        for thd in threads:
            thd.join()

        try:
            targets = { t.name : t for t in self._targets.targets }
            while True:
                data = pairs_to_do.get_nowait()
                their_counters = data[0]
                our_counters = self.counters._counts
                for key, value in their_counters.iteritems():
                    if key != "_counts":
                        our_counters[key] = our_counters.get(key, 0) + value
                for i in range(len(self._masks)):
                    m = self._masks[i]
                    d = data[1][i]
                    m.total += d[0]
                    m.kept += d[1]
                    m.update_with_count_data(d[2], targets)
        except Queue.Empty:
            pass

        self.counters.total_pairs += pair_db.count()
        if writeback:
            self.counters.unique_pairs = pair_db.unique_pairs()

        if not spats_config.quiet:
            self.report_counts(time.time() - start)


    def report_counts(self, delta = None):
        total = self.counters.total_pairs
        print "Successfully processed {} properly paired fragments:".format(self.counters.processed_pairs)
        warn_keys = [ "multiple_R1_match", ]
        countinfo = self.counters.counts_dict()
        for key in sorted(countinfo.keys(), key = lambda k : countinfo[k], reverse = True):
            print "  {}{} : {} ({:.1f}%)".format("*** " if key in warn_keys else "", key, countinfo[key], 100.0 * (float(countinfo[key])/float(total)) if total else 0)
        print "Masks:"
        tmap = { t.name : 0 for t in self._targets.targets }
        for m in self._masks:
            print "  {}: kept {}/{} ({:.1f}%)".format(m.chars, m.kept, m.total, (100.0 * float(m.kept)) / float(m.total) if m.total else 0)
        for t in self._targets.targets:
            tmap[t.name] += sum(m.counts(t))
        if 1 < len(self._targets.targets):
            print "Targets:"
            total = self.counters.processed_pairs
            for tgt in sorted(self._targets.targets, key = lambda t : tmap[t.name], reverse = True):
                if tmap[tgt.name] > 0:
                    print "  {}: {} ({:.1f}%)".format(tgt.name, tmap[tgt.name], (100.0 * float(tmap[tgt.name])) / float(total) if total else 0)
        if delta:
            print "Total time: ({:.1f}s)".format(delta)

    def compute_profiles(self):
        self._profiles = Profiles(self._targets, self._masks)
        self._profiles.compute()
        return self._profiles

    def write_reactivities(self, target_path):
        self._profiles.write(target_path)
