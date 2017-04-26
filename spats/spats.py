
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
# Target class, in particular Target.find_partial
#
# s_i, l_i, i_i are referred to as match_start, match_len, and
# match_index in the code (in Sequence class). a_i and b_i are
# referred to as left and right.
#
# right now, the code is a bit more cluttered than it should be, in an
# attempt to conform as much as possible to v1.0.2 behavior.
#


import math
import multiprocessing
import os
import Queue
import sys
import time

from config import spats_config
from mask import Mask, longest_match
from pair import Pair
from parse import fasta_parse
from target import Target
from util import _warn, _debug, reverse_complement, string_match_errors


class Spats(object):

    def __init__(self, target_path, output_folder):
        self.target_path = target_path
        self.output_folder = output_folder

        # user-configurable parameters
        self.masks = [ 'RRRY', 'YYYR' ]
        self.quiet = False
        self.adapter_t = spats_config.adapter_t
        self.adapter_b = spats_config.adapter_b

        # private vars
        self._targets = None
        self._masks = None

        self.total_pairs = self.processed_pairs = self.chucked_pairs = 0

    def setup(self):
        self._indexTargets()
        self._setupMasks()
        self.adapter_t_rc = reverse_complement(self.adapter_t)

    def _indexTargets(self):
        if self._targets:
            return
        targets = []
        for name, seq in fasta_parse(self.target_path):
            target = Target(name, seq)
            target.index()
            targets.append(target)
        self._targets = targets

        # TODO: handle multiple targets. for now just use the one
        self._target = targets[0]

    def _setupMasks(self):
        if self._masks:
            return
        self._masks = [ Mask(m, self._target.n) for m in self.masks ]
        self._maskSize = max([ len(m.chars) for m in self._masks ])

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
        pair.r1.find_in_target(self._target, reverse_complement = True)
        pair.r2.find_in_target(self._target)
        _debug([pair.r1.match_start, pair.r1.match_len, pair.r1.match_index, "--", pair.r2.match_start, pair.r2.match_len, pair.r2.match_index])

    def _trim_adapters(self, pair):

        # if we're here, we know that R2 hangs over the right edge. first, figure out by how much

        r2_seq = pair.r2.subsequence
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        r2_len_should_be = self._target.n - r2_start_in_target
        r2_length_to_trim = pair.r2.seq_len - r2_len_should_be
        pair.r2.trim(r2_length_to_trim)
        _debug("R2 trim: {}, {}, {}".format(r2_start_in_target, r2_len_should_be, r2_length_to_trim))

        if spats_config.minimum_adapter_len and r2_length_to_trim - 4 < spats_config.minimum_adapter_len:
            _debug("  !! v102 minimum adapter len {}".format(r2_length_to_trim - 4))
            return False

        if r2_length_to_trim < 4:
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
        pair.r1.trim(r1_length_to_trim, reverse_complement = True)
        pair.r1.adapter_errors = string_match_errors(r1_adapter_match, self.adapter_b)
        _debug("  R1 check = {}, errors = {}".format(r1_adapter_match, pair.r1.adapter_errors))
        if len(pair.r1.adapter_errors) > spats_config.allowed_adapter_errors:
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
            return

        self._match_mask(pair)
        if not pair.mask:
            pair.failure = "mask failure"
            return

        self._find_matches(pair)
        if not pair.matched:
            pair.failure = "no match"
            return

        # this is where R2 should start (if not a complete match, then r2.match_start will be > 0)
        r2_start_in_target = pair.r2.match_index - pair.r2.match_start
        if r2_start_in_target < 0:
            pair.failure = "R2 to left of site 0 failure on R2 for: {}".format(pair.identifier)
            return
        elif r2_start_in_target + pair.r2.original_len <= self._target.n:
            # we're in the middle -- no problem
            pass
        elif not self._trim_adapters(pair):
            # we're over the right edge, and adapter trimming failed
            pair.failure = pair.failure or "adapter trim failure"
            return
        else:
            # we're at the right and trimming went ok, cool
            pass

        if pair.r2.match_len == pair.r2.seq_len  and  pair.r1.match_len == pair.r1.seq_len:
            # everything that wasn't adapter trimmed was matched -- nothing to do
            pass
        else:
            # set the match to be the rest of the (possibly trimmed) sequence, and count errors
            pair.r1.match_to_seq(reverse_complement = True)
            pair.r2.match_to_seq()
            target = self._target.seq
            pair.r1.match_errors = string_match_errors(pair.r1.reverse_complement, target[pair.r1.match_index:])
            pair.r2.match_errors = string_match_errors(pair.r2.subsequence, target[pair.r2.match_index:])

            if max(len(pair.r1.match_errors), len(pair.r2.match_errors)) > spats_config.allowed_target_errors:
                if pair.r1.match_errors:
                    _debug("R1 errors: {}".format(pair.r1.match_errors))
                if pair.r2.match_errors:
                    _debug("R2 errors: {}".format(pair.r2.match_errors))
                pair.failure = "match errors failure"
                return

        n = self._target.n
        assert(pair.matched and pair.left >= 0 and pair.left <= n)

        # NOTE: this might change later due to "starts"
        if pair.right != n:
            pair.failure = "R1 right edge failure: {} - {}, n={}".format(pair.left, pair.right, n)
            return

        pair.register_count()

    #@profile
    def process_pair_data(self, data_r1_path, data_r2_path, max_pairs = 0):

        pairs_to_do = multiprocessing.Queue()

        #@profile
        def worker(x):
            chucked = 0
            processed = 0
            pair = Pair()
            while True:
                pairs = pairs_to_do.get()
                if not pairs:
                    break
                for lines in pairs:
                    pair.set_from_data(lines[0].split(' ')[0], lines[1].rstrip('\n'), lines[2].rstrip('\n'))
                    self.process_pair(pair)
                    if not pair.mask:
                        chucked += 1
                    elif pair.has_site:
                        processed +=1
            pairs_to_do.put((chucked, processed, [(m.total, m.kept, m.counts) for m in self._masks]))

        threads = []
        num_workers = max(1, spats_config.num_workers)
        for i in range(num_workers):
            thd = multiprocessing.Process(target = worker, args = (i,))
            threads.append(thd)
            thd.start()
        _debug("created {} workers".format(num_workers))

        total_pairs = 0
        start = time.time()

        with open(data_r1_path, 'rb') as r1_in:
            with open(data_r2_path, 'rb') as r2_in:
                r1_iter = iter(r1_in)
                r2_iter = iter(r2_in)

                while True:
                    pairs = []
                    for i in range(16384):
                        if max_pairs and total_pairs >= max_pairs:
                            break
                        try:
                            R1_id = r1_iter.next()
                            R1_seq = r1_iter.next()
                            r1_iter.next()
                            r1_iter.next()
                            r2_iter.next()
                            R2_seq = r2_iter.next()
                            r2_iter.next()
                            r2_iter.next()
                            pairs.append((R1_id, R1_seq, R2_seq))
                            total_pairs += 1
                        except StopIteration:
                            break
                    if 0 == len(pairs):
                        break
                    pairs_to_do.put(pairs)


        if self.quiet:
            _debug("Parsed {} records in {:.1f}s".format(total_pairs, time.time() - start))
        else:
            print "Parsed {} records in {:.1f}s".format(total_pairs, time.time() - start)
            

        for thd in threads:
            pairs_to_do.put(None) # just a dummy object to signal we're done
        for thd in threads:
            thd.join()

        try:
            while True:
                data = pairs_to_do.get_nowait()
                self.chucked_pairs += data[0]
                self.processed_pairs += data[1]
                for i in range(len(self._masks)):
                    m = self._masks[i]
                    d = data[2][i]
                    m.total += d[0]
                    m.kept += d[1]
                    counts = d[2]
                    for j in range(len(counts)):
                        m.counts[j] += counts[j]
        except Queue.Empty:
            pass
        self.total_pairs += total_pairs

        if not self.quiet:
            self.report_counts(time.time() - start)

    def report_counts(self, delta = None):
        m0 = self._masks[0]
        m1 = self._masks[1]
        format_str = "Processed {tot} properly paired fragments, " + \
                     "kept {kept0}/{tot0} ({pct0:.1f}%) treated, " + \
                     "{kept1}/{tot1} ({pct1:1f}%) untreated"
        if delta:
            format_str += " ({:.1f}s)".format(delta)
        print format_str.format(tot = m0.total + m1.total,
                                kept0 = m0.kept,
                                tot0 = m0.total,
                                pct0 = (100.0 * float(m0.kept)) / float(m0.total),
                                kept1 = m1.kept,
                                tot1 = m1.total,
                                pct1 = (100.0 * float(m1.kept)) / float(m1.total))

    def compute_profiles(self):
        # TODO: use numpy here ?
        n = self._target.n
        treated_counts = self._masks[0].counts
        untreated_counts = self._masks[1].counts
        betas = [ 0 for x in range(n+1) ]
        thetas = [ 0 for x in range(n+1) ]
        treated_sum = 0.0    # keep a running sum
        untreated_sum = 0.0  # for both channels
        running_c_sum = 0.0  # can also do it for c

        for k in range(n):
            X_k = float(treated_counts[k])
            Y_k = float(untreated_counts[k])
            treated_sum += X_k    #treated_sum = float(sum(treated_counts[:(k + 1)]))
            untreated_sum += Y_k  #untreated_sum = float(sum(untreated_counts[:(k + 1)]))
            if 0 == treated_sum  or  0 == untreated_sum:
                betas[k] = 0
                thetas[k] = 0
            else:
                Xbit = (X_k / treated_sum)
                Ybit = (Y_k / untreated_sum)
                if Ybit >= 1:
                    betas[k] = 0
                    thetas[k] = 0
                else:
                    betas[k] = max(0, (Xbit - Ybit) / (1 - Ybit))
                    thetas[k] = math.log(1.0 - Ybit) - math.log(1.0 - Xbit)
                    running_c_sum -= math.log(1.0 - betas[k])

        c = running_c_sum
        c_factor = 1.0 / c
        for k in range(n+1):
            thetas[k] = max(c_factor * thetas[k], 0)
        self.betas = betas
        self.thetas = thetas
        self.c = c

    def write_reactivities(self):
        out_path = os.path.join(self.output_folder, "rx.out")
        n = self._target.n
        treated_counts = self._masks[0].counts
        untreated_counts = self._masks[1].counts
        with open(out_path, 'wb') as outfile:
            outfile.write('sequence\trt_start\tfive_prime_offset\tnucleotide\ttreated_mods\tuntreated_mods\tbeta\ttheta\tc\n')
            format_str = "{name}\t{rt}\t".format(name = self._target.name, rt = n - 1) + "{i}\t{nuc}\t{tm}\t{um}\t{b}\t{th}" + "\t{c:.5f}\n".format(c = self.c)
            # TODO: xref https://trello.com/c/OtbxyiYt/23-3-nt-missing-from-reactivities-out
            # looks like we may want this to be range(n), chopping off was unintentional bug of previous version
            for i in range(n - 1):
                outfile.write(format_str.format(i = i,
                                                nuc = self._target.seq[i - 1] if i > 0 else '*',
                                                tm = treated_counts[i],
                                                um = untreated_counts[i],
                                                b = self.betas[i] if i > 0 else '-',
                                                th = self.thetas[i] if i > 0 else '-'))
