
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


import os
import time

from db import PairDB
from mask import Mask
from pair import Pair
from parse import fasta_parse, FastFastqParser
from processor import PairProcessor
from profiles import Profiles
from run import Run
from target import Targets
from util import _debug, _set_debug, _warn
from worker import SpatsWorker


class Spats(object):
    """The main SPATS driver.

       :param cotrans: pass `True` for cotrans-style experiments.

    """

    def __init__(self, cotrans = False):
        self.run = Run()
        self.run.cotrans = cotrans
        self.__processor = None
        self._targets = None
        self._masks = None
        self._profiles = None


    @property
    def _processor(self):
        if not self.__processor:
            self._addMasks()
            self.__processor = self.run._get_processor_class()(self.run, self._targets, self._masks)
        return self.__processor

    @property
    def targets(self):
        return self._targets

    def _addMasks(self):
        if not self._masks:
            self._masks = [ Mask(m) for m in self.run.masks ]


    def addTargets(self, *target_paths):
        """Used to add one or more target files for processing. Can be called multiple 
           times to add more targets. Inputs are expected to be in FASTA format with
           one or more targets per path. Must be called before processing.

           :param args: one or more filesystem paths to target files.
        """
        targets = []
        for path in target_paths:
            for name, seq in fasta_parse(path):
                targets.append((name, seq, 1 + len(targets)))
        self._addTargets(targets)

    def addTarget(self, name, seq, rowid = -1):
        self._addTargets( [ (name, seq, rowid if rowid != -1 else len(self._targets.targets)) ] )

    def loadTargets(self, pair_db):
        self._addTargets(pair_db.targets())

    def _addTargets(self, target_list):
        targets = self._targets or Targets()
        for name, seq, rowid in target_list:
            targets.addTarget(name, seq, rowid)
        if not targets.targets:
            raise Exception("didn't get any targets!")
        targets.minimum_match_length = self.run.minimum_target_match_length
        self._targets = targets


    def process_pair(self, pair):
        """Used process a single :class:`.pair.Pair`. Typically only used for debugging or analysis of specific cases.

           :param pair: a :class:`.pair.Pair` to process.
        """

        if not self.run.pair_length:
            self.run.pair_length = len(pair.r1.original_seq)

        _set_debug(self.run)
        _debug("> processing " + pair.identifier + "\n  --> " + pair.r1.original_seq + " , " + pair.r2.original_seq)
        _debug("  rc(R1): {}".format(pair.r1.reverse_complement))
        try:
            self._processor.process_pair(pair)
            if pair.failure:
                _debug(pair.failure)
            else:
                assert(pair.has_site)
                _debug("  ===> KEPT {}-{}".format(pair.site, pair.end))
        except:
            print("**** Error processing pair: {} / {}".format(pair.r1.original_seq, pair.r2.original_seq))
            raise


    def _memory_db_from_pairs(self, data_r1_path, data_r2_path):
        if not self.run.quiet:
            print("Parsing pair data...")
        start = time.time()
        db = PairDB()
        total_pairs = db.parse(data_r1_path, data_r2_path)
        report = "Parsed {} records in {:.1f}s".format(total_pairs, time.time() - start)

        # unclear if this helps, but potentially useful for further analysis later, and doesn't cost *too* much
        # but if it's holding things up, nuke it
        db.index()
        report += ", indexed in {:.1f}s".format(time.time() - start)

        if self.run.quiet:
            _debug(report)
        else:
            print(report)
        return db


    def process_pair_data(self, data_r1_path, data_r2_path):
        """Used to read and process a pair of FASTQ data files.

        Note that this parses the pair data into an in-memory SQLite
        database, which on most modern systems will be fine except for the
        largest input sets. If you hit memory issues, create a disk-based
        SQLite DB via :class:`.db.PairDB` and then use
        :meth:`.process_pair_db`.

        Note that this may be called multiple times to process more
        than one set of data files before computing profiles.

        :param data_r1_path: path to R1 fragments
        :param data_r2_path: path to matching R2 fragments.
        """
        self.run.apply_config_restrictions()
        use_quality = self.run._parse_quality
        if not self.run.skip_database and not use_quality:
            self.process_pair_db(self._memory_db_from_pairs(data_r1_path, data_r2_path))
        else:
            with FastFastqParser(data_r1_path, data_r2_path, use_quality) as parser:
                if not self.run.pair_length:
                    self.run.pair_length = parser.pair_length()
                self._process_pair_iter(parser.iterator(batch_size = 131072))

    def process_pair_db(self, pair_db, batch_size = 65536):
        """Processes pair data provided by a :class:`.db.PairDB`.

           Note that this may be called multiple times to process more
           than one set of inputs before computing profiles.

           :param pair_db: a :class:`.db.PairDB` of pairs to process.
        """

        self.run.apply_config_restrictions()
        if not self.run.pair_length:
            self.run.pair_length = pair_db.pair_length()

        if not self._targets:
            self.loadTargets(pair_db)

        result_set_id = pair_db.add_result_set(self.run.result_set_name or "default", self.run.resume_processing) if self.run.writeback_results else None
        if self._processor.uses_tags:
            self._processor.setup_tags(pair_db)

        if self.run.resume_processing:
            db_iter = pair_db.unique_pairs_with_counts_and_no_results(result_set_id, batch_size = batch_size)
        elif self.run._process_all_pairs:
            if not self.run.quiet:
                print("Using all_pairs...")
            db_iter = pair_db.all_pairs(batch_size = batch_size)
        else:
            db_iter = pair_db.unique_pairs_with_counts(batch_size = batch_size)

        self._process_pair_iter(db_iter, pair_db, result_set_id)

    #@profile
    def _process_pair_iter(self, pair_iter, pair_db = None, result_set_id = None):

        _set_debug(self.run)

        start = time.time()

        # force the processor to load and do whatever indexing/etc is required
        self._processor

        worker = SpatsWorker(self.run, self._processor, pair_db, result_set_id)

        if not self.run.quiet:
            print("Processing pairs...")

        worker.run(pair_iter)

        if not self.run.quiet:
            self._report_counts(time.time() - start)


    def _report_counts(self, delta = None):
        counters = self.counters
        total = counters.total_pairs
        print("Successfully processed {} properly paired fragments:".format(counters.registered_pairs))
        warn_keys = [ "multiple_R1_match", ]
        countinfo = counters.counts_dict()
        for key in sorted(countinfo.keys(), key = lambda k : countinfo[k], reverse = True):
            print("  {}{} : {} ({:.1f}%)".format("*** " if key in warn_keys else "", key, countinfo[key], 100.0 * (float(countinfo[key])/float(total)) if total else 0))
        print("Masks:")
        for m in self._masks:
            kept, total = counters.mask_kept(m), counters.mask_total(m)
            print("  {}: kept {}/{} ({:.1f}%)".format(m.chars, kept, total, (100.0 * float(kept)) / float(total) if total else 0))
        if 1 < len(self._targets.targets):
            print("Targets:")
            tmap = { t.name : counters.target_total(t) for t in self._targets.targets }
            total = counters.registered_pairs
            for tgt in sorted(self._targets.targets, key = lambda t : tmap[t.name], reverse = True):
                if tmap[tgt.name] > 0:
                    print("  {}: {} ({:.1f}%)".format(tgt.name, tmap[tgt.name], (100.0 * float(tmap[tgt.name])) / float(total) if total else 0))
        if delta:
            print("Total time: ({:.1f}s)".format(delta))

    @property
    def counters(self):
        """Returns the underlying :class:`.counters.Counters` object, which
        contains information about site and tag counts.
        """
        return self._processor.counters

    def compute_profiles(self):
        """Computes beta/theta/c reactivity values after pair data have been processed.

           :return: a :class:`.profiles.Profiles` object, which contains the reactivities for all targets.
        """

        self._profiles = Profiles(self._targets, self.run, self._processor.counters)
        self._profiles.compute()
        return self._profiles


    def write_reactivities(self, output_path):
        """Convenience function used to write the reactivities to an output
        file. Must be called after :meth:`.compute_profiles`.

           :param output_path: the path for writing the output.
        """

        self._profiles.write(output_path)

    def store(self, output_path):
        """Saves the state of the SPATS run for later processing.

           :param output_path: the path for writing the
                 output. Recommended file extension is `.spats`
        """

        if os.path.exists(output_path):
            os.remove(output_path)
        pair_db = PairDB(output_path)
        pair_db.store_run(self.run)
        pair_db.add_targets(self.targets)
        pair_db.store_counters("spats", self.counters)

    def load(self, input_path):
        """Loads SPATS state from a file.

           :param input_path: the path of a previously saved SPATS session.
        """

        pair_db = PairDB(input_path)
        pair_db.load_run(self.run)
        self.loadTargets(pair_db)
        pair_db.load_counters("spats", self.counters)

    def validate_results(self, data_r1_path, data_r2_path, algorithm = "find_partial", verbose = False):
        """Used to validate the results of the current run using against a
           different algorithm. Must be run after running
           :meth:`.process_pair_data`, or after loading the data
           (:meth:`.load`) from a previously-run session.

        :param data_r1_path: path to R1 fragments

        :param data_r2_path: path to matching R2 fragments.

        :param algorithm: Generally the default is correct, but you
           can select a particular algorithm for data validation (see
           :attr:`.run.Run.algorithm`).

        :param verbose: set to `True` for detailed output of mismatched sites.

        :return: `True` if results validate, `False` otherwise.
        """

        original_algorithm = self.run.algorithm
        if original_algorithm == algorithm:
            raise Exception("Validation cannot be run using the same algorithm.")

        if not self.counters.registered_dict():
            raise Exception("Normal SPATS run required first in order to validate the results.")

        other = Spats()
        other.run.load_from_config(self.run.config_dict())
        other.run.algorithm = algorithm
        other._targets = self._targets
        other.process_pair_data(data_r1_path, data_r2_path)

        match_count, total = self.compare_results(other, verbose = verbose)
        if match_count == total:
            print("Original results ({} algorithm) validated using {} algorithm, {} registered sites match.".format(original_algorithm, algorithm, match_count))
            return True
        else:
            print("Validation FAILURE: results ({} algorithm) only match {}/{} registered sites (when validated using {} algorithm).".format(original_algorithm, match_count, total, algorithm))
            return False


    def compare_results(self, other_spats, verbose = False):
        """Used to compare the results of the current run against another
           SPATS instance. Must be run after running
           :meth:`.process_pair_data`, or after loading the data
           (:meth:`.load`) from a previously-run session.

        :param other_spats: :class:`.Spats` instance to compare.

        :param verbose: set to `True` for detailed output of mismatched sites.

        :return: `(match_count, total)` : `match_count` indicates the
           number of sites matched, `total` indicates total number of
           sites.
        """

        our_counts = self.counters.registered_dict()
        their_counts = other_spats.counters.registered_dict()

        match_count = 0
        total = 0
        for key, value in our_counts.iteritems():
            total += 1
            if their_counts.get(key, 0) == value:
                match_count += 1
            elif verbose:
                print("Mismatch {}:  {} != {}".format(key, value, their_counts.get(key, 0)))
        return match_count, total
