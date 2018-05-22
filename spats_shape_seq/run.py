
import ast
import sys

from lookup import LookupProcessor, CotransLookupProcessor
from partial import PartialFindProcessor, CotransPartialFindProcessor
from native import CotransNativeProcessor
from tag import TagProcessor


class Run(object):
    """Encapsulates the inputs/config required for a Spats run.
    """

    def __init__(self):

        #: default ``[ 'RRRY', 'YYYR' ]``, treated mask is first.
        self.masks = [ 'RRRY', 'YYYR' ]

        #: default ``AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT``
        self.adapter_t = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"  # JJB: shows revcomped on the end of R2, from the back

        #: default ``AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC``
        self.adapter_b = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"                          # JJB: shows as-is on the end of R1, from the front

        #: default ``False``, set to true to suppress output messages
        self.quiet = False

        #: default ``False``, set to ``True`` to output detailed run information
        self.debug = False

        #: defaults to ``sys.stdout``, set to file-like object to gather debugging info
        self.log = sys.stdout

        #: Defaults to `10`. Runs will
        #: potentially be faster if you set it higher, but may miss
        #: some pairs that have only shorter matching
        #: subsequences. You can set it lower, but then
        #: there's some chance pairs will match the wrong place -- in
        #: which case, they will have too many errors and be discarded
        #: -- and it will allow shorter sequences at the end (which
        #: end up being adapter-trimmed) to be accepted. Might want to
        #: analyze your targets (xref
        #: :meth:`.target.Target.longest_target_self_matches`) to
        #: determine an appropriate value.
        self.minimum_target_match_length = 10


        #: Defaults to ``False``. If set to ``True``, will count both
        #: stops and muations, and incorporate the mutation information
        #: into the reactivity profile computations. Note that setting
        #: this will force ``allowed_target_errors`` to be ``1``.
        self.count_mutations = False

        #: Defaults to ``None``. If set to a phred-score ascii integer value (33 - 75), and
        #: ``count_mutations`` is ``True``, then this will require the
        #: quality score on any mutation to be greater than or equal
        #: to the indicated phred score to be counted.
        self.mutations_require_quality_score = None

        #: Defaults to ``None``. If set to ``stop_and_mut``, will count
        #: mutations that are at the site like any other mutation. If
        #: set to ``stop_only``, will count the stop but no mutation. In
        #: the default behavior, neither stops nor edge mutations are counted.
        self.count_edge_mutations = None


        #: Defaults to ``0``, set higher to require a minimal amount of
        #: adapter in order to do trimming. Generally not necessary since
        #: a positive match in the target is required before trimming.
        self.minimum_adapter_len = 0

        #: Default ``False``, set to ``True`` to allow indeterminate
        #: nucleotides to be processed as matches. For example, with this
        #: set to ``True``, the sequence ``ACNT`` in a pair will be
        #: considered a match (no error) to ``ACGT`` in the target. See
        #: also :attr:`.allowed_adapter_errors` and
        #: :attr:`.allowed_target_errors`.
        self.allow_indeterminate = False

        #: Default ``None``, which auto-detects the number of available cores
        #: (using ``multiprocessing.cpu_count()`` and creates that many
        #: workers. Set to an integer to force an explicit number. Set to ``1``
        #: to disable multiprocessing (sometimes useful for debugging). Only
        #: used when bulk processing input data from within
        #: :class:`.spats.Spats` or :meth:`.spats.Spats.process_pair_data` or ``Spats.process_pair_db``.
        self.num_workers = None

        #: Default ``0``, increase to allow the indicated number of errors
        #: when performing adapter trimming.
        self.allowed_adapter_errors = 0

        #: Default ``0``, increase to allow the indicated number of
        #: errors when matching to the target.
        self.allowed_target_errors = 0

        #: Default ``False``, set to ``True`` to write the results back to 
        #: the input database for further analysis or incremental/resumable
        #: processing. xref :attr:`.resume_processing`.
        self.writeback_results = False

        #: Defaults to ``"default"``, set to a string to choose the name for
        #: this result set. Also used for resuming processing, xref
        #: :attr:`.resume_processing`. Result sets can be compared using
        #: :meth:`.db.PairDB.differing_results`
        self.result_set_name = None

        #: Default ``False``, set to ``True`` to resume processing (if
        #: there's been a previous run using :attr:`.writeback_results`
        #: and the same :attr:`.result_set_name`).
        self.resume_processing = False

        #: Default ``True``, set to ``False`` to parse to a database and
        #: only process unique counts. Typically straight parsing is faster,
        #: but on some datasets it can be faster to determine
        #: and only process unique counts.
        self.skip_database = True

        #: Default ``None``, in which case the pair length is detected
        #: from input data. Otherwise, can be set explicitly.
        self.pair_length = None

        #: Default ``False``, in which the right edge must match the
        #: edge of the target. Set to ``True`` to allow other
        #: possibilities for the right edge.
        self.allow_multiple_rt_starts = False

        #: Default ``False``, set to ``True`` to run as a cotrans
        #: experiment. Pass a single target instead of a generated targets
        #: file.
        self.cotrans = False

        #: Default ``CTGACTCGGGCACCAAGGAC``, change as necessary for cotrans experiments.
        self.cotrans_linker = 'CTGACTCGGGCACCAAGGAC'

        #: Default ``20``, set to adjust the minimum number of bp to use from
        #: the cotrans target.
        self.cotrans_minimum_length = 20

        #: Default ``8``, set to adjust the minimum length for matching tags for the reads analyzer.
        self.minimum_tag_match_length = 8

        #: Default ``lookup``, set to ``find_partial`` to use the partial find algorithm.
        self.algorithm = "lookup"

        #: Default ``False``, set to ``True`` to allow beta, theta,
        #: and rho values to be negative (otherwise, negative values are
        #: set to ``0.0``).
        self.allow_negative_values = False

        #: Default ``False``, set to ``True`` to count and report information
        #: on pairs that align left of the 5' end. the count for each
        #: different prefix encountered will be reported. Setting this
        #: will force the usage of the ``find_partial`` algorithm.
        self.count_left_prefixes = False

        #: Default ``False``, set to ``True`` to treat any read with a
        #: 5' prefix as starting at site zero. Setting this
        #: will force the usage of the ``count_left_prefixes`` option.
        self.collapse_left_prefixes = False

        #: Default ``None``, which means collapse all prefixes. Set to
        #: a list of comma-separated strings to only collapse prefixes
        #: that appear in the list. Only applies with
        #: ``collapse_left_prefixes`` is ``True``
        self.collapse_only_prefixes = None

        #: Default ``None``, set to a string path to generate SAM outputs
        #: for the spats run.
        self.generate_sam = None

        # private config that should be persisted (use _p_ prefix)
        self._p_use_tag_processor = False
        self._p_processor_class = None
        self._p_v102_compat = False
        self._p_extra_tags = None # used by reads analyzer

        # private config that should not be persisted (use _ prefix)
        self._run_limit = 0 # for testing, only supported on num_workers=1
        self._process_all_pairs = False  # skip uniq'ing step, force all pairs to process (sometimes useful on large pair DB)
        self._parse_quality = False
        self._applied_restrictions = False


    def apply_config_restrictions(self):
        if self._applied_restrictions:
            return
        if self._p_use_tag_processor:
            self.count_mutations = False
        if self.count_mutations:
            self.allowed_target_errors = max(self.allowed_target_errors, 1)
        if self.collapse_left_prefixes:
            self.count_left_prefixes = True
        if self.count_left_prefixes or self.allow_multiple_rt_starts or self.allowed_target_errors > 1:
            self.algorithm = 'find_partial'
        if self.collapse_only_prefixes:
            self._p_collapse_only_prefix_list = [ x.strip() for x in self.collapse_only_prefixes.split(',') ]
        if self.count_mutations and self.mutations_require_quality_score:
            self._parse_quality = True
        if self.generate_sam:
            self.algorithm = 'find_partial'
            self.num_workers = 1
            self._parse_quality = True
        self._applied_restrictions = True
        self.validate_config()

    def validate_config(self):
        if self.mutations_require_quality_score and (self.mutations_require_quality_score < 33 or self.mutations_require_quality_score > 75):
            raise Exception('Invalid mutations_require_quality_score value: {}'.format(self.mutations_require_quality_score))
        if self.count_edge_mutations and (self.count_edge_mutations != 'stop_only' and self.count_edge_mutations != 'stop_and_mut'):
            raise Exception('Invalid count_edge_mutations value: {}'.format(self.count_edge_mutations))

    def _get_processor_class(self):
        self.apply_config_restrictions()
        if self._p_use_tag_processor:
            return TagProcessor
        else:
            return self._get_base_processor_class()

    def _get_base_processor_class(self):
        implementations = {
            "lookup" : LookupProcessor,
            "find_partial" : PartialFindProcessor,
            "cotrans_lookup" : CotransLookupProcessor,
            "cotrans_find_partial" : CotransPartialFindProcessor,
            "cotrans_native" : CotransNativeProcessor,
        }
        return implementations[("cotrans_" if self.cotrans else "") + self.algorithm]

    def config_dict(self):
        config = {}
        for attr in dir(self):
            if attr.startswith('_')  and  not attr.startswith('_p_'):
                continue
            val = getattr(self, attr)
            if callable(val) or attr == "log":
                continue
            config[attr] = val
        return config

    def config_string(self):
        config = self.config_dict()
        return "\n".join([ "{} = {}".format(attr, config[attr]) for attr in config.keys() ])

    def load_from_config(self, config_dict):
        for key in config_dict.keys():
            val = config_dict[key]
            try:
                val = ast.literal_eval(val)
            except:
                pass
            #print("run set {} = {} ({})".format(key, val, val.__class__))
            setattr(self, key, val)
