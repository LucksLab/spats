
import multiprocessing
import sys

class SpatsConfig(object):
    """A singleton is created and available as the
       :attr:`.spats.Spats.config` property on a :class:`.spats.Spats`
       instance. Change properties to control :class:`.spats.Spats`
       behavior.
    """

    def __init__(self):

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

        #: Set to ``None`` (default) to auto-set to ``k+1`` where
        #: ``k`` is the longest self-match in the target. Runs will
        #: potentially be faster if you set it higher, but may miss
        #: some pairs that have only shorter matching
        #: subsequences. You can set it to lower than ``k``, but then
        #: there's some chance pairs will match the wrong place -- in
        #: this case, they will have too many errors and be discarded
        #: -- and it will allow shorter sequences at the end (which
        #: end up being adapter-trimmed) to be accepted.
        self.minimum_target_match_length = None

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
        #: workers. Set to an integer to force an explicit number. Only
        #: used when bulk processing input data from within
        #: :class:`.spats.Spats` or :meth:`.spats.Spats.process_pair_data` or ``Spats.process_pair_db``.
        self.num_workers = None


        #: Default ``0``, increase to allow the indicated number of errors
        #: when performing adapter trimming.
        self.allowed_adapter_errors = 0


        #: Default ``0``, increase to allow the indicated number of
        #: errors when matching to the target.
        self.allowed_target_errors = 0


        # private config
        self._process_all_pairs = False  # skip uniq'ing step, force all pairs to process (sometimes useful on large pair DB)


spats_config = SpatsConfig()
