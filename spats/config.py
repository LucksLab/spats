
import multiprocessing

class SpatsConfig(object):

    def __init__(self):

        self.adapter_t = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"  # JJB: shows revcomped on the end of R2, from the back
        self.adapter_b = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"                          # JJB: shows as-is on the end of R1, from the front

        self.debug = False
        self.quiet = False

        # want to make this as large as possible for speed.
        # set to None to auto-set to (k+1) where k=longest self-match in the Target
        # runs will go faster if you set it higher, but may miss some pairs that
        # have only shorter matching subsequences
        # you can set it to lower than k, but then there's some chance pairs will
        # match the wrong place -- in this case, they will have too many errors
        # and be discarded, which is not bad; and it will allow shorter sequences
        # at the end (which end up being adapter-trimmed) to be accepted.
        self.minimum_target_match_length = None

        self.minimum_adapter_len = 0
        self.allow_indeterminate = False
        self.show_progress = True
        self.show_id_to_site = False
        self.num_workers = multiprocessing.cpu_count()

        # TODO: would it be better to base these off of some kind of acceptable %/freq, to control for error:length ratios?
        self.allowed_adapter_errors = 0
        self.allowed_target_errors = 0


spats_config = SpatsConfig()
