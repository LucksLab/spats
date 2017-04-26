
import multiprocessing

class SpatsConfig(object):

    def __init__(self):

        self.adapter_t = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"  # JJB: shows revcomped on the end of R2, from the back
        self.adapter_b = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"                          # JJB: shows as-is on the end of R1, from the front

        self.debug = False
        self.quiet = False
        self.minimum_target_match_length = 8 # want to make this as big as possible for speed.
        self.minimum_adapter_len = 0
        self.allow_indeterminate = False
        self.show_progress = True
        self.show_id_to_site = False
        self.num_workers = multiprocessing.cpu_count()

        # TODO: would it be better to base these off of some kind of acceptable %/freq, to control for error:length ratios?
        self.allowed_adapter_errors = 0
        self.allowed_target_errors = 0


spats_config = SpatsConfig()
