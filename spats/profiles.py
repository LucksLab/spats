
import math


class Profiles(object):

    def __init__(self, targets, masks):
        self._targets = targets
        self._masks = masks
        self._target = targets[0]

    def compute(self):
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
        self.treated_counts = treated_counts
        self.untreated_counts = untreated_counts
        self.betas = betas
        self.thetas = thetas
        self.c = c

    def write(self, target_path):
        n = self._target.n
        with open(target_path, 'wb') as outfile:
            outfile.write('sequence\trt_start\tfive_prime_offset\tnucleotide\ttreated_mods\tuntreated_mods\tbeta\ttheta\tc\n')
            format_str = "{name}\t{rt}\t".format(name = self._target.name, rt = n - 1) + "{i}\t{nuc}\t{tm}\t{um}\t{b}\t{th}" + "\t{c:.5f}\n".format(c = self.c)
            # TODO: xref https://trello.com/c/OtbxyiYt/23-3-nt-missing-from-reactivities-out
            # looks like we may want this to be range(n), chopping off was unintentional bug of previous version
            for i in range(n - 1):
                outfile.write(format_str.format(i = i,
                                                nuc = self._target.seq[i - 1] if i > 0 else '*',
                                                tm = self.treated_counts[i],
                                                um = self.untreated_counts[i],
                                                b = self.betas[i] if i > 0 else '-',
                                                th = self.thetas[i] if i > 0 else '-'))
