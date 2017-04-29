
import math


class Profiles(object):

    def __init__(self, targets, masks):
        self._targets = targets
        profiles = {}
        for target in self._targets.targets:
            # assume first mask is treated
            if len(masks) != 2:
                raise Exception("Expect exactly two masks to distinguish treated/untreated.")
            profiles[target.name] = TargetProfiles(target, masks[0].counts(target), masks[1].counts(target))
        self._profiles = profiles

    def profilesForTarget(self, target):
        return self._profiles[target.name]

    def profilesForTargetNamed(self, target_name):
        return self._profiles[target_name]

    def compute(self):
        for profile in self._profiles.values():
            profile.compute()

    def write(self, target_path):
        with open(target_path, 'wb') as outfile:
            outfile.write('sequence\trt_start\tfive_prime_offset\tnucleotide\ttreated_mods\tuntreated_mods\tbeta\ttheta\tc\n')
            for target in self._targets.targets:
                self.profilesForTarget(target).write(outfile)


class TargetProfiles(object):

    def __init__(self, target, treated_counts, untreated_counts):
        self._target = target
        self.treated_counts = treated_counts
        self.untreated_counts = untreated_counts

    def compute(self):
        n = self._target.n
        treated_counts = self.treated_counts
        untreated_counts = self.untreated_counts
        betas = [ 0 for x in range(n+1) ]
        thetas = [ 0 for x in range(n+1) ]
        treated_sum = 0.0    # keep a running sum
        untreated_sum = 0.0  # for both channels
        running_c_sum = 0.0  # can also do it for c

        for k in range(n):
            X_k = float(treated_counts[k])
            Y_k = float(untreated_counts[k])
            treated_sum += X_k    #running sum equivalent to: treated_sum = float(sum(treated_counts[:(k + 1)]))
            untreated_sum += Y_k  #running sum equivalent to: untreated_sum = float(sum(untreated_counts[:(k + 1)]))
            try:
                Xbit = (X_k / treated_sum)
                Ybit = (Y_k / untreated_sum)
                betas[k] = max(0, (Xbit - Ybit) / (1 - Ybit))
                thetas[k] = math.log(1.0 - Ybit) - math.log(1.0 - Xbit)
                running_c_sum -= math.log(1.0 - betas[k])
            except:
                #print "domain error: {} / {} / {} / {}".format(X_k, treated_sum, Y_k, untreated_sum)
                betas[k] = 0
                thetas[k] = 0

        c = running_c_sum
        c_factor = 1.0 / c if c else 1.0
        for k in range(n+1):
            thetas[k] = max(c_factor * thetas[k], 0)
        self.betas = betas
        self.thetas = thetas
        self.c = c

    def write(self, outfile):
        n = self._target.n
        format_str = "{name}\t{rt}\t".format(name = self._target.name, rt = n - 1) + "{i}\t{nuc}\t{tm}\t{um}\t{b}\t{th}" + "\t{c:.5f}\n".format(c = self.c)
        for i in range(n):
            outfile.write(format_str.format(i = i,
                                            nuc = self._target.seq[i - 1] if i > 0 else '*',
                                            tm = self.treated_counts[i],
                                            um = self.untreated_counts[i],
                                            b = self.betas[i] if i > 0 else '-',
                                            th = self.thetas[i] if i > 0 else '-'))
