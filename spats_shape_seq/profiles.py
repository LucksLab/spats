
import math


class Profiles(object):

    def __init__(self, targets, run, counters):
        self._targets = targets
        self._counters = counters
        self._cotrans = run.cotrans
        self._run = run
        count_muts = run.count_mutations
        count_indels = run.handle_indels
        masks = run.masks
        profiles = {}
        for target in self._targets.targets:
            n = len(target.seq)
            if run.cotrans:
                for end in xrange(run.cotrans_minimum_length, n + 1):
                    profiles["{}_{}".format(target.name, end)] = TargetProfiles(self, target,
                                                                                counters.mask_counts(target, masks[0], end),
                                                                                counters.mask_counts(target, masks[1], end),
                                                                                counters.mask_muts(target, masks[0], end) if count_muts else None,
                                                                                counters.mask_muts(target, masks[1], end) if count_muts else None,
                                                                                counters.mask_edge_muts(target, masks[0], end) if count_muts else None,
                                                                                counters.mask_edge_muts(target, masks[1], end) if count_muts else None,
                                                                                counters.mask_removed_muts(target, masks[0], end) if count_muts else None,
                                                                                counters.mask_removed_muts(target, masks[1], end) if count_muts else None,
                                                                                counters.mask_inserts(target, masks[0], end) if count_indels else None,
                                                                                counters.mask_inserts(target, masks[1], end) if count_indels else None,
                                                                                counters.mask_deletes(target, masks[0], end) if count_indels else None,
                                                                                counters.mask_deletes(target, masks[1], end) if count_indels else None)
            else:
                profiles[target.name] = TargetProfiles(self, target,
                                                       counters.mask_counts(target, masks[0], n),
                                                       counters.mask_counts(target, masks[1], n),
                                                       counters.mask_muts(target, masks[0], n) if count_muts else None,
                                                       counters.mask_muts(target, masks[1], n) if count_muts else None,
                                                       counters.mask_edge_muts(target, masks[0], n) if count_muts else None,
                                                       counters.mask_edge_muts(target, masks[1], n) if count_muts else None,
                                                       counters.mask_removed_muts(target, masks[0], n) if count_muts else None,
                                                       counters.mask_removed_muts(target, masks[1], n) if count_muts else None,
                                                       counters.mask_inserts(target, masks[0], n) if count_indels else None,
                                                       counters.mask_inserts(target, masks[1], n) if count_indels else None,
                                                       counters.mask_deletes(target, masks[0], n) if count_indels else None,
                                                       counters.mask_deletes(target, masks[1], n) if count_indels else None)

        self._profiles = profiles

    def profilesForTarget(self, target):
        return self._profiles[target.name]

    def profilesForTargetNamed(self, target_name):
        return self._profiles[target_name]

    def profilesForTargetAndEnd(self, target_name, end):
        if self._cotrans:
            return self._profiles["{}_{}".format(target_name, end)]
        else:
            return self._profiles[target_name]

    def compute(self):
        for profile in self._profiles.values():
            profile.compute()

    def write(self, target_path):
        with open(target_path, 'wb') as outfile:
            outfile.write('sequence\trt_start\tfive_prime_offset\tnucleotide\ttreated_mods\tuntreated_mods\tbeta\ttheta\tc\n')
            for key in sorted(self._profiles.keys()):
                self._profiles[key].write(outfile)

    def cotrans_data(self):
        if self._cotrans:
            return [ (int(key.split('_')[-1]), prof.data()) for key, prof in self._profiles.iteritems() ]
        else:
            return [ (len(prof.data()["t"]) - 1, prof.data()) for key, prof in self._profiles.iteritems() ]

    def cotrans_keys(self):
        return sorted(self._profiles.keys(), key = lambda x : int(x.split('_')[-1]))

    def data_range(self, data_type):
        vmin = None
        vmax = None
        for profile in self._profiles.values():
            vals = getattr(profile, data_type)
            vmin = min(vals) if vmin is None else min(min(vals), vmin)
            vmax = max(vals) if vmax is None else max(max(vals), vmax)
        return (vmin, vmax)


class TargetProfiles(object):

    def __init__(self, owner, target, treated_counts, untreated_counts, treated_muts, untreated_muts, treated_edge_muts, untreated_edge_muts, treated_removed_muts, untreated_removed_muts, treated_inserts, untreated_inserts, treated_deletes, untreated_deletes):
        self.owner = owner
        self._target = target
        self.treated_counts = treated_counts
        self.untreated_counts = untreated_counts
        self.treated_muts = treated_muts
        self.untreated_muts = untreated_muts
        self.treated_edge_muts = treated_edge_muts
        self.untreated_edge_muts = untreated_edge_muts
        self.treated_removed_muts = treated_removed_muts
        self.untreated_removed_muts = untreated_removed_muts
        self.treated_inserts = treated_inserts
        self.treated_deletes = treated_deletes
        self.untreated_inserts = untreated_inserts
        self.untreated_deletes = untreated_deletes

    @property
    def treated(self):
        return self.treated_counts

    @property
    def untreated(self):
        return self.untreated_counts

    @property
    def treated_mut(self):
        return self.treated_muts

    @property
    def untreated_mut(self):
        return self.untreated_muts

    @property
    def treated_insertions(self):
        return self.treated_inserts

    @property
    def untreated_insertion(self):
        return self.untreated_inserts

    @property
    def treated_deletions(self):
        return self.treated_deletes

    @property
    def untreated_deletions(self):
        return self.untreated_deletes

    @property
    def beta(self):
        return self.betas

    @property
    def theta(self):
        return self.thetas

    @property
    def rho(self):
        return self.rhos

    @property
    def r(self):
        return self.r_mut

    def compute(self):
        treated_counts = self.treated_counts
        untreated_counts = self.untreated_counts
        n = len(treated_counts) - 1
        betas = [ 0 for x in xrange(n+1) ]
        thetas = [ 0 for x in xrange(n+1) ]
        treated_sum = 0.0    # keep a running sum
        untreated_sum = 0.0  # for both channels
        running_c_sum = 0.0  # can also do it for c
        running_c_alt_sum = 0.0

        # NOTE: there is an index discrepancy here between indices
        # used in the code, and the indices used in the Aviran paper
        # where these formulae are derived: the indices are
        # reversed. so, where in the paper the formula uses
        # \sum_{i=k}^{n+1}, in the code we use \sum_{i=0}^{k+1}, and
        # this is intentional.
        #
        # for reference, here is the comment from the original SPATS code:
        #  // TargetProfile tracks an RNA of length n. arrays have n+1 entries, 
        #  // with index 1 corresponding to the 5'-most base, and index n 
        #  // corresponding to the 3'-most base.  Index 0 stores information about 
        #  // unmodified RNAs, where RT has fallen off the 5' end of the 
        #  // strand.  This convention differs from the paper, where we refer to 
        #  // the 5'-most base as index n, and the 3'-most base as index 1.

        for k in xrange(n):
            X_k = float(treated_counts[k])
            Y_k = float(untreated_counts[k])
            treated_sum += X_k    #running sum equivalent to: treated_sum = float(sum(treated_counts[:(k + 1)]))
            untreated_sum += Y_k  #running sum equivalent to: untreated_sum = float(sum(untreated_counts[:(k + 1)]))
            try:
                Xbit = (X_k / treated_sum)
                Ybit = (Y_k / untreated_sum)
                betas[k] = (Xbit - Ybit) / (1 - Ybit)
                thetas[k] = math.log(1.0 - Ybit) - math.log(1.0 - Xbit)
                running_c_alt_sum -= math.log(1.0 - betas[k])
                if not self.owner._run.allow_negative_values:
                    betas[k] = max(0.0, betas[k])
                    thetas[k] = max(0.0, thetas[k])
                running_c_sum -= math.log(1.0 - betas[k])
            except:
                #print("domain error: {} / {} / {} / {}".format(X_k, treated_sum, Y_k, untreated_sum))
                betas[k] = 0
                thetas[k] = 0

        c = running_c_sum
        c_factor = 1.0 / c if c else 1.0
        for k in xrange(n+1):
            thetas[k] = max(c_factor * thetas[k], 0)
        self.betas = betas
        self.thetas = thetas
        self.rhos = [ n * th for th in thetas ]
        self.c = c
        self.c_alt = running_c_alt_sum

        self.compute_mutated_profiles()

    def compute_mutated_profiles(self):
        if not self.treated_muts and not self.treated_inserts and not self.treated_deletes:
            return

        treated_counts = self.treated_counts
        untreated_counts = self.untreated_counts
        treated_muts = self.treated_muts
        untreated_muts = self.untreated_muts
        treated_removed_muts = self.treated_removed_muts
        untreated_removed_muts = self.untreated_removed_muts
        treated_inserts = self.treated_inserts
        untreated_inserts = self.untreated_inserts
        treated_deletes = self.treated_deletes
        untreated_deletes = self.untreated_deletes

        n = len(treated_counts) - 1
        mu = [ 0 for x in xrange(n+1) ]
        r_mut = [ 0 for x in xrange(n+1) ]
        depth_t = 0.0    # keep a running sum
        depth_u = 0.0  # for both channels
        running_c_sum = 0.0
        c_zero = False
        running_c_alt_sum = 0.0
        c_alt_zero = False

        # NOTE: there is an index discrepancy here between indices
        # used in the code, and the indices used in the derivation:
        # the indices are reversed. so, where formula uses
        # \sum_{i=j}^{n+1}, in the code we use \sum_{i=0}^{j+1}, and
        # this is intentional.

        for j in xrange(n):

            # xref https://trello.com/c/10pysbq7/261-mutation-depth-with-quality-filtering-when-calculating-mus
            # if we removed a mut due to low quality, then we want to remove the corresponding stop
            # from the analysis (even though it should still be counted in the non-mut analysis).
            s_j_t = float(treated_counts[j] - treated_removed_muts[j])     # s_j^+
            s_j_u = float(untreated_counts[j] - untreated_removed_muts[j]) # s_j^-

            # mut_j^+ - Only one of { mut, insert, delete } is currently possible at a site per pair
            mut_j_t = 0
            if treated_muts:
                mut_j_t += float(treated_muts[j])
            if treated_inserts:
                mut_j_t += float(treated_inserts[j])
            if treated_deletes:
                mut_j_t += float(treated_deletes[j])

            # mut_j^- - Only one of { mut, insert, delete } is currently possible at a site per pair
            mut_j_u = 0
            if untreated_muts:
                mut_j_u += float(untreated_muts[j])
            if untreated_inserts:
                mut_j_u += float(untreated_inserts[j])
            if untreated_deletes:
                mut_j_u += float(untreated_deletes[j])

            depth_t += s_j_t  #running sum equivalent to: depth_t = float(sum(treated_counts[:(j + 1)]))
            depth_u += s_j_u  #running sum equivalent to: depth_u = float(sum(untreated_counts[:(j + 1)]))
            curmu = 0.0
            try:
                Tbit = (mut_j_t / depth_t)
                Ubit = (mut_j_u / depth_u)
                curmu = mu[j] = (Tbit - Ubit) / (1 - Ubit)
                if not self.owner._run.allow_negative_values:
                    mu[j] = max(0.0, mu[j])
            except:
                #print("domain error: {} / {} / {} / {}".format(s_j_t, depth_t, s_j_u, depth_u))
                mu[j] = 0.0

            if mu[j] < 1.0:
                running_c_sum -= math.log(1.0 - mu[j]) # xref Yu_Estimating_Reactivities pdf, p24
            else:
                c_zero = True
            if curmu < 1.0:
                running_c_alt_sum -= math.log(1.0 - curmu)
            else:
                c_alt_zero = True

            r_mut[j] = self.betas[j] + mu[j]

        self.mu = mu
        self.r_mut = r_mut
        if c_zero:
            self.c = 0
        else:
            self.c += running_c_sum
        if c_alt_zero:
            self.c_alt = 0
        else:
            self.c_alt += running_c_alt_sum


    def write(self, outfile):
        n = len(self.treated_counts)
        format_str = "{name}\t{rt}\t".format(name = self._target.name, rt = n - 1) + "{i}\t{nuc}\t{tm}\t{um}\t{b}\t{th}" + "\t{c:.5f}\n".format(c = self.c)
        for i in xrange(n):
            outfile.write(format_str.format(i = i,
                                            nuc = self._target.seq[i - 1] if i > 0 else '*',
                                            tm = self.treated_counts[i],
                                            um = self.untreated_counts[i],
                                            b = self.betas[i] if i > 0 else '-',
                                            th = self.thetas[i] if i > 0 else '-'))

    def data(self):
        return { "t" : self.treated_counts,
                 "u" : self.untreated_counts }
