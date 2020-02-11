
def _dict_incr(d, key, m = 1):
    d[key] = d.get(key, 0) + m

class Counters(object):

    def __init__(self, run = None):
        self._run = run
        self.reset()

    def reset(self):
        self._counts = {}
        self._registered = {}
        self._depths = {}
        self._quality_depths = {}

    def __getattr__(self, key):
        if key.startswith('_'):
            return super.__getattr__(self, key)
        else:
            return self._counts.get(key, 0)

    def __setattr__(self, key, value):
        if key.startswith('_'):
            super.__setattr__(self, key, value)
        else:
            self._counts[key] = value

    def counts_dict(self):
        return { key : value for key, value in self._counts.iteritems() if not key.startswith('_') }

    def registered_dict(self):
        return self._registered

    def _count_key(self, pair):
        return "{}:{}:{}:{}".format(pair.target.rowid, pair.mask_label, pair.site, pair.end)

    def _depth_key(self, pair):
        return "{}:{}:{}".format(pair.target.rowid, pair.mask_label, pair.end)

    def _mut_key(self, pair, mut):
        return "{}:{}:M{}:{}".format(pair.target.rowid, pair.mask_label, mut, pair.end)

    def _low_quality_mut_key(self, pair):
        return "{}:{}:Mq{}:{}".format(pair.target.rowid, pair.mask_label, pair.site, pair.end)

    def _mut_edge_key(self, pair, mut):
        return "{}:{}:S{}M{}:{}".format(pair.target.rowid, pair.mask_label, pair.site, mut, pair.end)

    def _indel_key(self, pair, spot, indel):
        if indel.insert_type:
            return "{}:{}:I{}:{}".format(pair.target.rowid, pair.mask_label, spot + 1, pair.end)
        else:
            return "{}:{}:D{}:{}".format(pair.target.rowid, pair.mask_label, spot + 1, pair.end)

    def register_count(self, pair):
        if pair.mutations:
            for mut in pair.mutations:
                # TODO:  count mutations and indels at end as well as start
                # TODO   also count edge_muts in roi
                if pair.site == mut - 1:
                    # mutation on the edge
                    # xref https://trello.com/c/FulYfVjT/200-stop-map-mutation-on-edge-case
                    count_muts = self._run.count_edge_mutations
                    self.edge_muts += 1
                    if count_muts == 'stop_and_mut':
                        pair.edge_mut = 'stop_and_mut'
                        _dict_incr(self._registered, self._mut_key(pair, mut), pair.multiplicity)
                        self.mutations += pair.multiplicity
                        _dict_incr(self._counts, pair.mask_label + "_mut", pair.multiplicity)
                    elif count_muts == 'stop_only':
                        pair.edge_mut = 'stop_only'
                        pass
                    else:
                        # don't count this as a stop at all
                        pair.edge_mut = 'ignore'
                        return
                else:
                    _dict_incr(self._registered, self._mut_key(pair, mut), pair.multiplicity)
                    self.mutations += pair.multiplicity
                    _dict_incr(self._counts, pair.mask_label + "_mut", pair.multiplicity)
        if pair.removed_mutations:
            for mut in pair.removed_mutations:
                _dict_incr(self._registered, self._low_quality_mut_key(pair), pair.multiplicity)
        if pair.r1.indels or pair.r2.indels:
            # only count one indel at a spot per pair 
            for spot in set(pair.r1.indels.keys() + pair.r2.indels.keys()):
                indel = pair.r1.indels.get(spot, pair.r2.indels.get(spot))   # assumes types and length match at spot
                spot += 1    # treat indels as 1-based (like mutations) for reactivity computations
                _dict_incr(self._registered, self._indel_key(pair, spot, indel), pair.multiplicity)
                _dict_incr(self._counts, pair.mask_label + "_indels", pair.multiplicity)
                _dict_incr(self._counts, 'mapped_indel_len_{}'.format(len(indel.seq)), pair.multiplicity)
                self.indels += pair.multiplicity
                if spot <= pair.site  or  (not indel.insert_type and (spot - len(indel.seq)) < pair.site):
                    pair.edge_indel = True
                if indel.ambiguous:
                    pair.ambig_indel = True
                    # TODO:  check if it's next to a mut and count it
                    # TODO:  (do here or when finding amb?)
            if pair.edge_indel:
                self.edge_indel += pair.multiplicity
            if pair.ambig_indel:
                self.ambiguous_indels_pair += pair.multiplicity
            self.r1_indels += (pair.multiplicity * len(pair.r1.indels))
            self.r2_indels += (pair.multiplicity * len(pair.r2.indels))
        _dict_incr(self._registered, self._count_key(pair), pair.multiplicity)
        self._add_to_depth(pair)
        self.registered_pairs += pair.multiplicity
        _dict_incr(self._counts, pair.mask_label + "_kept", pair.multiplicity)
        # TODO: find and count complex indel where subst and repl lens may not match (will always be ambiguous?)
        # TODO:    will need to find like SHAPmapper2
        # TODO:    keep track of both lens

    def increment_mask(self, mask_label, multiplicity = 1):
        _dict_incr(self._counts, mask_label + "_total", multiplicity)

    def increment_key(self, counter_key, multiplicity = 1):
        _dict_incr(self._counts, counter_key, multiplicity)

    def _add_to_depth(self, pair):
        dk = self._depth_key(pair)
        n = min(pair.target.n, pair.end) + 1
        for spot in xrange(pair.site, n):
            self._depths.setdefault(dk, [0] * n)[spot] += pair.multiplicity
            if not pair.removed_mutations:
                self._quality_depths.setdefault(dk, [0] * n)[spot] += pair.multiplicity

    def register_prefix(self, prefix, pair):
        self.increment_key('prefix_{}_{}'.format(pair.mask_label, prefix), pair.multiplicity)

    def register_mapped_prefix(self, prefix, pair):
        self.increment_key('mapped_prefix_{}_{}'.format(pair.mask_label, prefix), pair.multiplicity)

    def register_mut_count(self, pair):
        self.increment_key('mut_count_{}'.format(len(pair.mutations) if pair.mutations else 0), pair.multiplicity)

    def register_mapped_mut_count(self, pair):
        self.increment_key('mapped_mut_count_{}'.format(len(pair.mutations) if pair.mutations else 0), pair.multiplicity)

    def count_data(self):
        return (self._counts, self._registered), (self._depths, self._quality_depths)

    def update_with_count_data(self, count_data, vect_data):
        their_counts, their_registered = count_data
        for key, their_value in their_counts.iteritems():
            _dict_incr(self._counts, key, their_value)
        for key, their_value in their_registered.iteritems():
            _dict_incr(self._registered, key, their_value)
        their_depths, their_quality_depths = vect_data
        for key, their_values in their_depths.iteritems():
            for i, their_value in enumerate(their_values):
                self._depths.setdefault(key, [0] * len(their_values))[i] += their_value
        for key, their_values in their_quality_depths.iteritems():
            for i, their_value in enumerate(their_values):
                self._quality_depths.setdefault(key, [0] * len(their_values))[i] += their_value

    def mask_total(self, mask):
        if mask.empty_place_holder:
            return self._counts.get(mask.empty_place_holder + "_total", 0)
        else:
            return self._counts.get(mask.chars + "_total", 0)

    def mask_kept(self, mask):
        if mask.empty_place_holder:
            return self._counts.get(mask.empty_place_holder + "_kept", 0)
        else:
            return self._counts.get(mask.chars + "_kept", 0)

    def target_total(self, target):
        total = 0
        for key, value in self._registered.items():
            if key.startswith("{}:".format(target.rowid)):
                total += value
        return total

    def mask_counts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:{}:{}".format(target.rowid, mask, site, end), 0) for site in xrange(end + 1) ]

    def mask_depths(self, target, mask, end):
        return self._depths.get("{}:{}:{}".format(target.rowid, mask, end), [0] * (end + 1))

    def mask_quality_depths(self, target, mask, end):
        return self._quality_depths.get("{}:{}:{}".format(target.rowid, mask, end), [0] * (end + 1))

    def mask_muts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:M{}:{}".format(target.rowid, mask, site, end), 0) for site in xrange(end + 1) ]

    def mask_edge_muts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:S{}M{}:{}".format(target.rowid, mask, site, site + 1, end), 0) for site in xrange(end + 1) ]

    def mask_removed_muts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:Mq{}:{}".format(target.rowid, mask, site, end), 0) for site in xrange(end + 1) ]

    def mask_inserts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:I{}:{}".format(target.rowid, mask, site, end), 0) for site in xrange(end + 1) ]

    def mask_deletes(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:D{}:{}".format(target.rowid, mask, site, end), 0) for site in xrange(end + 1) ]

    def site_count(self, target_id, mask, end, site):
        return self._registered.get("{}:{}:{}:{}".format(target_id, mask, site, end), 0)

    def site_mut_count(self, target_id, mask, end, site):
        return self._registered.get("{}:{}:M{}:{}".format(target_id, mask, site, end), 0)

    def load_from_db_data(self, data):
        c = self._registered
        for r in data:
            c["{}:{}:{}:{}".format(r[0], r[1], r[2], r[3])] = r[4]
            if r[2].isdigit():    # only add to depth if it's the count_key
                n = r[3]
                dk = "{}:{}:{}".format(r[0], r[1], r[3])
                for spot in xrange(r[2], r[3] + 1):
                    curdepths = self._depths.setdefault(dk, [0] * n)
                    if n > len(curdepths):
                        curdepths += [0] * (n - len(curdepths))
                    curdepths[spot] += r[4]
                # WARNING: The quality depths here are potentially
                # wrong, since we don't track low quality muts/pairs
                # in the Results or Pairs db tables.
                # Rather than leave uninitialized we consider all
                # pairs to have been of high enough quality (or that
                # the run was done with the run options of
                # `mutations_require_quality_score` set to None).
                # TAI:  If the viz tool ever gets more play, we should
                # probably track low_quality_muts in the Pairs db table.
                self._quality_depths[dk] = list(self._depths[dk])

