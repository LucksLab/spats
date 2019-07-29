
def _dict_incr(d, key, m = 1):
    d[key] = d.get(key, 0) + m

class Counters(object):

    def __init__(self, run = None):
        self._run = run
        self.reset()

    def reset(self):
        self._counts = {}
        self._registered = {}

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
        return "{}:{}:{}:{}".format(pair.target.rowid, pair.mask.chars, pair.site, pair.end)

    def _mut_key(self, pair, mut):
        return "{}:{}:M{}:{}".format(pair.target.rowid, pair.mask.chars, mut, pair.end)

    def _low_quality_mut_key(self, pair):
        return "{}:{}:Mq{}:{}".format(pair.target.rowid, pair.mask.chars, pair.site, pair.end)

    def _mut_edge_key(self, pair, mut):
        return "{}:{}:S{}M{}:{}".format(pair.target.rowid, pair.mask.chars, pair.site, mut, pair.end)

    def _indel_key(self, pair, spot, indel):
        if indel.insert_type:
            return "{}:{}:I{}:{}".format(pair.target.rowid, pair.mask.chars, spot + 1, pair.end)
        else:
            return "{}:{}:D{}:{}".format(pair.target.rowid, pair.mask.chars, spot + 1, pair.end)

    def register_count(self, pair):
        if pair.mutations:
            for mut in pair.mutations:
                if pair.site == mut - 1:
                    # mutation on the edge
                    # xref https://trello.com/c/FulYfVjT/200-stop-map-mutation-on-edge-case
                    count_muts = self._run.count_edge_mutations
                    self.edge_muts += 1
                    if count_muts == 'stop_and_mut':
                        pair.edge_mut = 'stop_and_mut'
                        _dict_incr(self._registered, self._mut_key(pair, mut), pair.multiplicity)
                        self.mutations += pair.multiplicity
                        _dict_incr(self._counts, pair.mask.chars + "_mut", pair.multiplicity)
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
                    _dict_incr(self._counts, pair.mask.chars + "_mut", pair.multiplicity)
        if pair.removed_mutations:
            for mut in pair.removed_mutations:
                _dict_incr(self._registered, self._low_quality_mut_key(pair), pair.multiplicity)
        if pair.r1.indels or pair.r2.indels:
            # only count one indel at a spot per pair 
            for spot in set(pair.r1.indels.keys() + pair.r2.indels.keys()):
                indel = pair.r1.indels.get(spot, pair.r2.indels.get(spot))   # assumes types and length match at spot
                spot += 1    # treat indels as 1-based (like mutations) for reactivity computations
                _dict_incr(self._registered, self._indel_key(pair, spot, indel), pair.multiplicity)
                _dict_incr(self._counts, pair.mask.chars + "_indels", pair.multiplicity)
                _dict_incr(self._counts, 'mapped_indel_len_{}'.format(len(indel.seq)), pair.multiplicity)
                self.indels += pair.multiplicity
                if spot <= pair.site  or  (not indel.insert_type and (spot - len(indel.seq)) < pair.site):
                    pair.edge_indel = True
            if pair.edge_indel:
                self.edge_indel += 1
            self.r1_indels += (pair.multiplicity * len(pair.r1.indels))
            self.r2_indels += (pair.multiplicity * len(pair.r2.indels))
        _dict_incr(self._registered, self._count_key(pair), pair.multiplicity)
        self.registered_pairs += pair.multiplicity
        _dict_incr(self._counts, pair.mask.chars + "_kept", pair.multiplicity)

    def increment_mask(self, mask, multiplicity = 1):
        _dict_incr(self._counts, mask.chars + "_total", multiplicity)

    def increment_key(self, counter_key, multiplicity = 1):
        _dict_incr(self._counts, counter_key, multiplicity)

    def register_prefix(self, prefix, pair):
        self.increment_key('prefix_{}_{}'.format(pair.mask.chars, prefix), pair.multiplicity)

    def register_mapped_prefix(self, prefix, pair):
        self.increment_key('mapped_prefix_{}_{}'.format(pair.mask.chars, prefix), pair.multiplicity)

    def register_mut_count(self, pair):
        self.increment_key('mut_count_{}'.format(len(pair.mutations) if pair.mutations else 0), pair.multiplicity)

    def register_mapped_mut_count(self, pair):
        self.increment_key('mapped_mut_count_{}'.format(len(pair.mutations) if pair.mutations else 0), pair.multiplicity)

    def count_data(self):
        return (self._counts, self._registered)

    def update_with_count_data(self, count_data):
        their_counts, their_registered = count_data
        for key, their_value in their_counts.iteritems():
            _dict_incr(self._counts, key, their_value)
        for key, their_value in their_registered.iteritems():
            _dict_incr(self._registered, key, their_value)

    def mask_total(self, mask):
        return self._counts.get(mask.chars + "_total", 0)

    def mask_kept(self, mask):
        return self._counts.get(mask.chars + "_kept", 0)

    def target_total(self, target):
        total = 0
        for key, value in self._registered.items():
            if key.startswith("{}:".format(target.rowid)):
                total += value
        return total

    def mask_counts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:{}:{}".format(target.rowid, mask, site, end), 0) for site in range(end + 1) ]

    def mask_muts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:M{}:{}".format(target.rowid, mask, site, end), 0) for site in range(end + 1) ]

    def mask_edge_muts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:S{}M{}:{}".format(target.rowid, mask, site, site + 1, end), 0) for site in range(end + 1) ]

    def mask_removed_muts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:Mq{}:{}".format(target.rowid, mask, site, end), 0) for site in range(end + 1) ]

    def mask_inserts(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:I{}:{}".format(target.rowid, mask, site, end), 0) for site in range(end + 1) ]

    def mask_deletes(self, target, mask, end):
        c = self._registered
        return [ c.get("{}:{}:D{}:{}".format(target.rowid, mask, site, end), 0) for site in range(end + 1) ]

    def site_count(self, target_id, mask, end, site):
        return self._registered.get("{}:{}:{}:{}".format(target_id, mask, site, end), 0)

    def site_mut_count(self, target_id, mask, end, site):
        return self._registered.get("{}:{}:M{}:{}".format(target_id, mask, site, end), 0)

    def load_from_db_data(self, data):
        c = self._registered
        for r in data:
            c["{}:{}:{}:{}".format(r[0], r[1], r[2], r[3])] = r[4]

