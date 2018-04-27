
def _dict_incr(d, key, m = 1):
    d[key] = d.get(key, 0) + m

class Counters(object):

    def __init__(self, run):
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

    def _mut_edge_key(self, pair, mut):
        return "{}:{}:S{}M{}:{}".format(pair.target.rowid, pair.mask.chars, pair.site, mut, pair.end)

    def register_count(self, pair):
        if pair.mutations:
            for mut in pair.mutations:
                if pair.site == mut + 1:
                    # mutation on the edge
                    # xref https://trello.com/c/FulYfVjT/200-stop-map-mutation-on-edge-case
                    count_muts = self._run.count_edge_mutations
                    if count_muts == 'stop_and_mut':
                        _dict_incr(self._registered, self._mut_key(pair, mut), pair.multiplicity)
                        self.mutations += pair.multiplicity
                        _dict_incr(self._counts, pair.mask.chars + "_mut", pair.multiplicity)
                    elif count_muts == 'stop_only':
                        pass
                    else:
                        # don't count this as a stop at all
                        return
                else:
                    _dict_incr(self._registered, self._mut_key(pair, mut), pair.multiplicity)
                    self.mutations += pair.multiplicity
                    _dict_incr(self._counts, pair.mask.chars + "_mut", pair.multiplicity)
        _dict_incr(self._registered, self._count_key(pair), pair.multiplicity)
        self.registered_pairs += pair.multiplicity
        _dict_incr(self._counts, pair.mask.chars + "_kept", pair.multiplicity)

    def increment_mask(self, mask, multiplicity = 1):
        _dict_incr(self._counts, mask.chars + "_total", multiplicity)

    def increment_key(self, counter_key, multiplicity = 1):
        _dict_incr(self._counts, counter_key, multiplicity)

    def register_prefix(self, prefix, pair):
        self.increment_key('prefix_{}_{}'.format(pair.mask.chars, prefix), pair.multiplicity)

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
        for key, value in self._registered:
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

    def site_count(self, target_id, mask, end, site):
        return self._registered.get("{}:{}:{}:{}".format(target_id, mask, site, end), 0)

    def site_mut_count(self, target_id, mask, end, site):
        return self._registered.get("{}:{}:M{}:{}".format(target_id, mask, site, end), 0)

    def load_from_db_data(self, data):
        c = self._registered
        for r in data:
            c["{}:{}:{}:{}".format(r[0], r[1], r[2], r[3])] = r[4]

