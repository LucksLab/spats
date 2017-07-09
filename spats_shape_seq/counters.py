
def _dict_incr(d, key, m = 1):
    d[key] = d.get(key, 0) + m

class Counters(object):

    def __init__(self):
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

    def _count_key(self, pair):
        return "{}:{}:{}:{}".format(pair.target.rowid, pair.mask.chars, pair.site, pair.right)

    def register_count(self, pair):
        _dict_incr(self._registered, self._count_key(pair), pair.multiplicity)
        self.registered_pairs += pair.multiplicity
        _dict_incr(self._counts, pair.mask.chars + "_kept", pair.multiplicity)

    def increment_mask(self, mask, multiplicity = 1):
        _dict_incr(self._counts, mask.chars + "_total", multiplicity)

    def count_data(self):
        return (self._counts, self._registered)

    def update_with_count_data(self, count_data):
        their_counts, their_registered = count_data
        for key, their_value in their_counts.iteritems():
            _dict_incr(self._counts, key, their_value)
        for key, their_value in their_registered.iteritems():
            _dict_incr(self._registered, key, their_value)

    def mask_total(self, mask):
        return self._counts.get(mask.chars + "_total")

    def mask_kept(self, mask):
        return self._counts.get(mask.chars + "_kept")

    def target_total(self, target):
        total = 0
        for key, value in self._registered:
            if key.startswith("{}:".format(target.rowid)):
                total += value
        return total
