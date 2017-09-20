
# utilities for working with jupyter notebooks
# meant to be leaded into the global notebook namespace, e.g.:
#  from spats_shape_seq.nb import *

import nbutil as _nbutil
import spats_shape_seq as _spats
import spats_shape_seq.db as _spats_db
import spats_shape_seq.util as _spats_util

colors = _spats_util.Colors()

def first_moment(data):
    return sum([(i + 1) * data[i] for i in range(len(data))])/sum(data)

def preseq_data():
    import json
    vals = json.loads(open('pre.spats', 'rb').read())
    res = _spats_util.SimpleObject()
    res.x_axis = range(len(vals[0]))
    res.treated = vals[0]
    res.untreated = vals[1]
    res.base = vals[2]
    return res


class _SpatsRunData(object):

    def __init__(self):
        self._spats = None
        self._db = None

    @property
    def processor(self):
        if not self._spats:
            self._loadDBAndModel()
        return self._spats._processor

    @property
    def db(self):
        if not self._db:
            self._loadDBAndModel()
        return self._db

    @property
    def spats(self):
        if not self._spats:
            self._loadDBAndModel()
        return self._spats

    def _loadDBAndModel(self):
        db = _spats_db.PairDB('run.spats')
        s = _spats.Spats()
        db.load_run(s.run)
        s.loadTargets(db)
        if db.has_counters():
            db.load_counters("spats", s.counters)
        s.compute_profiles()
        self._spats = s
        self._db = db

    @property
    def cotrans(self):
        return self.spats.run.cotrans
    @property
    def cotrans_target(self):
        return self.spats._targets.targets[0]

    @property
    def n(self):
        return self.cotrans_target.n

    @property
    def min_length(self):
        return self.spats.run.cotrans_minimum_length

    def _profile_for_end(self, end):
        return self.spats._profiles.profilesForTargetAndEnd(self.cotrans_target.name, end)

    def _total_counts(self, count_type):
        counts = []
        count_sum = 0.0
        for end in self.all_sites:
            x = sum(getattr(self._profile_for_end(end), count_type))
            counts.append(x)
            count_sum += x
        return counts

    @property
    def total_treated_counts(self):
        return self._total_counts("treated")

    @property
    def total_untreated_counts(self):
        return self._total_counts("untreated")

    @property
    def all_sites(self):
        return range(self.min_length, self.n + 1)

    @property
    def c_values(self):
        return [ self._profile_for_end(end).c for end in self.all_sites ]

def spats_run_data():
    return _SpatsRunData()
