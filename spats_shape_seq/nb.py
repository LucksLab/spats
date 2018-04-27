
# utilities for working with jupyter notebooks
# meant to be leaded into the global notebook namespace, e.g.:
#  from spats_shape_seq.nb import *

import spats_shape_seq as _spats
import db as _spats_db
import matrix as _matrix
import nbutil as _nbutil
import util as _spats_util

colors = _spats_util.Colors()

def first_moment(data):
    return sum([(i + 1) * data[i] for i in range(len(data))])/sum(data)

def normalize(data):
    s = float(sum(data))
    return map(lambda x : float(x) / s, data)

def preseq_data(key):
    import json
    vals = json.loads(open('pre_{}.spats'.format(key), 'rb').read())
    res = _spats_util.SimpleObject()
    res.x_axis = range(len(vals[0]))
    res.treated = vals[0]
    res.untreated = vals[1]
    res.base = vals[2]
    res.max_val = max(map(max, vals))
    return res


class _SpatsRunData(object):

    def __init__(self, path = None):
        self._path = path
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
        path = 'run.spats'
        if self._path:
            if self._path.endswith('.spats'):
                self._path = path
            else:
                import os
                path = os.path.join(self._path, path)
        db = _spats_db.PairDB(path)
        s = _spats.Spats()
        db.load_run(s.run)
        s.run.quiet = True
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
    def single_target(self):
        return self.spats._targets.targets[0]

    @property
    def cotrans_target(self):
        return self.single_target

    @property
    def n(self):
        return self.single_target.n

    @property
    def single_profile(self):
        res = self.spats._profiles.profilesForTarget(self.single_target)
        res.x_axis = range(self.n + 1)
        return res

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
    def total_treated_muts(self):
        return self._total_counts("treated_muts")

    @property
    def total_untreated_muts(self):
        return self._total_counts("untreated_muts")

    @property
    def all_sites(self):
        return range(self.min_length, self.n + 1)

    @property
    def c_values(self):
        return [ self._profile_for_end(end).c for end in self.all_sites ]

    def row(self, end):
        res = self._profile_for_end(end)
        res.x_axis = range(end + 1)
        return res

    def column(self, site):
        res = _spats_util.SimpleObject()
        res.x_axis = []
        plots = [ "treated", "untreated", "beta", "theta", "rho" ]
        for plot_type in plots:
            setattr(res, plot_type, [])
        for end in range(self.min_length, self.n + 1):
            if site <= end:
                profiles = self._profile_for_end(end)
                res.x_axis.append(end)
                for plot_type in plots:
                    data = getattr(profiles, plot_type)
                    getattr(res, plot_type).append(data[site])
        return res

def spats_run_data(path = None):
    return _SpatsRunData(path)

def cotrans_matrix_data(data_type, max_val = 0, x_range = None, y_range = None, path = None):
    run_data = _SpatsRunData(path)
    n = run_data.n
    min_length = run_data.min_length
    rows = [ [0] * (n + 1) for x in range(0, min_length) ]
    for end in range(min_length, n + 1):
        prof = run_data._profile_for_end(end)
        row = getattr(prof, data_type)
        if max_val > 0:
            row = [ min(x, max_val) for x in row ]
        row = row + ([0] * (n + 1 - len(row)))
        rows.append(row)
    if x_range:
        rows = map(lambda r: r[x_range[0]:x_range[1]], rows)
    if y_range:
        rows = rows[y_range[0]:y_range[1]]
    return rows

def cotrans_matrix(data, max_val = 0, flags = False):
    _nbutil.create_html_cell(_matrix.matrix_html(run_data.min_length, run_data.n, run_data.spats._profiles))
