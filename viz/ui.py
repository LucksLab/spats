
import os

import cjb.uif
import cjb.util.cfg
import viz.scenes
import viz.colorize
import viz.localizer

from spats_shape_seq import Spats
from spats_shape_seq.db import PairDB
from spats_shape_seq.pair import Pair
from spats_shape_seq.tag import TagProcessor
from spats_shape_seq.util import reverse_complement



class SpatsViz(cjb.uif.UIServer):

    def __init__(self):
        self.all_config = cjb.util.cfg.parse(os.path.expanduser("~/.spats_viz"))
        self.config = self.all_config["server"]
        cjb.uif.UIServer.__init__(self, self.config["host"], int(self.config["port"]), self.config["portfile"])
        self.colors = viz.colorize.Colorize(self.all_config.get("colors"))
        self.addFilter(self.colors)
        self.addFilter(viz.localizer.Localizer())
        self.reloadModel()

    def reloadModel(self):
        self._db = None
        self._spats = None

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
        data_config = self.all_config["data"]
        self._db = PairDB(data_config["dbfile"])
        self.result_set_id = self._db.result_set_id_for_name(data_config["result_set_name"])
        self._db.index_results()
        self.has_tags = bool(self.result_set_id)

        s = Spats()
        s.run._processor_class = TagProcessor
        s.run.allow_indeterminate = True
        s.run.allowed_target_errors = 2
        s.run.allowed_adapter_errors = 2
        if "run" in self.all_config:
            s.run.load_from_config(self.all_config["run"])

        s.loadTargets(self._db)
        p = s._processor
        for t in s._targets.targets:
            p.addTagTarget(t.name, t.seq)
            p.addTagTarget(t.name + "_rc", reverse_complement(t.seq))
            self.colors._colors[t.name.lower()] = self.colors.color("target")
        p.addTagTarget("adapter_t_rc", reverse_complement(s.run.adapter_t))
        p.addTagTarget("adapter_b", s.run.adapter_b)
        p.addTagTarget("linker_cotrans", s.run.cotrans_linker)
        p.addTagTarget("linker_cotrans_rc", reverse_complement(s.run.cotrans_linker))
        for tag, seq in self.all_config.get("extra_tags", {}).iteritems():
            p.addTagTarget(tag, seq)

        if self.has_tags:
            p.counters.load_from_db_data(self._db.counter_data_for_results(self.result_set_id))

        if "counters_key" in data_config:
            self._db.load_counters(data_config["counters_key"], p.counters)

        self._spats = s

    def home(self, message = None):
        self.reloadModel()
        sc = self.popToTop()
        if not sc:
            self.setScene(viz.scenes.Home(self))
