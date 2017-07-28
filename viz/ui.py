
import os

import cjb.uif
import cjb.util
import cjb.util.cfg
import viz.scenes
import viz.colorize
import viz.localizer
import viz.plotter

from spats_shape_seq import Spats
from spats_shape_seq.db import PairDB
from spats_shape_seq.pair import Pair
from spats_shape_seq.tag import TagProcessor
from spats_shape_seq.util import reverse_complement


class SpatsViz(cjb.uif.UIServer):

    def __init__(self):
        self.last_path = "/tmp/viz.info"
        self.db_name = "<No file loaded>"
        self.all_config = cjb.util.cfg.parse(os.path.expanduser("~/.spats_viz"))
        self.config = self.all_config.get("server", {})
        cjb.uif.UIServer.__init__(self, self.config.get("host", "0.0.0.0"), int(self.config.get("port", 17862)), self.config.get("portfile", "/tmp/uif.port"))
        self.colors = viz.colorize.Colorize(self.all_config.get("colors"))
        self.plotter = viz.plotter.Plotter()
        self.addFilter(self.colors)
        self.addFilter(viz.localizer.Localizer())
        self.reloadModel()

    def waitFor(self):
        cjb.uif.UIServer.waitFor(self)
        self.plotter.stop()

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

        rsnames = self._db.result_sets()
        if rsnames:
            self.result_set_id = self._db.result_set_id_for_name(rsnames[0])
            self._db.index_results()
            self.has_tags = bool(self.result_set_id)
        else:
            self.result_set_id = -1
            self.has_tags = False

        self.has_counters = self._db.has_counters()

        s = Spats()
        self._db.load_run(s.run)
        s.run._p_use_tag_processor = True
        s.loadTargets(self._db)
        if self.has_counters:
            self._db.load_counters("spats", s.counters)

        if self.has_tags:
            p = s._processor
            for t in s._targets.targets:
                p.addTagTarget(t.name, t.seq)
                p.addTagTarget(t.name + "_rc", reverse_complement(t.seq))
                self.colors._colors[t.name.lower()] = self.colors.color("target")
            p.addTagTarget("adapter_t_rc", reverse_complement(s.run.adapter_t))
            p.addTagTarget("adapter_b", s.run.adapter_b)
            if s.run.cotrans:
                p.addTagTarget("linker_cotrans", s.run.cotrans_linker)
                p.addTagTarget("linker_cotrans_rc", reverse_complement(s.run.cotrans_linker))
            if s.run._p_extra_tags:
                for tag, seq in s.run._p_extra_tags.iteritems():
                    p.addTagTarget(tag, seq)
            if not self.has_counters:
                p.counters.load_from_db_data(self._db.counter_data_for_results(self.result_set_id))

        self._spats = s

    def try_load_last(self):
        try:
            last_path = cjb.util.jsonAtPath(self.last_path)["last"]
            self.open_spats(last_path)
        except:
            pass

    def open_spats(self, path):
        self._db = PairDB(path)
        self._loadDBAndModel()
        self.db_name = os.path.basename(path)
        cjb.util.writeJsonToPath({ "last" : path}, self.last_path)

    def has_reads_data(self):
        return self._db and self.has_tags

    def has_counter_data(self):
        return self._db and (self.has_counters or self.has_tags)

    def home(self, message = None):
        sc = self.popToTop()
        if not sc:
            self.setScene(viz.scenes.Home(self))
