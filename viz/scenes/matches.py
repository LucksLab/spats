
import cjb.uif

from cjb.uif.layout import Size, Grid, Rect, layoutInScroller
from spats_shape_seq import Spats
from spats_shape_seq.pair import Pair
from spats_shape_seq.tag import TagProcessor
from spats_shape_seq.util import reverse_complement
from viz.scenes.base import BaseScene
from viz.scenes.pair import PairScene, RawPairScene
from viz.layout import buttonSize


class MatchedPair(object):

    def __init__(self, r1, r2, multiplicity, rowid, identifier):
        self.r1 = r1
        self.r2 = r2
        self.multiplicity = multiplicity
        self.rowid = rowid
        self.identifier = identifier

    @property
    def displayName(self):
        return self.r1 + "   " + self.r2


class Matches(BaseScene):

    def __init__(self, tagset_scene, tag):
        self.tagset_scene = tagset_scene
        self.include_tags = [tag]
        self.exclude_tags = []
        self._spats = None
        BaseScene.__init__(self, tagset_scene.ui, self.__class__.__name__)

    def addMatchView(self, pair):
        v = cjb.uif.views.Button(obj = pair)
        v.fontSize = 11
        self.addView(v)
        return v

    def build(self):
        BaseScene.build(self)
        pairs = self.tagset_scene.pair_db.results_matching(self.tagset_scene.result_set_id, self.include_tags, self.exclude_tags)
        matches = [ MatchedPair(p[2], p[3], p[4], p[0], p[1]) for p in pairs ]
        self.targetButtons([self.back])
        self.matchViews = [ self.addMatchView(m) for m in matches ]

    def layout(self, view):
        BaseScene.layout(self, view)
        self.buttonWithKey('back').frame = view.frame.topRightSubrect(size = buttonSize, margin = 20)
        cur = view.frame.centeredSubrect(w = 800, h = view.frame.size.h - 100)
        self.scroller = layoutInScroller(self.matchViews, cur, Size(600, 14), 2, self.scroller)
        return view

    @property
    def processor(self):
        if not self._spats:
            s = Spats()
            s.run._processor_class = TagProcessor
            s.run.allow_indeterminate = True
            s.run.allowed_target_errors = 2
            s.run.allowed_adapter_errors = 2
            s.loadTargets(self.tagset_scene.pair_db)
            p = s._processor
            for t in s._targets.targets:
                p.addTagTarget(t.name, t.seq)
                p.addTagTarget(t.name + "_rc", reverse_complement(t.seq))
            p.addTagTarget("adapter_t_rc", reverse_complement(s.run.adapter_t))
            p.addTagTarget("adapter_b", s.run.adapter_b)
            self._spats = s
        return self._spats._processor

    def handleViewMessage(self, scene, obj, message):
        if obj and isinstance(obj, MatchedPair):
            pair = Pair()
            pair.set_from_data(obj.identifier, obj.r1, obj.r2, obj.multiplicity)
            self.processor.process_pair_detail(pair)
            if pair.has_site:
                self.ui.setScene(PairScene(self, pair, expanded = True))
            else:
                self.ui.setScene(RawPairScene(self, pair, expanded = True))
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

    def back(self, message = None):
        self.ui.setScene(self.tagset_scene)

    def handleKeyEvent(self, keyInfo):
        handler = { "b" : self.back }.get(keyInfo["t"])
        if handler:
            handler()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

