
import cjb.uif

from cjb.uif.layout import Size, Grid, Rect, layoutInScroller
from cjb.uif.views import Label
from viz.scenes.base import BaseScene
from viz.scenes.matches import Matches
from viz.layout import buttonSize


class Tag(object):

    def __init__(self, tag, percent, query = None):
        self.tag = tag
        self.percent = percent
        self.query = query

    @property
    def displayName(self):
        return "{}: {}".format(self.query or self.tag, self.percent)


class Tagset(BaseScene):

    def __init__(self, ui, pair_db, result_set_id = 1):
        self.pair_db = pair_db
        self.result_set_id = result_set_id
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def addTagView(self, tag, bg):
        v = cjb.uif.views.Button(obj = tag)
        v.fontSize = 11
        v.bg = bg
        self.addView(v)
        return v

    def build(self):
        BaseScene.build(self)
        counts = self.pair_db.tag_counts(self.result_set_id)
        counts = sorted([ (key, counts[key]) for key in counts.keys() ], key = lambda x : x[1], reverse = True)
        total = float(self.pair_db.count())
        tags = [ Tag(tc[0], "{:.1f}%".format(float(tc[1]) * 100.0 / total)) for tc in counts ]
        self.tagViews = [ self.addTagView(t, [ 0.7, 0.7, 0.9 ]) for t in tags ]

    def layout(self, view):
        BaseScene.layout(self, view)
        cur = view.frame.centeredSubrect(w = 300, h = view.frame.size.h - 100)
        self.scroller = layoutInScroller(self.tagViews, cur, Size(200, 40), 20, self.scroller)
        return view

    def handleViewMessage(self, scene, obj, message):
        if obj and isinstance(obj, Tag):
            self.ui.setScene(Matches(self, obj.tag))
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

