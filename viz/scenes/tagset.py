
import copy
import cjb.uif

from cjb.uif.layout import Size, Grid, Rect, layoutInScroller
from cjb.uif.views import Label
from viz.scenes.base import BaseScene
from viz.scenes.matches import Matches
from viz.layout import buttonSize


class Tag(object):

    def __init__(self, tag, percent):
        self.tag = tag
        self.percent = percent

    @property
    def displayName(self):
        return "{}{} {}".format(self.tag, ":" if self.percent else "", self.percent)

class TagAction(object):

    def __init__(self, tag, incl_excl):
        self.tag = tag
        self.incl_excl = incl_excl

    @property
    def displayName(self):
        return "+" if self.incl_excl else "-"

    def new_tags(self, incl_excl_tags):
        tag = self.tag.tag
        if self.incl_excl:
            if tag in incl_excl_tags[0]:
                raise Exception("Tag already in query: {}".format(tag))
            elif tag in incl_excl_tags[1]:
                excl = copy.copy(incl_excl_tags[1])
                excl.remove(tag)
                return [ incl_excl_tags[0], excl ]
            else:
                return [ incl_excl_tags[0] + [ tag ], incl_excl_tags[1] ]
        else:
            if tag in incl_excl_tags[1]:
                raise Exception("Tag already excluded from query: {}".format(tag))
            elif tag in incl_excl_tags[0]:
                incl = copy.copy(incl_excl_tags[0])
                incl.remove(tag)
                return [ incl, incl_excl_tags[1] ]
            return [ incl_excl_tags[0], incl_excl_tags[1] + [ tag ] ]


COLOR_BG = [0.8, 0.8, 0.8]
COLOR_INCL = [0.6, 1.0, 0.6]
COLOR_EXCL = [1.0, 0.6, 0.6]

class Tagset(BaseScene):

    def __init__(self, ui, incl_excl_tags = [ [], [] ]):
        self.incl_excl_tags = incl_excl_tags
        self.query_scroller = None
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def addTagView(self, tag, bg, skip_incl = False, skip_excl = False):
        v = cjb.uif.views.View(obj = tag)
        v.bg = bg
        v.name_label = v.addSubview(cjb.uif.views.Label(tag.displayName, fontSize = 11))
        if skip_incl:
            v.incl = None
        else:
            v.incl = v.addSubview(cjb.uif.views.Button(obj = TagAction(tag, True)))
            self.addView(v.incl)
            v.incl.bg = COLOR_INCL
        if skip_excl:
            v.excl = None
        else:
            v.excl = v.addSubview(cjb.uif.views.Button(obj = TagAction(tag, False)))
            self.addView(v.excl)
            v.excl.bg = COLOR_EXCL
        self.addView(v)
        return v

    def addQueryTag(self, tag):
        if tag in self.incl_excl_tags[0]:
            tv = self.addTagView(Tag(tag, ""), COLOR_INCL, skip_incl = True)
        else:
            tv = self.addTagView(Tag(tag, ""), COLOR_EXCL, skip_excl = True)
        return tv

    def build(self):
        BaseScene.build(self)
        counts = self.ui.db.tag_counts(self.ui.result_set_id, self.incl_excl_tags[0], self.incl_excl_tags[1])
        counts = sorted([ (key, counts[key]) for key in counts.keys() ], key = lambda x : x[1], reverse = True)
        total = float(self.ui.db.count_matches(self.ui.result_set_id, self.incl_excl_tags[0], self.incl_excl_tags[1]))
        tags = [ Tag(tc[0], "{:.1f}%".format(float(tc[1]) * 100.0 / total)) for tc in counts ]
        query_tags = self.incl_excl_tags[0] + self.incl_excl_tags[1]
        self.tagViews = [ self.addTagView(t, COLOR_BG) for t in tags if t.tag not in query_tags ]
        self.targetButtons([self.showMatches])
        self.query_label = self.addView(cjb.uif.views.Label("Query: {} matches.".format(int(total)), fontSize = 14))
        self.queryTagViews = [ self.addQueryTag(t) for t in query_tags ]

    def layout(self, view):
        BaseScene.layout(self, view)
        cur = view.frame.centeredSubrect(w = 800, h = view.frame.size.h - 100)
        left = cur.leftSubrect(w = 400)
        self.scroller = layoutInScroller(self.tagViews, cur.leftover, Size(300, 40), 20, self.scroller)
        top = left.topSubrect(40)
        left = left.leftover
        self.query_label.frame = top.leftSubrect(200).centeredSubrect(w = 200, h = 16)
        self.buttonWithKey('showMatches').frame = top.leftover.leftCenteredSubrect(size = buttonSize, margin = 20)
        self.query_scroller = layoutInScroller(self.queryTagViews, left, Size(300, 40), 20, self.query_scroller)
        for tv in self.tagViews + self.queryTagViews:
            f = tv.frame.bounds()
            tv.name_label.frame = f.leftCenteredSubrect(w = 212, h = 16, margin = 8)
            if tv.incl:
                tv.incl.frame = f.leftover.leftCenteredSubrect(w = 36, h = 20, margin = 2)
            if tv.excl:
                tv.excl.frame = f.leftover.rightCenteredSubrect(w = 36, h = 20, margin = 2)
        return view

    def handleViewMessage(self, scene, obj, message):
        if obj and isinstance(obj, TagAction):
            self.ui.pushScene(Tagset(self.ui, obj.new_tags(self.incl_excl_tags)))
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

    def showMatches(self, message = None):
        self.ui.pushScene(Matches(self.ui, include_tags = self.incl_excl_tags[0], exclude_tags = self.incl_excl_tags[1]))
