
import math
import cjb.uif

from cjb.uif.layout import Size, Grid, Rect, layoutInScroller
from cjb.uif.views import Label
from viz.scenes.base import BaseScene
from viz.scenes.matches import Matches
from viz.layout import buttonSize

class Targets(BaseScene):

    def addTargetView(self, target, bg):
        v = cjb.uif.views.Button(obj = target, text = target[0])
        v.fontSize = 11
        v.bg = bg
        self.addView(v)
        return v

    def build(self):
        BaseScene.build(self)
        targets = sorted(self.ui.db.targets(), key = lambda x : x[0])
        self.targetViews = [ self.addTargetView(t, [ 1.0, 0.85, 0.7 ]) for t in targets ]

    def layout(self, view):
        BaseScene.layout(self, view)
        cur = view.frame.centeredSubrect(w = 300, h = view.frame.size.h - 100)
        self.scroller = layoutInScroller(self.targetViews, cur, Size(200, 40), 20, self.scroller)
        return view

    def handleViewMessage(self, scene, obj, message):
        if obj and 3 == len(obj):
            self.ui.pushScene(Target(self.ui, obj))
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

class Site(object):

    def __init__(self, target_id, index, nuc, treated_count, untreated_count):
        self.target_id = target_id
        self.index = index
        self.nuc = nuc
        self.treated_count = treated_count
        self.untreated_count = untreated_count

    @property
    def total(self):
        return self.treated_count + self.untreated_count


class Target(BaseScene):

    def __init__(self, ui, target):
        self.name = target[0]
        self.seq = target[1]
        self.target_id = target[2]
        self.scroller = None
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def addSiteView(self, site):
        v = cjb.uif.views.View(obj = site)
        v.site_label = v.addSubview(cjb.uif.views.Label(str(site.index), fontSize = 11))
        v.site_label.alignment = "right"
        v.nuc_label = v.addSubview(cjb.uif.views.Label(site.nuc, fontSize = 11, bg = [ 1.0, 0.85, 0.7 ]))
        v.bar = v.addSubview(cjb.uif.views.View())
        v.bar.bg = [ 0.8, 0.6, 1.0 ]
        v.treated_label = v.addSubview(cjb.uif.views.Label(str(site.treated_count), fontSize = 11))
        v.untreated_label = v.addSubview(cjb.uif.views.Label(str(site.untreated_count), fontSize = 11))
        v.target = lambda msg : self.handleViewMessage(None, site, msg)
        v.click = 1
        self.addView(v)
        return v

    def build(self):
        BaseScene.build(self)
        sitemap = { "{}_{}".format(s[0], s[1]) : s[2] for s in self.ui.db.result_sites(self.ui.result_set_id, self.target_id) }
        n = len(self.seq)
        total = 0
        treated = [0] * (n+1)
        untreated = [0] * (n+1)
        masks = self.ui.spats.run.masks
        self.siteViews = []
        for s in range(n + 1):
            site = Site(self.target_id,
                        s,
                        self.seq[s - 1] if s else "*",
                        sitemap.get("{}_{}".format(masks[0], s), 0),
                        sitemap.get("{}_{}".format(masks[1], s), 0))
            v = self.addSiteView(site)
            self.siteViews.append(v)
            total += site.untreated_count
            total += site.treated_count
        self.total = total

    def layoutSite(self, view):
        grid = Grid(frame = view.frame.bounds(), itemSize = Size(10, 16), columns = 60, rows = 1)
        view.site_label.frame = grid.frame(0, 3)
        view.nuc_label.frame = grid.frame(5)
        f = grid.frame(6, 20)
        # sqrt(sqrt(x)) for rescaling
        factor = math.sqrt(math.sqrt((float(view.obj.total) / float(self.total))))
        f.update(origin = f.origin, w = int(factor * float(f.size.width)), h = f.size.height)
        view.bar.frame = f
        view.treated_label.frame = grid.frame(40, 8)
        view.untreated_label.frame = grid.frame(50, 8)

    def layout(self, view):
        BaseScene.layout(self, view)
        cur = view.frame.centeredSubrect(w = 600, h = view.frame.size.h - 100)
        self.scroller = layoutInScroller(self.siteViews, cur, Size(600, 14), 2, self.scroller)
        for v in self.siteViews:
            self.layoutSite(v)
        return view

    def handleViewMessage(self, scene, obj, message):
        if obj and isinstance(obj, Site):
            self.ui.pushScene(Matches(self.ui, None, site = obj))
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

