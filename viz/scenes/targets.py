
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
            if self.ui.spats.run.cotrans:
                self.ui.pushScene(CotransTarget(self.ui, obj))
            else:
                self.ui.pushScene(Target(self.ui, obj))
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

    def handleKeyEvent(self, keyInfo):
        handler = { "0" : self.firstTarget }.get(keyInfo["t"])
        if handler:
            handler()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

    def firstTarget(self, message = None):
        self.handleViewMessage(None, self.targetViews[0].obj, None)


class Site(object):

    def __init__(self, target_id, site, end, nuc, treated_count, untreated_count):
        self.target_id = target_id
        self.site = site
        self.end = end
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
        v.site_label = v.addSubview(cjb.uif.views.Label(str(site.site), fontSize = 11))
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
        sitemap = { "{}_{}".format(s[0], s[2]) : s[3] for s in self.ui.db.result_sites(self.ui.result_set_id, self.target_id) }
        n = len(self.seq)
        total = 0
        treated = [0] * (n+1)
        untreated = [0] * (n+1)
        masks = self.ui.spats.run.masks
        self.siteViews = []
        for s in range(n + 1):
            site = Site(self.target_id,
                        s,
                        n,
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
        

class CotransTarget(BaseScene):

    def __init__(self, ui, target, data_type = "treated_counts"):
        self.name = target[0]
        self.seq = target[1]
        self.target_id = target[2]
        self.scroller = None
        self.data_type = data_type
        self.plot_type = "row"
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def build(self):
        BaseScene.build(self)
        self.targetButtons([self.treated, self.untreated, self.beta, self.theta, self.reactivity, self.togglePlotType])
        self.buttonWithKey('togglePlotType').text = "Plot Type: Row"
        self.type_label = self.addView(cjb.uif.views.Label("Query: {}".format(self.data_type), fontSize = 16))
        if self.ui.has_tags:
            sitemap = { "{}:{}:{}:{}".format(self.target_id, s[0], s[2], s[1]) : s[3] for s in self.ui.db.result_sites(self.ui.result_set_id, self.target_id) }
        else:
            sitemap = self.ui.spats.counters.registered_dict()
        max_count = float(max(sitemap.values()))
        seq = self.seq
        n = len(seq)
        total = 0
        spats = self.ui.spats
        self.profiles = spats.compute_profiles()
        masks = spats.run.masks
        self.siteViews = []
        self.scroller = cjb.uif.views.Scroller()
        self.addView(self.scroller)
        max_val = 0
        for end in range(spats.run.cotrans_minimum_length, n + 1):
            end_profiles = self.profiles.profilesForTargetAndEnd(self.name, end)
            data = getattr(end_profiles, self.data_type)
            for s in range(1, end + 1):
                site = Site(self.target_id,
                            s,
                            end,
                            seq[s - 1] if s else "*",
                            sitemap.get("{}:{}:{}:{}".format(self.target_id, masks[0], s, end), 0),
                            sitemap.get("{}:{}:{}:{}".format(self.target_id, masks[1], s, end), 0))
                v = cjb.uif.views.Button(obj = site)
                self.scroller.addSubview(v)
                self.addView(v)
                v.factor = float(data[s])
                max_val = max(max_val, v.factor)
                self.siteViews.append(v)
                total += site.total
        if max_val > 0:
            if self.data_type == "reactivities":
                max_val = 4.0 # TEMP
            for v in self.siteViews:
                bg = float("{:.2f}".format(min(v.factor, max_val) / max_val))
                v.bg = [ bg, bg, bg ]
        else:
            print "warning: data all zero"
        self.total = total

    def layout(self, view):
        BaseScene.layout(self, view)
        cur = view.frame.centeredSubrect(w = 1000, h = view.frame.size.h - 100)
        self.scroller.frame = cur
        n = len(self.seq)
        grid = Grid(frame = cur.bounds(), itemSize = Size(5, 5), columns = 200, rows = n - self.ui.spats.run.cotrans_minimum_length)
        spacer = 40
        for v in self.siteViews:
            site = v.obj
            v.frame = grid.frame(200 * (site.end - self.ui.spats.run.cotrans_minimum_length) + spacer + site.site)
        grid = Grid(frame = view.frame.leftCenteredSubrect(w = 120, h = 400, margin = 40), itemSize = Size(120, 40), columns = 1, rows = 5, spacing = Size(0, 40))
        grid.applyToViews(map(self.buttonWithKey, [ 'treated', 'untreated', 'beta', 'theta', 'reactivity' ]))
        self.type_label.frame = view.frame.topLeftSubrect(w = 240, h = 24, margins = Size(40, 100))
        self.buttonWithKey('togglePlotType').frame = view.frame.bottomCenteredSubrect(200, 40, margin = 24)
        return view

    def handleViewMessage(self, scene, obj, message):
        if obj and isinstance(obj, Site):
            self.show_plot(obj)
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

    def beta(self, message = None):
        self.ui.setScene(CotransTarget(self.ui, (self.name, self.seq, self.target_id), data_type = "betas"))

    def theta(self, message = None):
        self.ui.setScene(CotransTarget(self.ui, (self.name, self.seq, self.target_id), data_type = "thetas"))

    def reactivity(self, message = None):
        self.ui.setScene(CotransTarget(self.ui, (self.name, self.seq, self.target_id), data_type = "reactivities"))

    def treated(self, message = None):
        self.ui.setScene(CotransTarget(self.ui, (self.name, self.seq, self.target_id), data_type = "treated_counts"))

    def untreated(self, message = None):
        self.ui.setScene(CotransTarget(self.ui, (self.name, self.seq, self.target_id), data_type = "untreated_counts"))

    def togglePlotType(self, message = None):
        self.plot_type = ("column" if self.plot_type == "row" else "row")
        button = self.buttonWithKey('togglePlotType')
        self.sendViewMessage(button, "setText", "Plot Type: {}".format(self.plot_type.capitalize()))

    def handleKeyEvent(self, keyInfo):
        handler = { "b" : self.beta,
                    "h" : self.theta,
                    "r" : self.reactivity,
                    "t" : self.treated,
                    "p" : self.togglePlotType,
                    "u" : self.untreated }.get(keyInfo["t"])
        if handler:
            handler()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

    def show_plot(self, site):
        if self.plot_type == "row":
            profiles = self.profiles.profilesForTargetAndEnd(self.name, site.end)
            if self.data_type == "treated_counts" or self.data_type == "untreated_counts":
                plot = { "type" : "Treated/Untreated Counts, length = {}".format(site.end),
                         "data" : [ { "label" : "f+", "x" : range(site.end + 1), "y" : profiles.treated_counts, "m" : "r-" },
                                    { "label" : "f-", "x" : range(site.end + 1), "y" : profiles.untreated_counts, "m" : "b-" } ],
                         "x_axis" : "Site",
                         "y_axis" : "% of stops" }
            else:
                data = getattr(profiles, self.data_type)
                data = data[1:] # exclude 0
                plot = { "type" : "{}, length = {}".format(self.data_type, site.end),
                         "data" : [ { "x" : range(1, site.end + 1), "y" : data, "m" : "-" } ],
                         "x_axis" : "Site",
                         "y_axis" : self.data_type }
        else:
            plot_axis = []
            plot_data = []
            seq = self.seq
            n = len(seq)
            for end in range(self.ui.spats.run.cotrans_minimum_length, n + 1):
                if site.site <= end:
                    profiles = self.profiles.profilesForTargetAndEnd(self.name, end)
                    data = getattr(profiles, self.data_type)
                    plot_axis.append(end)
                    plot_data.append(data[site.site])
            plot = { "type" : "NT {}: {}".format(site.site, self.data_type),
                     "data" : [ { "x" : plot_axis, "y" : plot_data, "m" : "g-", "label": "Site {}".format(site.site) } ],
                     "x_axis" : "Length",
                     "y_axis" : self.data_type }
        self.ui.plotter.submit_plot(plot)
