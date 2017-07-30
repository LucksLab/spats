
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

    def __init__(self, ui, target, data_type = "treated"):
        self.name = target[0]
        self.seq = target[1]
        self.target_id = target[2]
        self.data_type = data_type
        self.plot_type = "row"
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def build(self):
        BaseScene.build(self)
        self.targetButtons([self.treated, self.untreated, self.beta, self.theta, self.rho, self.togglePlotType])
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
        self.matrix = cjb.uif.views.CustomView('SpatsMatrix')
        self.matrix.properties["d"] = self.profiles.cotrans_data()
        self.matrix.properties["plot"] = self.data_type
        self.matrix.properties["max"] = 5
        self.addView(self.matrix)

    def layout(self, view):
        BaseScene.layout(self, view)
        cur = view.frame.centeredSubrect(w = 800, h = 600)
        self.matrix.frame = cur
        self.type_label.frame = view.frame.topLeftSubrect(w = 240, h = 24, margins = Size(40, 100))
        self.buttonWithKey('togglePlotType').frame = view.frame.bottomCenteredSubrect(200, 40, margin = 24)
        grid = Grid(frame = view.frame.leftCenteredSubrect(w = 120, h = 400, margin = 40), itemSize = Size(120, 40), columns = 1, rows = 5, spacing = Size(0, 40))
        grid.applyToViews(map(self.buttonWithKey, [ 'treated', 'untreated', 'beta', 'theta', 'rho' ]))
        return view

    def handleViewMessage(self, scene, obj, message):
        if "site" == message.get("event"):
            self.show_plot(message["L"], message["site"])
        elif obj and isinstance(obj, Site):
            self.show_plot_site(obj)
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

    def change_plot(self, plot):
        self.data_type = plot
        if plot == "rho":
            # TODO
            max_val = 4
        else:
            max_val = 0
            n = len(self.seq)
            for end in range(self.ui.spats.run.cotrans_minimum_length, n + 1):
                profiles = self.profiles.profilesForTargetAndEnd(self.name, end)
                data = getattr(profiles, plot)
                max_val = max(max_val, max(data))
        self.sendViewMessage(self.matrix, "matrix_plot", { "plot" : plot, "max" : max_val })

    def beta(self, message = None):
        self.change_plot("beta")

    def theta(self, message = None):
        self.change_plot("theta")

    def rho(self, message = None):
        self.change_plot("rho")

    def treated(self, message = None):
        self.change_plot("treated")

    def untreated(self, message = None):
        self.change_plot("untreated")

    def togglePlotType(self, message = None):
        self.plot_type = ("column" if self.plot_type == "row" else "row")
        button = self.buttonWithKey('togglePlotType')
        self.sendViewMessage(button, "setText", "Plot Type: {}".format(self.plot_type.capitalize()))

    def handleKeyEvent(self, keyInfo):
        handler = None
        if "t" in keyInfo and keyInfo.get('ctrl'):
            handler = { "b" : self.beta,
                        "h" : self.theta,
                        "r" : self.rho,
                        "t" : self.treated,
                        "p" : self.togglePlotType,
                        "u" : self.untreated }.get(keyInfo["t"])
        if handler:
            handler()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

    def count_plot(self, profiles, L, site):
        return { "type" : "Treated/Untreated Counts, length = {}".format(L),
                 "data" : [ { "label" : "f+", "x" : range(L + 1), "y" : profiles.treated_counts, "m" : "r-" },
                            { "label" : "f-", "x" : range(L + 1), "y" : profiles.untreated_counts, "m" : "b-" } ],
                 "x_axis" : "Site",
                 "y_axis" : "% of stops" }

    def show_plot_site(self, site):
        self.show_plot(site.end, site.site)

    def show_plot(self, L, site):
        add_counts = False
        if self.plot_type == "row":
            profiles = self.profiles.profilesForTargetAndEnd(self.name, L)
            if self.data_type == "treated" or self.data_type == "untreated":
                plot = self.count_plot(profiles, L, site)
            else:
                data = getattr(profiles, self.data_type)
                data = data[1:] # exclude 0
                plot = { "type" : "{}, length = {}".format(self.data_type, L),
                         "data" : [ { "x" : range(1, L + 1), "y" : data, "m" : "-" } ],
                         "x_axis" : "Site",
                         "y_axis" : self.data_type }
                add_counts = True
        else:
            plot_axis = []
            plot_data = []
            seq = self.seq
            n = len(seq)
            for end in range(self.ui.spats.run.cotrans_minimum_length, n + 1):
                if site <= end:
                    profiles = self.profiles.profilesForTargetAndEnd(self.name, end)
                    data = getattr(profiles, self.data_type)
                    plot_axis.append(end)
                    plot_data.append(data[site])
            plot = { "type" : "NT {}: {}".format(site, self.data_type),
                     "data" : [ { "x" : plot_axis, "y" : plot_data, "m" : "g-", "label": "Site {}".format(site) } ],
                     "x_axis" : "Length",
                     "y_axis" : self.data_type }
        if add_counts:
            self.ui.plotter.submit_plots([plot, self.count_plot(profiles, L, site)])
        else:
            self.ui.plotter.submit_plot(plot)
