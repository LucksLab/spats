
import math
import cjb.uif

from cjb.uif.layout import Size, Grid, Rect, layoutInScroller
from cjb.uif.views import Button, Label, TextField, View
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
            self.ui.pushScene(CotransTarget(self.ui, obj))
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


class CotransTarget(BaseScene):

    def __init__(self, ui, target, data_type = "treated"):
        self.cotrans = ui.spats.run.cotrans
        self.name = target[0]
        self.seq = target[1]
        self.target_id = target[2]
        self.data_type = data_type
        self.plot_type = "row"
        self.selected_values = None
        self.show_counts = True
        self.data_types = [ "treated", "untreated", "beta", "theta", "rho" ]
        for data_type in self.data_types:
            setattr(self, data_type, self.changer(data_type))
            setattr(self, "show_" + data_type, self.shower(data_type))
            setattr(self, "update_" + data_type, self.updater(data_type))
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def changer(self, data_type):
        return lambda message = None: self.change_plot(data_type)

    def shower(self, data_type):
        return lambda message = None: self.show_max(data_type)

    def updater(self, data_type):
        return lambda message = None: self.update_plot_max(data_type, message)

    def show_max(self, data_type):
        tf = self.data_type_views[self.data_types.index(data_type)].textField
        self.sendViewMessage(tf, "show")
        self.sendViewMessage(tf, "focus")

    def update_plot_max(self, data_type, message):
        self.data_type = data_type
        dtv = self.data_type_views[self.data_types.index(data_type)]
        max_val = float(message["arg"]) or self.profiles.data_range(data_type)[1]
        self.maxes[self.data_type] = max_val
        self.sendViewMessage(dtv.textField, "hide")
        self.sendViewMessage(dtv.label, "setText", "max: {}".format(float(max_val)))
        self.sendViewMessage(self.matrix, "matrix_plot", { "plot" : data_type, "max" : max_val })

    def makeDataTypeView(self, data_type):
        dtv = self.addView(View())
        dtv.button = dtv.addSubview(Button(target = getattr(self, data_type), text = data_type))
        self.addView(dtv.button)
        dtv.label = dtv.addSubview(Label("max: {}".format(self.maxes[data_type]), fontSize = 11))
        dtv.label.click = 1
        dtv.label.target = getattr(self, "show_" + data_type)
        self.addView(dtv.label)
        dtv.textField = dtv.addSubview(TextField(target = getattr(self, "update_" + data_type)))
        dtv.textField.hidden = True
        self.addView(dtv.textField)
        return dtv

    def build(self):
        BaseScene.build(self)
        self.profiles = self.ui.spats.compute_profiles()
        self.maxes = { dt : (4 if dt == "rho" else self.profiles.data_range(dt)[1]) for dt in self.data_types }
        self.data_type_views = [ self.makeDataTypeView(data_type) for data_type in self.data_types ]
        self.targetButtons([self.togglePlotType, self.toggleCounts, self.totalCounts])
        self.buttonWithKey('togglePlotType').text = "Plot Type: Row"
        self.buttonWithKey('toggleCounts').text = "Show f+/f-: On"
        self.type_label = self.addView(cjb.uif.views.Label("Query: {}".format(self.data_type), fontSize = 16))
        if self.ui.has_tags:
            sitemap = { "{}:{}:{}:{}".format(self.target_id, s[0], s[2], s[1]) : s[3] for s in self.ui.db.result_sites(self.ui.result_set_id, self.target_id) }
        else:
            sitemap = self.ui.spats.counters.registered_dict()
        self.matrix = cjb.uif.views.CustomView('SpatsMatrix')
        self.matrix.properties["d"] = self.profiles.cotrans_data()
        self.matrix.properties["plot"] = self.data_type
        self.matrix.properties["max"] = self.profiles.data_range(self.data_type)[1]
        self.addView(self.matrix)

    def layout(self, view):
        BaseScene.layout(self, view)
        cur = view.frame.centeredSubrect(w = 800, h = 600)
        self.matrix.frame = cur
        self.type_label.frame = view.frame.topLeftSubrect(w = 240, h = 24, margins = Size(40, 100))
        grid = Grid(frame = view.frame.bottomCenteredSubrect(w = 600, h = 40, margin = 20), itemSize = Size(180, 40), columns = 3, rows = 1, spacing = Size(20, 0))
        grid.applyToViews(map(self.buttonWithKey, [ 'togglePlotType', 'toggleCounts', 'totalCounts' ]))
        grid = Grid(frame = view.frame.leftCenteredSubrect(w = 120, h = 400, margin = 40), itemSize = Size(120, 70), columns = 1, rows = 5, spacing = Size(0, 10))
        grid.applyToViews(self.data_type_views)
        for dtv in self.data_type_views:
            self.layoutDataTypeView(dtv)
        return view

    def layoutDataTypeView(self, dtv):
        f = dtv.frame.bounds()
        dtv.button.frame = f.topSubrect(40)
        f = f.leftover.topSubrect(24, 2)
        dtv.label.frame = f
        dtv.textField.frame = f
        return dtv

    def handleViewMessage(self, scene, obj, message):
        if "site" == message.get("event"):
            self.show_plot(message["L"], message["site"])
        elif "sel" == message.get("event"):
            self.selected_values = sorted(message["data"])
        elif "show_sel" == message.get("event"):
            self.show_selected_plot()
        elif obj and isinstance(obj, Site):
            self.show_plot_site(obj)
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

    def change_plot(self, data_type):
        self.data_type = data_type
        self.sendViewMessage(self.matrix, "matrix_plot", { "plot" : data_type, "max" : self.maxes[data_type] })
        self.sendViewMessage(self.type_label, "setText", "Query: {}".format(self.data_type))


    def togglePlotType(self, message = None):
        self.plot_type = ("column" if self.plot_type == "row" else "row")
        button = self.buttonWithKey('togglePlotType')
        self.sendViewMessage(button, "setText", "Plot Type: {}".format(self.plot_type.capitalize()))
        self.sendViewMessage(self.matrix, self.plot_type)

    def toggleCounts(self, message = None):
        self.show_counts = not self.show_counts
        button = self.buttonWithKey('toggleCounts')
        self.sendViewMessage(button, "setText", "Show f+/f-: {}".format("On" if self.show_counts else "Off"))

    def totalCounts(self, message = None):
        self.ui.plotter.submit_plot(self.total_reads_plot())

    def handleKeyEvent(self, keyInfo):
        handler = None
        if "t" in keyInfo and keyInfo.get('cmd'):
            handler = { "b" : self.beta,
                        "h" : self.theta,
                        "r" : self.rho,
                        "t" : self.treated,
                        "u" : self.untreated }.get(keyInfo["t"])
        elif "t" in keyInfo and keyInfo.get('ctrl'):
            handler = { "b" : self.show_beta,
                        "h" : self.show_theta,
                        "r" : self.show_rho,
                        "t" : self.show_treated,
                        "u" : self.show_untreated }.get(keyInfo["t"])
        elif "t" in keyInfo:
            handler = { "p" : self.togglePlotType,
                        "c" : self.toggleCounts,
                        "t" : self.totalCounts }.get(keyInfo["t"])
        if handler:
            handler()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

    def total_reads_plot(self):
        n = len(self.seq)
        min_length = self.ui.spats.run.cotrans_minimum_length
        treated = []
        untreated = []
        treated_sum = 0.0
        untreated_sum = 0.0
        for end in range(min_length, n + 1):
            profiles = self.profiles.profilesForTargetAndEnd(self.name, end)
            t = sum(profiles.treated)
            treated.append(t)
            treated_sum += t
            u = sum(profiles.untreated)
            untreated.append(u)
            untreated_sum += u
        return { "type" : "Total Treated/Untreated Counts",
                 "data" : [ { "label" : "f+", "x" : range(min_length, n + 1), "y" : treated, "m" : "r-" },
                            { "label" : "f-", "x" : range(min_length, n + 1), "y" : untreated, "m" : "b-" } ],
                 "xlim" : [ min_length, n + 1],
                 "x_axis" : "Length",
                 "y_axis" : "# of stops" }

    def count_plot(self, profiles, L):
        treated_sum = sum(profiles.treated_counts)
        untreated_sum = sum(profiles.untreated_counts)
        treated_pcts = map(lambda x: (float(x) / float(treated_sum)), profiles.treated_counts)
        untreated_pcts = map(lambda x: (float(x) / float(untreated_sum)), profiles.untreated_counts)
        return { "type" : "Treated/Untreated Counts, length = {}".format(L),
                 "data" : [ { "label" : "f+ ({})".format(int(treated_sum)), "x" : range(L + 1), "y" : treated_pcts, "m" : "r-" },
                            { "label" : "f- ({})".format(int(untreated_sum)), "x" : range(L + 1), "y" : untreated_pcts, "m" : "b-" } ],
                 "x_axis" : "Site",
                 "y_axis" : "% of stops" }

    def show_plot_site(self, site):
        self.show_plot(site.end, site.site)

    def show_plot(self, L, site):
        add_counts = False
        filname = None
        if self.plot_type == "row":
            profiles = self.profiles.profilesForTargetAndEnd(self.name, L)
            if self.data_type == "treated" or self.data_type == "untreated":
                plot = self.count_plot(profiles, L)
                filename = "counts_L{}".format(L)
            else:
                data = getattr(profiles, self.data_type)
                data = data[1:] # exclude 0
                plot = { "type" : "{}, length = {}".format(self.data_type, L),
                         "data" : [ { "x" : range(1, L + 1), "y" : data, "m" : "-" } ],
                         "ylim" : [ 0, self.maxes[self.data_type] ],
                         "x_axis" : "Site",
                         "y_axis" : self.data_type }
                filename = "{}_L{}".format(self.data_type, L)
                add_counts = self.show_counts
        else:
            plot_axis = []
            plot_data = []
            n = len(self.seq)
            for end in range(self.ui.spats.run.cotrans_minimum_length, n + 1):
                if site <= end:
                    profiles = self.profiles.profilesForTargetAndEnd(self.name, end)
                    data = getattr(profiles, self.data_type)
                    plot_axis.append(end)
                    plot_data.append(data[site])
            plot = { "type" : "NT {}: {}".format(site, self.data_type),
                     "data" : [ { "x" : plot_axis, "y" : plot_data, "m" : "g-", "label": "Site {}".format(site) } ],
                     "ylim" : [ 0, self.maxes[self.data_type] ],
                     "x_axis" : "Length",
                     "y_axis" : self.data_type }
            filename = "{}_nt{}".format(self.data_type, site)
        filename = "{}_{}".format(self.name, filename)
        if add_counts:
            self.ui.plotter.submit_plots([plot, self.count_plot(profiles, L)], filename)
        else:
            self.ui.plotter.submit_plot(plot, filename)


    def show_selected_plot(self):
        if not self.selected_values:
            return
        count_plots = []
        filname = None
        if self.plot_type == "row":
            plot_data = []
            useAll = (self.data_type == "treated" or self.data_type == "untreated")
            for L in self.selected_values:
                profiles = self.profiles.profilesForTargetAndEnd(self.name, L)
                data = getattr(profiles, self.data_type)
                if not useAll:
                    data = data[1:] # exclude 0
                plot_data.append({ "x" : range(0 if useAll else 1, L + 1), "y" : data, "m" : "-", "label" : "L={}".format(L) })
                if len(count_plots) < 3:
                    count_plots.append(self.count_plot(profiles, L))
            plot = { "type" : "{}, length = {}".format(self.data_type, self.selected_values),
                     "data" : plot_data,
                     "ylim" : [ 0, self.maxes[self.data_type] ],
                     "x_axis" : "Site",
                     "y_axis" : self.data_type }
            filename = "{}_L{}".format(self.data_type, "_".join(map(str,self.selected_values)))
        else:
            plot_axis = { val : [] for val in self.selected_values }
            plot_data = { val : [] for val in self.selected_values }
            n = len(self.seq)
            for site in self.selected_values:
                for end in range(self.ui.spats.run.cotrans_minimum_length, n + 1):
                    profiles = self.profiles.profilesForTargetAndEnd(self.name, end)
                    data = getattr(profiles, self.data_type)
                    if site <= end:
                        plot_axis[site].append(end)
                        plot_data[site].append(data[site])
            plot = { "type" : "{}, NTs: {}".format(self.data_type, self.selected_values),
                     "data" : [ { "x" : plot_axis[s], "y" : plot_data[s], "m" : "-", "label": "S={}".format(s) } for s in self.selected_values ],
                     "ylim" : [ 0, self.maxes[self.data_type] ],
                     "x_axis" : "Length",
                     "y_axis" : self.data_type }
            filename = "{}_nt{}".format(self.data_type, "_".join(map(str,self.selected_values)))
        if count_plots:
            self.ui.plotter.submit_plots([plot] + count_plots, filename)
        else:
            self.ui.plotter.submit_plot(plot, filename)
