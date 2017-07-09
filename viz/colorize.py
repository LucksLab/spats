
import ast

import cjb.uif


class Colorize(cjb.uif.Filter):

    def __init__(self, custom_colors = None):
        self._colors = {
            "target" : [ 1.0, 0.85, 0.7 ],
            "adapter_t" : [ 1.0, 0.5, 0.5 ],
            "adapter_b" : [ 0.5, 0.5, 1.0 ],
            "linker_cotrans" : [ 0.7, 0.9, 1.0 ],
            "rrry" : [ 0.2, 1.0, 0.1 ],
            "yyyr" : [ 0.2, 0.7, 0.1 ],
            "error" : [ 0.9, 0.1, 0.2 ],
            "bar" : [ 0.8, 0.6, 1.0 ],
            "grey" : [ 0.7, 0.7, 0.7 ],
            "light_grey" : [ 0.8, 0.8, 0.8 ],
        }
        if custom_colors:
            for key in custom_colors.keys():
                self._colors[key] = ast.literal_eval(custom_colors[key])

    def color(self, key):
        return self._colors[key.lower()]

    def filterView(self, view, scene):
        view.bg = [ 0.8, 0.85, 0.8 ]
        for v in view.recursiveSubviews():
            color = None
            if v.bg and isinstance(v.bg, basestring):
                v.bg = self.color(v.bg)
            elif isinstance(v, cjb.uif.views.Button):
                color = [ 0.7, 0.7, 0.7 ]
            if not v.bg and color:
                v.bg = color
        return view
