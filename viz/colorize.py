
import cjb.uif

class Colorize(cjb.uif.Filter):

    def __init__(self):
        self._colors = {
            "target" : [ 1.0, 0.85, 0.7 ],
            "adapter_t" : [ 1.0, 0.5, 0.5 ],
            "adapter_b" : [ 0.5, 0.5, 1.0 ],
            "RRRY" : [ 0.2, 1.0, 0.1 ],
            "YYYR" : [ 0.2, 0.7, 0.1 ],
            "error" : [ 0.9, 0.1, 0.2 ],
            "bar" : [ 0.8, 0.6, 1.0 ],
            "grey" : [ 0.7, 0.7, 0.7 ],
            "light_grey" : [ 0.8, 0.8, 0.8 ],
        }

    def color(self, key):
        return self._colors[key]

    def filterView(self, view, scene):
        view.bg = [ 0.8, 0.85, 0.8 ]
        for v in view.recursiveSubviews():
            color = None
            if v.bg and isinstance(v.bg, basestring):
                v.bg = self._colors[v.bg]
            elif isinstance(v, cjb.uif.views.Button):
                color = [ 0.7, 0.7, 0.7 ]
            if not v.bg and color:
                v.bg = color
        return view
