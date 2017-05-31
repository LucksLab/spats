
import cjb.uif

class Colorize(cjb.uif.Filter):

    def filterView(self, view, scene):
        view.bg = [ 0.8, 0.85, 0.8 ]
        for v in view.recursiveSubviews():
            color = None
            typeId = getattr(v.obj, 'typeId', None)
            if typeId:
                color = { 'rel' : [ 0.9, 0.6, 0.6 ],
                          'intf' : [ 0.6, 0.6, 0.9 ],
                          'med' : [ 0.9, 0.4, 0.1 ],
                          'method' : [ 0.7, 0.7, 1.0 ] }.get(typeId)
            elif isinstance(v, cjb.uif.views.Button):
                color = [ 0.7, 0.7, 0.7 ]
            if color:
                v.bg = color
        return view
