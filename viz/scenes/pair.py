
import cjb.uif
from cjb.uif.layout import Size, Grid
from viz.scenes.base import BaseScene
from viz.layout import buttonSize


class Nuc(object):

    def __init__(self, char, context):
        self.char = char
        self.context = context

    @property
    def displayName(self):
        return self.char

TAG_COLORS = [
    [ 0.7, 0.7, 0.7 ],
    [ 0.9, 0.6, 0.6 ],
    [ 0.6, 0.6, 0.9 ],
    [ 0.9, 0.4, 0.1 ],
    [ 0.7, 0.7, 1.0 ],
]

class PairScene(BaseScene):

    def __init__(self, ui, pair):
        self.pair = pair
        self.parts = {}
        BaseScene.__init__(self, ui, self.__class__.__name__)

    def addNucView(self, nuc, bg, part_name):
        v = cjb.uif.views.Button(obj = nuc)
        v.sideSpacing = 0
        v.bg = bg
        self.addView(v)
        if not self.parts.get(part_name):
            self.parts[part_name] = []
        self.parts[part_name].append(v)

    def build(self):
        BaseScene.build(self)
        for part_name in ( "r1", "r2" ):
            part = getattr(self.pair, part_name)
            seq = part.original_seq
            color_idx = 1
            color = TAG_COLORS[color_idx]
            idx = 0
            for tag in part.tags:
                while idx < tag[1]:
                    self.addNucView(Nuc(seq[idx], (self.pair, part_name, idx)), TAG_COLORS[0], part_name)
                    idx += 1
                while idx < tag[1] + tag[2]:
                    self.addNucView(Nuc(seq[idx], (self.pair, part_name, idx)), color, part_name)
                    idx += 1
                color_idx += 1
                color = TAG_COLORS[color_idx % len(TAG_COLORS)]
            while idx < len(seq):
                self.addNucView(Nuc(seq[idx], (self.pair, part_name, idx)), TAG_COLORS[0], part_name)
                idx += 1

    def layout(self, view):
        BaseScene.layout(self, view)
        frame = view.frame.centeredSubrect(600, 200)
        grid = Grid(frame = frame.topSubrect(100), itemSize = Size(12, 18), columns = 35, rows = 1)
        grid.applyToViews(self.parts["r1"])
        grid = Grid(frame = frame.leftover, itemSize = Size(12, 18), columns = 35, rows = 1)
        grid.applyToViews(self.parts["r2"])
        return view
