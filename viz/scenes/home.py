
from spats_shape_seq.pair import Pair

from cjb.uif.layout import Size, Grid
from viz.layout import buttonSize
from viz.scenes import BaseScene, PairScene


class Home(BaseScene):

    def test(self, message = None):
        pair = Pair()
        pair.set_from_data("21189", "TTTGGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACCAAAAGATCGGAAGA")
        pair.r1.tags = [('YYYR', 0, 4, 0), ('rc(5s)', 4, 21, 0), ('adapter_b', 25, 10, 1)]
        pair.r2.tags = [('5s', 0, 20, 123), ('adapter_b', 24, 11, 0)]
        self.ui.setScene(PairScene(self.ui, pair))

    def build(self):
        BaseScene.build(self)
        self.targetButtons([self.test]) #, self.interfaces, self.mediators, self.newRelationship, self.newInterface, self.newMediator])

    def handleKeyEvent(self, keyInfo):
        handlers = { "r" : [ self.relationships, self.newRelationship],
                     "i" : [ self.interfaces, self.newInterface],
                     "m" : [ self.mediators, self.newMediator]}.get(keyInfo["t"])
        if handlers:
            handlers[1 if "shift" in keyInfo else 0]()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

    def layout(self, view):
        grid = Grid(frame = view.frame, itemSize = buttonSize, columns = 3, rows = 2, spacing = Size(100, 300))
        grid.applyToViews([ v for v in view.subviews if getattr(v, "key", None) != 'home' ])
        return view
