
from spats_shape_seq.pair import Pair
from spats_shape_seq.target import _Target

from cjb.uif.layout import Size, Grid
from viz.layout import buttonSize
from viz.scenes import BaseScene, PairScene


class Home(BaseScene):

    def _pair(self, pairid):
        pair = Pair()

        if 0 == pairid:
            #pair.set_from_data("18333", "GAGTGTCCTTGGTGCCCGAGTCAGTGGTAGATCGG", "ACCACTGACTCGGGCACCAAGGACACTCAGATCGG")
            #pair.r1.tags = [('RRRY', 0, 4, 0), ('5s_rc', 4, 20, 0), ('5s', 24, 6, 76), ('adapter_b', 28, 7, 0)]
            #pair.r2.tags = [('5s', 4, 20, 123), ('RRRY', 24, 4, 0), ('5s_rc', 26, 6, 16), ('adapter_t_rc', 28, 7, 0)]
            pair.set_from_data("21189", "TTTGGTCCTTGGTGCCCGAGTCAGAGATCGGAAGA", "CTGACTCGGGCACCAAGGACCAAAAGATCGGAAGA")
            pair.r1.tags = [('YYYR', 0, 4, 0), ('5s_rc', 4, 21, 0), ('adapter_b', 25, 10, 1)]
            pair.r2.tags = [('5s', 0, 20, 123), ('YYYR', 20, 4, 0), ('adapter_t_rc', 24, 11, 0)]
        else:
            pair.set_from_data("1101:20069:1063", "TTTAGTCCTTGGTGCCCGAGTCAGATGCCTGGCAG", "TCCCACCTGACCCCATGCCGAACTCAGAAGTGAAA")
            pair.r1.tags = [('YYYR', 0, 4, 0), ('5s_rc', 4, 31, 0)]
            pair.r2.tags = [('5s', 0, 35, 27)]

        pair.target = _Target("5s", "GGATGCCTGGCGGCCGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCTGACTCGGGCACCAAGGAC", 0)
        return pair


    def test(self, message = None):
        self.ui.setScene(PairScene(self.ui, self._pair(0), expanded = False))

    def test2(self, message = None):
        self.ui.setScene(PairScene(self.ui, self._pair(1), expanded = False))

    def build(self):
        BaseScene.build(self)
        self.targetButtons([self.test, self.test2]) #, self.interfaces, self.mediators, self.newRelationship, self.newInterface, self.newMediator])

    def handleKeyEvent(self, keyInfo):
        handler = { "t" : self.test }.get(keyInfo["t"])
        if handler:
            handler()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

    def layout(self, view):
        grid = Grid(frame = view.frame, itemSize = buttonSize, columns = 3, rows = 2, spacing = Size(100, 300))
        grid.applyToViews([ v for v in view.subviews if getattr(v, "key", None) != 'home' ])
        return view
