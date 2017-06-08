
from cjb.uif.layout import Size, Grid
from viz.layout import buttonSize
from viz.scenes import BaseScene, Tagset, Targets


class Home(BaseScene):

    def targets(self, message = None):
        self.ui.pushScene(Targets(self.ui))

    def reads(self, message = None):
        self.ui.pushScene(Tagset(self.ui))

    def build(self):
        BaseScene.build(self)
        self.targetButtons([self.reads, self.targets])

    def handleKeyEvent(self, keyInfo):
        handler = { "r" : self.reads, "t" : self.targets }.get(keyInfo["t"])
        if handler:
            handler()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

    def layout(self, view):
        BaseScene.layout(self, view)
        grid = Grid(frame = view.frame, itemSize = buttonSize, columns = 2, rows = 2, spacing = Size(100, 300))
        grid.applyToViews([ v for v in view.subviews if getattr(v, "key", None) != 'home' ])
        return view
