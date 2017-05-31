
from cjb.uif.layout import Size, Grid
from viz.scenes.base import BaseScene
from viz.layout import buttonSize


class Home(BaseScene):

    def relationships(self, message = None):
        print "Rel"

    def build(self):
        BaseScene.build(self)
        self.targetButtons([self.relationships]) #, self.interfaces, self.mediators, self.newRelationship, self.newInterface, self.newMediator])

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
