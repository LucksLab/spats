
import cjb.uif

from cjb.uif.views import Label
from viz.layout import buttonSize


class BaseScene(cjb.uif.Scene):

    def __init__(self, ui, key = None):
        self.ui = ui
        self.scroller = None
        cjb.uif.Scene.__init__(self, ui.manager, key or self.__class__.__name__)
        self.container.properties['sendKeys'] = 1

    def build(self):
        # common UI
        if self.key != 'Home':
            self.targetButtons([self.ui.home, self.back])

    def layout(self, view):
        # common layout
        home = self.buttonWithKey('home')
        if home:
            home.frame = view.frame.bottomRightSubrect(size = buttonSize, margin = 10)
        back = self.buttonWithKey('back')
        if back:
            back.frame = view.frame.topRightSubrect(size = buttonSize, margin = 10)
        return view

    def back(self, message = None):
        self.ui.popScene()

    def addLabel(self, txt, bg = None):
        return self.addView(Label(txt, fontSize = 11, bg = bg))

    def addModelView(self, obj):
        return self.addView(cjb.uif.views.Button(obj = obj))

    def addModelViews(self, objs):
        map(self.addModelView, objs)

    def handleViewMessage(self, scene, obj, message):
        if obj:
            if isinstance(obj, Relationship):
                self.showRelationship(obj)
            #...
            else:
                print "App got message to " + str(obj) + ": " + str(message)
        elif message.get('event') == 'key':
            self.handleKeyEvent(message["arg"])
        else:
            print "App got general message: " + str(message)

    def handleKeyEvent(self, keyInfo):
        if keyInfo["t"] == "h" and 1 == len(keyInfo):
            self.ui.home()
        elif keyInfo["t"] == "b" and 1 == len(keyInfo):
            self.back()
        else:
            print "Unhandled key: " + str(keyInfo)
