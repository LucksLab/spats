
import cjb.uif

from viz.layout import buttonSize


class BaseScene(cjb.uif.Scene):

    def __init__(self, ui, key = None):
        self.ui = ui
        cjb.uif.Scene.__init__(self, ui.manager, key or self.__class__.__name__)
        self.container.properties['sendKeys'] = 1

    def build(self):
        # common UI
        if self.key != 'Home':
            self.targetButton(self.ui.home)

    def layout(self, view):
        # common layout
        home = self.buttonWithKey('home')
        if home:
            home.frame = view.frame.bottomRightSubrect(size = buttonSize, margin = 20)
        return view

    def addModelView(self, obj):
        self.addView(cjb.uif.views.Button(obj = obj))

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
        else:
            print "Unhandled key: " + str(keyInfo)
