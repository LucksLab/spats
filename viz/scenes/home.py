
from cjb.uif.layout import Size, Grid
from cjb.uif.views import Label
from viz.layout import buttonSize
from viz.scenes import BaseScene, Tagset, Targets


class Home(BaseScene):

    def targets(self, message = None):
        self.ui.pushScene(Targets(self.ui))

    def reads(self, message = None):
        self.ui.pushScene(Tagset(self.ui))

    def open_spats(self, message = None):
        self.sendViewMessage(self.container, "openFile", { "allow_multiple" : "0", "ext" : "spats", "allow_folders" : "0" })

    def handle_openFile(self, message):
        self.ui.open_spats(message["paths"][0])
        self.sendViewMessage(self.buttonWithKey('reads'), "show" if self.ui.has_reads_data() else "hide")
        self.sendViewMessage(self.buttonWithKey('targets'), "show" if self.ui.has_counter_data() else "hide")
        self.sendViewMessage(self.db_name, "setText", self.ui.db_name)

    def build(self):
        self.ui.try_load_last()
        BaseScene.build(self)
        self.targetButtons([self.open_spats, self.reads, self.targets])
        self.db_name = self.addView(Label(self.ui.db_name, fontSize = 12))
        if not self.ui.has_reads_data():
            self.buttonWithKey('reads').hidden = True
        if not self.ui.has_counter_data():
            self.buttonWithKey('targets').hidden = True

    def handleKeyEvent(self, keyInfo):
        handler = { "r" : self.reads, "t" : self.targets, "o" : self.open_spats }.get(keyInfo["t"])
        if handler:
            handler()
        else:
            BaseScene.handleKeyEvent(self, keyInfo)

    def handleViewMessage(self, scene, obj, message):
        if message.get("event") == "openFile":
            self.handle_openFile(message)
        else:
            BaseScene.handleViewMessage(self, scene, obj, message)

    def layout(self, view):
        BaseScene.layout(self, view)
        grid = Grid(frame = view.frame, itemSize = buttonSize, columns = 3, rows = 2, spacing = Size(100, 300))
        grid.applyToViews([ v for v in view.subviews if getattr(v, "key", None) != 'home' ])
        self.db_name.frame = view.frame.topLeftSubrect(200, 40, margin = 40);
        return view
