
import os

import cjb.uif
import cjb.util.cfg
import viz.scenes
import viz.colorize
import viz.localizer


class SpatsViz(cjb.uif.UIServer):

    def __init__(self):
        self.config = cjb.util.cfg.parse(os.path.expanduser("~/.spats_viz"))["ide"]
        cjb.uif.UIServer.__init__(self, self.config["host"], int(self.config["port"]), self.config["portfile"])
        self.addFilter(viz.colorize.Colorize())
        self.addFilter(viz.localizer.Localizer())

    def reloadModel(self):
        #parser = ide.parse.Parser(self.config["path"])
        #parser.parse()
        #self.model = parser.model
        pass

    def home(self, message = None):
        self.reloadModel()
        self.setScene(viz.scenes.Home(self))
