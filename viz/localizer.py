
import cjb.uif

strings = {
    'home' : 'Home',
    'back' : 'Back',
    'reads' : 'Reads',
    'targets' : 'Targets',
    'showAll' : 'Show All',
}

class Localizer(cjb.uif.Localizer):

    def __init__(self):
        cjb.uif.Localizer.__init__(self)
        self.addStringDictionary(strings)

