
import cjb.uif

strings = {
    'home' : 'Home',
    'back' : 'Back',
    'reads' : 'Reads',
    'targets' : 'Targets',
    'showMatches' : 'Show Matches',
    'dump' : 'Dump',
    'open_spats' : 'Open',
    'betas' : 'beta',
    'thetas' : 'theta',
    'reactivities' : 'rho',
}

class Localizer(cjb.uif.Localizer):

    def __init__(self):
        cjb.uif.Localizer.__init__(self)
        self.addStringDictionary(strings)

