
import cjb.uif

strings = {
    'home' : 'Home',
    'newInterface' : 'New Interface',
    'newInterfacePrompt' : 'New Interface Name:',
    'newMediator' : 'New Mediator',
    'relationships' : 'Relationships',
}

class Localizer(cjb.uif.Localizer):

    def __init__(self):
        cjb.uif.Localizer.__init__(self)
        self.addStringDictionary(strings)

