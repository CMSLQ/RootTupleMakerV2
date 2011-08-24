from FWCore.GuiBrowsers.ConfigToolBase import *

class AddPfMETType1Cor(ConfigToolBase):

    """ Add Type-I corrected pfMET collection to patEventContent
    """
    _label='addPfMETType1Cor'
    _defaultParameters=dicttypes.SortedKeysDict()

    def __init__(self):
        ConfigToolBase.__init__(self)
        self.addParameter(self._defaultParameters,'postfixLabel','PFType1Cor', '')
        self._parameters=copy.deepcopy(self._defaultParameters)
        self._comment = ''

    def getDefaultParameters(self):
        return self._defaultParameters

    def __call__(self,process,postfixLabel=None):
        if  postfixLabel is None:
            postfixLabel=self._defaultParameters['postfixLabel'].value
        self.setParameter('postfixLabel',postfixLabel)
        self.apply(process)

    def toolCode(self, process):
        postfixLabel=self._parameters['postfixLabel'].value


        ## add module as process to the default sequence
        def addAlso (label,value):
            existing = getattr(process, label)
            setattr( process, label+postfixLabel, value)
            process.patDefaultSequence.replace( existing, existing*value )

        ## clone and add a module as process to the
        ## default sequence
        def addClone(label,**replaceStatements):
            new = getattr(process, label).clone(**replaceStatements)
            addAlso(label, new)

        addClone('patMETs', metSource = cms.InputTag("metJESCorAK5PFJet"), addMuonCorrections = False)

        ## add new met collections output to the pat summary
        process.patCandidateSummary.candidates += [ cms.InputTag('patMETs'+postfixLabel) ]



class AddPfChargedMET(ConfigToolBase):

    """ Add pfChargedMET collection to patEventContent
    """
    _label='addPfChargedMET'
    _defaultParameters=dicttypes.SortedKeysDict()

    def __init__(self):
        ConfigToolBase.__init__(self)
        self.addParameter(self._defaultParameters,'postfixLabel','PFCharged', '')
        self._parameters=copy.deepcopy(self._defaultParameters)
        self._comment = ''

    def getDefaultParameters(self):
        return self._defaultParameters

    def __call__(self,process,postfixLabel=None):
        if  postfixLabel is None:
            postfixLabel=self._defaultParameters['postfixLabel'].value
        self.setParameter('postfixLabel',postfixLabel)
        self.apply(process)

    def toolCode(self, process):
        postfixLabel=self._parameters['postfixLabel'].value


        ## add module as process to the default sequence
        def addAlso (label,value):
            existing = getattr(process, label)
            setattr( process, label+postfixLabel, value)
            process.patDefaultSequence.replace( existing, existing*value )

        ## clone and add a module as process to the
        ## default sequence
        def addClone(label,**replaceStatements):
            new = getattr(process, label).clone(**replaceStatements)
            addAlso(label, new)

        addClone('patMETs', metSource = cms.InputTag("chargedPFMet"), addMuonCorrections = False)

        ## add new met collections output to the pat summary
        process.patCandidateSummary.candidates += [ cms.InputTag('patMETs'+postfixLabel) ]


##################################

addPfMETType1Cor=AddPfMETType1Cor()
addPfChargedMET=AddPfChargedMET()
