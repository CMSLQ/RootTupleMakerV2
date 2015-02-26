import FWCore.ParameterSet.Config as cms

rootTuplePFJets = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('slimmedJets'),
    ## InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
    #InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFresUp'),                                 
    #InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFresDown'),                                 
    #InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFenUpForCorrMEt'),                                 
    #InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFenDownForCorrMEt'),                                 
    Prefix = cms.string('PFJet'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(30),
    JECUncertainty = cms.string('AK5PF'),
    ReadJECuncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices')
)
