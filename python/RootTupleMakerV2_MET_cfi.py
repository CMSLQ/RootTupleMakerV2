import FWCore.ParameterSet.Config as cms

rootTupleCaloMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETs'),
    Prefix = cms.string('Calo'),
    Suffix = cms.string('')
)

rootTupleTCMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsTC'),
    Prefix = cms.string('TC'),
    Suffix = cms.string('')
)

rootTuplePFMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('')
)
