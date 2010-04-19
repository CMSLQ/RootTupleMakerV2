import FWCore.ParameterSet.Config as cms

rootTupleCaloMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETs'),
    Prefix = cms.string('Calo'),
    Suffix = cms.string('')
)

rootTupleTcMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsTC'),
    Prefix = cms.string('Tc'),
    Suffix = cms.string('')
)

rootTuplePfMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsPF'),
    Prefix = cms.string('Pf'),
    Suffix = cms.string('')
)
