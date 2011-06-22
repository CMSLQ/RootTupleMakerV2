import FWCore.ParameterSet.Config as cms

rootTupleTaus = cms.EDProducer("RootTupleMakerV2_Taus",
    InputTag = cms.InputTag('cleanPatTaus'),
    Prefix = cms.string('Tau'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10)
)
