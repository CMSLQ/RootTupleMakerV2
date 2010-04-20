import FWCore.ParameterSet.Config as cms

rootTupleMuons = cms.EDProducer("RootTupleMakerV2_Trigger",
    InputTagL1  = cms.InputTag('gtDigis'),
    InputTagHLT = cms.InputTag('TriggerResults','','HLT'),
    HLTPathsOfInterest = cms.vstring('')
)
