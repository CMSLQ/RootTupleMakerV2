import FWCore.ParameterSet.Config as cms

rootTupleEventSelection = cms.EDProducer("RootTupleMakerV2_EventSelection",
    L1InputTag  = cms.InputTag('gtDigis'),
    HcalNoiseInputTag = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
    FilterResultsInputTag = cms.InputTag('TriggerResults','','PAT')
)
