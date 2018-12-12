import FWCore.ParameterSet.Config as cms

rootTupleEventSelection = cms.EDProducer("RootTupleMakerV2_EventSelection",
    L1InputTag  = cms.InputTag('gtDigis'),
    FilterResultsInputTag = cms.InputTag('TriggerResults','','PAT'),
    # process label can sometimes be RECO and sometimes PAT
)

