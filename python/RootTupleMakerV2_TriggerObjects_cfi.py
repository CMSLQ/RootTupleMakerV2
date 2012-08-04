import FWCore.ParameterSet.Config as cms

rootTupleTriggerObjects= cms.EDProducer("RootTupleMakerV2_TriggerObjects",
    InputTag = cms.InputTag('hltTriggerSummaryAOD',"","HLT"),
    Prefix = cms.string("HLT"),
    Suffix = cms.string('')
)

