import FWCore.ParameterSet.Config as cms

rootTupleTriggerObjectsSingleEleIDISO = cms.EDProducer("RootTupleMakerV2_TriggerObjects",
    InputTag = cms.InputTag('hltTriggerSummaryAOD'),
    FilterID = cms.InputTag('hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter','','HLT'),
    Prefix = cms.string('HLTEle15IdIso'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10)
)

