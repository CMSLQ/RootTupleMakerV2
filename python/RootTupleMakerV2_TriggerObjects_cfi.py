import FWCore.ParameterSet.Config as cms

rootTupleTriggerObjectsHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilter = cms.EDProducer("RootTupleMakerV2_TriggerObjects",
    InputTag = cms.InputTag('hltTriggerSummaryAOD',"","HLT"),
    FilterID = cms.InputTag("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilter","","HLT"),
    Prefix = cms.string("HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilter"),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10)
)

rootTupleTriggerObjectsHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter = cms.EDProducer("RootTupleMakerV2_TriggerObjects",
    InputTag = cms.InputTag('hltTriggerSummaryAOD',"","HLT"),
    FilterID = cms.InputTag("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter","","HLT"),
    Prefix = cms.string("HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter"),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10)
)

rootTupleTriggerObjectsHLTEle27WP80TrackIsoFilter = cms.EDProducer("RootTupleMakerV2_TriggerObjects",
    InputTag = cms.InputTag('hltTriggerSummaryAOD',"","HLT"),
    FilterID = cms.InputTag("hltEle27WP80TrackIsoFilter","","HLT"),
    Prefix = cms.string("HLTEle27WP80TrackIsoFilter"),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10)
)





