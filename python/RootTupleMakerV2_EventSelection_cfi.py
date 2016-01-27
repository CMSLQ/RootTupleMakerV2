import FWCore.ParameterSet.Config as cms

rootTupleEventSelection = cms.EDProducer("RootTupleMakerV2_EventSelection",
    L1InputTag  = cms.InputTag('gtDigis'),
    HBHENoiseFilterResultInputTag = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
    HBHEIsoNoiseFilterResultInputTag = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
    FilterResultsInputTag = cms.InputTag('TriggerResults','','')
    # process label can sometimes be RECO and sometimes PAT
    # just get most recently produced collection
    # https://hypernews.cern.ch/HyperNews/CMS/get/physTools/3396/1/1/2.html
)
