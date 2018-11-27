import FWCore.ParameterSet.Config as cms

rootTupleEventSelection = cms.EDProducer("RootTupleMakerV2_EventSelection",
    L1InputTag  = cms.InputTag('gtDigis'),
    FilterResultsInputTag = cms.InputTag('TriggerResults','','HLT'),
    # process label can sometimes be RECO and sometimes PAT
    # just get most recently produced collection
    # https://hypernews.cern.ch/HyperNews/CMS/get/physTools/3396/1/1/2.html
    BadPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
    BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter")
)

from RecoMET.METFilters.BadPFMuonFilter_cfi import *
from RecoMET.METFilters.BadChargedCandidateFilter_cfi import *
BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

rootTupleEventSelectionSequence = cms.Sequence(
  BadPFMuonFilter*
  BadChargedCandidateFilter*
  rootTupleEventSelection
)

