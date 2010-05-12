import FWCore.ParameterSet.Config as cms

rootTupleGenJets = cms.EDProducer("RootTupleMakerV2_GenJets",
    InputTag = cms.InputTag('ak5GenJets'),
    Prefix = cms.string('GenJet'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10)
)

